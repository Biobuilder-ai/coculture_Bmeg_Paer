#!/usr/bin/env python3
"""COMETS 双菌共培养 cross-feeding 图.

用法（丢一个 log 进去就行，其余文件按 run id 自动找）:

    python crossfeeding_figure.py fluxlog_0x7371f0403640
    python crossfeeding_figure.py fluxlog_0x7371f0403640 -o crossfeed.pdf --top 12
    python crossfeeding_figure.py fluxlog_0x7371f0403640 --at-hour 6     # 见下

------------------------------------------------------------------------------
为什么必须是 flux log 而不是 media log
------------------------------------------------------------------------------
media log 只记录"某物质在格子里还剩多少 mmol"，是两个菌共同作用后的净结果。
如果 A 排 5、B 吃 5，media log 是一条水平线 —— 恰恰是最强的 cross-feeding，
在 media log 里完全看不见。要拿到"谁吃谁排、各多少"，必须开 flux log:

    params.set_param('writeFluxLog', True)
    params.set_param('FluxLogRate', 1)      # 5 也行，只是积分粗一点

------------------------------------------------------------------------------
脚本干了什么
------------------------------------------------------------------------------
1. 从传进来的文件名里抠出 COMETS 的 run id（那串 `0x....`），在同一个目录下
   找齐这次模拟的全套文件:
       .current_layout_0x...    -> 两个 .cmd 模型文件，以及它们在 flux log 里的顺序
       .current_global_0x...    -> FluxLogRate（flux log 每隔几个 cycle 记一次）
       .current_package_0x...   -> timeStep（每个 cycle 多少小时）
       totalbiomasslog_0x...    -> 每个 cycle 每个菌的生物量（g）
       medialog_0x...           -> 只用来做一次质量守恒自检，不参与画图
2. 解析 .cmd（cometspy 写给 COMETS 的模型文件）拿到:
       REACTION_NAMES      反应顺序 == flux log 里通量列的顺序
       EXCHANGE_REACTIONS  哪些列是交换反应
       SMATRIX             每个交换反应对应哪个胞外代谢物
3. 把通量积分成"这一整轮模拟总共交换了多少物质":

       amount = Σ_t  |flux(t)| * biomass(t) * Δt          [mmol]

   通量单位是 mmol/gDW/h，乘生物量(g)再乘时间(h)就是绝对量 mmol。
   uptake（flux<0）和 secretion（flux>0）分开累加 —— 同一个物质在不同时段
   可能先排后吃，两边都要留住，不能相减。
4. 按下面的规则筛物质（就是你说的"有交换或者都使用"）:
       cross-fed : 一个菌排出、另一个菌吃进    -> 真正的互喂
       shared    : 两个菌都在吃                -> 底物竞争
       exclusive : 只有一个菌碰它              -> 默认不画，--include-exclusive 打开
5. 三栏画图: 左菌 | 中间竖排代谢物 | 右菌。
   箭头指向菌 = 摄取，箭头背离菌 = 排出，线宽 ∝ 交换总量。

------------------------------------------------------------------------------
--at-hour: 画某一时刻的瞬时通量，而不是整轮累计量
------------------------------------------------------------------------------
默认画的是「这一整轮总共交换了多少物质」(mmol)。加上 --at-hour 6 就改画 t=6 h
那一刻的交换通量 (mmol/gDW/h) —— 也就是 COMETS 在那个 cycle 解出来的 FBA 通量
本身，不乘生物量也不乘时间。会自动取最接近的 cycle。

两者回答的是不同的问题:
    累计量   谁一共吃了多少 —— 决定最终生物量，图上大头是量多的物质
    瞬时通量 此刻谁在吃谁在排 —— 能看到「先排后吃」这种被累计量抹平的现象

挑时刻有讲究: t≈0 就是初始条件下的 FBA 解(可以和 notebook 预飞那格对照验证);
碳耗尽之后所有通量衰减到 0，图会空。指数期中段最有信息量。

注意分类(cross-fed / shared / ...)也跟着按那一刻重算，所以同一次模拟不同时刻
分类可能不同 —— 这正是瞬时图的用处，但图注里务必把时刻标清楚(脚本会自动标)。
和 media log 的质量守恒自检始终用累计量，不受 --at-hour 影响。

------------------------------------------------------------------------------
默认过滤
------------------------------------------------------------------------------
水、质子和一堆无机盐（O2/CO2/磷/硫/K/Ca/Mg/Fe/...）两个菌都在用，量还特别大，
不过滤掉的话前十名全是它们，有机物全被挤出去。所以默认扣掉，用
`--keep-inorganic` 可以全部留下，或者用 `--exclude` / `--only` 自己指定。
NH4 故意保留 —— 在这套培养基里 P. aeruginosa 排铵、B. megaterium 吃铵，
是真实的 cross-feeding，不是杂音。
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

# ----------------------------------------------------------------------------
# 外观：低饱和度的配色，投稿用不会太跳
# ----------------------------------------------------------------------------
C_SECRETION = "#C0553C"     # 排出（菌 -> 代谢物）
C_UPTAKE = "#2E6E8E"        # 摄取（代谢物 -> 菌）
C_ORG = ["#4C7A5D", "#7A5C9E"]           # 左菌、右菌各自的强调色
C_ORG_FILL = ["#EDF3EF", "#F2EFF6"]
C_PILL_CROSS_FILL, C_PILL_CROSS_EDGE = "#FBEAE1", "#C0553C"
C_PILL_SHARE_FILL, C_PILL_SHARE_EDGE = "#EEF1F3", "#7B8A93"
C_PILL_COSEC_FILL, C_PILL_COSEC_EDGE = "#F2F0E4", "#9A8F63"
C_PILL_EXCL_FILL, C_PILL_EXCL_EDGE = "#F7F7F5", "#B9B9B4"
C_TEXT, C_MUTED = "#1F2328", "#6B7178"

LW_MIN, LW_MAX = 0.8, 9.0   # 箭头线宽范围（pt）

# 默认扣掉的无机/大宗组分（NH4 不在其中，见上面说明）
INORGANIC = {
    "h2o_e", "h_e", "o2_e", "co2_e", "pi_e", "so4_e", "k_e", "na1_e", "cl_e",
    "ca2_e", "mg2_e", "fe2_e", "fe3_e", "mn2_e", "zn2_e", "cu2_e", "cobalt2_e",
    "mobd_e", "ni2_e", "sel_e", "slnt_e", "tungs_e", "cbl1_e", "no3_e", "no2_e",
}

# BiGG id -> 可读名。查不到就把 id 去掉 _e 直接用。
PRETTY = {
    "cit_e": "Citrate", "glu__L_e": "L-Glutamate", "asp__L_e": "L-Aspartate",
    "ser__L_e": "L-Serine", "glc__D_e": "D-Glucose", "sucr_e": "Sucrose",
    "fru_e": "D-Fructose", "sbt__D_e": "D-Sorbitol", "nh4_e": "Ammonium",
    "ac_e": "Acetate", "lac__D_e": "D-Lactate", "lac__L_e": "L-Lactate",
    "succ_e": "Succinate", "for_e": "Formate", "etoh_e": "Ethanol",
    "pyr_e": "Pyruvate", "akg_e": "2-Oxoglutarate", "mal__L_e": "L-Malate",
    "fum_e": "Fumarate", "co2_e": "CO2", "o2_e": "O2", "h2o_e": "H2O",
    "h_e": "H+", "pi_e": "Phosphate", "so4_e": "Sulfate", "ppi_e": "Diphosphate",
    "ala__L_e": "L-Alanine", "gly_e": "Glycine", "urea_e": "Urea",
    "23cgmp_e": "2',3'-cGMP", "2obut_e": "2-Oxobutanoate", "4abut_e": "GABA",
    "acald_e": "Acetaldehyde", "glyc_e": "Glycerol", "man_e": "D-Mannose",
    "gal_e": "D-Galactose", "xyl__D_e": "D-Xylose", "arab__L_e": "L-Arabinose",
    "mnl_e": "Mannitol", "tre_e": "Trehalose", "malt_e": "Maltose",
    "gthrd_e": "Glutathione", "gcald_e": "Glycolaldehyde", "glyclt_e": "Glycolate",
    "ade_e": "Adenine", "ins_e": "Inosine", "hxan_e": "Hypoxanthine",
    "thymd_e": "Thymidine", "ura_e": "Uracil", "orot_e": "Orotate",
}


# ============================================================================
# 一、找齐这次模拟的文件
# ============================================================================
def find_run_files(seed: Path) -> dict:
    """从任意一个 log 文件反推出这次 run 的全套文件。

    cometspy 给同一次模拟的所有文件加同一个后缀（对象内存地址），比如
    fluxlog_0x7371f0403640 / totalbiomasslog_0x7371f0403640 / .current_layout_0x...，
    所以给一个就能找到全部。
    """
    m = re.search(r"(0x[0-9a-fA-F]+)", seed.name)
    if not m:
        raise SystemExit(
            f"文件名里没有 COMETS 的 run id（形如 0x7371f0403640）: {seed.name}\n"
            "请直接传 fluxlog_0x... 这样的文件。"
        )
    rid, folder = m.group(1), seed.parent
    files = {"run_id": rid, "dir": folder}
    for key, prefix in [("flux", "fluxlog_"), ("biomass", "totalbiomasslog_"),
                        ("media", "medialog_"), ("layout", ".current_layout_"),
                        ("params", ".current_global_"), ("package", ".current_package_")]:
        path = folder / f"{prefix}{rid}"
        files[key] = path if path.is_file() else None
    if files["flux"] is None:
        raise SystemExit(
            f"没找到 {folder / ('fluxlog_' + rid)}。\n"
            "这次模拟没开 flux log，media log 无法区分是哪个菌在吃/排。\n"
            "在 params 里加上这两行再跑一次：\n"
            "    params.set_param('writeFluxLog', True)\n"
            "    params.set_param('FluxLogRate', 1)"
        )
    return files


def check_models_fresh(files: dict, model_files: list[str], force: bool) -> None:
    """确认 .cmd 模型文件是这次模拟当时的那一份。

    这是个真会踩的坑：cometspy 每跑一次都用同样的文件名重写 .cmd（./B_megaterium.cmd），
    所以后面再跑一次，就把之前那次 flux log 赖以解释的反应顺序覆盖掉了。反应个数往往
    还一样，于是通量列静悄悄地错位——图照画，全是错的。

    自检办法：同一个解里 EX_o2_e 和 O2tex/O2tpp 必须数值相等（氧气过了外膜多少就
    吸收多少）。对不上就说明列错位了。这里用更省事的时间戳版本：.cmd 比 flux log 新，
    就说明它是后来那次跑覆盖的。
    """
    if force:
        return
    flux_mtime = files["flux"].stat().st_mtime
    stale = [n for n in model_files
             if (files["dir"] / n).stat().st_mtime > flux_mtime + 60]
    if stale:
        raise SystemExit(
            f"模型文件 {', '.join(stale)} 比 {files['flux'].name} 还新。\n"
            "cometspy 每次模拟都会重写同名 .cmd，这几个已经被后来的运行覆盖了；\n"
            "反应顺序一变，flux log 的通量列就全部错位，画出来的图是错的。\n"
            "请重新跑一次（params 里带上 writeFluxLog=True），跑完立刻出图；\n"
            "或者把当时的 .cmd 找回来。确认没问题可以用 --force 跳过这个检查。"
        )


def read_keyvals(path: Path | None) -> dict:
    """读 `key = value` 形式的 .current_global / .current_package。"""
    out = {}
    if path and path.is_file():
        for line in path.read_text(errors="replace").splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
    return out


def read_layout_models(path: Path | None) -> list[str]:
    """layout 第一行 `model_file ./A.cmd ./B.cmd` —— 顺序就是 flux log 里的模型编号顺序。"""
    if not (path and path.is_file()):
        return []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("model_file"):
            return [Path(tok).name for tok in line.split()[1:]]
    return []


# ============================================================================
# 二、解析 .cmd 模型文件
# ============================================================================
_SECTION = re.compile(r"^([A-Z][A-Z_0-9]*)\b(.*)$")


def parse_cmd(path: Path) -> dict:
    """从 COMETS 的 .cmd 模型文件里拿出反应名、交换反应下标、交换反应对应的代谢物。

    .cmd 是分段的文本文件，段头顶格写、段体缩进、`//` 结束一段：
        SMATRIX  1349  1709        <- 稀疏计量矩阵，行=代谢物 列=反应（都是 1-based）
        METABOLITE_NAMES
        REACTION_NAMES             <- 顺序 == flux log 里通量列的顺序，这是全部关键
        EXCHANGE_REACTIONS         <- 一行下标，指向 REACTION_NAMES
    """
    sections, head, body = {}, None, []
    for raw in path.read_text(errors="replace").splitlines():
        if raw.strip() == "//":
            if head:
                sections[head] = body
            head, body = None, []
            continue
        m = _SECTION.match(raw)
        if m and not raw[:1].isspace():
            if head:
                sections[head] = body
            head, body = m.group(1), []
            if m.group(2).strip():
                body.append(m.group(2).strip())
            continue
        if head is not None and raw.strip():
            body.append(raw.strip())
    if head:
        sections[head] = body

    rxn_names = sections.get("REACTION_NAMES", [])
    met_names = sections.get("METABOLITE_NAMES", [])
    # EXCHANGE_REACTIONS 是一整行下标（1-based），转成 0-based
    exch = [int(t) - 1 for t in " ".join(sections.get("EXCHANGE_REACTIONS", [])).split()]

    # SMATRIX 第一行是 "1349 1709" 这种尺寸，跳过；其余是 "row col value"
    want, col2met = set(exch), {}
    for line in sections.get("SMATRIX", []):
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            row, col = int(parts[0]) - 1, int(parts[1]) - 1
        except ValueError:
            continue
        # 交换反应在 S 里只有一个非零项，就是它进出的那个胞外代谢物
        if col in want and col not in col2met and 0 <= row < len(met_names):
            col2met[col] = met_names[row]

    return {
        "name": path.stem,
        "n_reactions": len(rxn_names),
        # {反应列下标: 代谢物 id}，只保留能定位到代谢物的交换反应
        "exchanges": {c: col2met[c] for c in exch if c in col2met},
    }


# ============================================================================
# 三、读日志
# ============================================================================
def read_flux_log(path: Path, cols_by_model: dict[int, list[int]]) -> dict:
    """COMETS flux log: 每行 `cycle x y model_index flux1 flux2 ...`（model_index 从 1 开始）。

    只挑出各模型交换反应那几列，别的几千列不用解析。
    返回 {model_index: {cycle: {列下标: 通量}}}
    """
    out = {mi: {} for mi in cols_by_model}
    with path.open(errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 5:
                continue
            cycle, model = int(parts[0]), int(parts[3])
            cols = cols_by_model.get(model)
            if cols is None:
                continue
            fluxes = parts[4:]
            out[model][cycle] = {c: float(fluxes[c]) for c in cols if c < len(fluxes)}
    return out


def read_biomass_log(path: Path | None) -> dict[int, list[float]]:
    """totalbiomasslog: 每行 `cycle biomass_model1 biomass_model2 ...`（单位 g）。"""
    out = {}
    if not (path and path.is_file()):
        return out
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            try:
                out[int(parts[0])] = [float(v) for v in parts[1:]]
            except ValueError:
                pass
    return out


def read_media_log(path: Path | None) -> dict[str, tuple[float, float]]:
    """media log: 每行 `metabolite cycle x y mmol`。只取每个物质的首末值做自检。"""
    first, last = {}, {}
    if not (path and path.is_file()):
        return {}
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 5:
            met, cycle = parts[0], int(parts[1])
            try:
                val = float(parts[4])
            except ValueError:
                continue
            if met not in first or cycle < first[met][0]:
                first[met] = (cycle, val)
            if met not in last or cycle > last[met][0]:
                last[met] = (cycle, val)
    return {m: (first[m][1], last[m][1]) for m in first}


# ============================================================================
# 四、把通量积分成绝对交换量
# ============================================================================
def integrate(flux, models, biomass, dt_hours, tol=1e-9):
    """{代谢物: {模型下标: {'up': mmol, 'sec': mmol, 'peak_up':, 'peak_sec':}}}

    amount = Σ_t |flux| * biomass * Δt，通量 mmol/gDW/h × 生物量 g × 小时 = mmol。
    摄取和排出分开累加：同一个物质可能前期被排出、后期被吃回去，相减会把
    双向交换抹平，而那恰好是最该画出来的东西。
    """
    result = {}
    for mi, model in enumerate(models, start=1):
        series = flux.get(mi, {})
        for cycle, row in series.items():
            # 没有生物量日志就退化成"按通量本身加权"，等价于把生物量当作 1
            bm = biomass.get(cycle, [1.0] * len(models))
            b = bm[mi - 1] if mi - 1 < len(bm) else 1.0
            for col, v in row.items():
                if abs(v) < tol:
                    continue
                met = model["exchanges"][col]
                slot = result.setdefault(met, {}).setdefault(
                    mi, {"up": 0.0, "sec": 0.0, "peak_up": 0.0, "peak_sec": 0.0})
                if v < 0:                      # COMETS 约定：负 = 摄取
                    slot["up"] += -v * b * dt_hours
                    slot["peak_up"] = max(slot["peak_up"], -v)
                else:                          # 正 = 排出
                    slot["sec"] += v * b * dt_hours
                    slot["peak_sec"] = max(slot["peak_sec"], v)
    return result


def snapshot(flux, models, cycle, tol=1e-9):
    """{代谢物: {模型下标: {'up': 通量, 'sec': 通量}}}，某一个 cycle 的瞬时通量。

    形状和 integrate() 一模一样，所以 classify / 筛选 / 排序 / draw 全都不用改，
    换掉传下去的字典就行。区别只在值的含义：这里是 mmol/gDW/h，既不乘生物量也
    不乘时间，就是那一刻 FBA 解出来的交换通量本身。

    注意这份不能拿去和 media log 对账 —— 对账比的是整轮累计量，那得用 integrate()。
    """
    result = {}
    for mi, model in enumerate(models, start=1):
        for col, v in flux.get(mi, {}).get(cycle, {}).items():
            if abs(v) < tol:
                continue
            met = model["exchanges"][col]
            slot = result.setdefault(met, {}).setdefault(mi, {"up": 0.0, "sec": 0.0})
            if v < 0:                            # COMETS 约定：负 = 摄取
                slot["up"] = -v
            else:                                # 正 = 排出
                slot["sec"] = v
    return result


def classify(amounts: dict, n_models: int, floor: float) -> dict[str, str]:
    """给每个代谢物打标签。

    cross-fed    一个菌排出、另一个菌摄取 —— 真正的互喂
    shared       两个菌都在摄取 —— 底物竞争
    co-secreted  两个菌都在排出 —— 共同产物（CO2 这类），不是互喂
    exclusive    只有一个菌碰它 —— 独占底物 / 单方面产物

    floor 是「算不算真的动过」的下限，用来滤掉 1e-14 级别的求解器数值噪声。
    """
    labels = {}
    for met, per_model in amounts.items():
        up = {mi for mi, d in per_model.items() if d["up"] > floor}
        sec = {mi for mi, d in per_model.items() if d["sec"] > floor}
        if not up and not sec:
            continue
        if any(s != u for s in sec for u in up):
            labels[met] = "cross-fed"
        elif len(up) >= min(2, n_models):
            labels[met] = "shared"
        elif len(sec) >= min(2, n_models):
            labels[met] = "co-secreted"
        else:
            labels[met] = "exclusive"
    return labels


def check_against_media(amounts, media, factor=0.5, rel=0.2,
                        abs_tol=1e-8):
    """flux log 必须和 media log 对得上，否则图是假的。

    media log 记的是 COMETS 真正做了什么（池子里还剩多少），flux log 只是它
    写出来的一个通量向量。在 obj_style='MAX_OBJECTIVE_MIN_TOTAL' 下两者不相等：
    真正的碳源（glc/cit/sucr）在 flux log 里通量是 0，而一堆培养基里根本没有的
    物质（shikimate、quinate、valine）却有大通量——那是模型里的 infeasible loop。
    拿这种 flux log 画出来的 cross-feeding 全是伪迹，而且从图上看不出来。

    只查“确实被消耗掉了”的物质（变化超过初始量的 rel），static 的盐、O2、
    水只会微动，不会误报。
    """
    bad = []
    for met, (first, last) in media.items():
        delta = last - first
        if first <= 0 or abs(delta) < max(rel * first, abs_tol):
            continue
        net = sum(d["sec"] - d["up"] for d in amounts.get(met, {}).values())
        if abs(net) < factor * abs(delta):
            bad.append((met, delta, net))
    return bad


def pick_unit(max_mmol: float) -> tuple[float, str]:
    """自动选一个让数字落在 1~1000 的单位。"""
    for factor, name in [(1, "mmol"), (1e3, "µmol"), (1e6, "nmol"), (1e9, "pmol")]:
        if max_mmol * factor >= 1:
            return factor, name
    return 1e9, "pmol"


def pretty_met(met: str) -> str:
    return PRETTY.get(met, met[:-2] if met.endswith("_e") else met)


def pretty_species(name: str) -> str:
    """`B_megaterium` -> `B. megaterium`。"""
    parts = name.replace("-", "_").split("_")
    if len(parts) >= 2 and len(parts[0]) <= 2:
        return f"{parts[0]}. {' '.join(parts[1:])}"
    return name.replace("_", " ")


# ============================================================================
# 五、画图
# ============================================================================
def draw(rows, models, labels, amounts, unit, meta, out_paths, annotate, fig_w, row_h):
    """三栏图：左菌 | 中间竖排代谢物 | 右菌。

    rows 已经排好序，每项是 (代谢物 id, 类别)。
    箭头指向菌 = 摄取；箭头指向代谢物 = 排出；线宽 ∝ 交换总量（默认 sqrt 压缩，
    不然最大最小差三个数量级的时候细箭头会直接消失）。
    """
    factor, unit_name = unit
    n = len(rows)

    # --- 版面尺寸（先算英寸，再折算成 figure 比例，这样留白才可控）-------------
    pad_l, pad_r, pad_t, pad_b = 0.40, 0.40, 1.30, 1.15
    fig_h = pad_t + pad_b + n * row_h
    ax_w, ax_h = fig_w - pad_l - pad_r, fig_h - pad_t - pad_b
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([pad_l / fig_w, pad_b / fig_h, ax_w / fig_w, ax_h / fig_h])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    # 圆角在 x/y 数据单位上会被拉伸，用 mutation_aspect 抵消掉
    m_aspect = ax_w / ax_h

    # --- 横向坐标：菌盒 | 箭头区 | 代谢物 | 箭头区 | 菌盒 ---------------------
    box_l, box_r = (0.005, 0.145), (0.855, 0.995)
    pill_l, pill_r = 0.395, 0.605
    y_of = lambda i: 1.0 - (i + 0.5) / n       # 第 0 行在最上面
    half = 0.5 / n
    pill_h, offset = 0.62 * (1 / n), 0.19 * (1 / n)

    # --- 隔行浅底纹，帮眼睛对齐同一行的左右箭头 -----------------------------
    for i in range(n):
        if i % 2:
            ax.add_patch(plt.Rectangle((0.16, y_of(i) - half), 0.68, 2 * half,
                                       facecolor="#000000", alpha=0.028,
                                       edgecolor="none", zorder=0))

    # --- 两个菌：贯穿整列高度的圆角盒，名字竖排 -----------------------------
    for k, (model, (x0, x1)) in enumerate(zip(models, [box_l, box_r])):
        ax.add_patch(FancyBboxPatch(
            (x0, 0.004), x1 - x0, 0.992,
            boxstyle="round,pad=0,rounding_size=0.012", mutation_aspect=m_aspect,
            facecolor=C_ORG_FILL[k], edgecolor=C_ORG[k], linewidth=1.6, zorder=2))
        ax.text((x0 + x1) / 2 + 0.022, 0.5, pretty_species(model["name"]),
                rotation=90, ha="center", va="center", fontsize=15,
                fontstyle="italic", color=C_ORG[k], zorder=3)
        note = meta["biomass_note"][k]
        if note:
            ax.text((x0 + x1) / 2 - 0.030, 0.5, note, rotation=90, ha="center",
                    va="center", fontsize=8, color=C_MUTED, zorder=3)

    # --- 线宽映射：sqrt 压缩动态范围，最小值也保证画得出来 -------------------
    all_amounts = [v for met, _ in rows for d in amounts[met].values()
                   for v in (d["up"], d["sec"]) if v > 0]
    amax = max(all_amounts) if all_amounts else 1.0
    width = lambda a: LW_MIN + (LW_MAX - LW_MIN) * (min(a / amax, 1.0) ** 0.5)

    # --- 中间列：代谢物 pill --------------------------------------------------
    fills = {"cross-fed": (C_PILL_CROSS_FILL, C_PILL_CROSS_EDGE),
             "shared": (C_PILL_SHARE_FILL, C_PILL_SHARE_EDGE),
             "co-secreted": (C_PILL_COSEC_FILL, C_PILL_COSEC_EDGE),
             "exclusive": (C_PILL_EXCL_FILL, C_PILL_EXCL_EDGE)}
    for i, (met, cat) in enumerate(rows):
        y, (fc, ec) = y_of(i), fills[cat]
        ax.add_patch(FancyBboxPatch(
            (pill_l, y - pill_h / 2), pill_r - pill_l, pill_h,
            boxstyle="round,pad=0,rounding_size=0.010", mutation_aspect=m_aspect,
            facecolor=fc, edgecolor=ec, linewidth=1.2, zorder=4))
        ax.text(0.5, y + pill_h * 0.13, pretty_met(met), ha="center", va="center",
                fontsize=10.5, color=C_TEXT, zorder=5)
        ax.text(0.5, y - pill_h * 0.30, met, ha="center", va="center",
                fontsize=7, color=C_MUTED, family="monospace", zorder=5)

    # --- 箭头 ----------------------------------------------------------------
    def arrow(xa, xb, y, amount, color):
        lw = width(amount)
        ax.add_patch(FancyArrowPatch(
            (xa, y), (xb, y), arrowstyle="-|>",
            mutation_scale=9 + 1.7 * lw, linewidth=lw, color=color,
            shrinkA=1.5, shrinkB=1.5, joinstyle="miter", zorder=3))
        if annotate:
            ax.text((xa + xb) / 2, y + 0.30 * (1 / n), f"{amount * factor:.3g}",
                    ha="center", va="bottom", fontsize=7, color=color, zorder=6)

    for i, (met, _) in enumerate(rows):
        y = y_of(i)
        for mi in (1, 2):
            d = amounts[met].get(mi)
            if not d:
                continue
            x_box = box_l[1] if mi == 1 else box_r[0]
            x_pill = pill_l if mi == 1 else pill_r
            both = d["up"] > 0 and d["sec"] > 0     # 同一物质双向流动 -> 上下错开画
            if d["up"] > 0:                         # 代谢物 -> 菌
                arrow(x_pill, x_box, y + (offset if both else 0), d["up"], C_UPTAKE)
            if d["sec"] > 0:                        # 菌 -> 代谢物
                arrow(x_box, x_pill, y - (offset if both else 0), d["sec"], C_SECRETION)

    # --- 标题 -----------------------------------------------------------------
    fig.text(pad_l / fig_w, 1 - 0.42 / fig_h, meta["title"],
             fontsize=15.5, va="center", color=C_TEXT)
    fig.text(pad_l / fig_w, 1 - 0.76 / fig_h, meta["subtitle"],
             fontsize=9, va="center", color=C_MUTED)

    # --- 图例：方向配色 + pill 类别 + 线宽标尺 --------------------------------
    lax = fig.add_axes([pad_l / fig_w, 0.14 / fig_h, ax_w / fig_w, 0.80 / fig_h])
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    lax.legend(
        handles=[
            Line2D([], [], color=C_UPTAKE, lw=3, label="Uptake  (metabolite → cell)"),
            Line2D([], [], color=C_SECRETION, lw=3, label="Secretion  (cell → metabolite)"),
            Patch(facecolor=C_PILL_CROSS_FILL, edgecolor=C_PILL_CROSS_EDGE,
                  label="Cross-fed  (one secretes, other consumes)"),
            Patch(facecolor=C_PILL_SHARE_FILL, edgecolor=C_PILL_SHARE_EDGE,
                  label="Shared substrate  (both consume)"),
            Patch(facecolor=C_PILL_COSEC_FILL, edgecolor=C_PILL_COSEC_EDGE,
                  label="Co-secreted  (both produce)"),
        ],
        loc="upper left", bbox_to_anchor=(0, 1.10), frameon=False,
        fontsize=8.5, ncol=2, handlelength=2.2, columnspacing=1.8,
        labelspacing=0.45, borderaxespad=0)

    # 线宽标尺：三档参考值，读者才知道"粗"到底是多少
    x0 = 0.62
    lax.text(x0, 0.96, f"Arrow width = {meta['quantity']} ({unit_name})",
             fontsize=8.5, color=C_TEXT, va="top")
    for j, frac in enumerate([1.0, 0.4, 0.1]):
        a = amax * frac
        y = 0.62 - j * 0.26
        lax.add_line(Line2D([x0, x0 + 0.10], [y, y], color=C_MUTED,
                            lw=width(a), solid_capstyle="butt"))
        lax.text(x0 + 0.125, y, f"{a * factor:.3g}", fontsize=8,
                 color=C_MUTED, va="center")

    for p in out_paths:
        fig.savefig(p, dpi=600, bbox_inches=None,
                    facecolor="white", transparent=False)
    plt.close(fig)


# ============================================================================
# 六、主流程
# ============================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="从 COMETS flux log 画双菌 cross-feeding 图",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", type=Path,
                    help="这次模拟的任意一个 log（推荐 fluxlog_0x...），其余文件自动找")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="输出文件；不给扩展名就同时出 .pdf/.png/.svg（默认 crossfeeding_<runid>）")
    ap.add_argument("--top", type=int, default=14, help="最多画几个代谢物（按交换量排，默认 14）")
    ap.add_argument("--min-frac", type=float, default=0.004,
                    help="小于最大交换量这个比例的物质丢掉，滤数值噪声（默认 0.004）")
    ap.add_argument("--exclude", default="", help="额外排除的代谢物 id，逗号分隔")
    ap.add_argument("--only", default="", help="只画这些代谢物 id，逗号分隔（会覆盖其它筛选）")
    ap.add_argument("--keep-inorganic", action="store_true",
                    help="别自动扣掉水/质子/O2/CO2/无机盐")
    ap.add_argument("--include-exclusive", action="store_true",
                    help="连『只有一个菌用』的物质也画（默认只画互喂和共用）")
    ap.add_argument("--labels", default="", help="两个菌的显示名，逗号分隔，覆盖自动命名")
    ap.add_argument("--annotate", action="store_true", help="在箭头上标数值")
    ap.add_argument("--fig-width", type=float, default=11.0, help="图宽，英寸（默认 11）")
    ap.add_argument("--row-height", type=float, default=0.66, help="每行高，英寸（默认 0.66）")
    ap.add_argument("--at-hour", type=float, default=None,
                    help="画这个时刻的瞬时通量 mmol/gDW/h（会取最近的 cycle）；"
                         "不给就还是画整轮累计交换量")
    ap.add_argument("--csv", type=Path, default=None, help="把数值表也导出成 csv")
    ap.add_argument("--skip-check", action="store_true",
                    help="跳过『flux log vs media log』对账（确认过才用）")
    ap.add_argument("--force", action="store_true",
                    help="跳过『.cmd 比 flux log 新』的错位检查（确认过才用）")
    args = ap.parse_args(argv)

    # ---- 1. 找文件、读参数 --------------------------------------------------
    files = find_run_files(args.log)
    gparams, pparams = read_keyvals(files["params"]), read_keyvals(files["package"])
    time_step = float(pparams.get("timeStep", 0.05))
    flux_rate = int(float(gparams.get("FluxLogRate", 1)))
    dt = time_step * flux_rate                  # 每个 flux 采样点代表多少小时

    model_files = read_layout_models(files["layout"])
    if len(model_files) < 2:
        raise SystemExit(f"layout 里没读到两个模型: {files['layout']}")
    check_models_fresh(files, model_files[:2], args.force)
    models = [parse_cmd(files["dir"] / name) for name in model_files[:2]]
    if args.labels:
        for model, label in zip(models, args.labels.split(",")):
            model["name"] = label.strip()

    # ---- 2. 读 flux log，只取交换反应那些列 --------------------------------
    cols = {mi: sorted(m["exchanges"]) for mi, m in enumerate(models, start=1)}
    flux = read_flux_log(files["flux"], cols)
    if not any(flux.values()):
        raise SystemExit(f"flux log 是空的或者格式对不上: {files['flux']}")
    biomass = read_biomass_log(files["biomass"])
    if not biomass:
        print("! 没有 totalbiomasslog，箭头粗细退化为按通量加权（不是绝对物质量）",
              file=sys.stderr)

    cycles = sorted({c for s in flux.values() for c in s})
    sim_hours = max(cycles) * time_step

    # ---- 3. 积分 + 分类 ------------------------------------------------------
    # integrate() 永远要算：下面和 media log 对账靠它，那是唯一能识破坏 flux log
    # 的检查。至于图上画哪一份，另外决定。
    amounts = integrate(flux, models, biomass, dt)

    if args.at_hour is None:                     # 默认：整轮累计交换量
        plotted = amounts
        quantity = "total exchanged"                                  # 图例用，要短
        quantity_long = "total amount exchanged over the run"         # 副标题用
    else:                                        # 某一时刻的瞬时通量
        cycle = min(cycles, key=lambda c: abs(c * time_step - args.at_hour))
        got = cycle * time_step
        if abs(got - args.at_hour) > time_step:      # 要的时刻不在这次模拟里
            print(f"! --at-hour {args.at_hour} 超出了这次模拟的范围（0~{sim_hours:.1f} h），"
                  f"取最近的 t = {got:.2f} h", file=sys.stderr)
        plotted = snapshot(flux, models, cycle)
        quantity = quantity_long = f"flux at t = {got:.2f} h"
        if not plotted:
            raise SystemExit(
                f"cycle {cycle}（t = {cycle * time_step:.2f} h）的通量全是 0。\n"
                f"这次模拟跑了 {sim_hours:.1f} h，碳耗尽之后通量会衰减到 0，"
                "换一个更早的 --at-hour 试试。")

    biggest = max((v for d in plotted.values() for s in d.values()
                   for v in (s["up"], s["sec"])), default=0.0)
    # 分类只需要滤掉求解器数值噪声，门槛要低；显示门槛在筛完之后另算，
    # 否则 H2O/O2 那种大宗流量会把有机物的门槛顶到天上去，全被误判成 exclusive。
    labels = classify(plotted, len(models), floor=biggest * 1e-5)

    # ---- 3.5 与 media log 对账，对不上就不要画 ----------------------------
    media = read_media_log(files["media"])
    if media and not args.skip_check:
        bad = check_against_media(amounts, media)
        if bad:
            rows_bad = [f"    {m:<12} media 变化 {d * 1e6:>12.4g} nmol"
                        f"   flux 净交换 {n * 1e6:>12.4g} nmol" for m, d, n in bad]
            raise SystemExit("\n".join([
                "flux log 和 media log 对不上，下面这些物质在培养基里被消耗掉了，",
                "flux log 里却几乎没有对应的通量：",
                *rows_bad,
                "",
                "这说明 COMETS 写出的 flux log 不是它实际拿去更新培养基的那个解，",
                "常见原因是 obj_style='MAX_OBJECTIVE_MIN_TOTAL'（配上 GLOP）。改成",
                "    cm.obj_style = 'MAXIMIZE_OBJECTIVE_FLUX'",
                "重跑一次就好。想确认可以先跑 verify_fluxlog.py。",
                "确实知道自己在干什么可以用 --skip-check 跳过。"]))

    # ---- 4. 筛选 -------------------------------------------------------------
    if args.only:
        keep = {m.strip() for m in args.only.split(",") if m.strip()}
        rows_all = [(m, labels.get(m, "exclusive")) for m in plotted if m in keep]
    else:
        drop = set() if args.keep_inorganic else set(INORGANIC)
        drop |= {m.strip() for m in args.exclude.split(",") if m.strip()}
        wanted = {"cross-fed", "shared", "co-secreted"}
        if args.include_exclusive:
            wanted.add("exclusive")
        rows_all = [(m, c) for m, c in labels.items() if c in wanted and m not in drop]

    total_of = lambda m: sum(s["up"] + s["sec"] for s in plotted[m].values())
    # 显示门槛按「留下来的这批里最大的那个」算，而不是按全局最大
    shown_max = max((total_of(m) for m, _ in rows_all), default=0.0)
    rows_all = [r for r in rows_all if total_of(r[0]) > shown_max * args.min_frac]
    # 排序：互喂的排前面，同类内按交换量从大到小
    order = {"cross-fed": 0, "shared": 1, "co-secreted": 2, "exclusive": 3}
    rows_all.sort(key=lambda r: (order[r[1]], -total_of(r[0])))
    rows = rows_all[:args.top]
    if not rows:
        raise SystemExit("筛完没有任何代谢物。试试 --keep-inorganic 或 --include-exclusive。")

    # ---- 5. 打印数值表（图上看不清的数字在这里）-----------------------------
    # pick_unit() 是按「总量 mmol」设计的；瞬时通量量纲不同（mmol/gDW/h，本来就在
    # 0.1~20 量级），硬套会把它标成 nmol。
    unit = ((1.0, "mmol/gDW/h") if args.at_hour is not None
            else pick_unit(max(total_of(m) for m, _ in rows)))
    factor, unit_name = unit
    names = [pretty_species(m["name"]) for m in models]
    print(f"\nrun {files['run_id']}   {sim_hours:.1f} h   "
          f"timeStep {time_step} h   FluxLogRate {flux_rate}   "
          f"{len(cycles)} flux samples")
    print(f"models: {names[0]} (#1), {names[1]} (#2)\n")
    head = f"{'metabolite':<14}{'class':<11}"
    for nm in names:
        head += f"{nm + ' up':>16}{nm + ' sec':>17}"
    print(head + f"    [{unit_name}]")
    print("-" * len(head))
    for met, cat in rows_all:
        line = f"{met:<14}{cat:<11}"
        for mi in (1, 2):
            d = plotted[met].get(mi, {"up": 0.0, "sec": 0.0})
            line += f"{d['up'] * factor:>16.4g}{d['sec'] * factor:>17.4g}"
        print(line + ("" if (met, cat) in rows else "   (未画出)"))

    # 质量守恒自检：两菌净交换应当等于 media log 里池子的净变化。
    # 这里永远用 amounts（整轮累计）和固定的 nmol —— 比的是总量，不该跟着图的单位走。
    if media:
        print("\n自检  net exchange vs media log 净变化 "
              "（差得大通常是该物质被设成了 static，或 FluxLogRate>1 积分粗）  [nmol]")
        for met, _ in rows:
            net = sum(d["sec"] - d["up"] for d in amounts.get(met, {}).values())
            if met in media:
                obs = media[met][1] - media[met][0]
                print(f"  {met:<14}flux {net * 1e6:>12.4g}   "
                      f"media {obs * 1e6:>12.4g}")

    if args.csv:
        col_unit = "mmol_gDW_h" if args.at_hour is not None else "mmol"
        with args.csv.open("w", encoding="utf-8") as fh:
            fh.write("metabolite,class," + ",".join(
                f"{n}_uptake_{col_unit},{n}_secretion_{col_unit}" for n in names) + "\n")
            for met, cat in rows_all:
                vals = []
                for mi in (1, 2):
                    d = plotted[met].get(mi, {"up": 0.0, "sec": 0.0})
                    vals += [f"{d['up']:.6g}", f"{d['sec']:.6g}"]
                fh.write(f"{met},{cat}," + ",".join(vals) + "\n")
        print(f"\n数值表 -> {args.csv}")

    # ---- 6. 画 ---------------------------------------------------------------
    stem = args.out or (files["dir"] / f"crossfeeding_{files['run_id']}")
    out_paths = ([stem] if stem.suffix.lower() in {".pdf", ".png", ".svg"}
                 else [stem.with_suffix(e) for e in (".pdf", ".png", ".svg")])

    tally = ", ".join(f"{sum(1 for _, c in rows if c == cat)} {cat}"
                      for cat in ("cross-fed", "shared", "co-secreted", "exclusive")
                      if any(c == cat for _, c in rows))
    meta = {
        "title": f"{names[0]}  ↔  {names[1]}: metabolite exchange",
        "subtitle": (f"COMETS co-culture, {sim_hours:.0f} h · arrow width = "
                     f"{quantity_long} ({unit_name}) · "
                     f"{tally} · run {files['run_id']}"),
        "quantity": quantity,
        "biomass_note": [
            (f"final {biomass[max(biomass)][i]:.2e} g" if biomass else "")
            for i in range(len(models))],
    }
    draw(rows, models, labels, plotted, unit, meta, out_paths,
         args.annotate, args.fig_width, args.row_height)
    print("\n图 -> " + "  ".join(str(p) for p in out_paths))


if __name__ == "__main__":
    main()
