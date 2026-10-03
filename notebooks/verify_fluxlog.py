#!/usr/bin/env python3
"""跑一次短模拟，检查 flux log 是不是真的等于 COMETS 实际用的那个解.

为什么需要这个脚本
------------------
crossfeeding_figure.py 完全依赖 flux log 来判断"谁吃谁排"。但 flux log 只是
COMETS 写出来的一个通量向量，它不保证等于 COMETS 拿去更新培养基的那个解。

对不上的时候症状很隐蔽：图照画，箭头照给，只是画的全是模型里的
thermodynamically infeasible loop（凭空吃进 shikimate/quinate/valine，
真正的碳源 glc/cit/sucr 通量却是 0），看起来像 cross-feeding，其实是伪迹。

判据（唯一可靠的那个）
----------------------
media log 记的是每个物质在格子里还剩多少，这是 COMETS 真正做的事。
flux log 积分出来的净交换量必须等于 media log 的净变化：

    Sum_t flux(t) * biomass(t) * dt   ==   media(末) - media(初)

对每个非 static 的供给物质都成立，才说明 flux log 可用。

用法
----
    python verify_fluxlog.py                          # 默认两种 obj_style 都跑
    python verify_fluxlog.py --obj-style MAXIMIZE_OBJECTIVE_FLUX --hours 6
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
import tempfile
from pathlib import Path

import cobra
import cometspy as c

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# 配置：与 Coculture_Bmeg_Paer.ipynb 逐字一致，改了就不是同一个实验了
# ---------------------------------------------------------------------------
EXUDATE_UM = {'cit_e': 8.284, 'glu__L_e': 5.581, 'asp__L_e': 2.673, 'ser__L_e': 1.321,
              'glc__D_e': 1.065, 'sucr_e': 3.077, 'fru_e': 0.710, 'sbt__D_e': 0.355}
SALTS_MM = {'nh4_e': 4.5, 'k_e': 6.3, 'cl_e': 16.03, 'ca2_e': 5.0,
            'mg2_e': 2.0, 'so4_e': 2.087, 'pi_e': 0.3, 'mn2_e': 0.01438,
            'zn2_e': 0.00076, 'mobd_e': 0.00033, 'na1_e': 0.00066, 'fe2_e': 0.08669,
            'o2_e': 0.2, 'cu2_e': 0.001, 'cobalt2_e': 0.0001}
BULK = {'h2o_e': 1000.0, 'h_e': 1000.0}

DEFAULT_VMAX, DEFAULT_KM = 25.0, 3e-6
TIME_STEP, INITIAL_BIOMASS = 0.05, 1e-8
KM_PREFERRED, KM_REPRESSED = 0.5e-6, 20e-6

KINETICS_BM_RAW = {'glc__D_e': (6.0, KM_PREFERRED), 'fru_e': (5.0, KM_PREFERRED),
                   'sucr_e': (5.0, KM_PREFERRED), 'sbt__D_e': (4.0, KM_PREFERRED),
                   'cit_e': (1.2, KM_REPRESSED), 'glu__L_e': (1.2, KM_REPRESSED),
                   'asp__L_e': (1.2, KM_REPRESSED), 'ser__L_e': (1.2, KM_REPRESSED)}
KINETICS_PA_RAW = {'cit_e': (7.0, KM_PREFERRED), 'glu__L_e': (5.0, KM_PREFERRED),
                   'asp__L_e': (4.0, KM_PREFERRED), 'ser__L_e': (4.0, KM_PREFERRED),
                   'glc__D_e': (0.8, KM_REPRESSED), 'fru_e': (0.8, KM_REPRESSED)}
BLOCKED_PA = {'sbt__D_e'}
CARBON_CAP_BM, CARBON_CAP_PA, VMAX_O2_PA = 30.0, 60.0, 20.0
MODEL_FILE_BM = 'B_megaterium_DSM319_iJA1121_bigg.xml'
MODEL_FILE_PA = 'P_aeruginosa_iSD1509.xml'


# ---------------------------------------------------------------------------
# 建模（notebook 里的 helper 原样搬过来）
# ---------------------------------------------------------------------------
def exchange_map(model):
    return {list(r.metabolites)[0].id: r.id for r in model.exchanges}


def scale_to_carbon_cap(kinetics, cap, carbon_of):
    draw = {m: v * (EXUDATE_UM[m] * 1e-6) / (EXUDATE_UM[m] * 1e-6 + k) * carbon_of(m)
            for m, (v, k) in kinetics.items()}
    f = cap / sum(draw.values())
    return {m: (v * f, k) for m, (v, k) in kinetics.items()}


def make_comets_model(model, model_id, exchanges, kinetics, obj_style, blocked=()):
    cm = c.model(model)
    cm.id = model_id
    for rxn in model.boundary:                       # 胞内 sink 不当作交换反应
        mets = list(rxn.metabolites)
        if len(mets) == 1 and mets[0].compartment != 'e':
            cm.reactions.loc[cm.reactions.REACTION_NAMES == rxn.id, 'EXCH'] = False
            cm.reactions.loc[cm.reactions.REACTION_NAMES == rxn.id, 'EXCH_IND'] = 0
    for rxn in model.exchanges:
        cm.change_bounds(rxn.id, -1000, 1000)
    for met in blocked:
        if met in exchanges:
            cm.change_bounds(exchanges[met], 0, 1000)
    for met, (vmax, km) in kinetics.items():
        cm.change_vmax(exchanges[met], vmax)
        if hasattr(cm, 'change_km'):
            cm.change_km(exchanges[met], km)
        else:
            cm.reactions.loc[cm.reactions.REACTION_NAMES == exchanges[met], 'KM'] = km
    cm.initial_pop = [0, 0, INITIAL_BIOMASS]
    cm.optimizer = 'GLOP'
    cm.obj_style = obj_style
    return cm


def build(obj_style):
    """返回 (layout, 两个 cometspy model, static 集合)。"""
    bm = cobra.io.read_sbml_model(str(HERE / MODEL_FILE_BM))
    pa = cobra.io.read_sbml_model(str(HERE / MODEL_FILE_PA))
    pa.reactions.get_by_id('ATPM').bounds = (0.0, 8.39)
    ex_bm, ex_pa = exchange_map(bm), exchange_map(pa)

    def carbon(met):
        for mdl in (bm, pa):
            if mdl.metabolites.has_id(met):
                return (mdl.metabolites.get_by_id(met).elements or {}).get('C')

    medium = {m: v * 1e-6 for m, v in EXUDATE_UM.items()}
    medium.update({m: v * 1e-3 for m, v in SALTS_MM.items()})
    medium.update(BULK)
    static = set(medium) - set(EXUDATE_UM)          # 只有 8 个碳源会被消耗

    k_bm = scale_to_carbon_cap(KINETICS_BM_RAW, CARBON_CAP_BM, carbon)
    k_pa = scale_to_carbon_cap(KINETICS_PA_RAW, CARBON_CAP_PA, carbon)
    cm_bm = make_comets_model(bm, 'B_megaterium', ex_bm, k_bm, obj_style)
    cm_pa = make_comets_model(pa, 'P_aeruginosa', ex_pa,
                              {**k_pa, 'o2_e': (VMAX_O2_PA, DEFAULT_KM)},
                              obj_style, blocked=BLOCKED_PA)

    layout = c.layout()
    layout.grid = [1, 1]
    for cm in (cm_bm, cm_pa):
        layout.add_model(cm)
    for met, amount in medium.items():
        layout.set_specific_metabolite(met, amount, static=(met in static))
    return layout, [cm_bm, cm_pa], static


# ---------------------------------------------------------------------------
# 自检：flux log 积分 vs media log 净变化
# ---------------------------------------------------------------------------
def load_cf():
    """把 crossfeeding_figure.py 当模块用，复用它的 log 解析。"""
    spec = importlib.util.spec_from_file_location('cf', HERE / 'crossfeeding_figure.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def compare(run_dir: Path, static: set, obj_style: str) -> bool:
    """True 表示 flux log 可信。"""
    cf = load_cf()
    flux_file = next(run_dir.glob('fluxlog_0x*'))
    files = cf.find_run_files(flux_file)
    models = [cf.parse_cmd(files['dir'] / n)
              for n in cf.read_layout_models(files['layout'])[:2]]
    cols = {mi: sorted(m['exchanges']) for mi, m in enumerate(models, 1)}
    flux = cf.read_flux_log(files['flux'], cols)
    biomass = cf.read_biomass_log(files['biomass'])
    amounts = cf.integrate(flux, models, biomass, TIME_STEP)
    media = cf.read_media_log(files['media'])

    print(f"\n=== obj_style = {obj_style} ===")
    print(f"{'metabolite':12}{'media delta':>16}{'flux net':>16}   verdict   [nmol]")
    ok = True
    for met in sorted(EXUDATE_UM):                   # 只看会被消耗的碳源
        if met in static or met not in media:
            continue
        obs = (media[met][1] - media[met][0]) * 1e6
        net = sum(d['sec'] - d['up'] for d in amounts.get(met, {}).values()) * 1e6
        good = abs(net - obs) <= 0.10 * abs(obs) + 1e-6   # 10% 以内算过
        ok &= good
        print(f"{met:12}{obs:16.4g}{net:16.4g}   {'ok' if good else 'MISMATCH'}")
    print("flux log " + ("可信，可以拿去画 cross-feeding 图"
                         if ok else "不可信 -- 它不是 COMETS 实际用的那个解"))
    return ok


def run_one(obj_style: str, hours: float, keep: bool) -> bool:
    """在独立目录里跑，避免覆盖 notebook 那次运行的 .cmd / log。"""
    run_dir = Path(tempfile.mkdtemp(prefix='verify_', dir=str(HERE)))
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        layout, models, static = build(obj_style)
        params = c.params()
        for k, v in {'numRunThreads': 1, 'defaultVmax': DEFAULT_VMAX,
                     'defaultKm': DEFAULT_KM, 'maxCycles': int(hours / TIME_STEP),
                     'timeStep': TIME_STEP, 'spaceWidth': 1, 'maxSpaceBiomass': 10,
                     'minSpaceBiomass': 1e-11, 'deathRate': 0.0,
                     'writeMediaLog': True, 'MediaLogRate': 1,
                     'writeTotalBiomassLog': True,
                     'writeFluxLog': True, 'FluxLogRate': 1}.items():
            params.set_param(k, v)
        c.comets(layout, params).run(False)
        os.chdir(cwd)
        return compare(run_dir, static, obj_style)
    finally:
        os.chdir(cwd)
        if not keep:
            shutil.rmtree(run_dir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--obj-style', action='append', default=None,
                    help='可重复；默认两种都试')
    ap.add_argument('--hours', type=float, default=6.0)
    ap.add_argument('--keep', action='store_true', help='保留临时运行目录')
    args = ap.parse_args()
    styles = args.obj_style or ['MAX_OBJECTIVE_MIN_TOTAL', 'MAXIMIZE_OBJECTIVE_FLUX']
    results = {s: run_one(s, args.hours, args.keep) for s in styles}
    print('\n汇总: ' + '  '.join(f'{s}={"ok" if v else "BAD"}' for s, v in results.items()))
    return 0 if all(results.values()) else 1


if __name__ == '__main__':
    sys.exit(main())
