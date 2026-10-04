# coculture_Bmeg_Paer

**English** | [中文](#中文)

COMETS co-culture simulation of *Bacillus megaterium* (iJA1121) + *Pseudomonas aeruginosa* (iSD1509):
both strains share a single 1×1 well containing strawberry root exudate (SUGAR-SHIFTED profile,
14 carbon compounds totalling 142 µmol C/L). Spatial dynamic flux balance analysis (COMETS / cometspy)
is used to observe how they compete and cross-feed.

Main notebook: [`notebooks/Coculture_Bmeg_Paer_trans.ipynb`](notebooks/Coculture_Bmeg_Paer_trans.ipynb)

## Repository contents

| File | Description |
| --- | --- |
| `notebooks/Coculture_Bmeg_Paer_trans.ipynb` | Main notebook: modelling → pre-flight FBA → COMETS simulation → figures → run archiving |
| `notebooks/B_megaterium_DSM319_iJA1121_bigg.xml` | *B. megaterium* DSM319 genome-scale model (iJA1121, BiGG IDs) |
| `notebooks/P_aeruginosa_iSD1509.xml` | *P. aeruginosa* genome-scale model (iSD1509) |
| `notebooks/crossfeeding_figure.py` | Draws cross-feeding snapshot figures from flux/media/biomass logs (called automatically by the notebook) |
| `notebooks/verify_fluxlog.py` | Checks mass conservation between flux log and media log (optional standalone diagnostic) |
| `windows_run_comets_lab.bat` | Launch script for Windows |
| `linux_run_comets_lab.sh` | Launch script for Linux |
| `macos_run_comets_lab.sh` | Launch script for macOS |

## Runtime environment

The simulation runs inside the COMETS Docker environment (image name `comets-lab`), which contains:

- Python 3.10
- cometspy 0.6.3, cobra 0.30.0
- jupyterlab 4.5.4, matplotlib, numpy, pandas, scipy
- COMETS 2.12.3 (Java), with `COMETS_HOME=/opt/comets_linux` and `JAVA_TOOL_OPTIONS=-Dcomets.home=/opt/comets_linux`

> **Getting the Docker image**: the image file (`comets-lab.tar`, ~630 MB) exceeds GitHub's per-file
> limit, so it is distributed as a Release asset. Download it and place it in the **repository root**
> (next to the launch scripts):
>
> - Download URL: <https://github.com/Biobuilder-ai/coculture_Bmeg_Paer/releases/download/v1.0.0/comets-lab.tar>
> - Command line:
>   ```bash
>   curl -L -o comets-lab.tar https://github.com/Biobuilder-ai/coculture_Bmeg_Paer/releases/download/v1.0.0/comets-lab.tar
>   ```
> - Optional checksum: sha256 should be `d3e0e8f5bf6c2357da6eb8d6e9b02945ed91f82f870993735041e70a69161892`

## Quick start

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows/macOS) or Docker Engine (Linux). ≥ 8 GB RAM is recommended.
2. Put `comets-lab.tar` in the repository root.
3. Run the launch script for your platform:

   **Windows** (double-click or run in cmd):
   ```bat
   windows_run_comets_lab.bat
   ```

   **Linux**:
   ```bash
   chmod +x linux_run_comets_lab.sh
   ./linux_run_comets_lab.sh
   ```

   **macOS**:
   ```bash
   chmod +x macos_run_comets_lab.sh
   ./macos_run_comets_lab.sh
   ```

4. JupyterLab opens automatically in your browser (<http://localhost:8888>). The macOS script starts the
   container as a Jupyter server for VS Code to connect via "Existing Jupyter Server" (the script prints
   the URL at the end).
5. Open `notebooks/Coculture_Bmeg_Paer_trans.ipynb` and run all cells top to bottom.

### What the notebook does

1. **Helpers / Configuration**: loads both models, defines the strawberry root exudate medium and kinetic parameters.
2. **Model and layout**: builds the two-strain COMETS layout.
3. **Pre-flight check**: runs FBA on the medium COMETS will actually present, to confirm both strains grow before spending time on the simulation.
4. **Run**: runs the COMETS simulation, producing `fluxlog_*` / `medialog_*` / `totalbiomasslog_*` logs.
5. **Run archive + cross-feeding snapshots**: collects this run's logs and `.cmd` model files into a timestamped
   `comets_run_YYYYMMDD_HHMMSS/` directory, then calls `crossfeeding_figure.py` to emit instantaneous flux
   snapshot PDFs at given times (default 1 h, 5 h, 10 h). **Every rerun overwrites the `.cmd` files of the
   same name, so run this cell immediately after the simulation to archive the results.**

### Optional: verify the flux log

`crossfeeding_figure.py` relies entirely on the flux log to judge "who eats what, who excretes what".
If you suspect the arrows in the figure are artifacts of thermodynamically infeasible loops, run:

```bash
cd notebooks
python verify_fluxlog.py                          # both obj_styles by default
python verify_fluxlog.py --obj-style MAXIMIZE_OBJECTIVE_FLUX --hours 6
```

Criterion: the net exchange integrated from the flux log must equal the net change in the media log
(must hold for every non-static supplied metabolite).

## Stop the environment

```bash
docker stop comets-lab-container
```

## References / model sources

- Model iJA1121: *B. megaterium* DSM319 (Najar et al.)
- Model iSD1509: *P. aeruginosa* (Sánchez-García et al.)

---

# 中文

[English](#coculture_bmeg_paer) | **中文**

*Bacillus megaterium*（iJA1121）+ *Pseudomonas aeruginosa*（iSD1509）的 COMETS 共培养模拟 ——
两个菌株共享同一个 1×1 well 里的草莓根系分泌物（SUGAR-SHIFTED 配方，14 种碳化合物、共 142 µmol C/L），
用空间动态通量平衡（COMETS / cometspy）观察它们如何竞争与交叉喂养（cross-feeding）。

主 notebook：[`notebooks/Coculture_Bmeg_Paer_trans.ipynb`](notebooks/Coculture_Bmeg_Paer_trans.ipynb)

## 仓库内容

| 文件 | 说明 |
| --- | --- |
| `notebooks/Coculture_Bmeg_Paer_trans.ipynb` | 主 notebook：建模 → 预检 FBA → COMETS 模拟 → 出图 → run 归档 |
| `notebooks/B_megaterium_DSM319_iJA1121_bigg.xml` | *B. megaterium* DSM319 基因组尺度模型（iJA1121，BiGG 命名） |
| `notebooks/P_aeruginosa_iSD1509.xml` | *P. aeruginosa* 基因组尺度模型（iSD1509） |
| `notebooks/crossfeeding_figure.py` | 从 flux/media/biomass log 绘制 cross-feeding 快照图（notebook 会自动调用） |
| `notebooks/verify_fluxlog.py` | 校验 flux log 与 media log 物质守恒一致（可选的独立诊断工具） |
| `windows_run_comets_lab.bat` | Windows 启动脚本 |
| `linux_run_comets_lab.sh` | Linux 启动脚本 |
| `macos_run_comets_lab.sh` | macOS 启动脚本 |

## 运行环境

模拟依赖 COMETS 的 Docker 运行环境（镜像名 `comets-lab`），内含：

- Python 3.10
- cometspy 0.6.3、cobra 0.30.0
- jupyterlab 4.5.4、matplotlib、numpy、pandas、scipy
- COMETS 2.12.3（Java），环境变量 `COMETS_HOME=/opt/comets_linux`、`JAVA_TOOL_OPTIONS=-Dcomets.home=/opt/comets_linux`

> **获取 Docker 镜像**：镜像文件（`comets-lab.tar`，约 630 MB）体积超出 GitHub 仓库单文件限制，
> 已作为 Release 资产分发。下载后放到**仓库根目录**（与启动脚本同级）：
>
> - 下载地址：<https://github.com/Biobuilder-ai/coculture_Bmeg_Paer/releases/download/v1.0.0/comets-lab.tar>
> - 命令行下载：
>   ```bash
>   curl -L -o comets-lab.tar https://github.com/Biobuilder-ai/coculture_Bmeg_Paer/releases/download/v1.0.0/comets-lab.tar
>   ```
> - 校验（可选）：sha256 应为 `d3e0e8f5bf6c2357da6eb8d6e9b02945ed91f82f870993735041e70a69161892`

## 快速开始

1. 安装 [Docker Desktop](https://www.docker.com/products/docker-desktop/)（Windows/macOS）或 Docker Engine（Linux），建议 ≥ 8 GB 内存。
2. 把 `comets-lab.tar` 放到仓库根目录。
3. 按你的平台运行启动脚本：

   **Windows**（双击或在 cmd 中执行）：
   ```bat
   windows_run_comets_lab.bat
   ```

   **Linux**：
   ```bash
   chmod +x linux_run_comets_lab.sh
   ./linux_run_comets_lab.sh
   ```

   **macOS**：
   ```bash
   chmod +x macos_run_comets_lab.sh
   ./macos_run_comets_lab.sh
   ```

4. 浏览器会自动打开 JupyterLab（<http://localhost:8888>）。macOS 脚本会把容器作为
   Jupyter server 启动，供 VS Code 以 "Existing Jupyter Server" 方式连接（脚本末尾会打印 URL）。
5. 进入 `notebooks/`，打开 `Coculture_Bmeg_Paer_trans.ipynb`，从上到下执行全部 cell。

### notebook 里会发生什么

1. **Helpers / Configuration**：加载两个模型、定义草莓根系分泌物培养基与动力学参数。
2. **Model and layout**：构建 COMETS 双菌株 layout。
3. **Pre-flight check**：先在 COMETS 实际呈现的培养基上做 FBA，确认两株都能生长，再花时间跑模拟。
4. **Run**：执行 COMETS 模拟，生成 `fluxlog_*` / `medialog_*` / `totalbiomasslog_*` 等日志。
5. **Run archive + cross-feeding snapshots**：把本次运行的日志与 `.cmd` 模型文件收进带时间戳的
   `comets_run_YYYYMMDD_HHMMSS/` 目录，并调用 `crossfeeding_figure.py` 在指定时刻（默认 1 h、5 h、10 h）
   输出瞬时通量快照 PDF。**每次重新运行会重写同名 `.cmd` 文件，跑完请立刻执行这个 cell 归档。**

### 可选：校验 flux log

`crossfeeding_figure.py` 完全依赖 flux log 判断"谁吃谁排"。若怀疑图里的箭头是热力学不可行循环造成的伪迹，可运行：

```bash
cd notebooks
python verify_fluxlog.py                          # 默认两种 obj_style 都跑
python verify_fluxlog.py --obj-style MAXIMIZE_OBJECTIVE_FLUX --hours 6
```

判据：flux log 积分出的净交换量必须等于 media log 的净变化量（对每个非 static 供给物质都成立）。

## 停止环境

```bash
docker stop comets-lab-container
```

## 相关文献 / 模型出处

- 模型 iJA1121：*B. megaterium* DSM319（Najar et al.）
- 模型 iSD1509：*P. aeruginosa*（Sánchez-García et al.）
