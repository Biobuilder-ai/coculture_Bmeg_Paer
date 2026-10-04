# coculture_Bmeg_Paer

*Bacillus megaterium* (iJA1121) + *Pseudomonas aeruginosa* (iSD1509) 的 COMETS 共培养模拟 ——
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
