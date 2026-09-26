# Google Colab 运行入口

2026-09-26 启动修复：云端入口已更新为 `host_pip_v2`，使用 `venv --without-pip` 创建隔离环境，再由宿主 `pip --python` 安装依赖，绕过自动调用 ensurepip 的环节。环境创建、宿主 pip 和依赖安装分别记录到 `venv.log`、`host_pip.log`、`install.log`。原报错只确认 venv 创建失败，未保留底层 stderr，不能据此断言唯一原因。启动修复已通过本地隔离安装验证，云端 notebook 回读哈希一致；仍待用户重新运行验收。项目 ZIP、科学计算源码和参数保持原快照不变；旧标签页不要覆盖云端新版。

2026-09-25：以最新 `4daaf68` 进度为基准。默认 `adaptive_log_a`、严格联合初态、Prada12-HMF 和现有 K 设置；方案 B 的 gw_aged cohort 已纳入云端试算。旧方案文档中的 `e0ff270` 和“尚未实现 B”已过时。

运行 `python cloud/build_bundle.py` 生成 `cloud/build/` 中的项目 ZIP、逐文件 SHA-256 manifest、上传清单和 `colab_bootstrap.ipynb`。快照标识同时包含基准提交与实际文件清单哈希；云端续算扩展会明确记录为基准之上的改动。保留全部已跟踪研究文件、notebook、缓存及作者原图，排除 Git/编辑器缓存与其他临时目录。

将 ZIP 放入 `My Drive/PBH_Merger_Colab/snapshots/`，将入口 notebook 放入 `PBH_Merger_Colab/`。在 Colab 打开入口，选择 CPU 并全部运行。Drive 授权由用户完成；Google Drive 工具可上传文件，但不能直接启动 Colab 计算。

入口在 `/content` 解压并校验源码，在独立 venv 安装锁定版本依赖，防止 Colab 已加载的包污染计算。Python 需要 >=3.12。先比较 128 个已保存初态的 R1/R5、z=12/0、full/GW-only 共 8 组；逐元素事件分类要求一致，浮点量容差写入 `comparison.json`。比较不通过时停止扩大计算。

随后试跑 M0=1.5e6 的 B/full 与 B/GW-only：初始每壳 128、每入晕边界 8，seed=20260923。该规模只用于验证执行链，不是统计收敛研究。改变 notebook 参数可运行多种子、样本量、环境步长及 `gw_aged/pristine` 对照。每组输出含配置、源码哈希、人口账本、率 CSV/PNG、原始 NPZ 和最终存活轨道。

`run_cohort_shell_monte_carlo` 新增可选 `checkpoint_path`、`resume`：保存每个完整边界的账本、存活轨道/权重/cohort 身份和 RNG 状态。恢复时校验配置、源码、采样器和时间网格；拒绝混用旧检查点。已完成组直接读取 summary。此扩展不改变默认物理模型或无检查点调用的随机序列。目前仅 B 支持中途续算，A 保留原入口。

云端输出位于 `PBH_Merger_Colab/experiments/<snapshot>/`。`runtime_status.json` 区分 mounted、comparison passed、running、complete、failed/interrupted；`complete` 仅由真正完成的云端代码写出。安装/运行日志写回同目录。长运行中断后重新运行入口即可继续；新的代码/配置生成独立结果目录。

本地可复跑：

```powershell
python cloud/run_cloud.py reference --output cloud/reference
python cloud/run_cloud.py compare --output cloud/local_runs/comparison
python cloud/run_cloud.py cohort --output cloud/local_runs/pilot --mode full
python cloud/run_cloud.py cohort --output cloud/local_runs/pilot --mode gw_only
python cloud/build_bundle.py
```

本地运行与云端运行分别记录 platform；本地成功不代表 Colab 已实际执行。云端结果通过 Drive 回读后才验收云端计算。原 notebook 和旧缓存保留，云端入口产生独立结果。
