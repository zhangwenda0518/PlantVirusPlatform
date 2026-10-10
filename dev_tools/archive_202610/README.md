# dev_tools/archive_202610 —— 历史一次性脚本（勿直接运行）

2026-10-03 归档：以下脚本均已在当时**执行完毕**或属于一次性的数据校验/
修补代码，保留仅为可追溯（`git log -- dev_tools/archive_202610/<file>`）。
**不要直接运行**：它们大多直接改源码/改数据、依赖当时本机存在的数据目录
（如 `databases/_dl_staging/`、`run/fastq/GQMIX.*.fq.gz`），现在跑会失败
或被误改。

| 文件 | 性质（读文件头的结论） |
|---|---|
| `patch_ckpos.py` | 一次性 UI 补丁（工具②报告区复选框归位），已落入前端 |
| `patch_k2.py` | 文件头即写明「【已执行完毕，存档备查】」（Kraken2 端点/按钮） |
| `fix_246tax.py` | 一次性补丁：给 build_virus_dbs.py 打 taxonomy 准备逻辑（脚本内硬编码 `/home/USER/...` 服务器路径） |
| `fix_skipwrite.py` | 一次性补丁：按行号删除 build_virus_dbs.py 里的 taxid\|0 写入 |
| `_prune_i18n_20260927.py` | 文件头写明「一次性脚本」：按当时 run/_tmp_i18n_removable.json 删 i18n 孤儿键 |
| `_reorg_split.py` | 一次性目录整理（程序/数据库/示例分离布局），执行时用改名搬运，含本机 D: 绝对路径 |
| `_speedtest_kv.py` | 一次性基准测量（salmon_k31 vs minibwa），绑定当时 GQMIX 数据集与已建好的 kv_index |
| `_check_dl_ref.py` | 一次性校验当时下载的 Plant_Virus_Ref 数据集（databases/_dl_staging/ 已不在） |
| `_check_kv_discover.py` | 一次性验证鉴定库自动发现（起 app.test_client 手工核对） |
| `_check_plant_ref.py` | 一次性验证 plant_ref 库接入（同上） |
| `_fix_plant_virus_ref.py` | 一次性数据修补：Plant_Virus_Ref 表-序列对齐 + TSWV 注入 |
| `fill_species_gaps.py` | 2026-10-03 从 `docs/paper/benchmark/` 归档：会直接追加写平台参考库、依赖外部全量库路径（C-host_classify）与已消失的 kv_index；已改为默认 dry-run（写入需显式 `--write`） |

未归档（仍在使用，见 `dev_tools/`）：

- `pyz_probe.py` / `pyz_vs_source.py` —— 打包产物 PYZ 检查工具；后者被
  `tests/_selftest_pyz_cmp.py` 以 `from dev_tools.pyz_vs_source import _cmp` 引用。
- `rna_early_import_hook.py` —— 被 `dev_tools/VirusPlatform.spec:167` 作为
  PyInstaller runtime hook 引用，是打包链的活代码。
- `package.py`、`build_virus_db*.py`、`dist_*.py`、`pe_deps.py`、
  `audit_deps.py`、`cleanup_runs.py`、`convert_genus_lens.py`、
  `merge_genus_lens.py` —— 可重复使用的打包/审计/数据处理工具。
