# VirusPlatform 全面审查报告（2026-10-07）

> 审查对象：`Virus_Platform_Core/` 全部 127 个模块 + `main.py` / `app.py` / `webapp/` / `tests/`
> （核心包 88,631 行；含 10-06 新增、此前从未审查的 CRISPR 模块 2,217 行）
>
> 审查方式：全量语法编译 → 静态体检 → 主编排器精读 → 6 路并行只读代码审查 →
> 文档/代码量化对账 → **每条 P0/P1 亲自复核 + 造数据实测**。
>
> **本次未改动任何源文件**（`git status` 中 `Virus_Platform_Core/` 的 M 标记均为作者既有未提交改动；
> 会话内仅测试运行器自身写入 `tests/_run_all_report.txt`。探针临时脚本全部落 `%TEMP%` 并已删除，
> 探针产生的 `run/_probe2.log` 与空的 `run/tests_crispr/` 已清理）。
>
> 标注约定：✅ 已实测（附复现命令/输出） · 推断（深读代码得出，未跑）

---

## 0. 总体评价

**建议决定：小修（Minor Revision）。**

上一轮（2026-10-04）的 M1–M4 与三轮 P3 清零确实落地了，不是纸面承诺：我复跑了快速守卫子集 **9/9 通过**，
语法编译 **295 文件 0 错误**，`crispr_design` 单测 **18/18 通过**。工程纪律仍是同类项目里少见的：
`shell=False` 全程、`check_path` 集中防穿越、进程树统一终止（`utils._kill_proc_tree` 已替换三处裸 `kill`）、
`_run_all` 独立超时+互斥锁+陈旧锁回收、冻结分发"平台根 = exe 目录"这条硬约束被反复写进注释。

本轮发现的问题集中在**一个系统性模式**上，而不是散落的 bug：

> **同一语义存在多份实现时，修复只落在其中一处，同构路径被漏掉。**

三处独立证据都指向它：
- `.` 是否算 gap —— `phylo.py:605` 修了（注释明写"旧实现只认 '-'"），`sdt_exact.py:152` **没修**；
- `databases` 根该走 `DIRS` —— `tools_api.py:567` 知道（注释写了"分发版下两者不同"），同文件 `:630` **没改**；
- `crispr_design.py` 同一函数内 `:1851` 用了 `DIRS`，`:1874` **没用**。

这类"半修复"比完全没修更危险：它让 `grep 到一处正确写法` 变成不可信的通过信号。
建议把这三条（gap 字符集、数据库根、纯 ACGT 守卫）各自收敛成**单一权威实现**，而不是逐点补。

第二个模式是**守卫的覆盖盲区**：项目已经建立了很好的"发现一类问题就加一个守卫"的习惯，
但本轮实测出三个守卫的实际覆盖面小于它给人的印象（见 P1-2 / P1-3 / P1-4）。

---

## 1. P0 — 会导致结论错误或静默失效

### P0-1. SDT 的 gap 字符集与 identity 页不一致，最大差 9.72 个百分点 ✅已实测

三套口径并存：

| 位置 | 判据 | 说明 |
|---|---|---|
| `sdt_exact.py:152-157` `_sdt_formula` | `if x == '-' or y == '-'` | **只认 `-`** |
| `phylo.py:605-607` `pairwise_identity_matrix` | `if x not in '-.' and y not in '-.'` | 认 `-.`，注释明写"旧实现只认 `-`，把 `.` 当成碱基位点参与比较（2026-10-03 修复）" |
| `phylodyn_alignqc.py:35` `_GAP_CHARS = '-.?'` | 认 `-.?` | 又多一个 `?` |

实测（直接调 `_sdt_formula` 与 phylo 口径对拍）：

```
用例                 sdt_exact    phylo口径      差
纯 '-' gap          87.5000      87.5000      0.0000
'.' 作为 gap         88.8889      87.5000      1.3889
混合 . 与 -           87.5000      87.5000      0.0000
一侧 . 一侧碱基        77.7778      87.5000     -9.7222   ← 一侧 '.' 一侧有碱基
```

**为什么会改变结论**：SDT 是平台 ⑦ 步的核心输出、物种划界的直接依据（多数病毒属阈值 80%）。
9.72 个百分点足以把一个"应判同种"的对推到阈值之下。同一份含 `.` 的已比对 FASTA
（BioEdit / MEGA / Clustal 导出的常见形状）在 `/sdt` 页与 `/compare` 页会给出**两个答案**。

**触发条件（重要）**：平台的 MAFFT/trimAl 路径产出的 gap 是 `-`，**平台自产数据不触发**；
触发路径是用户上传的已比对 MSA（`tool_jobs.py:1237-1288` 的 `aligned=True` 分支，明确支持）。
`examples/` 下无此形状文件（已扫，见"误报排除"）。

**修法建议**：`_sdt_formula` 的 gap 判据改为 `x in '-.'`，与 `phylo.py` 对齐；
`?` 是否纳入属**判定语义**，请老师定方向后再动（`phylodyn_alignqc` 认它，`phylo` 不认）。

---

### P0-2. CRISPR 路线 B 漏出 IUPAC 兼并码引物 ✅已实测

`crispr_design.py:699-713` 的 `_ok_primer` 只挡字面 `N`：

```python
def _ok_primer(seq, c0):
    if 'N' in seq:            # ← 只挡 N
        return False
```

而同文件另两处守卫是完整的：路线 A（primer3）`:657-660` 用 `set(...) - _PURE_ACGT`，
探针 `:1667` 同款。作者在 `:352-355` **明确声明了契约**：

> 平票不猜：2 碱基平票记 IUPAC 兼并码（W/S/M/K/R/Y…）
> **crRNA/引物序列的纯 ACGT 约束由候选扫描与引物后置守卫另行把关**

即：路线 B 违反了作者自己写下的契约。

实测（造 600bp 窗口、侧翼埋 15 个 `W`，调 `_enumerate_primers`）：

```
路线 B 产出引物对数: 5
含非 ACGT 碱基的引物对: 5
  r=CAATGATTWWGAGCGCTAWWACATTATWWWTTAG (含['W'])
  ...
```

**为什么会改变结论**：候选池里混入含兼并码的"引物"，会一路流到 `ranked_sets.tsv` 与
`order_sheet.xlsx`（全仓 `_PURE_ACGT` 只在 `:388/465/477/487/495/659/1667` 出现，**无下游再过滤**），
用户照订购表下单会拿到无法合成的序列。

**触发条件**：共识序列在 2 碱基平票处记兼并码（`:352-356`，**与严格/放松档无关**，是默认行为）；
且 primer3 失败或返回 < `min_pairs_per_site`（默认 2）→ 走枚举兜底。两个条件都不罕见。

**修法建议**：`_ok_primer` 改 `if set(seq) - _PURE_ACGT: return False`。一行。

---

## 2. P1 — 健壮性 / 可移植性 / 可复现性

### P1-1. 「排队中取消」泄漏样品占位 → 样品永久不可再跑 ✅已实测

链路：
- `samples.py:571` 调 `_sample_busy_register(sample)` 登记占位（同理 `:201`、`:913`）；
- 释放只写在 job 函数的 `finally` 里（`:582-583`、`:213`、`:922-923`）；
- `tasks.py:229-236` 在"排队途中被取消"分支里直接 `return`，**从不调用 `fn`**。

于是 job 的 `finally` 永不执行 → `_RUNNING_SAMPLES` 里的样品名永不摘除。

实测（探针直调 `TaskManager`，heavy 名额压到 1 制造排队）：

```
排队中状态: running
queued_job 是否被调用 : False        ← fn 从未执行
release(finally) 是否执行: False      ← 清理从未执行
取消后状态: cancelled | error: 任务已停止（排队中取消）
```

**可达性**：heavy 名额默认 2（`tasks.py:44`），样品分析任务默认 heavy（`samples.py:220/590/932` 均未传 `weight`）。
跑 2 个样品 → 第 3 个排队 → 取消它 → 该样品被永久标记占用。

**后果**：`_sample_busy()`（`samples.py:966-978`）此后恒返回 True → 重新入队、删除样品、
清空结果全部被拒（提示"已有任务在运行"），且**无任何回收机制**（`_RUNNING_SAMPLES` 全仓只出现在
`323/329/337/970` 四处，无 TTL、无对账），必须重启进程。

**修法建议**：把释放挂到 TaskManager 的终态回调（而不是 job 的 `finally`），
或在 `_sample_busy()` 里对已终态任务做自愈清理。

---

### P1-2. 数据库根硬编码：分发版（数据库外置）下读错目录 ✅代码位置已实测 / 分发版后果为推断

| 位置 | 写法 | 同文件/同函数的正确对照 |
|---|---|---|
| `web/tools_api.py:630-631` | `os.path.join(PLATFORM_ROOT, 'databases', 'dsrna', 'panel')` | `:559-571` 用 `db_path()`/`DIRS['databases']`，`:567` 注释明写"分发版下两者不同（`DIRS['databases']` = `<数据库包>\databases`）" |
| `web/tools_api.py:778` | `.../databases/mirna_lib/plant_mirna_library.tsv` | 同上 |
| `crispr_design.py:1874` | `.../databases/crispr/near_neighbors` | 同函数 `:1851-1852` 已正确用 `DIRS.get('databases')` |
| `crispr_design.py:1905` | `.../databases/crispr/blast_cache` | 同上 |

**为什么改变结论**：分发版由 `config._detect_database_root()` 把数据库根重定向到同级
`VirusPlatform-Database`，两者不等价。后果分别是：
- `/api/tool/dsrna_panel`：代码把"目录缺失"当正常（docstring 明写"返回空清单，不当作错误"）
  → 面板清单恒为空，用户以为平台没有面板数据 —— **静默失效**；
- miRNA 库：`abort(500, '候选 miRNA 库缺失…')` —— 显式失败。

**守卫盲区**：`tests/_check_db_paths.py` 抓不到它。它只探 3 个端点
（`/api/dbs`、`/api/host_dbs`、`/api/kv_index_list`），且判据是"**返回给前端的路径必须绝对**"——
本问题是"**后端从错误的根读盘**"，性质不同，响应式判据原理上覆盖不到。

**修法建议**：统一改 `DIRS['databases']`；守卫加一条 **AST 扫描**（禁 `PLATFORM_ROOT, 'databases'`
字面组合，白名单：`config.py:58-61` 定义处、`db_migrate.py:139/145-147` 的迁移源目录）。

---

### P1-3. 86 个 pytest 用例不在任何自动入口 ✅已实测

`tests/_run_all.py` 是**硬编码清单**（`TASKS` 元组，`:105-312`），不是 glob。
它只收 `_it_*/_check_*/_smoke_*/_audit_*` 与两个脚本式 `test_logan_*`。

以下 6 个文件是 **pytest 风格**（用 `pytest.fixture`/`importorskip`，**无 `__main__` 块**，
因此也不能当脚本跑），既不在清单里、也没有 CI（无 `.github/workflows`）：

| 文件 | 用例数 |
|---|---|
| `tests/test_crispr_design.py` | 18 |
| `tests/test_host_genome_gget.py` | 16 |
| `tests/test_mirna_alignment.py` | 17 |
| `tests/test_utils_pure.py` | 14 |
| `tests/test_datasets_cli.py` | 11 |
| `tests/test_web_guard.py` | 5 |
| **合计** | **81（+ CRISPR 18 = 86）** |

实测：手动 `pytest` 跑这 6 个 → **86 passed**（CRISPR 18 单独跑 91.9s 全绿）。

**为什么是问题**：平台最新、最大的新增模块（`crispr_design.py`，2,217 行）的回归
**完全依赖人工记得敲 pytest**。项目自己在 `docs/review_VirusPlatform_20261004.md` 附录里写过
"一个没人跑的检查等于没有"——现在这个教训落在了新模块自己身上。

**修法建议**：`_run_all.py` 增加一个 pytest 分组（对 `test_*.py` 走 `python -m pytest`，
复用现有超时/互斥锁），或把这 6 个文件补上 `__main__` 入口改成脚本式。

---

### P1-4. SVG 良构守卫只覆盖 dfv_plot + `examples/` 快照 ✅调用点已实测 / "会复现"为推断

- `utils.verify_svg` 全仓只被 `dfv_plot.py:158` 与 `:184` 调用（生成时自检）。
- 其余 SVG 生产者**无生成时自检**：`gbdraw_plot.py`（圈图/线性图）、`viz.py`（圈图）、
  共线性比对图、`orf_annot` 的 `genome_diagrams/*.svg`。
- `tests/_check_svg_wellformed.py` 只 `os.walk(os.path.join(ROOT, 'examples'))` —— 这是**快照**判据，
  `run/tool_runs/**` 下新生成的产物不在覆盖内。

这正是上一轮 M3 的病灶：坏掉的 `examples/results/synteny/compare.svg`
（mtime 10-06 15:51，是**事后修补**的产物；18 处 `</path` 缺 `>`）的生成端没有加守卫。
即：**同一段生成逻辑再跑一次，仍可能产出坏 SVG，而守卫看不见它落在哪。**

**修法建议**：把 `verify_svg` 提到统一 savefig 出口；或把守卫扫描范围扩到 `run/tool_runs/**/*.svg`。

---

### P1-5. SDT 分布图在"全同输入"下崩溃，且崩溃点在产物半写之后 ✅已实测

`sdt_exact.py:508-512`：

```python
min_v = max(0.0, float(np.floor(np.nanmin(vals))))
counts, edges = np.histogram(vals, bins=np.arange(min_v, 101, 1))
plt.ylim(0, max(counts / len(vals)) * 1.1)      # ← 全同输入时这里炸
```

全部两两 identity == 100 → `min_v = 100` → `np.arange(100, 101, 1)` 只有一个边界
→ `counts` 为空数组 → `max(空)`。

实测：`counts = []  len = 0` → `ValueError: max() iterable argument is empty`。

**为什么是 P1 而不是 P2**：`:636` 的调用处**无 try/except**。而 `:627` 已写完
`sdt_matrix.csv`、`:634` 已写完热图、`:647` 的 `sdt_matrix.json` **还没写** ——
结果是"矩阵产物躺在目录里、任务却报失败"这一类最难排查的形态。

**触发条件**：输入含 ≥2 条完全相同的序列。这不是假想 —— 上一轮审查（M1 之前）
就在 `example_conserved_set.fasta` 里实测到过"其中 2 条完全相同"。

**修法建议**：`min_v` 与 `100` 重合时退化为单 bin（如 `np.arange(min_v - 1, 101, 1)`），
或给 `ylim` 加空数组保护。

---

## 3. P2 — 死代码 / 重复实现 / 遗留物

| 项 | 证据 | 说明 |
|---|---|---|
| **文档数字陈旧** | `wiki/Architecture.md:13,16,39` 写"70+ 模块 / 30+ 子命令"；实测 `main.py` 的 `add_parser` 计数 **41**，`Virus_Platform_Core` 顶层+子包模块 **127** | 数字不算错（41≥30、127≥70），但已明显脱离现实 |
| **CLI 文档漏 1/41** | 41 个子命令中，`crispr-design` 是**唯一**未出现在 `wiki/CLI-Reference.md` 的（`grep -c crispr` = 0）；`README.md`/`README_EN.md`/`wiki/Home.md` 完全未提 CRISPR，仅 `_Sidebar.md` 有 1 处 | 规划文档 §4.1 明确要求"对齐 CLI-Reference 页"，未完成 |
| **`pct_mean` 混源不可解释** | `public_meta/virome.py:416` 对每样本 pct 求均值；`:133` docstring 已警告"NCBI 用 analysed(已鉴定 reads)、GSA 用 analysed(总 reads) → 跨来源不可直接比" | 作用域**精确**：只影响 `by_virus` 的 `pct_mean` 列；`occ`/`n_samples` 不受影响；`sources` 计数已返回（`:395/419`），前端可自行分层。**不是**"全管线口径不统一" |
| **`.bak` 源码树残留 17 个** | `Virus_Platform_Core/public_meta/__init__.py.bak_20261001_prescreen`、`tests/_run_all.py.bak_20260916_geo移植前`、`tests/_it_platform.py.bak_*` ×6、`tests/_check_phylodyn_ui.py.bak_*` ×4 等 | 上轮归档了 4 个到 `archive/_bak_20261004/`，其余仍在原位。`_check_spec_datas` 现为 PASS（不影响打包），属卫生问题 |
| **`kv_engines.py` 15 处裸 subprocess** | 上轮已补 6 处 `timeout=7200`，管道 Popen 依赖外层 6h + 进程树兜底 | 属已知设计取舍，保留观察 |

---

## 4. 文档与代码不一致（量化对账）

| 文档承诺 | 代码实际 | 判定 |
|---|---|---|
| `wiki/Architecture.md` "30+ 子命令" | 41（`main.py` `add_parser`） | 陈旧但不算错 |
| `wiki/Architecture.md` "70+ 模块" | 127 | 陈旧 |
| `wiki/CLI-Reference.md` 应覆盖全部子命令 | 41 中漏 `crispr-design` | ❌ 不一致 |
| 规划文档 §5「二期全部落地：NucleaseProfile + Cas13a / EXO-NFO 探针 / intentional mismatch / 跨片段 panel / IUPAC 兼并码共识；测试 18 用例」 | `crispr_design.py` 确有对应实现；`pytest` 实测 **18 passed** | ✅ 一致（但该测试不在自动入口，见 P1-3） |
| `crispr_design.py:352-355` "引物序列的纯 ACGT 约束由引物后置守卫另行把关" | 路线 A ✅ / 探针 ✅ / **路线 B ❌** | ❌ 不一致（P0-2） |
| `phylo.py:605` 注释"gap 判定统一 `in '-.'`" | `phylo.py` ✅，`sdt_exact.py` 仍是 `'-'` | ❌ "统一"名不副实（P0-1） |
| `sdt_exact.py:5-13` 模块头 "口径与 SDT 原版一致：gaps = 任一序列为 gap 的列数" | 只认 `-`；`.` 不算 gap | ❌ 与自身声明不符（P0-1） |
| `tools_api.py:567` 注释"分发版下两者不同" | 同文件 `:630/778` 仍用 `PLATFORM_ROOT` | ❌ 不一致（P1-2） |

---

## 5. 误报排除（查证后**不报**，留档防止下次重犯）

按流程，子代理给的每条 P0/P1 都亲自复核过。以下 4 条经查证**不成立**，明确排除：

1. **「`crispr_design` 未登记 `spec hiddenimports` → 冻结分发缺模块」** —— 不成立。
   系统扫描发现**50 个顶层 Core 模块**同样是函数内惰性 import 且未登记
   （`assembly` / `consensus` / `sdt_exact` / `viz` / `pipeline` …），而分发版长期正常。
   说明作者在 spec 里显式登记的是**子包内相对导入**那类静态分析不到的情形
   （`spec:141-155` 注释即如此说明）。`from Virus_Platform_Core import crispr_design`
   是绝对导入，PyInstaller 能分析到。**结论：非缺陷。**
2. **「`treetime_ml.py:122` 硬编码 `C:\Python312`」** —— 不成立。
   `_PREFERRED_PY` 是**探测链的首选**，`:185` 明确 fallback 到 `sys.executable`；
   同理 `suvtk_submit.py:181` 是候选列表且有 `import suvtk` 校验。属**有意策略**。
3. **「`examples/` 含 `.` 作 gap 的比对文件」** —— 不成立。扫出的命中全部来自
   FASTQ（质量串含 `.`）与 `.fai` 索引，非比对序列。P0-1 的触发面**只在用户上传路径**。
4. **「日期口径全管线未统一 / `virphy_bridge.py` 保留年初」** —— 不成立。
   `Virus_Platform_Core/virphy_bridge.py` **不存在**（已重构），该说法来自过期的项目记忆。
   现存两处实现（`phylodyn_govern.decimal_year:120` 与 `phylodyn_local.to_decimal_year:142`）
   **默认 `mode='mid'` 一致**，且后者 docstring 明确声明"唯一权威实现见 `parse_date`"，
   并有 2026-09-18 / 2026-10-03 两次口径修复记录。**作用域清晰，不构成"不统一"。**

---

## 6. 修复优先级建议

**第一档 — 会改变结论（建议本轮修）**
1. **P0-1** `sdt_exact.py:152-157` gap 判据改 `x in '-.'`；`?` 是否纳入请老师定方向。
   → 配一个对拍守卫：同一份含 `.` 的 MSA 走 `/sdt` 与 `pairwise_identity_matrix`，断言逐对相等。
2. **P0-2** `crispr_design.py:699` 补 `set(seq) - _PURE_ACGT`。
   → 复用已实测的探针（600bp 窗口 + 侧翼 `W`）作为守卫用例。

**第二档 — 部署运维（分发版才会暴露，但一旦暴露很难查）**
3. **P1-2** `tools_api.py:630/631/778`、`crispr_design.py:1874/1905` 统一 `DIRS['databases']`；
   `_check_db_paths.py` 增一条 AST 判据。
4. **P1-1** 样品占位释放挂到任务终态回调。
5. **P1-5** SDT 分布图 bins/ylim 空数组保护。

**第三档 — 可维护性（防回归）**
6. **P1-3** `_run_all.py` 接入 pytest 分组（86 个用例重新回到回归面）。
7. **P1-4** `verify_svg` 提到统一 savefig 出口 / 扩大守卫扫描范围。
8. **P2** `wiki/Architecture.md` 与 `CLI-Reference.md` 补 `crispr-design` 与新模块数；
   清理 17 个 `.bak` 残留。

---

## 7. 值得保留的设计

- **冻结分发的"平台根"纪律**：`config.PLATFORM_ROOT` 单一来源 + `check_path` 集中防穿越 +
  `DIRS['databases']` 重定向，注释里逐条记着"已踩三处"的具体教训。这类约束在同类项目里极少见。
- **测试运行器的工程化**：独立超时（区分 TIMEOUT 与 FAIL）、互斥锁 + 陈旧锁回收、
  `VP_TEST_WAIT` 环境变量给批量满载留余量 —— 注释里记录了历史事故根因，不是空壳。
- **"示例即真实链路"**：`tests/make_example_results.py` 走与用户点击完全相同的 `/api/tool/run`，
  示例产物因此具备对账价值。
- **守卫驱动开发的习惯**：`_check_identity_aa` / `_check_genome_circos` / `_check_svg_wellformed`
  / `_check_primer` 都是"发现一类问题 → 加一个守卫"的产物，方向完全正确
  （本轮的三条覆盖盲区是在这个基础上的增量建议，不是否定）。
- **SDT 的 AA 锚定翻译**（`sdt_exact.py:707` `_anchored_aa_seqs`）：分簇 → 簇内最长 ORF 定 (链,框)
  → 整框宽容翻译（终止→X 不切断）→ 多轮重簇锚定 → 末位回退。修法比"加个告警"扎实得多，
  且把 NT 矩阵前置为分簇依据，逻辑自洽。
- **`_stage_execution_order` 的拓扑排序**：用 `STAGE_ORDER` 兜底定序，
  注释里给了"传 ['report','primer',…] 会得到报告先于组装"的具体反例。防御性写得很到位。

---

*审查产物：本报告。探针脚本均为临时文件（`%TEMP%`）已删除；
批量守卫子集结果见 `tests/_run_all_report.txt`（9/9 通过）。*
