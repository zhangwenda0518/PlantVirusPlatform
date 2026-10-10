# VirusPlatform 审查问题修复说明（2026-10-07）

> 配套文档：[审查报告 `review_VirusPlatform_20261007.md`](review_VirusPlatform_20261007.md)
>
> 本轮对审查报告 §6「修复优先级建议」的第一/二/三档全部条目执行了修复，
> 每条修复都配了**可独立复跑的验证脚本**，且带**对照臂**（从 `archive/_bak_20261007/`
> 原样抽取旧实现，同一输入 A/B 对拍）证明旧行为确实错。
>
> 标注约定：✅ 已实测（附脚本与实测数字） · ⏸ 有意未改（附理由与需老师定的方向）

---

## 0. 一句话总结

报告里 8 条待修项（P0-1 / P0-2 / P1-1…P1-5 / P2）**全部落地**；
其中 P1-4 与 P2 在实施时又各自暴露出**一处更深的问题**（守卫覆盖面小于字面、不可修复的残图仍被内嵌），
已一并修掉。本轮**新增守卫 1 个**（`tests/_it_svg_embed_chain.py`）、**扩面 2 个**
（`tests/_check_svg_wellformed.py` 扫到 `run/`；`tests/_check_spec_datas.py` 的 bak 断言扫到全树），
批量入口从 62 项扩到 **73 项**（`tests/_run_all.py --list | grep -c timeout=` 实测 73）。

**涉及判定语义的改动只有 1 处待定**：`?` 是否纳入 SDT 的 gap 字符集（见 §4.1），
本次**未动**，`_SDT_GAP` 常量处已就地标注待定。

**另有 1 项未完成**：`tests/_check_identity_aa.py` 与 `tests/_check_primer.py`
因 MAFFT 在本会话环境里不退出而**未能跑完，应视为未验证**（见 §7）——
该现象与本轮改动无关，但排查中发现了「`subprocess` 超时对 MAFFT 挂起无效」这个真实缺口。

---

## 1. P0 — 会改变结论的两条

### 1.1 P0-1 SDT gap 字符集与 identity 页对齐 ✅

**改动**：`Virus_Platform_Core/sdt_exact.py`

```python
_SDT_GAP = '-.'                      # 新增模块级常量（单一权威定义）

def _sdt_formula(a, b):
    ...
    if x in _SDT_GAP or y in _SDT_GAP:      # 原：if x == '-' or y == '-'
        gaps += 1
```

与 `phylo.py:605` 的 `pairwise_identity_matrix` 口径对齐（`-.` 都算 gap）。

**验证脚本**：`tests/_check_sdt_gap_parity.py`（三条判据）

| 判据 | 内容 | 实测 |
|---|---|---|
| A 一致性 | `_sdt_formula` 逐对 == `pairwise_identity_matrix` | 通过 |
| B **灵敏度（对照臂）** | 从备份抽取的旧实现，含 `.` 输入上**必须**与 phylo 不同 | `s1 vs s2: 旧 83.3333 → 新 90.9091（差 +7.5758 个百分点）` |
| C 改动面受限 | 无 `.` 的输入上新旧实现零差异 | 通过 |

判据 C 很关键：它把这次改动**限制在"含 `.` 的输入"**上，
证明平台自产数据（gap 恒为 `-`）的行为完全没变。

---

### 1.2 P0-2 CRISPR 路线 B 的引物纯度 ✅

**改动**：`Virus_Platform_Core/crispr_design.py`

```python
def _ok_primer(seq, c0):
    if set(seq) - _PURE_ACGT:        # 原：if 'N' in seq:
        return False
```

与同文件路线 A（`:657-660`）与探针（`:1667`）的口径一致，
即兑现 `:352-355` 自己写下的契约"纯 ACGT 约束由引物后置守卫另行把关"。

**验证脚本**：`tests/_check_crispr_primer_purity.py`

| 判据 | 实测 |
|---|---|
| A 新实现：600bp 窗口（侧翼埋 15 个 `W`）产出引物全为纯 ACGT | 通过 |
| B **对照臂**：旧实现漏出的含兼并码引物对 | **5 对**（应 >0，否则测试恒过） |
| C 纯 ACGT 窗口仍能正常产出（不误伤） | 通过 |

---

## 2. P1 — 健壮性 / 可移植性 / 可复现性

### 2.1 P1-1 「排队中取消」泄漏样品占位 ✅

**改动**：3 个文件

1. `Virus_Platform_Core/web/tasks.py` — `TaskManager.start()` 新增 `on_finish=None` 形参；
   `_run()` 内新增**幂等**的 `_call_on_finish()`，在「排队取消」分支与 `finally` 各调一次。
2. `Virus_Platform_Core/web/samples.py` — 三处 `tm.start(...)`（`_run_entry` / `api_analyze` /
   `api_pipeline_run`）都补 `on_finish=lambda: _sample_busy_release(s)`。
3. 幂等保护（`_fin = {'called': False}`）是为了防「排队取消分支调了、`finally` 又调一次」造成重复释放。

**为什么挂在终态回调而不是 job 的 `finally`**：`tasks.py` 在"排队途中被取消"时**从不调用 `fn`**，
所以 job 里的 `finally` 永不执行 —— 这正是泄漏的根因，把清理写回 `finally` 无法修好。

**验证脚本**：`tests/_check_task_on_finish.py`

| 判据 | 实测 |
|---|---|
| A **对照臂**：旧 TaskManager 在「排队取消」时 `fn` 的 `finally` 执行次数 | **0** |
| B 新实现：同一路径 `on_finish` 调用次数 | **1** |
| C 正常完成路径 | **1** |
| D 任务抛异常路径 | **1** |

脚本顶部置 `os.environ['VP_NO_RECOVER'] = '1'`，避免测试触发 TaskManager 的恢复扫描而污染 `run/`。

---

### 2.2 P1-2 数据库根硬编码 ✅

**改动**：4 个文件、6 处

| 文件 | 改动 |
|---|---|
| `web/tools_api.py` | `panel_dir` / `lethal_dir`（原 `:630-631`）、`_mirna_lib_path()`（原 `:778`）改 `DIRS['databases']`；删掉死回退的 `or os.path.join(PLATFORM_ROOT, 'databases')` |
| `crispr_design.py` | `index_cache` / `near_neighbors` / `blast_cache` 三处改 `DIRS['databases']`；导入行删掉 `PLATFORM_ROOT` |
| `web/examples_api.py` | `_ex_in()` 改 `DIRS['databases']` |
| `web/settings.py` | `'virus_src'` 回退分支改 `DIRS['databases']` |

后两个文件本不在报告的 P1-2 清单里 —— 改它们是为了让新守卫能定成**严格规则**
（"全仓不得出现 `os.path.join(PLATFORM_ROOT, 'databases', …)`"），
否则守卫必须开白名单，而白名单会掩盖同类新增。

**验证脚本**：`tests/_check_db_paths.py`（扩面）

新增 **AST 判据**（不是文本搜索 —— 注释里就写着这些路径）：

```
_is_platform_databases_join()   # 判 os.path.join(PLATFORM_ROOT, 'databases', ...)
_find_in_source() / _scan_platform_databases() / _selftest()
```

白名单只有 2 个（定义处）：`Virus_Platform_Core/config.py`、`db_migrate.py`。

| 判据 | 实测 |
|---|---|
| 全仓 AST 扫描 | **136 个文件，命中 0 处** |
| 自测①：注入代码形式（含 `'databases'` 字面量）必命中 | 1 处 ✅ |
| 自测②：注释 / 字符串形式必须 0 处 | 0 处 ✅ |
| 自测③：真实 `config.py` 必能被扫出（防"判据写错导致恒不匹配"） | >0 处 ✅ |

> ⚠️ 自测③ 是踩坑后补的：本守卫第一版把 AST 判据写成
> `node.func.value` 是 `Name('os')`，而 `os.path.join` 的真实 AST 是
> `Attribute(attr='join', value=Attribute(attr='path', value=Name('os')))`
> —— 于是判据**恒不匹配**，守卫**静默失灵**（假通过）。这正是 skill 铁律一
> "判据本身错了会带来假通过"的实例。

---

### 2.3 P1-3 86 个 pytest 用例回到回归面 ✅

**改动**：`tests/_run_all.py`

- `TASKS` 元组末尾接入 6 个 pytest 文件；
- 新增 `PYTEST_FILES` 集合 + `_pytest_available()`（探测一次并缓存）；
- `_run()` 按 `is_pytest` 分支走 `[sys.executable, '-m', 'pytest', p, '-q']`；
- 未装 pytest 时**整组 SKIP**（不是 FAIL —— 开发依赖未装的机器不该因此报红）；
- `n_skip` 计数并修正 total。

| 项 | 实测 |
|---|---|
| 批量入口条目数 | 62 → **73**（`--list` 实测） |
| 本轮新增/扩面的守卫 | 8 个（4 个上一轮新增 + 1 个本轮新增 + 3 个扩面）**全部登记在清单内** |
| pytest 分组单独跑 | **86 passed in 72.90s** |

---

### 2.4 P1-4 SVG 良构守卫扩面 ✅（实施中发现更深一层问题）

**改动**：4 个文件

1. `Virus_Platform_Core/gbdraw_plot.py` — 产物收集后对 `.svg` 调 `verify_svg`。
   （gbdraw 是外部工具，产物可能出 `</path` 缺 `>`；此前只有 `dfv_plot` 挂了守卫。）
2. `Virus_Platform_Core/genome_diag.py` — `render_contig_svg` 落盘后调 `verify_svg`。
3. `Virus_Platform_Core/viz.py` — **消费侧**兜底：内嵌进 `report.html` 之前先 `verify_svg`，
   **先修后读**；`verify_svg` 返回 False（修不动）时**跳过嵌入**并计 `n_broken` 告警。
4. `tests/_check_svg_wellformed.py` — 扫描范围从 `examples/` 扩到 `examples/ + run/`，
   并加 5 条判据。

**实施中发现的更深问题（原报告只提到"生成端没有守卫"）**：
只加生成端守卫**不够** —— 消费端（`viz.build_report`）会**静默读入**任何残缺 SVG 并内嵌，
浏览器容错渲染成残图，肉眼几乎看不出。所以"防静默读入"必须落在消费侧。
这正是上一轮 M3 病灶（`examples/results/synteny/compare.svg` 18 处残缺）之所以能进交付物的链路。

**验证脚本**：两个

`tests/_check_svg_wellformed.py`（5 条判据）

| 判据 | 内容 | 实测 |
|---|---|---|
| A 覆盖范围 | `examples/` + `run/` 全部 SVG 良构（只读，不改文件） | 扫描 **39 个**，非白名单残缺 **0** |
| B **灵敏度自测** | 残缺 / 良构 fixture 用同一判据必须被区分 | 残缺被判坏、良构未误判 ✅ |
| C 修复路径 | `verify_svg` 幂等自修仍可用 | 返回 True、修复后可解析、二次调用仍良构 ✅ |
| D 白名单不腐烂 | `KNOWN_BAD` 每条必须**仍然**残缺 | 2/2 ✅ |
| E 出口覆盖（AST） | `dfv_plot` / `gbdraw_plot` / `genome_diag` / `viz` 必调 `verify_svg`；且 `viz.py` 中 `verify_svg(svg_p)@840` 必须早于 `safe_open(svg_p)@841` | 全通过 ✅ |

`tests/_it_svg_embed_chain.py`（新增，3 条臂，端到端真跑 `build_report`）

| 臂 | 内容 | 实测 |
|---|---|---|
| 前置对照 | 两个坏样本读原文确实非良构（否则测试恒过） | ✅ |
| A 修得动 | 闭合标签缺 `>` → 就地修复、以良构片段内嵌（标记 `gfix` 可见） | ✅ |
| B 修不动 | 开标签也没闭合 → **不得内嵌**（标记 `gbad` 不可见） | ✅ |
| C 本来就良构 | 原样内嵌（标记 `gok` 可见） | ✅ |
| 汇总 | 报告内 SVG 片段非良构数 | **2 个片段 / 0 个非良构** |

灵敏度实证：把 `viz.verify_svg` 换成 no-op（等价于修复前行为）重跑该测试 →
**4 条 FAIL、退出码 1**（`报告内不含残缺字节串` / `非良构片段 = 2` /
`修不动的 SVG 已被跳过` / `落盘文件已修复`），证明它不是恒过。

> 附一条自己的踩坑：守卫第一版的"残缺 fixture"写成 `<path d="M0"</path` ——
> 开标签也没闭合，`verify_svg` 的正则**修不动**，于是判据 C 报红。
> 真实历史缺陷是"开标签合法、仅闭合 `</path` 缺 `>`"，fixture 必须同形才能测出修复路径。

---

### 2.5 P1-5 SDT 分布图全同输入崩溃 ✅

**改动**：`Virus_Platform_Core/sdt_exact.py` 的 `plot_distribution`

```python
min_v = max(0.0, float(np.floor(np.nanmin(vals))))
if min_v >= 100.0:                      # 全同输入：退化为单 bin
    edges = np.array([99.5, 100.5]); x_lo, x_hi = 99.5, 100.5
else:
    edges = np.arange(min_v, 101, 1); x_lo, x_hi = min_v, 100.0
counts, edges = np.histogram(vals, bins=edges)
prop = counts / len(vals)
plt.plot(0.5 * (edges[:-1] + edges[1:]), prop, color='#CC0000', lw=3.0)
plt.ylim(0, float(np.max(prop)) * 1.1)   # 原：max(counts/len(vals)) * 1.1
...
plt.xlim(x_lo, x_hi)
```

**验证脚本**：`tests/_check_sdt_distribution.py`

| 判据 | 实测 |
|---|---|
| A 新实现：全同矩阵（全 100）能正常出图（png + pdf） | 通过 |
| B **对照臂**：旧实现抛异常 | `ValueError: max() iterable argument is empty` |
| C 常规矩阵 Proportion 合计 == 1.0（钉住 2026-10-04 的 D3 修复，防回退） | 通过 |

---

## 3. P2 — 文档 / 遗留物

### 3.1 文档数字与 CLI 覆盖 ✅

| 文件 | 改动 |
|---|---|
| `wiki/Architecture.md:13` | `30+ 子命令` → **41 子命令** |
| `wiki/Architecture.md:16` | `70+ 模块` → **127 模块**；模块清单补 `crispr` |
| `wiki/Architecture.md:39` | `70+ 模块` → **127 模块** |
| `wiki/Architecture.md`（分层架构下方） | 新增「规模数字口径」条目：写明 127 的算法（顶层非 `__init__`）、41 的算法（AST 抽 `add_parser` 首参），并注明"不要目测" |
| `wiki/CLI-Reference.md` | 新增 `## CRISPR 检测设计` 段，补 `crispr-design` 全部参数（含默认值、档位语义） |

**核验口径**（可复跑）：

```
main.py add_parser 子命令数 = 41          （AST 抽首个常量实参）
Virus_Platform_Core 顶层非 __init__ = 127 （另含 7 个子包 __init__，共 134 个 .py / 88,730 行）
文档未提及的子命令 = 无（41/41 全覆盖）
```

> 说明：报告里写的"88,631 行"是审查当时的快照；本轮修复使行数增至 88,730。

### 3.2 `.bak` 残留归档 + **守卫扩面** ✅

**改动**：17 个 `.bak*` 从源码树迁移到 `archive/_bak_residue_20261007/`（保留相对路径）。

**实施中发现的更深问题**：`tests/_check_spec_datas.py` 里那条
`check('源码树已无 .bak 手工备份…')` 的标签写着"源码树"，
但它的 `bak_src` 来自 `os.walk(webapp)` —— **实际只扫 `webapp/`**。
所以源码树里 17 个 `.bak*`（`tests/` 12 个、`Virus_Platform_Core/` 2 个、`docs/` 1 个…）
长期存在而这条断言**始终绿灯**。这是"守卫覆盖面小于字面印象"的又一实例。

**改动**：该断言改为扫**整棵源码树**（排除 `archive/ 3rd/ dist/ run/ .git/ __pycache__/ …`），
并新增两条**防假绿**的伴随断言：

- **范围自证**：扫描必须真的覆盖 `Virus_Platform_Core` / `tests` / `docs` / `webapp` / `dev_tools`，
  否则"已无残留"是假绿；
- **灵敏度自测**：临时目录里放一个 `.bak`，必须被扫到。

| 判据 | 实测 |
|---|---|
| 源码树 `.bak` 残留 | **扫 1032 个文件，命中 0** |
| 范围自证 | 已覆盖 13 个顶层目录，未覆盖 `[]` |
| 灵敏度自测（临时目录） | 扫到 1 个 ✅ |
| **注入实证**：`docs/_sens_probe.md.bak_19700101_selftest` | **FAIL，出现 1 个** ✅（探针已删） |

**归档完整性核验**（独立脚本）：17 个文件全部

- 数量 17/17；
- 非空、**mtime 保留**（`2026-09-11` ~ `2026-10-01`，说明是 rename 而非无元数据拷贝）；
- 按类型可解析（`.py` 逐个 `py_compile`、`.json` 逐个 `json.load`、`.md` 可读 utf-8）；
- 基文件仍在原位（17/17）；
- 源码树残留 0。

---

## 4. ⏸ 有意未改（需老师定方向）

### 4.1 `?` 是否纳入 SDT 的 gap 字符集

现状三套并存：

| 位置 | 字符集 |
|---|---|
| `sdt_exact.py` `_SDT_GAP`（本轮） | `-.` |
| `phylo.py:605` | `-.` |
| `phylodyn_alignqc.py:35` `_GAP_CHARS` | `-.?` |

本轮只把 `sdt_exact` 对齐到 `-.`（与 identity 页一致，这是 P0-1 的要求），
**没有**顺手把 `?` 并进来 —— 因为 `?` 在 IUPAC 里是"未知碱基"而不是"缺失"，
算不算 gap 是**判定语义**问题（算 gap → 该列整体排除；不算 → 参与比较并计错配），
会直接改变 identity 数值。`_SDT_GAP` 常量处已就地标注"如需统一到 `-.?`，须先由老师定方向"。

### 4.2 `pct_mean` 混源不可解释

`public_meta/virome.py:416` 对每样本 `pct` 求均值，而 NCBI 用「已鉴定 reads」、
GSA 用「总 reads」作分母 → 跨来源不可直接比。报告已明确其**作用域精确**：
只影响 `by_virus` 的 `pct_mean` 列，`occ` / `n_samples` 不受影响，
且 `sources` 计数已返回（`:395/419`）供前端分层。

两种改法都有代价，需老师定：

- **A**：`pct_mean` 只在单一来源内计算，混源时返回 `null`（语义干净，但会让现存数字消失）；
- **B**：保留 `pct_mean`，额外返回按来源拆分的均值（不破坏既有前端，但字段变多）。

本次**未动**（属判定语义 + API 形状变更）。

### 4.3 `kv_engines.py` 15 处裸 subprocess

报告判定为"已知设计取舍，保留观察"（管道 Popen 依赖外层 6h + 进程树兜底），本轮未动。

### 4.4 MAFFT 挂起时 `subprocess` 超时无效（§7）

`sdt_exact.py:212` 的 `sp.run(..., timeout=1800)` 对"孙进程持管道"的挂起**不生效**。
改法需要老师定方向（用什么兜底进程树、超时值多少），本轮只记录。

---

## 5. 改动清单、备份与残留清理

### 5.1 改动文件（16 改 + 2 新增）

| 类别 | 文件 |
|---|---|
| 引擎（5） | `sdt_exact.py`（P0-1 / P1-5）、`crispr_design.py`（P0-2 / P1-2）、`gbdraw_plot.py`、`genome_diag.py`、`viz.py`（P1-4） |
| Web 层（5） | `web/tasks.py`、`web/samples.py`（P1-1）、`web/tools_api.py`、`web/examples_api.py`、`web/settings.py`（P1-2） |
| 测试（4） | `_run_all.py`（P1-3）、`_check_db_paths.py`、`_check_svg_wellformed.py`、`_check_spec_datas.py` |
| 文档（2） | `wiki/Architecture.md`、`wiki/CLI-Reference.md`（P2） |
| **新增** | `tests/_it_svg_embed_chain.py`（守卫）、`docs/review_VirusPlatform_20261007_FIXES.md`（本文） |

### 5.2 备份（⚠️ 有一处不完整，如实说明）

- `archive/_bak_20261007/`：修复开工时做的快照，**13 个文件**（含 `utils.py` ——
  该文件最终**没有改**，快照留着而已）。
- ⚠️ **有 4 个文件是在快照之后才决定要改的，因此没有"改前副本"**：
  `Virus_Platform_Core/viz.py`、`web/examples_api.py`、`web/settings.py`、
  `tests/_check_spec_datas.py`。
  这 4 处的**逐字 before/after 都写在本说明的对应小节里**（§2.2 / §2.4 / §3.2），
  改动本身很小（各 1–4 处），手工回退没有歧义；`git diff <file>` 也能看到。
- 顺带说明：本仓 `git status` 里 `Virus_Platform_Core/` 的若干 ` M` 标记是**作者既有的
  未提交改动**（实测：`gbdraw_plot.py` / `_check_db_paths.py` / `wiki/Architecture.md`
  的改前快照与 `HEAD` 一致，而 `sdt_exact.py` / `web/tasks.py` / `utils.py` 等与 `HEAD`
  不同）—— 所以**不能拿 `HEAD` 当"改前状态"**，这也是上面如实标注缺副本的原因。

### 5.3 残留清理

| 项 | 处理 |
|---|---|
| `tests/._run_all.lock` 陈旧锁 | 已确认持有者 pid 已死并删除（本次误用 `--help` 触发过一次全量跑，已终止） |
| 孤儿 MAFFT/Cygwin 进程 | 排查 §7 时累积了 36 个 `mafft.bat`/`bash.exe` 孤儿进程（我杀测试留下的），已按命令行含 `mafft-win` 精确终止，**复核残留 = 0** |
| 探针临时脚本 | 全部落 `%TEMP%` 或 `run/_tmp/`，已**逐个点名**删除（不用通配符） |
| 灵敏度探针注入的文件 | `examples/_sens_probe_broken.svg`、`docs/_sens_probe.md.bak_19700101_selftest` 均已删除并复查不存在 |
| 17 个 `.bak` 残留 | 迁到 `archive/_bak_residue_20261007/`（逐文件核验：数量 17/17、非空、**mtime 保留**、`.py` 可编译 / `.json` 可解析、基文件仍在原位） |
| 未清理项 | `run/tests_crispr/`（3.3MB，`test_crispr_design.py` 自己的产物目录，在可清理的 `run/` 下，留着可省下次重跑的准备时间） |

---

## 6. 回归执行记录（2026-10-07）

**全树语法编译**：`Virus_Platform_Core/ tests/ dev_tools/ scripts/ webapp/` 共 **298 个 `.py`，0 错误**。

**本轮改动直接相关的守卫（全部通过）**

| 脚本 | 结果 |
|---|---|
| `tests/_check_sdt_gap_parity.py` | PASS |
| `tests/_check_crispr_primer_purity.py` | PASS |
| `tests/_check_task_on_finish.py` | PASS |
| `tests/_check_db_paths.py` | PASS（扫 136 文件命中 0 + 3 条自测） |
| `tests/_check_sdt_distribution.py` | PASS |
| `tests/_check_svg_wellformed.py` | PASS（39 SVG / 0 非白名单残缺 / 5 判据） |
| `tests/_it_svg_embed_chain.py` | PASS（报告内 2 片段 / 0 非良构） |
| `tests/_check_spec_datas.py` | PASS（扫 1032 文件 / 0 `.bak` 残留） |

**同批回归（全部通过）**

| 脚本 | 结果 |
|---|---|
| `pytest tests/test_crispr_design.py test_host_genome_gget.py test_mirna_alignment.py test_utils_pure.py test_datasets_cli.py test_web_guard.py` | **86 passed in 72.90s** |
| `_check_enc_utils` / `_check_i18n` / `_check_js_syntax` / `_route_inventory` | PASS |
| `_check_examples_manifest` / `_check_ascii_stage` / `_check_start_modes` / `_selftest_pyz_cmp` | PASS |
| `_check_examples_inputs`（35 项）/ `_check_genome_circos` | PASS |

**批量入口**：`tests/_run_all.py --list` = **73 项**；本轮 8 个相关守卫均已登记在清单内。

### ⚠️ 两项未能完成（环境阻塞，非本轮改动所致）

`tests/_check_identity_aa.py` 与 `tests/_check_primer.py` **本轮未跑完**，
卡在 MAFFT 调用上。详见 §7。这两项**应视为未验证**。

---

## 7. 本轮额外发现：MAFFT 在带控制台的会话里不退出，且 `subprocess` 超时对它无效

排查上面那两项时挖出来的，**与本轮任何改动无关**（未改 MAFFT 调用、未改 `_sdt_pair_worker`），
但值得单独记一笔，因为它的第二半是**真实的健壮性缺口**。

### 7.1 现象（✅ 已实测）

`3rd/tools/mafft-win/mafft.bat --quiet --auto <2 条 6.6kb 序列>`：

```
rc=124（45s 超时）  输出 226 行
对齐结果本身是完整正确的：2 条 × 6694 列，格式无误
```

即 **MAFFT 干完活、写完正确输出，然后不退出**。复现条件已逐项排除：

| 变量 | 试过 | 结果 |
|---|---|---|
| 参数 | `--quiet --auto` / `--auto --quiet --inputorder --thread 4` / `--thread 1` | 四种组合**完全一样**（45s 超时、226 行） |
| stdin | 继承 / `< /dev/null` | 一样 |
| 沙箱 | 沙箱内 / `dangerouslyDisableSandbox` | 一样 |
| 孤儿进程 | 清空前 / 清空后（36 个 → 0 个） | 一样 |
| 控制台 | 默认 / `CREATE_NO_WINDOW` / `DETACHED_PROCESS` | 一样 |

**对照组（关键）**：同一台机器、同一时刻，**平台自己的调用是正常的** ——
`run/logs/crispr_design_20261007_154358.log` 16:04:17 记录：

```
[CMD] $ ...\mafft.bat --auto --quiet --inputorder --thread 4 ...\in.fa > aln.fa
[OUT] [退出码 0] 耗时 3.7s
```

⇒ **不是 MAFFT 坏了、也不是参数问题，而是"被谁调用"**。

**机制（推断，有代码依据）**：`mafft.bat` 是 Cygwin/MSYS 包装（`usr/bin/bash /usr/bin/mafft`），
`usr/bin/mafft` 脚本尾部有一段交互式 AWK：

```awk
printf( "Press Enter to exit." ) > "/dev/tty";
res = getline < "/dev/tty";
```

**带控制终端的会话**（本轮的 Git Bash pty）→ `/dev/tty` 可用 → 等按键，永不返回；
**平台的无终端调用**（GUI/子进程）→ `/dev/tty` 打不开 → `getline` 失败、`NF==0` 按"直接回车"处理 → 正常退出。
这也解释了为什么同一命令在 GUI 里 3.7s、在我这个会话里永久挂起。

### 7.2 真正的缺口：`subprocess.run(timeout=)` 在这种挂起下**不生效**（✅ 已实测 + 代码依据）

`Virus_Platform_Core/sdt_exact.py:212`（`_sdt_pair_worker`）：

```python
r = sp.run([_W['mafft'], '--quiet', '--auto', in_f],
           capture_output=True, text=True, timeout=1800,
           creationflags=getattr(sp, 'CREATE_NO_WINDOW', 0))
```

`faulthandler` 抓到的现场（4 个池工作进程全卡在同一处）：

```
File "...\subprocess.py", line 1628 in _communicate
File "...\subprocess.py", line 1209 in communicate
File "...\subprocess.py", line 550 in run
File "...\Virus_Platform_Core\sdt_exact.py", line 212 in _sdt_pair_worker
```

实测：`sp.run(timeout=40)` 在挂起输入上 **40s 后没有抛 `TimeoutExpired`**
（外层 `timeout 150` 才把它杀掉）。

**原因**（`subprocess.py` 实现）：超时分支是

```python
except TimeoutExpired as exc:
    process.kill()                      # 只杀直接子进程
    ...
    stdout, stderr = process.communicate()   # ← 卡在这里，异常还没抛出去
    raise
```

这里直接子进程是 `cmd.exe`（跑 `mafft.bat`），真正干活的是**孙进程** Cygwin `bash`/`mafft.exe`，
它们仍持有 stdout/stderr 管道句柄 → 第二次 `communicate()` 永远等不到 EOF →
**`TimeoutExpired` 永远不抛出**，`timeout=1800` 形同虚设。

**影响**：一旦 MAFFT 走到这个状态，SDT / identity / primer / CRISPR 等所有依赖 MAFFT 的
任务都不会被自己的超时救回来，只能靠外层进程树终止兜底。

**⏸ 本轮未改**：改法涉及"用什么手段兜底进程树 + 超时值"（`utils._kill_proc_tree` 已存在，
但 `_sdt_pair_worker` 是池工作进程，不能直接用），属需要老师定方向的改动，故只记录不动手。

**给老师的可复跑命令**（在**普通控制台**里跑应能通过；在带终端的自动化会话里会挂）：

```bat
python tests/_check_identity_aa.py
python tests/_check_primer.py
```

---

*本说明与审查报告 `review_VirusPlatform_20261007.md` 配套。*
*本轮改动 16 个文件（引擎 5 / Web 5 / 测试 4 / 文档 2）+ 新增 2 个；备份见 §5.2。*
