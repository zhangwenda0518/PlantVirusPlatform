# -*- coding: utf-8 -*-
r"""TreeTime —— 最大似然分子钟定年 / 离散性状地理重构 / 天际线 / 祖先序列。

## 定位

这是「本地时间与地理推断」的 **ML 引擎**（轨道 A）。它与平台另外两条定年通路
**并列而非替代**：

| 通路 | 钟模型 | 地理 | 群体动态 | 成本（476 叶实测） |
|---|---|---|---|---|
| `phylogeo.lsd2_dating` | 严格钟（最小二乘） | — | — | 秒级 |
| `treedater.py` | 非相关松弛钟 | — | LTT | 定年 12 s；parboot 约 8 min |
| **本模块** | ML 分子钟（`--relax` 可切松弛） | **mugration** | **skyline** | 定年 61.7 s + 地理 5.3 s |

⚠️ **三条通路的速率不可互相替代、不可混用**：实测同一份 476 叶 H3N2 数据，
treedater 3.7e-3 vs TreeTime 3.12e-3（差 7–9%），LSD2 2.10e-3。差异来自钟模型，
不是谁算错了。任何引用都必须带方法名。

## 五个子命令里平台用到四个

`clock`（定年 + RTT + 可选天际线）、`mugration`（地理迁移）、`ancestral`（祖先序列）、
`homoplasy`（替换饱和）。`arg` 未接（重组图，平台已有 RDP5）。

## ⚠️⚠️ 四个本机实测出来的坑（改代码前务必先读）

1. **`clock` 子命令是「缩水接口」，定年一律走「顶层无子命令」形式。**
   实测（TreeTime 0.12.1）：`treetime clock --confidence` 与
   `treetime clock --time-marginal only-final` 都直接报
   `error: unrecognized arguments`，`--coalescent` 同样不认。
   子命令只认 `--tree/--aln/--dates/--reroot/--keep-root/--clock-filter/
   --covariation/--plot-rtt/--allow-negative-rate` 那一小撮。
   而顶层解析器（不带子命令）才是老版完整接口，`--confidence`
   `--time-marginal` `--coalescent` `--n-skyline` `--relax` `--gen-per-year`
   `--max-iter` 全在它上面。**本模块 `run_clock` 一律用顶层形式**，
   调用方不用关心。
2. **dates / states CSV 的表头必须含 `name` / `strain` / `accession`。**
   否则 `treetime.utils.parse_dates` 抛
   `MissingDataError: need at least one column that contains the taxon labels`。
   平台内部一律用 `name`（并显式传 `--name-column name`），不吃这个亏。
3. **NEXUS 输出的注释里有逗号**，平台现成的 `phylogeo.parse_newick` 会在注释里断句。
   形如 `SWM_AF220109_1997.00000:106.378[&mutations="C131T,T367C,A384T"]` ——
   `parse_newick` 的 `_name_len()` 读到第一个 `,` 就收手，随后把 `T367C` 当成叶名。
   本模块自带 `parse_annotated_newick()`（状态机 + `[&…]` 整体吞掉），
   产出与 `phylogeo` 节点**同构**的 dict（`name`/`length`/`children` + 额外 `attrs`），
   因此 `phylogeo._iter_pre` / `tree_root_distances` / `time_tree_segments` 可直接复用。
4. **`confidence.csv` 里的状态是字母 `A`/`B`/`C`**，不是区划名。
   必须用 `GTR.txt` 的 `Character to attribute mapping:` 段翻译。
   拿不到 mapping 时**不许猜**（猜错会让整张走廊表静默错位），本模块返回空并告警。

## ⚠️⚠️⚠️ 第五个坑（最要命的一个）：**喂哪棵树，结果差一倍**

实测 RSV 209（时间信号 r²=0.04，接近噪声）同一份数据、同一套 states：

| 喂给 mugration 的树 | 迁移次数 | mu |
|---|---|---|
| `timetree.nwk`（枝长＝年） | **70** | 0.0255 /年 |
| `divergence_tree.nexus`（枝长＝替换/位点） | **139** | 79.9 |

拓扑完全相同（354 条边一一对应），差的是**枝长**。查清了原因：
`treetime clock` 输出的 `divergence_tree.nexus` **不是** `rate × 时间` —— 实测
div/time 比例中位数 3.12e-4，但离散度 `max/median = 201`，且零长枝 65 条
（时间树只有 5 条）。也就是说 divergence 树的枝长基本就是**输入树的原样**，
而时间树的枝长是严格钟下重新算的年代差（自洽性自检：0 条边不自洽）。

**时间信号强时两棵树近似成比例，两种喂法结果一致；时间信号弱时不成比例，
结果就会分叉。** 所以这不是"哪个对"的问题，而是"必须对拍并如实报出分歧"。

本模块的做法（也是平台的口径纪律）：
1. **默认喂时间树** —— 与四仓（phymap / Mugration-Analysis）一致，且只有时间树
   有定义好的时间轴，`events_year.tsv` 才算得出真年份；`mu` 的量纲也干净（/年）。
2. **`also_tree=` 自动对拍**：给了第二棵树就再跑一遍，把两边的迁移次数一起报出来；
   分歧 > 20% 时**显式告警**并说明"时间信号弱 → 地理重建对枝长口径敏感"。
   卡片上这两个数字并排显示，不许只留一个。

## 口径纪律（写死在返回值里）

- **R² 低 ≠ 没时间信号**（本模块实测更正过一次）：RSV 209 固定根位 R²=0.049，
  最大根位 R²=0.112，DRT 200 次置换 **p=0.005 通过**。判定看 DRT，不看 R²。
- **mugration 的 `confidence.csv` 是 ML 边际概率，不是贝叶斯后验。** 措辞不得简化；
  它给的是「在该树 + 该 GTR 下的似然边际」，没有先验、没有拓扑不确定性。
- **`skyline.tsv` 的上下界是 ±2 SD 的近似区间**（TreeTime 自己在表头写明
  `approximate confidence bounds (+/- 2.000000 standard deviations of the LH)`），
  **不是** BSP 的贝叶斯 95% HPD。且 `N_e` 的绝对标度取决于 `--gen-per-year`
  （默认 50），换假设就换数字，引用时必须一起写。
- **`mu` 的单位取决于喂进去的树**：喂替换/位点树 → 每位点每次换状态的概率；
  喂时间树（枝长＝年）→ **每年**换状态的概率。重建结果对枝长整体缩放不变，
  所以两种喂法的**走廊计数相同**，只有 `mu` 的量纲不同。

## 环境

TreeTime 装在平台首选解释器（`C:\Python312`）的 Scripts 下，与 `gbdraw` / `orfipy`
同一套约定。传递依赖 numpy / scipy / pandas / biopython **全部已在
`requirements.txt`**，零新增依赖。打包分发时 `treetime.exe` 属外部工具，
按 `platform.json` 的 `tools.treetime` 指定路径；未登记时 `find_treetime()` 会
按「环境变量 → platform.json → 3rd/ → 首选解释器同级 Scripts → PATH」逐级找。
"""

import io
import os
import re
import subprocess
import sys

from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.utils import check_path, task_check_cancel, \
    task_register_proc

# TreeTime 的五个子命令（0.12.1 实测；`arg` 未接）
TT_SUBCMDS = ('homoplasy', 'ancestral', 'mugration', 'clock', 'arg')

# TreeTime 会给「缺数据」补一个伪状态，在 GTR.txt 的 mapping 里显示为 `?`
# （实测 RSV 209：`F: ?`）。它不是区划，**绝不能**混进迁移矩阵 ——
# 否则会凭空多出一条「某区划 ↔ ?」的走廊，且计数被摊薄。
MISSING_STATES = ('?', 'Unknown', 'NA', 'N/A', '', 'nan')

# Python 侧 locale 净化：与 treedater.clean_env 同一理由（子进程别继承
# bash/MSYS 导出的 C.UTF-8，Windows 上的外部程序未必认）
_BAD_LOCALE_VARS = ('LC_ALL', 'LANG', 'LC_CTYPE', 'LC_COLLATE', 'LC_TIME',
                    'LC_MONETARY', 'LC_NUMERIC', 'LANGUAGE')

# 首选解释器同级 Scripts（platform.json 里 gbdraw / orfipy 也是这一套）
_PREFERRED_PY = r'C:\Python312'


def clean_env(extra=None):
    """剥掉外部程序不认的 locale 变量 + 强制 UTF-8 输出（公开便于测试复用）。"""
    env = dict(os.environ)
    for k in _BAD_LOCALE_VARS:
        env.pop(k, None)
    env['PYTHONIOENCODING'] = 'utf-8'
    env['PYTHONUTF8'] = '1'
    if extra:
        env.update(extra)
    return env


def find_treetime():
    """找 treetime 可执行入口。返回 (argv 前缀 list, 说明串)。

    逐级：环境变量 `TREETIME_EXE` → platform.json 的 `tools.treetime`
    → `3rd/tools/treetime/treetime.exe` → 首选解释器 Scripts → 当前解释器
    Scripts → PATH 上的 `treetime`。都找不到返回 `(None, 原因)` ——
    **不抛异常**，由调用方决定怎么如实报错（卡片要能显示"没装"而不是崩）。
    """
    tried = []

    def _ok(p, how):
        return ([p], f'{how}：{p}') if os.path.isfile(p) else None

    p = (os.environ.get('TREETIME_EXE') or '').strip()
    if p:
        tried.append(f'环境变量 TREETIME_EXE={p}')
        r = _ok(p, '环境变量')
        if r:
            return r

    try:
        from Virus_Platform_Core.web.state import cfg
        v = (cfg.tools or {}).get('treetime')
    except Exception:
        v = None
    if v and not str(v).startswith('bundled:'):
        cand = v if os.path.isabs(str(v)) else os.path.join(PLATFORM_ROOT, str(v))
        tried.append(f'platform.json tools.treetime={v}')
        r = _ok(cand, 'platform.json')
        if r:
            return r

    for cand in (os.path.join(PLATFORM_ROOT, '3rd', 'tools', 'treetime',
                              'treetime.exe'),
                 os.path.join(_PREFERRED_PY, 'Scripts', 'treetime.exe')):
        tried.append(cand)
        r = _ok(cand, '内置/首选解释器')
        if r:
            return r

    # 当前解释器同级 Scripts（服务可能就跑在带 treetime 的那个解释器上）
    cand = os.path.join(os.path.dirname(sys.executable), 'Scripts', 'treetime.exe')
    tried.append(cand)
    r = _ok(cand, '当前解释器')
    if r:
        return r

    # 最后：解释器本体 + `-m treetime`（不依赖 .exe 包装器）
    for py in (os.path.join(_PREFERRED_PY, 'python.exe'), sys.executable):
        if not py or not os.path.isfile(py):
            continue
        if _has_module(py, 'treetime'):
            return ([py, '-m', 'treetime'], f'python -m treetime：{py}')
        tried.append(f'{py} -m treetime（模块不可导入）')
    return (None, '找不到 TreeTime；试过：' + '；'.join(tried))


def _has_module(python_exe, mod, timeout=60):
    """探测某解释器能否 import 指定模块（TreeTime 的 .exe 是 console_script，
    发行版里可能只装了库没装脚本）。"""
    try:
        p = subprocess.run([python_exe, '-c', f'import {mod}'],
                           capture_output=True, env=clean_env(),
                           timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return p.returncode == 0


def check_tt_env(argv=None, timeout=120):
    """体检 TreeTime 环境：版本 + 可用子命令 + 首选解释器依赖。

    只跑两次轻量子进程（`version` 与顶层 `--help`），不跑真实分析 ——
    体检必须秒回，否则卡片一按"运行"就先卡住。
    """
    out = {'ok': False, 'exe': None, 'how': None, 'version': None,
           'subcommands': [], 'missing_subcommands': [], 'python': None,
           'problems': [], 'raw': ''}
    if argv is None:
        argv, how = find_treetime()
        out['how'] = how
    else:
        how = '（显式指定）'
    if not argv:
        out['problems'].append(how)
        return out
    out['exe'] = ' '.join(argv)

    rc, lines = _run_tt(argv + ['version'], timeout=timeout)
    out['raw'] = '\n'.join(lines[:20])
    m = re.search(r'treetime\s+([0-9][0-9.]*)', out['raw'], re.I)
    if m:
        out['version'] = m.group(1)
    elif rc != 0:
        out['problems'].append(
            f'`treetime version` 退出码 {rc}：{out["raw"][:200]}')

    # 子命令清单从顶层 --help 的 positional 行取；取不到就退到内置常量
    _rc2, h = _run_tt(argv + ['--help'], timeout=timeout)
    txt = '\n'.join(h)
    m = re.search(r'\{([a-z,]+)\}', txt)
    subs = [s.strip() for s in m.group(1).split(',')] if m else list(TT_SUBCMDS)
    out['subcommands'] = subs
    out['missing_subcommands'] = [s for s in ('clock', 'mugration')
                                  if s not in subs]
    if out['missing_subcommands']:
        out['problems'].append(
            'TreeTime 缺少必需子命令：' + '、'.join(out['missing_subcommands'])
            + '（版本过老？本模块按 0.12.x 的 CLI 契约写的）')

    # 依赖：TreeTime 本体在哪个解释器里，就查哪个
    py = None
    if len(argv) >= 3 and argv[1] == '-m':
        py = argv[0]
    elif argv[0].lower().endswith('.exe'):
        guess = os.path.join(os.path.dirname(os.path.dirname(argv[0])),
                             'python.exe')
        py = guess if os.path.isfile(guess) else None
    out['python'] = py
    if py:
        miss = [m for m in ('numpy', 'scipy', 'pandas', 'Bio')
                if not _has_module(py, m, timeout=90)]
        out['missing_python_deps'] = miss
        if miss:
            out['problems'].append(
                'TreeTime 所在解释器缺依赖：' + '、'.join(miss)
                + f'（{py}；按 requirements.txt 装齐即可）')
    out['ok'] = not out['problems']
    return out


def _run_tt(cmd, timeout=600, cwd=None, say=None, capture=True):
    """跑一条 treetime 命令：逐行流式写日志 + 同时收进内存供解析。

    与 `utils.run_cmd` 的区别：那个把 stdout/stderr 合并后只给日志，
    解析不了 `homoplasy` 这种「结果只打印到控制台」的子命令。
    返回 `(returncode, lines)`；**不抛异常**（超时/启动失败折成 rc=-1/-2），
    由调用方把原因写进 summary —— 卡片要能如实显示失败，不能 500。
    """
    cmd = [str(c) for c in cmd]
    if say:
        say('$ ' + subprocess.list2cmdline(cmd))
    try:
        task_check_cancel()
    except Exception:
        pass
    creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0) \
        if sys.platform == 'win32' else 0
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=False,
                                cwd=cwd, env=clean_env(),
                                creationflags=creationflags)
    except OSError as e:
        msg = f'[启动失败] {e.__class__.__name__}: {e}'
        if say:
            say(msg)
        return -2, [msg]
    task_register_proc(proc)
    lines = []
    import threading
    _to = {'hit': False}

    def _watchdog():
        try:
            proc.wait(timeout=timeout)
        except Exception:
            _to['hit'] = True
            try:
                proc.kill()
            except Exception:
                pass

    th = threading.Thread(target=_watchdog, daemon=True)
    th.start()
    try:
        for raw in proc.stdout:
            line = _decode(raw).rstrip()
            lines.append(line)
            if say and line:
                say(line)
    except Exception:
        pass
    th.join(timeout=5)
    rc = proc.returncode if proc.returncode is not None else -1
    if _to['hit']:
        if say:
            say(f'[超时] 命令运行超过 {timeout}s 被终止')
        return -1, lines + [f'[超时] >{timeout}s']
    if say:
        say(f'[退出码 {rc}]')
    return rc, lines


def _decode(raw):
    """字节 → 文本（TreeTime 输出一般 ASCII，中文路径可能按 GBK 打）。"""
    if isinstance(raw, str):
        return raw
    for enc in ('utf-8', 'gbk', 'latin-1'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', errors='replace')


# ---------------------------------------------------------------- 输入表

def write_dates_csv(pairs, path):
    """写 TreeTime 认的日期表：表头 **必须**含 `name`（坑 2）。

    pairs: `[(name, date_float_or_iso), …]`。date 直接给数值（十进制年）
    或 `YYYY-MM-DD` 字符串，TreeTime 两种都收。
    """
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        f.write('name,date\n')
        for n, d in pairs:
            f.write(f'{_csv_cell(n)},{_csv_cell(d)}\n')
    return path


def write_states_csv(pairs, path, attribute='region'):
    """写离散性状表：`name,<attribute>`（mugration 的 `--states` 输入）。"""
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        f.write(f'name,{attribute}\n')
        for n, s in pairs:
            f.write(f'{_csv_cell(n)},{_csv_cell(s)}\n')
    return path


def write_weights_csv(pairs, path, attribute='region'):
    """写均衡频率权重表（`--weights`）：`<attribute>,weight`。"""
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        f.write(f'{attribute},weight\n')
        for s, w in pairs:
            f.write(f'{_csv_cell(s)},{w}\n')
    return path


def _csv_cell(v):
    s = '' if v is None else str(v)
    if any(c in s for c in ',"\n'):
        return '"' + s.replace('"', '""') + '"'
    return s


# ---------------------------------------------------------------- NEXUS / 注释树

_ATTR_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s*=\s*"([^"]*)"')
_BRACKET_RE = re.compile(r'\[&([^\]]*)\]')


def nexus_tree_newick(text):
    """从 NEXUS 文本里抠出 `Tree tree1=…;` 那一行的 newick（保留注释）。

    抠不到抛 ValueError —— 静默返回空串会让下游把"没解析到"当成"树是空的"。
    """
    for ln in str(text).splitlines():
        s = ln.strip()
        if s.lower().startswith('tree ') and '=' in s:
            return s.split('=', 1)[1].strip().rstrip(';').strip()
    raise ValueError('NEXUS 里找不到 `Tree …=…;` 行（文件可能不是 TreeTime 的 NEXUS 输出）')


def parse_annotated_newick(text):
    """带 `[&k="v"]` 注释的 Newick 解析 → 节点树。

    节点 dict：`{'name','length','children','attrs'}` —— 前三个键与
    `phylogeo.parse_newick` **完全同构**，所以 `phylogeo._iter_pre` /
    `tree_root_distances` / `time_tree_segments` / `fitch_mugration` 可直接吃。

    ⚠️ 为什么不能用 `phylogeo.parse_newick`：`[&mutations="A,B,C"]` 里**有逗号**，
    那边的 `_name_len()` 见到 `,` 就收手，把注释后半截当成叶名（坑 3）。
    这里的状态机把 `[…]` 整体当"注释 token"吞掉，绝不进名字/枝长缓冲。
    """
    s = str(text).strip()
    if s.startswith('#NEXUS'):
        s = nexus_tree_newick(s)
    s = s.strip().rstrip(';')
    pos = [0]
    head = s[:60].replace('\n', ' ')

    def _attrs(tok):
        out = {}
        for blk in _BRACKET_RE.findall(tok):
            for k, v in _ATTR_RE.findall(blk):
                out[k] = v
            # 无 key 的纯数字注释（如 [&0.95]）留给调用方按需处理
        return out

    def _skip_comments():
        """跳过当前位置的 `[…]` 注释块（含嵌套）。返回拼起来的原文。"""
        buf = ''
        while pos[0] < len(s) and s[pos[0]] == '[':
            depth = 0
            start = pos[0]
            while pos[0] < len(s):
                c = s[pos[0]]
                if c == '[':
                    depth += 1
                elif c == ']':
                    depth -= 1
                    if depth == 0:
                        pos[0] += 1
                        break
                pos[0] += 1
            buf += s[start:pos[0]]
        return buf

    def _label():
        """读「名字 [注释] [:长度 [注释]]」，注释归到该节点。"""
        buf = ''
        while pos[0] < len(s) and s[pos[0]] not in ',():[]':
            buf += s[pos[0]]
            pos[0] += 1
        attrs = _attrs(_skip_comments())
        length = 0.0
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            lb = ''
            while pos[0] < len(s) and s[pos[0]] not in ',()[]':
                lb += s[pos[0]]
                pos[0] += 1
            attrs.update(_attrs(_skip_comments()))
            try:
                length = float(lb)
            except ValueError:
                length = 0.0
        return buf.strip(), length, attrs

    def _node():
        node = {'name': '', 'length': 0.0, 'children': [], 'attrs': {}}
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1
            while True:
                node['children'].append(_node())
                if pos[0] < len(s) and s[pos[0]] == ',':
                    pos[0] += 1
                    continue
                if pos[0] < len(s) and s[pos[0]] == ')':
                    pos[0] += 1
                    break
                raise ValueError(
                    'Newick 括号不配对：位置 %d 读到 %s，预期 "," 或 ")"。树开头：%s'
                    % (pos[0], '结尾' if pos[0] >= len(s) else repr(s[pos[0]]), head))
            node['name'], node['length'], node['attrs'] = _label()
        else:
            node['name'], node['length'], node['attrs'] = _label()
        return node

    root = _node()
    rest = s[pos[0]:].strip()
    if rest:
        raise ValueError('Newick 解析未走完：位置 %d 剩 %r（括号不配对？）。树开头：%s'
                         % (pos[0], rest[:40], head))
    return root


def iter_pre(root):
    """先序遍历（本地实现，避免为一个循环去 import phylogeo 的全部依赖）。"""
    yield root
    for c in root['children']:
        yield from iter_pre(c)


def tree_distances(root):
    """每节点到根的距离（枝长累加）→ {id(node): float}。"""
    out = {}

    def _rec(n, d):
        out[id(n)] = d
        for c in n['children']:
            try:
                ln = float(c.get('length') or 0.0)
            except (TypeError, ValueError):
                ln = 0.0
            _rec(c, d + ln)

    _rec(root, 0.0)
    return out


def to_newick(node, with_attrs=True):
    """节点树 → newick 文本（`parse_annotated_newick` 的逆向）。

    ⚠️ 枝长用 `%.10g` 而不是 `%.6f`：预比对树的枝长常见 1e-07 量级，
    固定 6 位小数会把它**四舍五入成 0**（等于把树抹平），10 位有效数字既保
    量级又不丢精度。

    ⚠️ `[&k=v]` 注释回写在**名字之后、冒号之前** —— TreeTime 只认这个位置
    （实测：写到枝长后面它读不到，节点属性全空）。解析时并进 `attrs` 的
    "枝长后注释"因此会挪到前面，键值不丢。
    """
    name = node.get('name') or ''
    lbl = ''
    if with_attrs and (node.get('attrs') or {}):
        lbl = '[&' + ','.join('%s=%s' % (k, v)
                             for k, v in node['attrs'].items()) + ']'
    ln = ''
    if node.get('length') is not None:
        try:
            ln = ':' + ('%.10g' % float(node['length']))
        except (TypeError, ValueError):
            ln = ''
    kids = node.get('children') or []
    if not kids:
        return '%s%s%s' % (name, lbl, ln)
    return '(%s)%s%s%s' % (','.join(to_newick(c, with_attrs) for c in kids),
                           name, lbl, ln)


def prune_tree(root, keep):
    """按叶名集合剪树 → `(新树, 统计)`。不在 `keep` 里的叶连同只服务于它们的
    内部节点一起去掉。

    `统计`：`removed_tips`（被删掉的叶名，调用方要能逐条说清剔了谁）/
    `n_removed` / `n_collapsed` / `n_kept_tips`。

    ⚠️ 只剩一个孩子的内部节点**必须压缩**（返回那个孩子、枝长相加）——
    留着 `(A:0.1):0.2` 这种单子节点的 newick 语法合法，但下游 TreeTime /
    LSD2 / ape 读进去会各自给出不同解释。根不压缩（根要保持是一棵树）。
    """
    keep = set(keep)
    removed = []
    collapsed = [0]

    def _rec(n, is_root=False):
        kids = n.get('children') or []
        attrs = dict(n.get('attrs') or {})
        if not kids:
            nm = n.get('name') or ''
            if nm in keep:
                return {'name': nm, 'length': n.get('length'),
                        'children': [], 'attrs': attrs}
            removed.append(nm)
            return None
        new_kids = []
        for c in kids:
            r = _rec(c)
            if r is not None:
                new_kids.append(r)
        if not new_kids:
            return None
        if len(new_kids) == 1 and not is_root:
            child = new_kids[0]
            try:
                child['length'] = (float(child.get('length') or 0.0)
                                   + float(n.get('length') or 0.0))
            except (TypeError, ValueError):
                pass
            collapsed[0] += 1
            return child
        return {'name': n.get('name') or '', 'length': n.get('length'),
                'children': new_kids, 'attrs': attrs}

    new_root = _rec(root, is_root=True)
    n_kept = 0 if new_root is None else len(
        [x for x in iter_pre(new_root) if not x['children']])
    return new_root, {'removed_tips': removed, 'n_removed': len(removed),
                      'n_collapsed': collapsed[0], 'n_kept_tips': n_kept}


def node_dates_from_branches(root, tip_dates):
    """时间树（枝长＝年）→ 每个节点的**日历年**。

    原理：TreeTime 的时间树把叶采样时刻钉死，所以
    `date_root = date_tip − dist(root→tip)` 对每个叶都成立（差异来自数值误差，
    取中位数最稳）。于是 `date_node = date_root + dist(root→node)`。

    ⚠️ 不用 `annotated_tree.nexus` 里的 `date=` 注释：mugration 重写标注树时
    **把 date 换成了区划**（实测），而且两边的 `NODE_xxxxxxx` 编号也对不上。
    从枝长算既不依赖命名也不依赖注释，最稳。

    返回 `({id(node): year}, root_year)`；叶日期缺失过多（< 半数）时 root_year
    为 None，调用方须据此拒绝输出年份轴（不许拿相对量冒充日历年）。
    """
    dist = tree_distances(root)
    tips = [n for n in iter_pre(root) if not n['children']]
    ests = []
    for t in tips:
        d = tip_dates.get(t['name'])
        if d is None:
            continue
        try:
            ests.append(float(d) - dist[id(t)])
        except (TypeError, ValueError):
            continue
    if not ests or len(ests) < max(3, len(tips) // 2):
        return {}, None
    ests.sort()
    root_year = ests[len(ests) // 2]
    return ({id(n): root_year + dist[id(n)] for n in iter_pre(root)}, root_year)


# ---------------------------------------------------------------- 结果解析

def parse_clock_txt(path):
    """读 `molecular_clock.txt` → {'rate','r2'}（TreeTime 的根到尾回归）。"""
    out = {}
    if not os.path.isfile(path):
        return out
    with io.open(path, encoding='utf-8', errors='replace') as f:
        txt = f.read()
    m = re.search(r'--rate:\s*([0-9.eE+-]+)', txt)
    if m:
        try:
            out['rate'] = float(m.group(1))
        except ValueError:
            pass
    m = re.search(r'--r\^2:\s*([0-9.eE+-]+)', txt)
    if m:
        try:
            out['r2'] = float(m.group(1))
        except ValueError:
            pass
    return out


def parse_gtr_txt(path):
    """读 `GTR.txt` → mapping / mu / pi / W_ij / Q_ij。

    返回 dict：
      `mapping` {字母: 状态名}、`states` [状态名]、`mu`、
      `pi` {状态: 频率}、`W` {(i,j): 对称速率}、`Q` {(i,j): 方向速率}
    """
    out = {'mapping': {}, 'states': [], 'mu': None, 'pi': {}, 'W': {}, 'Q': {}}
    if not os.path.isfile(path):
        return out
    with io.open(path, encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()

    mode = None
    hdr = []
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith('character to attribute mapping'):
            mode = 'map'
            continue
        if low.startswith('substitution rate'):
            m = re.search(r'([0-9.eE+-]+)\s*$', s)
            if m:
                try:
                    out['mu'] = float(m.group(1))
                except ValueError:
                    pass
            mode = None
            continue
        if low.startswith('equilibrium frequencies'):
            mode = 'pi'
            continue
        if low.startswith('symmetrized rates'):
            mode = 'W'
            hdr = []
            continue
        if low.startswith('actual rates'):
            mode = 'Q'
            hdr = []
            continue

        if mode == 'map':
            # ⚠️ 代号**不只有字母数字**：状态多于 36 个时 TreeTime 会用到
            #    `[` `\` `]` `^` `_` `` ` `` 等标点（实测 `]`=Russia、`[`=Peru）。
            #    原来写 `[A-Za-z0-9]+` → 这些行读不到 → 映射缺条 → 走廊表里
            #    真地名与未翻译代号混用、n_changes 虚高。
            m = re.match(r'^(\S)\s*:\s*(.*)$', s)
            if m:
                out['mapping'][m.group(1)] = m.group(2).strip()
        elif mode == 'pi':
            # 同 `map`：代号可能是标点（否则均衡频率表也缺条，采样偏倚校正会算歪）
            m = re.match(r'^(\S)\s*:\s*([0-9.eE+-]+)$', s)
            if m:
                try:
                    out['pi'][m.group(1)] = float(m.group(2))
                except ValueError:
                    pass
        elif mode in ('W', 'Q'):
            cells = s.split()
            if not hdr and len(cells) >= 2 and not _isnum(cells[0]):
                hdr = cells
                continue
            if len(cells) >= 2 and _isnum(cells[0]):
                # 首列是行标签（字母），随后是数值
                row = cells[0]
                for j, v in enumerate(cells[1:]):
                    try:
                        fv = float(v)
                    except ValueError:
                        continue
                    out[mode][(row, j)] = fv
    # 状态字母顺序：mapping 的插入序就是 TreeTime 的列序
    out['states'] = list(out['mapping'].keys())
    return out


def _isnum(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def parse_confidence_csv(path, mapping=None):
    """读 mugration 的 `confidence.csv` → (rows, summary)。

    rows: `[{'name','probs':{状态: p},'best','best_p','is_tip'}, …]`
    表头形如 `#name, A, B, C, D, E`（**字母**，坑 4）；给了 mapping 就翻译成区划名。

    summary 里放**必须一起看**的分布：`n_nodes` / `n_internal` /
    `median_max_p_internal` / `n_below_0_7` / `frac_below_0_7` / `min_max_p`。
    只报"成功率"而不报置信分布，是这类 ML 重构最常见的误用。
    """
    rows, summary = [], {}
    if not os.path.isfile(path):
        return rows, summary
    with io.open(path, encoding='utf-8', errors='replace') as f:
        lines = [ln.rstrip('\n').rstrip('\r') for ln in f if ln.strip()]
    if not lines:
        return rows, summary
    hdr = [c.strip().lstrip('#').strip() for c in lines[0].split(',')]
    if not hdr:
        return rows, summary
    keys = hdr[1:]
    labels = [(mapping or {}).get(k, k) for k in keys]
    summary['confidence_labels_raw'] = keys
    summary['confidence_labels'] = labels
    summary['confidence_labels_mapped'] = bool(mapping)
    # **部分**未映射也要能报数（原来只记"有没有 mapping"，全映射/半映射分不出来）
    summary['n_labels_unmapped'] = sum(1 for k in keys if k not in (mapping or {}))
    summary['confidence_labels_unmapped'] = [k for k in keys
                                             if k not in (mapping or {})]

    for ln in lines[1:]:
        cells = [c.strip() for c in ln.split(',')]
        if len(cells) < 2:
            continue
        nm = cells[0]
        probs = {}
        for k, v in zip(labels, cells[1:]):
            try:
                probs[k] = float(v)
            except ValueError:
                continue
        if not probs:
            continue
        best = max(probs, key=probs.get)
        rows.append({'name': nm, 'probs': probs, 'best': best,
                     'best_p': probs[best],
                     'is_tip': not str(nm).startswith('NODE_')})

    ints = [r['best_p'] for r in rows if not r['is_tip']]
    if ints:
        ints_sorted = sorted(ints)
        summary.update({
            'n_nodes': len(rows), 'n_internal': len(ints),
            'median_max_p_internal': ints_sorted[len(ints_sorted) // 2],
            'min_max_p': ints_sorted[0],
            'n_below_0_7': sum(1 for x in ints if x < 0.7),
            'n_below_0_5': sum(1 for x in ints if x < 0.5),
            'frac_below_0_7': round(sum(1 for x in ints if x < 0.7) / len(ints), 4),
        })
    return rows, summary


def parse_skyline_tsv(path):
    """读 `skyline.tsv` → (rows, meta)。

    表头注释里带 `assuming 50.0 gen/year` 与 `+/- 2.0 standard deviations`——
    这两个数字是**口径的一部分**，必须原样带出去，不能只留 `Ne` 列。
    rows: `[{'date','ne','lower','upper'}, …]`
    """
    rows, meta = [], {}
    if not os.path.isfile(path):
        return rows, meta
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s.startswith('#'):
                m = re.search(r'([0-9.]+)\s*gen/year', s)
                if m:
                    meta['gen_per_year'] = float(m.group(1))
                m = re.search(r'\+/-\s*([0-9.]+)\s*standard deviations', s)
                if m:
                    meta['sd_multiple'] = float(m.group(1))
                continue
            cells = s.split('\t') if '\t' in s else s.split()
            if len(cells) < 4:
                continue
            try:
                rows.append({'date': float(cells[0]), 'ne': float(cells[1]),
                             'lower': float(cells[2]), 'upper': float(cells[3])})
            except ValueError:
                continue
    return rows, meta


def parse_branch_mutations(path):
    """读 ancestral 的 `branch_mutations.txt` → rows（node/state1/pos/state2）。"""
    rows = []
    if not os.path.isfile(path):
        return rows
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for i, ln in enumerate(f):
            s = ln.rstrip('\n').rstrip('\r')
            if not s.strip():
                continue
            cells = s.split('\t') if '\t' in s else s.split()
            if i == 0 and cells and cells[0].lower() == 'node':
                continue
            if len(cells) >= 4:
                rows.append({'node': cells[0], 'from': cells[1],
                             'pos': cells[2], 'to': cells[3]})
    return rows


def parse_homoplasy_stdout(lines):
    """从 `treetime homoplasy` 的 stdout 里抠出同塑性表。

    ⚠️ 该子命令**不落文件**（实测 outdir 为空），结果只在控制台。
    表格形如 `\tC6T\t17`。解析失败返回空表 —— 不编数字。
    """
    rows, in_tab = [], False
    for ln in lines or []:
        s = ln.strip()
        if not s:
            continue
        if s.lower().startswith('mut') and 'multiplicity' in s.lower():
            in_tab = True
            continue
        if in_tab:
            cells = s.split()
            if len(cells) >= 2 and not cells[0].lower().startswith('mut'):
                try:
                    rows.append({'mut': cells[0], 'multiplicity': int(cells[1])})
                except ValueError:
                    pass
    return rows


# ---------------------------------------------------------------- 走廊计数

def count_corridors(root, state_of, node_year=None):
    """按分支统计迁移走廊（方向敏感）。

    `state_of`: {节点名: 状态}（mugration 用 `confidence.csv` 的 argmax）。
    `node_year`: {id(node): 日历年}（`node_dates_from_branches` 的产物）；
    给了就顺带把每条迁移的**发生年份**带上（借 Mugration-Analysis 的
    `Year, Origin, Destination` 三列），没有年份的列为 None。

    返回 dict：
      `n_changes` 迁移次数（**不含**涉及缺数据状态的边）
      `n_skipped_missing` 因缺数据被跳过的边数（必须报出来，否则"迁移少"
        会被误读成"迁移真的少"）
      `n_edges` 总边数 / `n_edges_missing` 两端状态至少一端缺失的边数
      `counts` {(from,to): 次数} / `rows` 逐条事件明细

    ⚠️ 这是 **ML 点估计**（每节点取边际最大者），不是后验。节点边际概率低的
    分支上，方向随时可能翻转 —— 置信分布必须一起看（`confidence_summary`）。
    """
    counts, rows = {}, []
    n_edges = n_missing = 0

    def _ok(s):
        return s is not None and str(s) not in MISSING_STATES

    for n in iter_pre(root):
        ps = state_of.get(n['name'])
        for c in n['children']:
            n_edges += 1
            cs = state_of.get(c['name'])
            if not (_ok(ps) and _ok(cs)):
                n_missing += 1
                continue
            if ps == cs:
                continue
            counts[(ps, cs)] = counts.get((ps, cs), 0) + 1
            rows.append({
                'from': ps, 'to': cs,
                'child': c['name'] or '(node)',
                'parent': n['name'] or '(root)',
                'year': (round(node_year[id(c)], 3)
                         if node_year and node_year.get(id(c)) is not None
                         else None),
                'branch_len': float(c.get('length') or 0.0),
            })
    return {'n_changes': sum(counts.values()), 'counts': counts, 'rows': rows,
            'n_edges': n_edges, 'n_edges_missing': n_missing,
            'n_skipped_missing': n_missing}


def corridors_table(ml_counts, fitch_counts=None):
    """把 ML 与 Fitch 两条口径的走廊**并列**成表（不混排名、不合并计数）。

    返回 `[{from,to,ml,fitch,delta}]`，按 ML 降序、Fitch 缺失置 None。
    ⚠️ 两列数字含义不同：`ml` = TreeTime 的 ML 点估计，`fitch` = 平台
    `phylogeo._fitch` 的最少变化数下的悬挂（tie-break 依赖树序）。
    它们**不是同一个统计量**，只能看趋势一致性，不能相加、不能平均。
    """
    keys = set(ml_counts or {}) | set(fitch_counts or {})
    rows = []
    for k in keys:
        a = (ml_counts or {}).get(k)
        b = (fitch_counts or {}).get(k)
        rows.append({'from': k[0], 'to': k[1], 'ml': a, 'fitch': b,
                     'delta': (a - b) if (a is not None and b is not None) else None})
    rows.sort(key=lambda r: (-(r['ml'] or 0), -(r['fitch'] or 0), r['from'], r['to']))
    return rows


def _state_counts_from_csv(path, col=1):
    """状态表 → 各区划条数（算采样偏倚用）。"""
    counts = {}
    try:
        with io.open(path, encoding='utf-8-sig', errors='replace') as f:
            f.readline()
            for ln in f:
                cells = [c.strip().strip('"') for c in ln.rstrip('\n').split(',')]
                if len(cells) > col and cells[col]:
                    counts[cells[col]] = counts.get(cells[col], 0) + 1
    except OSError:
        return counts
    return counts


# ---------------------------------------------------------------- 运行

def _base_args(argv, tree=None, aln=None, dates=None, out_dir=None,
               rng_seed=0, verbose=1, name_column='name', date_column='date'):
    """拼公共参数（顺序不影响 TreeTime 解析）。"""
    cmd = list(argv)
    if tree:
        cmd += ['--tree', tree]
    if aln:
        cmd += ['--aln', aln]
    if dates:
        cmd += ['--dates', dates, '--name-column', name_column,
                '--date-column', date_column]
    if rng_seed is not None:
        cmd += ['--rng-seed', str(int(rng_seed))]
    if out_dir:
        cmd += ['--outdir', out_dir]
    if verbose is not None:
        cmd += ['--verbose', str(int(verbose))]
    return cmd


def run_clock(argv, out_dir, tree=None, aln=None, dates=None, seq_len=None,
              reroot='least-squares', keep_root=False, relax=None,
              coalescent=None, n_skyline=None, gen_per_year=None,
              confidence=True, time_marginal=None, clock_filter=None,
              allow_negative_rate=False, gtr=None, max_iter=None,
              branch_length_mode=None, verbose=None, rng_seed=0, timeout=3600,
              say=None, log=None, clock_rate=None):
    """TreeTime 分子钟定年（＋可选天际线）。

    ⚠️ **一律走「顶层无子命令」形式**（坑 1）：`treetime clock` 是缩水接口，
    连 `--confidence` / `--time-marginal` 都不认。调用方不用关心。

    `relax`: `(sigma, coupling)`；`(1.0, 0)` ＝**非相关松弛钟**（可拿 CoV 的近亲），
    不给就是 TreeTime 默认的严格钟假设。
    `reroot`: `'least-squares'` / `'min_dev'` / `'oldest'` / `'best'` / None；
    `keep_root=True` 时加 `--keep-root`（原始根位有意义时用，如已有外群）。
    `time_marginal`: `'only-final'` 最省（与 TreeTime 推荐一致）。

    `gtr` / `max_iter` / `branch_length_mode`：**2026-09-18 为 `#t-rtt` 的
    TreeTime-RTT 补的三个口子**（默认 None ＝ 一个 flag 都不加，原有调用方
    行为**逐位不变**）。它们对应 `treetime.TreeTime` 构造与 `_run()` 的
    同名参数，用于逐条对齐 VirPhyKit 的调用：
      · `gtr='JC69'` ≡ VirPhyKit 的 `gtr='Jukes-Cantor'`
        （实测 `GTR.standard('JC69')` 与 `GTR.standard('Jukes-Cantor')` 的
        `W` / `mu` / `Pi` 全等）；
      · `branch_length_mode='joint'` ≡ `run(branch_length_mode='joint')`；
      · `max_iter=3` ≡ `run(max_iter=3)`。
    ⚠️ `--gtr` 不给时 treetime 默认是 `infer`（**推断** GTR 模型），
    不是 JC69 —— 要"固定成 JC69"必须显式传。
    ⚠️⚠️ 2026-10-03 实测并落定：本函数（**CLI 顶层通路**）的 `--gtr` **传了也不
    生效** —— treetime 0.12.1 的 `wrappers.py` 把 `infer_gtr` 算出来后没有传进
    `TreeTime.run`（详见 `API_DRIVER_SRC` 的注释），替换模型**一律被推断**。
    所以返回结果里**不再**写 `'gtr': <入参>`（那会被下游误读成"这批数字用的是
    固定模型"），改成如实的两项：`gtr_used=False`（固定模型确实没生效）+
    `gtr_requested=<入参>`（调用方原本请求的值，仅供追溯）。要真正固定替换模型
    必须走 `run_clock_api()`（API 旁路，`infer_gtr` 已接线）。

    `verbose`：`None` ＝ 沿用 1（原行为）；`4` ＝ 打**逐轮** `rate=` / `R^2=`。
    TreeTime-RTT 用 4 是为了留下收敛过程 —— `molecular_clock.txt` 只写收敛值，
    且 R² 被 treetime 写死成 `%1.2f`（`utils.py:32`），只有 2 位小数。
    `clock_rate`: **固定分子钟速率**（替换/位点/年）。给了就**关掉速率估计**
    （TreeTime 帮助原文："if specified, the rate of the molecular clock won't be
    [estimated]"），直接用该值定年 —— 这是**慢演化病原**（mpox 常设 `6e-5`）
    时间信号太弱、自动拟合不出钟时的**官方兜底**，对齐 phymap-workflow 的
    `treetime_clock_rate`（见 `git-repo/phymap-workflow/README.md` 的
    "Things to Know"）。
    ⚠️ **口径**：返回值里的 `clock_rate_fixed` 标明"这个速率是用户给的、
    不是 ML 估出来的"；引用时必须一起写，否则会被误读成"数据支持这个速率"。
    与 `relax` 同时给属于自相矛盾（固定速率 + 速率异质性），调用方自己负责。

    返回 dict（可 JSON 序列化）。**失败不抛异常**，原因写进 `error`。
    """
    say = say or (lambda *_a, **_k: None)
    os.makedirs(out_dir, exist_ok=True)
    res = {'ok': False, 'error': None, 'out_dir': out_dir, 'argv': None,
           'returncode': None, 'clock': {}, 'relax': list(relax) if relax else None,
           'coalescent': None, 'skyline': [], 'skyline_meta': {},
           'outputs': {}, 'warnings': [], 'stdout_tail': '',
           'tree_in': tree, 'n_skyline': n_skyline,
           # ⚠️ 替换模型（2026-10-03 修复）：CLI 通路的 `--gtr` 是**空转**的 ——
           # treetime 0.12.1 顶层 CLI 把 `infer_gtr=True` 写死（wrappers.py 算出的
           # infer_gtr 没传进 `TreeTime.run`，见 API_DRIVER_SRC 与 docstring），
           # 固定模型从未生效。旧实现把入参原样记成 `'gtr': 'JC69'`，下游会误读成
           # "这批数字用过 JC69"。现在拆成两项如实记录：
           #   gtr_used=False       —— 固定模型**没有**生效（模型一律被推断）
           #   gtr_requested=<入参> —— 调用方原本请求的值，仅作追溯
           # 真正要固定模型请走 run_clock_api()。
           'gtr_used': False, 'gtr_requested': gtr,
           'max_iter': max_iter,
           'branch_length_mode': branch_length_mode, 'verbose': verbose,
           # 定年口径的**实际生效值**（2026-10-03 补，供报告/结果文件追溯）：
           # clock_filter None ＝ 没传 flag（treetime 默认 4.0 剔叶）；
           # time_marginal None ＝ 没传；reroot 记 'keep-root' 或 reroot 串/None。
           'clock_filter': (int(clock_filter) if clock_filter is not None
                            else None),
           'time_marginal': (str(time_marginal) if time_marginal else None),
           'reroot': ('keep-root' if keep_root
                      else (str(reroot) if reroot else None)),
           # 用户指定的固定速率（None＝照常做 ML 估计）。见 docstring 的口径说明。
           'clock_rate_fixed': (float(clock_rate) if clock_rate else None)}
    if not argv:
        res['error'] = 'TreeTime 不可用（见 check_tt_env）'
        return res
    if not tree:
        res['error'] = '缺少输入树'
        return res

    # 顶层形式：不带子命令（带 `clock` 会丢掉 --confidence/--time-marginal/
    # --coalescent/--relax，见坑 1）
    cmd = list(argv)
    cmd += (['--keep-root'] if keep_root else
            (['--reroot', str(reroot)] if reroot else []))
    cmd = _base_args(cmd, tree=tree, aln=aln, dates=dates, out_dir=out_dir,
                     rng_seed=rng_seed,
                     # None ＝ 沿用 _base_args 的 verbose=1（原有调用方行为不变）。
                     # TreeTime-RTT 传 4：日志里才有**逐轮** `rate=` / `R^2=`
                     # （`molecular_clock.txt` 只写收敛值，且 R² 只有 2 位小数）。
                     **({} if verbose is None else {'verbose': int(verbose)}))
    if seq_len and not aln:
        cmd += ['--sequence-length', str(int(seq_len))]
    # 固定替换模型（不传时 treetime 默认 `--gtr infer`，会自己推断模型）
    if gtr:
        cmd += ['--gtr', str(gtr)]
    if branch_length_mode:
        cmd += ['--branch-length-mode', str(branch_length_mode)]
    # 固定速率：给了就跳过速率估计（慢演化病原的兜底；对齐 phymap 的
    # treetime_clock_rate）。用 str(float()) 而非 %g —— 6e-5 这类小数要原样传。
    if clock_rate:
        cmd += ['--clock-rate', str(float(clock_rate))]
    if relax:
        cmd += ['--relax', str(float(relax[0])), str(float(relax[1]))]
    if confidence:
        cmd += ['--confidence']
    if time_marginal:
        cmd += ['--time-marginal', str(time_marginal)]
    if max_iter is not None:
        cmd += ['--max-iter', str(int(max_iter))]
    if clock_filter is not None:
        cmd += ['--clock-filter', str(clock_filter)]
    if allow_negative_rate:
        cmd += ['--allow-negative-rate']
    if coalescent is not None:
        cmd += ['--coalescent', str(coalescent)]
        if n_skyline:
            cmd += ['--n-skyline', str(int(n_skyline))]
        if gen_per_year:
            cmd += ['--gen-per-year', str(float(gen_per_year))]
        res['coalescent'] = str(coalescent)
    res['argv'] = ' '.join(str(c) for c in cmd)

    rc, lines = _run_tt(cmd, timeout=timeout, say=say)
    res['returncode'] = rc
    res['stdout_tail'] = '\n'.join(lines[-60:])

    clock_txt = os.path.join(out_dir, 'molecular_clock.txt')
    res['clock'] = parse_clock_txt(clock_txt)
    # TreeTime 的输出文件**无树就什么都没写** —— 用「树文件在不在」判成败，
    # 而不是只看退出码（实测退出码 0 但树没写出的情况存在：日期全没匹配上时）
    tt_nex = os.path.join(out_dir, 'timetree.nexus')
    for tag, fn in (('timetree_nexus', 'timetree.nexus'),
                    ('timetree_nwk', 'timetree.nwk'),
                    ('divergence_nexus', 'divergence_tree.nexus'),
                    ('molecular_clock', 'molecular_clock.txt'),
                    ('sequence_model', 'sequence_evolution_model.txt'),
                    ('rtt_pdf', 'root_to_tip_regression.pdf'),
                    ('tree_pdf', 'timetree.pdf'),
                    ('dates_tsv', 'dates.tsv'),
                    ('auspice', 'auspice_tree.json'),
                    ('ancestral_fasta', 'ancestral_sequences.fasta'),
                    ('skyline_tsv', 'skyline.tsv'),
                    ('skyline_pdf', 'skyline.pdf'),
                    ('trace_log', 'trace_run.log')):
        fp = os.path.join(out_dir, fn)
        if os.path.isfile(fp):
            res['outputs'][tag] = fp
    if os.path.isfile(tt_nex):
        try:
            with io.open(tt_nex, encoding='utf-8', errors='replace') as f:
                nwk = nexus_tree_newick(f.read())
            p = os.path.join(out_dir, 'timetree.nwk')
            with io.open(p, 'w', encoding='utf-8', newline='\n') as f:
                f.write(nwk + '\n')
            res['outputs']['timetree_nwk'] = p
            # ⚠️ 真实叶数（2026-10-03 修复）：旧写
            # `len(re.findall(...)) and None` 是恒为 None 的死表达式（没有任何
            # 调用方读得到数）。这里用本模块已有的树解析器数"无子节点的节点"。
            _root = parse_annotated_newick(nwk)
            res['timetree_tips'] = sum(1 for _n in iter_pre(_root)
                                       if not _n['children'])
        except (OSError, ValueError) as e:
            res['warnings'].append(f'timetree.nwk 未导出：{type(e).__name__}: {e}')

    if res['coalescent'] and 'skyline_tsv' in res['outputs']:
        rows, meta = parse_skyline_tsv(res['outputs']['skyline_tsv'])
        res['skyline'] = rows
        res['skyline_meta'] = meta

    if not res['clock'] and not os.path.isfile(tt_nex):
        res['error'] = (f'TreeTime 没产出时间树（rc={rc}）。常见原因：'
                        '① 叶名与日期表对不上；② 树上枝长不是替换/位点；'
                        f'③ 序列数太少。末段输出：{res["stdout_tail"][-400:]}')
        return res
    if not res['clock']:
        res['warnings'].append(
            '没读到 `molecular_clock.txt` —— 根到尾回归（R²/斜率）缺项。'
            'R² 是判断"这棵树有没有时间信号"的第一道关，缺了必须人工补看')
    res['ok'] = True
    say(f'TreeTime 定年完成：rate={res["clock"].get("rate")}, '
        f'r²={res["clock"].get("r2")}')
    return res


# ------------------------------------------------ API 旁路（能固定替换模型）

API_DRIVER_NAME = '_treetime_api_driver.py'

# 驱动脚本源码：**只补 `infer_gtr` 一个参数**，参数解析 / 产物写出 / 日志格式 /
# 绘图 / `--verbose` 行为全部照走 treetime 自己的 CLI 代码路径 —— 这样产物与
# CLI 逐字节可比，差异只剩我们故意补的这一处。落到 `out_dir` 是为了让
# 「产生这批数字的代码」和数字待在一起（可单独重跑）。
API_DRIVER_SRC = r'''# -*- coding: utf-8 -*-
"""TreeTime CLI 的同一性驱动 —— 只补一个被 CLI 漏掉的参数。

为什么需要它
------------
treetime 0.12.1 的 CLI **顶层形式**（`treetime --tree … --aln …`）永远推断
替换模型，`--gtr` 是空转：

* `treetime/treetime.py:86`   `_run(self, …, infer_gtr=True, …)`，默认 **True**；
* `treetime/wrappers.py:427`  算出了 `infer_gtr = params.gtr == 'infer'`；
* `treetime/wrappers.py:491`  但 `myTree.run(root=root, …)` **没把它传进去**
  （该变量只喂给 520 行的输出判断）。

于是 `--gtr JC69` 只设了构造函数的初始模型，随即被 `infer_gtr=True` 覆盖。
实测（同数据、同 `--rng-seed 0`，只换 `--gtr JC69` / `--gtr infer`）：
`divergence_tree.nexus` / `timetree.nexus` / `auspice_tree.json` / `dates.tsv`
**逐字节相同**，日志计算行逐字节相同；且 `--gtr JC69` 时
`sequence_evolution_model.txt` 根本不写（520 行判 False）—— 模型被静默丢弃。

本脚本做什么
------------
**只补 `infer_gtr` 这一个参数**（外加一个可选开关，见下），其余一切照走
treetime 自己的 CLI 代码路径。

私有开关 `--api-pre-reroot`（在交给 treetime 解析前摘掉）：在 `run()` **之前**
先 `reroot()` 一次，照 VirPhyKit 的调用顺序（`tt.reroot('least-squares')` 然后
`tt.run(root='best')`）。⚠️ 这不是多余的：`_run()` 的第一件事是
`optimize_tree(max_iter=1)`（`treetime.py:243`），它按**当时的根位**优化枝长
—— 先 reroot 会换掉这个起点。实测（80 叶夹具，固定 JC69）：
预 reroot `3.0622e-03` vs 不预 `3.0877e-03`，**差 0.84%**。

用法（参数与 treetime CLI 顶层形式**完全一致**）：

    python _treetime_api_driver.py --tree … --aln … --dates … --gtr JC69 …

调用方：`Virus_Platform_Core/treetime_ml.py: run_clock_api()`。
"""
import sys

import treetime


def main(argv):
    from treetime import argument_parser

    # 本驱动私有的开关：treetime 不认，解析前先摘掉。
    pre_reroot = False
    if '--api-pre-reroot' in argv:
        argv = [a for a in argv if a != '--api-pre-reroot']
        pre_reroot = True

    params = argument_parser.make_parser().parse_args(argv)

    # 语义与 wrappers.py:427 完全一致：
    #   --gtr infer   → 推断（＝CLI 的原行为，逐位不变）
    #   --gtr <模型>  → 固定成该模型（CLI 做不到的那一半）
    infer_gtr = params.gtr == 'infer'
    reroot_to = None if params.keep_root else params.reroot

    orig_run = treetime.TreeTime.run

    def run_with_fixes(self, *a, **kw):
        kw.setdefault('infer_gtr', infer_gtr)
        if pre_reroot and reroot_to is not None:
            # 照上游顺序：先 reroot，再 run（见文件头「私有开关」）
            self.reroot(root=reroot_to)
        return orig_run(self, *a, **kw)

    treetime.TreeTime.run = run_with_fixes
    try:
        return params.func(params)
    finally:
        treetime.TreeTime.run = orig_run


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]) or 0)
'''


def tt_python(argv):
    """从 `find_treetime()` 的 argv 前缀反推「装了 TreeTime 的解释器」。

    `find_treetime()` 给 `[python.exe, '-m', 'treetime']` 时直接就是它；给
    `[treetime.exe]` 时取同级 `Scripts` 的上级 `python.exe`（与 `check_tt_env()`
    第 248-255 行的推断口径一致）。推不出来返回 None —— **不抛异常**。
    """
    if not argv:
        return None
    if len(argv) >= 3 and argv[1] == '-m':
        return argv[0]
    if str(argv[0]).lower().endswith('.exe'):
        cand = os.path.join(os.path.dirname(os.path.dirname(argv[0])), 'python.exe')
        return cand if os.path.isfile(cand) else None
    return None


def write_api_driver(out_dir):
    """把 API 驱动脚本写到 `out_dir`（与产物放一起，可事后核对 / 单独重跑）。"""
    p = os.path.join(out_dir, API_DRIVER_NAME)
    with io.open(p, 'w', encoding='utf-8', newline='\n') as f:
        f.write(API_DRIVER_SRC)
    return p


def run_clock_api(out_dir, tree=None, aln=None, dates=None, *,
                  gtr='JC69', infer_gtr=False, reroot='least-squares',
                  keep_root=False, max_iter=3, branch_length_mode='joint',
                  time_marginal='true', clock_filter=0, seq_len=None,
                  rng_seed=0, verbose=4, pre_reroot=False, timeout=3600,
                  say=None, log=None, py=None):
    """TreeTime 分子钟定年（**Python API 旁路**）—— 能真正固定替换模型。

    ⚠️ **为什么必须有这条路（已实测，不是推测）**：CLI 顶层形式**永远推断
    替换模型** —— `_run()` 的 `infer_gtr` 默认 True（`treetime.py:86`），而
    `wrappers.py:427` 算出的 `infer_gtr` 没被传进 `wrappers.py:491` 的
    `myTree.run()`。实测 `--gtr JC69` 与 `--gtr infer` 的
    `divergence_tree.nexus` / `timetree.nexus` / `auspice_tree.json`
    **逐字节相同** ⇒ **`run_clock(gtr=…)` 这个参数从来没生效过**。
    所以「固定替换模型」只能走 API。
    **`run_clock()` 本身一行未改**（老师 2026-09-18 定：只记录，先不改）。

    实现不是重写一遍调用，而是**驱动 treetime 自己的 CLI 代码路径**
    （见 `API_DRIVER_SRC`：给 `TreeTime.run` 补上 `infer_gtr` 后再调
    `params.func(params)`）—— 参数解析 / 产物 / 日志 / 绘图与原 CLI 同源，
    差异只剩我们故意补的这一处。驱动脚本会落盘到 `out_dir` 供事后核对。

    `infer_gtr=False` + `gtr='JC69'` ≡ VirPhyKit 的
    `gtr='Jukes-Cantor'` + `infer_gtr=False`（实测 `GTR.standard('JC69')` 与
    `GTR.standard('Jukes-Cantor')` 的 `W` / `mu` / `Pi` 全等）。

    `clock_filter=0`：treetime CLI 默认 `4.0`（剔偏离钟 >4 IQR 的叶），
    Python API 不传 `n_iqd` ＝ **不筛** —— 0 是假值，等价于不筛。

    `pre_reroot=True`：在 `run()` **之前**先 `reroot()` 一次，照 VirPhyKit 的
    调用顺序。⚠️ 实测这**不是**多余的（80 叶夹具、固定 JC69）：预 reroot
    `3.0622e-03` vs 不预 `3.0877e-03`，**差 0.84%** —— 因为 `_run()` 的第一件
    事是 `optimize_tree(max_iter=1)`（`treetime.py:243`），它按**当时的根位**
    优化枝长。默认 `False`（＝与 CLI 行为一致，只补 `infer_gtr`）。

    返回 dict，键与 `run_clock()` 对齐（`clock` / `outputs` / `argv` /
    `returncode` / `stdout_tail` / `warnings` / `error` / `ok`），另加
    `engine='treetime-api'` / `api_python` / `api_driver` / `infer_gtr` /
    `gtr` 等溯源字段。**失败不抛异常**，原因写进 `error`。
    """
    say = say or (lambda *_a, **_k: None)
    log = log or say
    os.makedirs(out_dir, exist_ok=True)
    res = {'ok': False, 'error': None, 'engine': 'treetime-api',
           'out_dir': out_dir, 'argv': None, 'returncode': None,
           'stdout_tail': '', 'clock': {}, 'outputs': {}, 'warnings': [],
           'api_python': None, 'api_driver': None,
           'infer_gtr': bool(infer_gtr), 'gtr': gtr,
           'reroot': None if keep_root else reroot,
           'keep_root': bool(keep_root), 'max_iter': max_iter,
           'branch_length_mode': branch_length_mode,
           'time_marginal': time_marginal, 'clock_filter': clock_filter,
           'rng_seed': rng_seed, 'verbose': verbose, 'seq_len': seq_len,
           'pre_reroot': bool(pre_reroot)}

    argv0, how = find_treetime()
    if not argv0:
        res['error'] = f'TreeTime 不可用：{how}'
        return res
    res['treetime'] = how
    res['api_python'] = py or tt_python(argv0)
    if not res['api_python'] or not os.path.isfile(res['api_python']):
        res['error'] = (f'推不出装了 TreeTime 的解释器（find_treetime 给的是 '
                        f'{argv0!r}）—— API 旁路要的是解释器本体，不是 .exe 包装。'
                        '可在 platform.json 的 tools.treetime 指到 python.exe。')
        return res

    if not tree:
        res['error'] = '缺少输入树'
        return res

    res['api_driver'] = write_api_driver(out_dir)

    # 参数与 `run_clock` 的 CLI 调用**逐条同源**（只取 RTT 需要的那几个），
    # 免得两处各写一份矩阵、日后漂移。
    cmd = [res['api_python'], res['api_driver']]
    cmd = _base_args(cmd, tree=tree, aln=aln, dates=dates, out_dir=out_dir,
                     rng_seed=rng_seed,
                     **({} if verbose is None else {'verbose': int(verbose)}))
    if seq_len and not aln:
        cmd += ['--sequence-length', str(int(seq_len))]
    if keep_root:
        cmd += ['--keep-root']
    elif reroot:
        cmd += ['--reroot', str(reroot)]
    if pre_reroot:
        cmd += ['--api-pre-reroot']
    if gtr:
        cmd += ['--gtr', str(gtr)]
    if branch_length_mode:
        cmd += ['--branch-length-mode', str(branch_length_mode)]
    if time_marginal:
        cmd += ['--time-marginal', str(time_marginal)]
    if max_iter is not None:
        cmd += ['--max-iter', str(int(max_iter))]
    if clock_filter is not None:
        cmd += ['--clock-filter', str(clock_filter)]
    res['argv'] = ' '.join(str(c) for c in cmd)

    rc, lines = _run_tt(cmd, timeout=timeout, say=say)
    res['returncode'] = rc
    res['stdout_tail'] = '\n'.join(lines[-60:])

    clock_txt = os.path.join(out_dir, 'molecular_clock.txt')
    res['clock'] = parse_clock_txt(clock_txt)
    tt_nex = os.path.join(out_dir, 'timetree.nexus')
    for tag, fn in (('timetree_nexus', 'timetree.nexus'),
                    ('timetree_nwk', 'timetree.nwk'),
                    ('divergence_nexus', 'divergence_tree.nexus'),
                    ('molecular_clock', 'molecular_clock.txt'),
                    ('sequence_model', 'sequence_evolution_model.txt'),
                    ('rtt_pdf', 'root_to_tip_regression.pdf'),
                    ('tree_pdf', 'timetree.pdf'),
                    ('dates_tsv', 'dates.tsv'),
                    ('auspice', 'auspice_tree.json'),
                    ('ancestral_fasta', 'ancestral_sequences.fasta'),
                    ('trace_log', 'trace_run.log')):
        fp = os.path.join(out_dir, fn)
        if os.path.isfile(fp):
            res['outputs'][tag] = fp
    if os.path.isfile(tt_nex):
        try:
            with io.open(tt_nex, encoding='utf-8', errors='replace') as f:
                nwk = nexus_tree_newick(f.read())
            p = os.path.join(out_dir, 'timetree.nwk')
            with io.open(p, 'w', encoding='utf-8', newline='\n') as f:
                f.write(nwk + '\n')
            res['outputs']['timetree_nwk'] = p
        except (OSError, ValueError) as e:
            res['warnings'].append(f'timetree.nwk 未导出：{type(e).__name__}: {e}')

    # 口径自检：既然要求固定模型，TreeTime 就**不该**写出推断出来的模型文件。
    # 它出现了说明补丁没生效（模型被推断覆盖）—— 必须报警，不能静默出数。
    if not infer_gtr and 'sequence_model' in res['outputs']:
        res['warnings'].append(
            '⚠️ infer_gtr=False 但 TreeTime 写出了 sequence_evolution_model.txt '
            '—— 替换模型被推断了，固定模型没生效，这批数字的模型口径存疑')

    if not res['clock'] and not os.path.isfile(tt_nex):
        res['error'] = (f'TreeTime 没产出时间树（rc={rc}）。末段输出：'
                        f'{res["stdout_tail"][-400:]}')
        return res
    if not res['clock']:
        res['warnings'].append(
            '没读到 `molecular_clock.txt` —— 根到尾回归（R²/斜率）缺项')
    res['ok'] = True
    say(f'TreeTime 定年完成（API 旁路，替换模型 {gtr}'
        + ('固定' if not infer_gtr else '推断')
        + f'）：rate={res["clock"].get("rate")}, r²={res["clock"].get("r2")}')
    return res


def run_mugration(argv, tree, states, out_dir, attribute='region',
                  confidence=True, sampling_bias_correction=None, pc=None,
                  weights=None, missing_data=None, rng_seed=0, timeout=3600,
                  say=None, also_tree=None, also_out_dir=None):
    """TreeTime mugration：离散性状（区划）最大似然祖先态重构 + GTR 迁移模型。

    `also_tree`：**第二棵树做对拍**（见模块头「第五个坑」）。通常传
    `divergence_tree.nexus`，默认树传 `timetree.nwk`。两棵树的迁移次数一起报出，
    分歧 > 20% 时置 `cross_check['divergent']=True` 并告警。

    返回 dict：`gtr`（mapping/mu/pi/W/Q）、`confidence`（rows）、
    `confidence_summary`、`state_of`（节点名 → 状态，argmax）、
    `annotated_nexus`、`outputs`、`warnings`、`cross_check`。失败原因进 `error`。
    """
    say = say or (lambda *_a, **_k: None)
    os.makedirs(out_dir, exist_ok=True)
    res = {'ok': False, 'error': None, 'out_dir': out_dir, 'argv': None,
           'returncode': None, 'attribute': attribute, 'gtr': {},
           'confidence': [], 'confidence_summary': {}, 'state_of': {},
           'outputs': {}, 'warnings': [], 'stdout_tail': '',
           'sampling_bias_correction': sampling_bias_correction,
           'tree_in': tree, 'cross_check': None}
    if not argv:
        res['error'] = 'TreeTime 不可用（见 check_tt_env）'
        return res
    if not tree or not states:
        res['error'] = '缺少输入树或状态表'
        return res

    def _cmd(sbc, tree_path, odir):
        c = list(argv) + ['mugration', '--tree', tree_path, '--states', states,
                          '--attribute', attribute,
                          '--name-column', 'name', '--outdir', odir]
        if confidence:
            c += ['--confidence']
        if sbc is not None:
            c += ['--sampling-bias-correction', str(float(sbc))]
        if pc is not None:
            c += ['--pc', str(float(pc))]
        if weights:
            c += ['--weights', weights]
        if missing_data:
            c += ['--missing-data', str(missing_data)]
        if rng_seed is not None:
            c += ['--rng-seed', str(int(rng_seed))]
        c += ['--verbose', '1']
        return c

    sbc = sampling_bias_correction
    if isinstance(sbc, str) and sbc.strip().lower() == 'auto':
        sbc = None                       # 第一遍不带校正，读完 GTR 再算
        res['sampling_bias_correction'] = 'auto'
        res['_auto_sbc'] = True

    cmd = _cmd(sbc, tree, out_dir)
    res['argv'] = ' '.join(str(c) for c in cmd)
    rc, lines = _run_tt(cmd, timeout=timeout, say=say)
    res['returncode'] = rc
    res['stdout_tail'] = '\n'.join(lines[-40:])

    gtr_path = os.path.join(out_dir, 'GTR.txt')
    res['gtr'] = parse_gtr_txt(gtr_path)
    for tag, fn in (('gtr_txt', 'GTR.txt'), ('confidence_csv', 'confidence.csv'),
                    ('annotated_nexus', 'annotated_tree.nexus')):
        fp = os.path.join(out_dir, fn)
        if os.path.isfile(fp):
            res['outputs'][tag] = fp

    if not res['gtr'].get('mapping'):
        res['error'] = (f'mugration 没产出可解析的 GTR.txt（rc={rc}）。'
                        '常见原因：① 状态表叶名与树对不上；'
                        '② 有效状态 < 2；③ 状态表缺 name 列。'
                        f'末段输出：{res["stdout_tail"][-400:]}')
        return res

    # ---- 采样偏倚校正：auto 的第二遍 ----
    if res.pop('_auto_sbc', False):
        sug = suggest_sampling_bias_correction(res['gtr'], states)
        res['sampling_bias_correction_suggested'] = sug
        if sug:
            say(f'采样偏倚校正系数自动估算：{sug:.4g}'
                '（(1−Σpᵢ²)/(1−Σtᵢ²)，pᵢ 为均衡频率、tᵢ 为叶上观测量）')
            cmd2 = _cmd(sug, tree, out_dir)
            rc2, lines2 = _run_tt(cmd2, timeout=timeout, say=say)
            res['returncode'] = rc2
            res['argv'] = ' '.join(str(c) for c in cmd2)
            res['sampling_bias_correction'] = sug
            res['gtr'] = parse_gtr_txt(gtr_path)
            res['stdout_tail'] = '\n'.join(lines2[-40:])
        else:
            res['warnings'].append('采样偏倚校正系数算不出（状态数 < 2）→ 本次未校正')
            res['sampling_bias_correction'] = None

    if confidence:
        rows, summ = parse_confidence_csv(
            os.path.join(out_dir, 'confidence.csv'), res['gtr'].get('mapping'))
        res['confidence'] = rows
        res['confidence_summary'] = summ
        # ⚠️ 条件必须是"**有**未映射的列"而不是"一条都没映射上"：
        #    部分未映射时（代号跨到标点）旧条件不触发 → 真地名与代号混在一张
        #    走廊表里也没人知道（2026-09-18 真跑踩到）。
        if (summ.get('n_labels_unmapped') or 0) and rows:
            res['warnings'].append(
                '%d/%d 个状态列没能用 GTR.txt 的字母映射翻译（未翻译的：%s）—— '
                '走廊表里这些**按原代号给**，不要当成区划名读'
                % (summ['n_labels_unmapped'], len(summ.get('confidence_labels_raw') or []),
                   ','.join((summ.get('confidence_labels_unmapped') or [])[:8])))
        # ⚠️ 原先这里还有一个 `elif not summ.get('confidence_labels_mapped') and rows:`
        #    分支（"状态列完全没能翻译"），其实**不可达**（2026-10-03 删除）：
        #    `confidence_labels_mapped = bool(mapping)`，mapping 为空 ⇒ 所有表头列
        #    都在 `n_labels_unmapped` 里；而 rows 非空 ⇒ 表头至少有 1 个状态列
        #    ⇒ 上面的 if 必中。留着只会误导后人以为它兜住了什么。
        res['state_of'] = {r['name']: r['best'] for r in rows}
        if summ.get('n_below_0_7'):
            res['warnings'].append(
                f'{summ["n_below_0_7"]}/{summ.get("n_internal")} 个内部节点的'
                f'ML 最大边际概率 < 0.7（最低 {summ.get("min_max_p"):.3g}）—— '
                '涉及这些节点的迁移方向不可靠，做方向性结论前先看置信分布')

    # ---- 第二棵树对拍（第五个坑的处置）----
    if also_tree and os.path.isfile(str(also_tree)):
        odir = also_out_dir or (out_dir.rstrip('/\\') + '_xcheck')
        os.makedirs(odir, exist_ok=True)
        say('对拍：换第二棵树再跑一遍 mugration（检验枝长口径敏感性）')
        _sbc2 = res['sampling_bias_correction']
        if isinstance(_sbc2, str):       # 'auto' 未解析成功 → 不带校正
            _sbc2 = None
        rc3, _l3 = _run_tt(_cmd(_sbc2, also_tree, odir),
                           timeout=timeout, say=say)
        g2 = parse_gtr_txt(os.path.join(odir, 'GTR.txt'))
        rows2, _s2 = parse_confidence_csv(os.path.join(odir, 'confidence.csv'),
                                          g2.get('mapping'))
        if g2.get('mapping') and rows2:
            so2 = {r['name']: r['best'] for r in rows2}
            res['cross_check'] = {'tree': str(also_tree), 'rc': rc3,
                                  'mu': g2.get('mu'), 'state_of': so2,
                                  'out_dir': odir}
            say(f'对拍第二棵树完成（rc={rc3}, mu={g2.get("mu")}）')
        else:
            res['cross_check'] = {'tree': str(also_tree), 'rc': rc3,
                                  'error': '第二棵树没跑出可解析结果'}
            res['warnings'].append(
                f'对拍的第二棵树（{os.path.basename(str(also_tree))}）没跑出结果，'
                '枝长口径敏感性未能检验')

    res['ok'] = True
    say(f'mugration 完成：mu={res["gtr"].get("mu")}, '
        f'状态 {len(res["gtr"].get("states") or [])} 个')
    return res


def leafset_key(node):
    """节点 → 其下全部叶名排序元组（树的**结构指纹**）。

    为什么用它当键：TreeTime 的 mugration 会把内部节点重命名成
    `NODE_0000123`（实测：输入树里叫 `NODE_0000142` 的节点，输出的标注树里
    变成 `NODE_0000005`），**两次运行的编号也不保证一致**。按节点名对齐会
    静默错位；按"该节点覆盖哪些叶"对齐则与命名无关，只依赖拓扑。
    """
    return tuple(sorted(t['name'] for t in iter_pre(node) if not t['children']))


def cross_check_corridors(root, cc, cc_tree_path=None, primary=None):
    """用对拍树的状态重算走廊，与主口径比较 → dict（供面板并排显示）。

    `root`：主树（时间树）根节点。
    `cc`：`run_mugration(also_tree=…)` 返回的 `cross_check`（含 `state_of`）。
    `cc_tree_path`：**对拍那棵树**的路径 —— 必须给，否则只能退化成按节点名对齐
    （见 `leafset_key` 的说明：mugration 会重命名内部节点，按名对齐会静默错位）。

    返回 dict：`n_changes` / `counts` / `n_skipped_missing` / `n_unmatched`
    / `rel_diff` / `divergent`（分歧 > 20%）。
    """
    if not cc or not cc.get('state_of'):
        return None
    so2 = cc['state_of']
    out = {'tree': os.path.basename(str(cc.get('tree') or '')),
           'mu': cc.get('mu'), 'n_changes': 0, 'counts': {},
           'n_skipped_missing': 0, 'n_unmatched': 0, 'aligned_by': 'node_name'}

    def _ok(s):
        return s is not None and str(s) not in MISSING_STATES

    # ---- 优先：按叶集合把对拍树的「结构 → 状态」建索引，再映射回主树 ----
    state_by_leafset = None
    if cc_tree_path and os.path.isfile(str(cc_tree_path)):
        try:
            with io.open(cc_tree_path, encoding='utf-8', errors='replace') as f:
                croot = parse_annotated_newick(f.read())
            state_by_leafset = {}
            for n in iter_pre(croot):
                st = so2.get(n['name'])
                if st is not None:
                    state_by_leafset[leafset_key(n)] = st
            out['aligned_by'] = 'leafset'
        except (OSError, ValueError):
            state_by_leafset = None

    for n in iter_pre(root):
        pn = n['name']
        if state_by_leafset is not None:
            ps = state_by_leafset.get(leafset_key(n))
        else:
            ps = so2.get(pn)
        for c in n['children']:
            if state_by_leafset is not None:
                cs = state_by_leafset.get(leafset_key(c))
            else:
                cs = so2.get(c['name'])
            if ps is None or cs is None:
                out['n_unmatched'] += 1
                continue
            if not (_ok(ps) and _ok(cs)):
                out['n_skipped_missing'] += 1
                continue
            if ps != cs:
                out['n_changes'] += 1
                # ⚠️ 键必须是字符串：这个 dict 会进 summary.json，元组键
                # 让 json.dump 直接抛 `keys must be str…, not tuple`
                k = '%s->%s' % (ps, cs)
                out['counts'][k] = out['counts'].get(k, 0) + 1
    if out['n_unmatched']:
        out['align_warning'] = (
            f'{out["n_unmatched"]} 个节点在对拍结果里找不到对应状态'
            f'（按 {out["aligned_by"]} 对齐）—— 对拍数字不可信，请人工核对')
    if primary:
        base = primary.get('n_changes') or 0
        out['primary_n_changes'] = base
        out['delta'] = out['n_changes'] - base
        out['rel_diff'] = (round(abs(out['n_changes'] - base) / base, 4)
                           if base else None)
        out['divergent'] = bool(base and out['rel_diff'] is not None
                                and out['rel_diff'] > 0.2)
    return out


def suggest_sampling_bias_correction(gtr, states_path):
    """按 TreeTime 官方文档的公式估算采样偏倚校正系数。

    公式：`(1 − Σpᵢ²) / (1 − Σtᵢ²)`，其中 pᵢ 是 GTR 的均衡频率、
    tᵢ 是叶上的**表观**频率（即实际采样条数占比）。
    """
    mapping = (gtr or {}).get('mapping') or {}
    pi = (gtr or {}).get('pi') or {}
    if len(mapping) < 2:
        return None
    counts = _state_counts_from_csv(states_path)
    total = sum(counts.values())
    if not total:
        return None
    t = [counts.get(s, 0) / total for s in mapping.values()]
    p = [pi.get(k) or 0.0 for k in mapping]
    denom = 1.0 - sum(x * x for x in t)
    num = 1.0 - sum(x * x for x in p)
    if denom <= 0:
        return None
    return round(num / denom, 6)


def run_ancestral(argv, aln, tree, out_dir, marginal=False, aa=False,
                  method_anc=None, gtr=None, rng_seed=0, timeout=3600, say=None):
    """TreeTime ancestral：祖先序列 + 逐枝替换。

    ⚠️ 产出文件名是 `ancestral_sequences.fasta`（**不是** TreeTime 文档里写的
    `ancestral.fasta`，0.12.1 实测），`branch_mutations.txt` 是制表符四列
    `node/state1/pos/state2`。
    """
    say = say or (lambda *_a, **_k: None)
    os.makedirs(out_dir, exist_ok=True)
    res = {'ok': False, 'error': None, 'out_dir': out_dir, 'argv': None,
           'returncode': None, 'n_mutations': 0, 'mutations': [],
           'outputs': {}, 'warnings': [], 'stdout_tail': ''}
    if not argv:
        res['error'] = 'TreeTime 不可用（见 check_tt_env）'
        return res
    if not aln or not tree:
        res['error'] = 'ancestral 需要比对 FASTA 与树两者'
        return res
    cmd = list(argv) + ['ancestral', '--aln', aln, '--tree', tree,
                        '--outdir', out_dir]
    if aa:
        cmd += ['--aa']
    if marginal:
        cmd += ['--marginal']
    if method_anc:
        cmd += ['--method-anc', str(method_anc)]
    if gtr:
        cmd += ['--gtr', str(gtr)]
    if rng_seed is not None:
        cmd += ['--rng-seed', str(int(rng_seed))]
    cmd += ['--verbose', '1']
    res['argv'] = ' '.join(str(c) for c in cmd)
    rc, lines = _run_tt(cmd, timeout=timeout, say=say)
    res['returncode'] = rc
    res['stdout_tail'] = '\n'.join(lines[-30:])
    for tag, fn in (('ancestral_fasta', 'ancestral_sequences.fasta'),
                    ('annotated_nexus', 'annotated_tree.nexus'),
                    ('branch_mutations', 'branch_mutations.txt'),
                    ('auspice', 'auspice_tree.json'),
                    ('sequence_model', 'sequence_evolution_model.txt')):
        fp = os.path.join(out_dir, fn)
        if os.path.isfile(fp):
            res['outputs'][tag] = fp
    if 'branch_mutations' in res['outputs']:
        rows = parse_branch_mutations(res['outputs']['branch_mutations'])
        res['mutations'] = rows
        res['n_mutations'] = len(rows)
    if 'ancestral_fasta' not in res['outputs']:
        res['error'] = (f'ancestral 没产出 ancestral_sequences.fasta（rc={rc}）。'
                        f'末段输出：{res["stdout_tail"][-400:]}')
        return res
    res['ok'] = True
    say(f'ancestral 完成：{res["n_mutations"]} 处逐枝替换')
    return res


def run_homoplasy(argv, aln, tree, out_dir, rng_seed=0, timeout=1800, say=None):
    """TreeTime homoplasy：替换饱和 / 同塑性检验（**结果只在 stdout**）。

    平台此前没有替换饱和检验；同塑性过高的位点集会让任何建树/定年结论打折，
    这张表是低成本的把关。
    """
    say = say or (lambda *_a, **_k: None)
    os.makedirs(out_dir, exist_ok=True)
    res = {'ok': False, 'error': None, 'outdir': out_dir, 'argv': None,
           'returncode': None, 'homoplasies': [], 'outputs': {},
           'warnings': [], 'stdout_tail': ''}
    if not argv:
        res['error'] = 'TreeTime 不可用（见 check_tt_env）'
        return res
    if not aln or not tree:
        res['error'] = 'homoplasy 需要比对 FASTA 与树两者'
        return res
    cmd = list(argv) + ['homoplasy', '--aln', aln, '--tree', tree,
                        '--outdir', out_dir]
    if rng_seed is not None:
        cmd += ['--rng-seed', str(int(rng_seed))]
    cmd += ['--verbose', '1']
    res['argv'] = ' '.join(str(c) for c in cmd)
    rc, lines = _run_tt(cmd, timeout=timeout, say=say)
    res['returncode'] = rc
    res['stdout_tail'] = '\n'.join(lines[-30:])
    res['homoplasies'] = parse_homoplasy_stdout(lines)
    if not res['homoplasies']:
        res['error'] = (f'homoplasy 没解析出同塑性表（rc={rc}）—— 该子命令不落文件，'
                        f'结果只在控制台。末段输出：{res["stdout_tail"][-300:]}')
        return res
    res['ok'] = True
    say(f'homoplasy 完成：同塑性最高的位点 {res["homoplasies"][0]["mut"]}'
        f'（{res["homoplasies"][0]["multiplicity"]} 次）')
    return res


# ---------------------------------------------------------------- 体检

def audit_clock(res, drt=None, n_tip=None):
    """定年结果体检 —— 面向页面讲清「能信到什么程度」。

    ⚠️ 红线：**时间信号没通过 DRT 检验时不给定年结论**。

    ⚠️ **R² 低 ≠ 没时间信号**（本模块实测更正过一次口径）：
    RSV 209 固定根位的根到尾 R² 只有 **0.049**，但把根位也纳入搜索后
    「最大根位 R²」= 0.1119，DRT 在 200 次置换下 **p=0.005 → 通过**。
    所以 R² 只是描述性的，**判定必须看 DRT**；反过来，R² 高也不等于通过。
    两者都给，别互相替代。
    """
    problems, warnings = [], list(res.get('warnings') or [])
    if not res.get('ok'):
        problems.append(res.get('error') or '未知失败')
        return {'ok': False, 'problems': problems, 'warnings': warnings,
                'verdict': '失败'}

    ck = res.get('clock') or {}
    rate, r2 = ck.get('rate'), ck.get('r2')
    if not (isinstance(rate, float) and rate > 0):
        problems.append(f'速率不是正数（{rate}）—— 结果不可用')

    if isinstance(r2, float) and r2 < 0.1:
        warnings.append(
            f'固定根位的根到尾回归 R²={r2:.3g} 偏低 —— 这**只说明根位可能没取好**，'
            '不等于没时间信号（RSV 209 实测：固定根位 R²=0.049，'
            '最大根位 R²=0.112，DRT p=0.005 **通过**）。'
            '请以 DRT 的 p 值为准，别用 R² 单独下结论')

    # DRT（日期随机化检验）：平台已有实现，这里只做**口径提示**，不重算
    if drt:
        p = drt.get('p_value', drt.get('p'))
        nv = drt.get('n_valid')
        if isinstance(p, float) and p >= 0.05:
            warnings.append(
                f'DRT 日期随机化检验 p={p:.3g} ≥ 0.05 —— **未通过**：'
                '观测到的"时间-遗传距离"关系与随机打乱日期后无法区分。'
                '本次定年数字仅作记录，不得作为定年结论引用')
        elif isinstance(p, float) and nv and p <= 1.0 / (nv + 1) + 1e-9:
            warnings.append(
                f'DRT p={p:.4g} 恰好等于置换次数的分辨率下限 '
                f'1/({nv}+1) —— 即"真实值超过全部置换"。'
                '想报更小的 p 必须加大置换次数（n_perm≥100 才能报 p<0.01）')
    elif drt is None:
        warnings.append('本次没有带上 DRT 日期随机化检验结果 —— 时间信号是否成立'
                        '尚未检验；R² 只是描述性的，不等价于 DRT')

    if res.get('coalescent') == 'skyline':
        meta = res.get('skyline_meta') or {}
        if not res.get('skyline'):
            warnings.append('指定了 skyline 但没解析出 skyline.tsv —— 天际线缺失')
        else:
            warnings.append(
                'skyline 的上下界是 TreeTime 的**近似区间**'
                f'（±{meta.get("sd_multiple", 2)} SD of the LH），**不是**贝叶斯 '
                '95% HPD；且 N_e 的绝对标度取决于 `--gen-per-year` 假设'
                f'（本次 {meta.get("gen_per_year", "?")} gen/year）。'
                '正式 95% 带仍需 BEAST 的 BSP（轨道 B）')
    if res.get('relax'):
        warnings.append(f'用了 `--relax {res["relax"][0]} {res["relax"][1]}`'
                        '（TreeTime 的松弛钟是**自相关/非相关先验下的枝端速率**，'
                        '与 treedater 的 CoV 不是同一个量，不可互相比较）')
    n = n_tip or res.get('n_tip')
    if n and n < 10:
        problems.append(f'定年树只有 {n} 个叶 —— 分子钟拟合的自由度不足，'
                        '速率/tMRCA 不可信')
    elif n and n < 20:
        warnings.append(f'只有 {n} 个叶，建议同时看 R² 与 DRT 的 p 值再下结论')
    return {'ok': not problems, 'problems': problems, 'warnings': warnings,
            'verdict': '通过' if not problems else '失败',
            'r2': r2, 'rate': rate}


def audit_mugration(res, n_regions=None, sampling_bias_correction=None):
    """地理重构体检。口径红线：ML 边际概率 ≠ 贝叶斯后验。"""
    problems, warnings = [], list(res.get('warnings') or [])
    if not res.get('ok'):
        problems.append(res.get('error') or '未知失败')
        return {'ok': False, 'problems': problems, 'warnings': warnings,
                'verdict': '失败'}
    n = n_regions if n_regions is not None else len(res.get('gtr', {}).get('states') or [])
    if n < 2:
        problems.append(f'有效区划只有 {n} 个 —— 迁移矩阵无意义')
    summ = res.get('confidence_summary') or {}
    if summ.get('frac_below_0_7'):
        f = summ['frac_below_0_7']
        if f > 0.5:
            warnings.append(
                f'{f:.0%} 的内部节点 ML 最大边际概率 < 0.7 —— '
                '多数节点的区划归属并不确定，走廊计数应视为**下界性的点估计**，'
                '做方向性结论前先做 RSPP bootstrap 或交 BEAST 求后验')
        elif f > 0.2:
            warnings.append(f'{f:.0%} 的内部节点 ML 最大边际概率 < 0.7 —— '
                            '部分迁移方向不可靠，面板里逐条标出')
    warnings.append('`confidence.csv` 给的是 **ML 边际概率**（该树 + 该 GTR 下的'
                    '似然边际），**不是贝叶斯后验**：没有先验、没有拓扑不确定性。'
                    '与 BEAST 的后验概率**不可比大小**')
    if sampling_bias_correction:
        warnings.append(f'已按采样偏倚系数 {sampling_bias_correction} 校正 —— '
                        '未校正时高置信度会虚高（RSV 实测未校正 max-p 中位数 '
                        '0.990，校正后 0.765）')
    else:
        warnings.append('**未做采样偏倚校正** —— 采样不均衡时置信度会系统性虚高。'
                        '把「采样偏倚校正」设为 auto 可让平台按 '
                        '(1−Σpᵢ²)/(1−Σtᵢ²) 自动估算')
    return {'ok': not problems, 'problems': problems, 'warnings': warnings,
            'verdict': '通过' if not problems else '失败'}
