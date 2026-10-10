# -*- coding: utf-8 -*-
r"""MJRM（Markov Jump Randomization Model）—— BEAST 1.x 的 **Markov jump 计数**。

## 先说清"做不做"这件事（2026-09-16 修正）

之前我把 **MJRM 生成器**与**跑 MCMC**混为一谈了，于是写了"不做 MatrixMJump（需 BEAST 类
MCMC）"。**这个说法不准确**：VirPhyKit 的 `src/MakovMJump/`（注意原作拼写）是**纯 Python**，
只做两件事 ——

1. `wrtcfg()` 生成 N×(N−1) 个 `<parameter id="A-to-B" value="{0/1 指示矩阵}"/>`
2. `wrt_rewards()` 生成 `<rewards>` 块
3. `process_xml()` 把上面两段插进一份**已有的** BEAST XML 模板

**它自己不跑 MCMC。** 真正需要 BEAST 的只是随后那次链（PVS 示例是 1×10⁸ 步）。
所以能在本地做的至少有三块：**生成配置 + 审计配置 + 解析结果**；只有"跑"这一步要 BEAST。

## 用真实 BEAST1 跑出来的实测结论（服务器 BEAST 1.10.4）

拿 VirPhyKit 自带示例 `Example/MJRM Generator/PVS/PVS_with_matrix.xml` 直接跑（链长压到
2×10⁶ 便于观察），主日志 **125 列里既没有任何 `A-to-B` 计数列，也没有 `_reward` 列** ——
也就是说 **插进 XML 的那一块是"哑"的**：指示矩阵被声明了，却没有任何 `<log>` 引用它们，
`<markovJumpsTreeLikelihood>` 开标签也没有 `rewards=` 属性。**跑完拿不到跳转计数。**

进一步对照实验（用 BEAST1 官方 `examples/TestXML/testMarkovJumps.xml` 的规范写法）：

| 写法 | 结果 |
|---|---|
| `<log><parameter idref="AS-to-OC"/></log>` | 能跑，但输出 `AS-to-OC1..36` 列，值恒为 **静态矩阵本身**（`0.0 1.0 0.0…`），**不是计数** |
| `<log><parameter idref="AS_reward"/></log>` | **BEAST1 直接解析失败** |
| `<log><treeLikelihood idref="…"/></log>` | BEAST1 官方示例的写法：把 likelihood 对象放进 log |

所以本模块的 `insert_into_template()` **不只照抄 VirPhyKit**，还顺手把"让计数真的出现"
的那一步补上（`with_log=True`），并给出 `audit_mjrm_xml()` —— **跑 BEAST 之前先检查
配置是否真的能产出结果**，别白等几小时。
"""

from __future__ import annotations

import io
import os
import re
from collections import OrderedDict

# VirPhyKit 的三个插入锚点（顺序即优先级）
_TARGET_MARKERS = (
    '<!--  END Ancestral state reconstruction',
    '<!-- END Ancestral state reconstruction',
    '</markovJumpsTreeLikelihood>',
)

_ROUTE_RE = re.compile(r'^([A-Za-z][\w.]*?)-to-([A-Za-z][\w.]*)$')
# Beauti 生成的 XML 会给「离散性状段」加这对注释标记 —— 拆「定年版」就靠它当锚点
_DT_START = '<!-- START Discrete Traits Model'
_DT_END = '<!-- END Discrete Traits Model'
_ROUTE_FLAT_RE = re.compile(r'^([A-Za-z][\w.]*?)-to-([A-Za-z][\w.]*?)(\d+)$')

# ⚠️ 实测（BEAST 1.10.4）真正的列名长这样（**带 `c_` 前缀和 `[1]` 下标**）：
#     c_Region.count[1]    总跳转次数
#     c_AS-to-OC[1]        逐路线期望跳转次数（**这就是我们要的计数**）
#     c_AS_reward[1]       各地点停留时间（reward = 1.0 在本地、0 在他处）
# 而 `<parameter idref="AS-to-OC">` 那种写法输出的是 `AS-to-OC1..AS-to-OC36`
# （展平的 6×6 常量矩阵，值恒为 0/1）—— 必须与计数区分开，否则会拿一堆 0/1
# 当成"迁移矩阵"，看起来还挺像。
_COL_BRACKET_ROUTE = re.compile(r'^(.+?)-to-(.+?)\[(\d+)\]$')
_COL_BRACKET_REWARD = re.compile(r'^(.+?)_reward\[(\d+)\]$')
_COL_BRACKET_TOTAL = re.compile(r'^(.+)\.count\[(\d+)\]$')
_COL_PLAIN_ROUTE = re.compile(r'^(.+?)-to-([A-Za-z][\w.]*?)$')
_COL_FLAT = re.compile(r'^(.+?)-to-([A-Za-z][\w.]*?)(\d+)$')


def _classify_columns(cols):
    """把日志列名分类成 ``route / reward / total / flat``。

    ``flat`` = 展平的静态矩阵列（`AS-to-OC1..36`）；它**不是计数**，
    分类出来是为了给调用方一个明确的拒绝理由，而不是悄悄按路线解析。
    """
    route, reward, total, flat = {}, {}, {}, {}
    for k, c in enumerate(cols):
        core = c[2:] if c.startswith('c_') else c
        m = _COL_BRACKET_ROUTE.match(core)
        if m:
            route[(m.group(1), m.group(2))] = k
            continue
        m = _COL_BRACKET_REWARD.match(core)
        if m:
            reward[m.group(1)] = k
            continue
        m = _COL_BRACKET_TOTAL.match(core)
        if m:
            total[m.group(1)] = k
            continue
        m = _COL_PLAIN_ROUTE.match(core)
        if m and not re.search(r'\d$', m.group(2)):
            route[(m.group(1), m.group(2))] = k
            continue
        m = _COL_FLAT.match(core)
        if m:
            flat.setdefault((m.group(1), m.group(2)), []).append(
                (int(m.group(3)), k))
    # 展平组如果凑不出"完整方阵"（组大小 = N²），多半是名字里本来就带数字
    # （如性状 `Line1`）→ 退回按普通路线解析，别把它当展平矩阵。
    # ⚠️ **恰好凑成 1×1 方阵判定**的单列（如 `X-to-Line1`）另行处理
    #    （2026-10-03 修复）：N=1 时方阵判定 `n²==len(items)` 恒成立 → 旧实现
    #    把它当"完整 1×1 展平矩阵" `continue` 掉，等于**静默丢弃**这一列。
    #    真正的展平矩阵至少 N²≥4 列，单列必然不是矩阵；但状态名可能真以数字
    #    结尾，故不擅自当路线，而是原样留在 flat 并登记到 `flat_single`，由调用方
    #    在结果里标注（不进矩阵）。**只动这一种**：单列且序号≠1 的组本来就走
    #    下面的回退路由（旧行为正确），不在此列、不受影响。
    flat_single = {}
    for (src, dst), items in list(flat.items()):
        if len(items) == 1:
            n = 1
            if sorted(i for i, _ in items) == list(range(1, n * n + 1)):
                flat_single[(src, dst)] = items[0][1]
                continue
        n = int(round(len(items) ** 0.5))
        if n * n == len(items) and sorted(i for i, _ in items) == list(range(1, n * n + 1)):
            continue
        for _, k in items:
            route[(src, dst + str(_))] = k
        flat.pop((src, dst))
    return {'route': route, 'reward': reward, 'total': total, 'flat': flat,
            'flat_single': flat_single}


# ---------------------------------------------------------------- 生成

def state_order_from_xml(xml_path):
    """读 BEAST 的**权威状态序号**：``<generalDataType><state code="…"/>``。

    这个顺序决定矩阵行列/奖励向量的下标含义 —— 实测 PVS 示例里是
    ``AS, EU, ME, NAm, OC, SAm``（字母序，Beauti 的写法），
    **不是** ``<attr name="Region">`` 在序列里首次出现的顺序。
    """
    with io.open(xml_path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    states = re.findall(r'<state\s+code="([^"]+)"', text)
    return [s.strip() for s in states if s.strip()]


def traits_from_xml(xml_path):
    """从 BEAST XML 里读离散性状的取值，并标明**顺序来源**。

    返回 ``{'trait', 'values', 'order_source', 'state_order', 'attr_order'}``

    ⚠️ **顺序必须用 `<state code=…>`（order_source='state_definition'）**，
    因为它就是 BEAST 内部的行列下标。``<attr name="…">`` 的首次出现顺序
    （order_source='attr_first_seen'）只是"序列在数据块里的排列"，
    与下标无关 —— 拿它去铺指示矩阵会让 `X-to-Y` **名不副实**
    （实测 VirPhyKit 的 PVS 示例就是这么错的：它按 attr 首现序
    `AS,OC,EU,ME,NAm,SAm` 铺矩阵，而 BEAST 的状态序是
    `AS,EU,ME,NAm,OC,SAm` → `AS-to-OC` 那个矩阵的 1 落在 (0,1)＝AS→EU）。
    """
    with io.open(xml_path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    attr_name, attr_vals = None, []
    for m in re.finditer(r'<attr\s+name="([^"]+)"\s*>([^<]*)</attr>', text):
        attr_name = m.group(1)
        v = m.group(2).strip()
        if v and v not in attr_vals:
            attr_vals.append(v)
    state_order = [s.strip() for s in re.findall(r'<state\s+code="([^"]+)"', text)
                   if s.strip()]
    if state_order:
        return {'trait': attr_name, 'values': state_order,
                'order_source': 'state_definition',
                'state_order': state_order, 'attr_order': attr_vals}
    return {'trait': attr_name, 'values': attr_vals,
            'order_source': 'attr_first_seen',
            'state_order': [], 'attr_order': attr_vals}


def build_matrices(traits):
    """生成 N×(N−1) 个 `<parameter id="A-to-B" value="{0/1 矩阵}"/>`。

    与 VirPhyKit `function_mmj.wrtcfg()` **逐字对齐**（含"值从下一行开始、
    每行空格分隔、结尾 `"/>`"这些格式细节），以便与参照实现逐字节对拍。
    """
    n = len(traits)
    zeros = ['0'] * n
    out = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            row = zeros[:]
            row[j] = '1'
            # ⚠️ id **不加宽度填充**（2026-10-03 修复）：旧写 `%2s-to-%2s`，单字符
            #    性状（如 'X'）会生成 `id=" X-to-AS"` 带前导空格的参数名 ——
            #    `audit_mjrm_xml` 的 declared/形状正则 `<parameter id="[A-Za-z]…`
            #    直接漏掉它，矩阵数、顺序、位置三项校验全部误判。除这一处外仍与
            #    VirPhyKit `wrtcfg()` 的格式对齐（值从下一行开始、每行空格分隔、
            #    结尾 `"/>`）。
            out.append('<parameter id="%s-to-%s" value="\n' % (traits[i], traits[j]))
            for _ in range(i):
                out.append(' '.join(zeros) + '\n')
            out.append(' '.join(row) + '\n')
            for _ in range(i + 1, n):
                out.append(' '.join(zeros) + '\n')
            out.append('"/>\n')
    return ''.join(out)


def build_rewards(traits, vykit_example_compat=False):
    """生成 `<rewards>` 块（每个性状一个 one-hot 向量）。

    默认与 VirPhyKit `function_mmj.wrt_rewards()` 的**代码**一致（4 空格缩进、
    `value=` 前 1 个空格、顺序跟随传入的 traits）。

    ``vykit_example_compat=True`` 则复刻它**示例文件** `PVS_with_matrix.xml` 的
    写法：`value=` 前**两个空格**、且按**字母序**排列。为什么要单独开这个开关 ——
    实测发现示例文件与它自己的代码不一致（示例是字母序 + 双空格，代码是跟随顺序 +
    单空格），说明那份示例是更早/手改过的产物。默认不迁就它，但保留开关以便
    "想跟示例逐字节对齐"的场合。
    """
    n = len(traits)
    order = sorted(traits) if vykit_example_compat else list(traits)
    gap = '  ' if vykit_example_compat else ' '
    # 示例文件里 rewards 那几行**顶格无缩进**（与它自己的代码也不同），
    # 所以"复刻示例"要连缩进一起去掉。
    ind = '' if vykit_example_compat else '    '
    out = ['<rewards>\n']
    for t in order:
        i = traits.index(t)
        vals = ['0.0'] * n
        vals[i] = '1.0'
        out.append('%s<parameter id="%s_reward"%svalue="%s" />\n'
                   % (ind, t, gap, ' '.join(vals)))
    out.append('</rewards>\n')
    return ''.join(out)


# likelihood 元素的两种方言：
#   · `markovJumpsTreeLikelihood`  —— VirPhyKit / Beauti（BEAST2 风格的祖先状态重建）
#   · `ancestralTreeLikelihood`    —— 参照管道的 `beast1_bridge.py` 从零生成的 BEAST1 写法
# 两者都能挂 `<parameter id="A-to-B">` + `<rewards>`，但**只有被 `<treeLikelihood idref>`
# 记录**时才会输出 `c_A-to-B[]` 计数（实测）。所以识别标签时必须把两种都认出来，
# 否则会误报"没有 likelihood 元素"。
#
# ── BEAST 1.10.4 的真实类体系（`javap` 自 beast.jar 读出，2026-09-16）──
#   XML 元素名由 **`<X>Parser` 去掉 Parser、首字母小写** 得到（所以 jar 里查不到字面量）。
#   能承载祖先状态重建/跳跃计数的likelihood 只有三档，且是**严格三级继承**：
#       BeagleTreeLikelihood                ← 纯序列似然，无性状
#         └─ AncestralStateBeagleTreeLikelihood   ＋TreeTraitProvider/AncestralStateTraitProvider
#              └─ MarkovJumpsBeagleTreeLikelihood ＋MarkovJumpsRegisterAcceptor/MarkovJumpsTraitProvider
#   对应 parser：BeagleTreeLikelihoodParser ← AncestralStateTreeLikelihoodParser
#              ← MarkovJumpsTreeLikelihoodParser
#   即 `markovJumpsTreeLikelihood` 是 `ancestralTreeLikelihood` 的**功能超集**；
#   而 ancestral 里**一个 Jump 相关方法都没有**（`javap` 实测），所以挂在它下面的
#   奖励矩阵/`<rewards>` 无处可用、被静默忽略。
#   另有 `markovJumpsLikelihoodLogger`（`MarkovJumpsLikelihoodLoggerParser`，按命名约定
#   派生自类名）是**专门的 jump 日志器**——官方还有这条接线途径，本模块走的是
#   "把 likelihood 放进 <log>" 那条（实测有效）。
_LIK_TAGS = (r'markovJumpsTreeLikelihood', r'ancestralTreeLikelihood')
# 不该挂 MJRM 参数的 likelihood（纯序列似然等）—— 只用于诊断，不做插入锚点
_OTHER_LIK_TAGS = (r'treeLikelihood', r'optimizedBeagleTreeLikelihood',
                   r'balancedBeagleTreeLikelihood',
                   r'approximateTreeLikelihood', r'alsTreeLikelihood')


def _find_likelihood(text):
    """找承载 MJRM 参数的 likelihood 元素 → ``(tag, id)`` 或 ``(None, None)``。"""
    for tag in _LIK_TAGS:
        m = re.search(r'<%s\b[^>]*\bid="([^"]+)"' % tag, text)
        if m:
            return tag, m.group(1)
    return None, None


def _insert_index(text):
    """该把 MJRM 块插到哪 —— **必须落在 likelihood 元素内部**。

    ⚠️ 这里踩过一个很隐蔽的坑（2026-09-16 实测）：模板里
    ``<!-- END Ancestral state reconstruction`` 这样的锚点**不止一处**
    （PVS 模板里有 5 处：1206 / 1263 / 1490 / 1541 / 1555 行），而
    ``<markovJumpsTreeLikelihood>`` 开标签在 1246 行 ——
    按"取第一个匹配行"插入（VirPhyKit `process_xml` 就是这么写的）会把矩阵插到
    **另一个 likelihood 里**。后果不会报错，而是 BEAST 跑完只吐出
    ``c_Region.count[1]`` **一列**、30 条路线全不见 —— 比报错难查得多。

    所以：先定位 likelihood 的开/闭标签，只在**其内部**取
    **最后一个**锚点（最靠近闭合标签的那个）。

    返回 ``(插入位置 或 None, (span_lo, span_hi) 或 None)``。
    """
    tag, _lid = _find_likelihood(text)
    if tag is None:
        return None, None
    open_m = re.search(r'<%s\b[^>]*>' % tag, text)
    close_m = re.search(r'</%s\s*>' % tag, text[open_m.end():])
    if not close_m:
        return None, None
    lo = open_m.end()
    hi = open_m.end() + close_m.start()
    inner = text[lo:hi]
    best = -1
    for mk in _TARGET_MARKERS:
        k = inner.rfind(mk)
        if k > best:
            best = k
    return (lo + best if best >= 0 else hi), (lo, hi)


def insert_into_template(template_path, traits=None, out_path=None,
                         trait_name=None, with_log=True, log_every=1000,
                         trait_order='state', fix_likelihood_tag=False):
    """把 MJRM 矩阵 + rewards 插进 BEAST XML 模板；顺带补上让计数输出的 log。

    ``trait_order``：矩阵的行列顺序。
      · ``'state'``（默认，**正确**）：按 XML 里 ``<generalDataType><state code=…>``
        的顺序 —— 那就是 BEAST 内部的行列下标。
      · ``'file'``：按 ``<attr>`` 首现顺序（＝VirPhyKit 的行为）。
        只在"要和它的旧产物逐字节对拍"时用；**会得到名不副实的 X-to-Y 标签**。
    ``traits=None`` 时自动探测（推荐）。传了 traits 会按目标顺序重排；
    集合与 XML 不一致时告警。

    ``fix_likelihood_tag=True``：若模板里已有的 MJRM 参数块挂在
    ``ancestralTreeLikelihood`` 下（参照管道 `beast1_bridge.py` 的产物就是这样），
    顺手把元素名改成 ``markovJumpsTreeLikelihood`` —— 实测原来的类**不支持**
    Markov jump，参数会被静默忽略、怎么跑都出不来计数。

    返回 dict：``{'ok', 'msg', 'out', 'traits', 'trait_name', 'n_matrices',
    'likelihood_id', 'likelihood_tag', 'logged', 'order_source', 'reordered',
    'reused_existing', 'warnings'}``
    """
    with io.open(template_path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    info = traits_from_xml(template_path)
    state_order = info['state_order']
    attr_order = info['attr_order']
    if trait_order not in ('state', 'file'):
        raise RuntimeError(f"trait_order 只支持 'state'/'file'，收到 {trait_order!r}")
    target = state_order if trait_order == 'state' else attr_order
    if not target:
        target = info['values']
    trait_name = trait_name or info['trait']
    warnings = []
    reordered = False

    if traits is None:
        traits = list(target)
    else:
        traits = [str(t).strip() for t in traits if str(t).strip()]
        if set(traits) != set(target) and target:
            warnings.append(
                f'传入的性状 {traits} 与 XML 里看到的 {target} 不是同一集合 —— '
                '矩阵可能错位')
        elif traits != target:
            # 按目标顺序重排（用户的输入顺序不决定下标，BEAST 的状态序才决定）
            traits = [t for t in target if t in traits]
            reordered = True
    if len(traits) < 2:
        return {'ok': False, 'msg': f'性状少于 2 个（拿到 {traits}）',
                'out': None, 'traits': traits, 'trait_name': trait_name,
                'n_matrices': 0, 'likelihood_id': None, 'logged': False,
                'order_source': info['order_source'], 'reordered': reordered,
                'warnings': warnings}
    if len(set(traits)) != len(traits):
        warnings.append('性状里有重复值 → 指示矩阵会重复/错位')
    bad = [t for t in traits if re.search(r'[^\w.\-]', t)]
    if bad:
        return {'ok': False, 'msg': f'性状名含非法字符（BEAST id 不允许）：{bad}',
                'out': None, 'traits': traits, 'trait_name': trait_name,
                'n_matrices': 0, 'likelihood_id': None, 'logged': False,
                'order_source': info['order_source'], 'reordered': reordered,
                'warnings': warnings}

    idx, span = _insert_index(text)
    if idx is None:
        return {'ok': False,
                'msg': '模板里找不到 `<markovJumpsTreeLikelihood id="…">` 及其闭合标签 '
                       '—— 请确认这份 XML 已经配好 Markov jumps（祖先状态重建）',
                'out': None, 'traits': traits, 'trait_name': trait_name,
                'n_matrices': 0, 'likelihood_id': None, 'logged': False,
                'order_source': info['order_source'], 'reordered': reordered,
                'warnings': warnings}
    # likelihood 的 id（两种方言都认：markovJumpsTreeLikelihood / ancestralTreeLikelihood）
    lik_tag, likelihood_id = _find_likelihood(text)

    block = build_matrices(traits) + build_rewards(traits)
    logged = False
    # 幂等保护：模板里**已经有** MJRM 块时不要再插一遍（否则 30 个矩阵翻倍、
    # BEAST 会因重复 id 报错）。参照管道 `beast1_bridge.py` 生成的 XML 就自带
    # 矩阵与 rewards，只是没把它们接进 log —— 那种情况只需要补 log。
    existing = re.findall(r'<parameter\s+id="\s*[A-Za-z][\w.]*-to-[A-Za-z][\w.]*"',
                          text)
    reused = bool(existing) and '<rewards>' in text
    out_text = text[:idx] + ('' if reused else block) + text[idx:]
    if reused:
        warnings.append(
            f'模板里已有 {len(existing)} 个 `A-to-B` 指示矩阵与 <rewards> '
            '→ 跳过矩阵插入，只按需补 jumpLog（避免重复 id）')
    # 一键修"类不对"：`ancestralTreeLikelihood` 不做 Markov jump（参数被静默忽略），
    # 改名成 `markovJumpsTreeLikelihood` 才生效。只动承载 id 的开/闭标签，idref 不动。
    fixed_tag = False
    if fix_likelihood_tag and reused and lik_tag == 'ancestralTreeLikelihood':
        pat_open = '<ancestralTreeLikelihood id="%s"' % likelihood_id
        if pat_open in out_text:
            out_text = out_text.replace(
                pat_open, '<markovJumpsTreeLikelihood id="%s"' % likelihood_id, 1)
            out_text = out_text.replace('</ancestralTreeLikelihood>',
                                        '</markovJumpsTreeLikelihood>', 1)
            fixed_tag = True
            lik_tag = 'markovJumpsTreeLikelihood'
            warnings.append(
                f'已把 `<ancestralTreeLikelihood id="{likelihood_id}">` 改名为 '
                '`<markovJumpsTreeLikelihood>` —— 实测原来的类不支持 Markov jump，'
                '改完计数才会真的产生')
    if fix_likelihood_tag and lik_tag == 'ancestralTreeLikelihood' and not reused:
        warnings.append('fix_likelihood_tag=True 但模板里没有现成的 MJRM 参数块，'
                        '未改名（先按普通生成流程处理）')
    if with_log and likelihood_id:
        # ⚠️ `<log>` **必须放在 `<mcmc>` 里面** —— 不能跟着矩阵一起插到锚点处，
        #    因为那个锚点（`<!-- END Ancestral state reconstruction` /
        #    `</markovJumpsTreeLikelihood>`）在 likelihood 元素**内部**。
        #    实测把它们放一起，BEAST1 直接报 XML 解析错误、连跑都跑不起来。
        #    官方写法（BEAST1 examples/TestXML/testMarkovJumps.xml）：
        #      <log logEvery="…" fileName="…"><treeLikelihood idref="…"/></log>
        m2 = re.search(r'</mcmc>', out_text)
        if not m2:
            warnings.append('XML 里没有 </mcmc>，无法把 jumpLog 放进 <mcmc> —— '
                            '跑完不会输出跳转计数')
        else:
            jumps = ('\t<log id="jumpLog" logEvery="%d" fileName="jumps.log">\n'
                     '\t\t<treeLikelihood idref="%s"/>\n'
                     '\t</log>\n' % (int(log_every), likelihood_id))
            out_text = out_text[:m2.start()] + jumps + out_text[m2.start():]
            logged = True
    elif with_log and not likelihood_id:
        warnings.append(
            '模板里没有 <markovJumpsTreeLikelihood id="…">，无法自动补 jumpLog；'
            '不加 log 的话 BEAST 跑完不会输出跳转计数（这正是 VirPhyKit 示例的情况）')
    if out_path is None:
        base, ext = os.path.splitext(template_path)
        out_path = f'{base}.mjrm{ext or ".xml"}'
    with io.open(out_path, 'w', encoding='utf-8', newline='') as f:
        f.write(out_text)
    msg = (f'已写入 {out_path}（{len(traits)} 个性状 → '
           f'{len(traits) * (len(traits) - 1)} 个指示矩阵'
           + ('，并补了 jumpLog' if logged else '，**未**补 jumpLog') + '）')
    if trait_order == 'file':
        warnings.append(
            "trait_order='file'：矩阵按 <attr> 首现顺序铺 —— 这是 VirPhyKit 的行为，"
            '当它与 <state code> 顺序不同时会得到**名不副实**的 X-to-Y 标签')
    elif state_order and attr_order and state_order != attr_order:
        warnings.append(
            f'矩阵顺序用了 BEAST 的状态序 {state_order}（正确）；'
            f'而 <attr> 首现序是 {attr_order} —— 两者不同，'
            'VirPhyKit 会用后者，这正是它示例错位的根因')
    return {'ok': True, 'msg': msg,
            'out': out_path, 'traits': traits, 'trait_name': trait_name,
            'n_matrices': len(traits) * (len(traits) - 1),
            'likelihood_id': likelihood_id, 'logged': logged,
            'likelihood_tag': lik_tag, 'fixed_likelihood_tag': fixed_tag,
            'order_source': info['order_source'], 'state_order': state_order,
            'attr_order': attr_order, 'reordered': reordered,
            'reused_existing': reused,
            'warnings': warnings}


# ---------------------------------------------------------------- 审计

def audit_mjrm_xml(xml_path):
    """**跑 BEAST 之前**检查这份 XML 到底能不能产出跳转计数。

    检查项（都是实测踩出来的）：
      1. 有没有 N×(N−1) 个 `A-to-B` 指示矩阵、形状是否是 N×N
      2. 有没有 `<rewards>`
      3. **这些参数有没有被任何 `<log>` 引用** —— 没有引用 = 跑完白跑
         （VirPhyKit 示例就是这种"哑"配置）
      4. log 里的引用方式是不是只会输出静态矩阵的那种
         （`<parameter idref="A-to-B">` 输出 `A-to-B1..N²` 的常量列，**不是计数**）

    返回 ``{'ok', 'traits', 'n_matrices', 'roots', 'logging': {...},
    'problems': [...], 'notes': [...]}``
    """
    with io.open(xml_path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    info = traits_from_xml(xml_path)
    traits = info['values']
    n = len(traits)
    # declared 正则带 `\s*`：旧版 VirPhyKit 风格产物用 `%2s` 填充，单字符性状会
    # 产出 `id=" X-to-AS"`；这里容忍前导空白（捕获组仍不含空格），新旧产物都能查。
    declared = re.findall(r'<parameter\s+id="\s*([A-Za-z][\w.]*-to-[A-Za-z][\w.]*)"',
                          text)
    has_rewards = '<rewards>' in text
    problems, notes = [], []

    # 形状检查：挨个看矩阵行数（同样容忍 id 前导空白）
    bad_shape = []
    for m in re.finditer(r'<parameter\s+id="\s*([A-Za-z][\w.]*-to-[A-Za-z][\w.]*)"'
                         r'\s+value="\s*\n(.*?)"\s*/>', text, re.S):
        rows = [r for r in m.group(2).splitlines() if r.strip()]
        if len(rows) != n or any(len(r.split()) != n for r in rows):
            bad_shape.append((m.group(1), len(rows)))
    if bad_shape:
        problems.append(f'{len(bad_shape)} 个矩阵形状与性状数 {n} 不符：'
                        f'{bad_shape[:3]}')

    # log 引用：**只看 likelihood 本身怎么被记录的** —— 光有 `<parameter idref="posterior">`
    # 这类引用完全不能说明 jump 计数会输出（VirPhyKit 示例就有 21 个参数引用，照样出不来计数）。
    lt_ids, alt_ids, param_refs = [], [], []
    for m in re.finditer(r'<log\b[^>]*>(.*?)</log>', text, re.S):
        seg = m.group(1)
        param_refs += re.findall(r'<parameter\s+idref="([^"]+)"', seg)
        lt_ids += re.findall(r'<(?<!ancestral)treeLikelihood\s+idref="([^"]+)"', seg)
        alt_ids += re.findall(r'<ancestralTreeLikelihood\s+idref="([^"]+)"', seg)
    flat_refs = [r for r in param_refs if _ROUTE_RE.match(r)]

    lik_tag, lik_id = _find_likelihood(text)
    if lik_id is None:
        # 退一步看看是不是挂到了"纯序列似然"上（`treeLikelihood` /
        # `optimizedBeagleTreeLikelihood` …）—— 那些类连性状都不管，更不做 jump
        other = None
        for tg in _OTHER_LIK_TAGS:
            m = re.search(r'<%s\b[^>]*\bid="([^"]+)"' % tg, text)
            if m:
                other = (tg, m.group(1))
                break
        if other and (declared or has_rewards):
            problems.append(
                f'MJRM 参数挂在了 `<{other[0]} id="{other[1]}">` 下 —— 这是**纯序列似然**'
                '（BeagleTreeLikelihood 系），既不管性状也不做 Markov jump，'
                '参数会被静默忽略。要挂到 `<markovJumpsTreeLikelihood>` 上'
                '（ancestral 系也不行：那个类没有 Jump 方法）')
        else:
            problems.append('模板里找不到承载 MJRM 参数的 likelihood 元素'
                            '（`<markovJumpsTreeLikelihood id="…">` 或 '
                            '`<ancestralTreeLikelihood id="…">`）—— '
                            '无法确认计数会输出（LSD2/RSPP 那类 likelihood 不含跳转计数）')
    elif lik_tag == 'ancestralTreeLikelihood' and declared and has_rewards:
        # ⚠️ 实测（BEAST 1.10.4，2026-09-16）：`ancestralTreeLikelihood` **不做** Markov
        #    jump —— 挂在它下面的 `A-to-B` 指示矩阵与 `<rewards>` 会被**静默忽略**。
        #    参照管道 `beast1_bridge.py` 从零生成的 XML 就是这种：真跑 2×10⁵ 步，
        #    主日志 84 列里 0 个 `-to-`、0 个 `_reward`；把元素名改成
        #    `markovJumpsTreeLikelihood`（+ 把 likelihood 接进 <log>）后，
        #    **37 列计数立刻出现**（30 路线 + count + 6 reward）。
        problems.append(
            f'承载 MJRM 参数的 likelihood 是 `<ancestralTreeLikelihood id="{lik_id}">` —— '
            '实测这个类**只做祖先状态重建、不做 Markov jump**，挂在它下面的 '
            '指示矩阵与 `<rewards>` 会被 BEAST **静默忽略**（真跑 2×10⁵ 步，主日志 '
            '84 列里 0 个计数列）。需要把元素名改成 `<markovJumpsTreeLikelihood '
            f'id="{lik_id}">`（并用 `insert_into_template(fix_likelihood_tag=True)` '
            '可一键改），再补 jumpLog')
    elif lik_id in lt_ids:
        pass                     # ✅ 唯一被实测证明能出计数的写法
    elif lik_id in alt_ids:
        # ⚠️ 实测（BEAST 1.10.4）：同一个 likelihood 用 `<ancestralTreeLikelihood idref>`
        #    记录时**只输出祖先状态概率**，不出 `c_A-to-B[]` 跳转计数。VirPhyKit 示例
        #    正是这种写法 —— 真跑了 2×10⁶ 步，主日志 125 列里 0 个计数列。
        problems.append(
            f'<log> 只通过 `<ancestralTreeLikelihood idref="{lik_id}">` 记录 —— '
            '实测这种写法只输出祖先状态概率，**不出** `c_A-to-B[]` 跳转计数'
            '（VirPhyKit 示例跑完 125 列里 0 个计数列就是这么来的）。'
            f'应改用 BEAST1 官方写法：`<treeLikelihood idref="{lik_id}"/>`')
    else:
        problems.append(
            f'没有任何 <log> 记录 likelihood `{lik_id}` —— 这份配置跑完不会输出'
            '跳转计数（指示矩阵被声明了也没用）。需要在 <mcmc> 里补一个 <log>，'
            '按 BEAST1 官方写法引用该 likelihood 对象')
    if flat_refs:
        problems.append(
            f'<log> 里直接引用了 {len(flat_refs)} 个 `A-to-B` 参数 —— 实测这只会输出'
            ' `A-to-B1..N²` 的**静态矩阵常量列**，不是跳转计数（值恒为 0/1）。'
            '应改为引用 likelihood 对象本身')

    if not declared:
        problems.append('没找到任何 `A-to-B` 指示矩阵')
    elif len(declared) != n * (n - 1):
        notes.append(f'指示矩阵 {len(declared)} 个，期望 N×(N−1)='
                     f'{n * (n - 1)} 个（N={n}）')
    if not has_rewards:
        problems.append('没有 <rewards> 块')

    # ── 位置检查：矩阵/rewards 必须在 `<markovJumpsTreeLikelihood>` **内部** ──
    _ltag, _lid = _find_likelihood(text)
    o_m = re.search(r'<%s\b[^>]*>' % _ltag, text) if _ltag else None
    c_m = (re.search(r'</%s\s*>' % _ltag, text[o_m.end():]) if o_m else None)
    if o_m and c_m:
        lo, hi = o_m.end(), o_m.end() + c_m.start()
        # 用正则而不是 str.find：旧版产物 id 可能带前导空白（`id=" X-to-AS"`）
        if declared:
            _pm = re.search(r'<parameter\s+id="\s*%s"' % re.escape(declared[0]),
                            text)
            first_p = _pm.start() if _pm else -1
        else:
            first_p = -1
        rw = text.find('<rewards>')
        if first_p >= 0 and not (lo < first_p < hi):
            problems.append(
                '**矩阵块不在 `<markovJumpsTreeLikelihood>` 内部**（在它之前或之后）'
                '—— 实测后果是 BEAST 照跑不误，但日志里**只有 '
                '`c_<trait>.count[]` 一列**、30 条路线全不见。模板里 '
                '`<!-- END Ancestral state reconstruction` 这类锚点往往**不止一处**，'
                '取"第一个匹配"就会插错元素；应取该 likelihood 内部的最后一个锚点')
        if rw >= 0 and not (lo < rw < hi):
            problems.append('`<rewards>` 不在 `<markovJumpsTreeLikelihood>` 内部')

    # ── 顺序检查（矩阵下标 vs BEAST 状态序）—— 实测最容易错、且错了看不出来 ──
    declared_order = []
    for d in declared:
        src = d.split('-to-', 1)[0].strip()
        if src not in declared_order:
            declared_order.append(src)
    state_order = [s.strip() for s in re.findall(r'<state\s+code="([^"]+)"', text)
                   if s.strip()]
    if state_order and declared_order and declared_order != state_order:
        problems.append(
            f'**矩阵下标与 BEAST 的状态序不一致**：矩阵按 {declared_order} 铺，'
            f'而 `<generalDataType><state code=…>` 的顺序是 {state_order}。'
            'BEAST 用后者当行列下标 → 每个 `X-to-Y` 的 1 都落到了别的状态对上，'
            '**标签名不副实**（例如标签 AS-to-OC 实际计的是 AS→EU）。'
            '生成时请用 trait_order="state"')
    # 内部自洽：rewards 的 one-hot 位置必须与矩阵顺序一致
    mismatch = []
    for m in re.finditer(r'<parameter\s+id="([\w.]+)_reward"[^>]*value="([^"]*)"',
                         text):
        name, vec = m.group(1), [x for x in m.group(2).split() if x]
        if len(vec) != n:
            continue
        try:
            one = vec.index('1.0')
        except ValueError:
            continue
        if declared_order and one < len(declared_order) \
                and declared_order[one] != name:
            mismatch.append(f'{name}_reward 的 1.0 在第 {one} 位 = '
                            f'矩阵顺序里的 {declared_order[one]}')
    if mismatch:
        problems.append(
            '**rewards 与矩阵顺序自相矛盾**：' + '；'.join(mismatch[:3])
            + '。同一个文件里两处必须用同一套下标，否则 reward（停留时间）'
              '与 jump 计数的标签都不可信')

    return {'ok': not problems, 'traits': traits, 'trait_name': info['trait'],
            'n_traits': n, 'n_matrices': len(declared), 'has_rewards': has_rewards,
            'log_refs': len(param_refs) + len(lt_ids) + len(alt_ids),
            'likelihood_id': lik_id, 'likelihood_tag': lik_tag,
            'route_refs': flat_refs[:5],
            'treeLikelihood_refs': lt_ids[:5], 'ancestral_refs': alt_ids[:5],
            'declared_order': declared_order, 'state_order': state_order,
            'problems': problems, 'notes': notes}


def make_dating_config(src_xml, out_xml=None, strip_logs=True, prefix=None):
    """从"含离散性状（系统地理）"的 XML 生成**只靠序列似然**的定年版。

    为什么需要它：定年与迁移路线常常**分开跑**。定年那一份应当是
    「序列似然 + 时钟 + 树先验」，**不该把位置性状模型也算进去** ——
    否则树会被迁移模型影响，而迁移那份又单独跑一次，两次的树不一致。

    做法（Beauti 生成的 XML 自带锚点，所以很安全）：
      · `<likelihood>` 与 `<prior>` 里被
        `<!-- START/END Discrete Traits Model -->` 包住的那些行 → 注释掉
      · `<log>` / `<logTree>` 里引用该性状（`<前缀>.xxx` / 性状 likelihood）的行 → 注释掉
        （否则会记录只由先验驱动的、无意义的数值）
      · **其它一律不动**：树先验（skyline/coalescent）、时钟（ucln/strict）、
        序列 site model、operators 全部保留

    ⚠️ 位置模型的 operators 仍然留着 → 那些参数会在先验下随机游走，
    对定年**无影响**（它们已不在 posterior 里），只是多花一点算力。
    """
    with io.open(src_xml, encoding='utf-8', errors='replace') as f:
        text = f.read()
    # 性状前缀（Beauti 用 `<generalDataType id="Region.dataType">`）
    if prefix is None:
        m = re.search(r'<generalDataType\s+id="([^."]+)\.dataType"', text)
        if m:
            prefix = m.group(1)
        else:
            m = re.search(r'<markovJumpsTreeLikelihood\s+id="([^.]+)\.', text)
            prefix = m.group(1) if m else None
    lik_id = None
    m = re.search(r'<markovJumpsTreeLikelihood\s+id="([^"]+)"', text)
    if m:
        lik_id = m.group(1)

    lines = text.splitlines(keepends=True)
    out = []
    n_lik = n_pri = n_log = 0
    section = None
    in_dt = False            # 是否处于 START/END Discrete Traits Model 之间
    for ln in lines:
        s = ln.strip()
        if s.startswith('<likelihood'):
            section = 'likelihood'
        elif s.startswith('<prior'):
            section = 'prior'
        elif s.startswith('<log') or s.startswith('<logTree'):
            section = 'log'
        elif s.startswith('</likelihood>') or s.startswith('</prior>') \
                or s.startswith('</log>'):
            section = None
        if s.startswith(_DT_START):
            in_dt = True
        if in_dt and section in ('likelihood', 'prior') and s and not s.startswith('<!'):
            out.append('<!-- MJRM-dating-strip --> ' + ln)
            if section == 'likelihood':
                n_lik += 1
            else:
                n_pri += 1
            continue
        if s.startswith(_DT_END):
            in_dt = False
        if strip_logs and section == 'log' and prefix:
            if (f'{prefix}.' in ln) or (lik_id and lik_id in ln) \
                    or (lik_id and f'idref="{lik_id}"' in ln):
                out.append('<!-- MJRM-dating-strip --> ' + ln)
                n_log += 1
                continue
        out.append(ln)

    new_text = ''.join(out)
    if out_xml is None:
        base, ext = os.path.splitext(src_xml)
        out_xml = f'{base}.dating{ext or ".xml"}'
    with io.open(out_xml, 'w', encoding='utf-8', newline='') as f:
        f.write(new_text)
    warnings = []
    # 诚实交代"没被清掉的部分"：`<prior>` 里还引用了该性状的参数（Beauti 把它们写在
    # Discrete Traits 标记**之外**）。它们只是给位置模型自己的先验项 ——
    # 不影响树/时钟/序列参数，但会让那几个参数在先验下空转（多几个无意义的采样列）。
    if prefix:
        pri = re.search(r'<prior\b[^>]*>([\s\S]*?)</prior>', new_text)
        resid = 0
        if pri:
            resid = sum(1 for ln in pri.group(1).splitlines()
                        if f'{prefix}.' in ln
                        and not ln.strip().startswith('<!--'))
        if resid:
            warnings.append(
                f'`<prior>` 里仍有 {resid} 行引用 `{prefix}.` 参数（Beauti 把它们放在'
                ' Discrete Traits 标记之外）—— 对定年参数**无影响**，只是那几个位置参数'
                '会在先验下空转；要彻底干净可手动删掉它们（本函数不动，免得误删树/时钟先验）')
    if n_lik == 0:
        warnings.append('`<likelihood>` 里没找到 Discrete Traits Model 段 —— '
                        '这份 XML 可能本来就不含性状（那就直接当定年版用）')
    return {'ok': True, 'out': out_xml, 'prefix': prefix, 'likelihood_id': lik_id,
            'n_likelihood': n_lik, 'n_prior': n_pri, 'n_log': n_log,
            'warnings': warnings}


def build_pair(template_path, out_dir=None, base_name=None, traits=None,
               trait_order='state', fix_likelihood_tag=True, log_every=1000):
    """一次生成**两份**配置（定年 / 迁移），因为这两件事常常分开跑。

    * ``<base>.routes.xml`` —— 迁移路线用：序列似然 **＋** `markovJumpsTreeLikelihood`
      （位置性状 + 跳跃计数）**＋** 让计数输出的 jumpLog。两者是**并列的两个
      likelihood**，不是替代关系（`markovJumpsTreeLikelihood` 继承自
      `ancestralTreeLikelihood`，但它只是"额外挂上去的第二份似然"）。
    * ``<base>.dating.xml`` —— 定年用：把 Discrete Traits Model 段注释掉，
      只留「序列似然 + 时钟 + 树先验」。

    返回 ``{'routes': {...}, 'dating': {...}, 'ok': bool}``。
    """
    out_dir = out_dir or os.path.dirname(os.path.abspath(template_path))
    os.makedirs(out_dir, exist_ok=True)
    base = base_name or os.path.splitext(os.path.basename(template_path))[0]
    routes_xml = os.path.join(out_dir, f'{base}.routes.xml')
    dating_xml = os.path.join(out_dir, f'{base}.dating.xml')
    r = insert_into_template(template_path, traits=traits, out_path=routes_xml,
                             trait_order=trait_order, with_log=True,
                             log_every=log_every,
                             fix_likelihood_tag=fix_likelihood_tag)
    d = make_dating_config(r['out'] if r.get('ok') else template_path,
                           out_xml=dating_xml)
    return {'ok': bool(r.get('ok')) and bool(d.get('ok')),
            'routes': r, 'dating': d,
            'routes_xml': routes_xml, 'dating_xml': dating_xml}


# ---------------------------------------------------------------- 解析

def parse_beast_log(log_path, burnin=0.1):
    """读 BEAST 制表符日志 → ``(columns, rows, n_skipped)``。

    BEAST1 的日志：若干 `#` 注释行 + 一行以 `state` 开头的**列名**行 + 数据行。
    烧弃按行数比例裁掉末尾 ``burnin`` 之外的开头部分。
    """
    with io.open(log_path, encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()
    hdr = None
    first_data = None
    for i, l in enumerate(lines):
        if l.startswith('state\t') or l.startswith('state '):
            hdr = l
            first_data = i + 1
            break
    if hdr is None:
        raise RuntimeError(f'{os.path.basename(log_path)} 里找不到以 state 开头的列名行'
                           '（不是 BEAST tab 日志？）')
    cols = [c.strip() for c in hdr.split('\t')]
    rows = []
    for l in lines[first_data:]:
        if not l.strip():
            continue
        parts = l.split('\t')
        if len(parts) < len(cols):
            parts = l.split()
            if len(parts) < len(cols):
                continue
        rows.append(parts)
    keep = int(len(rows) * (1.0 - float(burnin)))
    skipped = len(rows) - keep
    return cols, rows[max(0, skipped):], skipped


def jump_matrix_from_log(log_path, burnin=0.1, traits=None):
    """BEAST 日志 → 迁移（跳转）矩阵 + 停留时间，含逐路线均值 / 中位数 / 95% 区间。

    返回 ``{'routes', 'matrix', 'rewards', 'dwell', 'traits', 'n_samples',
    'burnin', 'total_mean', 'log', 'warnings'}``

    * ``routes``：``{'EU_to_OC': {'mean','median','lo95','hi95','n'}, …}``
      —— BEAST1 给的是**条件期望跳转次数**（在采样到的祖先状态/速率下），
      所以是小数而不是整数；这不是 bug。
    * ``rewards`` / ``dwell``：各地点**停留时间**（reward 定义＝在本地取 1.0），
      归一化后即"各区划占了多少进化时间"—— 判断 source/sink 的另一条独立证据。
    * ``total_mean``：优先取日志里的 `c_<trait>.count`（实测＝所有有序对的跳转之和，
      不是"事件数"），没有该列时才退化成逐路线求和。

    ⚠️ **两条防呆守卫**（都是实测踩出来的）：
      1. 只找到 `A-to-B1..N²` 展平列 → 直接报错（那是静态指示矩阵，不是计数）。
      2. 既没有 route 列也没有 reward 列 → 报错并列出可用列名。
    """
    cols, rows, skipped = parse_beast_log(log_path, burnin)
    warnings = []
    cls = _classify_columns(cols)
    route, reward, total, flat = (cls['route'], cls['reward'],
                                  cls['total'], cls['flat'])
    # 1×1 展平组单列（如 `X-to-Line1`）：不进矩阵也不静默丢 —— 在结果里标注
    flat_single = cls.get('flat_single') or {}
    if flat_single:
        warnings.append(
            f'有 {len(flat_single)} 个单列无法判定为 `A-to-B` 计数'
            f'（登记为 1×1 展平组、未参与矩阵统计）：'
            + ', '.join('%s-to-%s' % k for k in list(flat_single)[:5])
            + ' —— 状态名可能以数字结尾（如 Line1）；'
            '请核对 XML 的 <state code> 是否真是这个名字，再决定是否按路线统计')
    # 只有"真正的多列展平组"才按"log 了静态指示矩阵"报错（1×1 单列不算）
    flat_multi = {k: v for k, v in flat.items() if k not in flat_single}
    if not route and not reward and flat_multi:
        stems = list(flat_multi)[:2]
        raise RuntimeError(
            f'{os.path.basename(log_path)} 里只有 "A-to-B1..N²" 这种展平矩阵列'
            f'（如 {stems}）—— 这说明 XML 里 log 的是**静态指示矩阵本身**，'
            '不是跳转计数。请按 BEAST1 官方写法把 likelihood 对象放进 <log>'
            '（见 mjrm.insert_into_template(with_log=True)）；这个日志给不出迁移矩阵。')
    if not route and not reward:
        # 1×1 展平组的提示要跟着报错走：否则标注只进了 warnings，而这条路径
        # 直接 raise，用户看不到"列为什么没被统计"（2026-10-03）。
        hint = ''
        if flat_single:
            hint = ('；另有 %d 个单列被登记为 1×1 展平组（如 %s）——'
                    '状态名可能以数字结尾，请核对 XML 的 <state code>'
                    % (len(flat_single),
                       ', '.join('%s-to-%s' % k for k in list(flat_single)[:3])))
        raise RuntimeError(
            f'{os.path.basename(log_path)} 里既没有 `A-to-B` 计数列、也没有 '
            f'`_reward` 列。可用列（前 20）：{cols[:20]}' + hint)

    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return -1.0        # BEAST 用 '-' 表示未采样

    import numpy as np

    def _collect(idx_map):
        out = {}
        for key, k in idx_map.items():
            vals = []
            for r in rows:
                if k >= len(r):
                    continue
                v = _num(r[k])
                if v >= 0.0:
                    vals.append(v)
            out[key] = vals
        return out

    def _stats(vals):
        return {'n': len(vals), 'mean': round(float(np.mean(vals)), 4),
                'median': round(float(np.median(vals)), 4),
                'lo95': round(float(np.percentile(vals, 2.5)), 4),
                'hi95': round(float(np.percentile(vals, 97.5)), 4)}

    routes, matrix = OrderedDict(), OrderedDict()
    route_vals = _collect(route)
    for (src, dst), vals in route_vals.items():
        if len(vals) < 5:
            warnings.append(f'{src}-to-{dst} 有效采样仅 {len(vals)} 个，已跳过')
            continue
        st = _stats(vals)
        routes[f'{src}_to_{dst}'] = st
        matrix.setdefault(src, OrderedDict())[dst] = st['mean']

    rewards, dwell = OrderedDict(), OrderedDict()
    reward_vals = _collect(reward)
    for name, vals in reward_vals.items():
        if len(vals) < 5:
            continue
        rewards[name] = _stats(vals)
    tot_dwell = sum(v['mean'] for v in rewards.values())
    if tot_dwell > 0:
        for name, st in rewards.items():
            dwell[name] = round(st['mean'] / tot_dwell, 4)

    total_mean = None
    if total:
        name = next(iter(total))
        vals = _collect(total)[name]
        if len(vals) >= 5:
            total_mean = round(float(np.mean(vals)), 4)
    if total_mean is None:
        total_mean = round(sum(v['mean'] for v in routes.values()), 4)
        warnings.append('日志里没有 `c_<trait>.count` 列，总跳转数由逐路线求和得到')

    names = list(traits) if traits else sorted({s for s, _ in route}
                                               | {d for _, d in route})
    if skipped:
        warnings.append(f'按 burnin={burnin} 丢掉前 {skipped} 个采样点')
    if len(rows) < 200:
        warnings.append(f'有效采样点只有 {len(rows)} 个 —— 后验区间偏粗，'
                        '建议加长链或减小 logEvery')
    return {'routes': routes, 'matrix': matrix, 'rewards': rewards, 'dwell': dwell,
            'traits': names, 'n_samples': len(rows), 'burnin': burnin,
            'total_mean': total_mean, 'warnings': warnings,
            'log': os.path.basename(log_path)}


def write_jump_csvs(result, out_dir, prefix='mjrm'):
    """把跳转矩阵写成三张表（与 MOT 一致的落盘风格）。返回产物路径列表。"""
    os.makedirs(out_dir, exist_ok=True)
    outs = []
    p1 = os.path.join(out_dir, f'{prefix}_routes.tsv')
    with io.open(p1, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(['route', 'mean', 'median', 'lo95', 'hi95', 'n'])
                + '\n')
        for k, v in result['routes'].items():
            f.write('\t'.join([k, str(v['mean']), str(v['median']),
                               str(v['lo95']), str(v['hi95']), str(v['n'])])
                    + '\n')
    outs.append(p1)
    names = list(result['traits'])
    p2 = os.path.join(out_dir, f'{prefix}_matrix.tsv')
    with io.open(p2, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(['from\\to'] + names) + '\n')
        for a in names:
            row = [str(result['matrix'].get(a, {}).get(b, 0.0)) for b in names]
            f.write('\t'.join([a] + row) + '\n')
    outs.append(p2)
    if result.get('dwell'):
        p3 = os.path.join(out_dir, f'{prefix}_dwell.tsv')
        with io.open(p3, 'w', encoding='utf-8', newline='') as f:
            f.write('\t'.join(['region', 'time_mean', 'fraction']) + '\n')
            for k in names:
                if k in result['rewards']:
                    f.write('\t'.join([k, str(result['rewards'][k]['mean']),
                                       str(result['dwell'].get(k, ''))]) + '\n')
        outs.append(p3)
    return outs

