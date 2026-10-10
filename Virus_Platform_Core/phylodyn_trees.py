# -*- coding: utf-8 -*-
"""进化动力学·树类引擎（MCC 注释树解析 / TempMig / RRT / RSPP / BSP / TreeTime / LTT）。

与 VirPhyKit（Yin et al., 2025, Ecol Evol）方法学对齐、代码全部自研（GPL-3.0，
不复制其代码）。关键口径注明出处：

* MCC 注释解析：BEAST1 `[&K=v,...]` / BEAST2 `[&...]` 方括号注释，NEXUS 包装，
  根状态取 `K.set`+`K.set.prob`（或 BEAST2 的 `max.set`/`max`）。
* TempMig 迁移矩阵：对齐 VirPhyKit TempMig 底层（Brynildsrud global L4 MTB 脚本
  的分支-年计数法）：每条分支（父状态→子状态）对它跨越的**每个日历年** +1；
  年轴 = rootyear = `第二年新近叶年 - ceil(根高)`（复刻其口径，见下）。
* RRT（MCC 版）：真数据根状态 vs 区域随机化副本根状态；判据 = 真数据
  最大后验状态的根后验 > 所有随机副本中该状态的最大根后验。
* RSPP-Viz：单棵 MCC 树的根状态后验 → 条形/饼图数据。
* BSP-Viz：Tracer 式 TSV（Time/Median/Lower/Upper）为精确口径；BEAST .log
  （popSizes/groupSizes/rootHeight）为**近似口径**（组界按等区间数近似，
  前端明确标注）。
* TreeTime-RTT：本机装了 treetime 则走真引擎（最小二乘重根 + 分子钟），
  另补 AUDIT 记缺的 date-permutation 检验；treetime 不可用时退回自研回归。
* TreeDater-LTT：R/treedater 本机不可用，自研 LTT（谱系-时间曲线 +
  log-线性净多样化斜率）；可选先走 LSD2 定年（3rd/tools/lsd2，见 phylogeo）。
"""
import bisect
import io
import math
import os
import re

from Virus_Platform_Core.phylodyn_kit import (
    read_text, date_to_decimal, validate_date_str,
)

# ================================================================
# BEAST 注释 newick / NEXUS 解析
# ================================================================


def _strip_nexus(tree_text):
    """NEXUS/裸 newick → 纯树串。

    终点取**最后一个 ) 之后的注释放完整**（根节点的 [&...] 在最外层 ) 之后，
    到首个 ; 为止）——直接截到 ) 会把根注释切掉，TempMig/RSPP 就丢了根状态。
    """
    if '(' not in tree_text:
        raise ValueError('文件里没有找到树（缺少括号结构）')
    start = tree_text.index('(')
    end = tree_text.rindex(')')
    if end < start:
        raise ValueError('树串括号不配对')
    j = end + 1
    while j < len(tree_text) and tree_text[j] != ';':
        if tree_text[j] == '\n' and ')' not in tree_text[j:j + 200] \
                and 'End' in tree_text[j:j + 20]:
            break
        j += 1
    return tree_text[start:j].strip()


def _split_pairs(s):
    """按不在引号/花括号内的逗号切 key=value 片段。"""
    out, buf, depth, q = [], [], 0, None
    for ch in s:
        if q:
            buf.append(ch)
            if ch == q:
                q = None
            continue
        if ch in '"\'':
            q = ch
            buf.append(ch)
        elif ch == '{':
            depth += 1
            buf.append(ch)
        elif ch == '}':
            depth -= 1
            buf.append(ch)
        elif ch == ',' and depth == 0:
            out.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append(''.join(buf))
    return [x for x in (y.strip() for y in out) if x]


def _parse_ann(s):
    """`&k=v,k={a,b},k2="x"` → dict。值：字符串 / 列表 / 裸 token。"""
    ann = {}
    if s.startswith('&'):
        s = s[1:]
    for pair in _split_pairs(s):
        if '=' not in pair:
            continue
        k, v = pair.split('=', 1)
        k = k.strip()
        v = v.strip()
        if v.startswith('{') and v.endswith('}'):
            inner = v[1:-1]
            vals = [x.strip().strip('"').strip("'") for x in inner.split(',') if x.strip()]
            ann[k] = vals
        elif v.startswith('"') and v.endswith('"') and len(v) >= 2:
            ann[k] = v[1:-1]
        elif v.startswith("'") and v.endswith("'") and len(v) >= 2:
            ann[k] = v[1:-1]
        else:
            ann[k] = v
    return ann


def parse_beast_newick(text):
    """带 [&...] 注释的 newick → 根节点 dict {name,length,anns,children,cumdist}。"""
    s = _strip_nexus(text)
    pos = [0]

    def _anns():
        ann = {}
        while pos[0] < len(s) and s[pos[0]] == '[':
            depth, j, q = 0, pos[0], None
            while j < len(s):
                ch = s[j]
                if q:
                    if ch == q:
                        q = None
                elif ch in '"\'':
                    q = ch
                elif ch == '[':
                    depth += 1
                elif ch == ']':
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            ann.update(_parse_ann(s[pos[0] + 1:j]))
            pos[0] = j + 1
        return ann

    def _name():
        # 前导空白先跳过（旧实现靠末尾 strip，效果相同，但引号判定要在跳过之后）
        while pos[0] < len(s) and s[pos[0]] in ' \t':
            pos[0] += 1
        if pos[0] < len(s) and s[pos[0]] == "'":
            # Newick 引号名：`'…'`，内部 `''` = 一个字面单引号。
            # write_beast_newick 的 `_name_str` 会这样写含 `(),:;[]`/空白的名字
            # （2026-10-03）；不认引号的话 `'C,three'` 会被逗号切成两个节点。
            pos[0] += 1
            buf = ''
            while pos[0] < len(s):
                ch = s[pos[0]]
                pos[0] += 1
                if ch == "'":
                    if pos[0] < len(s) and s[pos[0]] == "'":   # '' → 字面单引号
                        buf += "'"
                        pos[0] += 1
                        continue
                    break                                       # 收尾引号
                buf += ch
            return buf
        buf = ''
        while pos[0] < len(s) and s[pos[0]] not in ',():[':
            buf += s[pos[0]]
            pos[0] += 1
        return buf.strip()

    def _len():
        length = 0.0
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            buf = ''
            while pos[0] < len(s) and s[pos[0]] not in ',()':
                buf += s[pos[0]]
                pos[0] += 1
            try:
                length = float(buf)
            except ValueError:
                length = 0.0
        return length

    def _node():
        node = {'name': '', 'length': 0.0, 'anns': {}, 'children': []}
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
                raise ValueError('树串在位置 %d 处意外终止' % pos[0])
            node['anns'] = _anns()
            node['name'] = _name()
            node['anns'].update(_anns())
            node['length'] = _len()
        else:
            node['name'] = _name()
            node['anns'] = _anns()
            node['length'] = _len()
        return node

    return _node()


def _walk(root):
    """先序遍历，带父指针与累计根距。yield (node, parent, cumdist)。"""
    stack = [(root, None, 0.0)]
    while stack:
        n, p, d = stack.pop()
        n['_parent'] = p
        yield n, p, d
        for c in n['children']:
            stack.append((c, n, d + (c['length'] or 0.0)))


def tree_height(root):
    """根高 = 根到最深叶的累计枝长。"""
    return max(d for _n, _p, d in _walk(root))


def node_ages(root, T=None):
    """每节点 age（距今）。优先用注释 height（距现在的年龄），否则 T-累计根距。

    兼容两种节点结构：本模块 parse_beast_newick（带 anns）与 phylogeo.parse_newick
    （裸 dict，无 anns）。
    """
    T = T if T is not None else tree_height(root)
    ages = {}
    for n, _p, d in _walk(root):
        h = (n.get('anns') or {}).get('height')
        try:
            ages[id(n)] = float(h)
        except (TypeError, ValueError):
            ages[id(n)] = T - d
    return ages, T


# ---------------------------------------------------------------- 状态提取

_ANN_NOISE = re.compile(r'\.(prob|set|set\.prob|rate|rate\.median|rate\.range|rate_95%_HPD|'
                        r'median|range|_95%_HPD|height|length|posterior|max\.prob)$')


def detect_trait(root):
    """自动探测离散性状键名：优先 max（BEAST2 MAP），其次 K.set / K+K.prob 成对。"""
    keys = set()
    for n, _p, _d in _walk(root):
        keys.update(n['anns'].keys())
    if 'max' in keys:
        return 'max'
    set_keys = {k[:-4] for k in keys if k.endswith('.set')}
    if set_keys:
        # 取在根注释里出现 set.set.prob 的那个优先
        root_keys = {k[:-4] for k in root['anns'] if k.endswith('.set')}
        for k in sorted(root_keys):
            return k
        return sorted(set_keys)[0]
    pair = {k for k in keys if '.' not in k and (k + '.prob') in keys}
    if pair:
        return sorted(pair)[0]
    return None


def node_state(n, trait):
    """单节点 MAP 状态：trait 直读 > set/set.prob argmax > max > 命名首 token 兜底。"""
    anns = n['anns']
    if trait:
        v = anns.get(trait)
        if isinstance(v, str) and v:
            return v
        s, p = anns.get(trait + '.set'), anns.get(trait + '.set.prob')
        if isinstance(s, list) and isinstance(p, list) and s and len(s) == len(p):
            try:
                return s[max(range(len(p)), key=lambda i: float(p[i]))]
            except (TypeError, ValueError):
                pass
        if isinstance(s, list) and s and (p is None):
            return s[0]
    if 'max' in anns and isinstance(anns['max'], str) and anns['max']:
        return anns['max']
    # 兜底：叶名首 token（VirPhyKit 示例的命名约定 CCD_AB622861_2008.58 -> CCD）
    if not n['children'] and n['name']:
        return n['name'].split('_')[0]
    return 'Unknown'


def root_state_posterior(root, trait_hint=None):
    """根注释 → [(state, prob)]（BEAST1 set/set.prob、BEAST2 max.set/max、单值三口径）。"""
    anns = root['anns']
    for k in ([trait_hint] if trait_hint else []):
        s, p = anns.get(k + '.set'), anns.get(k + '.set.prob')
        if isinstance(s, list) and isinstance(p, list):
            return _zip_sp(s, p)
        v = anns.get(k)
        if isinstance(v, str) and v:
            prob = anns.get(k + '.prob')
            try:
                return [(v, float(prob))]
            except (TypeError, ValueError):
                return [(v, 1.0)]
    if isinstance(anns.get('max.set'), list):
        return _zip_sp(anns['max.set'], anns.get('max.set.prob'))
    set_keys = sorted(k[:-4] for k in anns if k.endswith('.set'))
    for k in set_keys:
        s, p = anns.get(k + '.set'), anns.get(k + '.set.prob')
        if isinstance(s, list) and isinstance(p, list):
            return _zip_sp(s, p)
    singles = []
    for k, v in anns.items():
        if isinstance(v, str) and (k + '.prob') in anns and not _ANN_NOISE.search(k + '.'):
            try:
                singles.append((k, v, float(anns[k + '.prob'])))
            except (TypeError, ValueError):
                continue
    if singles:
        k, v, p = singles[0]
        return [(v, p)]
    raise ValueError('根节点没有可识别的状态后验注释'
                     '（需要 K.set/K.set.prob、max.set 或 K/K.prob）')


def _zip_sp(states, probs):
    """`set` 与 `set.prob` 两个后验列表 → [(状态, 概率), …]。

    ⚠️ 2026-10-03 补长度校验：原实现直接 `zip` —— 两个注解长度不一致时
    **静默按短的截断**（`K.set` 4 个状态、`K.set.prob` 3 个值 → 少一个状态，
    后验和就不是 1 了，还看不出来）。长度不一致直接报错，不猜。
    """
    states = list(states or [])
    probs = list(probs or [])
    if len(states) != len(probs):
        raise ValueError('set 与 set.prob 个数不一致（%d vs %d）—— '
                         '截断会把根状态后验静默算错，拒绝'
                         % (len(states), len(probs)))
    try:
        pairs = [(s, float(p)) for s, p in zip(states, probs)]
    except (TypeError, ValueError):
        raise ValueError('set 与 set.prob 数值无法解析')
    if not pairs:
        raise ValueError('根状态后验为空')
    return pairs


def _parse_translate(tree_text):
    """NEXUS Translate 块 → {数字: 真实叶名}（无块时 {}）。"""
    m = re.search(r'\btranslate\s+(.*?);', tree_text, re.S | re.I)
    if not m:
        return {}
    out = {}
    # 分隔符 [,;]? 必须可选：上面的截块 (.*?) 在块尾 ; 处停止，最后一条
    # 翻译项（"N name" 后直接跟块尾 ;）在 group(1) 里没有尾随分隔符，
    # 强制 [,;] 会把最后一条叶名映射静默丢掉。
    for mm in re.finditer(r'([A-Za-z0-9_.]+)\s+("[^"]*"|\'[^\']*\'|[^\s,;]+)\s*[,;]?',
                          m.group(1)):
        name = mm.group(2).strip().strip('"').strip("'")
        if name:
            out[mm.group(1)] = name
    return out


def load_mcc_tree(path):
    """读 MCC 树文件（NEXUS/newick，支持 Translate 块）→ (root, T)。"""
    text = read_text(path)
    trans = _parse_translate(text)
    root = parse_beast_newick(text)
    if trans:
        for n, _p, _d in _walk(root):
            if not n['children'] and n['name'] in trans:
                n['name'] = trans[n['name']]
    return root, tree_height(root)


# ================================================================
# TempMig —— 分支-年迁移矩阵（对齐 VirPhyKit TempMig 口径）
# ================================================================

# 叶名末段数字"像不像采样年"的两道闸（2026-10-03，见 `_tip_years`）：
#   ① 年份区间：病毒/类病毒的采样年都落在 1000–3000 内，之外必是编号/版本号；
#   ② 跨度上限：真实采样跨度不会超过 500 年（一百年都算长的了），
#      超过说明末段数字里混进了非年份 → 宁可报错让用户改写，不许算错年轴。
_TIP_YEAR_MIN, _TIP_YEAR_MAX = 1000.0, 3000.0
_TIP_YEAR_SPAN_MAX = 500.0


def _tip_years(root):
    """叶名末段 `_2008.58197` → 年；返回升序年列表（只收**可信**的采样年）。

    ⚠️ 2026-10-03 修两处：
      ① 原实现 `float(tail)` 一成功就当年份收下 —— `sample_1`、`|2`、`_3`
         这类编号也会进年列表，于是 `mostcurrent = distinct[-2]`（"第二大年"，
         VirPhyKit 口径）会被编号顶掉、日历整轴错位。现在候选必须落在
         [_TIP_YEAR_MIN, _TIP_YEAR_MAX)（先试整段、再试前 8 字符，与旧口径
         的两次尝试一致），且所有候选的**跨度**不超过 `_TIP_YEAR_SPAN_MAX`；
      ② 一个可信年份都没有时**抛中文错误**。原实现返回空表 → 调用方
         `distinct[-2]` 拿不到 → `mostcurrent = 0` → `rootyear = -ceil(T)`
         凭空产出负年轴（静默错）。
    """
    years = []
    for n, _p, _d in _walk(root):
        if n['children'] or not n['name']:
            continue
        tail = n['name'].rsplit('_', 1)[-1]
        for probe in (tail, tail[:8]):
            try:
                f = float(probe)
            except ValueError:
                continue
            if _TIP_YEAR_MIN <= f < _TIP_YEAR_MAX:
                years.append(f)
            break
    years.sort()
    if not years:
        raise ValueError(
            '树上一个可信的采样年都没有 —— 叶名末段需形如 `_2008.58197`'
            '（年份在 1000–3000 内）；采样年请写进 tip 名或改用元数据提供')
    if years[-1] - years[0] > _TIP_YEAR_SPAN_MAX:
        raise ValueError(
            '叶名末段的"年份"跨度 %.1f 年（%.1f ~ %.1f）不像真实采样年 ——'
            '多半混进了编号/版本号；采样年请写进 tip 名或改用元数据提供'
            % (years[-1] - years[0], years[0], years[-1]))
    return years


def migration_matrix_from_mcc(tree_path, out_dir=None, trait_hint=None,
                              n_resamples=None):
    """MCC 树 → 年 × 方向 的分支-年迁移矩阵。

    口径（对齐 VirPhyKit TempMig / Brynildsrud L4 脚本）：
      * 状态 = 各节点注释的 MAP 状态（max.set / K.set argmax / K 直读）
      * 每条分支（父→子）对 [ceil(T)-ceil(h_child)-floor(len), ceil(T)-ceil(h_child)]
        内的每个日历年 +1（区间端点含）
      * mostcurrentyear = 第二大叶年（复刻其 [-2] 口径，抗单个极端新近样本）
      * rootyear = mostcurrentyear - ceil(T)
    返回 summary dict；out_dir 给定时写出 TransposedMatrix.txt（Year 行 × 方向列，
    与 VirPhyKit TempMig 的 outfile_transposed.txt 同布局）。
    """
    root, _ = load_mcc_tree(tree_path)
    return migration_matrix_from_tree(root, out_dir=out_dir, trait_hint=trait_hint)


def migration_matrix_from_tree(root, out_dir=None, trait_hint=None):
    """同上，但直接吃已解析的树（供 .trees 后验集逐棵复用）。"""
    trait = trait_hint or detect_trait(root)
    if not trait:
        raise ValueError('树上没有找到离散性状注释（需要 BEAST MCC 的 K.set/max.set/'
                         'K+K.prob 注释）；叶名本身不能当状态，拒绝重建迁移矩阵')
    ages, T = node_ages(root)
    states = {id(n): node_state(n, trait) for n, _p, _d in _walk(root)}
    # 叶名兜底在 node_state 内做了；若仍全是 Unknown 说明注释无法解析
    if all(s == 'Unknown' for s in states.values()):
        raise ValueError('树上没有任何节点可提取离散状态，无法构建迁移矩阵')

    ty = _tip_years(root)
    distinct = sorted(set(ty))
    if len(distinct) >= 2:
        mostcurrent = int(distinct[-2])   # 复刻 VirPhyKit 第二大（去重）口径
    elif distinct:
        mostcurrent = int(distinct[-1])
    else:
        mostcurrent = 0
    ceilT = math.ceil(T)
    rootyear = int(mostcurrent - ceilT)

    uniq = sorted(set(states.values()))
    row_keys = [f'{a}_to_{b}' for a in uniq for b in uniq]
    calendar = {r: {} for r in row_keys}

    n_branch = 0
    for n, p, _d in _walk(root):
        if p is None:
            continue
        n_branch += 1
        h = ages[id(n)]
        end = int(math.ceil(T) - math.ceil(h)) + rootyear
        start = end - int(math.floor(n['length'] or 0.0))
        row = f'{states[id(p)]}_to_{states[id(n)]}'
        bucket = calendar[row]
        for z in range(start, end + 1):
            bucket[z] = bucket.get(z, 0) + 1

    years = list(range(rootyear, mostcurrent + 1))
    matrix = {r: [calendar[r].get(y, 0) for y in years] for r in row_keys}
    totals = {r: sum(v) for r, v in matrix.items()}

    out_file = None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, 'TransposedMatrix.txt')
        with open(out_file, 'w', encoding='utf-8', newline='\n') as f:
            f.write('Year\t' + '\t'.join(row_keys) + '\n')
            for i, y in enumerate(years):
                f.write(str(y) + '\t' + '\t'.join(str(matrix[r][i]) for r in row_keys) + '\n')

    diag = {r: totals[r] for r in row_keys if r.split('_to_')[0] == r.split('_to_')[1]}
    offd = {r: totals[r] for r in row_keys if r.split('_to_')[0] != r.split('_to_')[1]}
    return {
        'trait': trait or '(auto)',
        'n_branches': n_branch,
        'n_nodes': len(states),
        'n_tips': sum(1 for n, _p, _d in _walk(root) if not n['children']),
        'tree_height': T,
        'mostcurrent_year': mostcurrent,
        'root_year': rootyear,
        'years': years,
        'directions': row_keys,
        'matrix': matrix,
        'totals': totals,
        'top_migrations': sorted(offd.items(), key=lambda kv: -kv[1])[:12],
        'diag_lineage_years': diag,
        'out_file': out_file,
        'calendar_note': ('分支-年计数：每条分支（父状态→子状态）对其跨越的每个日历年 +1；'
                          '对角线 = 各状态内的谱系-年，非对角线 = 发生状态改变的分支所跨年份'),
    }


# ================================================================
# RRT（MCC 版）+ RSPP-Viz —— 根状态后验
# ================================================================

def _post_dict(pairs):
    d = {}
    for s, p in pairs:
        d[s] = d.get(s, 0.0) + p
    return d


def rspp_from_mcc(tree_path, trait_hint=None):
    """单棵 MCC 树的根状态后验（RSPP-Viz 同款输入）。"""
    root, _ = load_mcc_tree(tree_path)
    pairs = root_state_posterior(root, trait_hint)
    trait = trait_hint or detect_trait(root)
    d = _post_dict(pairs)
    total = sum(d.values())
    return {
        'trait': trait or '(auto)',
        'states': list(d.keys()),
        'probs': [d[s] for s in d],
        'prob_sum': round(total, 6),
        'map_state': max(d, key=d.get),
        'n_tip_annotated': None,
        'file': os.path.basename(tree_path),
    }


def rrt_from_mcc(original_path, randomized_paths, trait_hint=None):
    """区域随机化检验：真树 vs 随机副本的根状态后验表 + 判据。

    判据（对齐 VirPhyKit RRT）：真数据最大后验状态的根后验，大于所有随机副本
    中该状态的最大根后验 → PASS。
    """
    orig, _ = load_mcc_tree(original_path)
    op = _post_dict(root_state_posterior(orig, trait_hint))
    trait = trait_hint or detect_trait(orig)
    target = max(op, key=op.get)
    rows, rand_target = [], []
    for i, rp in enumerate(randomized_paths, 1):
        r, _ = load_mcc_tree(rp)
        d = _post_dict(root_state_posterior(r, trait_hint))
        rows.append({'replicate': f'Random{i}', 'file': os.path.basename(rp),
                     'probs': d})
        rand_target.append(d.get(target, 0.0))
    max_rand = max(rand_target) if rand_target else 0.0
    real_p = op.get(target, 0.0)
    return {
        'trait': trait or '(auto)',
        'target_state': target,
        'states': sorted({k for d in [op] + [r['probs'] for r in rows] for k in d}),
        'real': op,
        'randomized': rows,
        'min_rand_target': min(rand_target) if rand_target else None,
        'max_rand_target': max_rand,
        'real_prob': real_p,
        'passed': bool(rand_target) and real_p > max_rand,
        'n_random': len(rows),
        'verdict': ('RRT 通过：真数据根状态 %s 的后验 %.4f > 随机副本最大值 %.4f'
                    % (target, real_p, max_rand)) if (rand_target and real_p > max_rand)
                   else ('RRT 未通过：真数据根状态 %s 的后验 %.4f 未超过随机副本最大值 %.4f'
                         % (target, real_p, max_rand)),
    }


# ================================================================
# BSP-Viz —— 贝叶斯天际线
# ================================================================

_BSP_COL_ALIASES = {
    'time': ('time', 'year', 'x'),
    'median': ('median', 'median_pop', 'pop_median'),
    'lower': ('lower', 'lower_hpd', 'hpd_lower', 'lower_95_hpd', 'q025', '2.5%'),
    'upper': ('upper', 'upper_hpd', 'hpd_upper', 'upper_95_hpd', 'q975', '97.5%'),
}


def _bsp_col(header, kind):
    al = _BSP_COL_ALIASES[kind]
    for h in header:
        if h.strip().lower() in al:
            return h
    return None


def bsp_from_tsv(path):
    """Tracer 式天际线表（Time/Median/Lower/Upper）→ 数据（精确口径）。"""
    from Virus_Platform_Core.phylodyn_kit import _read_table
    rows = _read_table(path)
    hdr = {k: _bsp_col(rows[0], k) for k in _BSP_COL_ALIASES}
    missing = [k for k, v in hdr.items() if v is None]
    if missing:
        raise ValueError(f'天际线表缺列 {missing}；表头需含 Time/Median/Lower/Upper'
                         f'（实际: {list(rows[0])})')
    data = {'time': [], 'median': [], 'lower': [], 'upper': []}
    for r in rows:
        try:
            t = float(r[hdr['time']])
            me = float(r[hdr['median']])
            lo = float(r[hdr['lower']])
            up = float(r[hdr['upper']])
        except (TypeError, ValueError):
            continue
        data['time'].append(t)
        data['median'].append(me)
        data['lower'].append(lo)
        data['upper'].append(up)
    if len(data['time']) < 2:
        raise ValueError('天际线表可用数据行不足 2 行')
    return {'source': 'tsv', 'n': len(data['time']), 'data': data,
            'approx': False,
            'note': '输入为现成天际线表（Time/Median/Lower/Upper），精确口径'}


def bsp_from_beast_log(path, present_year=None, grid=120):
    """BEAST1 .log（bayesianSkyline.popSizesN/groupSizesN + treeModel.rootHeight）
    → 天际线数据。**近似口径**：组界按该样本根高的等区间数比例近似（逐样本
    step 曲线后在公共网格上取中位数与 95% 分位）；精确组界需要 .trees 文件，
    建议用 Tracer 导出的 TSV 走 bsp_from_tsv。
    """
    text = read_text(path)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    header = None
    samples = []
    for ln in lines:
        if ln.lstrip().startswith('#'):
            continue
        cols = ln.split('\t')
        if len(cols) == 1:
            cols = ln.split()
        if header is None:
            header = [c.strip() for c in cols]
            continue
        if len(cols) != len(header):
            continue
        samples.append(dict(zip(header, cols)))
    if not samples:
        raise ValueError('.log 里没有数据行')
    pop_cols = sorted([c for c in header if re.fullmatch(r'bayesianSkyline\.popSizes\d+', c)],
                      key=lambda c: int(re.search(r'(\d+)$', c).group(1)))
    gs_cols = sorted([c for c in header if re.fullmatch(r'bayesianSkyline\.groupSizes\d+', c)],
                     key=lambda c: int(re.search(r'(\d+)$', c).group(1)))
    rh_col = next((c for c in header if c == 'treeModel.rootHeight'), None)
    if not (pop_cols and gs_cols and rh_col):
        # 退一步：单列空格分隔向量
        pop_vec = next((c for c in header if c.endswith('.popSizes')), None)
        gs_vec = next((c for c in header if c.endswith('.groupSizes')), None)
        if pop_vec and gs_vec and rh_col:
            pop_cols, gs_cols = [pop_vec], [gs_vec]
        else:
            raise ValueError('.log 里找不到 bayesianSkyline.popSizes/groupSizes/'
                             'treeModel.rootHeight 列（BEAST2 的 log 口径不同，'
                             '请用 Tracer 导出 TSV）')

    def _vec(srow, cols):
        if len(cols) == 1:
            return [float(x) for x in srow[cols[0]].replace(',', ' ').split()]
        return [float(srow[c]) for c in cols]

    per_grid = []
    t_max = 0.0
    for srow in samples:
        try:
            ps = _vec(srow, pop_cols)
            gs = _vec(srow, gs_cols)
            T = float(srow[rh_col])
        except (TypeError, ValueError):
            continue
        if not ps or len(ps) != len(gs) or T <= 0:
            continue
        total = sum(gs)
        if total <= 0:
            continue
        bounds, cum = [0.0], 0
        for g in gs:
            cum += g
            bounds.append(T * cum / total)
        per_grid.append((bounds, ps))
        t_max = max(t_max, T)
    if len(per_grid) < 2:
        raise ValueError('.log 里可解析的天际线样本不足 2 个')
    # ⚠️ 2026-10-03 修：present 不给**不再编造 2026.0** —— 硬编码会把"相对时间"
    #    冒充成日历年。缺参时横轴改为相对时间（根 = 0，向现代递增），与
    #    `ltt_from_tree`「无定年 → 横轴 = 根起算距离」同一口径，并在 note 里说明。
    if present_year:
        present = float(present_year)
    else:
        present = None
    ages = [t_max * (i + 1) / grid for i in range(grid)]
    med, lo, up, times = [], [], [], []
    for a in ages:
        vals = []
        for bounds, ps in per_grid:
            gi = bisect.bisect_right(bounds, a) - 1
            gi = min(max(gi, 0), len(ps) - 1)
            vals.append(ps[gi])
        vals.sort()
        times.append(round((present - a) if present is not None else a, 4))
        med.append(_quant(vals, 0.5))
        lo.append(_quant(vals, 0.025))
        up.append(_quant(vals, 0.975))
    note = ('BEAST log 近似口径：组界按「等区间数 × 逐样本根高」近似；'
            '精确组界需要 .trees，建议改用 Tracer 导出的 TSV')
    if present is None:
        note += ('；**调用方未提供 present 年** → 横轴为相对时间'
                 '（根 = 0，单位年，向现代递增），不是日历年')
    return {'source': 'beast_log', 'n_samples': len(per_grid), 'n': len(times),
            'data': {'time': times, 'median': med, 'lower': lo, 'upper': up},
            'approx': True,
            'axis': 'calendar_year' if present is not None else 'relative',
            'present': present,
            'note': note}


def _quant(sorted_vals, q):
    if not sorted_vals:
        return 0.0
    i = q * (len(sorted_vals) - 1)
    lo = int(math.floor(i))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = i - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


# ================================================================
# TreeTime-RTT —— 真引擎 + date-permutation；不可用时自研回归兜底

# ================================================================
# TreeTime-RTT —— 真引擎 + date-permutation；不可用时自研回归兜底
# ================================================================

def treetime_available():
    try:
        import treetime  # noqa: F401
        return True, 'treetime 可用'
    except Exception as e:  # ImportError 及其它
        return False, f'本机未安装 treetime（pip install treetime）：{e}'


def _read_dates_meta(meta_path, name_col='name', date_col='date'):
    """元数据 CSV/TSV → ({name: 小数年}, {name: location}, [(行号, name, 原因)])。"""
    from Virus_Platform_Core.phylodyn_kit import _read_table
    rows = _read_table(meta_path)
    if name_col not in rows[0] or date_col not in rows[0]:
        raise ValueError(f'元数据表需要 {name_col},{date_col} 列；实际表头: {list(rows[0])}')
    dates, locs, bad = {}, {}, []
    for i, r in enumerate(rows, 2):
        nm = (r.get(name_col) or '').strip()
        raw = (r.get(date_col) or '').strip()
        if not nm:
            continue
        try:
            dates[nm] = float(raw)
        except ValueError:
            ok, v, err = validate_date_str(raw)
            if ok:
                dates[nm] = date_to_decimal(v)
            else:
                bad.append((i, nm, err))
                continue
        loc = (r.get('location') or '').strip()
        if loc:
            locs[nm] = loc
    if len(dates) < 3:
        raise ValueError('有效日期行不足 3 条' + (f'；首条错误: {bad[0]}' if bad else ''))
    return dates, locs, bad


def _region_lookup(mapping_path):
    """Mapping.txt：组<TAB>地点 → {地点: 组}。"""
    if not mapping_path or not os.path.exists(mapping_path):
        return {}
    out = {}
    for line in read_text(mapping_path).splitlines():
        parts = line.rstrip('\r\n').split('\t')
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            out[parts[1].strip()] = parts[0].strip()
    return out


def _rtt_points_and_fit(dates, dists, regions=None):
    """根到尾散点 + 最小二乘拟合。

    返回 `(points, reg, r2, n_meta_not_on_tree)`。

    ⚠️ 2026-10-03 修：只回归**两边都有**的叶（`set(dates) & set(dists)`）——
    原实现 `names = sorted(dates)` 后直接 `dists[n]`，元数据里存在、但树上没有
    的叶（常见：元数据整表导入、树只含子集）会 **KeyError 整条链崩掉**。
    落在树外的元数据条数由第 4 个返回值如实上报，调用方写进结果（不静默）。
    """
    from scipy.stats import linregress
    names = sorted(set(dates) & set(dists))
    if len(names) < 2:
        raise ValueError(
            '可比对的（树 × 元数据）叶不足 2 个：%d 个元数据名字里只有 %d 个能在'
            '树上找到 —— 请核对叶名与元数据 name 列是否逐字一致'
            % (len(dates), len(names)))
    xs = [dates[n] for n in names]
    ys = [dists[n] for n in names]
    reg = linregress(xs, ys)
    r2 = reg.rvalue ** 2
    points = [{'name': n, 'date': dates[n], 'dist': dists[n],
               'region': (regions or {}).get(n)} for n in names]
    return points, reg, r2, len(set(dates) - set(dists))


def _perm_test_r2(xs, ys, n_perm=100, seed=0):
    """打乱日期重算 R² 的零分布 → 经验 p（≥ 观测 R² 的比例）。"""
    import random
    rng = random.Random(seed)
    from scipy.stats import linregress
    obs = linregress(xs, ys).rvalue ** 2
    ge = 0
    for _ in range(int(n_perm)):
        perm = xs[:]
        rng.shuffle(perm)
        if linregress(perm, ys).rvalue ** 2 >= obs:
            ge += 1
    return obs, (ge + 1) / (n_perm + 1)


def treetime_rtt(aln_path, tree_path, meta_path, out_dir,
                 mapping_path=None, n_perm=100, time_marginal=False,
                 prog=None):
    """TreeTime 根到尾回归 + 时间树。

    三道前置精确校验（缺日期 / 缺序列 直接中文报错，不进引擎）；
    距离用 mutation_length 累加的遗传距离（时间枝长会得到 R²≡1 的假象）。
    """
    ok, why = treetime_available()
    if not ok:
        raise RuntimeError(why)
    from treetime import TreeTime
    from treetime.utils import parse_dates
    from Bio import AlignIO, Phylo
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO as _SeqIO
    import tempfile

    def _log(m):
        if prog:
            prog(m)

    dates, locs, bad = _read_dates_meta(meta_path)
    if bad:
        _log(f'元数据 {len(bad)} 行日期无法解析，已剔除（首条: 行{bad[0][0]} {bad[0][1]} {bad[0][2]}）')

    aln = AlignIO.read(aln_path, 'fasta')
    recs = []
    for rec in aln:
        recs.append(SeqRecord(Seq(str(rec.seq).replace('U', 'T')), id=rec.id, description=''))
    seq_len = len(aln[0])
    tree_tips = {c.name for c in Phylo.read(tree_path, 'newick').get_terminals()}
    missing = tree_tips - set(dates)
    if missing:
        raise ValueError(f'树上 {len(missing)} 个叶缺日期（例如 {sorted(missing)[:3]}）；'
                         'TreeTime 遇到缺日期会中止，请先补齐元数据或剔除这些叶')
    no_seq = tree_tips - {r.id for r in recs}
    if no_seq:
        raise ValueError(f'比对缺少 {len(no_seq)} 个树叶的序列（例如 {sorted(no_seq)[:3]}）；'
                         '比对与树必须是同一批序列（名字逐字一致）')
    with tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False,
                                     encoding='utf-8', newline='\n') as tf:
        _SeqIO.write(recs, tf, 'fasta')
        tmp_fas = tf.name
    # 只留「树上有日期」的叶，写干净的小数年 CSV（treetime parse_dates 吃这个最稳）
    with tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False,
                                     encoding='utf-8', newline='\n') as td2:
        td2.write('name,date\n')
        for k in sorted(tree_tips & set(dates)):
            td2.write(f'{k},{dates[k]}\n')
        tmp_dates = td2.name
    try:
        _log(f'TreeTime 载入 {len(recs)} 条 × {seq_len} bp ...')
        tt = TreeTime(aln=tmp_fas, tree=tree_path, gtr='Jukes-Cantor',
                      dates=parse_dates(tmp_dates),
                      seq_len=seq_len, verbose=1)
        tt.reroot(root='least-squares')
        _log('最小二乘重根完成，推断分子钟 ...')
        tt.run(root='best', branch_length_mode='joint', infer_gtr=False,
               infer_clock=True, resolve_polytomies=True,
               time_marginal=bool(time_marginal), max_iter=3)
    finally:
        for tmp in (tmp_fas, tmp_dates):
            try:
                os.unlink(tmp)
            except OSError:
                pass

    rate = float(tt.date2dist.clock_rate) if hasattr(tt, 'date2dist') and tt.date2dist else None
    engine_r2 = None
    try:
        engine_r2 = float(tt.clock_model.get('r_val'))
    except (AttributeError, KeyError, TypeError):
        pass
    # 根到尾遗传距离：沿 mutation_length 累加（branch_length 在定年后是时间单位，
    # 用它回归会得到恒等 R²=1 的假象）
    def _genetic_dist(tip):
        d = 0.0
        for cl in tt.tree.get_path(tip):
            ml = getattr(cl, 'mutation_length', None)
            d += ml if ml is not None else (cl.branch_length or 0.0)
        return d
    dists = {c.name: _genetic_dist(c) for c in tt.tree.get_terminals()
             if c.name and not getattr(c, 'bad_branch', False)}
    regions = {}
    loc2reg = _region_lookup(mapping_path)
    if locs and loc2reg:
        for nm, loc in locs.items():
            regions[nm] = loc2reg.get(loc, 'Unknown')
    elif locs:
        regions = {nm: loc for nm, loc in locs.items()}

    points, reg, r2, n_meta_not_on_tree = _rtt_points_and_fit(dates, dists, regions)
    if n_meta_not_on_tree:
        _log('元数据有 %d 条名字不在树上 → 本次回归**已忽略**（只回归树 × 元数据'
             '的交集；原实现会直接 KeyError）' % n_meta_not_on_tree)
    xs = [p['date'] for p in points]
    ys = [p['dist'] for p in points]
    obs_r2, perm_p = _perm_test_r2(xs, ys, n_perm=n_perm)

    os.makedirs(out_dir, exist_ok=True)
    nwk = os.path.join(out_dir, 'timetree_inferred.nwk')
    Phylo.write(tt.tree, nwk, 'newick')
    tsv = os.path.join(out_dir, 'rtt_points.tsv')
    with open(tsv, 'w', encoding='utf-8', newline='\n') as f:
        f.write('name\tdate\tdist\tregion\n')
        for p in points:
            f.write(f"{p['name']}\t{p['date']}\t{p['dist']}\t{p['region'] or ''}\n")
    pdf = None
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5.5))
        regs = sorted({p['region'] for p in points if p['region']} or {'(全部)'})
        cmap = plt.get_cmap('tab20', max(len(regs), 1))
        for i, rg in enumerate(regs):
            seg = [p for p in points if (p['region'] or '(全部)') == rg]
            ax.scatter([p['date'] for p in seg], [p['dist'] for p in seg],
                       s=26, color=cmap(i), label=rg, alpha=0.85)
        slope = reg.slope
        inter = reg.intercept
        xs2 = [min(xs), max(xs)]
        ax.plot(xs2, [slope * x + inter for x in xs2], color='#00A0E9', lw=2)
        ax.set_xlabel('Sampling date (decimal year)')
        ax.set_ylabel('Distance from root (subs/site)')
        ax.set_title(f'Root-to-tip (TreeTime)  rate={rate:.3e}  R²={r2:.3f}'
                     if rate else f'Root-to-tip  R²={r2:.3f}')
        if 1 <= len(regs) <= 12:
            ax.legend(fontsize=7, loc='best')
        fig.tight_layout()
        pdf = os.path.join(out_dir, 'RootToTip.pdf')
        fig.savefig(pdf, format='pdf')
        plt.close(fig)
    except Exception:
        pdf = None

    return {
        'engine': 'treetime',
        'n_tips': len(points),
        'n_meta_not_on_tree': n_meta_not_on_tree,
        'clock_rate': rate,
        'engine_r2': engine_r2,
        'r2': r2,
        'slope': reg.slope,
        'p_value': reg.pvalue,
        'perm_n': n_perm,
        'perm_r2': obs_r2,
        'perm_p': perm_p,
        'root_date': min(dates.values()),
        'timetree': os.path.basename(nwk),
        'points_file': os.path.basename(tsv),
        'rtt_pdf': os.path.basename(pdf) if pdf else None,
        'n_regions': len({p['region'] for p in points if p['region']}),
        'points': points,
    }


def inhouse_rtt(tree_path, meta_path, n_perm=100):
    """treetime 不可用时的自研根到尾回归（免外部引擎；只做回归与置换检验）。"""
    from Virus_Platform_Core import phylogeo as pg
    dates, locs, bad = _read_dates_meta(meta_path)
    root = pg.parse_newick(read_text(tree_path))
    dists = {}
    for n, _p, d in _walk(root):
        if not n['children'] and n['name']:
            nm = n['name'].split('|')[0].strip()
            if nm in dates:
                dists[nm] = d
    if len(dists) < 3:
        raise ValueError('树上可配对日期的叶不足 3 个')
    points, reg, r2, n_meta_not_on_tree = _rtt_points_and_fit(dates, dists)
    obs_r2, perm_p = _perm_test_r2([p['date'] for p in points],
                                   [p['dist'] for p in points], n_perm=n_perm)
    return {'engine': 'inhouse', 'n_tips': len(points), 'clock_rate': reg.slope,
            'n_meta_not_on_tree': n_meta_not_on_tree,
            'r2': r2, 'slope': reg.slope, 'p_value': reg.pvalue,
            'perm_n': n_perm, 'perm_r2': obs_r2, 'perm_p': perm_p,
            'root_date': min(dates.values()), 'timetree': None,
            'points_file': None, 'rtt_pdf': None, 'n_regions': 0, 'points': points}



# ================================================================
# BEAST .trees 后验集解析 —— 修「单棵 MCC 的 MAP 状态」局限
# ================================================================

def iter_beast_trees(path, burnin=0.1, max_trees=None):
    """BEAST .trees 后验树集 → 逐棵 (样本序号, root, T)。

    * burnin：按比例丢弃最前段（默认 10%）
    * max_trees：超过时**等步长抽样**（横跨整段，不只取前段）
    * 支持 Translate 块；每棵树必须是单行（BEAST 标准输出即如此）
    * ⚠️ 2026-10-03 改**流式**：原实现 `read_text` 整文件 + 把所有 tree 行
      缓冲进列表 —— 大 .trees（数千棵 × 长枝长/叶名）光文本就可达 GB 级，
      再加解析出的树对象内存翻倍。现在两遍**逐行**读：
        第一遍只数 `tree ` 行数（顺带在首个 tree 行前抓 Translate 块，它总在
          文件头），第二遍按等步长选出的行号**只解析被选中的行**。
      峰值内存 ≈ 被选中的 `max_trees` 棵（`max_trees` 因此同时封顶解析量与内存），
      且不把整份文本读进内存。行号口径与旧实现完全一致（同一批等步长下标）。
    """
    n_tree, header = 0, []
    with open(path, encoding='utf-8-sig', errors='replace') as f:
        for ln in f:
            if '(' in ln and ln.lstrip().lower().startswith('tree '):
                n_tree += 1
            elif n_tree == 0:
                header.append(ln)     # 首个 tree 行之前的头部（Translate 在这里）
    if not n_tree:
        raise ValueError(f'{path} 里没有 tree 行（不是 BEAST .trees 后验树集？）')
    trans = _parse_translate(''.join(header))
    skip = int(n_tree * float(burnin))
    n_pick = n_tree - skip
    if max_trees and n_pick > max_trees:
        step = n_pick / float(max_trees)
        idx = [skip + int(i * step) for i in range(max_trees)]
    else:
        idx = list(range(skip, n_tree))
    sel = set(idx)
    out = []
    k = -1
    with open(path, encoding='utf-8-sig', errors='replace') as f:
        for ln in f:
            if '(' not in ln or not ln.lstrip().lower().startswith('tree '):
                continue
            k += 1
            if k not in sel:
                continue
            root = parse_beast_newick(ln)
            if trans:
                for n, _p, _d in _walk(root):
                    if not n['children'] and n['name'] in trans:
                        n['name'] = trans[n['name']]
            out.append((k + 1, root, tree_height(root)))
    return out


def rspp_posterior(trees_path, trait_hint=None, burnin=0.1, max_trees=400):
    """根状态后验的**后验集口径**：逐棵树读根状态分布，
    输出均值后验向量 + 每棵树 MAP 状态的频率 + 每状态后验的范围。"""
    trees = iter_beast_trees(trees_path, burnin=burnin, max_trees=max_trees)
    if not trees:
        raise ValueError('后验集里没有可解析的树')
    sums, map_counts, per_tree = {}, {}, {}
    trait = None
    for _i, root, _T in trees:
        trait = trait or detect_trait(root)
        pairs = root_state_posterior(root, trait_hint)
        d = _post_dict(pairs)
        for st, p in d.items():
            sums[st] = sums.get(st, 0.0) + p
            per_tree.setdefault(st, []).append(p)
        mapst = max(d, key=d.get)
        map_counts[mapst] = map_counts.get(mapst, 0) + 1
    n = len(trees)
    states = sorted(sums)
    return {
        'mode': 'posterior',
        'trait': trait_hint or trait or '(auto)',
        'n_trees': n,
        'states': states,
        'mean_probs': {st: sums[st] / n for st in states},
        'map_freq': {st: map_counts.get(st, 0) / n for st in states},
        'map_consensus': max(map_counts, key=map_counts.get) if map_counts else None,
        'range': {st: [min(per_tree[st]), max(per_tree[st])] for st in states},
        'note': '后验集口径：跨后验树取均值，不再是单棵 MCC 的点估计',
    }


def tempmig_posterior(trees_path, trait_hint=None, burnin=0.1, max_trees=100,
                      out_dir=None):
    """迁移矩阵的**后验集口径**：逐棵树跑分支-年矩阵后取均值。

    各树年轴略有差异 → 取交集轴 [max(根年), min(现年)]；
    输出每方向「分支-年总数」的均值与标准差 + 交集轴上的均值矩阵。
    """
    trees = iter_beast_trees(trees_path, burnin=burnin, max_trees=max_trees)
    if not trees:
        raise ValueError('后验集里没有可解析的树')
    trait = None
    per_totals, per_year, rootyears, highs = {}, {}, [], []
    n = 0
    for _i, root, _T in trees:
        trait = trait_hint or detect_trait(root)
        m = migration_matrix_from_tree(root, trait_hint=trait)
        n += 1
        rootyears.append(m['root_year'])
        highs.append(m['mostcurrent_year'])
        for d, tot in m['totals'].items():
            per_totals.setdefault(d, []).append(tot)
        for d, row in m['matrix'].items():
            for y, v in zip(m['years'], row):
                per_year.setdefault((d, y), []).append(v)
    lo, hi = max(rootyears), min(highs)
    dirs = sorted(per_totals)
    mean_year = {d: [] for d in dirs}
    if lo <= hi:
        for d in dirs:
            for y in range(lo, hi + 1):
                vals = per_year.get((d, y))
                mean_year[d].append(round(sum(vals) / len(vals), 4) if vals else 0.0)
    years = list(range(lo, hi + 1)) if lo <= hi else []
    import math as _math
    mean_tot = {d: sum(per_totals[d]) / n for d in dirs}
    sd_tot = {d: _math.sqrt(sum((x - mean_tot[d]) ** 2 for x in per_totals[d]) / n)
              for d in dirs}
    res = {
        'mode': 'posterior',
        'trait': trait_hint or trait or '(auto)',
        'n_trees': n,
        'year_axis_note': ('交集年轴 [%d, %d]（各后验树根年/现年略有差异）' % (lo, hi)
                           if lo <= hi else '各树年轴无交集'),
        'years': years,
        'directions': dirs,
        'matrix_mean': mean_year,
        'totals_mean': {d: round(mean_tot[d], 3) for d in dirs},
        'totals_sd': {d: round(sd_tot[d], 3) for d in dirs},
        'top_migrations': sorted(mean_tot.items(), key=lambda kv: -kv[1])[:12],
        'out_file': None,
        'note': '后验集口径：跨后验树取均值 ± 标准差，不再是单棵 MCC 的点估计',
    }
    if out_dir and years:
        os.makedirs(out_dir, exist_ok=True)
        fp = os.path.join(out_dir, 'TransposedMatrix_posterior.txt')
        with io.open(fp, 'w', encoding='utf-8', newline='\n') as f:
            f.write('Year\t' + '\t'.join(dirs) + '\n')
            for i, y in enumerate(years):
                f.write(str(y) + '\t' + '\t'.join(
                    ('%.3f' % mean_year[d][i]) for d in dirs) + '\n')
        res['out_file'] = 'pdtempmig/TransposedMatrix_posterior.txt'
    return res


# ================================================================
# MCC 区域随机化（Slatkin–Maddison/BaTS 式 Fitch 置换，接通 t-pdrrt 断点）
# ================================================================

_TRAIT_FAMILY = ('', '.prob', '.set', '.set.prob')


def _ann_str(v):
    """注解值 → newick 注解文本（字符串加引号，列表加花括号）。"""
    if isinstance(v, list):
        return '{' + ','.join(('"%s"' % x) if re.search(r'[^A-Za-z0-9_.]', str(x)) else str(x)
                              for x in v) + '}'
    s = str(v)
    if s == '' or re.search(r'[^A-Za-z0-9_.+-]', s):
        return '"%s"' % s
    return s


def _name_str(name):
    """Newick 名字 → 可安全写回的字符串（2026-10-03）。

    含 `(),:;[]`、空白或单引号的名字必须用**单引号**包裹（内部的 `'` 按 Newick
    标准翻倍转义），否则 `A B` / `x:1` / `(a)` 这类名字回读（或被 BEAST 读）时
    会被解析成别的东西 —— 严重时整棵树解析失败。不含特殊字符的名字原样输出，
    与既有产物保持逐字一致。
    """
    s = '' if name is None else str(name)
    if not s:
        return ''                       # 无名字：保持"什么都不写"（不写空引号）
    if any(c in s for c in "(),:;[]' \t\r\n"):
        return "'" + s.replace("'", "''") + "'"
    return s


def write_beast_newick(root, trait=None, states=None, root_set=None):
    """带注解 newick 序列化。

    trait 给定时：抹掉该性状族（trait/.prob/.set/.set.prob）原注解，按 states
    （{id(node): 状态}）重写（叶补 trait.prob=1.0）；root_set 给定时在根上写
    trait.set/trait.set.prob（等概率）。其余注解（height/posterior/…）原样保留。
    """
    strip = {trait + suf for suf in _TRAIT_FAMILY} if trait else set()

    def _anns(n):
        anns = {k: v for k, v in n['anns'].items() if k not in strip}
        if trait and states:
            st = states.get(id(n))
            if st:
                anns[trait] = st
                if not n['children']:
                    anns[trait + '.prob'] = 1.0
        if trait and root_set is not None and n is root and root_set:
            anns[trait + '.set'] = sorted(root_set)
            anns[trait + '.set.prob'] = [round(1.0 / len(root_set), 6)] * len(root_set)
        if not anns:
            return ''
        return '[&' + ','.join('%s=%s' % (k, _ann_str(v)) for k, v in anns.items()) + ']'

    def _rec(n):
        # 叶名/内部名经 `_name_str` 引号化：含 `(),:;[]`/空白的名字不引号
        # 会写出**非法 newick**（回读或喂 BEAST 时被解析成别的东西）
        head = _name_str(n['name']) + _anns(n)
        if n['children']:
            return '(' + ','.join(_rec(c) for c in n['children']) + ')' + head + ':' + repr(float(n['length'] or 0.0))
        return head + ':' + repr(float(n['length'] or 0.0))

    return _rec(root) + ';'


def _root_candidate_set(root, states):
    """Fitch up-pass 只算根候选集（states 按**叶名**键）。"""
    from Virus_Platform_Core import phylogeo as pg
    sets = {}
    for n in pg._iter_post(root):
        if not n['children']:
            st = states.get(n['name'], '')
            sets[id(n)] = {st} if st else {'Unknown'}
            continue
        inter, union = None, set()
        for c in n['children']:
            cs = sets[id(c)]
            union |= cs
            inter = cs if inter is None else (inter & cs)
        sets[id(n)] = inter if inter else union
    return sets[id(root)]


def mcc_randomize_test(tree_path, n=20, out_dir=None, seed=None, trait_hint=None):
    """MCC 树的区域随机化检验（Fitch/简约法口径，Slatkin–Maddison 式）。

    统计量 = Fitch 最少迁移步数；零分布 = 打乱叶区划标注后重跑 Fitch 的步数；
    p = (1 + #{随机步数 ≤ 真实步数}) / (N + 1)。随机树同时落盘（其余注解原样
    保留、性状族按重算状态重写、根写候选集等概率）——可直接喂
    「区域随机化检验（MCC）」卡或 VirPhyKit 的 RRT。

    与 t-pdrrt 的分工：t-pdrrt = BEAST 后验口径（用户已有随机化 MCC 树）；
    本工具 = 简约法口径（只有一棵 MCC 树就能跑），回答同一个科学问题。
    """
    import random as _random
    from Virus_Platform_Core import phylogeo as pg
    root, _ = load_mcc_tree(tree_path)
    trait = trait_hint or detect_trait(root)
    if not trait:
        raise ValueError('树上没有离散性状注释（需要 K.set/max.set/K+K.prob），无法随机化')
    tip_names, tip_states = [], []
    for nd, _p, _d in _walk(root):
        if not nd['children'] and nd['name']:
            tip_names.append(nd['name'])
            tip_states.append(node_state(nd, trait))
    if len(tip_names) < 3:
        raise ValueError('可用叶不足 3 条')
    if len(set(tip_states)) < 2:
        raise ValueError('全部叶状态相同，随机化检验无意义')
    states_by_name = dict(zip(tip_names, tip_states))
    _, trans_real, _ = pg._fitch(root, states_by_name)
    steps_real = len(trans_real)

    rng = _random.Random(seed)
    reps, trees_written = [], []
    census = {}
    for st in tip_states:
        census[st] = census.get(st, 0) + 1
    for i in range(1, int(n) + 1):
        perm = tip_states[:]
        rng.shuffle(perm)
        perm_states = dict(zip(tip_names, perm))
        by_id, trans_i, _ = pg._fitch(root, perm_states)
        reps.append(len(trans_i))
        if out_dir:
            states_out = {}
            for nd, _p, _d in _walk(root):
                states_out[id(nd)] = by_id.get(id(nd)) or perm_states.get(nd['name'])
            txt = write_beast_newick(root, trait=trait, states=states_out,
                                     root_set=_root_candidate_set(root, perm_states))
            fp = os.path.join(out_dir, 'random%d.tre' % i)
            with io.open(fp, 'w', encoding='utf-8', newline='\n') as f:
                f.write(txt + '\n')
            trees_written.append(fp)
    p_val = (1 + sum(1 for x in reps if x <= steps_real)) / (len(reps) + 1)
    sig = p_val <= 0.05
    return {
        'trait': trait or '(auto)',
        'n_tips': len(tip_names),
        'n_random': len(reps),
        'seed': seed,
        'steps_real': steps_real,
        'steps_random': reps,
        'steps_random_min': min(reps) if reps else None,
        'steps_random_max': max(reps) if reps else None,
        'steps_random_mean': (sum(reps) / len(reps)) if reps else None,
        'p_value': p_val,
        'significant': sig,
        'verdict': ('显著（p=%.4f ≤ 0.05）：区划与树的结构协同强于随机' % p_val
                    if sig else
                    '不显著（p=%.4f > 0.05）：随机化下迁移步数与真实相当' % p_val),
        'census': census,
        'random_trees': [os.path.basename(x) for x in trees_written],
        'note': ('简约法（Fitch）口径：统计量 = 最少迁移步数；随机树已按重算状态'
                 '重写性状注解，可直接喂「区域随机化检验（MCC）」卡'),
    }


# ================================================================
# TreeDater 真引擎（便携 R + treedater 包，随 3rd/ 分发、目标机免装免联网）
# ================================================================

def find_rscript():
    """定位 Rscript：platform.json tools.rscript 优先，其次 3rd/R/*/bin（打包即用）。"""
    from Virus_Platform_Core.config import get_config, PLATFORM_ROOT
    try:
        p = get_config().tool('rscript')
        if p and os.path.isfile(p):
            return p
    except Exception:
        pass
    import glob as _g
    for pat in ('3rd/R/*/bin/Rscript.exe', '3rd/R/*/bin/Rscript'):
        hits = sorted(_g.glob(os.path.join(PLATFORM_ROOT, *pat.split('/'))))
        if hits:
            return hits[-1]
    return None


def _r_env(exe=None):
    """R 调用的干净环境。

    * 剥掉 LC_*/LANG：Git Bash/某些终端的 LC_CTYPE=C.UTF-8 会让 R 的 Windows
      运行时把基础包 namespace 搞挂（实测）。
    * 剥掉 R_LIBS*/R_ENVIRON/R_PROFILE/R_HOME：免疫宿主机残留配置。
    * R_HOME 按 exe 位置显式回填（bin/..），目标机无注册表也能跑（便携分发）。
    """
    env = dict(os.environ)
    for k in list(env):
        if k.startswith('LC_') or k in ('LANG', 'LANGUAGE', 'R_LIBS', 'R_LIBS_USER',
                                        'R_LIBS_SITE', 'R_ENVIRON', 'R_PROFILE',
                                        'R_HOME'):
            env.pop(k)
    if exe:
        env['R_HOME'] = os.path.dirname(os.path.dirname(exe))
    return env


def treedater_available():
    """(可用?, Rscript 路径或人话原因)——要求 treedater 包可加载。"""
    import subprocess
    exe = find_rscript()
    if not exe:
        return False, '未找到 Rscript（3rd/R/*/bin/Rscript.exe 或 platform.json tools.rscript）'
    try:
        # encoding/errors 显式给：Windows 上 text=True 默认用 ANSI 代码页（GBK），
        # R 的英文输出里混进非 GBK 字节会抛 UnicodeDecodeError（2026-10-03 修）
        r = subprocess.run(
            [exe, '-e', "cat(requireNamespace('treedater', quietly=TRUE))"],
            capture_output=True, text=True, encoding='utf-8', errors='replace',
            timeout=180, env=_r_env(exe))
        ok = r.stdout.strip().endswith('TRUE')
        msg = ('R 在但缺 treedater 包：Rscript -e '
               + chr(34) + 'install.packages(' + chr(34) + 'treedater' + chr(34) + ')'
               + chr(34))
        return (ok, exe if ok else msg)
    except (OSError, subprocess.SubprocessError, UnicodeDecodeError) as e:
        return False, f'Rscript 调用失败: {e}'


def treedater_ltt(tree_path, meta_path, out_dir, seq_len=None, ncpu=1, prog=None):
    """真 TreeDater：R treedater::dater() 定年 + parboot LTT 置信带 PDF。

    元数据 name,date（ISO 或小数年，内部统一小数年）。速率在 R 内用
    「根到端遗传距离 ~ 采样时间」回归直接算出并写入 TD_summary.txt
    （不依赖 stdout 解析），同时落 `rtt_r2`（root-to-tip 回归 R²，常用时间
    信号诊断）。返回 summary；曲线由自研 LTT 在定年树上重算（横轴=日历年）。
    .engine='treedater'。

    ⚠️ 未定根的树（ML 树常态：根三分叉）其「根到端距离」**没有意义**，直接
    回归会得到负斜率。所以 R 侧先判断 `is.rooted()`，未定根时以**最老样本**
    为外类群定根再回归——实测 Example 的 h3n2_na_500：未定根 → −3.49e−3
    (R²=0.85)、定根后 → +3.82e−3 (R²=0.99)，与同树 TreeTime 的 +3.11e−3
    (R²=0.986) 同量级。**已定根的输入树不做任何改动**（行为不变）。
    """
    import subprocess
    exe = find_rscript()
    if not exe:
        raise RuntimeError('未找到 Rscript（3rd/R/*/bin/Rscript.exe 或 platform.json '
                           'tools.rscript）；可改用自研 LTT 口径')
    dates = _read_dates_meta(meta_path)[0]
    seq_len = int(seq_len) if seq_len else 10000
    os.makedirs(out_dir, exist_ok=True)
    meta_csv = os.path.join(out_dir, 'treedater_meta.csv')
    with io.open(meta_csv, 'w', encoding='utf-8', newline='\n') as f:
        f.write('name,date\n')
        for k in sorted(dates):
            f.write('%s,%s\n' % (k, dates[k]))
    BS = chr(92)
    TAB = BS + 't'          # NL 原与 TAB 一起赋值但从未使用（F841），2026-10-03 删
    r_lines = [
        'require(treedater); require(ape)',
        'args <- commandArgs(trailingOnly = TRUE)',
        'tree_file <- args[1]; meta_file <- args[2]',
        'seqlen <- as.numeric(args[3]); out_dir <- args[4]; ncpu <- as.numeric(args[5])',
        'tre <- read.tree(tree_file)',
        'Times <- read.csv(meta_file)',
        'sts <- setNames(Times[,2], Times[,1])',
        # 未定根的 ML 树（根三分叉）→ 根到端距离无意义，先以最老样本为外类群定根
        'if (!is.rooted(tre)) tre <- root(tre, outgroup = names(which.min(sts)),',
        '                                 resolve.root = FALSE)',
        'dtr <- dater(tre, sts, seqlen, clock = "uncorrelated", ncpu = ncpu)',
        'write.tree(dtr, file.path(out_dir, "dated_tree.nwk"))',
        'pdf(file.path(out_dir, "Phylogeny.pdf"), width = 10, height = 8)',
        'plot(dtr, no.mar = TRUE, cex = 0.5); dev.off()',
        'root <- length(tre$tip.label) + 1',
        'tips <- seq_along(tre$tip.label)',
        'gdist <- dist.nodes(tre)[root, tips]',
        'tipdates <- as.numeric(sts[tre$tip.label])',
        'fit <- lm(gdist ~ tipdates)',
        'rate <- as.numeric(fit$coefficients[[2]])',
        'rtt_r2 <- as.numeric(summary(fit)$r.squared)',
        'rootage <- as.numeric(max(node.depth.edgelength(dtr)))',
        'sum_file <- file.path(out_dir, "TD_summary.txt")',
        'writeLines(c(paste0("rate' + TAB + '", ifelse(is.na(rate), "NA",',
        '  format(rate, digits = 6))), paste0("rtt_r2' + TAB + '",',
        '  ifelse(is.na(rtt_r2), "NA", format(rtt_r2, digits = 5))),',
        '  paste0("root_age' + TAB + '",',
        '  format(rootage, digits = 8)), paste0("ok' + TAB + '1")), sum_file)',
        'pb <- tryCatch(parboot(dtr, ncpu = ncpu), error = function(e) NULL)',
        'if (!is.null(pb)) { pdf(file.path(out_dir, "LTT.pdf"), width = 10, height = 8)',
        '  plot(pb); dev.off() }',
    ]
    r_script = os.path.join(out_dir, 'run_treedater.R')
    with io.open(r_script, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(r_lines) + '\n')
    if prog:
        # 实测（476 叶 h3n2、seq_len=1400、ncpu=1）：dater() 约 0.5-1 分钟出定年树与
        # 摘要，但随后的 parboot 置信带还要 ~35 分钟 → 旧文案「1-3 分钟」会让人以为卡死。
        prog('R treedater 定年中（定年约 1 分钟；紧随其后的 parboot 置信带在单核下'
             '可能要 10-40 分钟，请勿中断）...')
    # ⚠️ 不能用 capture_output=True（PIPE 管道）：R 偶发在 parboot 中途退出时，
    # 其 stdout/stderr 管道写句柄可能被孤儿孙进程继承持有，EOF 永不到来 →
    # communicate() 无限阻塞，连 timeout=3600 也救不了（timeout 杀掉直接
    # 子进程后还要 drain 管道，照样卡死；2026-09-18 实测
    # _check_phylodyn_example 卡死 23 分钟即此因）。改为**文件重定向**：
    # 子进程直接写日志文件，run() 只等进程句柄，与管道 EOF 彻底解耦。
    # （shell=False 的列表参数调用，无任何字符串拼壳。）
    _r_out = os.path.join(out_dir, 'treedater.stdout.log')
    _r_err = os.path.join(out_dir, 'treedater.stderr.log')
    with io.open(_r_out, 'wb') as _fo, io.open(_r_err, 'wb') as _fe:
        _r = subprocess.run([exe, r_script, tree_path, meta_csv,
                             str(seq_len), out_dir, str(ncpu)],
                            stdout=_fo, stderr=_fe, shell=False,
                            timeout=3600, env=_r_env(exe))
    with io.open(_r_out, 'r', encoding='utf-8', errors='replace') as f:
        _out_tail = f.read()[-2000:]
    with io.open(_r_err, 'r', encoding='utf-8', errors='replace') as f:
        _err_tail = f.read()[-2000:]
    r = subprocess.CompletedProcess(_r.args, _r.returncode,
                                    stdout=_out_tail, stderr=_err_tail)
    dated = os.path.join(out_dir, 'dated_tree.nwk')
    sum_file = os.path.join(out_dir, 'TD_summary.txt')
    if not os.path.isfile(dated):
        raise RuntimeError('treedater 运行失败：' + (r.stderr or r.stdout or '')[-500:])
    tags = {}
    if os.path.isfile(sum_file):
        for ln in read_text(sum_file).splitlines():
            k, _, v = ln.partition('\t')
            tags[k.strip()] = v.strip()
    if not tags or tags.get('ok') != '1':
        raise RuntimeError('treedater 摘要缺失：' + (r.stderr or r.stdout or '')[-300:])

    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None
    rate = _num(tags.get('rate'))
    root_age_years = _num(tags.get('root_age'))
    rtt_r2 = _num(tags.get('rtt_r2'))
    res = ltt_from_tree(dated, meta_path=meta_path, out_dir=out_dir)
    res.update({
        'engine': 'treedater',
        'rscript': exe,
        'rate': rate,
        # root-to-tip 回归的 R²（时间信号诊断）；未定根树定根后才有意义
        'rtt_r2': rtt_r2,
        'root_age_years': root_age_years,
        'dated_tree': 'pdltt/dated_tree.nwk',
        'phylo_pdf': 'pdltt/Phylogeny.pdf' if os.path.isfile(os.path.join(out_dir, 'Phylogeny.pdf')) else None,
        'ltt_pdf': 'pdltt/LTT.pdf' if os.path.isfile(os.path.join(out_dir, 'LTT.pdf')) else None,
        'r_tail': (r.stderr or r.stdout or '')[-300:],
    })
    if res['axis'] != 'calendar_year':
        res['dated_note'] = (res.get('dated_note') or '') + \
            '（定年树横轴应为日历年；若此提示出现请检查树单位）'
    return res


def ltt_from_tree(tree_path, meta_path=None, use_lsd2=False, seq_len=None,
                  out_dir=None, prog=None):
    """谱系-时间（LTT）曲线。

    meta 给定时：叶年 = 元数据 date（ISO 或小数年），横轴为日历年；
    否则横轴为距今的遗传距离（若树本身已是时间树亦可）。
    use_lsd2=True 且本机有 LSD2 时：先用 LSD2 定年（需要 meta + seq_len）。
    """
    from Virus_Platform_Core import phylogeo as pg

    def _log(m):
        if prog:
            prog(m)

    dated_note = None
    if use_lsd2:
        if not meta_path:
            raise ValueError('LSD2 定年需要日期元数据 CSV（name,date）')
        exe = pg.find_lsd2()
        if not exe:
            raise RuntimeError('未找到 LSD2（platform.json tools.lsd2 或 3rd/tools/lsd2/lsd2.exe）；'
                               '取消勾选 LSD2 也可直接画 LTT（横轴=遗传距离）')
        dates = _read_dates_meta(meta_path)[0]
        seq_len = seq_len or 10000
        os.makedirs(out_dir or '.', exist_ok=True)
        _log('LSD2 定年中（引擎: %s）...' % exe)
        res = pg.lsd2_dating(tree_path, dates, seq_len, out_dir or '.', exe=exe)
        dated_path = res.get('dated_tree') if isinstance(res, dict) else None
        if dated_path and os.path.exists(dated_path):
            tree_path = dated_path
            dated_note = '已先用 LSD2 定年（.result.nwk），LTT 横轴为定年时间'
        else:
            dated_note = 'LSD2 未返回定年树，退回输入树'

    root = pg.parse_newick(read_text(tree_path))
    T = tree_height(root)
    tip_dates = None
    if meta_path:
        dates, _locs, _bad = _read_dates_meta(meta_path)
        m = {}
        for n, _p, _d in _walk(root):
            if not n['children'] and n['name']:
                nm = n['name']
                for cand in (nm, nm.split('|')[0].strip(), nm.split()[0] if nm.split() else nm):
                    if cand in dates:
                        m[nm] = dates[cand]
                        break
        if len(m) >= 3:
            # 时间轴单位检测：时间树的根龄不可能早于最早样本的跨度之下
            # （T >= span 才能把枝长当「年」；遗传距离树的树高 ≪ 年跨度）
            yrs = sorted(m.values())
            span = yrs[-1] - yrs[0]
            if T >= span - 0.01:
                tip_dates = m
            else:
                dated_note = ('输入树的枝长不是「年」（树高 %.4g < 元数据年跨度 %.4g）：'
                              'LTT 横轴为遗传距离；要日历年轴请先 LSD2 定年'
                              % (T, span))
    present = max(tip_dates.values()) if tip_dates else 0.0

    # LTT 曲线（事件法，支持多分叉）：起点 = 根的子数；每经过一个非根内节点
    # 增加「子数 − 1」条谱系（二叉 +1，三叉 +2 …）；终点 = 叶数
    ages, _ = node_ages(root)
    root_age = ages[id(root)]
    events = {}
    for n, p, _d in _walk(root):
        if p is None or not n['children']:
            continue
        a = ages[id(n)]
        if a >= root_age - 1e-12:
            continue
        events[a] = events.get(a, 0) + (len(n['children']) - 1)
    n_lin = max(len(root['children']), 1)
    curve = [(root_age, n_lin)]
    for a in sorted(events, reverse=True):
        n_lin += events[a]
        curve.append((a, n_lin))
    curve.append((0.0, n_lin))
    n_tips = sum(1 for n, _p, _d in _walk(root) if not n['children'])
    T = root_age if root_age > 0 else T

    def _x(a):
        return (present - a) if tip_dates else a

    pts = [{'x': round(_x(a), 6), 'n': n_lin} for a, n_lin in curve]
    #（原 `xs = [p['x'] …]` 从未被使用 → F841，2026-10-03 删；下面用的是 xs2）
    ys = [math.log(p['n']) for p in pts if p['n'] > 0]
    xs2 = [_x(a) for a, n in curve if n > 0]
    slope = None
    if len(xs2) >= 3:
        from scipy.stats import linregress
        slope = linregress(xs2, ys).slope
    csv_out = None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        csv_out = os.path.join(out_dir, 'ltt_curve.tsv')
        with open(csv_out, 'w', encoding='utf-8', newline='\n') as f:
            f.write('x\tn_lineages\n')
            for p in pts:
                f.write(f"{p['x']}\t{p['n']}\n")
    return {
        'n_tips': n_tips,
        'root_age': T,
        'axis': 'calendar_year' if tip_dates else 'subst_distance',
        'present': present if tip_dates else None,
        'n_dated_tips': len(tip_dates) if tip_dates else 0,
        'log_linear_slope': slope,
        'slope_note': ('log N 对时间的回归斜率 ≈ 净多样化率（birth−death 的近似；'
                       '横轴为日历年时单位 = 每年）' if slope is not None else None),
        'dated_note': dated_note,
        'curve': pts,
        'curve_file': os.path.basename(csv_out) if csv_out else None,
    }
