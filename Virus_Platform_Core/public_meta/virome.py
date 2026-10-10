# -*- coding: utf-8 -*-
"""病毒谱：把预筛产物从「样本视角」翻成「病毒视角」。

预筛的 `stat_run_summary.tsv` 回答的是「哪些样本含病毒」（一行一个样本）。
本模块回答另一半问题：**这批样本里到底检出了哪些病毒、各自属于什么分类、
占比多少、在多少个样本里出现过**。

输入（预筛产物，纯读，不联网）
------------------------------
* `<项目>/prescreen/stat_viral_taxa.tsv` —— 每样本的病毒 taxon 明细
* `<项目>/prescreen/stat_run_summary.tsv` —— 样本清单（取分母）

输出
----
见 `build_virome()` / `build_virome_all()` 的返回结构，前端拿它直接渲染表格与图。

口径（使用者已定，改动前必须确认）
----------------------------------
* **发生率分母 = 全部样本数**（含 STAT/GSA 本来就没做分类的样本）。
  「无数据」按未检出计 —— 最保守口径。`n_with_data` 一并返回，
  需要「有数据样本内发生率」时前端自行换算，不改变主口径。
* **默认口径「最特异检出」**：每样本只保留**没有检出后代**的病毒节点。

为什么必须去重
--------------
STAT 的 tax_table 给的是**完整谱系**，父子节点计数递增。实测同一 run：

    Betarhabdovirinae   total=307899   （亚科，含全部后代）
    ├─ Cytorhabdovirus  total= 88283
    │   └─ Cytorhabdovirus sp. 'lycii'  total= 88283   ← 最特异
    └─ Betacytorhabdovirus  total= 45155
        └─ Betacytorhabdovirus lycii  total= 45155
            └─ Goji cytorhabdovirus A  total= 45155    ← 最特异

按名字直接汇总，会把**一个物种同时算进它的属、亚科、科、目……**
发生率整体虚高。故默认只取叶子（无检出后代者）；
需要看某个层级时用 `level=<rank>` 显式切换，此时取该 rank 的节点各计一次。

踩过的三个坑（都有实测依据，别改回去）
--------------------------------------
1. **taxid 必须过 `merged.dmp` 归一。** STAT 的 tax_table 里同一个名字会挂在
   **多个新旧 taxid** 上（实测 `Caudoviricetes` → 10662 / 10744 / 2731619）。
   不归一，同一个节点会被拆成好几行（实测 831 行 vs 归一后 784 行）。
2. **占比要用 `total_count / analysed` 重算，不能用产物里的 `pct`。**
   产物 `pct` 只有 2 位小数，而病毒占比普遍是痕量级 —— 实测 121 个样本检出
   `Citrus exocitosis viroid`，`pct` 全是 `0.00`，直接用这一列就废了。
   （已实测：15643 行里 0 行与重算值不符。）
3. **名字反查优先用数据自建**（NCBI 行自带 name→taxid），只在剩余名字足够多时
   才通读 `names.dmp`。实测 1422 个名字里 1421 个能由数据自身覆盖，
   剩 1 个 —— 为 1 个名字通读 470 万行不值得。

性能
----
单次计算约 15~25 秒，**瓶颈是 `TaxonomyGraph.load()`（实测 13.4s，
要建 300 万条 nodes + 300 万条 names）**，与数据量无关。
故本模块做了两层缓存：
* 进程内 `_TG_CACHE`：同一个 tax_dir 只加载一次；
* 调用方（Web 层）另有按项目的磁盘缓存。
`build_virome_all()` 一次加载算完所有层级 —— 层级切换不要再重算。
"""
import csv
import os
import sys
from collections import Counter, defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
try:                                                      # 独立脚本 / --run-engine
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.taxonomy import TaxonomyGraph, taxonomy_ready
except ImportError:                                       # 作为包被 import 时
    from ..config import DIRS
    from ..taxonomy import TaxonomyGraph, taxonomy_ready

# 可选层级。'specific' = 最特异检出（默认）；其余为 taxonomy rank 名。
LEVELS = ('specific', 'species', 'genus', 'family', 'order', 'class',
          'phylum', 'kingdom', 'realm')

# names.dmp 里可当作「同一个东西」的名字类型（scientific name 优先）
_ALT_NAME_TYPES = ('equivalent name', 'synonym', 'genbank synonym', 'includes',
                   'genbank common name', 'common name')

# GSA path 的根是字面量 '{}'（空树根），不是分类单元
_PATH_ROOT_NOISE = ('{}', '', 'Viruses', 'root')

# 剩余名字少于此数就不通读 names.dmp（一次约 9 秒，为几个名字不值）
_RESIDUAL_SCAN_MIN = 20

# tax_dir -> 已 load() 的 TaxonomyGraph（进程内复用，避免每次 13 秒）
_TG_CACHE = {}


def _tax_graph(tax_dir=None):
    key = str(tax_dir or DIRS['taxonomy'])
    tg = _TG_CACHE.get(key)
    if tg is None:
        tg = TaxonomyGraph(key)
        tg.load()                # resolve() 自己不触发 load()，必须先显式加载
        _TG_CACHE[key] = tg
    return tg


# ------------------------------------------------------------------
# 读
# ------------------------------------------------------------------
def _read_tsv(path):
    if not os.path.isfile(path):
        return []
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def _to_int(v, default=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _to_float(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _precise_pct(row, denom):
    """用精确计数重算占比，而不是直接用产物里的 `pct`（只有 2 位小数）。

    分母口径与引擎一致：NCBI 用 analysed（已鉴定 reads），
    GSA 用 analysed（总 reads）—— 故跨来源仍不可直接比，需用 `virus_pct_ident`。
    """
    if denom and denom > 0:
        return round(100.0 * _to_int(row.get('total_count')) / denom, 6)
    return _to_float(row.get('pct'))


# ------------------------------------------------------------------
# 名字反查
# ------------------------------------------------------------------
def _scan_names_dmp(names, tax_dir=None):
    """{小写名: taxid}，scientific name 优先，异名兜底。只收 `names` 里的。"""
    want = {str(n).strip().lower() for n in names if str(n).strip()}
    if not want:
        return {}
    try:
        base = tax_dir or DIRS['taxonomy']
        p = os.path.join(base, 'names.dmp')
        if not os.path.isfile(p):
            return {}
        sci, alt = {}, {}
        with open(p, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                parts = line.split('|')
                if len(parts) < 4:
                    continue
                low = parts[1].strip().lower()
                if low not in want:
                    continue
                try:
                    tid = int(parts[0])
                except ValueError:
                    continue
                ntype = parts[3].strip()
                if ntype == 'scientific name':
                    sci.setdefault(low, tid)
                elif ntype in _ALT_NAME_TYPES:
                    alt.setdefault(low, tid)
        out = dict(alt)
        out.update(sci)                                   # 学名覆盖异名
        return out
    except OSError:
        return {}


def build_name_index(rows, tax_dir=None):
    """name(小写) -> taxid。**优先用数据自建**，再按需回退扫 names.dmp。

    NCBI 行自带 tax_id，所以「名字→taxid」大多能从数据本身推出来 ——
    实测 1422 个名字里 1421 个可自建，只剩 1 个要查库。
    为 1 个名字通读 470 万行 names.dmp（约 9 秒）不值，故设了阈值。
    """
    idx, residual = {}, set()
    for r in rows:
        org = str(r.get('org') or '').strip().lower()
        if not org:
            continue
        tid = str(r.get('tax_id') or '').strip()
        if tid:
            idx.setdefault(org, _to_int(tid))
        else:
            residual.add(org)
    residual -= set(idx)
    if len(residual) >= _RESIDUAL_SCAN_MIN:
        idx.update(_scan_names_dmp(residual, tax_dir))
    return idx


# ------------------------------------------------------------------
# 单行 → 统一表示
# ------------------------------------------------------------------
class _Ctx:
    """把一行 taxon 记录解析成可比较的 ident / 祖先集 / rank / 谱系。

    两种来源用不同的 ident 空间（NCBI 用 taxid，GSA 用名字小写），
    但**同一 acc 只会来自一个来源**（实测 311 个样本里 0 个多来源），
    所以单样本内父子比较始终在同一空间内。
    """

    def __init__(self, tax_dir=None, name_idx=None):
        self.tg = _tax_graph(tax_dir) if taxonomy_ready(tax_dir) else None
        self.name_idx = name_idx or {}
        self._lin = {}

    def lineage(self, taxid):
        if taxid not in self._lin:
            self._lin[taxid] = (self.tg.lineage(taxid) if self.tg else [])
        return self._lin[taxid]

    def rank_of(self, taxid):
        return self.tg.rank_of(taxid) if self.tg else ''

    def resolve(self, row):
        """→ dict(ident, ancestors, rank, taxid, lineage) 或 None（无信息）"""
        raw_tid = str(row.get('tax_id') or '').strip()
        org = str(row.get('org') or '').strip()
        path = str(row.get('path') or '').strip()

        if raw_tid:                                       # NCBI：taxid 为准
            taxid = _to_int(raw_tid, 0)
            if not taxid:
                return None
            # ⚠ 必须过 merged.dmp 归一（见模块 docstring 坑 1）
            rt = self.tg.resolve(taxid) if self.tg else taxid
            chain = self.lineage(rt)
            return {
                'ident': ('t', rt),
                'ancestors': {('t', t) for t, _r, _n in chain[1:]},
                'rank': self.rank_of(rt),
                'taxid': rt,
                'lineage': {r: n for t, r, n in chain},
            }

        if not org:                                       # GSA：只有名字
            return None
        parts = [p.strip() for p in path.split('>')] if path else []
        anc = {('n', p.lower()) for p in parts[:-1]
               if p and p not in _PATH_ROOT_NOISE}
        anc.discard(('n', org.lower()))
        taxid = self.name_idx.get(org.lower(), 0)
        chain = self.lineage(taxid) if taxid else []
        return {
            'ident': ('n', org.lower()),
            'ancestors': anc,
            'rank': self.rank_of(taxid) if taxid else '',
            'taxid': taxid,
            'lineage': {r: n for t, r, n in chain},
        }

    def canon(self, row, info):
        """跨样本聚合用的规范键 + 展示名。"""
        taxid = info.get('taxid') or 0
        if not taxid:                                     # 退回名字键
            taxid = self.name_idx.get(str(row.get('org') or '').lower(), 0)
        if taxid:
            return ('t', taxid), self._disp_name(taxid, row)
        org = str(row.get('org') or '')
        return ('n', org.lower()), org

    def _disp_name(self, taxid, row):
        """优先用 taxdump 里的学名；查不到就退回 STAT/GSA 自己给的名字。

        （STAT 的分类库与本地 taxdump 不同版本，实测有 7 个名字本地查不到 ——
        退回原始名比显示 `taxid:3403061` 有用得多。）
        """
        org = str(row.get('org') or '')
        if self.tg:
            n = self.tg.name(taxid)
            if n and not n.startswith('taxid:'):
                return n
        return org or f'taxid:{taxid}'


# ------------------------------------------------------------------
# 第一步：解析（贵，只做一次）
# ------------------------------------------------------------------
def _prepare(taxa_tsv, summary_tsv, tax_dir=None):
    summary = _read_tsv(summary_tsv)
    taxa = _read_tsv(taxa_tsv)
    denom_by_acc = {str(r.get('acc') or '').strip(): _to_int(r.get('analysed'))
                    for r in summary}

    ctx = _Ctx(tax_dir, build_name_index(taxa, tax_dir))

    by_acc = defaultdict(list)
    for r in taxa:
        acc = str(r.get('acc') or '').strip()
        if acc:
            by_acc[acc].append(r)

    rank_hist = Counter()
    rank_samples = defaultdict(set)
    per_sample = {}
    per_sample_n = {}

    for acc, rows in by_acc.items():
        parsed = []
        for r in rows:
            info = ctx.resolve(r)
            if not info:
                continue
            key, disp = ctx.canon(r, info)
            parsed.append({
                'key': key, 'name': disp,
                'rank': info['rank'], 'taxid': info['taxid'],
                'lineage': info['lineage'],
                'ident': info['ident'], 'ancestors': info['ancestors'],
                'source': str(r.get('source') or ''),
                'pct': _precise_pct(r, denom_by_acc.get(acc, 0)),
            })
            if info['rank']:
                rank_hist[info['rank']] += 1
                rank_samples[info['rank']].add(acc)

        if not parsed:
            continue

        # 叶子判定：该行不是任何**其它**检出行的祖先
        all_anc = set()
        for p in parsed:
            all_anc |= p['ancestors']
        for p in parsed:
            p['is_leaf'] = p['ident'] not in all_anc

        per_sample[acc] = parsed
        per_sample_n[acc] = len({p['key'] for p in parsed if p['is_leaf']})

    return {
        'n_total': len(summary),
        'n_with_data': sum(1 for r in summary
                           if (r.get('status') or '') == 'OK'),
        'per_sample': per_sample,
        'per_sample_n': per_sample_n,
        'rank_hist': rank_hist,
        'rank_samples': rank_samples,
    }


# ------------------------------------------------------------------
# 第二步：按层级聚合（便宜，可反复调）
# ------------------------------------------------------------------
def _aggregate(prep, level='specific', top_n=30, matrix_top=40,
               matrix_max_samples=400):
    if level not in LEVELS:
        level = 'specific'
    per_sample = prep['per_sample']
    per_sample_n = prep['per_sample_n']
    rank_hist = prep['rank_hist']
    rank_samples = prep['rank_samples']
    n_total = prep['n_total']

    # 每样本：选出该层级的行，按规范键去重（同一键取占比最大者）
    slots = {}
    for acc, parsed in per_sample.items():
        if level == 'specific':
            picked = [p for p in parsed if p['is_leaf']]
        else:
            picked = [p for p in parsed if p['rank'] == level]
        slot = {}
        for p in picked:
            cur = slot.get(p['key'])
            if cur is None or p['pct'] > cur['pct']:
                slot[p['key']] = p
        if slot:
            slots[acc] = slot

    # 跨样本聚合
    agg = {}
    for acc, slot in slots.items():
        for k, p in slot.items():
            a = agg.get(k)
            if a is None:
                a = agg[k] = {
                    'key': k,
                    'name': p['name'], 'rank': p['rank'], 'taxid': p['taxid'],
                    'lineage': p['lineage'],
                    'n_samples': 0, 'sum_pct': 0.0, 'max_pct': 0.0,
                    'sources': Counter(), 'samples': [],
                }
            a['n_samples'] += 1
            a['sum_pct'] += p['pct']
            a['max_pct'] = max(a['max_pct'], p['pct'])
            a['sources'][p['source']] += 1
            a['samples'].append(acc)
            if not a['rank'] and p['rank']:
                a['rank'] = p['rank']
            if not a['lineage'] and p['lineage']:
                a['lineage'] = p['lineage']
            if not a['taxid'] and p['taxid']:
                a['taxid'] = p['taxid']

    denom = n_total or 1
    ordered = sorted(agg.values(),
                     key=lambda a: (-a['n_samples'], -a['sum_pct'], a['name']))
    viruses = []
    for a in ordered:
        lin = a['lineage'] or {}
        viruses.append({
            'name': a['name'],
            'taxid': a['taxid'],
            'rank': a['rank'] or '未标注',
            'n_samples': a['n_samples'],
            'occ': round(100.0 * a['n_samples'] / denom, 6),
            'pct_mean': round(a['sum_pct'] / a['n_samples'], 6),
            'pct_max': round(a['max_pct'], 6),
            'pct_mean_all': round(a['sum_pct'] / denom, 6),
            'sources': dict(a['sources']),
            'realm': lin.get('realm', ''),
            'family': lin.get('family', ''),
            'genus': lin.get('genus', ''),
        })

    # 分类构成（按科）。
    # ⚠ `n_hits` 是「样本×病毒」对数，一个样本可能在同一科里有多个病毒，
    #   所以它 ≥ 去重样本数；两个都返回，别把 n_hits 当成样本数用。
    fam = {}
    fam_samples = defaultdict(set)
    for a in agg.values():
        f = (a['lineage'] or {}).get('family') or '未归类到科'
        d = fam.setdefault(f, {'family': f, 'n_hits': 0, 'viruses': 0,
                               'n_samples': 0})
        d['n_hits'] += a['n_samples']
        d['viruses'] += 1
        fam_samples[f].update(a['samples'])
    for f, s in fam_samples.items():
        fam[f]['n_samples'] = len(s)
    by_family = sorted(fam.values(),
                       key=lambda x: (-x['n_samples'], -x['viruses'], x['family']))

    # 热图矩阵：top N 病毒 × 有检出的样本（按检出数降序，密集的排前面）
    top = ordered[:matrix_top]
    hit_accs = sorted(slots, key=lambda a: (-per_sample_n.get(a, 0), a))
    hit_accs = hit_accs[:matrix_max_samples]
    z = [[round(slots.get(acc, {}).get(a['key'], {}).get('pct', 0.0), 6)
          for acc in hit_accs] for a in top]

    return {
        'level': level,
        'n_total': n_total,
        'n_with_data': prep['n_with_data'],
        'n_with_virus': len(slots),
        'n_distinct': len(viruses),
        'levels_available': [lv for lv in LEVELS
                             if lv == 'specific' or rank_hist.get(lv)],
        'rank_hist': [{'rank': r, 'rows': rank_hist[r],
                       'samples': len(rank_samples[r])}
                      for r in LEVELS if r != 'specific' and rank_hist.get(r)],
        'viruses': viruses,
        'top': viruses[:top_n],
        'by_family': by_family,
        'matrix': {'viruses': [a['name'] for a in top],
                   'keys': [list(a['key']) for a in top],
                   'samples': hit_accs, 'z': z},
        'richness': [{'n': k, 'samples': v}
                     for k, v in sorted(Counter(per_sample_n.values()).items())],
    }


# ------------------------------------------------------------------
# 对外入口
# ------------------------------------------------------------------
def build_virome(taxa_tsv, summary_tsv, level='specific', tax_dir=None,
                 top_n=30, matrix_top=40, matrix_max_samples=400):
    """单层级。CLI 与简单调用用；Web 层请用 `build_virome_all`。"""
    prep = _prepare(taxa_tsv, summary_tsv, tax_dir)
    return _aggregate(prep, level, top_n, matrix_top, matrix_max_samples)


def build_virome_all(taxa_tsv, summary_tsv, tax_dir=None, levels=LEVELS,
                     top_n=30, matrix_top=40, matrix_max_samples=400):
    """一次解析、算完全部层级 → {level: result}。

    层级切换不该重算 —— 贵的是 taxonomy 加载，与层级无关。
    """
    prep = _prepare(taxa_tsv, summary_tsv, tax_dir)
    avail = [lv for lv in levels
             if lv == 'specific' or prep['rank_hist'].get(lv)]
    return {lv: _aggregate(prep, lv, top_n, matrix_top, matrix_max_samples)
            for lv in avail}


# ------------------------------------------------------------------
# CLI（排查用）
# ------------------------------------------------------------------
def _fmt_pct(v):
    """占比普遍极小（痕量级），定点 3 位会把它们全抹成 0 —— 按量级换精度。"""
    v = float(v or 0)
    if v == 0:
        return '0'
    if v >= 1:
        return f'{v:.3f}'
    if v >= 0.001:
        return f'{v:.4f}'
    return f'{v:.2e}'


def _main(argv=None):
    import argparse
    import json
    p = argparse.ArgumentParser(description='预筛产物 → 病毒谱（只读，不联网）')
    p.add_argument('project_dir', help='<项目>/prescreen 目录，或项目目录')
    p.add_argument('--level', default='specific', choices=list(LEVELS))
    p.add_argument('--top', type=int, default=30)
    p.add_argument('--json', action='store_true', help='输出完整 JSON')
    a = p.parse_args(argv)

    d = a.project_dir
    if not os.path.isfile(os.path.join(d, 'stat_viral_taxa.tsv')):
        d = os.path.join(d, 'prescreen')
    res = build_virome(os.path.join(d, 'stat_viral_taxa.tsv'),
                       os.path.join(d, 'stat_run_summary.tsv'),
                       level=a.level, top_n=a.top)
    if a.json:
        print(json.dumps(res, ensure_ascii=False))
        return 0

    print(f"样本总数 {res['n_total']}（有分类数据 {res['n_with_data']}）"
          f" · 检出病毒样本 {res['n_with_virus']}"
          f" · 层级 {res['level']} · 去重病毒 {res['n_distinct']}")
    print(f"{'#':>4} {'病毒':<44}{'rank':<12}{'样本':>6}{'发生率%':>9}"
          f"{'均值%':>11}{'最大%':>11}")
    for i, v in enumerate(res['top'], 1):
        print(f"{i:>4} {v['name'][:43]:<44}{v['rank']:<12}{v['n_samples']:>6}"
              f"{v['occ']:>9.2f}{_fmt_pct(v['pct_mean']):>11}"
              f"{_fmt_pct(v['pct_max']):>11}")
    print()
    print('按科（前 10，n_hits 为「样本×病毒」对，n_samples 为去重样本数）:')
    for f in res['by_family'][:10]:
        print(f"   {f['family'][:40]:<42} 样本 {f['n_samples']:>5}"
              f"  命中 {f['n_hits']:>5}  病毒 {f['viruses']:>4}")
    return 0


if __name__ == '__main__':
    raise SystemExit(_main())
