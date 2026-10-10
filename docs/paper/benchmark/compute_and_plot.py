# -*- coding: utf-8 -*-
"""评测指标计算 + 中文图绘制。

输入: <runs>/results.json（run_bench.py 产物）
输出: 同目录 metrics_*.csv 与 figs/fig_b1..b5（PNG + SVG）
指标口径与参考文献一致（Massart 2019 性能测试 / PhytoPipe / Viroscope）:
  灵敏度 Sensitivity = TP/(TP+FN)
  精确率 Precision  = TP/(TP+FP)
  假发现率 FDR      = FP/(TP+FP)
  F1 = 2·P·S/(P+S)
"""
import csv
import glob
import json
import math
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
# 输入目录（results.json + 指标/性能 CSV 的输出目录）：
#   优先级 argv[1] > 环境变量 VP_BENCH_RUNS > 仓库内 results/（归档数据）。
#   旧默认 C:\vp_srna_test\bench 为评测机路径，2026-10-03 改为仓库内默认。
RUNS = (sys.argv[1] if len(sys.argv) > 1
        else os.environ.get('VP_BENCH_RUNS')
        or os.path.join(HERE, 'results'))
FIGS = os.path.join(HERE, 'figs')
os.makedirs(FIGS, exist_ok=True)

# D11–D18 真值：默认取仓库内 truth_viromock.json（与旧外部路径
# C:\vp_srna_test\bench_truth_viromock.json 逐字节相同，2026-10-03 核实）；
# 可用环境变量 VP_BENCH_TRUTH 覆盖。
VM_TRUTH = json.load(open(os.environ.get('VP_BENCH_TRUTH')
                          or os.path.join(HERE, 'truth_viromock.json'),
                          encoding='utf-8'))
VM_D110 = json.load(open(os.path.join(HERE, 'truth_viromock_d110.json'),
                         encoding='utf-8'))

# D1-D10 物种分组真值：每物种一名多别名（含 ICTV 改名/常见近缘写法）。
# 2026-10-03 起从脚本常量外置为数据文件 truth_viromock_d110_groups.json
# （审计核实与旧常量逐字相同）；与同目录 truth_viromock_d110.json（原始
# 解析、别名较简）物种集合一致，差异仅在别名完整度（本文件含 ICTV 改名
# 与常见近缘写法），评估采用本文件。可用环境变量 VP_BENCH_D110_GROUPS
# 覆盖。保留原注记：D9：PiVB 库内无该物种；Pistacia emaravirus 参考
# 100% 覆盖疑似误标，判 TP 但注记存疑。D7：TSWV 在评测快照库内完全缺失
# （真实数据库缺口，FN；2026-09-19 起 plant_ref 库已注入 TSWV）。
D110_GROUPS = json.load(open(os.environ.get('VP_BENCH_D110_GROUPS')
                             or os.path.join(HERE,
                                             'truth_viromock_d110_groups.json'),
                             encoding='utf-8'))

# Bari 真值（Massart et al. 2019 Phytopathology 109:488, Table 1，人工转录）。
# 已知缺口：S9 的 Grapevine virus A 在源表有检出记录，但本转录未含其
# reads 数（expected 缺该项）→ 定量点不含 GVA；S9 检出判定（species 列表）
# 仍含 GVA，不受影响。缺口已在 BENCHMARK_zh.md 数据说明节披露。
BARI_TRUTH = {
    'S3': {'species': ['Potato virus X', 'Potato virus B'],
           'expected': {'Potato virus X': 58204,
                        # RNA1 8,333 + RNA2 4,363（两段合并到物种级）
                        'Potato virus B': 12696}},
    'S9': {'species': ['Grapevine leafroll-associated virus 1', 'Grapevine virus A',
                       'Grapevine virus B', 'Grapevine rupestris stem pitting-associated virus',
                       'Grapevine red globe virus', 'Grapevine rupestris vein feathering virus',
                       'Grapevine Syrah virus 1', 'Hop stunt viroid',
                       'Grapevine yellow speckle viroid 1'],
           'expected': {'Grapevine leafroll-associated virus 1': 14292,
                        'Hop stunt viroid': 1181, 'Grapevine yellow speckle viroid 1': 341,
                        'Grapevine virus B': 2984, 'Grapevine rupestris stem pitting-associated virus': 1002,
                        'Grapevine red globe virus': 78, 'Grapevine rupestris vein feathering virus': 593,
                        'Grapevine Syrah virus 1': 60}},
}
# VIROMOCK 真值物种 → 平台库物种名匹配别名（规范化后双向包含）。
# 注意 ICTV 改名：Pepino mosaic virus 在库里亦作 Potexvirus pepini。
VM_ALIAS = {
    '11': ['pepino mosaic virus', 'potexvirus pepini'],
    '12': ['african cassava mosaic virus', 'cassava mosaic virus',
           'begomovirus manihotis'],
    '13': ['banana streak'],
    '14': ['potato virus y', 'potyvirus yituberosi'],
    '15': ['eggplant mottled dwarf', 'alphanucleorhabdovirus melongenae'],
    '16': ['bell pepper endornavirus', 'bell pepper alphaendornavirus'],
    '17': ['little cherry virus 1'],
    '18': ['barley yellow dwarf', 'luteovirus pavhordei'],
}


def norm(s):
    return re.sub(r'[^a-z0-9 ]', '', (s or '').lower()).strip()


def species_match(reported, aliases):
    r = norm(reported)
    return any(a in r or r in a for a in map(norm, aliases) if a)


def load():
    rows = json.load(open(os.path.join(RUNS, 'results.json'), encoding='utf-8'))
    return {r['exp']: r for r in rows}


def kv_species(res):
    return [(c['species'], c) for c in res.get('confirmed', [])]


def eval_agent_set(confirmed, groups):
    """物种级判定（分组别名版）：返回 (tp, fp_names, fn)。

    groups: [{'name': 真值物种, 'al': [可匹配别名(规范化前任意写法)]}]
    TP 按真值物种去重（含 ICTV 改名：命中任一别名即该物种 TP）；
    fp_names 为去重后的假阳性物种名。
    """
    hit, fp_names = set(), {}
    for sp, _ in confirmed:
        nsp = norm(sp)
        matched = None
        for g in groups:
            for a in [norm(g['name'])] + [norm(x) for x in g.get('al', [])]:
                if a and (a in nsp or nsp in a):
                    matched = g['name']
                    break
            if matched:
                break
        if matched:
            hit.add(matched)
        else:
            fp_names[sp] = fp_names.get(sp, 0) + 1
    tp = len(hit)
    fn = len(groups) - tp
    return tp, sorted(fp_names, key=fp_names.get, reverse=True), max(fn, 0)


def prf(tp, fp, fn):
    s = tp / (tp + fn) if tp + fn else 0.0
    p = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * p * s / (p + s) if p + s else 0.0
    fdr = fp / (tp + fp) if tp + fp else 0.0
    return {'sens': s, 'prec': p, 'fdr': fdr, 'f1': f1,
            'tp': tp, 'fp': fp, 'fn': fn}


def main():
    R = load()
    detail, summary, perf, quant = [], [], [], []

    # —— 作废运行判定（数据驱动）——
    # minibwa 引擎记录应带 mapped/total（平台运行日志的 mapped=X/Y 汇总行）。
    # 缺失 = 运行异常提前退出（本轮 = vm_d5_bwa：0.9s / 2MB，无 mapped/total，
    # 监测失真），该记录作废、不计入任何指标与图；记录仍在 perf 表中保留并
    # 标注（见 BENCHMARK_zh.md §4）。判定条件与数据字段绑定：若记录补齐
    # mapped/total 则自动恢复为有效。
    void = {}
    for k, res in R.items():
        if (k.endswith('_bwa') and res.get('confirmed') is not None
                and res.get('mapped') is None and res.get('total') is None):
            void[k] = 'minibwa 运行异常提前退出（无 mapped/total，监测失真），记录作废'
            print(f'[warn] {k}：{void[k]}，不计入指标')

    # —— VIROMOCK（salmon k31 常规路线，物种级）——
    for d, t in VM_TRUTH.items():
        if str(d).startswith('_'):          # 真值文件的 _notes 元数据段
            continue
        res = R.get(f'vm_d{d}')
        if not res or 'confirmed' not in res:
            continue
        conf = kv_species(res)
        groups = [{'name': t['species'], 'al': VM_ALIAS[d]}]
        tp, fp_names, fn = eval_agent_set(conf, groups)
        m = prf(tp, len(fp_names), fn)
        exp_reads = sum(x['reads'] for x in t['isolates'])
        obs = sum(v for sp, v in (res.get('species_em') or {}).items()
                  if species_match(sp, VM_ALIAS[d]))
        summary.append({'组': 'VIROMOCK-8', '方法': 'salmon-k31(常规)',
                        '数据集': f"D{d} {t['species']}", **m,
                        'FP物种': '; '.join(fp_names[:3])})
        detail.append({'dataset': f'D{d}', 'truth': t['species'], 'tp': tp,
                       'fp': fp_names, 'fn': fn, 'expected_reads': exp_reads,
                       'observed_em': obs,
                       'mapped': res.get('mapped'), 'total': res.get('total')})
        quant.append({'点': f"D{d} {t['species'][:14]}", '预期': exp_reads,
                      '观测': obs, '组': 'VIROMOCK'})
        perf.append({'实验': f"VIROMOCK D{d}", '方法': 'salmon-k31',
                     'wall_s': res.get('wall_s'), 'rss_mb': res.get('peak_rss_mb'),
                     'mapped': res.get('mapped'), 'total': res.get('total'),
                     'reads_in': (res.get('total') or 0) or None})

    # —— Bari sRNA：鉴定（k15 适配）——
    for s in ('S3', 'S9'):
        res = R.get(f'bari_{s.lower()}_k15')
        if not res or 'confirmed' not in res:
            continue
        conf = kv_species(res)
        groups = [{'name': t, 'al': [t]} for t in BARI_TRUTH[s]['species']]
        tp, fp_names, fn = eval_agent_set(conf, groups)
        m = prf(tp, len(fp_names), fn)
        summary.append({'组': f'Bari-{s}', '方法': 'salmon-k15(小RNA适配)',
                        '数据集': f'{s} (250K sRNA)', **m,
                        'FP物种': '; '.join(fp_names[:3])})
        em_by_sp = res.get('species_em') or {}
        for hit, exp in BARI_TRUTH[s]['expected'].items():
            obs = sum(v for sp, v in em_by_sp.items()
                      if norm(sp) == norm(hit) or norm(hit) in norm(sp))
            quant.append({'点': f"{s} {hit[:16]}", '预期': exp,
                          '观测': obs, '组': 'Bari-sRNA'})
        perf.append({'实验': f'Bari {s}', '方法': 'salmon-k15',
                     'wall_s': res.get('wall_s'), 'rss_mb': res.get('peak_rss_mb'),
                     'mapped': res.get('mapped'), 'total': res.get('total'),
                     'reads_in': res.get('total'),
                     'index_build_s': res.get('index_build_s')})

    # —— Bari sRNA：组装（srna 模式 + BLAST ≥60bp）——
    for s in ('S3', 'S9'):
        res = R.get(f'bari_{s.lower()}_asm')
        if not res:
            continue
        sp_hits = res.get('species_contigs', {})
        conf = [(sp, None) for sp in sp_hits]
        groups = [{'name': t, 'al': [t]} for t in BARI_TRUTH[s]['species']]
        tp, fp_names, fn = eval_agent_set(conf, groups)
        m = prf(tp, len(fp_names), fn)
        summary.append({'组': f'Bari-{s}', '方法': 'srna组装+BLAST(≥60bp)',
                        '数据集': f'{s} (250K sRNA)', **m,
                        'FP物种': '; '.join(fp_names[:3])})
        perf.append({'实验': f'Bari {s}', '方法': 'srna组装',
                     'wall_s': res.get('wall_s'), 'rss_mb': res.get('peak_rss_mb'),
                     'mapped': None, 'total': None,
                     'reads_in': 250000, 'n_contigs': res.get('n_contigs'),
                     'n_ge60': res.get('n_ge60')})

    # —— VIROMOCK D1-D10 半人工/真实集（salmon + minibwa 双引擎，物种级）——
    for d, t in VM_D110.items():
        for eng, tag in (('salmon', ''), ('minibwa', '_bwa')):
            exp = f'vm_d{d}{tag}'
            res = R.get(exp)
            if not res or 'confirmed' not in res:
                continue
            if exp in void:      # 作废记录：只保留 perf 行（带标注），不进指标
                perf.append({'实验': f'VIROMOCK D{d}', '方法': f'{eng}',
                             'wall_s': res.get('wall_s'),
                             'rss_mb': res.get('peak_rss_mb'),
                             'mapped': res.get('mapped'),
                             'total': res.get('total'),
                             'reads_in': (res.get('total') or 0) or None,
                             '备注': void[exp]})
                continue
            conf = kv_species(res)
            groups = D110_GROUPS[d]
            # P1-8：真值标 ambiguous 的组不进主指标（TP/FN 受不确定数据库
            # 注释驱动）；另出一行敏感性口径（计入 ambiguous）供对比。
            main_groups = [g for g in groups if not g.get('ambiguous')]
            amb_groups = [g for g in groups if g.get('ambiguous')]
            tp, fp_names, fn = eval_agent_set(conf, main_groups)
            m = prf(tp, len(fp_names), fn)
            summary.append({'组': 'VIROMOCK-D110', '方法': f'{eng}引擎',
                            '数据集': f"D{d} {t['species'][0]}", **m,
                            'FP物种': '; '.join(fp_names[:3])})
            if amb_groups:
                tp_a, fp_a, fn_a = eval_agent_set(conf, main_groups
                                                   + amb_groups)
                m_a = prf(tp_a, len(fp_a), fn_a)
                summary.append({'组': 'VIROMOCK-D110(敏感性)',
                                '方法': f'{eng}引擎+存疑真值',
                                '数据集': f"D{d} 含ambiguous",
                                **m_a, 'FP物种': '; '.join(fp_a[:3]),
                                '备注': '含 D9 存疑别名真值（详见'
                                        'BENCHMARK_zh.md），仅作敏感性对照'})
            obs = sum(v for sp, v in (res.get('species_em') or {}).items()
                      if species_match(sp, [a for g in groups for a in
                                            [g['name']] + g.get('al', [])]))
            perf.append({'实验': f'VIROMOCK D{d}', '方法': f'{eng}',
                         'wall_s': res.get('wall_s'),
                         'rss_mb': res.get('peak_rss_mb'),
                         'mapped': res.get('mapped'),
                         'total': res.get('total'),
                         'reads_in': (res.get('total') or 0) or None})

    # —— D11-D18 minibwa 引擎 ——
    for d, t in VM_TRUTH.items():
        if str(d).startswith('_'):
            continue
        res = R.get(f'vm_d{d}_bwa')
        if not res or 'confirmed' not in res or f'vm_d{d}_bwa' in void:
            continue
        conf = kv_species(res)
        groups = [{'name': t['species'], 'al': VM_ALIAS[d]}]
        tp, fp_names, fn = eval_agent_set(conf, groups)
        m = prf(tp, len(fp_names), fn)
        summary.append({'组': 'VIROMOCK-8', '方法': 'minibwa引擎',
                        '数据集': f"D{d} {t['species']}", **m,
                        'FP物种': '; '.join(fp_names[:3])})
        perf.append({'实验': f'VIROMOCK D{d}', '方法': 'minibwa',
                     'wall_s': res.get('wall_s'),
                     'rss_mb': res.get('peak_rss_mb'),
                     'mapped': res.get('mapped'), 'total': res.get('total'),
                     'reads_in': (res.get('total') or 0) or None})

    # —— Bari 50K 深度稀释 ——
    for s in ('S3', 'S9'):
        res = R.get(f'bari_{s.lower()}_50k')
        if not res or 'confirmed' not in res:
            continue
        conf = kv_species(res)
        groups = [{'name': t, 'al': [t]} for t in BARI_TRUTH[s]['species']]
        tp, fp_names, fn = eval_agent_set(conf, groups)
        m = prf(tp, len(fp_names), fn)
        summary.append({'组': f'Bari-{s}', '方法': 'salmon-k15(50K稀释)',
                        '数据集': f'{s} (50K sRNA)', **m,
                        'FP物种': '; '.join(fp_names[:3])})
        perf.append({'实验': f'Bari {s}', '方法': 'salmon-k15·50K',
                     'wall_s': res.get('wall_s'), 'rss_mb': res.get('peak_rss_mb'),
                     'mapped': res.get('mapped'), 'total': res.get('total'),
                     'reads_in': res.get('total')})

    # —— Bari 全链路（fastp→鉴定→组装）——
    for s in ('S3', 'S9'):
        res = R.get(f'chain_{s.lower()}')
        if not res:
            continue
        conf = kv_species(res) + [(sp, None) for sp in
                                  (res.get('species_contigs') or {})]
        groups = [{'name': t, 'al': [t]} for t in BARI_TRUTH[s]['species']]
        tp, fp_names, fn = eval_agent_set(conf, groups)
        m = prf(tp, len(fp_names), fn)
        summary.append({'组': f'Bari-{s}', '方法': '全链路(fastp→鉴定+组装)',
                        '数据集': f'{s} (250K sRNA)', **m,
                        'FP物种': '; '.join(fp_names[:3])})
        st = res.get('stages') or {}
        perf.append({'实验': f'Bari {s}', '方法': '全链路',
                     'wall_s': res.get('driver_wall_s'),
                     'rss_mb': res.get('peak_rss_mb'),
                     'mapped': res.get('mapped'), 'total': res.get('total'),
                     'reads_in': res.get('total'),
                     'stages': st, 'qc_survived': res.get('qc_survived_reads')})

    # —— 对照组：Bari 小RNA 上的引擎原始配置（S3 + S9；数据取自 results.json）——
    # salmon-k31 对照：输入 S3/S9 250K 21–24nt sRNA；所用索引为 3 参考小库
    #   idx_k31（NC_011620.1 PVX + NC_043447.1/NC_043448.1 PVB，num_refs=3，
    #   2026-09-15 构建），非平台全库——已在 CSV 备注与 BENCHMARK_zh.md 披露。
    # minibwa 对照：平台入口对 21–24nt 明确拒绝（rejected=true，含 message）。
    # 指标不再写死：confirmed 为空时按 0 命中计（TP0/FP0/FN=真值物种数），
    # 若记录将来带 confirmed 则自动按其评估。
    def _ctrl_rows(exp, name, method, extra_note):
        rec = R.get(exp) or {}
        n_t = len(BARI_TRUTH[name]['species'])
        conf = [(c['species'], c) for c in (rec.get('confirmed') or [])]
        if conf:
            groups = [{'name': t, 'al': [t]}
                      for t in BARI_TRUTH[name]['species']]
            tp, fp_names, fn = eval_agent_set(conf, groups)
            m = prf(tp, len(fp_names), fn)
        else:
            m = prf(0, 0, n_t)
        summary.append({'组': f'Bari-{name}', '方法': method,
                        '数据集': f'{name} (250K sRNA)', **m,
                        'FP物种': '', '备注': extra_note})

    for exp, name in (('ctrl_s3_k31', 'S3'), ('ctrl_s9_k31', 'S9')):
        rec = R.get(exp) or {}
        note = (f"对照：3 参考小库 idx_k31（PVX+PVB，num_refs=3）；"
                f"num_mapped={rec.get('num_mapped')}，"
                f"percent_mapped={rec.get('percent_mapped')}")
        _ctrl_rows(exp, name, 'salmon-k31(原始配置,对照)', note)
        perf.append({'实验': f'Bari {name}', '方法': 'salmon-k31对照',
                     'wall_s': rec.get('wall_s'),
                     'rss_mb': rec.get('peak_rss_mb'),
                     'mapped': rec.get('num_mapped'), 'total': 250000,
                     'reads_in': 250000, '备注': note})
    for exp, name in (('ctrl_s3_bwa', 'S3'), ('ctrl_s9_bwa', 'S9')):
        rec = R.get(exp) or {}
        msg = (rec.get('message') or '').strip()
        note = ('对照：minibwa 引擎'
                + ('（rejected=true，平台拒绝 21–24nt reads）'
                   if rec.get('rejected') else '')
                + (f'；message：{msg}' if msg else ''))
        _ctrl_rows(exp, name, 'minibwa(对照)', note)
        perf.append({'实验': f'Bari {name}', '方法': 'minibwa对照',
                     'wall_s': rec.get('wall_s'),
                     'rss_mb': rec.get('peak_rss_mb'),
                     'mapped': None, 'total': 250000,
                     'reads_in': 250000, '备注': note})

    # —— 写表 ——
    def write_csv(name, rows):
        # 列 = 首行键序 + 后续行新增键（备注等可选列不丢字段）
        cols = list(rows[0].keys()) if rows else []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        p = os.path.join(RUNS, name)
        with open(p, 'w', newline='', encoding='utf-8-sig') as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        print('written', p)

    write_csv('metrics_summary.csv', summary)
    write_csv('perf.csv', perf)

    # —— 聚合（画图用） ——
    agg = {}
    for r in summary:
        key = (r['组'], r['方法'])
        a = agg.setdefault(key, {'tp': 0, 'fp': 0, 'fn': 0})
        for k in a:
            a[k] += r[k]
    def _wilson(k, n, z=1.96):
        """Wilson 95% CI（比例的小样本友好区间；审查 P2-3）。"""
        if not n:
            return None, None
        p_ = k / n
        d = 1 + z * z / n
        c = (p_ + z * z / (2 * n)) / d
        h = z * math.sqrt(p_ * (1 - p_) / n + z * z / (4 * n * n)) / d
        return round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)

    agg_rows = []
    for (grp, meth), a in agg.items():
        m = prf(a['tp'], a['fp'], a['fn'])
        sl, sh = _wilson(a['tp'], a['tp'] + a['fn'])
        pl, ph = _wilson(a['tp'], a['tp'] + a['fp'])
        agg_rows.append({'组': grp, '方法': meth, **m,
                         'sens_lo95': sl, 'sens_hi95': sh,
                         'prec_lo95': pl, 'prec_hi95': ph,
                         'n_truth': a['tp'] + a['fn'], 'n_call': a['tp'] + a['fp']})
    # 聚合行（含 Wilson 95% CI）导出成表：点估计与区间同表可见（P2-3）
    write_csv('metrics_agg_ci.csv', agg_rows)

    draw(summary, agg_rows, detail, quant, perf)
    return summary, agg_rows, perf, quant


def draw(summary, agg_rows, detail, quant, perf):
    C = {'ok': '#5B9A57', 'bad': '#C0504D', 'mid': '#4C72B0', 'warn': '#D9A441',
         'gray': '#AAAAAA'}

    # ── 图1：Bari sRNA 检出灵敏度（方法对比 + 论文参考线）──
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    methods = [('salmon-k15(小RNA适配)', '#4C72B0'),
               ('srna组装+BLAST(≥60bp)', '#5B9A57'),
               ('salmon-k31(原始配置,对照)', '#C0504D'),
               ('minibwa(对照)', '#AAAAAA')]
    labels, vals = [], []
    for meth, _ in methods:
        sub = [r for r in summary if r['方法'] == meth]
        # 合并（micro）口径：柱高与柱顶 TP/总数标注同源，均按 TP/(TP+FN) 汇总，
        # 不再用"样本均值"（与标注口径矛盾，2026-10-03 修正）。
        tp = sum(r['tp'] for r in sub)
        fn = sum(r['fn'] for r in sub)
        tot = tp + fn
        sens = tp / tot if tot else 0.0
        labels.append(meth)
        vals.append((sens, tp, tot))
    xs = range(len(labels))
    bars = ax.bar(xs, [v[0] * 100 for v in vals],
                  color=[m[1] for m in methods], width=0.62)
    for x, (sens, tp, tot) in zip(xs, vals):
        ax.text(x, sens * 100 + 2, f'{tp}/{tot}', ha='center', fontsize=10,
                weight='bold')
    ax.axhline(76, ls='--', lw=1.2, color='#D9A441')
    ax.text(len(labels) - 0.45, 77.5, 'Massart 2019：21 实验室 250K 均值 76%',
            ha='right', fontsize=8, color='#8a6414')
    ax.axhline(100, ls=':', lw=1.0, color='#888888')
    ax.text(len(labels) - 0.45, 93, '最佳管线 100%', ha='right', fontsize=8,
            color='#666666')
    ax.set_xticks(list(xs))
    ax.set_xticklabels(['salmon k15\n(小RNA适配)', 'srna 组装\n+BLAST≥60bp',
                        'salmon k31\n(原始配置)', 'minibwa\n(对照)'], fontsize=9)
    ax.set_ylabel('物种级检出灵敏度 (%)')
    ax.set_ylim(0, 108)
    ax.set_title('Bari 小RNA数据（2 样品 × 11 真值病原）检出性能', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'fig_b1_detection.png'), dpi=200)
    fig.savefig(os.path.join(FIGS, 'fig_b1_detection.svg'))
    plt.close(fig)

    # ── 图2：阳性判定汇总（TP/FP/FN 堆叠 + 指标）──
    _SHORT = {'salmon-k15(小RNA适配)': 'salmon-k15',
              'srna组装+BLAST(≥60bp)': 'srna组装',
              'salmon-k31(原始配置,对照)': 'salmon-k31对照',
              'minibwa(对照)': 'minibwa对照',
              'salmon-k15(50K稀释)': 'salmon-k15·50K',
              '全链路(fastp→鉴定+组装)': '全链路',
              'salmon引擎': 'salmon', 'minibwa引擎': 'minibwa',
              'salmon-k31(常规)': 'salmon-k31'}

    def _short(meth):
        return _SHORT.get(meth, meth.split('(')[0])

    rows = [r for r in agg_rows if r['tp'] + r['fp'] + r['fn'] > 0]
    rows.sort(key=lambda r: -(r['tp'] + r['fp'] + r['fn']))
    # 行多时按行数自适应加高，行标签单行短名（此前两行标签字号 8 且行距固定，
    # 图 b2 左图 y 轴标签重叠；2026-10-03 修正）
    fig_h2 = max(4.6, 0.30 * len(rows))
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(12.6, fig_h2),
                                   gridspec_kw={'width_ratios': [1.05, 1]})
    ys = range(len(rows))
    axl.barh(ys, [r['tp'] for r in rows], color=C['ok'], label='真阳性 TP')
    axl.barh(ys, [r['fp'] for r in rows], left=[r['tp'] for r in rows],
             color=C['bad'], label='假阳性 FP')
    axl.barh(ys, [r['fn'] for r in rows],
             left=[r['tp'] + r['fp'] for r in rows], color=C['warn'],
             label='假阴性 FN')
    for y, r in zip(ys, rows):
        axl.text(r['tp'] + r['fp'] + r['fn'] + 0.3, y,
                 f"TP{r['tp']} FP{r['fp']} FN{r['fn']}", va='center', fontsize=7.5)
    axl.set_yticks(list(ys))
    axl.set_yticklabels([f"{r['组']}·{_short(r['方法'])}" for r in rows],
                        fontsize=7.5)
    axl.set_xlabel('判定数（个）')
    axl.set_title('阳性/阴性判定构成', fontsize=11)
    axl.legend(fontsize=8, loc='upper right')
    axl.spines[['top', 'right']].set_visible(False)

    xs = range(len(rows))
    w = 0.27
    axr.bar([x - w for x in xs], [r['sens'] * 100 for r in rows], w,
            color=C['mid'], label='灵敏度')
    axr.bar(list(xs), [r['prec'] * 100 for r in rows], w,
            color=C['ok'], label='精确率')
    axr.bar([x + w for x in xs], [r['f1'] * 100 for r in rows], w,
            color=C['warn'], label='F1')
    axr.set_xticks(list(xs))
    # 完整标签（不再 [:10] 截断），旋转 45° 容纳
    axr.set_xticklabels([f"{r['组']}·{_short(r['方法'])}" for r in rows],
                        rotation=45, ha='right', fontsize=7)
    axr.set_ylabel('%')
    axr.set_ylim(0, 112)
    axr.set_title('灵敏度 / 精确率 / F1', fontsize=11)
    axr.legend(fontsize=8)
    axr.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'fig_b2_confusion.png'), dpi=200)
    fig.savefig(os.path.join(FIGS, 'fig_b2_confusion.svg'))
    plt.close(fig)

    # ── 图3：VIROMOCK 检出与定量（期望 vs 观测 reads）──
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    det = sorted(detail, key=lambda d: d['dataset'])
    xs = range(len(det))
    w = 0.38
    ax.bar([x - w / 2 for x in xs], [d['expected_reads'] for d in det], w,
           color='#B9CCE4', label='期望人工病毒 reads（真值）')
    ax.bar([x + w / 2 for x in xs], [d['observed_em'] or 0 for d in det], w,
           color=C['mid'], label='平台定量 EM reads（salmon k31）')
    ymax = max(d['expected_reads'] for d in det)
    for x, d in zip(xs, det):
        st = '√ 已确诊' if d['tp'] else '× 未确诊'
        y = max(d['expected_reads'], d['observed_em'] or 0)
        ax.text(x, y * 1.03, st, ha='center', fontsize=9,
                color=C['ok'] if d['tp'] else C['bad'], weight='bold')
    ax.set_xticks(list(xs))
    ax.set_xticklabels([f"D{d['dataset'][1:]}·{d['truth'].replace('*', '')[:10]}"
                        for d in det], fontsize=7.5, rotation=18, ha='right')
    ax.set_ylabel('reads 数')
    ax.set_title('VIROMOCK 人工数据集（8 种病毒）检出与定量'
                 '（柱顶 √/× = 物种是否确诊）', fontsize=10.5)
    ax.legend(fontsize=8.5, ncol=1, loc='upper right', framealpha=0.95)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylim(0, ymax * 1.30)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'fig_b3_viromock.png'), dpi=200)
    fig.savefig(os.path.join(FIGS, 'fig_b3_viromock.svg'))
    plt.close(fig)

    # ── 图4：定量准确性散点（预期 vs 观测，log-log）──
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    pts = [(q['预期'], q['观测'], q['组']) for q in quant if q['预期'] > 0]
    for grp, color, mk in (('VIROMOCK', C['mid'], 'o'), ('Bari-sRNA', C['ok'], 's')):
        px = [p[0] for p in pts if p[2] == grp]
        py = [p[1] for p in pts if p[2] == grp]
        ax.scatter(px, py, c=color, s=52, label=grp, alpha=0.85, marker=mk,
                   edgecolors='white', linewidths=0.6)
    lim = [30, max(max(p[0] for p in pts), max(p[1] for p in pts)) * 2]
    ax.plot(lim, lim, '--', color='#999999', lw=1, label='y = x（完全一致）')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('预期病毒 reads（真值）')
    ax.set_ylabel('平台定量 reads（EM）')
    ok = [(a, b) for a, b, _ in pts if b > 0]
    if len(ok) >= 3:
        import math
        xs = [math.log10(a) for a, _ in ok]
        ys = [math.log10(b) for _, b in ok]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        r = (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
             / math.sqrt(sum((x - mx) ** 2 for x in xs)
                         * sum((y - my) ** 2 for y in ys)))
        ax.text(0.05, 0.93, f'Pearson r = {r:.3f}（log 空间, n={len(ok)}）',
                transform=ax.transAxes, fontsize=9)
    # 口径注记（2026-10-03）：真值按 reads 计，PE 集 reads=2×fragments，
    # 平台 EM 计的是片段 → 整体位于 y=x 下方（系统偏移）；r 只表征相关、
    # 不表征准确度。详见 BENCHMARK_zh.md §3 定量小节。
    ax.text(0.05, 0.735,
            '注：PE 数据集真值为 reads（=2×片段）\n观测为平台 EM 片段计数，'
            '存在系统偏移\nr 仅表征相关强度，不表征准确度\n'
            '（后验核算，见 BENCHMARK_zh.md §3）',
            transform=ax.transAxes, fontsize=6.8, color='#555555', va='top')
    ax.set_title('定量准确性：预期 vs 观测', fontsize=11)
    ax.legend(fontsize=8, loc='lower right')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'fig_b4_quant.png'), dpi=200)
    fig.savefig(os.path.join(FIGS, 'fig_b4_quant.svg'))
    plt.close(fig)

    # ── 图5：耗时 + 峰值内存 ──
    perf_ok = [p for p in perf if p.get('wall_s') is not None]
    fig_h = max(4.4, 0.33 * len(perf_ok))
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, fig_h))
    ys = range(len(perf_ok))
    labels = [f"{p['实验']}\n{p['方法']}" for p in perf_ok]
    colors = ['#4C72B0' if 'salmon' in p['方法'] else
              '#5B9A57' if '组装' in p['方法'] else '#AAAAAA' for p in perf_ok]
    axl.barh(ys, [p['wall_s'] for p in perf_ok], color=colors)
    for y, p in zip(ys, perf_ok):
        extra = (f"（含索引 {p['index_build_s']:.0f}s）"
                 if p.get('index_build_s') else '')
        axl.text(p['wall_s'] + 1, y, f"{p['wall_s']:.0f}s{extra}",
                 va='center', fontsize=7.5)
    axl.set_yticks(list(ys))
    axl.set_yticklabels(labels, fontsize=7.5)
    axl.invert_yaxis()
    axl.set_xlabel('墙钟时间 (s)')
    axl.set_title('单样品耗时', fontsize=11)
    axl.spines[['top', 'right']].set_visible(False)

    rss = [p.get('rss_mb') or 0 for p in perf_ok]
    axr.barh(ys, rss, color=colors)
    for y, v in zip(ys, rss):
        axr.text(v + max(rss) * 0.01, y, f'{v:,.0f}', va='center', fontsize=7.5)
    axr.set_yticks(list(ys))
    axr.set_yticklabels([])
    axr.invert_yaxis()
    axr.set_xlabel('峰值内存 (MB，进程树合计)')
    axr.set_title('峰值内存', fontsize=11)
    axr.spines[['top', 'right']].set_visible(False)
    fig.suptitle('运行开销（4 线程，本机 Windows）', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'fig_b5_perf.png'), dpi=200)
    fig.savefig(os.path.join(FIGS, 'fig_b5_perf.svg'))
    plt.close(fig)

    # ── 图6：salmon vs minibwa 双引擎对比（半人工集 D1-D10）──
    # 数据集集合从真值文件键名解析（D1..D10 精确匹配），此前用
    # `'D1' in p['实验']` 子串过滤只匹配到 D1/D10-D18，丢失 D2-D9 的 8 条
    # minibwa 柱（2026-10-03 修正）。
    ds_order = [f'VIROMOCK D{d}' for d in sorted(VM_D110, key=int)]
    _ds_set = set(ds_order)
    sa = {p['实验']: p for p in perf
          if p['方法'] == 'salmon' and p['实验'] in _ds_set}
    mb = {p['实验']: p for p in perf
          if p['方法'] == 'minibwa' and p['实验'] in _ds_set}
    sa_ts = [sa.get(k, {}).get('wall_s') for k in ds_order]
    mb_ts = [mb.get(k, {}).get('wall_s') for k in ds_order]
    mb_notes = [mb.get(k, {}).get('备注') for k in ds_order]
    if any(v is not None for v in mb_ts):
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, 4.4))
        xs = range(len(ds_order))
        w = 0.38
        axl.bar([x - w / 2 for x in xs], [t or 0 for t in sa_ts], w,
                color=C['mid'], label='salmon（伪比对 EM 定量）')
        axl.bar([x + w / 2 for x in xs], [t or 0 for t in mb_ts], w,
                color=C['warn'], label='minibwa（真比对计数）')
        for x, (a, b) in zip(xs, zip(sa_ts, mb_ts)):
            for v, off in ((a, -w / 2), (b, w / 2)):
                if v:
                    axl.text(x + off, v * 1.05, f'{v:.0f}', ha='center',
                             fontsize=7.5)
        # 缺失/作废记录标注由数据字段驱动（记录缺 wall_s 或带"作废"备注），
        # 不再写死"超时"字符串（results.json 无 timed_out 字段；2026-10-03 修正）
        for x, note in zip(xs, mb_notes):
            if note and '作废' in note:
                axl.text(x + w / 2, 12, '作废', ha='center', fontsize=7,
                         color=C['bad'], rotation=90)
        axl.set_yscale('log')
        axl.set_xticks(list(xs))
        axl.set_xticklabels([k.replace('VIROMOCK ', '') for k in ds_order],
                            fontsize=8)
        axl.set_ylabel('墙钟时间 (s，log)')
        axl.set_title('单样品耗时：双引擎（4 线程）', fontsize=11)
        axl.legend(fontsize=8.5)
        axl.spines[['top', 'right']].set_visible(False)

        eng_rows = [r for r in summary if r['组'] == 'VIROMOCK-D110']
        labs, s_vals, m_vals = [], [], []
        for d in sorted(VM_D110, key=int):
            sub = {r['方法']: r for r in eng_rows if r['数据集'].startswith(f'D{d} ')}
            labs.append(f'D{d}')
            s_vals.append(sub.get('salmon引擎', {}).get('sens', 0) * 100)
            # D5 minibwa 记录作废 → 无 minibwa 指标：留空（nan）并标注
            if 'minibwa引擎' in sub:
                m_vals.append(sub['minibwa引擎']['sens'] * 100)
            else:
                m_vals.append(float('nan'))
        xs2 = range(len(labs))
        axr.bar([x - w / 2 for x in xs2], s_vals, w, color=C['mid'],
                label='salmon')
        axr.bar([x + w / 2 for x in xs2], m_vals, w, color=C['warn'],
                label='minibwa')
        for x, v in zip(xs2, m_vals):
            if v != v:      # nan：该集 minibwa 记录作废
                axr.text(x + w / 2, 3, '作废', ha='center', fontsize=6.5,
                         color=C['bad'], rotation=90)
        axr.axhline(100, ls=':', lw=0.8, color='#999999')
        axr.set_xticks(list(xs2))
        axr.set_xticklabels(labs, fontsize=8)
        axr.set_ylabel('物种级检出灵敏度 (%)')
        axr.set_ylim(0, 112)
        axr.set_title('半人工集 D1-D10 检出灵敏度', fontsize=11)
        axr.legend(fontsize=8.5)
        axr.spines[['top', 'right']].set_visible(False)
        fig.suptitle('双引擎对比（VIROMOCK D1-D10 半人工/真实集）', fontsize=11.5)
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, 'fig_b6_engine.png'), dpi=200)
        fig.savefig(os.path.join(FIGS, 'fig_b6_engine.svg'))
        plt.close(fig)

    # ── 图7：测序深度效应（Bari 50K vs 250K）──
    depth_rows = [r for r in summary if '50K稀释' in r['方法']
                  or r['方法'] == 'salmon-k15(小RNA适配)']
    if depth_rows:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        groups = {}
        for r in depth_rows:
            depth = '50K' if '50K' in r['方法'] else '250K'
            groups.setdefault((r['组'], depth), {'tp': 0, 'fn': 0})
            groups[(r['组'], depth)]['tp'] += r['tp']
            groups[(r['组'], depth)]['fn'] += r['fn']
        labs = [f'{k[0]} {k[1]}' for k in groups]
        bars = ax.bar(range(len(groups)),
                      [g['tp'] / (g['tp'] + g['fn']) * 100
                       for _, g in groups.items()], 0.55,
                      color=['#4C72B0', '#7FA6D9', '#5B9A57', '#82BC80'][:len(groups)])
        for i, (k, g) in enumerate(groups.items()):
            ax.text(i, g['tp'] / (g['tp'] + g['fn']) * 100 + 2,
                    f"{g['tp']}/{g['tp'] + g['fn']}", ha='center', fontsize=9,
                    weight='bold')
        ax.axhline(76, ls='--', lw=1.1, color='#D9A441')
        ax.text(len(groups) - 0.4, 77.5, '论文 250K 均值 76%', ha='right',
                fontsize=8, color='#8a6414')
        ax.axhline(46, ls='--', lw=1.1, color='#C0504D')
        ax.text(len(groups) - 0.4, 47.5, '论文 50K 均值 46%', ha='right',
                fontsize=8, color='#8a3535')
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(labs, fontsize=8.5)
        ax.set_ylabel('物种级检出灵敏度 (%)')
        ax.set_ylim(0, 108)
        ax.set_title('测序深度效应：Bari 小RNA 检出灵敏度', fontsize=11)
        ax.spines[['top', 'right']].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, 'fig_b7_depth.png'), dpi=200)
        fig.savefig(os.path.join(FIGS, 'fig_b7_depth.svg'))
        plt.close(fig)

    # fig5 高度随实验数自适应
    print('figures ->', FIGS)


if __name__ == '__main__':
    main()
