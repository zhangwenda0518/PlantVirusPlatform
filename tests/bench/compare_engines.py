# -*- coding: utf-8 -*-
"""对照 VIROMOCK 官方基准（基准测试.xlsx）评估双引擎定量/检出。

归一化口径：观测占比 = 该株定量值 / 「与任一目标株同物种的全部参考」定量值之和
（消除同物种额外参考分摊 reads 的影响，如 D14/D5 中的 PVY 参考池）。
"""
import csv
import json
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
# 平台根：默认自动推导（tests/bench/ 上溯 3 层 = 平台根），
# 可用环境变量 VP_ROOT 覆盖（换机/换目录不用改代码）。
PLATFORM = (os.environ.get('VP_ROOT')
            or os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))))
OUT_ROOT = os.path.join(PLATFORM, '_bench_kvs2')
LIB_INFO = os.path.join(PLATFORM, 'databases', 'virusref_db', 'viromock_kv',
                        'reference.ref_info.tsv')

# 官方 ground truth（基准测试.xlsx：人工添加株预期比例 %）
TRUTH = {
    'D1':  [('JQ911663', 69.75), ('KU883267', 29.22), ('MH323442', 1.03)],
    'D5':  [('EF026076', 80), ('AY884983', 20)],
    'D6':  [('PVY_artificial_strain', 100)],
    'D10': [('PPV_MK387313_artificial', 100)],
    'D11': [('DQ000985', 39.9), ('AJ606359', 15.32), ('MF422616', 14.96),
            ('JQ314460', 9.94), ('MK133092', 9.94), ('HG313807', 9.94)],
    'D12': [('HE979770', 50), ('HE979758', 19.98), ('AJ314739', 15.01), ('KR611579', 15.01)],
    'D13': [('DQ451009', 37.93), ('KT895259', 15.14), ('DQ092436', 15.76),
            ('AY750155', 10.73), ('KT895258', 10.53), ('AY493509', 9.91)],
    'D14': [('AB711147', 37.79), ('KC634004', 19.63), ('MF176828', 18.89),
            ('JQ969039', 19.63), ('FJ214726', 9.56)],
    'D15': [('LN680656', 49.95), ('FR751552', 25.02), ('KJ082087', 25.02)],
    'D16': [('KX977568', 50), ('KR080326', 20), ('JN019858', 15), ('JQ951943', 15)],
    'D17': [('MH300061', 40.28), ('KX192366', 19.91), ('EU715989', 14.93),
            ('LN794218', 14.93), ('MG934545', 9.95)],
    'D18': [('EF521843', 31.11), ('KF523382', 15.58), ('KY593456', 15.58),
            ('D11028', 14.73), ('KC559092', 12.63), ('EU332308', 10.37)],
}


def load_species_map():
    m = {}
    with open(LIB_INFO, encoding='utf-8') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            acc = (row['Accession'] or '').split('.')[0]
            m[acc] = row.get('Species_NCBI') or ''
    return m


ACC2SP = load_species_map()


def quant_value(row, engine):
    for col in (('EM_Reads', 'Uniq_Reads') if engine == 'salmon'
                else ('Uniq_Reads', 'EM_Reads')):
        v = (row.get(col) or '').strip()
        if v:
            try:
                return float(v)
            except ValueError:
                pass
    return 0.0


def eval_run(ds, engine):
    pj = os.path.join(OUT_ROOT, f'{ds}_{engine}', 'quant_rows.json')
    if not os.path.isfile(pj):
        return None
    rows = json.load(open(pj, encoding='utf-8'))
    targets = TRUTH[ds]

    def row_of(acc):
        for r in rows:
            if (r.get('Accession') or '').split('.')[0] == acc:
                return r
        return None

    truth_species = {ACC2SP.get(a, '') for a, _ in targets}

    def in_pool(r):
        sp = ACC2SP.get((r.get('Accession') or '').split('.')[0], '')
        return any(sp == s or (s and sp.startswith(s)) for s in truth_species)

    pool = [r for r in rows if in_pool(r)]
    tot_pool = sum(quant_value(r, engine) for r in pool)
    detail = []
    n_det = 0
    errs = []
    for acc, exp in targets:
        r = row_of(acc)
        v = quant_value(r, engine) if r else 0.0
        obs = (v / tot_pool * 100) if tot_pool else 0.0
        det = v > 0
        n_det += det
        d = {'acc': acc, 'exp': exp, 'obs': round(obs, 2), 'detected': det,
             'cov': (r or {}).get('Coverage(%)', '')}
        if det:
            d['abs_err_pp'] = round(abs(obs - exp), 2)
            errs.append(d['abs_err_pp'])
        detail.append(d)
    return {'ds': ds, 'engine': engine, 'targets': len(targets),
            'detected': n_det,
            'mae_pp': round(sum(errs) / len(errs), 2) if errs else None,
            'max_err_pp': round(max(errs), 2) if errs else None,
            'detail': detail}


def main():
    print('===== 引擎 × 数据集 对照官方基准 =====')
    summary = []
    for ds in TRUTH:
        for eng in ('salmon', 'minibwa'):
            r = eval_run(ds, eng)
            if r:
                summary.append(r)
                print(f"{ds:>4}/{eng:<8} 检出 {r['detected']}/{r['targets']}"
                      f"  MAE={r['mae_pp']}pp  最大={r['max_err_pp']}pp")
            else:
                print(f"{ds:>4}/{eng:<8} —— 无结果")
    print('\n===== 引擎级汇总 =====')
    for eng in ('salmon', 'minibwa'):
        rs = [r for r in summary if r['engine'] == eng]
        if not rs:
            continue
        n_t = sum(r['targets'] for r in rs)
        n_d = sum(r['detected'] for r in rs)
        errs = [d['abs_err_pp'] for r in rs for d in r['detail'] if d['detected']]
        print(f"{eng}: 检出 {n_d}/{n_t} ({n_d/n_t*100:.0f}%)"
              f"  定量MAE={sum(errs)/len(errs):.2f}pp  最大={max(errs):.1f}pp"
              f"  (n={len(errs)})")
    out_p = os.path.join(OUT_ROOT, 'compare_summary.json')
    json.dump(summary, open(out_p, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    print(f'\n明细 -> {out_p}')


if __name__ == '__main__':
    main()
