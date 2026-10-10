# -*- coding: utf-8 -*-
"""共识序列完整性/准确性基准：真值参考子集 + 全量 reads 回贴 → 共识 vs 真值。"""
import csv
import json
import os
import sys
import time

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
# 平台根：默认自动推导（tests/bench/ 上溯 3 层 = 平台根），
# 可用环境变量 VP_ROOT 覆盖（换机/换目录不用改代码）。
PLATFORM = (os.environ.get('VP_ROOT')
            or os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, PLATFORM)
# 基准数据集（作者本机路径）：可用环境变量 VP_BENCH_DATA 覆盖
DATA = os.environ.get('VP_BENCH_DATA', r'E:\谷歌下载\测试数据验证')
OUT_ROOT = os.path.join(PLATFORM, '_bench_cons')
LIB_FA = os.path.join(PLATFORM, 'databases', 'virusref_db', 'viromock_kv',
                      'reference.fasta')


def _safe_path(base, *parts):
    """规范化并确保落在 base 目录内（防拼接越界）。"""
    base = os.path.normpath(base)
    p = os.path.normpath(os.path.join(base, *parts))
    if os.path.commonpath([base, p]) != base:
        raise ValueError(f'路径越界: {p}')
    return p


D = lambda n: os.path.join(DATA, n)
# 数据集 -> (reads, 真值参考 accession 列表)
TARGETS = {
    'D1':  (D('Dataset_1_R1.fastq.gz'), D('Dataset_1_R2.fastq.gz'),
            ['JQ911663', 'KU883267', 'MH323442']),
    'D4':  (D('Dataset_4_R1.fastq.gz'), D('Dataset_4_R2.fastq.gz'),
            ['DQ377131']),
    'D5':  (D('Dataset_5.fastq.gz'), None,
            ['EF026076', 'AY884983']),
    'D6':  (D('Dataset_6.fastq.gz'), None,
            ['PVY_artificial_strain']),
    'D7':  (D('Dataset_7_R1.fastq.gz'), D('Dataset_7_R2.fastq.gz'),
            ['NC_002050', 'NC_002051', 'NC_002052']),
    'D8':  (D('Dataset_8_R1.fastq.gz'), D('Dataset_8_R2.fastq.gz'),
            ['NC_005286', 'NC_040543']),
    'D9':  (D('Dataset_9_R1.fastq.gz'), D('Dataset_9_R2.fastq.gz'),
            ['NC_078016', 'NC_078017', 'NC_078018', 'NC_078019',
             'NC_078020', 'NC_078021', 'NC_078022', 'NC_078023']),
    'D10': (D('Dataset_10.fastq.gz'), None,
            ['PPV_MK387313_artificial']),
    'D11': (D('Dataset_11_R1.fastq.gz'), D('Dataset_11_R2.fastq.gz'),
            ['DQ000985', 'AJ606359', 'MF422616', 'JQ314460', 'MK133092', 'HG313807']),
    'D12': (D('Dataset_12_R1.fastq.gz'), D('Dataset_12_R2.fastq.gz'),
            ['HE979770', 'HE979758', 'AJ314739', 'KR611579']),
    'D13': (D('Dataset_13_R1.fastq.gz'), D('Dataset_13_R2.fastq.gz'),
            ['DQ451009', 'KT895259', 'DQ092436', 'AY750155', 'KT895258', 'AY493509']),
    'D14': (D('Dataset_14_R1.fastq.gz'), D('Dataset_14_R2.fastq.gz'),
            ['AB711147', 'KC634004', 'MF176828', 'JQ969039', 'FJ214726']),
    'D15': (D('Dataset_15_R1.fastq.gz'), D('Dataset_15_R2.fastq.gz'),
            ['LN680656', 'FR751552', 'KJ082087']),
    'D16': (D('Dataset_16_R1.fastq.gz'), D('Dataset_16_R2.fastq.gz'),
            ['KX977568', 'KR080326', 'JN019858', 'JQ951943']),
    'D17': (D('Dataset_17_R1.fastq.gz'), D('Dataset_17_R2.fastq.gz'),
            ['MH300061', 'KX192366', 'EU715989', 'LN794218', 'MG934545']),
    'D18': (D('Dataset_18_R1.fastq.gz'), D('Dataset_18_R2.fastq.gz'),
            ['EF521843', 'KF523382', 'KY593456', 'D11028', 'KC559092', 'EU332308']),
}


def load_lib():
    recs, name, seq = {}, None, []
    for line in open(LIB_FA, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            if name:
                recs[name] = ''.join(seq)
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        recs[name] = ''.join(seq)
    return recs


def subset_fasta(lib, accs, path):
    found = {}
    with open(path, 'w', encoding='utf-8') as f:
        for acc in accs:
            for k, s in lib.items():
                if k.split('.')[0] == acc:
                    f.write(f'>{k}\n{s}\n')
                    found[acc] = k
                    break
    return found


def main():
    only = sys.argv[1:] or None
    from Virus_Platform_Core.consensus import consensus_and_variants
    lib = load_lib()
    os.makedirs(OUT_ROOT, exist_ok=True)
    results = []
    for ds, (r1, r2, accs) in TARGETS.items():
        if only and ds not in only:
            continue
        run_dir = _safe_path(OUT_ROOT, ds)
        os.makedirs(run_dir, exist_ok=True)
        fa = os.path.join(run_dir, 'targets.fasta')
        found = subset_fasta(lib, accs, fa)
        missing = [a for a in accs if a not in found]
        t0 = time.time()
        err = None
        try:
            consensus_and_variants(run_dir, fa, [(r1, r2)], threads=19,
                                   logger=type('L', (), {'log': lambda self, *a: None})())
        except Exception as e:
            err = repr(e)
        wall = round(time.time() - t0, 1)

        # 解析 coverage.tsv + consensus.fa + 真值比对
        cov = {}
        cov_p = os.path.join(run_dir, 'consensus', 'coverage.tsv')
        if os.path.isfile(cov_p):
            with open(cov_p, encoding='utf-8') as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    cov[(row.get('contig') or '').split('.')[0]] = row
        cons, cname, cseq = {}, None, []
        for line in open(os.path.join(run_dir, 'consensus', 'consensus.fa'),
                         encoding='utf-8'):
            if line.startswith('>'):
                if cname:
                    cons[cname] = ''.join(cseq)
                cname, cseq = line[1:].split()[0], []
            else:
                cseq.append(line.strip())
        if cname:
            cons[cname] = ''.join(cseq)

        detail = []
        for acc, full in found.items():
            ref = lib.get(full) or lib.get(acc)
            c = cons.get(full) or cons.get(acc) or ''
            row = cov.get(acc, {})
            n_pos = sum(1 for a, b in zip(ref, c) if b != 'N')
            match = sum(1 for a, b in zip(ref, c) if b == a and b != 'N')
            detail.append({
                'acc': acc, 'ref_len': len(ref),
                'coverage_pct': row.get('coverage_pct'),
                'called_pct': row.get('called_pct'),
                'identity_pct': round(match * 100.0 / n_pos, 3) if n_pos else None,
                'n_snv': row.get('n_snv'), 'n_isnv': row.get('n_isnv'),
                'mean_depth': row.get('mean_depth'),
            })
        res = {'ds': ds, 'wall_sec': wall, 'error': err, 'detail': detail,
               'missing_refs': missing}
        json.dump(res, open(os.path.join(run_dir, 'bench_cons.json'), 'w',
                            encoding='utf-8'), ensure_ascii=False, indent=1)
        results.append(res)
        print(f">>> {ds}: {wall}s err={err}", flush=True)
        for d in detail:
            print(f"    {d['acc']:<24} 覆盖={d['coverage_pct']}% "
                  f"可信位点={d['called_pct']}% 一致率={d['identity_pct']}% "
                  f"SNV={d['n_snv']} iSNV={d['n_isnv']} depth={d['mean_depth']}",
                  flush=True)
    json.dump(results, open(os.path.join(OUT_ROOT, 'bench_cons_summary.json'),
                            'w', encoding='utf-8'), ensure_ascii=False, indent=1)


if __name__ == '__main__':
    main()
