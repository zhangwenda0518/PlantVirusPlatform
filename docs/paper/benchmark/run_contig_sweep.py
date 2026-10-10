# -*- coding: utf-8 -*-
"""min contig 门槛扫描：同一套 S3 srna 组装 contigs，按不同最短长度阈值
过滤后对 PVB RNA1/RNA2 真值做覆盖/一致率评估，量化「放宽 min contig 的收益」。

用法: python run_contig_sweep.py <out_dir>
结束打印 RESULT:{...}
"""
import glob
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, REPO)

# 评测快照参考库（2026-09 的 kv_index；2026-09-29 起被 plant_ref 取代）。
# 优先级：环境变量 VP_BENCH_DB > 旧默认路径；缺库时报清晰中文错误。
_DB_DIR = os.environ.get('VP_BENCH_DB') or os.path.join(
    REPO, 'databases', 'virusref_db', 'kv_index')
KV_REF = os.path.join(_DB_DIR, 'reference.fasta')
if not os.path.isfile(KV_REF):
    raise SystemExit(
        f'[基准库缺失] 找不到参考库 FASTA：{KV_REF}\n'
        f'说明：本评测快照绑定 2026-09 的 kv_index 库（manifest 2026-09-19），'
        f'2026-09-29 起平台默认库已改为 virusref_db/plant_ref/，旧目录已移除。\n'
        f'复现请先重建等价库，或用环境变量 VP_BENCH_DB=<库目录> 显式指定。')
# Bari 原始数据路径：评测机外部目录，可用环境变量 VP_BENCH_BARI 覆盖
BARI_S3 = os.environ.get(
    'VP_BENCH_BARI',
    r'E:\谷歌下载\测试数据验证\Raw-data-Bari.Workshop\Sample_3.fastq')
CUTOFFS = [30, 40, 50, 60, 80, 100, 150, 200]
TARGETS = ['NC_043447.1', 'NC_043448.1']
OUT = (sys.argv[1] if len(sys.argv) > 1
       else os.environ.get('VP_BENCH_SWEEP_OUT')
       or r'C:\vp_srna_test\bench\contig_sweep')


def read_fa(path):
    seqs, name, cur = {}, None, []
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('>'):
            if name:
                seqs[name] = ''.join(cur)
            name = ln[1:].split()[0]
            cur = []
        else:
            cur.append(ln.strip())
    if name:
        seqs[name] = ''.join(cur)
    return seqs


def merge(iv):
    iv.sort()
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def blast_cov(contig_fa, subject_fa, blastn):
    """返回 (覆盖 bases, HSP 总长, 一致 bases) —— 对 subject 合并区间。"""
    r = subprocess.run(
        [blastn, '-query', contig_fa, '-subject', subject_fa,
         '-outfmt', '6 sstart send length nident',
         '-evalue', '1e-3', '-num_threads', '4'],
        capture_output=True, text=True, encoding='utf-8', errors='replace')
    iv, aln, ident = [], 0, 0
    for ln in r.stdout.splitlines():
        p = ln.split('\t')
        ss, se, ln_, ni = map(int, p)
        iv.append((min(ss, se), max(ss, se)))
        aln += ln_
        ident += ni
    cbases = sum(b - a for a, b in merge(iv))
    return cbases, aln, ident


def main():
    os.makedirs(OUT, exist_ok=True)
    peak = {'rss': 0}

    def _mon():
        import psutil
        p = psutil.Process()
        while not stop.is_set():
            try:
                t = p.memory_info().rss
                for c in p.children(recursive=True):
                    try:
                        t += c.memory_info().rss
                    except (psutil.Error, OSError):
                        pass
                peak['rss'] = max(peak['rss'], t)
            except (psutil.Error, OSError):
                pass
            time.sleep(1.0)

    stop = threading.Event()
    threading.Thread(target=_mon, daemon=True).start()

    from Virus_Platform_Core.assembly import run_spades
    scratch = os.path.join(REPO, '_zz_bench_sweep')
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        sp = run_spades(BARI_S3, None, scratch, mode='srna', threads=4,
                        memory_gb=8)
        asm_s = round(time.time() - t0, 1)
        contigs = []
        nm, cur = None, []
        for ln in open(sp, encoding='utf-8'):
            if ln.startswith('>'):
                if nm:
                    contigs.append((nm, ''.join(cur)))
                nm, cur = ln[1:].split()[0], []
            else:
                cur.append(ln.strip())
        if nm:
            contigs.append((nm, ''.join(cur)))

        truth = read_fa(KV_REF)
        blastn = glob.glob(os.path.join(REPO, '3rd', 'tools', 'Blast', 'bin',
                                        'blastn.exe'))[0]

        sweep = []
        for cut in CUTOFFS:
            sub = [(n, s) for n, s in contigs if len(s) >= cut]
            qfa = os.path.join(OUT, f'contigs_ge{cut}.fasta')
            with open(qfa, 'w', encoding='utf-8') as f:
                for n2, s2 in sub:
                    f.write(f'>{n2}\n')
                    for i in range(0, len(s2), 70):
                        f.write(s2[i:i + 70] + '\n')
            row = {'cutoff': cut, 'n_contigs': len(sub),
                   'total_bp': sum(len(s) for _, s in sub), 'targets': {}}
            for acc in TARGETS:
                truth_seq = truth[acc]
                sfa = os.path.join(OUT, f'truth_{acc}.fasta')
                if not os.path.isfile(sfa):
                    with open(sfa, 'w', encoding='utf-8') as f:
                        f.write(f'>{acc}\n')
                        for i in range(0, len(truth_seq), 70):
                            f.write(truth_seq[i:i + 70] + '\n')
                cb, aln, ident = blast_cov(qfa, sfa, blastn)
                row['targets'][acc] = {
                    'coverage_pct': round(cb / len(truth_seq) * 100, 2),
                    'identity_pct': round(ident / aln * 100, 2) if aln else None}
            sweep.append(row)
            tstr = ' / '.join(f"{row['targets'][a]['coverage_pct']:.1f}%"
                              for a in TARGETS)
            print(f"[cut≥{cut:>3}] n={len(sub):>3} bp={row['total_bp']:>6} "
                  f"PVB覆盖(RNA1/RNA2): {tstr}", flush=True)
            os.remove(qfa)

        res = {'asm_wall_s': asm_s, 'n_contigs_total': len(contigs),
               'sweep': sweep, 'peak_rss_mb': round(peak['rss'] / 1048576, 1)}
        print('RESULT:' + json.dumps(res, ensure_ascii=False))
        sys.stdout.flush()
        os._exit(0)
    finally:
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == '__main__':
    main()
