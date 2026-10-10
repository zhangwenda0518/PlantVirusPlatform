# -*- coding: utf-8 -*-
"""共识序列比较评测：
  A. 模拟自洽（注意：**不是**准确度测验）：库内参考 → 模拟 150bp PE **无错误**
     reads（200× 平铺、无测序错误、无污染）→ 平台共识(minibwa) → 与源参考
     逐碱基比对 → 覆盖率/一致率。所得"零引入错误"是无错输入下的自洽性测试，
     仅证明管线自身不引入变异，不表征真实数据（含测序错误/低质量碱基）的
     处理能力。（2026-10-03 披露）
  B. Bari 小RNA 实测：S3 真实 sRNA reads → 平台共识（0 映射 → 全 N 静默产物，
     发现=记录）
  C. 组装替代路线：S3 srna 组装 contigs 对 PVB 两段的覆盖/一致率

参数披露（2026-10-03）：A/B/C 均显式用 vc_depth=5、vc_qual=20、vc_freq=0.5。
评测当时（2026-09）平台默认 vc_depth 实为 1（帮助文本误写 5，2026-10 已修正
默认=5），故本结果不代表平台当时默认口径；**平台默认 vc_depth=1（旧默认）下的
表现未测**。

用法: python run_cons_bench.py <out_dir>
结束打印 RESULT:{...}
"""
import argparse
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
BARI_S3 = (os.environ.get('VP_BENCH_BARI')
           or r'E:\谷歌下载\测试数据验证\Raw-data-Bari.Workshop\Sample_3.fastq')
OUT = (sys.argv[1] if len(sys.argv) > 1
       else os.environ.get('VP_BENCH_CONS_OUT')
       or r'C:\vp_srna_test\bench\cons')

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(COMP)[::-1]


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


def sim_pe_reads(seq, acc, out_dir, read_len=150, depth=200):
    """无错误 PE reads：步进平铺 + 每层错位 37bp，全覆盖 depth×。"""
    fq1 = os.path.join(out_dir, f'{acc}_R1.fastq')
    fq2 = os.path.join(out_dir, f'{acc}_R2.fastq')
    L = len(seq)
    n = 0
    with open(fq1, 'w') as f1, open(fq2, 'w') as f2:
        for rep in range(depth):
            off = (rep * 37) % (read_len * 2)
            pos = -off
            while pos + read_len * 2 <= L:
                a = seq[max(pos, 0):pos + read_len]
                b = seq[pos + read_len:pos + read_len * 2]
                if len(a) == read_len and len(b) == read_len and pos >= 0:
                    n += 1
                    f1.write(f'@{acc}_{n}/1\n{a}\n+\n{"I" * read_len}\n')
                    f2.write(f'@{acc}_{n}/2\n{revcomp(b)}\n+\n{"I" * read_len}\n')
                pos += read_len * 2
    return n


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
            time.sleep(0.5)

    stop = threading.Event()
    threading.Thread(target=_mon, daemon=True).start()

    from Virus_Platform_Core.known_virus_suite.kv_common import ToolRegistry
    from Virus_Platform_Core.known_virus_suite.kv_consensus import ConsensusStage
    import polars as pl

    tools = ToolRegistry()
    tools.probe()
    tools.require('minibwa', 'samtools', 'viral_consensus')

    def make_args():
        return argparse.Namespace(
            reference=KV_REF, threads=4, resume=False, jobs=2,
            vc_depth=5, vc_qual=20, vc_freq=0.5, vc_ambig='N')

    seqs = read_fa(KV_REF)
    res = {'targets': []}

    # ── A. 模拟自洽（3 个参考：PVB RNA1/RNA2 + GLRaV-1）──
    targets = ['NC_043447.1', 'NC_043448.1', 'NC_016509.1']
    sim_dir = os.path.join(OUT, 'sim')
    os.makedirs(sim_dir, exist_ok=True)
    t_all = time.time()
    for acc in targets:
        t0 = time.time()
        truth = seqs[acc]
        n = sim_pe_reads(truth, acc, sim_dir, depth=200)
        df = pl.DataFrame({'Sample': [f'sim_{acc.split(".")[0]}'],
                           'Accession': [acc]})
        # 每目标独立 stage 目录：harmonized 索引按目标隔离（共用会缓存
        # 第一个目标的参考，后续目标比对到错误索引上）
        stage_dir = os.path.join(sim_dir, 'stage_' + acc.split('.')[0])
        stage = ConsensusStage(make_args(), tools, _Log(), stage_dir)
        stage.run([{'name': f'sim_{acc.split(".")[0]}', 'r1': f'{sim_dir}/{acc}_R1.fastq',
                    'r2': f'{sim_dir}/{acc}_R2.fastq'}], df)
        wall = time.time() - t0
        # 产物在 consensus/Unannotated_<acc>/ 下（无扩展名），rglob 收集
        cons_files = [str(p) for p in Path(stage_dir, 'consensus').rglob('*')
                      if p.is_file()]
        ev = {'acc': acc, 'ref_len': len(truth), 'sim_pairs': n,
              'wall_s': round(wall, 1), 'n_cons_files': len(cons_files)}
        best = _best_eval(cons_files, acc, {acc: truth})
        ev.update(best)
        res['targets'].append(ev)
        print(f"[A] {acc}: {ev}", flush=True)
    res['total_wall_s'] = round(time.time() - t_all, 1)

    # ── B. Bari sRNA 共识实测（PVB RNA1/RNA2 两目标）──
    t0 = time.time()
    df = pl.DataFrame({'Sample': ['S3', 'S3'],
                       'Accession': ['NC_043447.1', 'NC_043448.1']})
    stage_dir_b = os.path.join(OUT, 'srna')
    stage = ConsensusStage(make_args(), tools, _Log(), stage_dir_b)
    try:
        stage.run([{'name': 'S3', 'r1': BARI_S3, 'r2': None}], df)
        srna = {'ran': True}
    except Exception as e:
        srna = {'ran': False, 'error': str(e)[:200]}
    srna['wall_s'] = round(time.time() - t0, 1)
    cons = [str(p) for p in Path(stage_dir_b, 'consensus').rglob('*')
            if p.is_file()]
    truth = {a: seqs[a] for a in ('NC_043447.1', 'NC_043448.1')}
    srna['eval'] = _best_eval(cons, None, truth) if cons else {
        'coverage_pct': 0.0, 'identity_pct': None, 'note': '无共识产物'}
    res['srna'] = srna
    print(f"[B] sRNA: {json.dumps(srna, ensure_ascii=False)[:300]}", flush=True)

    # ── C. 组装替代路线（S3 srna 组装 contigs 对 PVB 两段）──
    t0 = time.time()
    from Virus_Platform_Core.assembly import run_spades
    scratch = os.path.join(REPO, '_zz_bench_cons_asm')
    os.makedirs(scratch, exist_ok=True)
    try:
        sp = run_spades(BARI_S3, None, scratch, mode='srna', threads=4,
                        memory_gb=8)
        contigs = []
        nm, cur = None, []
        for ln in open(sp, encoding='utf-8'):
            if ln.startswith('>'):
                if nm:
                    contigs.append((nm, ''.join(cur)))
                nm, cur = ln[1:].strip(), []
            else:
                cur.append(ln.strip())
        if nm:
            contigs.append((nm, ''.join(cur)))
        wall = time.time() - t0
        tf = os.path.join(OUT, 's3_contigs.fasta')
        with open(tf, 'w', encoding='utf-8') as f:
            for nm2, s2 in contigs:
                f.write(f'>{nm2}\n')
                for i in range(0, len(s2), 70):
                    f.write(s2[i:i + 70] + '\n')
        blastn = glob.glob(os.path.join(REPO, '3rd', 'tools', 'Blast', 'bin',
                                        'blastn.exe'))[0]
        route = {}
        for acc, truth in truth.items():
            tfa = os.path.join(OUT, f'truth_{acc}.fasta')
            with open(tfa, 'w', encoding='utf-8') as f:
                f.write(f'>{acc}\n')
                for i in range(0, len(truth), 70):
                    f.write(truth[i:i + 70] + '\n')
            r = subprocess.run(
                [blastn, '-query', tf, '-subject', tfa,
                 '-outfmt', '6 qstart qend sstart send length nident gaps',
                 '-evalue', '1e-5', '-num_threads', '4'],
                capture_output=True, text=True, encoding='utf-8',
                errors='replace')
            covered, ident, aln = 0, 0, 0
            iv = []
            for ln in r.stdout.splitlines():
                p = ln.split('\t')
                qs, qe, ss, se, ln_, ni, gp = map(int, p)
                iv.append((min(ss, se), max(ss, se)))
                covered += ln_
                ident += ni
                aln += ln_
            merged = _merge(iv)
            cbases = sum(b - a for a, b in merged)
            route[acc] = {
                'coverage_pct': round(cbases / len(truth) * 100, 2),
                'identity_pct': round(ident / aln * 100, 2) if aln else None,
                'aln_bases': aln}
        res['assembly_route'] = {'wall_s': round(wall, 1), 'n_contigs': len(contigs),
                                 'targets': route}
        print(f"[C] 组装路线: {json.dumps(route)[:300]}", flush=True)
    finally:
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)

    stop.set()
    time.sleep(0.6)
    res['peak_rss_mb'] = round(peak['rss'] / 1048576, 1)
    print('RESULT:' + json.dumps(res, ensure_ascii=False))
    sys.stdout.flush()
    os._exit(0)


def _merge(iv):
    iv.sort()
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _best_eval(cons_files, acc, truth):
    """共识产物 vs 真值参考：返回覆盖率/一致率最优的一组。

    acc=None 时（sRNA 多目标）逐真值尝试，取覆盖最高者。
    """
    best = None
    for cf in cons_files:
        seq = read_fa(cf)
        if not seq:
            continue
        cname, cseq = max(seq.items(), key=lambda kv: len(kv[1]))
        if acc is None:      # sRNA：两目标都试，取覆盖最高者
            candidates = sorted(truth.items())
        else:
            key = acc.split('.')[0]
            candidates = [(k, v) for k, v in truth.items()
                          if k.startswith(key)]
        if not candidates:
            continue
        cand_best = None
        for ck, target in candidates:
            ev = _eval_pair(cname, cseq, ck, target)
            if ev and (cand_best is None or
                       ev['coverage_pct'] > cand_best['coverage_pct']):
                cand_best = ev
        if cand_best is None:
            continue
        if best is None or cand_best['coverage_pct'] > best['coverage_pct']:
            best = cand_best
    return best or {'coverage_pct': 0.0, 'identity_pct': None,
                    'note': '共识与真值无比对'}


def _eval_pair(cname, cseq, acc, target):
    import subprocess as _sp
    blastn = glob.glob(os.path.join(REPO, '3rd', 'tools', 'Blast', 'bin',
                                    'blastn.exe'))[0]
    tfa = os.path.join(OUT, '_t.fasta')
    with open(tfa, 'w', encoding='utf-8') as f:
        f.write('>' + cname + '\n')
        for i in range(0, len(cseq), 70):
            f.write(cseq[i:i + 70] + '\n')
    t2 = os.path.join(OUT, '_t2.fasta')
    with open(t2, 'w', encoding='utf-8') as f:
        f.write('>' + acc + '\n')
        for i in range(0, len(target), 70):
            f.write(target[i:i + 70] + '\n')
    r = _sp.run(
        [blastn, '-query', tfa, '-subject', t2,
         '-outfmt', '6 qstart qend sstart send length nident gaps',
         '-evalue', '1e-10', '-num_threads', '4'],
        capture_output=True, text=True, encoding='utf-8', errors='replace')
    cov, ident, aln = 0, 0, 0
    iv = []
    for ln in r.stdout.splitlines():
        p = ln.split('\t')
        qs, qe, ss, se, ln_, ni, gp = map(int, p)
        iv.append((min(ss, se), max(ss, se)))
        cov += ln_
        ident += ni
        aln += ln_
    merged = _merge(iv)
    cbases = sum(b - a for a, b in merged)
    return {'cons_len': len(cseq),
            'coverage_pct': round(cbases / len(target) * 100, 2),
            'identity_pct': round(ident / aln * 100, 2) if aln else None,
            'mismatch_gaps': aln - ident, 'target': acc}


class _Log:
    def _p(self, *a, **k):
        print(*a, **k, flush=True)

    def __getattr__(self, attr):
        # info/warning/error/debug/exception 等任意日志级别都兜住
        return self._p


if __name__ == '__main__':
    main()
