#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_speedtest_kv.py — 用平台已建好的 kv_index 实测 salmon_k31 vs minibwa 病毒定量

输入：run/fastq/GQMIX.R1.fq.gz + GQMIX.R2.fq.gz（509,325 对双端片段）
索引：databases/virusref_db/kv_index/{salmon_k31, minibwa}
      （参考 final.cluster.ref.fasta，8,464 条 / 32.43 Mb）

口径（唯一重要的事）：
  salmon 端  = quant（伪比对定量，索引已建好）
  minibwa 端 = map（真比对 SAM）+ samtools view/sort/index
               + pandepth 覆盖度 + idxstats + ANI
  后处理必须算进 minibwa，否则是拿半程对全程。

指标：墙钟耗时（N 轮取中位）+ 子进程内存峰值（psutil 采样 RSS）

用法：
  python _speedtest_kv.py                  # 默认 3 轮
  python _speedtest_kv.py --repeat 5
  python _speedtest_kv.py --threads 16
"""

import argparse
import os
import shutil
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'Virus_Platform_Core'))
import config  # noqa: E402

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'run' / 'fastq'
R1 = DATA / 'GQMIX.R1.fq.gz'
R2 = DATA / 'GQMIX.R2.fq.gz'

KV_IDX = ROOT / 'databases' / 'virusref_db' / 'kv_index'
SALMON_IDX = KV_IDX / 'salmon_k31'
MBW_IDX = KV_IDX / 'minibwa' / 'minibwa'
REF = ROOT / 'databases' / 'virusref_db' / 'final.cluster.ref.fasta'
WORK = ROOT / '_speedtest_work'

TOOLS = config.detect_tools()
SALMON = TOOLS['salmon']
MINIBWA = TOOLS['minibwa']
SAMTOOLS = TOOLS['samtools']
PANDEPTH = TOOLS['pandepth']

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


def dir_size_mb(p):
    p = Path(p)
    if p.is_file():
        return p.stat().st_size / 1024 / 1024
    tot = 0
    for f in p.rglob('*'):
        if f.is_file():
            try:
                tot += f.stat().st_size
            except OSError:
                pass
    return tot / 1024 / 1024


def run(cmd, log_path, env=None):
    """跑命令，返回 (耗时s, 峰值RSS_MB)。psutil 轮询进程树。"""
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    peak = [0.0]
    stop = threading.Event()

    def _sample(p):
        root = None
        while not stop.is_set():
            try:
                if root is None:
                    root = psutil.Process(p.pid)
                tot = 0
                for c in [root] + root.children(recursive=True):
                    try:
                        tot += c.memory_info().rss
                    except Exception:  # noqa: BLE001
                        pass
                if tot > peak[0]:
                    peak[0] = tot
            except Exception:  # noqa: BLE001
                pass
            time.sleep(0.01)

    t0 = time.perf_counter()
    with open(log_path, 'w', encoding='utf-8') as fl:
        p = subprocess.Popen(cmd, stdout=fl, stderr=subprocess.STDOUT, env=env)
        th = None
        if _HAS_PSUTIL:
            th = threading.Thread(target=_sample, args=(p,), daemon=True)
            th.start()
        p.wait()
        stop.set()
        if th:
            th.join(timeout=0.5)
        if p.returncode != 0:
            raise RuntimeError(f'rc={p.returncode}: {" ".join(str(c) for c in cmd[:4])}')
    return time.perf_counter() - t0, peak[0] / 1024 / 1024


# ══════════════════════════════════════════════════════════
def salmon_run(base, logs, threads):
    """salmon quant（索引已建好）"""
    rounds = []
    for phase in ('r1', 'r2', 'r3'):
        qd = base / f'quant_{phase}'
        shutil.rmtree(qd, ignore_errors=True)
        t, m = run([SALMON, 'quant', '-i', str(SALMON_IDX), '-p', str(threads),
                    '-l', 'A', '-1', str(R1), '-2', str(R2), '-o', str(qd)],
                   logs / f'salmon_quant_{phase}.log')
        hits = 0
        qf = qd / 'quant.sf'
        if qf.exists():
            with open(qf, encoding='utf-8') as f:
                f.readline()
                for line in f:
                    c = line.split('\t')
                    if len(c) > 3 and float(c[3]) > 0:
                        hits += 1
        rounds.append({'t': t, 'mem': m, 'hits': hits})
        print(f'    salmon quant  #{phase}   {t:>7.2f} s  {m:>7.0f} MB   '
              f'命中参考 {hits}', flush=True)
    return rounds


def _pd_env():
    env = dict(os.environ)
    env['PATH'] = os.pathsep.join(
        [str(Path(PANDEPTH).parent), str(Path(SAMTOOLS).parent),
         env.get('PATH', '')])
    return env


def minibwa_run(base, logs, threads):
    """minibwa map + 完整后处理"""
    rounds = []
    for i in range(3):
        phase = f'r{i+1}'
        ad = base / f'align_{phase}'
        shutil.rmtree(ad, ignore_errors=True)
        ad.mkdir(parents=True, exist_ok=True)
        sam = ad / 'sample.sam'
        plog = logs / f'minibwa_post_{phase}.log'

        t_map, m_map = run([MINIBWA, 'map', f'-t{threads}', '-u', '-o', str(sam),
                            str(MBW_IDX), str(R1), str(R2)],
                           logs / f'minibwa_map_{phase}.log')
        sam_mb = sam.stat().st_size / 1024 / 1024
        nrows = n_mapped = 0
        with open(sam, encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                nrows += 1
                if not (int(line.split('\t', 2)[1]) & 4):
                    n_mapped += 1

        raw = ad / 'sample.bam'
        sbam = ad / 'sample.sorted.bam'
        post = {}
        for nm, cmd, env in [
            ('view', [SAMTOOLS, 'view', '-b', '-F', '0x04', '-o', str(raw), str(sam)], None),
            ('sort', [SAMTOOLS, 'sort', '-@', str(threads), '-o', str(sbam), str(raw)], None),
            ('index', [SAMTOOLS, 'index', str(sbam)], None),
            ('pandepth', [PANDEPTH, '-a', '-i', str(sbam), '-o', str(ad / 'sample.pd'),
                          '-t', str(threads)], _pd_env()),
            ('idxstats', [SAMTOOLS, 'idxstats', str(sbam)], None),
            ('ANI', [SAMTOOLS, 'view', '-F', '0x904', str(sbam)], None),
        ]:
            tt, _ = run(cmd, plog, env=env)
            post[nm] = tt

        t_post = sum(post.values())
        rounds.append({'map': t_map, 'map_mem': m_map, 'post': t_post,
                       'postdetail': post, 'sam_mb': sam_mb,
                       'nrows': nrows, 'n_mapped': n_mapped})
        print(f'    minibwa map  #{phase}   {t_map:>7.2f} s  {m_map:>7.0f} MB   '
              f'SAM {sam_mb:.1f} MB / {nrows} 行', flush=True)
        print(f'      └ 后处理 {t_post:.2f} s '
              f'({" ".join(f"{k}:{v:.2f}" for k, v in post.items())})', flush=True)
    return rounds


# ══════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--repeat', type=int, default=3)
    a = ap.parse_args()

    print('=' * 72)
    print('salmon_k31 vs minibwa — 病毒定量实测（平台现成索引）')
    print('=' * 72)
    print(f'  输入    {R1.name} + {R2.name}  (509,325 对双端片段)')
    print(f'  参考    {REF.name}  (8,464 条 / 32.43 Mb)')
    print(f'  索引    salmon_k31 {dir_size_mb(SALMON_IDX):.0f} MB   '
          f'minibwa {dir_size_mb(MBW_IDX):.0f} MB')
    for nm, p in [('salmon', SALMON), ('minibwa', MINIBWA),
                  ('samtools', SAMTOOLS), ('pandepth', PANDEPTH)]:
        if not (p and Path(p).exists()):
            print(f'  工具缺失: {nm} -> {p}')
            return 1
    for p in (SALMON_IDX, MBW_IDX.with_suffix('.mbw'), REF):
        if not p.exists():
            print(f'  路径缺失: {p}')
            return 1
    print(f'  线程 {a.threads}   轮数 {a.repeat}')

    base = WORK
    shutil.rmtree(base, ignore_errors=True)
    logs = base / 'logs'
    logs.mkdir(parents=True, exist_ok=True)

    print('\n[1] salmon_k31（伪比对定量）', flush=True)
    sr = salmon_run(base, logs, a.threads)

    print('\n[2] minibwa（真比对 + 覆盖度 + ANI）', flush=True)
    mr = minibwa_run(base, logs, a.threads)

    s_t = statistics.median([r['t'] for r in sr])
    s_m = statistics.median([r['mem'] for r in sr])
    m_t = statistics.median([r['map'] for r in mr])
    m_p = statistics.median([r['post'] for r in mr])
    m_m = statistics.median([r['map_mem'] for r in mr])

    print(f'\n{"=" * 72}')
    print('结果（中位数）')
    print(f'{"=" * 72}')
    print(f'  {"阶段":<34}{"salmon":>13}{"minibwa":>13}')
    print(f'  {"-" * 60}')
    print(f'  {"定量 / 比对 map":<32}{s_t:>11.2f}s{m_t:>11.2f}s')
    print(f'  {"后处理（覆盖度/ANI 等）":<30}{"—":>13}{m_p:>11.2f}s')
    print(f'  {"-" * 60}')
    total_m = m_t + m_p
    print(f'  {"合计":<34}{s_t:>11.2f}s{total_m:>11.2f}s')
    print(f'  {"内存峰值":<32}{s_m:>10.0f}MB{m_m:>10.0f}MB')
    ratio = total_m / s_t if s_t else 0
    print(f'\n  minibwa / salmon = {ratio:.2f}x')
    print(f'  侧证：一次运行的实际数据量 {statistics.median([r["nrows"] for r in mr])} 行 SAM')

    pd_ = mr[0]['postdetail']
    print('\n  minibwa 后处理拆解:')
    for k, v in pd_.items():
        print(f'    {k:<16}{v:>8.2f} s')

    print(f'\n  产物: {WORK}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
