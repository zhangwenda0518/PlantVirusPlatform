# -*- coding: utf-8 -*-
"""双引擎（salmon/minibwa）× 18 数据集跑批：直接驱动引擎入口 identify+filter。

不抽病毒 reads（本基准只对比定量），输出目录 _bench_kvs2/<数据集>_<引擎>/。
所有进入命令行/路径的变量均来自本文件常量表（DATASETS/ENGINES），不接受外部输入。
"""
import csv
import json
import os
import subprocess
import sys
import threading
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
OUT_ROOT = os.path.join(PLATFORM, '_bench_kvs2')
LIB = os.path.join(PLATFORM, 'databases', 'virusref_db', 'viromock_kv')
REF = os.path.join(LIB, 'reference.fasta')
INFO = os.path.join(LIB, 'reference.ref_info.tsv')


def _safe_path(base, *parts):
    """规范化并确保落在 base 目录内（防拼接越界）。"""
    base = os.path.normpath(base)
    p = os.path.normpath(os.path.join(base, *parts))
    if os.path.commonpath([base, p]) != base:
        raise ValueError(f'路径越界: {p}')
    return p


D = lambda n: os.path.join(DATA, n)
DATASETS = [
    ('D1',  D('Dataset_1_R1.fastq.gz'),  D('Dataset_1_R2.fastq.gz')),
    ('D2',  D('Dataset_2_R1.fastq.gz'),  D('Dataset_2_R2.fastq.gz')),
    ('D3',  D('Dataset_3_R1.fastq.gz'),  D('Dataset_3_R2.fastq.gz')),
    ('D4',  D('Dataset_4_R1.fastq.gz'),  D('Dataset_4_R2.fastq.gz')),
    ('D5',  D('Dataset_5.fastq.gz'),     None),
    ('D6',  D('Dataset_6.fastq.gz'),     None),
    ('D7',  D('Dataset_7_R1.fastq.gz'),  D('Dataset_7_R2.fastq.gz')),
    ('D8',  D('Dataset_8_R1.fastq.gz'),  D('Dataset_8_R2.fastq.gz')),
    ('D9',  D('Dataset_9_R1.fastq.gz'),  D('Dataset_9_R2.fastq.gz')),
    ('D10', D('Dataset_10.fastq.gz'),    None),
    ('D11', D('Dataset_11_R1.fastq.gz'), D('Dataset_11_R2.fastq.gz')),
    ('D12', D('Dataset_12_R1.fastq.gz'), D('Dataset_12_R2.fastq.gz')),
    ('D13', D('Dataset_13_R1.fastq.gz'), D('Dataset_13_R2.fastq.gz')),
    ('D14', D('Dataset_14_R1.fastq.gz'), D('Dataset_14_R2.fastq.gz')),
    ('D15', D('Dataset_15_R1.fastq.gz'), D('Dataset_15_R2.fastq.gz')),
    ('D16', D('Dataset_16_R1.fastq.gz'), D('Dataset_16_R2.fastq.gz')),
    ('D17', D('Dataset_17_R1.fastq.gz'), D('Dataset_17_R2.fastq.gz')),
    ('D18', D('Dataset_18_R1.fastq.gz'), D('Dataset_18_R2.fastq.gz')),
]
ENGINES = ['salmon', 'minibwa']
_DATASET_NAMES = {name for name, _r1, _r2 in DATASETS}
_ENGINE_NAMES = set(ENGINES)


def tree_peak(proc, stop):
    import psutil
    peak = 0
    while not stop.is_set():
        try:
            total = proc.memory_info().rss
            for c in proc.children(recursive=True):
                try:
                    total += c.memory_info().rss
                except Exception:
                    pass
            peak = max(peak, total)
        except Exception:
            pass
        stop.wait(1.0)
    return peak


def run_engine(entry, ds_name, r1, r2, engine, threads=19):
    # 白名单校验：命令行/路径只允许常量表里的数据集与引擎
    if ds_name not in _DATASET_NAMES:
        raise ValueError(f'未知数据集: {ds_name}')
    if engine not in _ENGINE_NAMES:
        raise ValueError(f'未知引擎: {engine}')
    out_dir = _safe_path(OUT_ROOT, f'{ds_name}_{engine}')
    os.makedirs(out_dir, exist_ok=True)
    log_p = os.path.join(out_dir, 'engine.log')
    base = ['--out', out_dir, '--reference', REF, '--ref-info', INFO,
            '--engine', engine, '--threads', str(threads),
            '--align-threads', str(threads), '--sample-name', ds_name.lower(),
            '--index-dir', LIB]
    sub = ['identify', '-1', r1] + (['-2', r2] if r2 else [])
    t_id = t_fl = None
    import psutil
    me = psutil.Process()
    stop = threading.Event()
    peak_box = {}

    def mon():
        peak_box['v'] = tree_peak(me, stop)

    th = threading.Thread(target=mon, daemon=True)
    th.start()
    t0 = time.time()
    err = None
    try:
        with open(log_p, 'w', encoding='utf-8') as lf:
            t1 = time.time()
            subprocess.run(entry + sub + base, stdout=lf,
                           stderr=subprocess.STDOUT, check=True, shell=False)
            t_id = round(time.time() - t1, 1)
            t2 = time.time()
            subprocess.run(entry + ['filter'] + base, stdout=lf,
                           stderr=subprocess.STDOUT, check=True, shell=False)
            t_fl = round(time.time() - t2, 1)
    except Exception as e:
        err = f'{type(e).__name__}: {e}'
    stop.set()
    th.join(timeout=3)
    wall = round(time.time() - t0, 1)

    # 汇总过滤表
    n_pass = 0
    rows = []
    fp = os.path.join(out_dir, 'filter', 'filtered.tsv')
    if os.path.isfile(fp):
        with open(fp, encoding='utf-8') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        n_pass = len(rows)
    res = {'ds': ds_name, 'engine': engine, 'wall_sec': wall,
           'identify_sec': t_id, 'filter_sec': t_fl,
           'peak_rss_gb': round(peak_box.get('v', 0) / 2**30, 2),
           'passed': n_pass, 'error': err}
    with open(os.path.join(out_dir, 'bench.json'), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=2)
    # 保存定量明细供对比
    with open(os.path.join(out_dir, 'quant_rows.json'), 'w', encoding='utf-8') as fh:
        json.dump(rows, fh, ensure_ascii=False)
    print(f">>> {ds_name}/{engine}: {wall}s (identify {t_id}s + filter {t_fl}s) "
          f"peak {res['peak_rss_gb']}GB pass {n_pass}" + (f" ERR {err}" if err else ''), flush=True)
    return res


def main():
    from Virus_Platform_Core.kv_stage import engine_cmd
    entry_exe, entry_pre = engine_cmd()
    cmd_prefix = entry_pre if entry_pre else []
    # entry_exe 是 python 解释器；cmd_prefix 是 -c 引导
    full_entry = [entry_exe] + cmd_prefix
    only = sys.argv[1:] or None
    results = []
    for ds, r1, r2 in DATASETS:
        if only and ds not in only:
            continue
        for eng in ENGINES:
            try:
                results.append(run_engine(full_entry, ds, r1, r2, eng))
            except Exception as e:
                print(f"!!! {ds}/{eng} 异常: {e}", flush=True)
    print('\n===== 汇总 =====')
    for r in results:
        print(json.dumps(r, ensure_ascii=False))


if __name__ == '__main__':
    main()
