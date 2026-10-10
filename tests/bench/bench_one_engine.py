# -*- coding: utf-8 -*-
"""单引擎跑批（带已完成跳过）：用法 python bench_one_engine.py <engine> <threads>"""
import json
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
# 平台根：默认自动推导（tests/bench/ 上溯 3 层 = 平台根），
# 可用环境变量 VP_ROOT 覆盖（换机/换目录不用改代码）。
PLATFORM = (os.environ.get('VP_ROOT')
            or os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, PLATFORM)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bench_engines import DATASETS, ENGINES, run_engine  # noqa: E402

engine = sys.argv[1] if len(sys.argv) > 1 else 'salmon'
threads = int(sys.argv[2]) if len(sys.argv) > 2 else 19
assert engine in ENGINES, f'引擎仅限 {ENGINES}'

from Virus_Platform_Core.kv_stage import engine_cmd  # noqa: E402
entry_exe, entry_pre = engine_cmd()
full_entry = [entry_exe] + entry_pre

for ds, r1, r2 in DATASETS:
    bench_p = os.path.join(PLATFORM, '_bench_kvs2', f'{ds}_{engine}', 'bench.json')
    bench_p = os.path.normpath(bench_p)
    if not bench_p.startswith(os.path.normpath(PLATFORM)):
        raise ValueError(f'路径越界: {bench_p}')
    if os.path.isfile(bench_p):
        try:
            prev = json.load(open(bench_p, encoding='utf-8'))
            if not prev.get('error'):
                print(f'-- 跳过 {ds}/{engine}（已完成 {prev["wall_sec"]}s）', flush=True)
                continue
        except Exception:
            pass
    try:
        run_engine(full_entry, ds, r1, r2, engine, threads=threads)
    except Exception as e:
        print(f'!!! {ds}/{engine} 异常: {e}', flush=True)
print(f'[{engine}] 跑批结束', flush=True)
