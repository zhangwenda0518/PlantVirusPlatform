# -*- coding: utf-8 -*-
"""为「同源比对（BLASTN · BLASTX）」生成示例结果。

默认走**本地引擎**（离线：blastn vs 病毒参考核酸库 /
DIAMOND vs RefSeq 病毒蛋白库），秒级完成；`--engine online` 可改用 NCBI 在线
（分钟级、需联网）。产物固化到 examples/results/hom/，
与真实运行结果隔离。
（原 cdd 模块已随 CDD 下线移除，2026-09-29。）

用法（需平台已在 http://127.0.0.1:8765 运行）：
    python tests/make_example_annot.py [--engine local|online]
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import sys
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.environ.get('VP_BASE', 'http://127.0.0.1:8765')
EX = os.path.join(ROOT, 'examples')
OUT = os.path.join(EX, 'results')
CONTIG_FA = os.path.join(EX, 'example_contig_1.fasta')

TITLES = {
    'hom': '同源比对（BLASTN / BLASTX）',
}


def _post(path: str, payload: dict) -> dict:
    req = urllib.request.Request(
        BASE + path, data=json.dumps(payload).encode('utf-8'),
        headers={'Content-Type': 'application/json'}, method='POST')
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode('utf-8'))


def _get(path: str):
    with urllib.request.urlopen(BASE + path, timeout=60) as r:
        return json.loads(r.read().decode('utf-8'))


def first_record(path: str):
    name, buf = None, []
    for ln in io.open(path, encoding='utf-8'):
        ln = ln.strip()
        if ln.startswith('>'):
            if name:
                break
            name = ln[1:].split()[0]
        elif name:
            buf.append(ln)
    return name, ''.join(buf)


def run_one(contig: str, seq: str, action: str, engine: str) -> str:
    """提交分析并等待完成，返回缓存 JSON 路径。"""
    d = _post('/api/tool/analyze',
              {'run': '_seq_input', 'contig': contig, 'action': action,
               'seq': seq, 'engine': engine})
    if d.get('error'):
        raise RuntimeError(f'{action}: 启动失败 {d["error"]}')
    tid = d.get('task')
    if not tid:
        print(f'    {action}: 已有缓存，跳过运行')
    else:
        t0 = time.time()
        while True:
            snap = _get('/api/task/' + tid)
            st = snap.get('status')
            if st == 'done':
                break
            if st in ('failed', 'cancelled'):
                raise RuntimeError(f'{action}: 任务 {st} {snap.get("error")}')
            if time.time() - t0 > 3600:
                raise RuntimeError(f'{action}: 超时')
            time.sleep(3)
        print(f'    {action}: 完成（{time.time() - t0:.0f}s）')
    key = hashlib.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
    safe = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
    tag = '' if engine == 'online' else '_' + engine
    stem = f'{safe}_{key}{tag}'
    p = os.path.join(ROOT, 'run', 'tool_runs', '_seq_input', stem, 'analysis',
                     f'{stem}_{action}.json')
    if not os.path.isfile(p):
        raise RuntimeError(f'{action}: 缓存文件不存在 {p}')
    return p


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine', choices=('local', 'online'), default='local')
    args = ap.parse_args()
    if not os.path.isfile(CONTIG_FA):
        print(f'缺少示例 contig 文件: {CONTIG_FA}')
        return 2
    contig, seq = first_record(CONTIG_FA)
    print(f'示例 contig: {contig} ({len(seq)} bp) · 引擎: {args.engine}')

    jobs = [('blastn', 'hom'), ('blastx', 'hom')]
    made = {}
    for action, module in jobs:
        try:
            src = run_one(contig, seq, action, args.engine)
        except Exception as e:                       # noqa: BLE001
            print(f'  ! {action} 失败: {e}')
            continue
        dst_dir = os.path.join(OUT, module)
        os.makedirs(dst_dir, exist_ok=True)
        name = f'{contig}_{action}_{args.engine}.json'
        shutil.copy2(src, os.path.join(dst_dir, name))
        made.setdefault(module, []).append(name)
        print(f'  → examples/results/{module}/{name}'
              f'  ({os.path.getsize(src) / 1024:.1f} KB)')

    mp = os.path.join(OUT, 'manifest.json')
    man = json.load(io.open(mp, encoding='utf-8')) if os.path.isfile(mp) else {}
    for module, files in made.items():
        d = os.path.join(OUT, module)
        allf = []
        for dirpath, _dirs, names in os.walk(d):
            for n in names:
                allf.append(os.path.relpath(os.path.join(dirpath, n), d)
                            .replace('\\', '/'))
        man[module] = {
            'title': TITLES.get(module, module),
            'files': sorted(allf),
            'generated_at': time.strftime('%Y-%m-%d %H:%M:%S'),
            'source': 'examples/example_contig_1.fasta',
            'note': ('本地引擎（离线）：blastn vs 病毒参考核酸库 / '
                     'DIAMOND vs RefSeq 病毒蛋白库'
                     if args.engine == 'local' else
                     'NCBI 在线（BLASTN / BLASTX）'),
        }
    with io.open(mp, 'w', encoding='utf-8') as f:
        json.dump(man, f, ensure_ascii=False, indent=1, sort_keys=True)
    print(f'manifest 已更新：{len(man)} 个模块')
    return 0


if __name__ == '__main__':
    sys.exit(main())
