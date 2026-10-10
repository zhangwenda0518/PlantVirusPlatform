# -*- coding: utf-8 -*-
"""引物设计守卫：保守区模式全链（建树组→保守区→primer3→热力学评分/分级）。

背景（2026-10-04 审查 + 2026-10-06 复核）：Primer 模块 wiki 强调的热力学
评分 / recommendation 分级 / 宿主特异性此前无任何专项测试（仅
_it_annotate 顺带断言 n_primers>0）；且示例快照一度 n_primers=0
（均匀点突变把连续保守段打散到 233bp < 400 阈值，2026-10-06 示例数据
改为块状突变后恢复 489bp/9 对）。本守卫经 /api/tool/run 与用户点击完全
同链路跑一遍，钉住全链不退化。

用法: python tests/_check_primer.py
"""
import csv
import os
import shutil
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

FAILS = []


def check(cond, msg):
    tag = 'ok  ' if cond else 'FAIL'
    print(f'  {tag} {msg}')
    if not cond:
        FAILS.append(msg)


def main():
    import app as appmod
    c = appmod.app.test_client()
    r = c.post('/api/tool/run', json={
        'tool': 'primer',
        'params': {'fasta': os.path.join(ROOT, 'examples',
                                         'example_conserved_set.fasta'),
                   'mode': 'conserved'}})
    j = r.get_json() or {}
    tid, out_dir = j.get('task'), j.get('out_dir')
    check(r.status_code == 200 and bool(tid and out_dir),
          f'任务提交 200（task={tid}）')

    status = None
    for _ in range(120):
        time.sleep(1.0)
        snap = c.get('/api/tasks?logs=1').get_json()
        tasks = snap if isinstance(snap, list) else (snap or {}).get('tasks') or []
        rec = next((t for t in tasks if t.get('id') == tid), None)
        if rec and rec.get('status') in ('done', 'failed', 'cancelled'):
            status = rec.get('status')
            break
    check(status == 'done', f'任务完成（{status}）')
    if status != 'done':
        print('✘ 任务未完成，终止')
        return 1

    summary = os.path.join(out_dir, '06_primer', 'summary.json')
    import json
    with open(summary, encoding='utf-8') as f:
        n = json.load(f).get('n_primers', 0)
    check(n >= 3, f'设计出引物对 ≥3（实际 {n}）')

    tsv = os.path.join(out_dir, '06_primer', 'primers.tsv')
    check(os.path.isfile(tsv), 'primers.tsv 产出')
    if os.path.isfile(tsv):
        rows = list(csv.DictReader(open(tsv, encoding='utf-8-sig'),
                                   delimiter='\t'))
        check(len(rows) == n, f'表行数 = n_primers（{len(rows)} vs {n}）')
        scores = [float(r['score']) for r in rows
                  if r.get('score') not in (None, '', 'NA')]
        check(len(scores) == len(rows), '每对引物都有热力学 score')
        if scores:
            check(all(0 <= s <= 100 for s in scores),
                  f'score 值域 0-100（{min(scores):.1f}-{max(scores):.1f}）')
        recs = {r.get('recommendation') for r in rows}
        check(recs <= {'RECOMMENDED', 'USABLE', 'RISKY', 'REJECT', 'NA'},
              f'recommendation 只取已知分级（{sorted(recs)}）')
        check(all(18 <= len(r['F_primer(5\'-3\')']) <= 30
                  and 18 <= len(r['R_primer(5\'-3\')']) <= 30
                  for r in rows), '引物长度 18-30nt')
        check(all(float(r['F_Tm']) >= 50 and float(r['R_Tm']) >= 50
                  for r in rows), 'Tm ≥ 50℃')
    shutil.rmtree(out_dir, ignore_errors=True)   # 守卫产物即跑即清

    if FAILS:
        print(f'✘ {len(FAILS)} 项未通过')
        return 1
    print('✓ 引物设计守卫全部通过')
    return 0


if __name__ == '__main__':
    sys.exit(main())
