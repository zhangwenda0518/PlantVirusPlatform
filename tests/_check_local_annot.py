# -*- coding: utf-8 -*-
"""本地引擎端到端检查：走 Flask test client 提交 blastn/blastx（engine=local）。

不依赖运行中的服务；结果与前端「命中表」同源。
（原 cdd 分支已随 CDD 下线移除，2026-09-29。）
"""
import io
import os
import sys
import time
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

FA = os.path.join(ROOT, 'examples', 'example_contig_1.fasta')
CONTIG = 'example_contig_1'


def read_body(path):
    out = []
    for ln in io.open(path, encoding='utf-8'):
        ln = ln.strip()
        if ln and not ln.startswith('>'):
            out.append(ln)
    return ''.join(out)


def main() -> int:
    import app as appmod
    c = appmod.app.test_client()
    seq = read_body(FA)
    print(f'示例 contig: {CONTIG} ({len(seq)} bp)')
    fails = []
    for action in ('blastn', 'blastx'):
        r = c.post('/api/tool/analyze', json={
            'run': '_seq_input', 'contig': CONTIG, 'action': action,
            'seq': seq, 'engine': 'local'})
        if r.status_code != 200:
            print(f'  ✗ {action}: HTTP {r.status_code} {r.get_data()[:200]}')
            fails.append(action)
            continue
        d = r.get_json()
        if d.get('cached'):
            print(f'  {action}: 命中缓存')
        else:
            tid = d['task']
            t0 = time.time()
            while True:
                s = c.get(f'/api/task/{tid}').get_json()
                if s['status'] in ('done', 'failed', 'cancelled'):
                    break
                if time.time() - t0 > 300:
                    print(f'  ✗ {action}: 超时')
                    fails.append(action)
                    break
                time.sleep(1)
            if s['status'] != 'done':
                print(f'  ✗ {action}: {s["status"]} {s.get("error")}')
                fails.append(action)
                continue
            print(f'  {action}: 任务完成 {time.time() - t0:.1f}s')

        rr = c.get('/api/tool/analysis', query_string={
            'contig': CONTIG, 'action': action, 'engine': 'local', 'seq': seq})
        if rr.status_code != 200:
            print(f'  ✗ {action}: 读结果 HTTP {rr.status_code}')
            fails.append(action)
            continue
        res = rr.get_json()
        hits = res.get('hits') or []
        print(f'  {action}: engine={res.get("engine")} hits={len(hits)}')
        for h in hits[:3]:
            print(f"        {h['accession']:16s} {h['identity']:6.1f}% "
                  f"{str(h['title'])[:46]:46s} {h['sciname'][:28]}")
        if not hits:
            fails.append(f'{action}(无命中)')
    print('ALL OK' if not fails else 'FAILED: ' + ', '.join(fails))
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
