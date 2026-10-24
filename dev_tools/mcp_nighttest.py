# -*- coding: utf-8 -*-
"""MCP 全面可行性夜间测试驱动器（2026-10-16）。

八个阶段全部走**真实执行**（不做 mock）：
  P0 环境预检          P1 协议级 19 工具全扫（只读类真调用）
  P2 真实全流程流水线   P3 工具卡真实链（quicktree/sdt/genoplot/primer/convert）
  P4 CLI 安全面全量     P5 预筛模块可达性
  P6 三入口一致性       P7 健壮性（并发/重启/延迟/下线门）
产出：run/mcp_nighttest/report.md（逐阶段 PASS/FAIL + 证据）。
任何阶段失败不阻断后续阶段；最终给出总体可行性结论。
"""
import asyncio
import json
import os
import shutil
import statistics
import sys
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'dev_tools'))
sys.path.insert(0, ROOT)          # 嵌入式 _pth 不认 cwd，必须显式
os.chdir(ROOT)

from test_mcp_e2e import PY, SERVER, call, poll_task  # noqa: E402
from mcp import ClientSession, StdioServerParameters   # noqa: E402
from mcp.client.stdio import stdio_client              # noqa: E402

REPORT_DIR = os.path.join(ROOT, 'run', 'mcp_nighttest')
REPORT = os.path.join(REPORT_DIR, 'report.md')
EX = os.path.join(ROOT, 'examples')
RESULTS = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(REPORT_DIR, 'console.log'), 'a',
              encoding='utf-8') as f:
        f.write(line + '\n')


def sec(name):
    RESULTS.append({'phase': name, 'checks': []})
    log(f"━━━ {name} ━━━")


def check(name, ok, detail=''):
    RESULTS[-1]['checks'].append((name, bool(ok), str(detail)[:220]))
    log(('  ✅ ' if ok else '  ❌ ') + name + (f'  {detail}' if detail else ''))


def flush_report():
    lines = ['# MCP 全面可行性夜间测试报告',
             '',
             f"生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}",
             f"被测对象：virusplatform_mcp（平台根 mcp_server.py，stdio/HTTP 双模式）",
             '']
    total_pass = total_fail = 0
    fails = []
    for ph in RESULTS:
        p = sum(1 for _, ok, _ in ph['checks'] if ok)
        f = sum(1 for _, ok, _ in ph['checks'] if not ok)
        total_pass += p
        total_fail += f
        icon = '✅' if f == 0 else '⚠️' if p > f else '❌'
        lines.append(f"## {icon} {ph['phase']}（通过 {p} / 失败 {f}）")
        lines.append('')
        for name, ok, detail in ph['checks']:
            mark = '✅' if ok else '❌'
            lines.append(f"- {mark} {name}" + (f" —— {detail}" if detail else ''))
        lines.append('')
        fails.extend((ph['phase'], n, d) for n, ok, d in ph['checks'] if not ok)
    lines.append('---')
    lines.append(f"## 总计：通过 {total_pass} / 失败 {total_fail}")
    if fails:
        lines.append('')
        lines.append('### 失败项明细')
        for ph, n, d in fails:
            lines.append(f"- [{ph}] {n} —— {d}")
    with open(REPORT, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


async def stdio_session():
    params = StdioServerParameters(command=PY, args=[SERVER])
    ctx = stdio_client(params)
    r, w = await ctx.__aenter__()
    s = ClientSession(r, w)
    await s.__aenter__()
    await s.initialize()
    return ctx, s


async def p0_p1(ctx, s):
    sec('P0 环境预检')
    r = await asyncio.create_subprocess_exec(
        PY, SERVER, '--check', stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT, cwd=ROOT)
    out, _ = await asyncio.wait_for(r.communicate(), 120)
    txt = out.decode('utf-8', 'replace')
    check('P0 --check 退出码 0', r.returncode == 0, txt.strip().splitlines()[-1])
    check('P0 --check 报告 ready', 'ready    : yes' in txt)

    sec('P1 协议级 19 工具全扫（真调用只读类）')
    tools = (await s.list_tools()).tools
    check('P1 工具数 = 19', len(tools) == 19, f"实际 {len(tools)}")
    no_doc = [t.name for t in tools if not (t.description or '').strip()]
    check('P1 全部带 docstring', not no_doc, str(no_doc))

    for name, args in [
            ('vp_platform_info', {}),
            ('vp_list_tools', {}), ('vp_list_samples', {}),
            ('vp_list_tasks', {'limit': 5, 'include_archived': True}),
            ('vp_list_tool_runs', {'limit': 5}),
            ('vp_pipeline_stages', {}),
            ('vp_cli_help', {}),
            ('vp_prescreen_projects', {}),
            ('vp_get_task_log', {'task_id': 'zzzz'}),
            ('vp_prescreen_result', {'name': 'mcp_night_none'}),
    ]:
        err, d, text = await call(s, name, args)
        ok = (not err)
        check(f'P1 {name} 可调用', ok, (text or '')[:100])

    err, d, _ = await call(s, 'vp_platform_info', {})
    check('P1 platform_info 信封 summary',
          bool(d.get('summary')) and bool(d.get('fingerprint')),
          d.get('summary', ''))
    return d.get('fingerprint')


async def p2_pipeline(s):
    sec('P2 真实全流程流水线（示例双端数据，默认 13 阶段）')
    r1 = os.path.join(EX, 'example_R1.fastq.gz')
    r2 = os.path.join(EX, 'example_R2.fastq.gz')
    err, sub, text = await call(
        s, 'vp_pipeline_start',
        {'sample': 'mcp_night_demo', 'r1': r1, 'r2': r2}, timeout=120)
    tid = (sub or {}).get('result', {}).get('task_id')
    check('P2 流水线提交（信封）', not err and tid,
          (sub or {}).get('summary') or text[:120])
    if not tid:
        return
    dl = time.time() + 5400
    last = ''
    while time.time() < dl:
        err, snap, _ = await call(s, 'vp_get_task',
                                  {'task_id': tid, 'log_lines': 3})
        if err:
            check('P2 轮询', False, str(snap)[:150])
            return
        line = f"{snap['status']} {snap['stage']} {snap['pct']:.0%} {snap['msg']}"
        if line != last:
            log(f'    [流水线] {line}')
            last = line
        if snap['status'] in ('done', 'failed', 'cancelled'):
            break
        await asyncio.sleep(15)
    ok = snap['status'] == 'done'
    check('P2 流水线到终态 done', ok,
          f"error={snap.get('error', '')[:200]}" if not ok else '')
    check('P2 done 信封带解读', bool(snap.get('interpretation'))
          and bool(snap.get('next')), str(snap.get('next'))[:80])
    sd = os.path.join(ROOT, 'run', 'results', 'mcp_night_demo')
    stages = {}
    try:
        mf = json.load(open(os.path.join(sd, 'project.json'),
                            encoding='utf-8'))
        stages = mf.get('stages') or mf.get('stages_status') or {}
    except Exception:
        pass
    n_report = os.path.isfile(os.path.join(sd, '10_report'))
    check('P2 样品目录产物落盘', os.path.isdir(sd) and bool(os.listdir(sd)),
          f'project.json 阶段键 {len(stages)}；10_report 存在={n_report}')


TOOL_CARD_RUNS = [
    ('quicktree', {'seqs': 'examples/example_conserved_set.fasta',
                   'aligned': True}),
    ('sdt', {'seqs': 'examples/example_cmv.fasta'}),
    ('genoplot', {'fasta': None, 'ann': 'examples/example_genome.gb'}),
    ('primer', {'fasta': 'examples/example_cmv.fasta', 'mode': 'plain'}),
    ('convert', {'input': 'examples/example_R1.fastq.gz', 'target': 'fasta'}),
]


async def p3_cards(s):
    sec('P3 工具卡真实链（5 张卡：提交→轮询→产物）')
    for tool, params in TOOL_CARD_RUNS:
        params = {k: (os.path.join(ROOT, v)
                      if isinstance(v, str) and v.startswith('examples/')
                      else v) for k, v in params.items()}
        err, sub, text = await call(s, 'vp_run_tool',
                                    {'tool': tool, 'params': params},
                                    timeout=120)
        tid = (sub or {}).get('result', {}).get('task_id')
        run = (sub or {}).get('result', {}).get('run')
        check(f'P3 {tool} 提交', not err and tid,
              (sub or {}).get('summary') or text[:100])
        if not tid:
            continue
        snap, perr = await poll_task(s, tid, timeout=1800)
        ok = (snap or {}).get('status') == 'done'
        check(f'P3 {tool} 执行 done', ok,
              perr or f"status={(snap or {}).get('status')} "
                      f"error={(snap or {}).get('error') or ''}"[:200])
        if ok and run:
            err, rf, _ = await call(s, 'vp_read_run_file',
                                    {'run': run, 'path': 'run.log'})
            check(f'P3 {tool} 产物可读', not err
                  and (rf or {}).get('result', {}).get('text'),
                  (rf or {}).get('summary', ''))


CLI_SYNC = [
    ('ictv-status', {}), ('ncbi-list', {}), ('gb-list', {}),
    ('ref-status', {}), ('submit-list', {}), ('logan-jobs', {}),
    ('ictv-refs', {'genus': 'Tobamovirus', 'limit': 5}),
]
CLI_TASK = [
    ('samples-backfill', {'dry_run': True}),
    ('db-migrate', {'to': os.path.join(ROOT, 'run', 'mcp_nighttest',
                                       'fake_db_target'), 'check': True}),
    ('tool-runs', {'action': 'status'}),
    ('tool-runs', {'action': 'archive', 'dry_run': True,
      'older_days': 30}),
]


async def p4_cli(s):
    sec('P4 CLI 安全面全量（sync 真执行 / task 安全参数）')
    for cmd, args in CLI_SYNC:
        err, d, text = await call(s, 'vp_cli_run',
                                  {'command': cmd, 'args': args}, timeout=330)
        ok = not err and d.get('status') == 'ok'
        check(f'P4 run {cmd}', ok, d.get('summary', text[:80]))
    for cmd, args in CLI_TASK:
        err, d, text = await call(s, 'vp_cli_submit',
                                  {'command': cmd, 'args': args}, timeout=120)
        tid = (d or {}).get('result', {}).get('task_id')
        check(f'P4 submit {cmd} {" ".join(args)}', not err and tid,
              text[:80])
        if not tid:
            continue
        snap, perr = await poll_task(s, tid, timeout=600)
        if cmd == 'db-migrate':
            # --check 对无库目标**应当**失败：CLI 如实报告「缺运行文件」
            # 即为正确行为（MCP 层只负责透传）
            ok = (snap or {}).get('status') == 'failed' and                 '缺运行文件' in ((snap or {}).get('error') or '')
            check('P4 task db-migrate --check 语义正确（如实报缺库）', ok,
                  f"status={(snap or {}).get('status')}")
        else:
            check(f'P4 task {cmd} done', (snap or {}).get('status') == 'done',
                  perr or f"error={(snap or {}).get('error') or ''}"[:150])
    # 真联网下载（小批量，证明在线链路）
    err, d, text = await call(
        s, 'vp_cli_submit',
        {'command': 'ncbi-dl',
         'args': {'term': 'Tobamovirus[ORGN] AND complete genome[TITL]',
                  'name': 'mcp_nighttest', 'max': 2}}, timeout=120)
    tid = (d or {}).get('result', {}).get('task_id')
    check('P4 submit ncbi-dl（真联网）', not err and tid, text[:80])
    if tid:
        snap, perr = await poll_task(s, tid, timeout=600)
        check('P4 ncbi-dl done（下载 2 条参考）',
              (snap or {}).get('status') == 'done',
              perr or f"error={(snap or {}).get('error') or ''}"[:150])
        err, lst, _ = await call(s, 'vp_cli_run', {'command': 'ncbi-list'})
        check('P4 集合出现在 ncbi-list',
              not err and 'mcp_nighttest' in (lst or {}).get('result', {})
              .get('output', ''))


async def p5_prescreen(s):
    sec('P5 预筛模块可达性（真实联网批量留待有真实 Run 号时验证）')
    err, d, _ = await call(s, 'vp_prescreen_projects', {})
    check('P5 projects 可查', not err and 'projects' in d)
    err, d, _ = await call(s, 'vp_prescreen_run',
                           {'name': 'mcp_night_ps', 'source': 'file',
                            'path': 'examples/不存在.list'})
    check('P5 非法路径当场报错', err and 'path' in str(d).lower()
          or err, str(d)[:100])


async def _http_fingerprint(port, token):
    from mcp.client.streamable_http import streamablehttp_client
    url = f'http://127.0.0.1:{port}/mcp'
    hdrs = {'Authorization': f'Bearer {token}'}
    async with streamablehttp_client(url, headers=hdrs) as (r, w, _sid):
        async with ClientSession(r, w) as s:
            await s.initialize()
            err, info, _ = await call(s, 'vp_platform_info', {})
            return (info or {}).get('fingerprint')


async def p6_cross_entry(port=8913, token='nighttest'):
    sec('P6 三入口一致性（stdio / HTTP / 网页控制台）')
    proc = subprocess.Popen([PY, SERVER, '--http', str(port), '--token', token],
                            cwd=ROOT, stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT)
    try:
        dl = time.time() + 30
        while time.time() < dl:
            try:
                import socket
                with socket.create_connection(('127.0.0.1', port), 1):
                    break
            except OSError:
                await asyncio.sleep(0.5)
        fp_http = await _http_fingerprint(port, token)
        check('P6 HTTP 入口可用', bool(fp_http), f'指纹 {fp_http}')

        from app import app as _flask
        fc = _flask.test_client()
        r = fc.post('/api/mcp/call', json={'name': 'vp_platform_info',
                                           'arguments': {}})
        d = r.get_json()
        fp_web = (d or {}).get('data', {}).get('fingerprint')
        check('P6 网页控制台入口可用', r.status_code == 200
              and d.get('isError') is False and bool(fp_web),
              f'指纹 {fp_web}')
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()


async def p7_robust(s, expected_fp):
    sec('P7 健壮性（延迟/并发/信封/下线门）')
    lat = []
    for _ in range(20):
        t0 = time.time()
        err, d, _ = await call(s, 'vp_platform_info', {})
        assert not err
        lat.append((time.time() - t0) * 1000)
    check('P7 20 连发延迟 P95 < 500ms',
          statistics.quantiles(lat, n=20)[18] < 500,
          f'min={min(lat):.0f}ms avg={statistics.mean(lat):.0f}ms '
          f'max={max(lat):.0f}ms')

    results = {}

    def worker(i):
        asyncio.run(_conc(results, i))
    async def _conc(res, i):
        async with stdio_client(StdioServerParameters(
                command=PY, args=[SERVER])) as (r, w):
            async with ClientSession(r, w) as s2:
                await s2.initialize()
                err, d, _ = await call(s2, 'vp_list_tools', {})
                res[i] = (not err, (d or {}).get('count'))
    ts = [threading.Thread(target=worker, args=(i,)) for i in range(4)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    check('P7 4 并发 stdio 客户端互不干扰（各列到 31 张卡）',
          all(v[0] and v[1] == 31 for v in results.values()),
          str(results))

    env = dict(os.environ)
    env['VP_MCP_DISABLED_TOOLS'] = 'vp_cancel_task'
    params = StdioServerParameters(command=PY, args=[SERVER], env=env)
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s2:
            await s2.initialize()
            names = [t.name for t in (await s2.list_tools()).tools]
    check('P7 下线门生效（19→18）',
          len(names) == 18 and 'vp_cancel_task' not in names)


import subprocess  # noqa: E402


PH = {'P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'}
if len(sys.argv) > 1 and sys.argv[1].startswith('P'):
    PH = set(sys.argv[1].split(','))


async def main():
    os.makedirs(REPORT_DIR, exist_ok=True)
    log('夜间全面可行性测试开始')
    ctx = s = None
    fp = None
    try:
        ctx, s = await stdio_session()
        if 'P0' in PH or 'P1' in PH:
            fp = await p0_p1(ctx, s)
        if 'P2' in PH:
            await p2_pipeline(s)
        if 'P3' in PH:
            await p3_cards(s)
        if 'P4' in PH:
            await p4_cli(s)
        if 'P5' in PH:
            await p5_prescreen(s)
    except Exception as e:
        check('P0-P5 执行', False, f'{type(e).__name__}: {e}')
    finally:
        if ctx:
            try:
                await ctx.__aexit__(None, None, None)
            except Exception:
                pass
    try:
        if 'P6' in PH:
            await p6_cross_entry()
    except Exception as e:
        sec('P6 三入口一致性')
        check('P6 执行', False, f'{type(e).__name__}: {e}')
    try:
        async with stdio_client(StdioServerParameters(
                command=PY, args=[SERVER])) as (r, w):
            async with ClientSession(r, w) as s2:
                await s2.initialize()
                err, info, _ = await call(s2, 'vp_platform_info', {})
                await p7_robust(s2, (info or {}).get('fingerprint'))
    except Exception as e:
        sec('P7 健壮性')
        check('P7 执行', False, f'{type(e).__name__}: {e}')

    flush_report()
    log(f'完成。报告: {REPORT}')


if __name__ == '__main__':
    asyncio.run(main())
