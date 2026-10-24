# -*- coding: utf-8 -*-
"""mcp_server.py stdio 全链路自测（不依赖平台包，仅 mcp SDK + 标准库）。

用法（用平台内置 Python 跑，客户端与服务端同为 mcp v1）：
    3rd\\python\\python.exe dev_tools\\test_mcp_e2e.py

覆盖：握手 → 工具清单 → 环境总览 → 工具箱/阶段/样品/运行清单 →
未知工具与缺参错误路径 → convert 真实运行（提交→轮询→读产物）。
加 --pipeline 追加流水线冒烟；加 --http [port] 测 streamable-HTTP 模式
（客户端断开后任务照跑、新连接重连轮询到 done）。
"""
import asyncio
import json
import os
import socket
import subprocess
import sys
import tempfile
import time

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, '3rd', 'python', 'python.exe')
SERVER = os.path.join(ROOT, 'mcp_server.py')

PASS, FAIL = [], []


def check(name, ok, detail=''):
    (PASS if ok else FAIL).append(name)
    print(('  ✅ ' if ok else '  ❌ ') + name + (f'  {detail}' if detail else ''))


async def call(s, tool, args, timeout=60):
    res = await asyncio.wait_for(s.call_tool(tool, args or {}), timeout)
    text = ''.join(getattr(c, 'text', '') for c in res.content)
    try:
        data = json.loads(text)
    except (ValueError, TypeError):
        data = None
    return bool(getattr(res, 'isError', False)), data, text


def _make_fastq():
    tmp = tempfile.mkdtemp(prefix='vp_mcp_e2e_')
    fq = os.path.join(tmp, 'tiny_R1.fastq')
    with open(fq, 'w', encoding='ascii') as f:
        for i in range(1, 5):
            f.write(f'@read{i}\nACGTACGTACGTACGTACGT\n+\n{"I" * 20}\n')
    return fq


async def test_http(port, token=None):
    """streamable-HTTP 模式：断开客户端 → 任务照跑 → 重连轮询到 done。
    给 token 时另验：无 Authorization 被拦 401，带 token 正常。"""
    from mcp.client.streamable_http import streamablehttp_client

    url = f'http://127.0.0.1:{port}/mcp'
    hdrs = {'Authorization': f'Bearer {token}'} if token else None
    proc = subprocess.Popen(
        [PY, SERVER, '--http', str(port)] +
        (['--token', token] if token else []), cwd=ROOT,
        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    try:
        dl = time.time() + 30
        while time.time() < dl:
            try:
                with socket.create_connection(('127.0.0.1', port), 1):
                    break
            except OSError:
                await asyncio.sleep(0.5)
        else:
            check('HTTP 服务端口就绪', False, '30s 内未监听')
            return

        print('== HTTP-0. 鉴权 ==')
        if token:
            rejected = False
            try:
                async with streamablehttp_client(url) as (r, w, _sid):
                    async with ClientSession(r, w) as s:
                        await s.initialize()
            except Exception:
                rejected = True
            check('无 token 被拦', rejected)
            ok_auth = False
            try:
                async with streamablehttp_client(url, headers=hdrs) as \
                        (r, w, _sid):
                    async with ClientSession(r, w) as s:
                        await s.initialize()
                ok_auth = True
            except Exception as e:
                print('    带token失败:', repr(e)[:120])
            check('带 token 通过', ok_auth)
            if not ok_auth:
                return
        print('== HTTP-1. 握手与提交 ==')
        async with streamablehttp_client(url, headers=hdrs) as (r, w, _sid):
            async with ClientSession(r, w) as s:
                await s.initialize()
                tools = (await s.list_tools()).tools
                check('HTTP list_tools 22 个工具', len(tools) == 22,
                      ','.join(sorted(t.name for t in tools)))
                err, info, _ = await call(s, 'vp_platform_info', {})
                check('HTTP vp_platform_info ready',
                      not err and info.get('ready'))
                err, sub, text = await call(
                    s, 'vp_run_tool',
                    {'tool': 'convert',
                     'params': {'input': _make_fastq(), 'target': 'fasta'},
                     'threads': 2})
                check('HTTP 提交 convert',
                      not err and sub.get('result', {}).get('task_id'),
                      text[:100])
                tid = sub['result']['task_id']
                run = sub['result']['run']

        print('== HTTP-2. 客户端已断开，等待任务在服务端跑完 ==')
        status_a = None
        async with streamablehttp_client(url, headers=hdrs) as (r, w, _sid):
            async with ClientSession(r, w) as s:
                await s.initialize()
                err, snap, _ = await call(
                    s, 'vp_get_task', {'task_id': tid, 'log_lines': 0})
                status_a = (snap or {}).get('status')
            # 立刻再断开，只在服务端轮询
        deadline = time.time() + 180
        final = None
        while time.time() < deadline:
            async with streamablehttp_client(url, headers=hdrs) as (r, w, _sid):
                async with ClientSession(r, w) as s:
                    await s.initialize()
                    err, snap, _ = await call(
                        s, 'vp_get_task', {'task_id': tid, 'log_lines': 0})
            final = (snap or {}).get('status')
            if final in ('done', 'failed', 'cancelled'):
                break
            await asyncio.sleep(2)
        check('客户端断开后任务到终态 done', final == 'done',
              f'断开时={status_a} 终态={final}')

        print('== HTTP-3. 重连读产物 ==')
        async with streamablehttp_client(url, headers=hdrs) as (r, w, _sid):
            async with ClientSession(r, w) as s:
                await s.initialize()
                err, rf, _ = await call(s, 'vp_read_run_file',
                                        {'run': run, 'path': 'run.log'})
                check('重连读 run.log',
                      not err and 'fq2fa'
                      in (rf or {}).get('result', {}).get('text', ''))
                err, data, _ = await call(s, 'vp_list_tool_runs', {'limit': 5})
                check('重连清单含本次运行',
                      not err and any(x['run'] == run
                                      for x in data.get('runs', [])))
                err, del_res, _ = await call(s, 'vp_delete_tool_run',
                                             {'run': run, 'confirm': True})
                check('重连删除测试运行', not err and del_res.get('deleted'))
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
    print(f'\nHTTP 模式：通过 {len(PASS)} / 失败 {len(FAIL)}')
    if FAIL:
        print('失败项:', ', '.join(FAIL))


async def poll_task(s, tid, timeout=180):
    dl = time.time() + timeout
    last = ''
    while time.time() < dl:
        err, snap, _ = await call(s, 'vp_get_task',
                                  {'task_id': tid, 'log_lines': 5})
        if err:
            return None, f'vp_get_task 报错: {snap}'
        line = f"{snap.get('status')} {snap.get('stage')} " \
               f"{snap.get('pct', 0):.0%} {snap.get('msg') or ''}"
        if line != last:
            print(f'    [{time.strftime("%H:%M:%S")}] {line.strip()}')
            last = line
        if snap.get('status') in ('done', 'failed', 'cancelled'):
            return snap, None
        await asyncio.sleep(1.5)
    return None, '轮询超时'


async def test_gate():
    """启动 gate：env 摘除工具、无效 deferred.json fail-closed 拒启。"""
    from mcp.client.stdio import get_default_environment
    print('== GATE. 下线机制 ==')
    env = get_default_environment()
    env['VP_MCP_DISABLED_TOOLS'] = 'vp_cancel_task,not_a_tool'
    params = StdioServerParameters(command=PY, args=[SERVER], env=env)
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            names = sorted(t.name for t in (await s.list_tools()).tools)
            check('env 下线摘除 vp_cancel_task', 'vp_cancel_task' not in names,
                  f'{len(names)} 个工具')
            err, _, text = await call(s, 'vp_cancel_task', {'task_id': 'x'})
            check('被下线工具不可调用',
                  err and ('未知工具' in text or 'Unknown tool' in text),
                  text[:60])

    dj = os.path.join(ROOT, 'deferred.json')
    with open(dj, 'w', encoding='utf-8') as f:
        f.write('{bad json!!')
    try:
        p = subprocess.run([PY, SERVER, '--check'], cwd=ROOT,
                           capture_output=True, text=True, timeout=120,
                           errors='replace')
        check('坏 deferred.json 拒启（fail-closed）',
              p.returncode != 0 and 'fail-closed' in (p.stderr or '')
              + (p.stdout or ''))
    finally:
        os.remove(dj)


async def main():
    params = StdioServerParameters(command=PY, args=[SERVER])
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            print('== 1. 握手与清单 ==')
            tools = (await s.list_tools()).tools
            names = sorted(t.name for t in tools)
            check('list_tools 返回 22 个工具', len(names) == 22, ','.join(names))
            no_desc = [t.name for t in tools if not (t.description or '').strip()]
            check('全部工具带 docstring 描述', not no_desc, str(no_desc))

            err, info, _ = await call(s, 'vp_platform_info', {})
            check('vp_platform_info ready', not err and info.get('ready'),
                  f"host_db={info.get('host_db_ready')} "
                  f"virus_db={info.get('virus_db_ready')}")
            check('版本+能力指纹齐备',
                  bool(info.get('version')) and bool(info.get('fingerprint'))
                  and info.get('tools_count') == 22
                  and len(info.get('tools_enabled') or []) == 22,
                  f"v={info and info.get('version')} "
                  f"fp={info and info.get('fingerprint')}")

            err, tl, _ = await call(s, 'vp_list_tools', {})
            # 2026-10-20：+1 = conan（CONAN-SWIFT Cas3 crRNA 设计）→ 32 张卡
            check('vp_list_tools 32 张卡', not err and tl.get('count') == 32,
                  f"count={tl and tl.get('count')}")
            tmap = {t['tool']: t for t in tl['tools']}
            no_schema = [k for k, t in tmap.items() if not t.get('params')]
            check('32 张卡全带参数 schema', not no_schema, str(no_schema))
            n_params = sum(len(t.get('params') or []) for t in tmap.values())
            check('参数 schema 覆盖 ≥120 项', n_params >= 120,
                  f'{n_params} 项')
            ks = next(p['desc'] for p in tmap['kvsuite']['params']
                      if p['key'] == 'stage')
            check('kvsuite.stage 带枚举口径', 'identify+filter' in ks)
            check('convert.required_params 含 input',
                  'input' in tmap['convert']['required_params'])
            check('conan.required_params 含 input，负组可选',
                  'input' in tmap.get('conan', {}).get('required_params', [])
                  and 'negative' not in tmap.get('conan', {})
                  .get('required_params', []))
            check('未下线时全部卡 gated=false',
                  all(t.get('gated') is False for t in tmap.values()))

            err, tl2, _ = await call(s, 'vp_list_tasks',
                                     {'limit': 5, 'include_archived': True})
            check('归档任务并入列表',
                  not err and isinstance(tl2.get('archived_count'), int)
                  and isinstance(tl2.get('archived'), list),
                  f"archived_count={tl2 and tl2.get('archived_count')}")

            err, st, _ = await call(s, 'vp_pipeline_stages', {})
            check('vp_pipeline_stages 15 阶段',
                  not err and len(st.get('stage_order', [])) == 15)

            for tname in ('vp_list_samples', 'vp_list_tool_runs',
                          'vp_list_tasks'):
                err, data, _ = await call(s, tname, {})
                check(tname, not err,
                      f"total={data and (data.get('total') or data.get('count'))}")

            print('== 2. 错误路径 ==')
            err, data, text = await call(s, 'vp_run_tool', {'tool': 'nope'})
            check('未知工具被拒', err and '未知工具' in text, text[:80])
            err, data, text = await call(s, 'vp_run_tool',
                                         {'tool': 'convert', 'params': {}})
            check('缺必填参数被拒', err and '缺少必填参数' in text, text[:80])
            err, data, text = await call(s, 'vp_get_task', {'task_id': 'zzzz'})
            check('未知任务被拒', err and '任务不存在' in text, text[:80])

            print('== 3. 真实运行（convert FASTQ→FASTA）==')
            tmp = tempfile.mkdtemp(prefix='vp_mcp_e2e_')
            fq = os.path.join(tmp, 'tiny_R1.fastq')
            with open(fq, 'w', encoding='ascii') as f:
                for i in range(1, 5):
                    f.write(f'@read{i}\nACGTACGTACGTACGTACGT\n+\n'
                            f'{"I" * 20}\n')
            err, sub, _ = await call(
                s, 'vp_run_tool',
                {'tool': 'convert',
                 'params': {'input': fq, 'target': 'fasta'}, 'threads': 2})
            check('提交成功返回 task_id/run（信封）',
                  not err and sub.get('result', {}).get('task_id'),
                  f"run={sub and sub.get('run')} weight={sub and sub.get('weight')}")
            if err or not sub:
                return
            snap, perr = await poll_task(s, sub['result']['task_id'])
            check('任务到终态 done', bool(snap) and
                  snap.get('status') == 'done',
                  perr or f"status={(snap or {}).get('status')} "
                          f"error={(snap or {}).get('error')}")
            expect_keys = {'id', 'name', 'status', 'stage', 'pct', 'msg',
                           'error', 'eta', 'out_dir', 'link', 'started_h',
                           'finished_h', 'log', 'result', 'summary',
                           'interpretation', 'next'}
            check('快照稳定投影（键集合恒定）',
                  set((snap or {}).keys()) == expect_keys,
                  str(sorted(set((snap or {}).keys()) ^ expect_keys)))
            out_dir = (snap or {}).get('result', {}).get('out_dir') or \
                sub['result']['out_dir']
            check('结果注入 out_dir', bool(out_dir), str(out_dir))

            err, rf, _ = await call(s, 'vp_read_run_file',
                                    {'run': sub['result']['run'],
                                     'path': 'run.log'})
            check('vp_read_run_file 读 run.log',
                  not err and 'fq2fa'
                  in (rf or {}).get('result', {}).get('text', ''),
                  f"{rf and rf.get('summary')}")

            err, rf2, _ = await call(s, 'vp_read_run_file',
                                     {'run': sub['result']['run'],
                                      'path': '../../etc/passwd'})
            check('路径越界被拒', err)

            err, runs, _ = await call(s, 'vp_list_tool_runs',
                                      {'limit': 5})
            check('运行出现在清单里',
                  not err and any(r['run'] == sub['result']['run']
                                  for r in runs.get('runs', [])))

            print('== 3b. 删除守卫与执行 ==')
            err, _, text = await call(s, 'vp_delete_tool_run',
                                      {'run': sub['result']['run'],
                                       'confirm': False})
            check('未 confirm 拒删', err and 'confirm' in text, text[:80])
            err, _, text = await call(s, 'vp_delete_tool_run',
                                      {'run': 'no_such_run_xyz',
                                       'confirm': True})
            check('不存在的运行拒删', err, text[:60])
            err, del_res, text = await call(s, 'vp_delete_tool_run',
                                            {'run': sub['result']['run'],
                                             'confirm': True})
            check('confirm 后删除成功', not err and del_res.get('deleted'),
                  text[:80])
            err, runs2, _ = await call(s, 'vp_list_tool_runs', {'limit': 5})
            check('删除后清单不再含该运行',
                  not err and not any(r['run'] == sub['result']['run']
                                      for r in runs2.get('runs', [])))

            print('== 3c. conan 真实运行（CONAN-SWIFT Cas3 crRNA 设计）==')
            # 合成已知答案：目标 2 株（位点 AAG+32nt）、负组 1 条（位点
            # 原样保留 → leak 率 100%）；产物断言 summary.json 的 top 行
            rng_c = __import__('random').Random(20261020)
            _r = lambda n: ''.join(rng_c.choice('ACGT') for _ in range(n))
            sp32 = 'ACGTTGCACTGGATCAAGTCAAGGTTACGATA'
            site_c = _r(199) + 'AAG' + sp32 + _r(166)
            conan_dir = tempfile.mkdtemp(prefix='vp_mcp_conan_')
            tfa = os.path.join(conan_dir, 'targets.fa')
            nfa = os.path.join(conan_dir, 'negatives.fa')
            with open(tfa, 'w', encoding='ascii') as f:
                f.write(f'>isoA\n{site_c}\n>isoB\n{site_c}\n')
            with open(nfa, 'w', encoding='ascii') as f:
                f.write(f'>near_leak\n{site_c}\n')
            err, sub3, _ = await call(
                s, 'vp_run_tool',
                {'tool': 'conan',
                 'params': {'input': tfa, 'negative': nfa}})
            check('conan 提交成功返回 task_id',
                  not err and sub3.get('result', {}).get('task_id'),
                  f"run={sub3 and sub3.get('run')}")
            if not err and sub3:
                snap3, perr3 = await poll_task(s, sub3['result']['task_id'])
                check('conan 任务到终态 done',
                      bool(snap3) and snap3.get('status') == 'done',
                      perr3 or f"status={(snap3 or {}).get('status')} "
                               f"error={(snap3 or {}).get('error')}")
                err4, rf3, _ = await call(
                    s, 'vp_read_run_file',
                    {'run': sub3['result']['run'],
                     'path': 'conan/summary.json'})
                ok_sm = False
                try:
                    sm3 = json.loads(
                        (rf3 or {}).get('result', {}).get('text', ''))
                    ok_sm = (sm3.get('n_candidates', 0) > 0
                             and bool(sm3.get('top'))
                             and sm3['top'][0].get('leak') == 100.0
                             and sm3['top'][0].get('all_targets_pass') is True)
                except Exception:
                    ok_sm = False
                check('conan summary 已知答案命中（leak=100，全目标通过）',
                      not err4 and ok_sm,
                      ((rf3 or {}).get('result', {}).get('text') or '')[:120])

            print('== 3d. 任务日志与取消路径 ==')
            err, lg, _ = await call(s, 'vp_get_task_log',
                                    {'task_id': sub['result']['task_id'],
                                     'lines': 50})
            check('vp_get_task_log 读已完成任务',
                  not err and (lg or {}).get('lines', 0) > 0
                  and any('fq2fa' in x for x in lg.get('log', [])),
                  f"lines={lg and lg.get('lines')}")
            err, sub2, text = await call(
                s, 'vp_run_tool',
                {'tool': 'convert',
                 'params': {'input': _make_fastq(), 'target': 'fasta'}})
            if not err and sub2:
                err2, cx, _ = await call(s, 'vp_cancel_task',
                                         {'task_id': sub2['result']['task_id']})
                ok_cancel = (not err2 and cx.get('accepted')) or \
                    '已结束' in str(cx)
                check('vp_cancel_task 受理或任务已完', bool(ok_cancel),
                      f"accepted={cx and cx.get('accepted')} "
                      f"status={cx and cx.get('status')}")
                csnap, cerr = await poll_task(
                    s, sub2['result']['task_id'], timeout=120)
                check('取消后到合法终态',
                      bool(csnap) and csnap.get('status') in
                      ('done', 'failed', 'cancelled'),
                      f"status={(csnap or {}).get('status')}")
                err3, _, _ = await call(s, 'vp_delete_tool_run',
                                        {'run': sub2['result']['run'],
                                         'confirm': True})

            print('== 3d. CLI 能力面 ==')
            err, ch, _ = await call(s, 'vp_cli_help', {})
            check('vp_cli_help 37 命令', not err and ch.get('count') == 37,
                  f"count={ch and ch.get('count')}")
            bhd = ch['commands']['build-host-db']
            req = [a['key'] for a in bhd['args'] if a['required']]
            check('build-host-db 必填 genome/taxid',
                  bhd['mode'] == 'task' and req == ['genome', 'taxid'])
            err, st, _ = await call(s, 'vp_cli_run', {'command': 'ictv-status'})
            check('vp_cli_run 真跑 ictv-status',
                  not err and st.get('status') == 'ok'
                  and st.get('exit_code') == 0,
                  (st or {}).get('summary', ''))
            err, _, text = await call(s, 'vp_cli_run',
                                      {'command': 'build-host-db',
                                       'args': {'genome': 'x', 'taxid': 1}})
            check('task 命令被 run 拒收',
                  err and 'vp_cli_submit' in text, text[:80])
            err, _, text = await call(s, 'vp_cli_run', {'command': 'nope'})
            check('未知 CLI 命令被拒', err and '未知 CLI' in text, text[:80])
            err, _, text = await call(s, 'vp_cli_submit',
                                      {'command': 'db-migrate',
                                       'args': {'to': 'D:\\nonexist_target'}})
            check('db-migrate 无 confirm 拒提', err and 'confirm' in text,
                  text[:80])
            err, sb, text = await call(
                s, 'vp_cli_submit',
                {'command': 'samples-backfill',
                 'args': {'dry_run': True}})
            check('vp_cli_submit 提交成功',
                  not err and sb.get('result', {}).get('task_id'),
                  text[:80])
            if not err and sb:
                ssnap, serr = await poll_task(s, sb['result']['task_id'],
                                              timeout=180)
                check('CLI 任务到终态 done', bool(ssnap) and
                      ssnap.get('status') == 'done',
                      serr or f"status={(ssnap or {}).get('status')} "
                              f"error={((ssnap or {}).get('error') or '')[:150]}")

                print('== 3e. 预筛模块 ==')
            err, prj, _ = await call(s, 'vp_prescreen_projects', {})
            check('vp_prescreen_projects 同步可查',
                  not err and isinstance(prj.get('projects'), list),
                  f"count={prj and prj.get('count')}")
            err, _, text = await call(
                s, 'vp_prescreen_run',
                {'name': 'bad name!', 'runs': ['SRR1']})
            check('预筛非法项目名被拒',
                  err and '项目名' in text, text[:80])
            err, ps, text = await call(
                s, 'vp_prescreen_run', {'name': 'mcp_ps_empty', 'runs': []})
            check('预筛空 runs 被拒', err and 'Run' in text, text[:80])
            err, pr, _ = await call(s, 'vp_prescreen_result',
                                    {'name': 'mcp_ps_empty'})
            check('预筛结果不存在时 exists=false',
                  not err and pr.get('exists') is False)

            if '--pipeline' in sys.argv[1:]:
                print('== 4. 流水线驱动（fq2fa 单阶段冒烟）==')
                err, psub, text = await call(
                    s, 'vp_pipeline_start',
                    {'sample': 'mcp_smoke', 'r1': fq,
                     'stages': ['fq2fa'], 'threads': 2})
                check('流水线提交成功',
                      not err and psub.get('result', {}).get('task_id'),
                      text[:120])
                if not err and psub:
                    psnap, perr = await poll_task(
                        s, psub['result']['task_id'], timeout=300)
                    check('流水线到终态 done', bool(psnap) and
                          psnap.get('status') == 'done',
                          perr or f"status={(psnap or {}).get('status')} "
                                  f"error={((psnap or {}).get('error') or '')[:200]}")
                    import shutil as _sh
                    _sh.rmtree(psub['result']['sample_dir'],
                               ignore_errors=True)

    print(f'\n通过 {len(PASS)} / 失败 {len(FAIL)}')
    if FAIL:
        print('失败项:', ', '.join(FAIL))
        return 1
    return 0


if __name__ == '__main__':
    if '--gate' in sys.argv[1:]:
        sys.exit(asyncio.run(test_gate()))
    elif '--http' in sys.argv[1:]:
        _a = sys.argv[1:]
        _i = _a.index('--http')
        _port = (int(_a[_i + 1])
                 if len(_a) > _i + 1 and _a[_i + 1].isdigit() else 8899)
        _tok = _a[_a.index('--token') + 1] if '--token' in _a else None
        sys.exit(asyncio.run(test_http(_port, _tok)))
    sys.exit(asyncio.run(main()))
