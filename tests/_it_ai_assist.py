# -*- coding: utf-8 -*-
"""AI 助手（对话驱动代理）集成测试：Flask test_client + 本地模拟端点。

不依赖任何真实 AI 服务：起一个线程内的模拟 OpenAI 兼容 SSE 服务器，
第一轮返回工具调用、第二轮流式返回正文，完整覆盖「上下文注入 → 工具
循环 → 流式下发」链路。

用法: python tests/_it_ai_assist.py
退出码: 0 = 全部通过；1 = 存在失败。
"""
import json
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

for _s in (sys.stdout, sys.stderr):
    try:
        _stream = _s
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault('VP_NO_RECOVER', '1')
# 本守卫会发布大量事件 / 建执行单；关掉审计落盘，避免往 run/ai_audit/
# 写测试垃圾（审计本身的正确性由 tests/_check_ai_audit.py 在临时目录专测）。
os.environ.setdefault('VP_AI_AUDIT', '0')

FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


# ------------------------------------------------------------------
# 模拟 OpenAI 兼容端点：/v1/chat/completions
#   - 请求消息里没有 role=tool → 回 tool_calls（get_sample_digest）
#   - 有 role=tool        → 流式吐正文分片
# ------------------------------------------------------------------
class _MockAI(BaseHTTPRequestHandler):
    def log_message(self, *a):        # 静默访问日志
        pass

    def do_POST(self):
        n = int(self.headers.get('Content-Length') or 0)
        body = json.loads(self.rfile.read(n) or b'{}')
        msgs = body.get('messages') or []
        has_tool_result = any(m.get('role') == 'tool' for m in msgs)
        self.send_response(200)
        self.send_header('Content-Type', 'text/event-stream')
        self.end_headers()
        chunks = []
        if not has_tool_result:
            chunks.append({'choices': [{'delta': {'tool_calls': [
                {'index': 0, 'id': 'call_1', 'type': 'function',
                 'function': {'name': 'get_sample_digest', 'arguments': ''}}]}}]})
            chunks.append({'choices': [{'delta': {'tool_calls': [
                {'index': 0, 'function': {'arguments': '{"sample":'}}]}}]})
            chunks.append({'choices': [{'delta': {'tool_calls': [
                {'index': 0, 'function': {'arguments': ' "smoke_e2e"}'}}]}}]})
            chunks.append({'choices': [{'delta': {},
                                        'finish_reason': 'tool_calls'}]})
        else:
            for piece in ('模拟回答：', '检出 3 种病毒。'):
                chunks.append({'choices': [{'delta': {'content': piece}}]})
            chunks.append({'choices': [{'delta': {},
                                        'finish_reason': 'stop'}]})
        for c in chunks:
            self.wfile.write(f"data: {json.dumps(c)}\n\n".encode())
        self.wfile.write(b'data: [DONE]\n\n')
        self.wfile.flush()


def _sse_events(resp):
    """把测试客户端的流式响应解析成事件列表。"""
    out = []
    for frame in resp.get_data(as_text=True).split('\n\n'):
        frame = frame.strip()
        if frame.startswith('data:'):
            try:
                out.append(json.loads(frame[5:].strip()))
            except ValueError:
                pass
    return out


def main():
    from app import app as flask_app
    c = flask_app.test_client()

    # 原配置留档，结束时恢复
    orig_ai = dict(flask_app.config.get('cfg').ai) if False else None
    from Virus_Platform_Core.web.state import cfg
    orig_ai = dict(cfg.ai)

    print('== 1) 配置接口 ==')
    r = c.get('/api/ai/config')
    d = r.get_json()
    check('GET /api/ai/config 200', r.status_code == 200)
    check('presets 存在', len(d.get('presets') or []) >= 5)
    # 掩码格式：''（未配）或 '••••+尾4位'（本体绝不出后端）。
    # 本机可能已配置真实 Key（run/ai_local.json），这里做环境无关断言。
    mk = d['config'].get('api_key') or ''
    check('api_key 掩码格式（本体不出后端）',
          (not mk.strip()) or (mk.startswith('••••') and len(mk) <= 8
                               and not mk[4:].startswith('sk-')), mk)

    r = c.post('/api/ai/config', json={})
    check('空体 POST config → 400', r.status_code == 400)
    r = c.post('/api/ai/test', json={})
    check('空体 POST test → 400（无网络请求）', r.status_code == 400)
    r = c.post('/api/ai/chat', json={})
    check('空体 POST chat → 400', r.status_code == 400)

    print('== 2) 未配置时的降级 ==')
    cfg.ai = {'enabled': False, 'provider': 'custom', 'base_url': '',
              'api_key': '', 'model': '', 'temperature': 0.3,
              'max_tokens': 2048}
    r = c.post('/api/ai/chat', json={'messages': [{'role': 'user',
                                                   'content': '你好'}]})
    check('未配置 chat 仍 200（SSE 引导）', r.status_code == 200)
    evs = _sse_events(r)
    check('事件含 error 引导配置',
          any(e.get('type') == 'error' and '配置' in (e.get('message') or '')
              for e in evs), str(evs))
    check('事件含 done', any(e.get('type') == 'done' for e in evs))

    r = c.get('/api/ai/context')
    d = r.get_json()
    check('context 含样品清单', isinstance(d.get('samples'), list))
    check('context ready=False', d.get('ready') is False)

    print('== 3) 模拟端点全链路（工具循环 + 流式） ==')
    srv = ThreadingHTTPServer(('127.0.0.1', 0), _MockAI)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    port = srv.server_address[1]
    try:
        cfg.ai = {'enabled': True, 'provider': 'custom',
                  'base_url': f'http://127.0.0.1:{port}/v1',
                  'api_key': 'test-key-123456', 'model': 'mock-1',
                  'temperature': 0.2, 'max_tokens': 512}
        r = c.post('/api/ai/chat', json={
            'messages': [{'role': 'user', 'content': '解读样品'}],
            'sample': 'smoke_e2e'})
        check('chat 200', r.status_code == 200)
        evs = _sse_events(r)
        types = [e.get('type') for e in evs]
        check('出现 tool 事件', 'tool' in types, str(types))
        check('工具事件名为 get_sample_digest',
              any(e.get('type') == 'tool' and
                  e.get('name') == 'get_sample_digest' for e in evs))
        body = ''.join(e.get('content') or '' for e in evs
                       if e.get('type') == 'delta')
        check('正文增量拼接完整', '模拟回答：检出 3 种病毒。' == body, body)
        check('done 收尾', types[-1] == 'done', str(types))

        # 工具直调：路径穿越必须被拦
        from Virus_Platform_Core.ai_agent import digest, tools as ai_tools
        bad = ai_tools.dispatch('read_result_table',
                                {'sample': '..', 'relpath': 'platform.json'})
        check('read_result_table 拦越界', '"error"' in bad, bad)
        ok = ai_tools.dispatch('read_result_table', {
            'sample': 'smoke_e2e',
            'relpath': '02b_kvsuite/identify/all_viruses.summary.tsv'})
        check('read_result_table 正常读表', 'header' in ok, ok[:120])

        # 2026-10-08 新增通道：工具箱任务目录（run/tool_runs/<run>/）。
        # 与 read_result_table 共用 _safe_read_table_under，但根白名单不同，
        # 所以越界与正常读两条都要各自验一遍（防"只修一处"）。
        bad2 = ai_tools.dispatch('read_tool_run_file',
                                 {'run': '..', 'relpath': 'platform.json'})
        check('read_tool_run_file 拦越界', '"error"' in bad2, bad2[:120])
        bad3 = ai_tools.dispatch('read_tool_run_file',
                                 {'run': 'nosuchrun', 'relpath': 'x.tsv'})
        check('read_tool_run_file 拦不存在任务', '"error"' in bad3, bad3[:120])
        runs = digest.tool_runs_brief(limit=1)
        if runs:
            good = ai_tools.dispatch('read_tool_run_file',
                                     {'run': runs[0]['run'],
                                      'relpath': 'run.log'})
            check('read_tool_run_file 正常读任务文件', '"error"' not in good,
                  good[:120])
        else:
            print('  [SKIP] run/tool_runs 下暂无任务，跳过正常读分支')

        # 2026-10-08 新增通道：预筛产物目录（run/meta_search/<项目>/prescreen/）。
        # 预筛**不写** tool_runs（web/prescreen.py 只写 DIRS['meta_search']），
        # 故必须另立根白名单；越界 / 不存在 / 正常读三条各自验一遍。
        bad4 = ai_tools.dispatch('read_prescreen_file',
                                 {'project': '..', 'relpath': 'platform.json'})
        check('read_prescreen_file 拦非法项目名', '"error"' in bad4, bad4[:120])
        bad5 = ai_tools.dispatch('read_prescreen_file',
                                 {'project': 'nosuchproj',
                                  'relpath': 'stat_run_summary.tsv'})
        check('read_prescreen_file 拦不存在项目', '"error"' in bad5, bad5[:120])
        # 遍历真实项目：brief 的第一个项目不一定带 stat_run_summary.tsv
        # （残目录只有别的文件时此前会误报失败），找到有该表的为止
        pres = digest.prescreen_brief(limit=10)
        good2, used = None, None
        for _p in pres:
            t = ai_tools.dispatch('read_prescreen_file', {
                'project': _p['project'],
                'relpath': 'stat_run_summary.tsv'})
            if '"error"' not in t:
                good2, used = t, _p['project']
                break
        if good2 is not None:
            check(f'read_prescreen_file 正常读预筛表（{used}）',
                  '"error"' not in good2, good2[:120])
        else:
            print('  [SKIP] run/meta_search 下暂无预筛项目，跳过正常读分支')

        # digest 覆盖面（2026-10-08 补齐）：主链 4 段 + 扩展分析知情层。
        # 用**最小夹具**做功能性验证——只断言"键存在"是弱断言，无法区分
        # "代码在但读错文件"；这里造一个真样品目录让 digest 去读。
        import json as _json
        import shutil
        import tempfile
        from Virus_Platform_Core.config import DIRS as _DIRS
        tmp = tempfile.mkdtemp(prefix='vp_ai_digest_')
        try:
            sd = os.path.join(tmp, 'probe_sample')
            for sub in ('00_prep', '01_host_removal', '09_genome_plots'):
                os.makedirs(os.path.join(sd, sub))
            with open(os.path.join(sd, '01_host_removal', 'stats.json'),
                      'w', encoding='utf-8') as f:
                _json.dump({'ok': True, 'host_ratio': 0.9871,
                            'total_pairs': 1000000, 'kept_pairs': 12900,
                            'dropped_pairs': 987100,
                            'host_db': {'name': 'tomato', 'taxid': 4081},
                            'classify_params': {'confidence': 0.1}}, f)
            with open(os.path.join(sd, '00_prep', 'subsample.json'),
                      'w', encoding='utf-8') as f:
                _json.dump({'n_pairs': 200000}, f)
            with open(os.path.join(sd, '09_genome_plots', 'summary.json'),
                      'w', encoding='utf-8') as f:
                _json.dump({'n_seqs': 3, 'plots': ['a.svg', 'b.svg']}, f)
            _orig = digest.DIRS
            digest.DIRS = dict(_DIRS, results=tmp)
            try:
                det = (digest.sample_digest('probe_sample') or {}).get(
                    'detail') or {}
            finally:
                digest.DIRS = _orig
            check('digest 读到 ①宿主去除 host_ratio',
                  (det.get('host') or {}).get('host_ratio') == 0.9871,
                  str(det.get('host')))
            check('digest 读到 子采样 n_pairs',
                  (det.get('subsample') or {}).get('n_pairs') == 200000,
                  str(det.get('subsample')))
            check('digest 读到 ⑨基因组图 n_plots',
                  (det.get('gbdraw') or {}).get('n_plots') == 2,
                  str(det.get('gbdraw')))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        pb = digest.platform_brief()
        check('platform_brief 含 tool_runs/extensions/prescreen',
              all(k in pb for k in ('tool_runs', 'extensions', 'prescreen')),
              str(sorted(pb.keys())))

        # 扩展分析「稳定契约」通道（2026-10-08）：用最小夹具验证 digest 真的
        # 按契约读到数字，而不是只看到目录存在。四个契约工具各造一条。
        tmp2 = tempfile.mkdtemp(prefix='vp_ai_ext_')
        try:
            tr = os.path.join(tmp2, 'tool_runs')
            ms = os.path.join(tmp2, 'meta_search')
            def _w(rel, obj):
                p = os.path.join(tmp2, rel)
                os.makedirs(os.path.dirname(p), exist_ok=True)
                with open(p, 'w', encoding='utf-8') as f:
                    _json.dump(obj, f)
            _w('tool_runs/phylodyn_20261008_120000_abcdef/phylodyn_summary.json',
               {'clock': {'timetree_source': 'TreeTime', 'root_year': 1979.5,
                          'n_methods_ok': 1, 'rtt_fit': {'r2': 0.44}},
                'phylogeo': {'n_changes_ml': 4, 'n_corridors': 3},
                'warnings': ['w1', 'w2'], 'errors': []})
            _w('tool_runs/rtt_20261008_120000_abcdef/rtt/summary.json',
               {'slope': 0.008, 'r2': 0.619, 'n': 15, 'method': 'nj',
                'strong_signal': True, 'warnings': [],
                'aln': 'D:\\x\\example_phylogeo.fasta'})
            srs = os.path.join(ms, 'probe_proj', 'prescreen',
                               'stat_run_summary.tsv')
            os.makedirs(os.path.dirname(srs), exist_ok=True)
            with open(srs, 'w', encoding='utf-8', newline='\n') as f:
                f.write('acc\tstatus\tpass\n')
                f.write('R1\tOK\tPASS\n')
                f.write('R2\tOK\t\n')
                f.write('R3\tNO_ANALYSIS\t\n')
            _orig2 = digest.DIRS
            digest.DIRS = dict(_orig2, tool_runs=tr, meta_search=ms)
            try:
                rows = {r['tool']: r for r in digest.tool_runs_brief(limit=10)}
                ps = digest.prescreen_brief(limit=5)
            finally:
                digest.DIRS = _orig2
            pd_ = rows.get('phylodyn') or {}
            check('扩展契约：phylodyn 读到 root_year/n_corridors',
                  (pd_.get('metrics') or {}).get('root_year') == 1979.5
                  and (pd_.get('metrics') or {}).get('n_corridors') == 3,
                  str(pd_))
            check('扩展契约：phylodyn 派生 n_warn=2',
                  (pd_.get('metrics') or {}).get('n_warn') == 2, str(pd_))
            rt = rows.get('rtt') or {}
            check('扩展契约：rtt 读到 r2 与输入文件名',
                  (rt.get('metrics') or {}).get('r2') == 0.619
                  and rt.get('input') == 'example_phylogeo.fasta', str(rt))
            check('预筛计数：n_runs=3 / n_pass=1 / n_ok=2',
                  bool(ps) and ps[0].get('n_runs') == 3
                  and ps[0].get('n_pass') == 1 and ps[0].get('n_ok') == 2,
                  str(ps))
        finally:
            shutil.rmtree(tmp2, ignore_errors=True)

        print('== 4) 配置保存往返 ==')
        r = c.post('/api/ai/config', json={'enabled': True,
                                           'provider': 'deepseek',
                                           'temperature': 9})
        d = r.get_json()
        check('保存 200', r.status_code == 200)
        check('温度越界被钳到 2.0',
              d['config']['temperature'] == 2.0, str(d['config']))
        check('掩码含尾 4 位', d['config']['api_key'].endswith('3456'),
              d['config']['api_key'])
        r2 = c.get('/api/ai/config')
        check('GET 掩码一致', r2.get_json()['config']['api_key']
              == d['config']['api_key'])

        # 密钥分离铁律：platform.json（进 git）永不落 Key，
        # 真身在 run/ai_local.json（gitignore 覆盖）
        import json as _json
        from Virus_Platform_Core.config import PLATFORM_ROOT, _ai_local_path
        with open(os.path.join(PLATFORM_ROOT, 'platform.json'),
                  encoding='utf-8') as f:
            disk_ai = (_json.load(f).get('ai') or {})
        check('platform.json 落盘无 api_key',
              not str(disk_ai.get('api_key') or '').strip(), str(disk_ai))
        with open(_ai_local_path(), encoding='utf-8') as f:
            local_key = (_json.load(f) or {}).get('api_key') or ''
        check('run/ai_local.json 持有 Key 真身',
              local_key == 'test-key-123456')
        # 冷加载合并：新 Config 实例（模拟重启）能从密钥文件读回
        from Virus_Platform_Core.config import Config
        c2 = Config()
        check('重启后 Key 从密钥文件合并', c2.ai.get('api_key')
              == 'test-key-123456')
    finally:
        srv.shutdown()
        cfg.ai = orig_ai
        cfg.save()

    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
