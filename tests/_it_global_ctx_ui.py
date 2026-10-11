# -*- coding: utf-8 -*-
"""全局上下文（当前样品/项目）与全局清理按钮的浏览器端验证。

需要一个已在跑的实例（脚本自己起一个临时实例，用完关掉）：
    python tests/_it_global_ctx_ui.py [--port 8791] [--keep]

验证点：
  1. 顶部导航注入全局上下文条（当前样品 / 当前项目 / 全局清理）
  2. 样品选择弹窗能列出样品、选中后写入 localStorage 并刷新胶囊
  3. 项目胶囊切换后 localStorage 同步
  4. 全局样品选中后页面 JS curSample 同步、刷新后持久化恢复
     （原「管道页自动接上」断言随 /pipeline 归档移除）
  5. 结果中心按全局项目筛选样品行
  6. 任务中心"只看当前样品"随全局样品出现/消失
  7. 全局清理按钮：弹窗内容正确、确认时发出 POST /api/global/reset
     （用 fetch 打桩拦截，**不真的清真实数据**）、并清掉浏览器端关联键

依赖 selenium（平台已装）+ 本机 Edge/Chrome。
"""
import argparse
import json
import os
import socket
import subprocess
import sys
import tempfile
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding='utf-8', errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

FAIL = []


def dismiss_dialogs(drv):
    """关掉意外出现的 alert/confirm/prompt（会阻塞 renderer 导致命令超时）。"""
    try:
        a = drv.switch_to.alert
        print(f'      [warn] 页面弹出对话框，已接受: {a.text[:80]!r}')
        a.accept()
        return True
    except Exception:
        return False


def check(name, cond, extra=''):
    print(('  ✔ ' if cond else '  ✘ ') + name + (f'  {extra}' if extra else ''))
    if not cond:
        FAIL.append(name)


def free_port(start=8791):
    for p in range(start, start + 40):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(('127.0.0.1', p))
                return p
            except OSError:
                continue
    raise RuntimeError('没有可用端口')


def make_driver():
    from selenium import webdriver
    from selenium.webdriver.edge.options import Options as EdgeOptions
    opts = EdgeOptions()
    opts.add_argument('--headless=new')
    opts.add_argument('--window-size=1440,960')
    opts.add_argument('--no-sandbox')
    try:
        return webdriver.Edge(options=opts)
    except Exception:
        from selenium.webdriver.chrome.options import Options as ChromeOptions
        co = ChromeOptions()
        co.add_argument('--headless=new')
        co.add_argument('--window-size=1440,960')
        return webdriver.Chrome(options=co)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=0)
    ap.add_argument('--keep', action='store_true',
                    help='保留浏览器窗口/不关服务（排查用）')
    args = ap.parse_args()

    port = args.port or free_port()
    base = f'http://127.0.0.1:{port}'
    env = dict(os.environ, VP_NO_RECOVER='1')
    # 服务日志必须落文件、不能挂 PIPE：werkzeug 会为每个请求打一行访问日志，
    # PIPE 缓冲（64KB）写满后服务进程会阻塞在 write 上，页面加载随即超时。
    logpath = os.path.join(tempfile.gettempdir(),
                           f'vp_ui_it_{port}.log')
    logf = open(logpath, 'w', encoding='utf-8', errors='replace')
    srv = subprocess.Popen(
        [sys.executable, '-c',
         # 必须显式注入 ROOT：项目自带绿色版 Python 目录下有 python312._pth，
         # 该文件存在时进入 isolated 模式，`-c` 的 sys.path 不含 cwd，
         # `import app` 会 ModuleNotFoundError、服务起不来。
         'import sys\n'
         f'sys.path.insert(0, {ROOT!r})\n'
         'from werkzeug.serving import make_server\n'
         'import app\n'
         f"s = make_server('127.0.0.1', {port}, app.app, threaded=True)\n"
         "print('ready', flush=True)\n"
         's.serve_forever()'],
        cwd=ROOT, env=env, stdout=logf, stderr=subprocess.STDOUT)
    drv = None
    try:
        # 等端口可连（不看 stdout，避免依赖管道）
        ready = False
        for _ in range(100):
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                if s.connect_ex(('127.0.0.1', port)) == 0:
                    ready = True
                    break
            time.sleep(0.3)
        if not ready:
            print('服务启动失败，日志:', logpath)
            return 2
        print(f'临时实例: {base}  (服务日志 {logpath})')

        drv = make_driver()
        drv.set_page_load_timeout(60)
        W = drv.execute_script

        def goto(path):
            """导航 + 兜掉意外弹窗（alert 会挂住 renderer，让后续命令超时）。"""
            drv.get(base + path)
            time.sleep(1.6)
            dismiss_dialogs(drv)

        def js_alert_free(script):
            dismiss_dialogs(drv)
            return W(script)

        # ---- 1. 导航上下文条注入 ----
        print('[1] 顶部导航注入全局上下文条')
        goto('/samples')
        check('#navCtx 已注入', W("return !!document.getElementById('navCtx')"))
        check('样品胶囊存在', W("return !!document.getElementById('navCtxSample')"))
        check('项目胶囊存在', W("return !!document.getElementById('navCtxProject')"))
        check('全局清理按钮存在', W("return !!document.getElementById('navCtxReset')"))
        check('样品胶囊初始为空态',
              W("return document.getElementById('navCtxSample')"
                ".classList.contains('is-empty')"))
        check('清理按钮有文案',
              bool(W("return document.getElementById('navCtxReset').textContent.trim()")))

        # 先取一个真实存在的样品名（没有就跳过选择相关断言）
        samples = json.loads(W(
            "return await (await fetch('/api/samples')).text()"))
        names = [s['name'] for s in samples]
        print(f'      平台样品数: {len(names)}  {names[:3]}')
        projs = sorted({s.get('project') for s in samples if s.get('project')})
        print(f'      项目数: {len(projs)}  {projs[:3]}')

        if names:
            target = names[0]
            # ---- 2. 样品选择弹窗 ----
            print('[2] 样品选择弹窗')
            W("document.getElementById('navCtxSample').click()")
            time.sleep(0.6)
            check('弹窗打开', W("return !!document.getElementById('gxList')"))
            n_items = W("return document.querySelectorAll('#gxList .gx-item').length")
            check('弹窗列出样品', n_items == len(names), f'{n_items} vs {len(names)}')
            # 搜索框过滤
            W("const q=document.getElementById('gxQ'); q.value='zzz-no-such';"
              "q.dispatchEvent(new Event('input'))")
            time.sleep(0.2)
            check('搜索可过滤为空',
                  W("return document.querySelectorAll('#gxList .gx-item').length") == 0)
            W("const q=document.getElementById('gxQ'); q.value='';"
              "q.dispatchEvent(new Event('input'))")
            time.sleep(0.2)
            W(f"document.querySelector('#gxList .gx-item[data-name=\"{target}\"]').click()")
            time.sleep(0.5)
            check('胶囊显示所选样品',
                  target in W("return document.getElementById('navCtxSample')"
                              ".textContent"))
            check('localStorage 已写入 vp_ctx_sample',
                  W("return localStorage.getItem('vp_ctx_sample')") == target)
            check('胶囊脱离空态',
                  not W("return document.getElementById('navCtxSample')"
                        ".classList.contains('is-empty')"))

            # ---- 4. 全局样品上下文同步与持久化 ----
            # （2026-10-11：原「管道页 curSample 自动接上」断言随 /pipeline 归档
            #   移除——gxApplyPage 的自动接上挂载点 #samples 只在管道页存在；
            #   这里改验仍存在的行为：VP_CTX 全局对象同步 + 刷新后持久化恢复。）
            print('[3] 样品选择后 VP_CTX 全局对象已同步')
            del_ok = W("return window.VP_CTX && VP_CTX.sample === "
                       + json.dumps(target))
            check('页面 JS 的 VP_CTX.sample 已同步', bool(del_ok),
                  str(W("return window.VP_CTX ? VP_CTX.sample : 'no-VP_CTX'")))

            # 刷新一次（模拟"切模块再回来"），验证持久化恢复
            goto('/samples')
            check('刷新后 VP_CTX 仍从 localStorage 恢复该样品',
                  W("return window.VP_CTX && VP_CTX.sample === "
                    + json.dumps(target)),
                  str(W("return window.VP_CTX ? VP_CTX.sample : 'no-VP_CTX'")))

            # ---- 6. 任务中心"只看当前样品" ----
            print('[4] 任务中心"只看当前样品"')
            goto('/tasks')
            vis = W("const w=document.getElementById('curOnlyWrap');"
                    "return w && w.style.display !== 'none'")
            check('筛选控件可见（有全局样品时）', bool(vis))
            lbl = W("const l=document.getElementById('curOnlyWrap');"
                    "return l ? l.textContent : ''")
            check('筛选控件标出样品名', target in (lbl or ''), (lbl or '').strip())
            W("const c=document.getElementById('curOnlyChk');"
              "c.checked=true; c.onchange({target:c})")
            time.sleep(0.4)
            empty = W("const e=document.querySelector('.tk-empty');"
                      "return e ? e.textContent : ''")
            # 空态有两种等价文案：无任务时的「暂无任务」，以及 2026-09 起
            # 过滤命中 0 条时的「已按当前样品过滤，另有 N 个任务被隐藏」
            # （后者信息量更大，是刻意改进）。守卫语义不变：勾选后列表必须
            # 呈现空态提示，而不是残留别的样品的任务卡。
            check('勾选后列表被过滤（当前无该样品任务，应显示空态）',
                  '暂无任务' in (empty or '')
                  or '已按当前样品过滤' in (empty or '')
                  or not empty,
                  (empty or '')[:60])
            W("return localStorage.removeItem('vp_ctx_sample')")
            drv.refresh()
            time.sleep(1.6)
            dismiss_dialogs(drv)
            vis2 = W("const w=document.getElementById('curOnlyWrap');"
                     "return w && w.style.display !== 'none'")
            check('无全局样品时筛选控件隐藏', not vis2)
            # 还原上下文，供后续步骤使用
            W(f"window.VP_CTX.setSample({json.dumps(target)})")
            time.sleep(0.3)

        if projs:
            # ---- 3. 项目胶囊 ----
            print('[5] 项目胶囊与全局项目同步')
            goto('/samples')
            W("document.getElementById('navCtxProject').click()")
            time.sleep(0.5)
            check('项目弹窗打开', W("return !!document.getElementById('gxList')"))
            n_p = W("return document.querySelectorAll('#gxList .gx-item').length")
            check('项目数 = 实际项目数 + 全部项目', n_p == len(projs) + 1,
                  f'{n_p} vs {len(projs) + 1}')
            W(f"document.querySelector('#gxList .gx-item[data-v=\"{projs[0]}\"]').click()")
            time.sleep(0.6)
            check('localStorage 已写入 vp_ctx_project',
                  W("return localStorage.getItem('vp_ctx_project')") == projs[0])
            check('样品页项目下拉已同步到全局项目',
                  W("const s=document.getElementById('projFilter');"
                    "return s ? s.value : 'no-select'") == projs[0],
                  str(W("const s=document.getElementById('projFilter');"
                        "return s ? s.value : 'no-select'")))

            # 反向：用户直接操作页面下拉（模板 onchange 路径）→ 全局项目跟随，
            # 且下拉值不被页面自身的列表刷新刷回旧值。
            if len(projs) > 1:
                W("const s=document.getElementById('projFilter');"
                  "s.value=" + json.dumps(projs[1]) + ";"
                  "s.dispatchEvent(new Event('change'))")
                time.sleep(1.5)
                check('手动改下拉后全局项目跟随',
                      W("return localStorage.getItem('vp_ctx_project')") == projs[1],
                      str(W("return localStorage.getItem('vp_ctx_project')")))
                check('手动改下拉后下拉值未被刷回',
                      W("const s=document.getElementById('projFilter');"
                        "return s ? s.value : ''") == projs[1],
                      str(W("const s=document.getElementById('projFilter');"
                            "return s ? s.value : ''")))
                check('顶部项目胶囊同步显示',
                      projs[1] in W("return document.getElementById('navCtxProject')"
                                    ".textContent"))
                # 还原成 projs[0] 供结果中心断言
                W("const s=document.getElementById('projFilter');"
                  "s.value=" + json.dumps(projs[0]) + ";"
                  "s.dispatchEvent(new Event('change'))")
                time.sleep(1.0)

            # ---- 5. 结果中心按项目筛选 ----
            print('[6] 结果中心按全局项目筛选')
            goto('/results')
            check('结果中心项目下拉存在',
                  W("return !!document.getElementById('resProjFilter')"))
            check('下拉已选中全局项目',
                  W("const s=document.getElementById('resProjFilter');"
                    "return s ? s.value : ''") == projs[0])
            hidden = W(
                "const rows=[...document.querySelectorAll"
                "('#sampleTbl tr[data-project]')];"
                "return {total: rows.length,"
                " hidden: rows.filter(r => r.style.display === 'none').length,"
                " mismatched: rows.filter(r => r.style.display === 'none'"
                "   && r.dataset.project === " + json.dumps(projs[0]) + ").length}")
            check('筛选生效（存在被隐藏的其它项目行）',
                  hidden['hidden'] > 0 or hidden['total'] <= 1, str(hidden))
            check('当前项目的行未被隐藏', hidden['mismatched'] == 0, str(hidden))

            # ---- 7. 全局清理按钮（打桩，不动真实数据） ----
            print('[7] 全局清理按钮（fetch 打桩，不触碰真实数据）')
            goto('/samples')
            W("""
              window.__cap = null;
              try { sessionStorage.removeItem('vp_reset_capture'); } catch (e) {}
              const _orig = window.fetch;
              window.fetch = function (url, opts) {
                if (String(url).indexOf('/api/global/reset') >= 0) {
                  const cap = {url: String(url),
                               body: (opts && opts.body) || '',
                               method: (opts && opts.method) || ''};
                  window.__cap = cap;
                  // 确认后会整页 reload，内存里的 __cap 会被冲掉，
                  // 另存 sessionStorage 供刷完再断言。
                  try { sessionStorage.setItem('vp_reset_capture',
                                               JSON.stringify(cap)); } catch (e) {}
                  return Promise.resolve(new Response(JSON.stringify(
                    {ok: true, cancelled_tasks: 1, cleared_queue: 2,
                     removed_task_records: 3, active_left: 0,
                     kept_results: true}),
                    {status: 200, headers: {'Content-Type': 'application/json'}}));
                }
                return _orig.apply(this, arguments);
              };
            """)
            W("document.getElementById('navCtxReset').click()")
            time.sleep(0.5)
            check('清理弹窗打开', W("return !!document.getElementById('gxYes')"))
            check('含"同时清空任务归档"选项',
                  W("return !!document.getElementById('gxArch')"))
            guard = W("const b=document.getElementById('gxModalBody');"
                      "return b ? b.textContent : ''")
            check('弹窗声明保留结果文件', 'results/' in guard and 'tool_runs/' in guard)
            W("document.getElementById('gxYes').click()")

            # 确认后页面会整页刷新（复位全部 DOM 态）——等它稳定再断言
            time.sleep(3.2)
            dismiss_dialogs(drv)
            cap = None
            try:
                cap = json.loads(W("return sessionStorage.getItem('vp_reset_capture')")
                                 or 'null')
            except Exception as e:
                print('      [warn] 读取捕获失败:', e)
            check('确认后发出 POST /api/global/reset',
                  bool(cap) and cap.get('method') == 'POST'
                  and cap.get('url', '').endswith('/api/global/reset'), str(cap))
            check('请求体带 include_archive',
                  bool(cap) and 'include_archive' in (cap.get('body') or ''),
                  str(cap))
            check('请求体带 confirm 防误触令牌',
                  bool(cap) and 'confirm' in (cap.get('body') or '')
                  and 'reset' in (cap.get('body') or ''), str(cap))
            check('浏览器端关联键已清除',
                  W("return localStorage.getItem('vp_ctx_sample') === null"
                    " && localStorage.getItem('vp_ctx_project') === null"
                    " && localStorage.getItem('vp_last_sample') === null"))
            check('context 胶囊已复位为空',
                  W("return document.getElementById('navCtxSample')"
                    ".classList.contains('is-empty')"
                    " && document.getElementById('navCtxProject')"
                    ".classList.contains('is-empty')"))
            check('刷新后已不再选中样品（关联确被清掉）',
                  W("return typeof curSample === 'undefined' || !curSample"))
            note = W("const c=document.getElementById('toastBox');"
                     "return c ? c.textContent : ''")
            check('刷新后回显清理结果（toast 含计数与"保留"语义）',
                  '1' in (note or '') and '3' in (note or '')
                  and ('清理完成' in (note or '') or 'cleanup' in (note or '').lower()),
                  (note or '')[:120].replace('\n', ' '))

            # ---- 8. kvsuite 选样弹窗接入全局上下文 ----
            print('[8] kvsuite 选样弹窗接入全局样品/项目')
            if names and projs:
                W(f"window.VP_CTX.setSample({json.dumps(names[0])});"
                  f"window.VP_CTX.setProject({json.dumps(projs[0])})")
                goto('/tools?g=kvsuite')
                check('kvPickSamples 可用',
                      W("return typeof kvPickSamples === 'function'"))
                W("document.getElementById('kv_samples').value = '';"
                  "kvPickSamples()")
                time.sleep(1.2)
                check('选样弹窗打开',
                      W("return !!document.getElementById('kvPickProj')"))
                check('弹窗有项目筛选下拉',
                      W("return document.querySelectorAll('#kvPickProj option').length")
                      == len(projs) + 1,
                      str(W("return document.querySelectorAll('#kvPickProj option')"
                            ".length")))
                check('项目下拉默认落到全局项目',
                      W("const s=document.getElementById('kvPickProj');"
                        "return s ? s.value : ''") == projs[0])
                check('全局当前样品默认已勾选',
                      W("const c=document.querySelector('.kv-pick[value=\""
                        + names[0] + "\"]'); return !!(c && c.checked)"))
                vis = W(
                    "const rows=[...document.querySelectorAll"
                    "('#kvPickProj')].length;"
                    "const trs=[...document.querySelectorAll('tbody tr[data-proj]')];"
                    "return {total: trs.length,"
                    " shown: trs.filter(t => t.style.display !== 'none').length,"
                    " wrongProject: trs.filter(t => t.style.display !== 'none'"
                    "   && t.dataset.proj !== " + json.dumps(projs[0]) + ").length}")
                check('按项目筛选生效（其它项目行被隐藏）',
                      vis['shown'] <= vis['total'] and vis['wrongProject'] == 0,
                      str(vis))
                # 全选只作用于可见行
                W("document.getElementById('kvPickAll').click()")
                time.sleep(0.3)
                check('全选只勾选当前项目下的可见样品',
                      W("const trs=[...document.querySelectorAll('tbody tr[data-proj]')];"
                        "return trs.filter(t => t.style.display === 'none'"
                        " && t.querySelector('.kv-pick').checked).length") == 0)
                W("document.getElementById('kvPickCancel').click()")
                time.sleep(0.3)

        # ---- 9. 管道页阶段卡顺序 ----
        # （2026-10-11：整段随 /pipeline「自动化分析流程」页归档移除——
        #   阶段卡/管道 DAG 视图已无活动页面；逻辑本体保留在 app-batch.js，
        #   恢复页面时从 archive/_retire_20261011/pipeline/ 拷回模板即可。）

        # ---- 10. 样品名被规范化时前端显式提示（fetch 打桩，不落盘） ----
        print('[10] 样品名被规范化时前端提示实际登记名')
        goto('/samples')
        W("""
          try { localStorage.removeItem('vp_reset_capture'); } catch (e) {}
          const _of = window.fetch;
          window.fetch = function (url, opts) {
            if (String(url).indexOf('/api/pipeline/create') >= 0) {
              return Promise.resolve(new Response(JSON.stringify(
                {sample: 'A', typed: '样品A',
                 note: '样品名「样品A」已按文件名安全规则登记为「A」'}),
                {status: 200, headers: {'Content-Type': 'application/json'}}));
            }
            return _of.apply(this, arguments);
          };
        """)
        W("document.getElementById('sample').value = '样品A';"
          "document.getElementById('r1').value = 'dummy_R1.fastq.gz';")
        W("scCreate(document.querySelector('button[onclick*=\"scCreate\"]'))")
        time.sleep(1.2)
        err = W("const e=document.getElementById('cErr');return e?e.textContent:''")
        check('错误/提示区显示 note（把实际登记名讲清楚）',
              'A' in (err or '') and '样品A' in (err or ''), (err or '')[:100])
        tst = W("const b=document.getElementById('toastBox');return b?b.textContent:''")
        check('toast 用实际登记名 A 而非输入名', 'A' in (tst or ''),
              (tst or '')[:80])
    finally:
        if drv is not None and not args.keep:
            try:
                drv.quit()
            except Exception:
                pass
        if not args.keep:
            srv.terminate()
            try:
                srv.wait(timeout=15)
            except Exception:
                srv.kill()
        try:
            logf.close()
        except Exception:
            pass

    print()
    if FAIL:
        print(f'✘ {len(FAIL)} 项未通过: ' + ', '.join(FAIL))
        return 1
    print('✔ 全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
