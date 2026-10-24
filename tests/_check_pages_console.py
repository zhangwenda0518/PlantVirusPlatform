# -*- coding: utf-8 -*-
"""全页面控制台巡检（Playwright）：每个页面真实加载，收集 JS 报错 / 页面异常 /
失败的请求，另验证关键全局函数存在。2026-09-13 app.js 拆分时引入，防拆分回归。

运行: python tests/_check_pages_console.py [--base http://127.0.0.1:8765]
不带 --base 时自起 Flask 服务（werkzeug 线程），跑完即关。

两种驱动方式：
  · 默认 —— 单进程内依次访问全部页面（最快）；
  · `--isolated` —— 每页起一个独立子进程（`--page <url>` 单页模式）。
    本机实测：同一 chromium 进程内连续访问会在第 3~6 页之间随机段错误
    （headless shell 与完整 chromium 均如此，单页独立访问则始终正常），
    属于累积状态问题。要一次跑完 21 页时用 `--isolated` 规避，代价是
    每页多一次服务启动开销。
"""
import os
import subprocess
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

BASE = None
if '--base' in sys.argv:
    BASE = sys.argv[sys.argv.index('--base') + 1]

# 单页模式：只巡检这一个 URL（供 --isolated 的父进程调用）
ONLY_PAGE = None
if '--page' in sys.argv:
    ONLY_PAGE = sys.argv[sys.argv.index('--page') + 1]

if not BASE and ('--page' in sys.argv or '--fast' in sys.argv):
    # 父进程（默认逐页独立进程模式）只负责派发子进程，服务由各子进程自己起
    # （串行，端口不冲突）。18766 被残留进程占用时自动换临时端口——
    # 否则子进程在这里崩掉，父进程只会看到「无结论输出」。
    from werkzeug.serving import make_server
    import threading
    import app as app_module
    try:
        _srv = make_server('127.0.0.1', 18766, app_module.app, threaded=True)
        _port = 18766
    except OSError:
        _srv = make_server('127.0.0.1', 0, app_module.app, threaded=True)
        _port = _srv.server_port
    threading.Thread(target=_srv.serve_forever, daemon=True).start()
    BASE = f'http://127.0.0.1:{_port}'
    import time
    time.sleep(1.0)


def pages_to_check():
    """全部无参数的 GET 页面路由（排除 API / 静态 / virome SPA）。"""
    import app as app_module
    out = []
    for rule in app_module.app.url_map.iter_rules():
        if 'GET' not in rule.methods:
            continue
        r = rule.rule
        if (r.startswith('/api/') or r.startswith('/static')
                or '<' in r or r.rstrip('/').endswith('/virome')
                or r == '/virome'):
            continue
        out.append(r)
    return sorted(set(out))


def main():
    from playwright.sync_api import sync_playwright
    pages = [ONLY_PAGE] if ONLY_PAGE else pages_to_check()
    if not ONLY_PAGE:
        print(f'待巡检页面 {len(pages)} 个')
    bad = []
    with sync_playwright() as pw:
        # 受限环境下 chromium 的沙箱与 /dev/shm 会引发随机段错误
        browser = pw.chromium.launch(
            headless=True,
            args=['--no-sandbox', '--disable-dev-shm-usage', '--disable-gpu'])
        ctx = browser.new_context(viewport={'width': 1440, 'height': 900})
        page = ctx.new_page()
        issues = []

        page.on('pageerror', lambda e: issues.append(('pageerror', str(e)[:200])))
        page.on('console', lambda m: issues.append(('console.' + m.type, m.text[:200]))
                if m.type == 'error' else None)
        page.on('response', lambda r: issues.append(('http ' + str(r.status), r.url[:120]))
                if r.status >= 400 else None)

        for p in pages:
            issues.clear()
            try:
                page.goto(BASE + p, wait_until='networkidle', timeout=20000)
            except Exception as e:
                # networkidle 对轮询页可能超时：退回 load
                try:
                    issues.clear()
                    page.goto(BASE + p, wait_until='load', timeout=20000)
                    page.wait_for_timeout(800)
                except Exception as e2:
                    # 截 400 字符：截 120 时 URL 恰好被切在 "http://127.0.0.1:"
                    # 处，看起来像「invalid URL 是端口丢了」，实际是别的错
                    bad.append((p, 'goto失败: ' + str(e2)[:400]))
                    continue
            page.wait_for_timeout(400)
            if issues:
                bad.append((p, ' | '.join(f'{k}:{v}' for k, v in issues[:4])))

        # 关键全局函数抽查（拆分后 app.js / app-compare.js 的函数仍为全局）
        checks = [
            ('/tools', 'jfetch', 'function'),
            ('/tools', 'rhRefresh', 'function'),
            # 比较基因组（⑤⑥ + ⑦MSA/⑧树/⑨SDT 工作区）内嵌在工具页
            ('/tools', 'loadGbCollections', 'function'),
            ('/tools', 'defaultCollName', 'function'),
            ('/results', 'pagerHtml', 'function'),
            ('/results', 'rhRefresh', 'function'),
            # 2026-10-11：/pipeline 归档，全局函数抽查换到 /samples（同一批 app.js 函数）
            ('/samples', 'selectSample', 'function'),
            ('/samples', 'pasteSeq', 'function'),
        ]
        if ONLY_PAGE:                      # 单页模式：全局函数抽查只在 /tools /results 上做
            checks = [c for c in checks if c[0] == ONLY_PAGE]
        for route, fn, want in checks:
            issues.clear()
            try:
                page.goto(BASE + route, wait_until='load', timeout=20000)
                page.wait_for_timeout(300)
            except Exception:
                pass
            got = page.evaluate(f'typeof {fn}')
            if got != want:
                bad.append((route, f'全局 {fn} = {got}（期望 {want}）'))

        browser.close()

    if ONLY_PAGE:
        # 单页模式：父进程只解析这一行（末行 = 结论）
        print(('OK ' if not bad else 'BAD ')
              + (bad[0][1] if bad else 'no issues'))
        return 0 if not bad else 1

    print('=' * 64)
    if bad:
        print(f'✘ {len(bad)} 个问题:')
        for p, msg in bad:
            print(f'  {p}\n      {msg}')
        sys.exit(1)
    print(f'✔ 全部 {len(pages)} 个页面加载无 JS 报错，全局函数抽查通过')


def run_isolated():
    """逐页独立进程巡检：规避同进程内连续访问导致的 chromium 随机段错误。"""
    try:
        import playwright  # noqa: F401
    except ImportError:
        # playwright 属开发依赖（requirements-dev.txt），未装时跳过而非失败，
        # 与其他浏览器探针的约定一致。
        print('[SKIP] 未安装 playwright，跳过真浏览器页面巡检')
        print('       安装: python -m pip install -r requirements-dev.txt')
        return 0
    pages = pages_to_check()
    print(f'逐页独立进程巡检，共 {len(pages)} 页（每页一个新进程）')
    me = os.path.abspath(__file__)
    bad = []
    import tempfile
    for i, p in enumerate(pages, 1):
        # 子进程输出必须落**文件**而不是管道（2026-10-24 修）：playwright
        # 拉起的 chrome 是孙进程，会继承管道句柄——子进程超时被 kill 后
        # 管道写端仍被 chrome 握着，subprocess.run 要等 EOF 就**永久挂死**
        # （实测批量回归卡在这里 6 分钟无输出）。文件句柄没有这个问题。
        log = os.path.join(tempfile.gettempdir(), f'_vp_page_{i:02d}.log')
        msg = ''
        rc = None
        with open(log, 'w', encoding='utf-8', errors='replace') as fh:
            proc = subprocess.Popen(
                [sys.executable, me, '--page', p],
                stdout=fh, stderr=fh, stdin=subprocess.DEVNULL)
            try:
                rc = proc.wait(timeout=180)
            except subprocess.TimeoutExpired:
                # taskkill /T 连 chrome 孙进程一起收掉，不留僵尸
                subprocess.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                               capture_output=True)
                try:
                    proc.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    proc.kill()
                rc, msg = proc.returncode, 'TIMEOUT（180s）'
        with open(log, encoding='utf-8', errors='replace') as fh:
            lines = [ln for ln in fh.read().splitlines() if ln.strip()]
        last = lines[-1] if lines else ''
        # 判据是子进程自报的末行结论（`OK ` / `BAD `），不是 returncode：
        # `--page` 模式曾经把 main() 的返回值丢掉（见文件末尾修注），
        # 于是无论页面坏成什么样 returncode 都是 0。现在两者都查。
        if msg == '' and last.startswith('OK ') and rc == 0:
            print(f'[{i:2d}/{len(pages)}] {p:28s} ✔', flush=True)
        else:
            if msg == '':
                msg = last[4:] if last.startswith('BAD ') else (
                    last or f'子进程无结论输出（rc={rc}）')
            print(f'[{i:2d}/{len(pages)}] {p:28s} ✘ {msg[:240]}', flush=True)
            bad.append((p, msg))
    print('=' * 64)
    if bad:
        print(f'✘ {len(bad)}/{len(pages)} 个页面有问题:')
        for p, msg in bad:
            print(f'  {p}\n      {msg}')
        return 1
    print(f'✔ 全部 {len(pages)} 个页面加载无 JS 报错，全局函数抽查通过')
    return 0


if __name__ == '__main__':
    # 默认走逐页独立进程：单进程内连续访问在本机 chromium 下会随机段错误
    # （第 3~6 页之间不等），会让「页面巡检」变成不可信的红灯。要快可显式
    # 加 --fast（单进程连续访问）。
    # 2026-09-16 修：单页模式必须 sys.exit(main()) 把结论传出去 —— 原来只写
    # `main()`，返回值被丢掉 ⇒ 子进程永远 exit 0 ⇒ 父进程的 returncode 判据
    # 恒真，页面全坏也报「全部页面 ✔」。
    if '--page' in sys.argv or '--fast' in sys.argv:
        sys.exit(main())
    else:
        sys.exit(run_isolated())
