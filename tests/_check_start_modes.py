# -*- coding: utf-8 -*-
"""启动模式检查：桌面窗口（gui）与网页（web）两条路径的选择与回退。

背景：原先只能用环境变量 `VP_GUI=browser` 切模式，双击 .bat 的人没有入口；
现在支持 `--gui` / `--web`（命令行优先于环境变量），并新增两个启动脚本。
模式选择是纯用户可见行为，且很容易在改动 main() 时被弄坏，故固化在此。

做法：在独立子进程里跑 app.main()，把副作用全换成桩——
  · _pick_port / _start_server：假端口与假 server
  · _open_window：与真函数同语义（web 模式立即返回 False）
  · _open_browser：只记录
  · threading.Thread：同步执行（消除"浏览器调用尚未发生"的竞态）
  · time.sleep：抛 KeyboardInterrupt，从浏览器模式的死循环里跳出
不真的开窗口、不真的开浏览器、不占端口。
"""
import json
import os
import subprocess
import sys

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fails = []

CHILD = r'''
import json, sys, types
sys.path.insert(0, %r)
import app

mode = app._gui_mode()
calls = []

def fake_window(url, m='gui'):
    calls.append(('window', m))
    return m != 'web'          # 与真实现同语义

app._open_window = fake_window
app._open_browser = lambda url: calls.append(('browser', url))
app._pick_port = lambda: 8799
# 测试子进程不与真实例互斥：撞上环境里遗留的平台实例锁会让
# main() 走 busy 退出分支，测试输出无从产生（2026-10-04 回归实测）
app._acquire_instance_lock = lambda: ('ok', None)
app._release_instance_lock = lambda p: None
app._start_server = lambda port: types.SimpleNamespace(shutdown=lambda: None)

class SyncThread:
    def __init__(self, target=None, args=(), **kw):
        self._t, self._a = target, args
    def start(self):
        self._t(*self._a)

app.threading = types.SimpleNamespace(Thread=SyncThread, get_ident=lambda: 0)

# app.time 只把 sleep 换成"抛 KeyboardInterrupt"（跳出浏览器模式死循环），
# 其余属性（strftime/time/monotonic…）透传真 time —— app.main() 的其它动作
# （如单实例锁写 pid+时间戳）也要用；只塞一个 sleep 的 SimpleNamespace 会
# 让 add 新功能时报 AttributeError（2026-10-03 实例锁 strftime 即此）。
import time as _real_time

class TimeStub:
    def __getattr__(self, name):
        return getattr(_real_time, name)
    def sleep(self, s):
        raise KeyboardInterrupt()

app.time = TimeStub()
try:
    app.main()
except KeyboardInterrupt:
    pass
print('@@' + json.dumps({'mode': mode, 'calls': calls}))
''' % ROOT


def run(args, env_extra):
    env = dict(os.environ)
    env.pop('VP_GUI', None)
    env.update(env_extra)
    out = subprocess.run([sys.executable, '-c', CHILD] + args, cwd=ROOT, env=env,
                         capture_output=True, text=True, encoding='utf-8')
    line = [l for l in (out.stdout or '').splitlines() if l.startswith('@@')]
    if not line:
        return {'mode': '?', 'calls': [],
                'raw': (out.stdout or '') + (out.stderr or '')}
    return json.loads(line[-1][2:])


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        fails.append(msg)


def kinds(r):
    return [c[0] for c in r['calls']]


def main():
    print('[默认] 无参数、无环境变量 → 桌面窗口，且不回退浏览器')
    r = run([], {})
    check(r['mode'] == 'gui', f"模式判定 = {r['mode']}")
    check(r['calls'] and r['calls'][0] == ['window', 'gui'],
          f"以 gui 模式调用 _open_window: {r['calls']}")
    check('browser' not in kinds(r), f"窗口成功时不回退浏览器: {kinds(r)}")

    print('\n[--web] 网页模式：不尝试开窗口，直接进浏览器分支')
    r = run(['--web'], {})
    check(r['mode'] == 'web', f"模式判定 = {r['mode']}")
    check(r['calls'] and r['calls'][0] == ['window', 'web'],
          f"传给 _open_window 的是 web（它会立即返回 False）: {r['calls']}")
    check('browser' in kinds(r), f"确实走浏览器分支: {kinds(r)}")

    print('\n[--gui] 强制桌面窗口')
    r = run(['--gui'], {})
    check(r['mode'] == 'gui' and r['calls'][0] == ['window', 'gui'],
          f"模式={r['mode']} 首次调用={r['calls'][:1]}")
    check('browser' not in kinds(r), '未回退浏览器')

    print('\n[VP_GUI=browser] 旧写法（环境变量）仍生效')
    r = run([], {'VP_GUI': 'browser'})
    check(r['mode'] == 'web' and 'browser' in kinds(r),
          f"模式={r['mode']} 调用={kinds(r)}")

    print('\n[优先级] 命令行 --gui 覆盖 VP_GUI=browser')
    r = run(['--gui'], {'VP_GUI': 'browser'})
    check(r['mode'] == 'gui' and 'browser' not in kinds(r),
          f"模式={r['mode']} 调用={kinds(r)}")

    print('\n[别名] --browser / --no-gui / --webview / VP_GUI=gui')
    check(run(['--browser'], {})['mode'] == 'web', '--browser → web')
    check(run(['--no-gui'], {})['mode'] == 'web', '--no-gui → web')
    check(run(['--webview'], {})['mode'] == 'gui', '--webview → gui')
    check(run([], {'VP_GUI': 'gui'})['mode'] == 'gui', 'VP_GUI=gui → gui')

    print('\n[启动脚本] 两个显式 .bat 存在且带对应对参数')
    for fn, flag in (('启动平台-桌面窗口.bat', '--gui'),
                     ('启动平台-网页.bat', '--web')):
        p = os.path.join(ROOT, fn)
        ok = os.path.isfile(p)
        body = ''
        if ok:
            import io
            body = io.open(p, encoding='utf-8', errors='replace').read()
        check(ok and f'app.py {flag}' in body,
              f'{fn} 存在且调用 app.py {flag}')

    print('\n' + '=' * 52)
    if fails:
        print(f'FAILED: {len(fails)} 项')
        for f in fails:
            print('  - ' + f)
        return 1
    print('启动模式检查通过 ✔')
    return 0


if __name__ == '__main__':
    sys.exit(main())
