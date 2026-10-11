# -*- coding: utf-8 -*-
"""顶部导航宽度自适应验证（真浏览器实测，不是看 CSS 猜）。

背景（2026-09-11 用户截图报障）：
`.topnav` 原先 `height: 54px` 固定高度，而 `.navlinks` 是 `flex-wrap: wrap`。
一级导航当时 13 个链接（约占 1150px，现按 NAV_GROUPS 派生），窗口一窄
`.navlinks` 内部就换行，但绿色横条仍只有 54px —— 换行出来的第二行
（样品/项目/全局清理/EN/📂）直接溢出到横条外的白底上。`height: auto` 的兜底
只在 ≤720px 生效，远低于用户窗口宽度，所以平时必然踩到。

⚠️ 新增一级导航项会加宽这条横条（2026-09-18 加 detect 组即 +1 项），
宽屏"必须单行"的断言（1920px 那条）是随之最可能先红的项 —— 加组后务必跑本脚本。

本脚本在多个窗口宽度下断言：
  1. 横条不溢出：`scrollHeight <= clientHeight + 1`
  2. 横条内每个子元素都完整落在横条矩形内（第二行没有跑到条外）
  3. 页面没有横向溢出：`documentElement.scrollWidth <= innerWidth + 1`
  4. 上下文胶囊 / 语言按钮 / 📂 都仍可见可点（自适应不能靠"藏掉功能"实现）
  5. 窄窗口下允许换行（row 数 ≥ 1），但必须被横条包住

用法：python tests/_it_nav_responsive.py [--keep]
依赖 selenium + 本机 Edge/Chrome。
"""
import argparse
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

# 一级导航链接数 = NAV_GROUPS 每个分组一条（.nav-drop 里的 <a>）
# + 4 条非分组项（/、/tasks、/settings、/help）。
# 别写死数字：09-17 恢复 compare 组、09-18 新增 detect 组、10-11 归档
# /pipeline，三次都会让它陈旧；而"数字对不上"这种红是噪声红，会掩盖真的溢出问题。
from Virus_Platform_Core.web.pages import NAV_GROUPS as _NAV_GROUPS  # noqa: E402

NON_GROUP_LINKS = 4
EXPECT_LINKS = len(_NAV_GROUPS) + NON_GROUP_LINKS

# 覆盖：宽屏单行 → 逐步变窄到必须换行 → 很窄
WIDTHS = [1920, 1680, 1560, 1440, 1366, 1280, 1180, 1024, 960, 820, 720, 640]
PAGES = ['/samples', '/tools?g=kvsuite', '/results']

FAIL = []


def check(name, cond, extra=''):
    print(('  ✔ ' if cond else '  ✘ ') + name + (f'  {extra}' if extra else ''))
    if not cond:
        FAIL.append(name)


def free_port(start=8791):
    for p in range(start, start + 60):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(('127.0.0.1', p))
                return p
            except OSError:
                continue
    raise RuntimeError('无可用端口')


def make_driver():
    from selenium import webdriver
    from selenium.webdriver.edge.options import Options as EdgeOptions
    o = EdgeOptions()
    o.add_argument('--headless=new')
    o.add_argument('--no-sandbox')
    o.add_argument('--window-size=1920,1000')
    try:
        return webdriver.Edge(options=o)
    except Exception:
        from selenium.webdriver.chrome.options import Options as ChromeOptions
        c = ChromeOptions()
        c.add_argument('--headless=new')
        c.add_argument('--window-size=1920,1000')
        return webdriver.Chrome(options=c)


MEASURE = """
const nav = document.querySelector('.topnav');
if (!nav) return {error: 'no .topnav'};
const nr = nav.getBoundingClientRect();
const linksBox = document.querySelector('.navlinks');
// 视觉行数不能用 round(top) 去重：.topnav 是 align-items:center，同一行里
// 不同高度的元素 top 本来就不同（实测把单行数成 5 行）。改为按"中心 Y 分簇"：
// 中心 Y 相差 > 10px 才算新的一行。
function countRows(rects) {
  const cs = rects.map(r => (r.top + r.bottom) / 2).sort((a, b) => a - b);
  let rows = 0, last = -1e9;
  for (const c of cs) { if (c - last > 10) { rows++; last = c; } }
  return rows;
}
const kids = [...nav.children].flatMap(c =>
  c.classList.contains('navlinks') ? [...c.children] : [c]);
const vis = kids.filter(k => k.getClientRects().length);
const out = vis.filter(k => {
  const r = k.getBoundingClientRect();
  return r.bottom > nr.bottom + 1 || r.top < nr.top - 1 ||
         r.right > nr.right + 1 || r.left < nr.left - 1;
}).map(k => (k.id || k.className || k.tagName) + '');
// 一级导航链接 = NAV_GROUPS 分组 + 非分组项（本脚本的 EXPECT_LINKS 与之同源）
const lv1 = [...document.querySelectorAll(
    '.navlinks > a, .navlinks > .nav-drop > a')]
  .filter(a => a.getClientRects().length);
const r = el => { if (!el || !el.getClientRects().length) return null;
                  const b = el.getBoundingClientRect();
                  return {w: Math.round(b.width), h: Math.round(b.height)}; };
const chip = document.getElementById('navCtxSample');
const reset = document.getElementById('navCtxReset');
const lang = document.getElementById('langBtn');
const icon = document.querySelector('.navlinks .iconbtn');
return {
  overflowY: nav.scrollHeight - nav.clientHeight,
  overflowX: document.documentElement.scrollWidth - window.innerWidth,
  outside: out,
  navH: Math.round(nr.height),
  navRows: countRows(vis.map(k => k.getBoundingClientRect())),
  linkRows: countRows(lv1.map(a => a.getBoundingClientRect())),
  links: lv1.length,
  linksBoxW: linksBox ? Math.round(linksBox.getBoundingClientRect().width) : 0,
  chip: r(chip), reset: r(reset), lang: r(lang), icon: r(icon),
  vw: window.innerWidth,
};
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--keep', action='store_true')
    args = ap.parse_args()

    port = free_port()
    base = f'http://127.0.0.1:{port}'
    env = dict(os.environ, VP_NO_RECOVER='1')
    # 服务日志必须落文件：挂 PIPE 时 werkzeug 的访问日志写满 64KB 后服务会
    # 阻塞在 write 上，页面加载随之超时。
    logpath = os.path.join(tempfile.gettempdir(), f'vp_nav_it_{port}.log')
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
        ok = False
        for _ in range(100):
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                if s.connect_ex(('127.0.0.1', port)) == 0:
                    ok = True
                    break
            time.sleep(0.3)
        if not ok:
            print('服务启动失败:', logpath)
            return 2
        print(f'临时实例 {base}（日志 {logpath}）')

        drv = make_driver()
        drv.set_page_load_timeout(60)
        W = drv.execute_script

        def viewport(w, h=900):
            """用 CDP 覆写视口尺寸。

            不要用 drv.set_window_size()：headless 下每次真正改窗口尺寸极慢
            （实测 12 个宽度 × 3 个页面跑不完 600s），而 CDP 的
            Emulation.setDeviceMetricsOverride 是纯渲染层覆写，毫秒级，
            且同样会触发媒体查询与 innerWidth 变化。
            """
            drv.execute_cdp_cmd('Emulation.setDeviceMetricsOverride',
                                {'width': w, 'height': h,
                                 'deviceScaleFactor': 1, 'mobile': False})

        # 先造出"最坏情况"的上下文：选中一个较长样品名 + 一个项目，
        # 让胶囊内容尽量宽（样品名最长 50 字符）。
        for page in PAGES:
            print(f'\n=== {page} ===')
            viewport(1920)
            drv.get(base + page)
            time.sleep(1.4)
            names = W("return await (await fetch('/api/samples')).text()")
            import json as _json
            try:
                lst = _json.loads(names)
            except Exception:
                lst = []
            if lst:
                longest = max((x['name'] for x in lst), key=len)
                projs = [x.get('project') for x in lst if x.get('project')]
                W(f"window.VP_CTX.setSample({_json.dumps(longest)});"
                  f"window.VP_CTX.setProject({_json.dumps(projs[0] if projs else '')})")
                time.sleep(0.3)
                # 加宽到 50 字符，压满胶囊
                W("document.getElementById('navCtxSample').querySelector('.v')"
                  ".textContent = 'X'.repeat(50)")
            for w in WIDTHS:
                t0 = time.time()
                viewport(w)
                m = W(MEASURE)
                dt = time.time() - t0
                tag = f'{w}px'
                if dt > 3:
                    print(f'      [warn] {tag} 测量耗时 {dt:.1f}s')
                if not m or m.get('error'):
                    check(f'{page} {tag} 能取到导航', False, str(m))
                    continue
                det = (f"navH={m['navH']} navRows={m['navRows']} "
                       f"linkRows={m['linkRows']} ovY={m['overflowY']} "
                       f"ovX={m['overflowX']}")
                check(f'{page} {tag} 横条不纵向溢出', m['overflowY'] <= 1, det)
                check(f'{page} {tag} 无子元素跑到横条外',
                      not m['outside'], f"{det} outside={m['outside'][:3]}")
                check(f'{page} {tag} 页面无横向溢出',
                      m['overflowX'] <= 1, det)
                check(f'{page} {tag} 上下文/语言/目录入口都还在且可见',
                      all(m[k] and m[k]['w'] > 0 and m[k]['h'] > 0
                          for k in ('chip', 'reset', 'lang', 'icon')),
                      f"chip={m['chip']} reset={m['reset']} "
                      f"lang={m['lang']} icon={m['icon']}")
                check(f'{page} {tag} {EXPECT_LINKS} 个一级链接全在',
                      m['links'] == EXPECT_LINKS,
                      f"{m['links']}（期望 {EXPECT_LINKS}）")
                # 行数预算按宽度分档：640px 下两行物理放不下（14 个一级项
                # 实测约 1285px，两行上限只有 ~1224px，纯几何），允许 3 行，
                # 但"被横条完整包住"由上面的 outside 断言另行守住。
                row_limit = 3 if w <= 640 else 2
                check(f'{page} {tag} 一级链接最多 {row_limit} 行（不允许再碎）',
                      m['linkRows'] <= row_limit, det)            # 宽屏必须回到"品牌与导航同一行、横条 54px"
            viewport(1920)
            m = W(MEASURE)
            if m and not m.get('error'):
                check(f'{page} 1920px 单行且横条高度≈54px',
                      m['navRows'] == 1 and 52 <= m['navH'] <= 58,
                      f"navRows={m['navRows']} navH={m['navH']}")
                check(f'{page} 1920px 一级链接单行',
                      m['linkRows'] == 1, f"linkRows={m['linkRows']}")
        drv.execute_cdp_cmd('Emulation.clearDeviceMetricsOverride', {})
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
        print(f'✘ {len(FAIL)} 项未通过:')
        for f in FAIL[:25]:
            print('   - ' + f)
        if len(FAIL) > 25:
            print(f'   … 另有 {len(FAIL) - 25} 项')
        return 1
    print('✔ 全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
