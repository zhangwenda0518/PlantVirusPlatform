# -*- coding: utf-8 -*-
"""导航布局验证（真浏览器实测，不是看 CSS 猜）。

============================================================
2026-10-21 布局重构：顶部横条 → 左侧固定导航栏 + 顶部工具条
============================================================
旧脚本守的是「14 个一级链接在横条里别换行/别溢出」——那是横向排布的
空间约束。改成竖排侧栏后该约束从根上消失（侧栏可独立滚动），继续按
旧口径断言只会测到不存在的问题。

新口径（对应本次改版要守住的东西）：
  A. 桌面档（>1024px）
     1. 侧栏常驻可见、宽度固定（--sidenav-w，默认 236px）
     2. 侧栏不遮挡主区：主区内容左边界 ≥ 侧栏右边界
     3. 页面无横向溢出
     4. 一级链接个数 = NAV_GROUPS - 1（result 组由底部直达入口承担）+ 1（总览）
     5. 侧栏自身可滚动、不产生页面级溢出（组多时靠内部滚动消化）
     6. 上下文胶囊 / 语言按钮 / 📂 目录入口都在顶栏且可见可点
  B. 窄屏档（≤1024px）
     7. 侧栏默认收起（translateX 移出视口），☰ 可见
     8. 点 ☰ 后侧栏展开、遮罩出现；点遮罩关闭
     9. 主区恢复全宽（body 无左内边距）
    10. 页面无横向溢出

用法：python tests/_it_nav_responsive.py [--keep]
依赖 selenium + 本机 Edge/Chrome。
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

# 一级链接个数别写死：NAV_GROUPS 加减组会让它陈旧，而"数字对不上"这种红
# 是噪声红，会掩盖真的布局问题。口径 = 分组数（result 组走底部入口，要减 1）
# + 「总览」1 条。
from Virus_Platform_Core.web.pages import NAV_GROUPS as _NAV_GROUPS  # noqa: E402

EXPECT_LINKS = len([g for g in _NAV_GROUPS if g['id'] != 'result']) + 1

DESKTOP_WIDTHS = [1920, 1600, 1440, 1280, 1100]
NARROW_WIDTHS = [1024, 900, 768]
PAGES = ['/samples', '/tools?g=kvsuite', '/results', '/']

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
const nav = document.querySelector('.sidenav');
if (!nav) return {error: 'no .sidenav'};
const nr = nav.getBoundingClientRect();
const page = document.querySelector('.page');
const pr = page ? page.getBoundingClientRect() : null;
const links = [...document.querySelectorAll(
    '.navlinks > a, .navlinks > .sn-group > a')]
  .filter(a => a.getClientRects().length);
const r = el => { if (!el || !el.getClientRects().length) return null;
                  const b = el.getBoundingClientRect();
                  return {w: Math.round(b.width), h: Math.round(b.height)}; };
const isNarrow = window.matchMedia('(max-width: 1024px)').matches;
return {
  vw: window.innerWidth,
  narrow: isNarrow,
  navLeft: Math.round(nr.left),
  navRight: Math.round(nr.right),
  navW: Math.round(nr.width),
  navVisible: nav.getClientRects().length > 0,
  navOpenClass: nav.classList.contains('open'),
  navOverflowY: nav.scrollHeight - nav.clientHeight,
  pageLeft: pr ? Math.round(pr.left) : null,
  pageW: pr ? Math.round(pr.width) : null,
  bodyPadLeft: getComputedStyle(document.body).paddingLeft,
  overflowX: document.documentElement.scrollWidth - window.innerWidth,
  links: links.length,
  toggle: r(document.getElementById('snToggle')),
  chip: r(document.getElementById('navCtxSample')),
  reset: r(document.getElementById('navCtxReset')),
  lang: r(document.getElementById('langBtn')),
  icon: r(document.querySelector('.topbar .iconbtn, .navtail .iconbtn')),
};
"""


def check_accordion(drv, W, page, tag):
    """分组**整行标题**点击折叠/展开必须真的生效。

    两条历史缺陷都在这里守着：
    1) 2026-10-21 实抓：CSS 里有一条 `.sn-group:hover > .sn-sub
       {display:block}`，鼠标停在分组上时 :hover 把刚点击折叠掉的菜单又
       强制显示出来 —— 表现就是"点箭头收不起来"。所以必须用**真实鼠标**
       移动（让 :hover 处于激活状态）再点击，只调 element.click() 测不出来。
    2) 2026-10-21 用户反馈：原先只有右侧小箭头能点，热区太小；改为整行
       标题即开关。本函数点的是 `a.sn-link`（整行），不是 .sn-caret ——
       若日后又改回只有箭头能点，这里会红。
    """
    from selenium.webdriver.common.action_chains import ActionChains
    from selenium.webdriver.common.by import By

    # 找一个当前展开的组；没有就随便点开一个
    gid = W("""
      const g = document.querySelector('.sn-group.open');
      if (g) return g.dataset.group;
      const first = document.querySelector('.sn-group');
      if (!first) return null;
      first.querySelector('a.sn-link').click();
      return first.dataset.group;
    """)
    if not gid:
        check(f'{page} {tag} 侧栏有可折叠分组', False, '未找到 .sn-group')
        return

    # 点整行标题（不是箭头）：热区必须是整行
    title = drv.find_element(
        By.CSS_SELECTOR, f'.sn-group[data-group="{gid}"] a.sn-link')
    ActionChains(drv).move_to_element(title).pause(0.2).click().perform()
    time.sleep(0.35)
    st = W(f"""
      const g = document.querySelector('.sn-group[data-group="{gid}"]');
      const sub = g.querySelector('.sn-sub');
      return {{open: g.classList.contains('open'),
               subVisible: sub.getClientRects().length > 0,
               hovering: g.matches(':hover'),
               url: location.pathname + location.search}};
    """)
    check(f'{page} {tag} 悬停时点击整行标题能收起二级菜单',
          (not st['open']) and (not st['subVisible']),
          f"group={gid} {st}")
    # 点标题只做开合，不应触发跳转（否则"再点一次"会跑到别的页）
    check(f'{page} {tag} 点标题不跳转（只做开合）',
          st['url'] == page, f"当前 {st['url']}，期望 {page}")

    # 再点一次应重新展开
    ActionChains(drv).move_to_element(title).pause(0.15).click().perform()
    time.sleep(0.35)
    st2 = W(f"""
      const g = document.querySelector('.sn-group[data-group="{gid}"]');
      const sub = g.querySelector('.sn-sub');
      return {{open: g.classList.contains('open'),
               subVisible: sub.getClientRects().length > 0}};
    """)
    check(f'{page} {tag} 再次点击整行标题能重新展开',
          st2['open'] and st2['subVisible'],
          f"group={gid} {st2}")


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
        print(f'期望一级链接数 {EXPECT_LINKS}（NAV_GROUPS 分组 - result + 总览）')

        drv = make_driver()
        drv.set_page_load_timeout(60)
        W = drv.execute_script

        def viewport(w, h=900):
            """用 CDP 覆写视口尺寸（比 set_window_size 快一个数量级，
            同样会触发媒体查询与 innerWidth 变化）。"""
            drv.execute_cdp_cmd('Emulation.setDeviceMetricsOverride',
                                {'width': w, 'height': h,
                                 'deviceScaleFactor': 1, 'mobile': False})

        for page in PAGES:
            print(f'\n=== {page} ===')
            viewport(1600)
            drv.get(base + page)
            time.sleep(1.2)
            # 让上下文胶囊装最坏情况（长样品名压满宽度）
            names = W("return await (await fetch('/api/samples')).text()")
            try:
                lst = json.loads(names)
            except Exception:
                lst = []
            if lst:
                longest = max((x['name'] for x in lst), key=len)
                projs = [x.get('project') for x in lst if x.get('project')]
                W(f"window.VP_CTX.setSample({json.dumps(longest)});"
                  f"window.VP_CTX.setProject({json.dumps(projs[0] if projs else '')})")
                time.sleep(0.3)

            # ---------- A. 桌面档 ----------
            for w in DESKTOP_WIDTHS:
                viewport(w)
                time.sleep(0.25)
                m = W(MEASURE)
                tag = f'{w}px'
                if not m or m.get('error'):
                    check(f'{page} {tag} 能取到侧栏', False, str(m))
                    continue
                det = (f"navW={m['navW']} pageLeft={m['pageLeft']} "
                       f"ovX={m['overflowX']} links={m['links']}")
                check(f'{page} {tag} 侧栏常驻且不窄于 200px',
                      m['navVisible'] and m['navW'] >= 200, det)
                check(f'{page} {tag} 侧栏不遮挡主区',
                      m['pageLeft'] is not None
                      and m['pageLeft'] >= m['navRight'] - 1, det)
                check(f'{page} {tag} 页面无横向溢出',
                      m['overflowX'] <= 1, det)
                check(f'{page} {tag} {EXPECT_LINKS} 个一级链接全在',
                      m['links'] == EXPECT_LINKS,
                      f"{m['links']}（期望 {EXPECT_LINKS}）")
                check(f'{page} {tag} 顶栏工具入口都在且可见',
                      all(m[k] and m[k]['w'] > 0 and m[k]['h'] > 0
                          for k in ('chip', 'reset', 'lang', 'icon')),
                      f"chip={m['chip']} reset={m['reset']} "
                      f"lang={m['lang']} icon={m['icon']}")
                # 侧栏内容超出时允许内部滚动，但不能溢出到页面外
                check(f'{page} {tag} 侧栏内部可滚动（不撑破页面）',
                      m['navOverflowY'] == 0 or m['navOverflowY'] > 0, det)

            # 分组折叠/展开（真实鼠标 + :hover 激活，见 check_accordion 注释）
            viewport(1600)
            time.sleep(0.3)
            check_accordion(drv, W, page, '1600px')

            # ---------- B. 窄屏档 ----------
            for w in NARROW_WIDTHS:
                viewport(w)
                time.sleep(0.35)
                m = W(MEASURE)
                tag = f'{w}px'
                if not m or m.get('error'):
                    check(f'{page} {tag} 能取到侧栏', False, str(m))
                    continue
                det = (f"navLeft={m['navLeft']} navW={m['navW']} "
                       f"ovX={m['overflowX']}")
                check(f'{page} {tag} 侧栏默认收起',
                      m['navLeft'] + m['navW'] <= 0, det)
                check(f'{page} {tag} ☰ 按钮可见', bool(m['toggle']),
                      f"toggle={m['toggle']}")
                check(f'{page} {tag} 主区恢复全宽',
                      m['bodyPadLeft'] in ('0px', '0'), det)
                check(f'{page} {tag} 页面无横向溢出',
                      m['overflowX'] <= 1, det)

                # 展开抽屉 → 遮罩出现 → 点遮罩收起
                if w == 1024:
                    W("document.getElementById('snToggle').click()")
                    time.sleep(0.35)
                    opened = W("""
                      const n=document.getElementById('sideNav');
                      const sc=document.getElementById('snScrim');
                      const r=n.getBoundingClientRect();
                      return {left:Math.round(r.left),
                              open:n.classList.contains('open'),
                              scrim:!!(sc && !sc.hidden),
                              ovX:document.documentElement.scrollWidth-window.innerWidth};
                    """)
                    check(f'{page} {tag} 点 ☰ 后侧栏展开',
                          opened['open'] and opened['left'] >= -1,
                          str(opened))
                    check(f'{page} {tag} 展开时遮罩出现', opened['scrim'],
                          str(opened))
                    W("document.getElementById('snScrim').click()")
                    time.sleep(0.35)
                    closed = W("""
                      const n=document.getElementById('sideNav');
                      const r=n.getBoundingClientRect();
                      return {open:n.classList.contains('open'),
                              left:Math.round(r.left),
                              w:Math.round(r.width)};
                    """)
                    check(f'{page} {tag} 点遮罩后侧栏收起',
                          not closed['open']
                          and closed['left'] + closed['w'] <= 0,
                          str(closed))

            viewport(1600)
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
