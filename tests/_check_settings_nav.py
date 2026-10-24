# -*- coding: utf-8 -*-
"""设置页侧栏锚点对齐巡检（Playwright，真浏览器）。

============================================================
2026-10-21 修的问题（用户报障："点左侧功能跳转位置很乱、没有对齐"）
============================================================
三个独立缺陷叠在一起，看起来都像"跳转乱"：

1. **多列网格**：`.settings-grid` 是三列 auto-fit，而侧栏按文档顺序逐个
   列出 —— 「输出/输入/数据库目录」三张卡在同一行、top 完全相同（实测
   都是 167），点哪个都滚到同一处。
2. **锚点被顶栏压住**：`scroll-margin-top:14px` 而吸顶顶栏高 54px，
   滚动后目标卡片顶边贴在视口 0 处，标题被吃掉 40px（实测 top=14）。
3. **侧栏顺序 ≠ 文档顺序**：`sec-storage` / `sec-paths` 在 DOM 里排在
   「通用设置 / 数据库 / 分析默认参数」之后，侧栏却把它们归在
   「目录与存储」组前面 —— 点「平台目录」时按位置算出的当前区块是
   夹在中间的 `sec-defaults`，高亮跳到别人身上。

修复：设置页改单列流式；`scroll-margin-top` 引用 `--topbar-h`；
侧栏按 DOM 实际位置排序并同步调整了模板物理顺序。

本脚本断言（每个侧栏条目）：
  A. 点击后目标区块顶边落在顶栏下方 12~30px 区间（不被压住、也不偏下）
  B. 点击后高亮的是被点的那一项（不是相邻项）
  C. 不同条目滚动到**不同**位置（排除"多列同一行"回归）
  D. 侧栏顺序与文档顺序一致
  E. 全程无 JS 报错

用法：python tests/_check_settings_nav.py
"""
import os
import sys
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from playwright.sync_api import sync_playwright
except ImportError:                     # pragma: no cover
    sync_playwright = None

FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name
          + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


def main():
    if sync_playwright is None:
        print('[SKIP] 未安装 playwright，跳过真浏览器验证')
        return 0

    from werkzeug.serving import make_server
    import threading
    os.environ.setdefault('VP_NO_RECOVER', '1')
    import app as app_module
    srv = make_server('127.0.0.1', 18788, app_module.app, threaded=True)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = 'http://127.0.0.1:18788'
    time.sleep(1.0)

    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(
                headless=True,
                args=['--no-sandbox', '--disable-dev-shm-usage', '--disable-gpu'])
            ctx = browser.new_context(viewport={'width': 1600, 'height': 900})
            page = ctx.new_page()
            errors = []
            page.on('pageerror', lambda e: errors.append(str(e)))
            page.on('console',
                    lambda m: errors.append(m.text) if m.type == 'error' else None)

            page.goto(base + '/settings', wait_until='load')
            page.wait_for_timeout(1500)

            # ---- 布局前提：设置卡片必须是单列 ----
            cols = page.evaluate("""
              () => {
                const g = document.querySelector('.se-main > .settings-grid');
                if (!g) return null;
                const tops = [...g.querySelectorAll(':scope > .card')]
                  .map(c => Math.round(c.getBoundingClientRect().top + window.scrollY));
                return {count: tops.length, unique: new Set(tops).size};
              }""")
            check('设置卡片为单列（每张卡独占一行）',
                  cols and cols['count'] == cols['unique'],
                  f'卡片数={cols and cols["count"]} 不同 top 数={cols and cols["unique"]}')

            # ---- 侧栏顺序 vs 文档顺序 ----
            order = page.evaluate("""
              () => {
                const docTop = el => { let y=0, n=el;
                  while (n) { y += n.offsetTop||0; n = n.offsetParent; } return y; };
                const side = [...document.querySelectorAll('#seSide a')]
                  .map(a => a.dataset.for);
                const doc = [...side].sort((a, b) =>
                  docTop(document.getElementById(a)) - docTop(document.getElementById(b)));
                return {side, doc};
              }""")
            check('侧栏顺序与文档顺序一致',
                  order['side'] == order['doc'],
                  f"侧栏={order['side']}")

            ids = order['side']
            print(f'  侧栏条目 {len(ids)} 项')

            positions = {}
            for sid in ids:
                page.evaluate(
                    f"document.querySelector('#seSide a[data-for=\"{sid}\"]').click()")
                # 平滑滚动最长 700ms 锁 + 余量
                page.wait_for_timeout(1100)
                r = page.evaluate(f"""
                  () => {{
                    const el = document.getElementById('{sid}');
                    const b = el.getBoundingClientRect();
                    const tb = document.querySelector('.topbar').getBoundingClientRect();
                    const on = document.querySelector('#seSide a.on');
                    return {{top: Math.round(b.top),
                             topbarBottom: Math.round(tb.bottom),
                             active: on ? on.dataset.for : null,
                             y: Math.round(window.scrollY)}};
                  }}""")
                positions[sid] = r['y']
                offset = r['top'] - r['topbarBottom']
                check(f'{sid} 落点不被顶栏遮挡且贴合（偏移 12~30px）',
                      12 <= offset <= 30,
                      f"top={r['top']} topbarBottom={r['topbarBottom']} 偏移={offset}")
                check(f'{sid} 点击后高亮的就是该项',
                      r['active'] == sid, f"高亮={r['active']}")

            # ---- 不同条目必须落在不同滚动位置 ----
            dup = len(positions) - len(set(positions.values()))
            check('各条目滚动位置互不相同（无"同排卡片挤在一处"）',
                  dup == 0, f'{dup} 组重复: {positions}')

            # ---- 深链 #hash 也要对齐 ----
            page.goto(base + '/settings#aiCard', wait_until='load')
            page.wait_for_timeout(1800)
            h = page.evaluate("""
              () => {
                const el = document.getElementById('aiCard');
                const b = el.getBoundingClientRect();
                const tb = document.querySelector('.topbar').getBoundingClientRect();
                return {top: Math.round(b.top), topbarBottom: Math.round(tb.bottom)};
              }""")
            check('带 #hash 直接进入也落在顶栏下方',
                  h['top'] >= h['topbarBottom'],
                  f"top={h['top']} topbarBottom={h['topbarBottom']}")

            check('全程无 JS 报错', not errors, f'{errors[:3]}')
            ctx.close()
            browser.close()
    finally:
        srv.shutdown()

    print()
    if FAIL:
        print(f'✘ {len(FAIL)} 项未通过:')
        for f in FAIL[:20]:
            print('   - ' + f)
        return 1
    print('✔ 设置页侧栏锚点全部对齐')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
