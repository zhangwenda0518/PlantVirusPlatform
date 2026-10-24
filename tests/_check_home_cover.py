# -*- coding: utf-8 -*-
"""启动页封面回归：真浏览器打开 8765 的首页，截浅色/深色两版封面。

用法: python tests/_check_home_cover.py [out_dir]

判据（任一不满足就打印 FAIL，退出码 1）:
  1. 首页 200 且 .hero.cover 存在；
  2. 封面已无插图（.cover-art 不应再存在，主视觉位图已按需求删除）；
  3. 封面文字非空（i18n 键漏配时这里会空）；
  4. 已按需求移除四张叙事卡与平台架构总览（.cover-card / .arch 不应存在）；
  5. 封面文字没溢出容器（右边界不出视口，横向无滚动条）。
"""
import io
import sys
from pathlib import Path

BASE = 'http://127.0.0.1:8765/'
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else None
if OUT is None:
    OUT = Path(__file__).resolve().parent.parent / 'run' / 'cover_shots'


def main() -> int:
    import urllib.request
    with urllib.request.urlopen(BASE, timeout=10) as r:
        if r.status != 200:
            print(f'FAIL 首页 HTTP {r.status}')
            return 1
    OUT.mkdir(parents=True, exist_ok=True)

    from playwright.sync_api import sync_playwright
    bad = []
    with sync_playwright() as pw:
        b = pw.chromium.launch(headless=True)
        for theme in ('light', 'dark'):
            ctx = b.new_context(viewport={'width': 1440, 'height': 1000},
                                device_scale_factor=2)
            pg = ctx.new_page()
            if theme == 'dark':
                # 与 home.html 内联脚本同一条键：先置位再进页面
                pg.add_init_script(
                    "try{localStorage.setItem('vp_theme','dark')}catch(e){}")
            pg.goto(BASE, # 对话坞常驻 EventSource 使 networkidle 永不触发（与其它巡检一致改 load）
                wait_until='load')
            pg.wait_for_timeout(700)

            if not pg.locator('.hero.cover').count():
                bad.append('缺少 .hero.cover')
            # 封面插图已按需求删除：留着说明 home.html 或 CSS 没清干净
            if pg.locator('.cover-art').count():
                bad.append('封面仍有 .cover-art 插图（应已删除）')
            lead = (pg.locator('.hero.cover').inner_text() or '').strip()
            if len(lead) < 12:
                bad.append(f'封面文案过短/为空: {lead!r}')
            # 2026-10-23 按需求移除：静态叙事卡与重复的架构总览
            if pg.locator('.cover-card').count():
                bad.append('仍有 .cover-card 叙事卡（应已移除）')
            if pg.locator('.arch').count():
                bad.append('仍有 .arch 架构总览（应已移除）')

            # 溢出检查：封面的实际右边界不能越过视口
            box = pg.locator('.hero.cover').bounding_box()
            if box and box['x'] + box['width'] > pg.viewport_size['width'] + 1:
                bad.append(f'.hero.cover 横向溢出 {box["x"] + box["width"]:.0f}px')
            scroll = pg.evaluate('document.documentElement.scrollWidth')
            if scroll > pg.viewport_size['width'] + 1:
                bad.append(f'页面出现横向滚动 {scroll}px')

            out = OUT / f'home_cover_{theme}.png'
            pg.screenshot(path=str(out), clip={'x': 0, 'y': 0,
                                               'width': pg.viewport_size['width'],
                                               'height': 760})
            print(f'  截图 {out.name}')
            ctx.close()
        b.close()

    if bad:
        print('FAIL')
        for x in bad:
            print('  -', x)
        return 1
    print('OK 首页封面通过（浅色/深色）')
    return 0


if __name__ == '__main__':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')
    sys.exit(main())
