"""真浏览器验证 /genome「打开即出图」（Playwright，需平台已在 8765 运行）。

两段检查：
  A) 纯默认路径：打开页面后不做任何操作，等示例自动填充触发出图；
  B) 交互模式开关：取消勾选 → 回到 standby；重新勾选 → 再次出图。

用法：python tests/_genome_default_check.py
"""
import os
import sys
import tempfile
import time

from playwright.sync_api import sync_playwright
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


URL = 'http://127.0.0.1:8765/genome'
# 截图落盘到系统临时目录（此前硬编码了本机用户名路径，换机即失效）
SHOT = os.path.join(tempfile.gettempdir(), 'vp_genome_default.png')


def main() -> int:
    fails = []
    with sync_playwright() as pw:
        b = pw.chromium.launch(headless=True)
        pg = b.new_page(viewport={'width': 1600, 'height': 1000})
        errs = []
        pg.on('console', lambda m: errs.append(m.text) if m.type == 'error' else None)
        pg.on('pageerror', lambda e: errs.append(str(e)))

        pg.goto(URL, wait_until='domcontentloaded')
        pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")
        pg.reload(wait_until='domcontentloaded')
        pg.wait_for_selector('#gp_interactive', timeout=15000)

        # ---------- A. 纯默认路径 ----------
        t0 = time.time()
        try:
            pg.wait_for_selector('#gp_preview_body svg', timeout=90000)
            print(f'  [A] 打开即出图，耗时 {time.time() - t0:.1f}s（无任何人工操作）')
        except Exception as e:                     # noqa: BLE001
            fails.append('纯默认路径未出图')
            print(f'  [A] ✗ 未出图: {str(e)[:120]}')

        for name, cond in (
            ('交互模式默认勾选', pg.is_checked('#gp_interactive')),
            ('预览区默认可见', pg.is_visible('#gp_preview')),
            ('standby 提示默认隐藏', not pg.is_visible('#gp_preview_standby')),
            ('缩放条可见', pg.is_visible('#gp_zoom')),
        ):
            print(f'  [A] {"ok " if cond else "FAIL "}{name}')
            if not cond:
                fails.append(name)

        box = pg.eval_on_selector(
            '#gp_preview_body svg',
            'el=>{const r=el.getBoundingClientRect();'
            'return {w:Math.round(r.width),h:Math.round(r.height)}}')
        print(f'  [A] SVG 渲染尺寸 {box["w"]}x{box["h"]}，'
              f'输入={pg.input_value("#gp_ann")!r}')
        pg.screenshot(path=SHOT)
        print(f'  [A] 截图: {SHOT}')

        # ---------- B. 开关往返 ----------
        pg.uncheck('#gp_interactive')
        pg.dispatch_event('#gp_interactive', 'change')
        pg.wait_for_timeout(300)
        off_ok = pg.is_visible('#gp_preview_standby') and not pg.is_visible('#gp_preview')
        print(f'  [B] {"ok " if off_ok else "FAIL "}取消勾选 → 回到 standby')
        pg.check('#gp_interactive')
        pg.dispatch_event('#gp_interactive', 'change')
        try:
            pg.wait_for_selector('#gp_preview_body svg', timeout=90000)
            on_ok = True
        except Exception:                          # noqa: BLE001
            on_ok = False
        print(f'  [B] {"ok " if on_ok else "FAIL "}重新勾选 → 再次出图')
        if not (off_ok and on_ok):
            fails.append('开关往返异常')

        if errs:
            print('  控制台错误:')
            for e in errs[:8]:
                print('   ', e[:160])
        b.close()

    print(('ALL OK' if not fails else 'FAILED: ' + ', '.join(fails)))
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
