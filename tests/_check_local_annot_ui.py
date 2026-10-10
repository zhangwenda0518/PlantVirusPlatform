# -*- coding: utf-8 -*-
"""真浏览器验证 t-hom 的「本地引擎」按钮链路（Playwright）。

/tools 一次只显示一个二级模块，所以用 ?g=annotate#t-hom 进入。
覆盖：示例自动填入 → 引擎默认 local → 点 ▶ BLASTN →
命中表出现且标注「本地库（离线）」。
（原 t-cdd 段已随 CDD 卡下线移除，2026-09-29。）
用法：python tests/_check_local_annot_ui.py
"""
import sys

from playwright.sync_api import sync_playwright
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


BASE = 'http://127.0.0.1:8765/tools?g=annotate'


def _rows(pg, body_id):
    pg.wait_for_selector(f'#{body_id} table tbody tr', timeout=300000)
    return (pg.eval_on_selector_all(f'#{body_id} table tbody tr',
                                    'els=>els.length'),
            pg.inner_text(f'#{body_id}'))


def main() -> int:
    fails = []
    with sync_playwright() as pw:
        b = pw.chromium.launch(headless=True)
        pg = b.new_page(viewport={'width': 1600, 'height': 1000})
        errs = []
        pg.on('console', lambda m: errs.append(m.text) if m.type == 'error' else None)
        pg.on('pageerror', lambda e: errs.append(str(e)))

        pg.goto(BASE, wait_until='domcontentloaded')
        pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")

        # ---------- t-hom（BLASTN 本地） ----------
        pg.goto(BASE + '#t-hom', wait_until='domcontentloaded')
        pg.wait_for_selector('#hom_engine', state='visible', timeout=30000)
        v = pg.input_value('#hom_engine')
        print(f'  hom_engine 默认 = {v}')
        if v != 'local':
            fails.append('hom_engine 默认不是 local')
        pg.wait_for_function(
            "() => { const t = document.getElementById('homSeqText');"
            " return t && t.value.trim().length > 0; }", timeout=30000)
        txt = pg.input_value('#homSeqText')
        print(f'  homSeqText 自动填入 {len(txt)} 字符，首行 {txt.splitlines()[0][:40]!r}')
        pg.click('#t-hom button:has-text("▶ BLASTN")')
        rows, head = _rows(pg, 'homBody')
        print(f'  BLASTN 命中 {rows} 行 · 表头含「本地库」: {"本地库" in head}')
        if rows < 1 or '本地库' not in head:
            fails.append('BLASTN 本地结果异常')
        # 切到在线再切回（只验证选择器，不真跑在线）
        pg.select_option('#hom_engine', 'online')
        print(f'  切换在线后 = {pg.input_value("#hom_engine")}')
        pg.select_option('#hom_engine', 'local')

        if errs:
            print('  控制台错误:')
            for e in errs[:6]:
                print('   ', e[:150])
        b.close()
    print('ALL OK' if not fails else 'FAILED: ' + ', '.join(fails))
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
