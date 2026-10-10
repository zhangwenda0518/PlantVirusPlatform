# -*- coding: utf-8 -*-
"""验证「粘贴序列」框改为**独立内联大文本框**（尺寸对齐 BLAST 卡）。

检查点：
1. 点 📋 后，文件输入行下方出现独立 textarea，高度 ≥ 220px（与 BLAST 卡一致）；
2. 粘贴内容后自动写入临时文件，输入框回填 uploads/paste_*.fasta 路径；
3. 再点一次 📋 收起；Esc 收起。
（参照卡原为 t-cdd，CDD 卡下线后改用 t-hom 的文本框，2026-09-29。）
用法：python tests/_check_paste_box.py
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


ORF = 'http://127.0.0.1:8765/orf'
HOM = 'http://127.0.0.1:8765/tools?g=annotate#t-hom'
SEQ = '>demo_seq\nATGGCTAGCTAGCTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGCTAA\n'


def main() -> int:
    fails = []
    with sync_playwright() as pw:
        b = pw.chromium.launch(headless=True)
        pg = b.new_page(viewport={'width': 1600, 'height': 1200})
        errs = []
        pg.on('pageerror', lambda e: errs.append(str(e)))

        # 参考尺寸：BLAST 卡的序列文本框
        pg.goto(HOM, wait_until='domcontentloaded')
        pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")
        pg.wait_for_selector('#homSeqText', state='visible', timeout=30000)
        hom_h = pg.eval_on_selector('#homSeqText',
                                    'e=>Math.round(e.getBoundingClientRect().height)')
        print(f'  参照 BLAST 文本框高度: {hom_h}px')

        # /orf 卡的 📋
        pg.goto(ORF, wait_until='domcontentloaded')
        pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")
        pg.wait_for_selector('#of_fa', timeout=30000)
        pg.click('#of_fa ~ button:has-text("📋"), button[onclick*="pasteSeq(\'of_fa\'"]')
        pg.wait_for_selector('#pastebox-of_fa textarea', timeout=15000)
        h = pg.eval_on_selector('#pastebox-of_fa textarea',
                                'e=>Math.round(e.getBoundingClientRect().height)')
        print(f'  内联粘贴框高度: {h}px（≥220 且与 BLAST 卡相当: '
              f'{h >= 220 and abs(h - hom_h) <= 60}）')
        if h < 220 or abs(h - hom_h) > 60:
            fails.append(f'粘贴框高度不合适 ({h} vs BLAST 卡 {hom_h})')

        # 粘贴 → 自动写入临时文件
        pg.fill('#pastebox-of_fa textarea', SEQ)
        pg.wait_for_function(
            "() => { const v = document.getElementById('of_fa').value;"
            " return v && /uploads[\\\\/]paste_/.test(v); }", timeout=20000)
        val = pg.input_value('#of_fa')
        print(f'  粘贴后输入框回填: {val}')
        if 'uploads/paste_' not in val.replace('\\', '/'):
            fails.append('未回填临时文件路径')

        # 再点 📋 收起
        pg.click('button[onclick*="pasteSeq(\'of_fa\'"]')
        pg.wait_for_timeout(300)
        gone = pg.query_selector('#pastebox-of_fa') is None
        print(f'  再点 📋 收起: {gone}')
        if not gone:
            fails.append('再次点击未收起')

        if errs:
            print('  页面错误:')
            for e in errs[:6]:
                print('   ', e[:150])
        b.close()
    print('ALL OK' if not fails else 'FAILED: ' + ', '.join(fails))
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
