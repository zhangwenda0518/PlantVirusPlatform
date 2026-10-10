# -*- coding: utf-8 -*-
"""真浏览器验证「功能注释」页的注释层选择（prot / pfam / bfvd，Playwright）。

覆盖：
  1. 三个勾选框默认全选；引擎/库下拉可用，「全部层」按钮存在
  2. 取消层1 → 引擎/参考库下拉禁用 + 出现提示；重新勾上 → 恢复
  3. 三个全取消 → 点运行只弹提示、不创建运行目录
  4. 只勾 pfam 用方式 B 实跑 → run.log 出现「启用层 pfam」与
     「层1 序列同源搜索按选择跳过」，summary.json 的 libs==['pfam']、engine 为空
（bfvd = BFVD 病毒域层，2026-09-29 起 CDD 的替代。）
用法：python tests/_check_annot_libs.py
"""
import json
import sys
import urllib.request

from playwright.sync_api import sync_playwright
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


ORIGIN = 'http://127.0.0.1:8765'
BASE = ORIGIN + '/annotation'
LIB_IDS = ('oa_lib_prot', 'oa_lib_pfam', 'oa_lib_bfvd')


def _get(path):
    with urllib.request.urlopen(ORIGIN + path, timeout=20) as r:
        return json.load(r)


def _text(path):
    with urllib.request.urlopen(ORIGIN + path, timeout=20) as r:
        return r.read().decode('utf-8', 'replace')


def main() -> int:
    fails = []

    def ck(cond, msg):
        print(('  ok  ' if cond else '  FAIL ') + msg)
        if not cond:
            fails.append(msg)

    with sync_playwright() as pw:
        b = pw.chromium.launch(headless=True)
        pg = b.new_page(viewport={'width': 1600, 'height': 1100})
        errs = []
        pg.on('console', lambda m: errs.append(m.text) if m.type == 'error' else None)
        pg.on('pageerror', lambda e: errs.append(str(e)))
        dialogs = []
        pg.on('dialog', lambda d: (dialogs.append(d.message), d.accept()))
        # 首次访问时 examples.js 会自动弹出「示例结果总览」遮罩（class=dlgmask），
        # 会挡住点击；置位 vp_example_seen 跳过（与 _check_local_annot_ui.py 一致）。
        pg.add_init_script(
            "try{localStorage.setItem('vp_example_seen','1')}catch(e){}")

        pg.goto(BASE, wait_until='domcontentloaded')
        pg.wait_for_selector('#oa_lib_prot', state='visible', timeout=30000)
        pg.evaluate("() => document.querySelectorAll('.dlgmask')"
                    ".forEach(m => { if (m.id !== 'dlgMask') m.remove(); })")

        print('[1] 默认状态')
        for i in LIB_IDS:
            ck(pg.is_checked('#' + i), f'{i} 默认勾选')
        ck(pg.is_enabled('#oa_engine'), '层1 开启时引擎下拉可用')
        ck(pg.is_enabled('#oa_db'), '层1 开启时参考库下拉可用')
        ck(not pg.is_visible('#oa_libs_warn'), '层1 开启时无警告')
        ck(pg.is_visible('button:has-text("全部层")'), '「全部层」按钮存在')

        print('[2] 取消层1')
        pg.uncheck('#oa_lib_prot')
        pg.wait_for_timeout(150)
        ck(pg.is_disabled('#oa_engine'), '层1 关闭后引擎下拉禁用')
        ck(pg.is_disabled('#oa_db'), '层1 关闭后参考库下拉禁用')
        ck(pg.is_visible('#oa_libs_warn'), '层1 关闭后出现提示')
        ck('层1' in pg.inner_text('#oa_libs_warn'), '提示文案提到层1')
        ck(not pg.is_visible('#oa_db_custom'), '层1 关闭时自备库行隐藏')
        pg.click('button:has-text("全部层")')
        pg.wait_for_timeout(150)
        ck(all(pg.is_checked('#' + i) for i in LIB_IDS), '「全部层」重新勾满三层')
        ck(pg.is_enabled('#oa_engine'), '恢复后引擎下拉可用')

        print('[3] 全不勾 → 应拦截')
        for i in LIB_IDS:
            pg.uncheck('#' + i)
        pg.wait_for_timeout(120)
        pg.fill('#oa_fa', 'examples/example_viral_contigs.fasta')
        before = [x['name'] for x in _get('/api/tool/runs')][:1]
        pg.click('button:has-text("▶ 预测并注释")')
        pg.wait_for_timeout(600)
        ck(any('至少勾选' in m or '注释层' in m for m in dialogs),
           f'全不勾时弹出提示（{dialogs[-1][:32] if dialogs else "无 dialog"}）')
        ck([x['name'] for x in _get('/api/tool/runs')][:1] == before,
           '全不勾时未创建运行目录')

        print('[4] 只勾 pfam 实跑（方式 B）')
        pg.check('#oa_lib_pfam')
        pg.evaluate("() => document.querySelectorAll('.dlgmask')"
                    ".forEach(m => { if (m.id !== 'dlgMask') m.remove(); })")
        pg.click('button:has-text("▶ 预测并注释")')
        pg.wait_for_selector('#toolrun-orfa', timeout=30000)
        page_txt = ''
        for _ in range(90):
            pg.wait_for_timeout(2000)
            page_txt = pg.inner_text('body')      # 含任务日志面板
            if '任务完成' in page_txt or '失败' in page_txt:
                break
        ck('任务完成' in page_txt, '页面出现「任务完成」')
        # 页面日志面板只渲染尾部若干行且长行被 CSS 截断，完整日志以 run.log 为准
        newest = [x['name'] for x in _get('/api/tool/runs')
                  if x['name'].startswith('orfa_')][0]
        log = _text(f'/tool_runs/{newest}/run.log')
        for key in ('启用层 pfam', '层1 序列同源搜索按选择跳过', 'Pfam-A-Viruses'):
            ck(key in log, f'run.log 出现「{key}」')
        for line in log.strip().splitlines():
            if '启用层' in line or '层1 序列同源' in line:
                print('  run.log | ' + line.strip()[:110])
        s = _get(f'/tool_runs/{newest}/04b_orf_annot/summary.json')
        print(f'  {newest}: libs={s.get("libs")} engine={s.get("engine")!r} '
              f'n_annotated={s.get("n_annotated")}/{s.get("n_orfs")}')
        ck(s.get('libs') == ['pfam'], "summary.json libs == ['pfam']")
        ck(not s.get('engine'), 'summary.json engine 为空（层1 未跑）')
        ck('层1 序列同源（RefSeq 病毒蛋白）' not in (s.get('libs_labels') or ['x']),
           'summary.json libs_labels 不含层1')
        ck(s.get('evidence') and 'seq' not in s['evidence'],
           'evidence 中无 seq（层1 未跑）')

        if errs:
            print('  控制台错误:')
            for e in errs[:6]:
                print('   ', e[:150])
            fails.append('浏览器控制台有错误')
        b.close()
    print('ALL OK' if not fails else 'FAILED: ' + '; '.join(fails))
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
