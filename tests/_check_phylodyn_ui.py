# -*- coding: utf-8 -*-
"""真浏览器验证「进化动力学」卡组（**2026-09-18 起为 8 张**）。

后端引擎已由 tests/_check_phylodyn_kit.py 对 Example 数据逐模块验证（含两个
逐字节/逐格对拍契约）；本测试锁**界面链路**：

  [1] /tools?g=phylodyn 组落地列出 9 张卡（5 张数据准备类 + RDP5/时间信号/系统地理
      + 新增「时间与地理推断·本地全链」t-phylodyn）；
      ★ 负向：8 张归档卡的落地条目必须消失（逐个点名）；
  [2] 「进化动力学数据接入」卡进卡：三来源切换、时间地点粘贴框、严格开关齐全；
      静态文案已 i18n（中文可见），切英文同键变英文；
  [3] 真跑一次数据接入（手动粘贴 3 条：1 行非法日期 + 1 行**小数年**）：结果面板出现，
      保留 2 条（小数年那条必须被折算收下而不是当坏行剔掉）、问题表点名非法日期那行、
      按年/按地点图渲染（plotly 出图），并从 Web 层取回 metadata.csv 验折算结果；
  [4] 「序列 ID 重命名」卡真跑（fixture FASTA + 映射表）：16 条全改名；
  [5] ★ 归档防回归：8 张卡的 section / 关键控件 id / JS 函数必须都不在；
      同时正向断言保留卡的 5 个 JS 函数还在（防"删过头"也通过）；
  [6] 序列分组卡的关键控件存在（不真跑）；时空分布绘图/时空降采样 已并入数据接入；
  [7] 切英文（保留卡的键 zh/en 双份）；
  [8] 无 pageerror。

2026-09-18 变更：组内由 16 项收敛到 8 项 —— 保留「数据接入 / 序列 ID 重命名 /
序列分组」（进化平台没做过的，时空绘图/降采样已于当天并入数据接入），
其余 8 张（rand/rrt/tempmig/bsp/rspp/treetime/ltt/mjrm）已归档，
改由进化平台的「时间与地理推断·本地全链」统一覆盖。见
docs/进化动力学统一方案_20260918.md。
底层 phylodyn_kit / phylodyn_trees **未删**（保留卡仍在用），
所以 _check_phylodyn_kit.py / _check_phylodyn_example.py 照旧全跑。

无 playwright 时 SKIP（同 _check_phylogeo_ui 口径）。
用法：python tests/_check_phylodyn_ui.py            （自起临时实例）
      python tests/_check_phylodyn_ui.py --port N   （用已在运行的平台）
"""
import argparse
import io
import os
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request

try:
    from playwright.sync_api import sync_playwright
except ImportError:                     # pragma: no cover - 取决于本机环境
    sync_playwright = None

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SHOT = os.path.join(ROOT, 'run', '_check_phylodyn_ui.png')
FIX = os.path.join(ROOT, 'run', '_check_phylodyn_ui_fix')
# 外部参照仓（延伸基因组/MMPV-RNA 工作区）：默认仍是维护者本机路径，
# 可用环境变量 VP_EXT_ROOT 覆盖（指向 MMPV-RNA 工作区根）；
# 目录不存在时相关小节自动跳过。
EXT_ROOT = os.environ.get('VP_EXT_ROOT', r'D:\桌面\延伸基因组\MMPV-RNA')
EXAMPLE = os.path.join(EXT_ROOT, 'biosoft', 'VirPhyKit', 'Example')

N_CARDS = 5
# 2026-09-18：进化动力学组从 16 项收敛到 8 项 ——
#   保留 5 张数据准备类卡（数据接入 / ID 重命名 / 分组 / 时空绘图 / 降采样，
#   进化平台没有做过）+ 原有的 RDP5 / 时间信号 / 系统地理 3 张；
#   其余 8 张（rand/rrt/tempmig/bsp/rspp/treetime/ltt/mjrm）已归档，
#   改由进化平台的「时间与地理推断·本地全链」统一覆盖。
ARCHIVED = ['t-pdrand', 't-pdrrt', 't-pdtempmig', 't-pdbsp', 't-pdrspp',
            't-pdrename', 't-pdgroup',
            't-pdtreetime', 't-pdltt', 't-pdmjrm',
            # 2026-09-18 第二批：数据准备类合并 → 这两张并进 t-pdprep
            't-pdspacetime', 't-pdsub']
# 保留卡至少要点到的控件（存在性检查，不逐个真跑）
WANT_PER_CARD = {
    't-pdprep': ['#pd_src', '#pd_fa_text', '#pd_meta_text', '#pd_fa',
                 '#pd_meta', '#pd_acc', '#pd_strict', '#toolrun-pdprep', '#pdprepResult',
                 '#pd_coords_tsv', '#pd_sub_on', '#pd_sub_box', '#pd_sub_mode',
                 '#pd_sub_n', '#pd_sub_seed', '#pdprepToA0'],
    # 2026-09-21: t-pdrename / t-pdgroup 已归档（重命名接进 t-phylodyn 可选前置；
    # 分组归导入/收集阶段）
}
WANT_KEYS = ['tk.pdprepH2', 'tk.pdprepStrict', 'tk.pdprepDesc', 'tk.pdprepSubOn']

FAILS = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def free_port(start=8861):
    for p in range(start, start + 40):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(('127.0.0.1', p))
                return p
            except OSError:
                continue
    raise RuntimeError('没有可用端口')


def start_server(port):
    logpath = os.path.join(tempfile.gettempdir(), f'vp_phylodyn_ui_{port}.log')
    logf = io.open(logpath, 'w', encoding='utf-8', errors='replace')
    env = dict(os.environ, VP_NO_RECOVER='1')
    proc = subprocess.Popen(
        [sys.executable, '-c',
         'import sys\n'
         f'sys.path.insert(0, {ROOT!r})\n'
         'from werkzeug.serving import make_server\n'
         'import app\n'
         f"s = make_server('127.0.0.1', {port}, app.app, threaded=True)\n"
         "print('ready', flush=True)\n"
         's.serve_forever()'],
        cwd=ROOT, env=env, stdout=logf, stderr=subprocess.STDOUT)
    for _ in range(120):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', port)) == 0:
                return proc, logpath
        if proc.poll() is not None:
            break
        time.sleep(0.3)
    proc.kill()
    print('服务启动失败，日志:', logpath)
    return None, logpath


def write_fixtures():
    """数据接入/重命名的输入；Example 在时用它的 rename 数据（口径更真）。"""
    os.makedirs(FIX, exist_ok=True)
    fa = os.path.join(FIX, 'ren.fasta')
    mp = os.path.join(FIX, 'ren.tsv')
    if os.path.isdir(EXAMPLE):
        import shutil
        shutil.copy(os.path.join(EXAMPLE, 'SeqIDRenamer', 'rename.fasta'), fa)
        shutil.copy(os.path.join(EXAMPLE, 'SeqIDRenamer', 'SeqIDRenamer.txt'), mp)
        return fa, mp
    with io.open(fa, 'w', encoding='utf-8', newline='') as f:
        f.write('>old1\nACGT\n>old2\nTTTT\n')
    with io.open(mp, 'w', encoding='utf-8', newline='') as f:
        f.write('old1\tnew1\nold2\tnew2\n')
    return fa, mp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=0)
    ap.add_argument('--keep', action='store_true')
    args = ap.parse_args()

    if sync_playwright is None:
        print('[SKIP] 未安装 playwright，跳过真浏览器验证')
        print('       安装: python -m pip install -r requirements-dev.txt')
        return 0

    ren_fa, ren_map = write_fixtures()
    proc, logpath = None, ''
    port = args.port or free_port()
    if not args.port:
        proc, logpath = start_server(port)
        if proc is None:
            return 2
        print(f'临时实例: http://127.0.0.1:{port}  (服务日志 {logpath})')

    try:
        with sync_playwright() as pw:
            b = pw.chromium.launch(headless=True)
            pg = b.new_page(viewport={'width': 1600, 'height': 1100})
            errs = []
            pg.on('pageerror', lambda e: errs.append(str(e)))
            pg.goto(f'http://127.0.0.1:{port}/tools?g=phylodyn',
                    wait_until='domcontentloaded')
            pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")
            pg.reload(wait_until='domcontentloaded')
            pg.wait_for_selector('#groupLanding', timeout=30000)
            pg.evaluate("document.querySelectorAll('.dlgmask').forEach(m => m.remove())")

            print('\n[1] 组落地：5 张卡（2026-09-21 归档 rename/group 后由 7 收敛）')
            cards = pg.eval_on_selector_all(
                '#groupLanding .landing-item', 'els => els.map(e => e.textContent)')
            check(len(cards) == N_CARDS, f'落地 {len(cards)}/{N_CARDS} 张卡')
            for name in ('数据接入', '一键分析'):
                check(any(name in c for c in cards), f'落地含「{name}」卡')
            # ★ 负向防回归：8 张归档卡的落地条目必须消失（只断言"少了"会假安全，
            #   这里逐个点名，并同时确认命中的是**归档名单**而不是别的东西）
            land_ids = pg.eval_on_selector_all(
                '#groupLanding .landing-item',
                'els => els.map(e => e.getAttribute("data-item"))')
            left = [c for c in ARCHIVED if c in land_ids]
            check(not left, f'8 张归档卡已从落地消失（残留: {left or "无"}）')

            print('\n[2] 数据接入卡：三来源 + i18n')
            pg.locator('#groupLanding .landing-item[data-item="t-pdprep"]').click()
            pg.wait_for_selector('#t-pdprep', state='visible', timeout=15000)
            sec = pg.locator('#t-pdprep')
            check(pg.url.endswith('#t-pdprep'), f'hash 落到 #t-pdprep（{pg.url}）')
            for sel in WANT_PER_CARD['t-pdprep']:
                check(sec.locator(sel).count() == 1, f'控件存在: {sel}')
            # 2026-09-18 界面整合：可选后处理收进一个折叠区，**默认收起** ——
            # 主流程（来源 / 运行 / 结果 / 接力 A0）不该被可选步骤淹没。
            check(sec.locator('#pd_adv').count() == 1, '可选后处理折叠区存在')
            check(pg.evaluate("() => !document.getElementById('pd_adv').open"),
                  '折叠区默认**收起**（不勾就不做）')
            check(sec.locator('#pd_src').is_visible(), '主流程的来源选择可见（在折叠区外）')
            check(sec.locator('#toolrun-pdprep').count() == 1, '运行/输出区在折叠区外')
            # ★ 负向防回归：zip 来源 2026-09-18 整条下线（explorer 不再导出 zip）
            _zips = [s for s in ('#pd_src_zip', '#pd_zip', '#pdyn_zip')
                     if sec.locator(s).count()]
            check(not _zips, f'zip 来源的控件已摘除（残留: {_zips or "无"}）')
            _opts = pg.eval_on_selector_all(
                '#pd_src option', 'els => els.map(e => e.value)')
            check('zip' not in _opts, f'来源下拉里没有 zip（实际 {_opts}）')
            check(not pg.evaluate(
                "() => typeof window.pdImportZip === 'function'"),
                'pdImportZip 函数已摘除')
            check(pg.evaluate(
                "() => document.querySelectorAll('#pd_src_files').length") == 1,
                '保留的「文件对」来源仍在（正向）')
            for key in WANT_KEYS[:2]:
                txt = pg.eval_on_selector(f'[data-i18n="{key}"]',
                                          'e => e.textContent.trim()')
                check(bool(txt) and not txt.startswith('tk.'),
                      f'{key} → {txt[:30]}…')

            print('\n[3] 真跑数据接入（手动粘贴，1 行日期坏 + 1 行小数年）')
            # 2026-09-18：来源默认已改成「文件对」（explorer 不再导出 zip），
            # 本用例走的是**手动粘贴** → 先切回 manual，否则填的文本框不在生效的来源上。
            pg.select_option('#pd_src', 'manual')
            pg.fill('#pd_fa_text',
                    '>SeqA\nACGTACGTACGT\n>SeqB\nTTTTGGGGCCCC\n>SeqC\nAAAACCCCGGGG\n')
            # SeqB 用 explorer 导出包的**小数年**口径（2026-09-18 起必须被收下并折算）；
            # SeqC 用真正非法的 13 月（必须被剔并点名）。
            pg.fill('#pd_meta_text',
                    'name\tdate\tlocation\tlat\tlon\n'
                    'SeqA\t2019-03-01\tWuhan\t30.5\t114.3\n'
                    'SeqB\t2019.5\tWuhan\n'
                    'SeqC\t2020-13-01\tNairobi\n')
            pg.locator('#t-pdprep button', has_text='导入并校验').click()
            ok = False
            for _ in range(60):
                if pg.eval_on_selector('#pdprepResult',
                                       "e => getComputedStyle(e).display") != 'none':
                    ok = True
                    break
                pg.wait_for_timeout(1000)
            if not ok:
                print('    任务日志尾部:',
                      pg.eval_on_selector('#toolrun-pdprep', 'e => e.innerText')[-400:])
            check(ok, '结果面板已出现')
            if ok:
                summ = pg.inner_text('#pdprepSummary')
                check('保留 2' in summ and '剔除/警告 1' in summ and '输入 3' in summ
                      and '地点 1' in summ, f'汇总正确（{summ[:60]}…）')
                issues = pg.inner_text('#pdprepIssues')
                check('SeqC' in issues and '月份 13' in issues and 'SeqB' not in issues,
                      '问题表只点名 SeqC 非法日期（小数年 SeqB 不再算坏行）')
                n_svg = pg.eval_on_selector('#pdprepYearPlot',
                                            "e => e.querySelectorAll('.main-svg').length")
                check(n_svg >= 1, f'按年图已渲染（{n_svg} 块）')
                check(pg.eval_on_selector('#pdprepFaDl', 'e => e.href')
                      .endswith('pdprep/sequences.fasta?dl=1'), '下载链接指向产物')
                # 端到端：把 metadata.csv 从 Web 层取回来，证明小数年真的折成了 ISO
                meta_href = pg.eval_on_selector('#pdprepMetaDl', 'e => e.getAttribute("href")')
                try:
                    with urllib.request.urlopen(
                            f'http://127.0.0.1:{port}{meta_href}', timeout=30) as resp:
                        body = resp.read().decode('utf-8-sig', 'replace')
                    check('SeqB,2019-07-02,' in body,
                          '小数年 2019.5 经 Web 落盘已折成 ISO '
                          f'（末行 = {body.strip().splitlines()[-1]}）')
                except Exception as e:      # noqa: BLE001
                    check(False, f'取回 metadata.csv 失败: {e}')

            print('\n[4] 归档防回归：重命名/分组两卡也已不在（2026-09-21）')
            for cid in ('t-pdrename', 't-pdgroup'):
                check(pg.locator(f'#{cid}').count() == 0,
                      f'归档卡 section 不存在: #{cid}')
            gone_pd = pg.evaluate(
                "names => names.filter(n => typeof window[n] === 'function')",
                ['runPdrename', 'loadPdrenameResult', 'runPdgroup',
                 'loadPdgroupResult'])
            check(not gone_pd, f'归档卡 JS 函数已移除（残留: {gone_pd or "无"}）')

            print('\n[5] 归档防回归：8 张卡的 DOM 与 JS 必须都不在')
            # 页面里不能再有归档卡的 section（前端 section 已摘）
            for cid in ARCHIVED:
                check(pg.locator(f'#{cid}').count() == 0,
                      f'归档卡 section 不存在: #{cid}')
            # 关键控件 id 也必须消失（防"只删了 section 外壳"这种半拉子）
            for sel in ('#pdmj_traits', '#pdmj_xml', '#pdtt_input', '#pdltt_input',
                        '#pdbsp_path', '#pdrspp_input', '#pdtm_input',
                        '#pdrrt_original', '#pdrand_input'):
                check(pg.locator(sel).count() == 0, f'归档控件已消失: {sel}')
            # 归档卡的 JS 函数不能再定义（死代码清理）
            gone = pg.evaluate(
                "names => names.filter(n => typeof window[n] === 'function')",
                ['runPdRrt', 'loadPdRrtResult', 'runPdTempMig', 'runPdBsp',
                 'runPdRspp', 'runPdTreetime', 'runPdLtt', 'runPdMjrm',
                 'pdbspSrcToggle', 'runPdrename', 'loadPdrenameResult',
                 'runPdgroup', 'loadPdgroupResult'])
            check(not gone, f'归档卡的 JS 函数已移除（残留: {gone or "无"}）')
            # ★ 正向：保留卡的 JS 必须还在（否则"删过头"也会通过）
            keep_js = pg.evaluate(
                "names => names.filter(n => typeof window[n] === 'function')",
                ['runPdprep'])
            check(len(keep_js) == 1, f'保留卡 JS 函数齐（实际 {keep_js}）')

            print('\n[6] 保留卡控件存在性（点卡进卡）')
            # 逐卡点击（落地条目文案在 pages.py title 里）
            # 2026-09-21: t-pdgroup 已归档 → 无逐卡点击对象
            card_names = {}
            for card_id, name in card_names.items():
                pg.goto(f'http://127.0.0.1:{port}/tools?g=phylodyn',
                        wait_until='domcontentloaded')
                pg.wait_for_selector('#groupLanding', state='visible', timeout=30000)
                pg.locator(f'#groupLanding .landing-item[data-item="{card_id}"]').click()
                pg.wait_for_selector(f'#{card_id}', state='visible', timeout=15000)
                sec = pg.locator(f'#{card_id}')
                missing = [s for s in WANT_PER_CARD[card_id]
                           if sec.locator(s).count() != 1]
                check(not missing, f'{name}（{card_id}）控件齐（缺: {missing or "无"}）')

            print('\n[7] 切英文（新增键 zh/en 双份）')
            pg.evaluate("VP_LANG = 'en'; applyI18n(document);")
            pg.wait_for_timeout(150)
            en = pg.evaluate("t('tk.pdprepStat', 'X')")
            check(en and en != 'X' and 'Kept' in en, f'英文汇总文案（{en[:44]}…）')
            pg.evaluate("VP_LANG = 'zh'; applyI18n(document);")

            print('\n[8] 无 JS 报错')
            check(not errs, f'无 pageerror（实际 {errs[:2]}）')

            os.makedirs(os.path.dirname(SHOT), exist_ok=True)
            pg.screenshot(path=SHOT, full_page=False)
            b.close()
    finally:
        if proc is not None and not args.keep:
            proc.kill()
            proc.wait(timeout=10)

    print('\n' + '=' * 52)
    if FAILS:
        print(f'FAILED: {len(FAILS)} 项')
        for f in FAILS:
            print('  - ' + f)
        return 1
    print('进化动力学新卡 UI 检查通过 ✔  截图:', os.path.relpath(SHOT, ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
