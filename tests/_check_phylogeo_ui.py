# -*- coding: utf-8 -*-
"""真浏览器验证主平台「系统地理分析」卡（2026-09-16 裁剪版移植）。

背景：09-15 拆分把 `compare` / `phylodyn` 两组从一级导航摘掉后，`#t-rdp` /
`#t-rtt` / `#t-phylogeo` 三张卡**路由还在、卡片还在，但从界面上到不了**
（组落地路由要求组存在于 `NAV_GROUPS`）。2026-09-16 移植了自研系统地理能力
（haversine 距离三分类 / MOTP 时间分箱 / 权重分带 / 迁移 GIF）后恢复了最小
分组与一级导航项，本测试锁住「真的能点到、参数真的在、跑完真的出面板」。

检查项：
  [1] 一级导航有「进化动力学分析」→ /tools?g=phylodyn 组落地列出 9 张卡
      （不是空白页，也不是掉到别的组）；
  [2] 点「系统地理分析」卡片进卡：本卡输入/参数/输出容器齐全；
  [3] 静态方法学提醒与新增参数标签已 i18n 成中文（不是裸键名）；
  [4] 语言开关切英文后同一批文案变英文（zh/en 双份都在）；
  [5] 真跑一次（8 条合成比对 + 元数据 + 坐标表，勾 RRT/RSSP/MOTP）：
      结果面板出现，三分类表 / 权重分带表 / GIF 区都渲染出内容，
      且分带口径写明「非贝叶斯」（不许被叫成贝叶斯后验）；
  [6] GIF 区在「有 MOTP 但坐标不足」时给明确指引（纯函数直接调）；
  [6b] 迁移弧线地图：断言 scattergeo 图元与弧线真的画出来（读 plotly
       `_fullData`，不是断言 DOM 存在）；负控＝可用坐标区划 <2 时必须
       显式列出缺口，不许画成空白（空白会被读成"没有迁移"，
       而真相是"没有坐标"）；
  [7] 无 pageerror。

用法：python tests/_check_phylogeo_ui.py            （自起临时实例，随机空闲端口）
      python tests/_check_phylogeo_ui.py --port 8765（用已在运行的平台）
"""
import argparse
import io
import json
import os
import socket
import subprocess
import sys
import tempfile
import time

# playwright 属开发依赖（requirements-dev.txt），未装时本测试**跳过**而非失败：
# 它是「真浏览器验证」的加强项，不该让没装开发依赖的环境整批测试变红。
try:
    from playwright.sync_api import sync_playwright
except ImportError:                     # pragma: no cover - 取决于本机环境
    sync_playwright = None

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SHOT = os.path.join(ROOT, 'run', '_check_phylogeo_ui.png')
FIX = os.path.join(ROOT, 'run', '_check_phylogeo_ui_fix')

# 卡内定位的控件（pg_* 全部在本卡内，不会误命中别处）
WANT_CTRLS = [
    ('#pg_input', '比对 FASTA'),
    ('#pg_meta', '元数据 CSV'),
    ('#pg_trait', '区划列名'),
    ('#pg_method', '建树方法'),
    ('#pg_rrt', 'RRT 开关'),
    ('#pg_rrt_perm', '置换次数'),
    ('#pg_rssp', 'RSSP 开关'),
    ('#pg_rssp_n', 'Bootstrap 复本数'),
    ('#pg_motp', 'MOTP 开关'),
    ('#pg_motp_bin', '分箱宽度'),
    ('#pg_dated_tree', 'LSD2 定年树'),
    ('#pg_date_trait', '日期列名'),
    ('#pg_coords', '区域坐标表'),
    ('#pg_threads', '线程数'),
    # 轨道 B（BEAST 打包 / 产物导入，2026-10-08 接线）
    ('#pg_beast_log', 'BEAST 后验日志'),
    ('#pg_beast_trees', 'BEAST 树'),
    ('#pg_beast_burnin', 'BEAST burn-in'),
    ('#pg_beast_min_prob', 'BEAST 后验阈值'),
    ('#pg_beast_max_trees', 'BEAST 最多树数'),
    ('#pg_export_tag', 'BEAST 作业名'),
]
WANT_BOXES = ['#toolrun-phylogeo', '#pgResult', '#pgTitle', '#pgSummary',
              '#pgGeo', '#pgMatrix', '#pgTransClasses', '#pgBands', '#pgBeast',
              '#pgGif', '#pgExport']
# 必须被 i18n 覆盖的静态文案（裸键名 = 语言开关失灵）
WANT_KEYS = ['tk.pgWarnStatic', 'tk.pgRrtLbl', 'tk.pgRsspLbl', 'tk.pgMotpLbl',
             'tk.pgDatedTreeLbl', 'tk.pgDateTraitLbl', 'tk.pgCoordsLbl',
             'tk.pgCoordsHint', 'tk.pgBeastLbl', 'tk.pgExportLbl',
             'tk.pgExportBtn', 'tk.pgExportHint']

FAILS = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def free_port(start=8821):
    for p in range(start, start + 40):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(('127.0.0.1', p))
                return p
            except OSError:
                continue
    raise RuntimeError('没有可用端口')


def start_server(port):
    """起临时实例（日志落文件；`-c` 启动必须显式注入 ROOT 到 sys.path，
    绿色版 Python 的 python312._pth 会让 `-c` 的 sys.path 不含 cwd）。"""
    logpath = os.path.join(tempfile.gettempdir(), f'vp_phylogeo_ui_{port}.log')
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


def write_fixture():
    """8 条 60nt 两个簇（Asia / Europe）+ 元数据 + 坐标表；跑得秒级。"""
    os.makedirs(FIX, exist_ok=True)
    base = 'ACGT' * 15
    nxt = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
    muts = {'s1': {3: 1}, 's2': {3: 1, 7: 1}, 's3': {3: 1, 12: 1},
            's4': {3: 1, 20: 1}, 's5': {30: 1}, 's6': {30: 1, 35: 1},
            's7': {30: 1, 41: 1}, 's8': {30: 1, 50: 1}}
    region = {f's{i}': ('Asia' if i <= 4 else 'Europe') for i in range(1, 9)}
    lines = []
    for nm in sorted(muts):
        s = list(base)
        for pos in muts[nm]:
            s[pos] = nxt[s[pos]]
        lines += ['>' + nm, ''.join(s)]
    aln = os.path.join(FIX, 'geo.fasta')
    with io.open(aln, 'w', encoding='utf-8', newline='') as f:
        f.write('\n'.join(lines) + '\n')
    meta = os.path.join(FIX, 'geo_meta.csv')
    rows = ['seq_id,region,year,latitude,longitude']
    for i, nm in enumerate(sorted(muts)):
        rows.append(f'{nm},{region[nm]},{2000 + i},{35.0 if i < 4 else 51.0},'
                    f'{100.0 if i < 4 else 12.0 + i}')
    with io.open(meta, 'w', encoding='utf-8', newline='') as f:
        f.write('\n'.join(rows) + '\n')
    coords = os.path.join(FIX, 'geo_coords.tsv')
    with io.open(coords, 'w', encoding='utf-8', newline='') as f:
        f.write('Region\tLatitude\tLongitude\nAsia\t35\t105\nEurope\t51\t12\n')
    return aln, meta, coords


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=0,
                    help='用已在运行的平台端口（默认自起临时实例）')
    ap.add_argument('--keep', action='store_true', help='保留临时实例（排查用）')
    args = ap.parse_args()

    if sync_playwright is None:
        print('[SKIP] 未安装 playwright，跳过真浏览器验证')
        print('       安装: python -m pip install -r requirements-dev.txt')
        return 0

    aln, meta, coords = write_fixture()
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
            # 首次访问的「示例结果」浮层会挡住卡片；标记已看过再刷新
            pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")
            pg.reload(wait_until='domcontentloaded')
            pg.wait_for_selector('#groupLanding', timeout=30000)
            pg.evaluate("document.querySelectorAll('.dlgmask').forEach(m => m.remove())")

            print('\n[1] 一级导航 + 组落地')
            nav_ok = pg.evaluate(
                "Array.from(document.querySelectorAll('.navlinks a[data-group]'))"
                ".some(a => a.dataset.group === 'phylodyn'"
                " && a.getAttribute('href').includes('g=phylodyn'))")
            check(nav_ok, '一级导航有「进化动力学分析」且指向 ?g=phylodyn')
            cards = pg.eval_on_selector_all(
                '#groupLanding .landing-item',
                'els => els.map(e => e.textContent.trim())')
            # 2026-09-21 起组内 5 张卡：3 张原有（RDP5 / 时间信号 / 系统地理）
            # + 1 张数据接入 + 1 张「⚡ 一键分析（定年 → 地理迁移 → 汇总溯源）」。
            # （ID 重命名 / 分组 两卡归档：重命名接进全链可选前置，分组归导入收集）
            # 其余 8 张（rand/rrt/tempmig/bsp/rspp/treetime/ltt/mjrm）已于 09-18 归档。
            check(len(cards) == 5,
                  f'组落地列出 5 张卡（实际 {len(cards)}：'
                  f'{[c.split(chr(10))[0] for c in cards][:4]}…）')
            check(any('系统地理' in c for c in cards), '落地含「系统地理分析」卡')
            check(pg.eval_on_selector('#groupLanding',
                                      "e => getComputedStyle(e).display") != 'none',
                  '显示的是组落地（不是空白页）')

            print('\n[2] 点卡片进卡：参数 / 输出容器齐全')
            pg.locator('#groupLanding .landing-item', has_text='系统地理').click()
            pg.wait_for_selector('#t-phylogeo', state='visible', timeout=15000)
            sec = pg.locator('#t-phylogeo')
            check(pg.url.endswith('#t-phylogeo'), f'hash 落到 #t-phylogeo（{pg.url}）')
            for sel, what in WANT_CTRLS:
                check(sec.locator(sel).is_visible(), f'参数可见: {what} {sel}')
            for sel in WANT_BOXES:
                check(sec.locator(sel).count() == 1, f'容器存在: {sel}')
            check(sec.locator('.tool-out').count() >= 1, '本卡有输出区 .tool-out')
            check(pg.locator('#pg_rrt').is_checked() is False
                  and pg.locator('#pg_motp').is_checked() is False,
                  '统计开关默认关（不勾不拖慢常规运行）')

            print('\n[3] 静态文案已 i18n（中文可见，不是裸键名）')
            for key in WANT_KEYS:
                txt = pg.eval_on_selector(
                    f'[data-i18n="{key}"]', 'e => e.textContent.trim()')
                check(bool(txt) and not txt.startswith('tk.'),
                      f'{key} → {txt[:34]}…')

            print('\n[4] 语言开关切英文（同一批键有英文）')
            # 直接改 VP_LANG + 重刷文案：走 setLang() 会 POST /api/settings 改平台
            # 语言并整页 reload —— 用 --port 8765 打真实例时那是用户设置的副作用。
            zh_hint = pg.eval_on_selector('[data-i18n="tk.pgCoordsHint"]',
                                          'e => e.textContent.trim()')
            pg.evaluate("VP_LANG = 'en'; applyI18n(document);")
            pg.wait_for_timeout(150)
            en_hint = pg.eval_on_selector('[data-i18n="tk.pgCoordsHint"]',
                                          'e => e.textContent.trim()')
            check(en_hint and en_hint != zh_hint and not en_hint.startswith('tk.'),
                  f'切英文生效（{en_hint[:40]}…）')
            en_title = pg.evaluate("t('tk.pgRrtLbl', 'X')")
            check(en_title and en_title != 'X' and 'tk.' not in en_title,
                  f'JS 文案同源取到英文（{en_title}）')
            pg.evaluate("VP_LANG = 'zh'; applyI18n(document);")
            pg.wait_for_timeout(150)

            print('\n[5] 真跑一次（RRT/RSSP/MOTP 全开）')
            pg.fill('#pg_input', aln)
            pg.fill('#pg_meta', meta)
            pg.fill('#pg_coords', coords)
            pg.fill('#pg_rrt_perm', '60')
            pg.fill('#pg_rssp_n', '5')
            pg.check('#pg_rrt')
            pg.check('#pg_rssp')
            pg.check('#pg_motp')
            sec.locator('button', has_text='运行系统地理分析').first.click()
            ok = False
            for _ in range(120):            # 最多等 120s（真建树 + 60 次置换）
                if pg.eval_on_selector('#pgResult',
                                       "e => getComputedStyle(e).display") != 'none':
                    ok = True
                    break
                pg.wait_for_timeout(1000)
            if not ok:
                tlog = pg.eval_on_selector('#toolrun-phylogeo',
                                           'e => e.innerText').strip()
                print('    任务日志尾部:', tlog[-400:] or '(空)')
            check(ok, '结果面板已出现（页面轮询到任务完成）')
            if ok:
                summary = pg.inner_text('#pgSummary')
                check('迁移事件' in summary, '汇总行含迁移事件数')
                tc = pg.inner_text('#pgTransClasses')
                check('迁移事件分类' in tc, f'三分类表已渲染（{tc[:40]}…）')
                for lbl in ('Direct', 'Indirect', 'Distant'):
                    check(lbl in tc, f'三分类表含 {lbl}')
                check(pg.locator('#pgTransClasses table tbody tr').count() >= 1,
                      '三分类表有数据行')
                wb = pg.inner_text('#pgBands')
                check('后验权重分带' in wb and '非贝叶斯' in wb,
                      f'分带表已渲染且口径写明非贝叶斯（{wb[:40]}…）')
                check(pg.locator('#pgBands table tbody tr').count() == 4,
                      '分带表四档齐全')
                gif = pg.inner_text('#pgGif')
                check('motp' in pg.inner_text('#pgSummary').lower()
                      or 'MOTP' in summary, '汇总含 MOTP 段')
                check(pg.locator('#pgGifBtn').count() == 1,
                      f'GIF 导出按钮已出（GIF 区: {gif[:40] or "按钮"}…）')
                check('贝叶斯' in wb and 'BEAST' not in wb and 'BSSVS' not in wb,
                      '分带措辞不冒充 BEAST/BSSVS')

                # 迁移弧线地图（2026-09-18 新增）：验「图真的画出来了」——
                # 不是断言 DOM 元素存在。plotly 6.x 的 to_plotly_json() 会把
                # numpy 编成 {dtype,bdata}，所以读 el._fullData 而不是 el.data。
                check(pg.locator('#pgGeoPlot').count() == 1, '弧线地图容器已出')
                gd = pg.eval_on_selector('#pgGeoPlot', """e => {
                  const t = (e._fullData || []);
                  const geos = t.filter(x => x.type === 'scattergeo');
                  return {n: t.length, geo: geos.length,
                          lines: geos.filter(x => x.mode === 'lines').length,
                          note: (document.getElementById('pgGeoNote') || {}).innerText || ''};
                }""")
                check(gd['geo'] >= 2,
                      f"真画出 scattergeo 图元（{gd['geo']} 个 / 共 {gd['n']} trace）")
                check(gd['lines'] >= 1,
                      f"至少一条迁移弧线（{gd['lines']} 条）")
                check('弧线' in gd['note'], '说明行写明弧线口径')
                check('坐标出处' in gd['note'],
                      '说明行标出坐标出处（不许把近似质心当准坐标）')
                # 缺坐标必须显式交代：本 fixture 的坐标表覆盖全部区划，
                # 所以这里应当**没有**缺项；一旦有，说明行必须写出「未落点」。
                if '缺坐标未落点' in gd['note']:
                    check('未落点' in gd['note'],
                          '有缺坐标时必须显式列出，不许静默少画')
                else:
                    check('已画弧线' in gd['note'], '说明行给出已画弧线条数')

            print('\n[6] GIF 区：有 MOTP 但坐标不足时给指引（纯函数）')
            msg = pg.evaluate(
                "pgGifHtml({motp: {bins: [1], rows: [1]}, coords: {A: [0, 0]},"
                " coord_missing: ['B']}, 'x')")
            check('不足 2 个' in msg and 'B' in msg,
                  f'坐标不足时显式列出缺口（{msg[:60]}…）')
            msg2 = pg.evaluate("pgGifHtml({motp: {bins: [], rows: []}}, 'x')")
            check('MOTP' in msg2, f'没勾 MOTP 时给「先勾 MOTP」指引（{msg2[:40]}…）')

            print('\n[6b] 弧线地图负控：坐标不足 2 个时必须显式列出缺口（纯函数）')
            # 负控的意义：面板在「画不出来」时不能变成一块空白 —— 空白会被
            # 读成"没有迁移"，而真相是"没有坐标"。所以缺口必须写出来。
            msg = pg.evaluate("""() => {
              const d = document.createElement('div');
              d.id = 'wbGeoNeg';
              document.body.appendChild(d);
              drawPgGeo('wbGeoNeg', {coords: {A: [0, 0]},
                                     coord_missing: ['B', 'C']}, null);
              const s = d.innerText;
              d.remove();
              return s;
            }""")
            check('不足 2 个' in msg and 'B' in msg and 'C' in msg,
                  f'坐标不足时列出缺口（{msg[:60]}…）')

            print('\n[7] 页面无 JS 报错')
            check(not errs, f'无 pageerror（实际 {errs[:2]}）')

            os.makedirs(os.path.dirname(SHOT), exist_ok=True)
            pg.screenshot(path=SHOT, full_page=True)
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
    print('系统地理卡 UI 检查通过 ✔  截图:', os.path.relpath(SHOT, ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
