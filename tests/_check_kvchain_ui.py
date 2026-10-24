# -*- coding: utf-8 -*-
"""真浏览器验证「⚡一键分析（全流程）」卡：输入 / 参数 / 输出都在本卡内。

背景（2026-09-11 用户截图报障）：「一键分析（全流程）」页面没有输入、没有输出、
没有参数——参数读的是另一个二级模块 t-kvsuite 卡的 kv_* 字段（点开本页看不到也
改不了），本卡也没有 .toolrun 输出容器，跑起来日志/产物无处落。

检查项：
  1. 本卡可见输入 >= 15 个，样品/参考/注释/鉴定库/引擎/阈值/开关/线程数都在；
  2. 本卡有 .tool-out 输出区 + #toolrun-kvchain 容器；
  3. 结果面板容器齐全（徽标/流程各段/汇总目录/鉴定表/过滤表/变异图/演化图/深度图）；
  4. 「✨ 示例」填本卡直接 reads（#kvc_reads_r1/r2 → 内置示例 reads 绝对路径），
     并清空样品框——不再填样品名 EXAMPLE（那个样品无人创建，点了必报
     「路径不存在: run\\results\\EXAMPLE」，2026-10-13 用户换测试数据后复现）；
  5. 「⤓ 取该卡参数」把 t-kvsuite 卡当前表单值搬过来；
  6. 未选样品点运行时给明确提示，不静默；
  7. 内联结果链路：用 page.route 桩住 chain_result / kvsuite_result，
     验证「汇总清单 → 子运行 → 内联 kvc* 表格/图」渲染正确，
     且不写入 t-kvsuite 卡的 kv* 容器（两卡互不覆盖）。

用法：python tests/_check_kvchain_ui.py            （自起临时实例，随机空闲端口）
      python tests/_check_kvchain_ui.py --port 8765（用已在运行的平台）
"""
import argparse
import io
import os
import socket
import subprocess
import sys
import tempfile
import time

# playwright 属开发依赖（requirements-dev.txt），未装时本测试**跳过**而非失败：
# 它是「真浏览器验证」的加强项，不该让没装开发依赖的环境整批测试变红。
# 与其他浏览器探针（_it_download_ui / _it_async_ui 等）的约定保持一致。
try:
    from playwright.sync_api import sync_playwright
except ImportError:                     # pragma: no cover - 取决于本机环境
    sync_playwright = None

# 控制台编码兜底：Windows 默认代码页是 GBK，✔/✘ 等字符会让 print 抛
# UnicodeEncodeError。只改错误处理为 replace（编码不动，中文照常可读）。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SHOT = os.path.join(ROOT, 'run', '_check_kvchain_ui.png')

# 期望在本卡内可见的参数控件（卡内定位，避免误命中另一张卡的同名字段）
WANT_CTRLS = [
    ('#kvc_samples', '样品'),
    ('#kvc_reads_r1', '直填测序数据 R1'),
    ('#kvc_reads_r2', '直填测序数据 R2'),
    ('#kvc_min_cov', '覆盖率下限'),
    ('#kvc_min_depth', '平均深度下限'),
    ('#kvc_min_reads', 'reads 下限'),
    ('#kvc_min_poisson', '泊松置信度'),
    ('#kvc_variant_qual', '变异 QUAL 下限'),
    ('#kvc_min_freq', '次等位频率下限'),
    ('#kvc_aa_label_cutoff', '氨基酸标签 AF 下限'),
    ('#kvc_max_aa_labels', '标签数量上限'),
    ('#kvc_threads', '线程数'),
    ('#kvc_evo', '扩展变异分析'),
    ('#kvc_snpgenie', 'SNPGenie'),
    ('#kvc_all_variants', '画全部变异'),
]
# 期望存在的容器（输出区 + 结果面板）
WANT_BOXES = ['#toolrun-kvchain', '#kvChainResult', '#kvcSummary', '#kvcSteps',
              '#kvcEntries', '#kvcIdentifyTable', '#kvcFilteredTable',
              '#kvcConsensusTable', '#kvcVariantTable', '#kvcVariantPlots',
              '#kvcEvoSummary', '#kvcEvoPlots', '#kvcEvoTables', '#kvcPlots']

# 桩数据：_chain.json 形态（子运行 = kvsuite_ui_chk）与 kvsuite_result 形态
_CHAIN_JSON = """{
 "run": "kvchain_ui_chk",
 "out_dir": "D:\\\\tmp\\\\run\\\\tool_runs\\\\kvchain_ui_chk",
 "meta": {"chain": "kvchain", "run": "kvchain_ui_chk", "done": true, "steps": [
   {"key": "kvsuite", "run": "kvsuite_ui_chk",
    "title": "已知病毒识别与定量（全流程）",
    "summary": {"stage": "all", "engine": "minibwa", "n_confirmed": 2}}]},
 "entries": [
   {"name": "step_kvsuite", "is_dir": true, "is_link": true, "size": null},
   {"name": "kvsuite__run.log", "is_dir": false, "is_link": true, "size": "1.2 KB"}]}"""

_KV_JSON = """{
 "run": "kvsuite_ui_chk",
 "summary": {"n_identify": 3, "n_filtered": 2, "n_discarded": 1, "n_consensus": 2,
             "n_variants": 5, "n_variant_plots": 0, "n_evo_plots": 0},
 "plots": [], "variant_plots": [], "evo_plots": [], "evo_manifest": {},
 "popgen_tables": {}, "snpgenie_products": {},
 "identify": [{"Sample": "EXAMPLE", "Accession": "NC_001", "Species": "TMV",
               "Molecule_Type2": "ssRNA(+)", "Length": "6395",
               "Coverage(%)": "98.2", "MeanDepth": "120.5",
               "Uniq_Reads": "3210", "EM_Reads": "3200", "taxid": "12242"}],
 "filtered": [{"Sample": "EXAMPLE", "Accession": "NC_001", "Species": "TMV",
               "Segment": "", "is_segmented": "False", "Length": "6395",
               "Coverage(%)": "98.2", "MeanDepth": "120.5",
               "Uniq_Reads": "3210", "EM_Reads": "3200"}],
 "discarded": [],
 "consensus_qc": [{"Sample": "EXAMPLE", "Reference": "NC_001", "Feature": "ORF1",
                   "Gene": "RdRp", "Length": "3000", "N_ratio(%)": "0.1",
                   "Stops": "0", "QC": "PASS"}],
 "variants": [{"genome": "NC_001", "POS": "123", "REF": "A", "ALT": "G",
               "Gene_Name": "RdRp", "Annotation": "missense",
               "Impact": "MODERATE", "HGVS_c": "c.1A>G", "HGVS_p": "p.K1E",
               "AA": "K/E", "CDS": "cds1", "QUAL": "30.2", "DP": "100",
               "AF": "0.42"}]}"""

FAILS = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def free_port(start=8791):
    for p in range(start, start + 40):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(('127.0.0.1', p))
                return p
            except OSError:
                continue
    raise RuntimeError('没有可用端口')


def start_server(port):
    """起一个临时平台实例（服务日志落文件：PIPE 写满会让服务阻塞）。

    注意 `-c` 代码里那行 sys.path 注入：项目自带的绿色版 Python 目录下有
    `python312._pth`，存在该文件时 Python 进入 isolated 模式，`-c` 的
    sys.path **不含当前工作目录**（只有 python.exe 所在目录），于是
    `import app` 直接 ModuleNotFoundError、服务永远起不来。用脚本文件跑
    时 sys.path[0] 是脚本所在目录，所以只有这种 `-c` 起子进程的写法会中招。
    """
    logpath = os.path.join(tempfile.gettempdir(), f'vp_kvchain_ui_{port}.log')
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

    proc, logpath = None, ''
    port = args.port or free_port()
    if not args.port:
        proc, logpath = start_server(port)
        if proc is None:
            return 2
        print(f'临时实例: http://127.0.0.1:{port}  (服务日志 {logpath})')
    url = f'http://127.0.0.1:{port}/tools?g=kvsuite#t-kvchain'

    try:
        with sync_playwright() as pw:
            b = pw.chromium.launch(headless=True)
            pg = b.new_page(viewport={'width': 1600, 'height': 1100})
            errs = []
            pg.on('pageerror', lambda e: errs.append(str(e)))

            pg.goto(url, wait_until='domcontentloaded')
            # 首次访问的「示例结果」浮层会挡住卡片按钮；标记已看过再刷新
            pg.evaluate("try{localStorage.setItem('vp_example_seen','1')}catch(e){}")
            pg.reload(wait_until='domcontentloaded')
            pg.wait_for_selector('#kvc_samples', state='visible', timeout=30000)
            pg.evaluate("document.querySelectorAll('.dlgmask').forEach(m => m.remove())")
            sec = pg.locator('#t-kvchain')

            print('\n[1] 本卡输入与参数')
            n_in = sec.locator('input:visible').count()
            check(n_in >= 13, f'本卡可见输入 {n_in} 个（>=13）')
            for sel, what in WANT_CTRLS:
                check(sec.locator(sel).is_visible(), f'参数可见: {what} {sel}')

            print('\n[2] 输出区与结果面板容器')
            check(sec.locator('.tool-out').count() >= 1, '本卡有输出区 .tool-out')
            for sel in WANT_BOXES:
                check(sec.locator(sel).count() == 1, f'容器存在: {sel}')

            print('\n[3] 示例填入落到本卡')
            ex = sec.locator('.vp-ex-btn', has_text='示例').first
            check(ex.count() > 0, '标题栏有「✨ 示例」按钮')
            pg.evaluate("document.querySelector('#kvc_samples').value='STALE_S'")
            pg.evaluate("document.querySelector('#kvc_reads_r1').value=''")
            pg.evaluate("document.querySelector('#kvc_reads_r2').value=''")
            ex.click()
            pg.wait_for_timeout(800)   # 等异步取 /api/example_paths 绝对路径
            v = pg.eval_on_selector('#kvc_samples', 'e => e.value')
            check(v == '', f'「示例」清空样品框（不依赖 run/results/<样品> 存在）'
                           f'（实际 {v!r}）')
            r1 = pg.eval_on_selector('#kvc_reads_r1', 'e => e.value')
            r2 = pg.eval_on_selector('#kvc_reads_r2', 'e => e.value')
            check(r1.endswith('example_R1.fastq.gz')
                  and r2.endswith('example_R2.fastq.gz'),
                  f'「示例」填本卡直接 reads R1/R2（实际 {r1!r} / {r2!r}）')

            print('\n[4] 「取该卡参数」从 t-kvsuite 卡搬运')
            pg.evaluate("document.querySelector('#kv_samples').value='PULLED_A,B'")
            pg.evaluate("document.querySelector('#kv_min_cov').value='33'")
            pg.evaluate("document.querySelector('#kv_evo').checked=false")
            sec.locator('button', has_text='取该卡参数').first.click()
            pg.wait_for_timeout(300)
            check(pg.eval_on_selector('#kvc_samples', 'e => e.value') == 'PULLED_A,B',
                  '搬来样品名')
            check(pg.eval_on_selector('#kvc_min_cov', 'e => e.value') == '33',
                  '搬来覆盖率阈值')
            check(pg.eval_on_selector('#kvc_evo', 'e => e.checked') is False,
                  '搬来开关状态')

            print('\n[5] 未选样品点运行的提示')
            dialogs = []
            pg.on('dialog', lambda d: (dialogs.append(d.message), d.dismiss()))
            pg.evaluate("document.querySelector('#kvc_samples').value=''")
            # [3] 的示例填充已把直接 reads 填上；「无输入」场景须连 reads 一起清
            pg.evaluate("document.querySelector('#kvc_reads_r1').value=''")
            pg.evaluate("document.querySelector('#kvc_reads_r2').value=''")
            sec.locator('button', has_text='一键运行').first.click()
            pg.wait_for_timeout(800)
            check(bool(dialogs) and '样品' in dialogs[0],
                  f'有明确提示（实际 {dialogs[:1]}）')

            print('\n[6] 内联结果渲染（桩数据：chain_result + kvsuite_result）')
            # 用桩响应打通「汇总清单 → 子运行 → 内联 kvc* 容器」这条链路，
            # 不依赖本机真有 kvchain 运行记录。
            pg.route('**/api/tool/chain_result*', lambda route: route.fulfill(
                status=200, content_type='application/json', body=_CHAIN_JSON))
            pg.route('**/api/tool/kvsuite_result*', lambda route: route.fulfill(
                status=200, content_type='application/json', body=_KV_JSON))
            pg.evaluate("loadChainResult('kvchain', 'kvchain_ui_chk')")
            pg.wait_for_timeout(1200)
            check(pg.eval_on_selector('#kvChainResult',
                                      "e => getComputedStyle(e).display") != 'none',
                  '汇总面板已展开')
            check('kvchain_ui_chk' in pg.inner_text('#kvcTitle'), '标题含汇总运行名')
            check('kvsuite_ui_chk' in pg.inner_text('#kvcSubRun'), '标出定量子运行名')
            check('已知病毒识别与定量' in pg.inner_text('#kvcSteps'), '流程各段表已渲染')
            check(pg.locator('#kvcEntries tr').count() >= 2, '汇总目录清单已渲染')
            check('确认存在 2' in pg.inner_text('#kvcSummary'), '计量徽标来自子运行结果')
            check(pg.locator('#kvcIdentifyTable tbody tr').count() == 1, '鉴定表已渲染')
            check(pg.locator('#kvcFilteredTable tbody tr').count() == 1, '过滤表已渲染')
            check(pg.locator('#kvcConsensusTable tbody tr').count() == 1, '共识 QC 表已渲染')
            check(pg.locator('#kvcVariantTable tbody tr').count() == 1, '变异注释表已渲染')
            # 内联渲染不得污染另一张卡（t-kvsuite 的结果容器仍为空）
            check(pg.eval_on_selector('#kvIdentifyTable',
                                      'e => e.innerHTML.trim()') == '',
                  '内联渲染未写入 t-kvsuite 卡容器')

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
    print('一键全流程卡 UI 检查通过 ✔  截图:', os.path.relpath(SHOT, ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
