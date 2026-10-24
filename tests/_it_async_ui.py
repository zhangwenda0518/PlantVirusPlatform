# -*- coding: utf-8 -*-
"""四个「同步阻塞请求 → 后台任务」改造的**真实浏览器**验证
（Selenium + 本机 Edge/Chrome，headless）。

被验的是阶段一第 ② 项的**前端**部分：请求已改成异步，但页面必须真的
能跑起轮询、把中间进度显示出来、并在结束时把结果渲染到原来的位置。

  1. /submit  「🚀 生成 .sqn」     —— 原最坏 45 分钟阻塞 HTTP，无进度无取消
  2. /submit  「🌐 物种名校验」    —— 原每物种一次 esearch（每次最长 120s）串行
  3. /meta    「📝 AI 总结」      —— 原单次调用客户端超时 180s
  4. /cds-export「导出」          —— 原全量重解析集合 + 现场翻译

外部依赖（suvtk / DeepSeek / NCBI）全部在本进程内用桩替换，离线可跑；
临时实例起在 127.0.0.1 的独立端口，不碰用户的 8989。
无浏览器/无 selenium 时自动 SKIP 并返回 0。
用法: python tests/_it_async_ui.py [--port 8798] [--keep] [--headed]
"""
import argparse
import json
import os
import shutil
import socket
import sys
import threading
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
os.environ['VP_NO_RECOVER'] = '1'

FAIL = []
_CLEAN = []


def check(label, ok, extra=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {label}' + (f'  {extra}' if extra else ''),
          flush=True)
    if not ok:
        FAIL.append(label)


def skip(label, why):
    print(f'  [SKIP] {label}  {why}', flush=True)


def free_port(preferred):
    for p in (preferred, 0):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', p))
                return s.getsockname()[1]
        except OSError:
            continue
    raise RuntimeError('无可用端口')


def _run_once():
    global _CLEAN
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=8798)
    ap.add_argument('--keep', action='store_true')
    ap.add_argument('--headed', action='store_true')
    a = ap.parse_args()

    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.edge.options import Options as EdgeOptions
    except ImportError as e:
        skip('selenium', f'未安装（{e}）')
        return 0

    import pandas as pd
    import app as appmod
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.ncbi_submit import store as sstore
    from Virus_Platform_Core import suvtk_submit as ss
    from Virus_Platform_Core import ncbi_download
    from Virus_Platform_Core.web import meta as vpm
    from Virus_Platform_Core import cds_export
    from Virus_Platform_Core.gb_collection import gb_collection_dir

    PROJ = 'uitest_async_proj'
    RUN_NAME = 'orf_uitest_async'
    COLL = 'uitest_async_cds'
    MPROJ = 'uitest_async_meta'
    tool_runs = DIRS.get('tool_runs')
    run_dir = os.path.join(tool_runs, RUN_NAME)
    _CLEAN += [sstore.table_dir(PROJ), run_dir, gb_collection_dir(COLL),
               vpm._meta_proj(MPROJ)]

    # ── 夹具 ────────────────────────────────────────────────────
    shutil.rmtree(sstore.table_dir(PROJ), ignore_errors=True)
    shutil.rmtree(run_dir, ignore_errors=True)
    shutil.rmtree(gb_collection_dir(COLL), ignore_errors=True)
    shutil.rmtree(vpm._meta_proj(MPROJ), ignore_errors=True)

    os.makedirs(os.path.join(run_dir, '03_assembly'), exist_ok=True)
    with open(os.path.join(run_dir, '03_assembly', 'viral_contigs.fasta'),
              'w', encoding='utf-8') as f:
        f.write('>contig1 taxid=11001 taxon=Test_virus\n' + 'ACGT' * 120 + '\n')
    sstore.create_table(PROJ, sample='demo')
    df = sstore.load_table(PROJ)
    df.loc[df.index[0], 'sequence_name'] = 'contig1'
    df.loc[df.index[0], 'organism'] = 'Test virus'
    if len(df) > 1:
        df.loc[df.index[1], 'organism'] = 'Notarealspecies xyz'
    df = df.loc[df.index[:2]]
    sstore.save_table(PROJ, df)
    sstore._write_link(PROJ, 'orf', {'run': RUN_NAME})

    cdir = gb_collection_dir(COLL)
    os.makedirs(cdir, exist_ok=True)
    with open(os.path.join(cdir, 'x.gb'), 'w', encoding='utf-8') as f:
        f.write("""LOCUS       ITASYNC                   60 bp    RNA     linear
ACCESSION   ITASYNC
ORGANISM  Test virus
  Viruses.
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="Test virus"
     CDS             1..60
                     /gene="cp"
                     /product="coat protein"
                     /translation="MKLT"
ORIGIN
        1 atggctaaat tgacttaata a
//
""")
    cds_export.build_cds_table(COLL)          # 预热缓存，导出会更快

    vpm._meta_store_df(MPROJ, 'core14', pd.DataFrame({
        'Run': ['SRR1', 'SRR2'], 'ScientificName': ['Test virus'] * 2,
        'Database': ['SRA', 'GSA'], 'FileSize_GB': ['1.0', '2.0'],
        'CollectionDate': ['2024-01-01', '2024-02-01']}))

    # ── 桩：外部依赖全部替换 ────────────────────────────────────
    CANNED = 'CANNED SUMMARY: public runs collected for this study.'
    _orig = (ss.build_sqn, ncbi_download.esearch, vpm._meta_deepseek_key,
             vpm._meta_ai_chat, vpm._meta_ai_client, vpm._meta_load_df)
    stages_seen = []

    def stub_build_sqn(*args, **kw):
        log = kw.get('log') or (lambda m: None)
        out_dir = kw.get('out_dir') or sstore.table_dir(PROJ)
        for kw_word in ('suvtk features', 'suvtk comments', 'table2asn'):
            stages_seen.append(kw_word)
            log(f'  [stub] {kw_word} ...')
            time.sleep(0.6)
        with open(os.path.join(out_dir, 'sqn.sqn'), 'wb') as f:
            f.write(b'\x00' * 64)
        with open(os.path.join(out_dir, 'sqn.val'), 'w', encoding='utf-8') as f:
            f.write('Warning: stub warning line\nInfo: stub info line\n')
        return {'sqn': os.path.join(out_dir, 'sqn.sqn'),
                'stats': {'error': 0, 'warning': 1, 'info': 1},
                'val_head': ['Warning: stub warning line']}

    tax_calls = []

    def stub_esearch(term, db='taxonomy', retmax=1):
        tax_calls.append(term)
        if 'Test virus' in term:
            return {'count': 1, 'ids': ['11001']}
        return {'count': 0, 'ids': []}

    ss.build_sqn = stub_build_sqn
    ncbi_download.esearch = stub_esearch
    vpm._meta_deepseek_key = lambda: 'stub-key'
    vpm._meta_ai_chat = lambda *x, **k: CANNED
    vpm._meta_ai_client = lambda *x, **k: object()

    _real_load = vpm._meta_load_df

    def stub_load_df(name, which):
        if name == MPROJ:
            return pd.DataFrame({
                'ScientificName': ['Test virus', 'Test virus'],
                'Database': ['SRA', 'GSA'], 'FileSize_GB': ['1.0', '2.0'],
                'CollectionDate': ['2024-01-01', '2024-02-01']})
        return _real_load(name, which)

    vpm._meta_load_df = stub_load_df

    # ── 起临时实例 ─────────────────────────────────────────────
    from werkzeug.serving import make_server
    port = free_port(a.port)
    srv = make_server('127.0.0.1', port, appmod.app, threaded=True)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f'http://127.0.0.1:{port}'
    print(f'临时实例: {base}  （用户的 8989 未被触碰）', flush=True)

    drv = None
    try:
        opts = EdgeOptions()
        if not a.headed:
            opts.add_argument('--headless=new')
        opts.add_argument('--no-sandbox')
        opts.add_argument('--disable-gpu')
        opts.add_argument('--window-size=1440,1000')

        def _new_drv(retries=3, backoff=1.5):
            # Selenium Manager 偶发拿不到驱动（网络抖动/正在自更新），重试
            _last = None
            for _att in range(retries):
                try:
                    return webdriver.Edge(options=opts)
                except Exception as e:                   # noqa: BLE001
                    _last = e
                    time.sleep(backoff)
            raise RuntimeError(f'Edge 连续 {retries} 次启动失败: '
                               f'{type(_last).__name__}: {str(_last)[:120]}')

        def _session():
            # 每个流程一个全新浏览器会话（2026-10-24）：同一会话连续加载
            # 4 个重页面时 headless Edge 会在流程 4 原生崩溃（WebDriver 会话
            # 中断），与 _check_pages_console 的 chromium 段错误同因——
            # 累积状态问题；分session后与「逐页独立进程」同口径。
            d = _new_drv()
            d.set_page_load_timeout(40)
            # 批量回归满载时元素等待偶发不够（2026-10-06 实测 40s 超时翻车、
            # 单跑 67s 全过）；批量运行器注入 VP_TEST_WAIT 放宽，单跑不变
            return d, WebDriverWait(d, int(os.environ.get('VP_TEST_WAIT', 40)))

        def _recycle(old):
            # 会话切换（2026-10-24 稳健版）：旧会话 quit **容错**——批量满载下
            # 第 4 次切换时 quit 会抛空消息 WebDriverException（会话/驱动已
            # 死，放弃即可）；新会话起不来时递增退避重试 5 次，仍败才冒泡。
            try:
                old.quit()
            except Exception:                            # noqa: BLE001
                pass
            _last = None
            for _i in range(5):
                try:
                    d = webdriver.Edge(options=opts)
                    d.set_page_load_timeout(40)
                    return d, WebDriverWait(
                        d, int(os.environ.get('VP_TEST_WAIT', 40)))
                except Exception as e:                   # noqa: BLE001
                    _last = e
                    time.sleep(1.0 + _i)
            raise RuntimeError(f'Edge 重建会话连续 5 次失败: '
                               f'{type(_last).__name__}: {str(_last)[:160]}')

        # 首个会话起不来 = 无浏览器环境，SKIP 不阻塞套件（与原口径一致）；
        # 后续会话起不来则让 WebDriverException 冒泡给 main() 重试。
        try:
            drv, wait = _session()
        except Exception as e:                           # noqa: BLE001
            skip('Edge 驱动', f'{type(e).__name__}: {str(e)[:80]}'
                              f'（已重试 3 次）')
            return 0
        from Virus_Platform_Core.web.tasks import tm

        def running_tasks():
            return [x for x in tm.list_all(logs=False)
                    if x['status'] == 'running']

        # ═══════ 1) .sqn 异步 ═══════
        print('\n== 1) /submit 「生成 .sqn」异步 ==', flush=True)
        drv.get(base + '/submit')
        wait.until(lambda d: d.find_elements(By.ID, 'sqnOut'))
        drv.execute_script("CUR = arguments[0];", PROJ)
        btn = drv.find_element(By.XPATH, "//button[contains(., '生成 .sqn')]")
        btn.click()
        # 轮询抓中间态：证明真的显示了进度/日志，而不是干等
        mids, seen_final = [], None
        t0 = time.time()
        while time.time() - t0 < 60:
            txt = drv.find_element(By.ID, 'sqnOut').text
            if '✅' in txt:
                seen_final = txt
                break
            if txt and txt not in mids:
                mids.append(txt)
            time.sleep(0.4)
        check('（1）.sqn 任务已在后台跑过（suvtk 三阶段都被调用）',
              stages_seen[:3] == ['suvtk features', 'suvtk comments',
                                  'table2asn'], str(stages_seen))
        check('（1）页面显示了中间进度/日志（不是干等）',
              any('进度' in m or 'suvtk' in m for m in mids),
              f'中间态 {len(mids)} 帧：{mids[:1]}')
        check('（1）结束时渲染出 .sqn 产物与校验摘要',
              bool(seen_final) and 'Error 0' in seen_final
              and 'Warning 1' in seen_final,
              f'最终文本={str(seen_final)[:150]!r}')
        check('（1）任务确实进了任务中心',
              any(x['name'].startswith('提交 .sqn 生成') or '.sqn' in x['name']
                  for x in tm.list_all(logs=False)),
              str([x['name'] for x in tm.list_all(logs=False)][:4]))

        # ═══════ 2) 物种名校验异步 ═══════
        drv, wait = _recycle(drv)        # 换全新会话，规避累积状态崩溃
        print('\n== 2) /submit 「物种名校验」异步 ==', flush=True)
        # 新会话要重新进入 /submit 并恢复页面级选中状态 CUR（流程 1 的
        # 会话里由 execute_script 设置，换会话即失效）
        drv.get(base + '/submit')
        wait.until(lambda d: d.find_elements(
            By.XPATH, "//button[contains(., '物种名校验')]"))
        drv.execute_script("CUR = arguments[0];", PROJ)
        drv.find_element(By.XPATH, "//button[contains(., '物种名校验')]").click()
        t0 = time.time()
        val_txt = ''
        while time.time() - t0 < 60:
            val_txt = drv.find_element(By.ID, 'valOut').text
            if '物种名校验（' in val_txt:
                break
            time.sleep(0.4)
        check('（2）校验结果渲染出来（读回 taxonomy_check.json）',
              '物种名校验（' in val_txt, f'{val_txt[:150]!r}')
        check('（2）✓ 已收录物种带 taxid',
              '✓' in val_txt and '11001' in val_txt, f'{val_txt[:150]!r}')
        check('（2）✗ 未收录物种被标出',
              '✗' in val_txt or '未收录' in val_txt, f'{val_txt[:150]!r}')
        check('（2）服务端确实逐个查过（stub 被调用）', len(tax_calls) >= 1,
              str(tax_calls))
        check('（2）结果文件已落项目目录',
              os.path.isfile(os.path.join(sstore.table_dir(PROJ),
                                          'taxonomy_check.json')))

        # ═══════ 3) AI 总结异步 ═══════
        drv, wait = _recycle(drv)        # 换全新会话，规避累积状态崩溃
        print('\n== 3) /meta 「AI 总结」异步 ==', flush=True)
        drv.get(base + '/meta')
        wait.until(lambda d: d.find_elements(By.ID, 'sumText'))
        drv.execute_script("setProject(arguments[0]);", MPROJ)
        # 用 JS 直接调 aiSummary：core14 的按钮挂在未激活的页签里，
        # 原生 click 会 ElementNotInteractable（元素不可见）。JS 调用
        # 同样走真实前端逻辑（轮询 + 结果渲染），只是绕开可见性限制。
        drv.execute_script(
            "const b=[...document.querySelectorAll('button')].find("
            "  x=>(x.getAttribute('onclick')||'').includes(\"aiSummary('core14'\"));"
            "if(!b) throw new Error('未找到 AI 总结按钮');"
            "aiSummary('core14', b);")
        # 抓「生成中」中间态
        saw_busy = False
        t0 = time.time()
        while time.time() - t0 < 60:
            if drv.find_elements(By.XPATH,
                                 "//button[contains(., 'AI 生成中')]"):
                saw_busy = True
            if not drv.find_elements(By.ID, 'sumModal'):
                break
            disp = drv.find_element(By.ID, 'sumModal').value_of_css_property(
                'display')
            if disp != 'none':
                break
            time.sleep(0.3)
        modal = drv.find_element(By.ID, 'sumModal')
        val = drv.find_element(By.ID, 'sumText').get_attribute('value')
        check('（3）AI 任务完成后弹窗展示摘要',
              modal.value_of_css_property('display') != 'none'
              and val.strip() == CANNED,
              f'display={modal.value_of_css_property("display")} value={val[:60]!r}')
        check('（3）过程中出现过「AI 生成中」忙碌态（不是同步卡死）',
              saw_busy, f'saw_busy={saw_busy}')

        # ═══════ 4) CDS/PEP 导出异步 ═══════
        drv, wait = _recycle(drv)        # 换全新会话，规避累积状态崩溃
        print('\n== 4) /cds-export 「导出」异步 ==', flush=True)
        drv.get(base + '/cds-export')
        wait.until(lambda d: d.find_elements(By.ID, 'cpColl'))
        drv.execute_script(
            "const s=document.getElementById('cpColl');"
            "s.value=arguments[0];"
            "if (typeof cdsPickOpen==='function') cdsPickOpen();", COLL)
        wait.until(lambda d: d.find_elements(
            By.CSS_SELECTOR, '#cpTable input[type=checkbox]'))
        boxes = drv.find_elements(By.CSS_SELECTOR,
                                  '#cpTable input[type=checkbox]')
        check('（4）挑选表里有可勾选的 CDS', bool(boxes), f'{len(boxes)} 个')
        if boxes:
            drv.execute_script("arguments[0].click();", boxes[0])
        exported = drv.find_elements(By.XPATH,
                                     "//button[contains(., '导出')]")
        check('（4）找到导出按钮', bool(exported),
              f'{len(exported)} 个候选')
        if exported:
            drv.execute_script("arguments[0].click();", exported[0])
            t0 = time.time()
            msg = ''
            while time.time() - t0 < 90:
                msg = drv.find_element(By.ID, 'cpMsg').text
                if '✔ 导出' in msg or '失败' in msg:
                    break
                time.sleep(0.5)
            check('（4）导出完成并显示条数/基因数',
                  '✔ 导出' in msg and '条 CDS' in msg, f'{msg[:160]!r}')
            check('（4）导出任务进了任务中心',
                  any(x['name'].startswith('CDS/PEP 导出')
                      for x in tm.list_all(logs=False)),
                  str([x['name'] for x in tm.list_all(logs=False)][:5]))

        print('\n== 控制台无 JS 报错 ==', flush=True)
        try:
            logs = drv.get_log('browser')
            js_err = [l for l in logs if l.get('level') == 'SEVERE'
                      and l.get('source') == 'javascript']
            check('无 SEVERE 级 JS 错误', not js_err,
                  str([l.get('message', '')[:120] for l in js_err])[:250])
        except Exception as e:
            skip('控制台日志', f'{type(e).__name__}')

    finally:
        (ss.build_sqn, ncbi_download.esearch, vpm._meta_deepseek_key,
         vpm._meta_ai_chat, vpm._meta_ai_client, vpm._meta_load_df) = _orig
        if drv is not None:
            try:
                drv.quit()
            except Exception:
                pass
        try:
            srv.shutdown()
        except Exception:
            pass
        for p in _CLEAN:
            shutil.rmtree(p, ignore_errors=True)

    print('\n' + ('ALL PASS' if not FAIL
                  else f'{len(FAIL)} FAIL: {FAIL}'), flush=True)
    return 0 if not FAIL else 1


def main():
    """跑一次；WebDriver 会话偶发原生崩溃时重试一次（见 _it_download_ui）。"""
    try:
        from selenium.common.exceptions import WebDriverException
    except ImportError:
        return _run_once()
    for attempt in (1, 2):
        try:
            return _run_once()
        except WebDriverException as e:
            if attempt == 2:
                print(f'  [FAIL] WebDriver 会话连续两次中断：'
                      f'{type(e).__name__}: {str(e)[:120]}',
                      flush=True)
                return 1
            print(f'  [RETRY] WebDriver 会话中断（{type(e).__name__}: '
                  f'{str(e)[:80]}），重试一次…',
                  flush=True)
            FAIL.clear()
    return 1


if __name__ == '__main__':
    sys.exit(main())
