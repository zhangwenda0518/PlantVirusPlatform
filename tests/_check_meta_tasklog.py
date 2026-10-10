# -*- coding: utf-8 -*-
"""meta 页「④ 运行任务」卡片的**真实浏览器**验收（Selenium + Edge headless）。

背景：④ 的任务卡原本只有状态/最后一条消息，「已完成」的任务没有任何入口
回看它做过什么（历史日志在任务中心 /tasks 有，但④没有）。本次改动给每张
卡加了「日志」<details>（懒加载 /api/task/<id>/log）与「任务中心」入口。

锁定两点：
  1. 点开「日志」后能看到任务日志文本（而非永远「加载中…」）；
  2. 卡片列表底部有通往 /tasks（任务中心·历史归档）的链接。
  3. 展开状态跨 10s 周期重绘保持（_mtLogOpen）——展开后等到下一轮
     重绘再断言 details 仍为 open 且文本仍在。

用法: python tests/_check_meta_tasklog.py [--headed]
"""
import argparse
import os
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


def check(label, ok, extra=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {label}' + (f'  {extra}' if extra else ''),
          flush=True)
    if not ok:
        FAIL.append(label)


def free_port(preferred):
    for p in (preferred, 0):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', p))
                return s.getsockname()[1]
        except OSError:
            continue
    raise RuntimeError('无可用端口')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=8807)
    ap.add_argument('--headed', action='store_true')
    a = ap.parse_args()

    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.edge.options import Options as EdgeOptions
        from selenium.webdriver.support.ui import WebDriverWait
    except ImportError as e:
        print(f'[SKIP] selenium 未安装: {e}')
        return 0

    import app as appmod
    from Virus_Platform_Core.web.tasks import tm

    # ── 1) 种一个「公共检索」任务（名字前缀决定它出现在④） ──────
    MARK = 'UITEST-LOG-MARKER-6F3A'

    def job(emit, progress, cancel):
        emit('[00:00:01] 检索 Lycium barbarum（验收用）')
        emit(MARK)
        return {'n': 1}

    tid = tm.start('公共检索 UITEST日志验收', job, weight='light')
    deadline = time.time() + 15
    while time.time() < deadline:
        snap = tm.snapshot(tid, log_lines=5)
        if snap and snap.get('status') != 'running':
            break
        time.sleep(0.2)

    # ── 2) 起服务（不碰用户的 8989） ─────────────────────────────
    from werkzeug.serving import make_server
    port = free_port(a.port)
    srv = make_server('127.0.0.1', port, appmod.app, threaded=True)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f'http://127.0.0.1:{port}'
    print(f'临时实例: {base}\n', flush=True)

    opts = EdgeOptions()
    if not a.headed:
        opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--disable-gpu')
    opts.add_argument('--window-size=1440,1000')
    drv = webdriver.Edge(options=opts)
    try:
        drv.set_page_load_timeout(30)
        wait = WebDriverWait(drv, 20)

        drv.get(base + '/meta')
        wait.until(lambda d: d.find_elements(
            By.CSS_SELECTOR, '#taskBox .card'))
        cards = drv.find_elements(By.CSS_SELECTOR, '#taskBox .card')
        check('④ 渲染出任务卡（含验收任务）', len(cards) >= 1,
              f'{len(cards)} 张')
        card = next((c for c in cards
                     if 'UITEST' in c.find_element(By.TAG_NAME, 'b').text),
                    cards[0])

        details = card.find_element(By.CSS_SELECTOR, 'details.tk-log-wrap')
        check('任务卡带「日志」折叠块', details.is_displayed())
        drv.execute_script('arguments[0].scrollIntoView({block:"center"});',
                           details)
        time.sleep(0.3)
        details.find_element(By.TAG_NAME, 'summary').click()
        wait.until(lambda d: MARK in d.find_element(
            By.CSS_SELECTOR, '[data-tid] .tk-log').text)
        check('展开日志能看到任务日志文本',
              MARK in details.find_element(By.CSS_SELECTOR, '.tk-log').text)

        # 展开状态跨周期重绘保持（done 任务 10s 一轮；等到下一轮再查）
        time.sleep(11)
        cards2 = drv.find_elements(By.CSS_SELECTOR, '#taskBox .card')
        card2 = next((c for c in cards2
                      if 'UITEST' in c.find_element(By.TAG_NAME, 'b').text),
                     None)
        det2 = card2.find_element(By.CSS_SELECTOR, 'details.tk-log-wrap')
        check('重绘后折叠块仍展开且日志文本保留',
              det2.get_attribute('open') and MARK in det2.text)

        link = drv.find_element(By.CSS_SELECTOR, '#taskBox a[href="/tasks"]')
        check('卡片列表底部有「任务中心」入口', link.is_displayed(),
              link.text.strip())
    finally:
        drv.quit()
        srv.shutdown()

    print('\n' + ('ALL PASS' if not FAIL else f'FAILED: {FAIL}'))
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
