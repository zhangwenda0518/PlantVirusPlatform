# -*- coding: utf-8 -*-
"""格式转换「双端输入」的真实浏览器验收（Edge headless）。

2026-09-30 需求：格式转换支持 R1/R2 双端，双端分别转换、分别输出。

锁定：cv_input2 字段渲染；双端提交后输出区列出 R1/R2 两个产物、
日志含两条 seqkit 命令。
用法: python tests/_check_convert_paired.py
"""
import os
import sys
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


def main():
    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.edge.options import Options as EdgeOptions
        from selenium.webdriver.support.ui import WebDriverWait
    except ImportError as e:
        print(f'[SKIP] selenium 未安装: {e}')
        return 0

    import app as appmod

    from werkzeug.serving import make_server
    srv = make_server('127.0.0.1', 0, appmod.app, threaded=True)
    port = srv.server_address[1]
    import threading
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f'http://127.0.0.1:{port}'
    print(f'临时实例: {base}\n', flush=True)

    opts = EdgeOptions()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--disable-gpu')
    opts.add_argument('--window-size=1440,1000')
    drv = webdriver.Edge(options=opts)
    try:
        drv.set_page_load_timeout(30)
        wait = WebDriverWait(drv, 60)
        drv.get(base + '/tools?g=sample')
        wait.until(lambda d: d.find_elements(By.ID, 'cv_input'))
        drv.execute_script("showModule && showModule('t-convert')")
        time.sleep(0.5)

        inp2 = drv.find_element(By.ID, 'cv_input2')
        check('R2 输入框渲染', inp2.is_displayed())

        r1 = os.path.join(ROOT, 'examples', 'example_R1.fastq.gz')
        r2 = os.path.join(ROOT, 'examples', 'example_R2.fastq.gz')
        drv.execute_script('arguments[0].value = arguments[1];',
                           drv.find_element(By.ID, 'cv_input'), r1)
        drv.execute_script('arguments[0].value = arguments[1];', inp2, r2)
        btn = drv.find_element(
            By.XPATH, "//button[contains(@onclick,'runTool')]")
        drv.execute_script("runTool('convert', arguments[0]);", btn)

        def done(d):
            txt = d.find_element(By.ID, 'toolrun-convert').text
            return ('已完成' in txt or 'failed' in txt.lower()) and 'fa.gz' in txt
        wait.until(done)
        box = drv.find_element(By.ID, 'toolrun-convert').text
        check('输出区列出 R1 产物', 'example_R1.fa.gz' in box)
        check('输出区列出 R2 产物', 'example_R2.fa.gz' in box)
        check('运行日志含两条转换命令（R1 与 R2）',
              box.count('fq2fa') >= 2 and 'example_R1' in box
              and 'example_R2' in box,
              f"fq2fa 出现 {box.count('fq2fa')} 次")
        check('任务状态为已完成', '已完成' in box)
    finally:
        drv.quit()
        srv.shutdown()

    print('\n' + ('ALL PASS' if not FAIL else f'FAILED: {FAIL}'))
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
