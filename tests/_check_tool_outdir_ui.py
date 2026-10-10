# -*- coding: utf-8 -*-
"""样本处理四卡「输出目录（可选）」的真实浏览器验收（Edge headless）。

锁定：t-convert / t-fastp / t-prepchain / hostremoval 页的输出目录输入框
渲染；convert 端到端（填自定义目录 → 运行 → 输出区显示外部 out_dir）。
用法: python tests/_check_tool_outdir_ui.py
"""
import os
import shutil
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

    OUT = os.path.join(ROOT, 'run', '_uitest_outdir_ui')
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)

    opts = EdgeOptions()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--disable-gpu')
    opts.add_argument('--window-size=1440,1000')
    drv = webdriver.Edge(options=opts)
    try:
        drv.set_page_load_timeout(30)
        wait = WebDriverWait(drv, 90)

        drv.get(base + '/tools?g=sample')
        wait.until(lambda d: d.find_elements(By.ID, 'cv_input'))
        for mod, fid in (('t-convert', 'cv_outdir'),
                         ('t-fastp', 'f_outdir'),
                         ('t-prepchain', 'pc_outdir')):
            drv.execute_script("showModule && showModule(arguments[0]);", mod)
            time.sleep(0.3)
            el = drv.find_element(By.ID, fid)
            check(f'{mod} 输出目录输入框渲染',
                  el.is_displayed()
                  and 'browseDir(' in (el.find_element(By.XPATH, '..')
                                       .find_element(By.TAG_NAME, 'button')
                                       .get_attribute('onclick')))

        drv.get(base + '/hostremoval')
        wait.until(lambda d: d.find_elements(By.ID, 'hr_outdir'))
        check('hostremoval 页输出目录输入框渲染',
              drv.find_element(By.ID, 'hr_outdir').is_displayed())

        # 端到端：convert + 自定义输出目录
        drv.get(base + '/tools?g=sample')
        wait.until(lambda d: d.find_elements(By.ID, 'cv_input'))
        drv.execute_script("showModule && showModule('t-convert')")
        time.sleep(0.3)
        r1 = os.path.join(ROOT, 'examples', 'example_R1.fastq.gz')
        drv.execute_script('arguments[0].value = arguments[1];',
                           drv.find_element(By.ID, 'cv_input'), r1)
        drv.execute_script('arguments[0].value = arguments[1];',
                           drv.find_element(By.ID, 'cv_outdir'), OUT)
        btn = drv.find_element(
            By.XPATH, "//button[contains(@onclick,'runTool')]")
        drv.execute_script("runTool('convert', arguments[0]);", btn)

        def done(d):
            t = d.find_element(By.ID, 'toolrun-convert').text
            return 'fa.gz' in t and ('已完成' in t or 'failed' in t.lower())
        wait.until(done)
        box = drv.find_element(By.ID, 'toolrun-convert').text
        check('输出区显示外部输出目录', '_uitest_outdir_ui' in box,
              box[:120].replace('\n', ' | '))
        check('任务已完成', '已完成' in box)
        # 产物确实落在外部目录
        subs = [d for d in os.listdir(OUT) if d.startswith('convert_')]
        check('产物落盘 <out_dir>/<运行名>/', bool(subs) and any(
            f.endswith('.fa.gz')
            for f in os.listdir(os.path.join(OUT, subs[0]))), subs)
    finally:
        drv.quit()
        srv.shutdown()
        shutil.rmtree(OUT, ignore_errors=True)

    print('\n' + ('ALL PASS' if not FAIL else f'FAILED: {FAIL}'))
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
