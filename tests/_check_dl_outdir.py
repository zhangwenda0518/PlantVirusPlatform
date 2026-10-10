# -*- coding: utf-8 -*-
"""下载页「输出位置」的真实浏览器验收（Selenium + Edge headless）。

锁定：
  1. 创建表单渲染出「输出位置」输入框与浏览按钮；
  2. 带 out_dir 创建批次后，批次卡显示输出位置路径 + 打开按钮；
  3. 页面无 SEVERE 级 JS 错误。
用法: python tests/_check_dl_outdir.py
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
    from Virus_Platform_Core import public_data
    from Virus_Platform_Core.config import DIRS

    DL = os.path.join(ROOT, 'run', '_uitest_outdir_dl')
    CUSTOM = os.path.join(ROOT, 'run', '_uitest_outdir_custom', 'GSA数据')
    shutil.rmtree(DL, ignore_errors=True)
    shutil.rmtree(os.path.join(ROOT, 'run', '_uitest_outdir_custom'),
                  ignore_errors=True)
    os.makedirs(DL, exist_ok=True)
    DIRS['downloads'] = DL
    public_data._manager = None

    from werkzeug.serving import make_server
    srv = make_server('127.0.0.1', 0, appmod.app, threaded=True)
    port = srv.server_address[1]
    threading_srv = __import__('threading')
    threading_srv.Thread(target=srv.serve_forever, daemon=True).start()
    base = f'http://127.0.0.1:{port}'
    print(f'临时实例: {base}\n', flush=True)

    opts = EdgeOptions()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--disable-gpu')
    opts.add_argument('--window-size=1440,1000')
    drv = webdriver.Edge(options=opts)
    logs = []
    try:
        drv.set_page_load_timeout(30)
        wait = WebDriverWait(drv, 20)
        drv.get(base + '/download')
        wait.until(lambda d: d.find_elements(By.ID, 'dlOutDir'))
        check('表单渲染「输出位置」输入框', True)
        btn = drv.find_element(
            By.XPATH, "//button[contains(@onclick,\"browseDir('dlOutDir')\")]")
        check('浏览按钮存在且绑定 browseDir', btn.is_displayed())

        # 带 out_dir 走真实 API 创建（items 全注释 → 秒完成，不真下载）
        import json as _json
        import urllib.request as _ur
        req = _ur.Request(
            base + '/api/dl/create',
            data=_json.dumps({'name': 'uitest输出位置', 'items': ['# 注释'],
                              'convert_sra': False, 'out_dir': CUSTOM}
                             ).encode('utf-8'),
            headers={'Content-Type': 'application/json'}, method='POST')
        _ur.urlopen(req, timeout=15).read()
        wait.until(lambda d: d.find_elements(
            By.CSS_SELECTOR, '#batchList .card'))
        time.sleep(1.0)
        card = drv.find_element(By.CSS_SELECTOR, '#batchList .card')
        txt = card.text
        check('批次卡显示输出位置路径',
              '_uitest_outdir_custom' in txt.replace('\n', ' '),
              txt[:80].replace('\n', ' | '))
        ok_btn = card.find_elements(By.CSS_SELECTOR, 'button[data-d]')
        check('批次卡有「打开目录」按钮', bool(ok_btn))
        if ok_btn:
            check('打开按钮携带正确目录',
                  os.path.normpath(ok_btn[0].get_attribute('data-d'))
                  == os.path.normpath(CUSTOM))

        for entry in drv.get_log('browser'):
            if entry['level'] == 'SEVERE':
                logs.append(entry['message'][:120])
        check('无 SEVERE 级 JS 错误', not logs, str(logs[:2]))
    finally:
        drv.quit()
        srv.shutdown()
        shutil.rmtree(DL, ignore_errors=True)
        shutil.rmtree(os.path.join(ROOT, 'run', '_uitest_outdir_custom'),
                      ignore_errors=True)

    print('\n' + ('ALL PASS' if not FAIL else f'FAILED: {FAIL}'))
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
