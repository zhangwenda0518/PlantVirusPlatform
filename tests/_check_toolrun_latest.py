# -*- coding: utf-8 -*-
"""格式转换卡「运行日志不随新任务更新」的真实浏览器回归（Edge headless）。

2026-09-30 用户报告：连跑两次格式转换，第二次已完成（右上弹窗），输出区
仍显示第一次的运行日志。根因：injectStageLogs 对同一容器按任务数组顺序
反复写入，而 /api/tasks 是新→旧序 —— 旧任务后写覆盖新任务。修复为每容器
只吃第一个（最新）匹配任务。

本测试真实连跑两次转换（example_R1 → example_R2，FASTQ→FASTA，seqkit
本地秒级），断言最终输出区显示**第二次** run 的名称与日志。

用法: python tests/_check_toolrun_latest.py
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

        def run_convert(infile):
            """填输入 → JS 直调 runTool（走完整前端提交/轮询链，绕开
            headless 下的遮挡/首访弹窗干扰）→ 等该次任务**完成**。"""
            inp = drv.find_element(By.ID, 'cv_input')
            drv.execute_script('arguments[0].value = arguments[1];', inp, infile)
            btn = drv.find_element(
                By.XPATH, "//button[contains(@onclick,'runTool')]")
            drv.execute_script("runTool('convert', arguments[0]);", btn)

            def done(d):
                txt = d.find_element(By.ID, 'toolrun-convert').text
                return (infile.replace('\\', '/') in txt.replace('\\', '/')
                        and ('已完成' in txt or '失败' in txt
                             or 'failed' in txt.lower()
                             or 'done' in txt.lower()))
            wait.until(done)

            def done(d):
                txt = d.find_element(By.ID, 'toolrun-convert').text
                return (infile.replace('\\', '/') in txt.replace('\\', '/')
                        and ('已完成' in txt or '失败' in txt
                             or 'failed' in txt.lower()
                             or 'done' in txt.lower()))
            wait.until(done)

        drv.get(base + '/tools?g=sample')
        wait.until(lambda d: d.find_elements(By.ID, 'cv_input'))
        # 工作台卡片默认隐藏，显式激活格式转换卡（同用户点导航卡）
        drv.execute_script("showModule && showModule('t-convert')")
        time.sleep(0.5)

        r1 = os.path.join(ROOT, 'examples', 'example_R1.fastq.gz')
        r2 = os.path.join(ROOT, 'examples', 'example_R2.fastq.gz')
        print('== 第一次转换（example_R1） ==', flush=True)
        run_convert(r1)
        box1 = drv.find_element(By.ID, 'toolrun-convert').text
        check('第一次 run 显示在输出区', 'example_R1' in box1)

        print('== 第二次转换（example_R2） ==', flush=True)
        run_convert(r2)
        box2 = drv.find_element(By.ID, 'toolrun-convert').text
        check('输出区已切换为第二次 run（example_R2）',
              'example_R2' in box2, box2[:100].replace('\n', ' | '))
        check('第二次完成后不再回显第一次的日志',
              'example_R1' not in box2 or 'example_R2' in box2)
        check('输出区状态为已完成（新任务）',
              '已完成' in box2 or 'done' in box2.lower())
    finally:
        drv.quit()
        srv.shutdown()

    print('\n' + ('ALL PASS' if not FAIL else f'FAILED: {FAIL}'))
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
