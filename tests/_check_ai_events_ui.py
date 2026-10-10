# -*- coding: utf-8 -*-
"""主动播报事件流的 UI 巡检（Playwright）：真实浏览器渲染 + 闸门按钮闭环。

运行: python tests/_check_ai_events_ui.py
自起 Flask（端口 18769），跑完即关。

流程：
  1. 页面加载（对话坞关闭）→ 依次发布 stage_start / ai_read / ai_gate(pending)
  2. FAB 未读红点计数；展开后红点清零、三类卡片渲染
  3. 点击闸门「继续」→ 待决闸门被 API 真实放行（线程返回 continue）
  4. 全程无 JS 控制台错误；截图供目检
"""
import os
import sys
import threading
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('VP_NO_RECOVER', '1')
# 本守卫会发布大量事件 / 建执行单；关掉审计落盘，避免往 run/ai_audit/
# 写测试垃圾（审计本身的正确性由 tests/_check_ai_audit.py 在临时目录专测）。
os.environ.setdefault('VP_AI_AUDIT', '0')

# 与其他浏览器探针（_check_kvchain_ui / _it_download_ui 等）约定一致：
# 未装 playwright（开发依赖缺失）时 SKIP 并返回 0，不让加强项把整批测试变红。
try:
    from playwright.sync_api import sync_playwright
except ImportError:                     # pragma: no cover - 取决于本机环境
    sync_playwright = None

FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


def main():
    if sync_playwright is None:
        print('[SKIP] 未安装 playwright，跳过真浏览器验证')
        print('       安装: python -m pip install -r requirements-dev.txt')
        return 0
    from werkzeug.serving import make_server
    import app as app_module
    from Virus_Platform_Core.ai_agent import broadcast
    from Virus_Platform_Core.web.state import cfg
    srv = make_server('127.0.0.1', 18769, app_module.app, threaded=True)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = 'http://127.0.0.1:18769'
    time.sleep(1.0)

    orig_ai = dict(cfg.ai)
    # enabled 必须为 True：闸门 arming 条件是 enabled && gate（评估已打桩，
    # 不依赖真实 AI 端点）
    cfg.ai = {'enabled': True, 'provider': 'custom', 'base_url': 'x',
              'api_key': 'k', 'model': 'm', 'temperature': 0.3,
              'max_tokens': 2048, 'explain': True, 'broadcast': True,
              'gate': True, 'gate_timeout_min': 15,
              'broadcast_stages': ['kvsuite'], 'gate_stages': ['kvsuite']}
    # 闸门评估打桩（不依赖 AI 配置）
    broadcast._ai_assess = (lambda s, st:
                            '建议继续：检出 3 种病毒，定量与组装质量正常。')

    out = {}
    def _gate_thread():
        out['r'] = broadcast.gate_point('smoke_e2e', 'kvsuite', timeout_s=60)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(
                headless=True,
                args=['--no-sandbox', '--disable-dev-shm-usage',
                      '--disable-gpu'])
            page = browser.new_context(
                viewport={'width': 1440, 'height': 900}).new_page()
            errors = []
            page.on('console', lambda m: errors.append(m.text)
                    if m.type == 'error' else None)
            page.on('pageerror', lambda e: errors.append(str(e)))
            # 常驻 EventSource 会让 networkidle 永不触发，用 load（与
            # _check_pages_console 的回退策略一致）
            page.goto(base + '/', wait_until='load')
            time.sleep(1.0)                       # 等 EventSource 建连

            # 发布三类事件（闸门事件由真实 gate_point 线程产生）
            broadcast.note_stage_start('smoke_e2e', 'kvsuite')
            broadcast.publish('ai_read', sample='smoke_e2e', stage='kvsuite',
                              name='②b 已知病毒识别与定量',
                              text='**速读**：检出 3 种病毒，3600 fragments。')
            time.sleep(2.0)                       # 等 SSE 推达
            check('坞关闭时 FAB 红点出现',
                  page.locator('#aiFabBadge').is_visible())
            check('标题带未读前缀', page.title().startswith('('), page.title())

            # 闸门待决（由真实 gate_point 线程产生）→ 应自动弹出对话坞
            # （2026-10-02 起新行为：闸门阻塞流水线，必须有人拍板）
            th = threading.Thread(target=_gate_thread, daemon=True)
            th.start()
            deadline = time.time() + 10
            while time.time() < deadline and not broadcast.pending_gate():
                time.sleep(0.1)
            deadline = time.time() + 10
            while time.time() < deadline and not page.locator('#aiDock').is_visible():
                page.wait_for_timeout(300)
            check('闸门待决 → 对话坞自动展开',
                  page.locator('#aiDock').is_visible())
            page.wait_for_timeout(600)
            check('展开后红点清零', not page.locator('#aiFabBadge').is_visible())
            check('讲解卡片渲染',
                  page.locator('.ai-ev-start').count() >= 1)
            check('AI 播报卡片渲染',
                  page.locator('.ai-read').count() >= 1)
            check('闸门卡片渲染且按钮可用',
                  page.locator('.ai-gate .ai-gate-btns button').count() == 2
                  and page.locator('.ai-gate-st.pending').count() >= 1)
            page.screenshot(path='run/logs/ai_events_feed.png')

            page.locator('.ai-gate-btns button').first.click()
            th.join(timeout=8)
            check('点继续 → 闸门线程放行', out.get('r') == 'continue',
                  str(out))
            page.wait_for_timeout(1500)
            check('闸门卡片原地翻转为已继续',
                  page.locator('.ai-gate-st.done').count() == 1
                  and page.locator('.ai-gate').count() == 1)
            page.screenshot(path='run/logs/ai_events_gate_done.png')

            real_errors = [e for e in errors
                           if 'favicon' not in e and '404' not in e]
            check(f'无 JS 控制台错误（{len(real_errors)} 条）',
                  not real_errors, '; '.join(real_errors[:3]))
            browser.close()
    finally:
        cfg.ai = orig_ai
        cfg.save()
        srv.shutdown()

    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 播报事件流 UI 巡检全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
