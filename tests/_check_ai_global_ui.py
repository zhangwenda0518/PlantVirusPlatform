# -*- coding: utf-8 -*-
"""AI 对话坞全局化与易用性巡检（Playwright）。

运行: python tests/_check_ai_global_ui.py
自起 Flask（端口 18771），跑完即关。

验证点：
  1. 对话坞随 _nav 全局化：/results 等非工具页也有 FAB 与坞
     （原列举的 /pipeline 已于 2026-10-11 归档）
  2. 闸门待决事件自动弹出对话坞（唯一自动开坞的事件）
  3. 标题未读计数（(N) 前缀）随红点出现/清零
  4. 对话历史 localStorage 持久化：写入历史 → 刷新 → 恢复
  5. /api/ai/config 带 stage_options；设置页渲染两组阶段勾选框
  6. 全程无 JS 控制台错误
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _fixture_smoke_e2e import ensure as _ensure_smoke_e2e  # noqa: E402
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
    _ensure_smoke_e2e()
    from werkzeug.serving import make_server
    import app as app_module
    from Virus_Platform_Core.ai_agent import broadcast
    from Virus_Platform_Core.web.state import cfg
    srv = make_server('127.0.0.1', 18771, app_module.app, threaded=True)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = 'http://127.0.0.1:18771'
    time.sleep(1.0)

    orig_ai = dict(cfg.ai)
    cfg.ai = {'enabled': True, 'provider': 'custom', 'base_url': 'x',
              'api_key': 'k', 'model': 'm', 'temperature': 0.3,
              'max_tokens': 2048, 'explain': True, 'broadcast': True,
              'gate': True, 'gate_timeout_min': 15,
              'broadcast_stages': ['kvsuite'], 'gate_stages': ['kvsuite']}
    broadcast._ai_assess = (lambda s, st: '建议继续：数据正常。')

    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(
                headless=True,
                args=['--no-sandbox', '--disable-dev-shm-usage',
                      '--disable-gpu'])
            ctx = browser.new_context(viewport={'width': 1440, 'height': 900})
            page = ctx.new_page()
            errors = []
            page.on('console', lambda m: errors.append(m.text)
                    if m.type == 'error' else None)
            page.on('pageerror', lambda e: errors.append(str(e)))

            # -- 1) 全局坞：非总览页 --
            page.goto(base + '/results', wait_until='load')
            page.wait_for_timeout(1200)
            check('results 页有 FAB', page.locator('#aiFab').count() == 1)
            check('results 页有对话坞', page.locator('#aiDock').count() == 1)

            # -- 2/3) 闸门待决自动弹坞 + 标题计数 --
            # 先造一条普通事件（坞关着 → 未读+标题）
            broadcast.publish('ai_read', sample='smoke_e2e', stage='kvsuite',
                              name='②b', text='**速读**：正常。')
            page.wait_for_timeout(2000)
            t0 = page.title()
            check('标题带未读前缀', t0.startswith('('), t0)
            th = threading.Thread(
                target=lambda: broadcast.gate_point('smoke_e2e', 'kvsuite',
                                                    timeout_s=30),
                daemon=True)
            th.start()
            deadline = time.time() + 10
            while time.time() < deadline and not broadcast.pending_gate():
                time.sleep(0.1)
            page.wait_for_timeout(2500)
            check('闸门待决 → 对话坞自动弹出',
                  page.locator('#aiDock').is_visible())
            page.locator('.ai-gate-btns button').first.click()
            th.join(timeout=8)
            page.wait_for_timeout(1200)
            check('标题未读随展开清零', not page.title().startswith('('),
                  page.title())

            # -- 4) 历史持久化 --
            page.locator('#aiInput').fill('持久化测试消息')
            page.click('#aiSendBtn')
            page.wait_for_timeout(1500)          # 请求（无真实配置→降级错误也算内容）
            page.goto(base + '/results', wait_until='load')
            page.wait_for_timeout(800)
            page.click('#aiFab')
            page.wait_for_timeout(600)
            bubbles = page.locator('#aiMsgs .ai-user').count()
            check('刷新后对话历史恢复', bubbles >= 1, str(bubbles))
            page.locator('#aiDock .ai-ico[title*="清空"], #aiDock .ai-ico').first.click()
            page.evaluate("localStorage.removeItem('vp_ai_msgs')")

            # -- 4b) AI 执行单卡片：真实创建待决单 → 渲染 → 驳回闭环 --
            from Virus_Platform_Core.ai_agent import executors
            rec, err = executors.create_action(
                'retry_failed', 'smoke_e2e',
                reason='UI 巡检测试执行单', source='chat')
            check('执行单创建成功', rec is not None, str(err))
            deadline = time.time() + 10
            while time.time() < deadline and not page.locator(
                    f'.ai-action[data-aid="{rec["id"]}"]').count():
                page.wait_for_timeout(300)
            check('执行单卡片渲染（含批准/驳回按钮）',
                  page.locator(f'.ai-action[data-aid="{rec["id"]}"] '
                               '.ai-gate-btns button').count() == 2)
            page.locator(f'.ai-action[data-aid="{rec["id"]}"] '
                         '.ai-gate-btns button').last.click()
            deadline = time.time() + 8
            while time.time() < deadline:
                a = executors.get_action(rec['id'])
                if a and a['status'] == 'rejected':
                    break
                time.sleep(0.2)
            check('点驳回 → 执行单状态翻转 rejected',
                  (executors.get_action(rec['id']) or {}).get('status')
                  == 'rejected')
            page.wait_for_timeout(1200)
            check('执行单卡片原地更新为已驳回',
                  page.locator(f'.ai-action[data-aid="{rec["id"]}"] '
                               '.ai-gate-st.done').count() == 1)

            # -- 5) stage_options 与设置页勾选组 --
            r = page.evaluate(
                "fetch('/api/ai/config').then(r => r.json())")
            check('config 带 stage_options',
                  len(r.get('stage_options') or []) >= 15)
            page.goto(base + '/settings', wait_until='load')
            page.locator('#aiCard').scroll_into_view_if_needed()
            page.wait_for_timeout(600)
            nb = page.locator('#aiBcastStages input').count()
            ng = page.locator('#aiGateStages input').count()
            check(f'设置页阶段勾选组渲染（播报 {nb} / 闸门 {ng}）',
                  nb >= 15 and ng >= 15)
            page.screenshot(path='run/logs/ai_settings_stages.png')

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
    print('✔ 全局坞与易用性巡检全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
