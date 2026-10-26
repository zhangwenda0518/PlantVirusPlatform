# -*- coding: utf-8 -*-
"""MCP 控制台：页面可达 + 破坏面工具拦截（不拉起 MCP 服务器进程）。

背景（2026-10-26）：
  · MCP 控制台从首页正文迁到独立 /mcp 页 —— 它是开发者调试面板（手填 JSON、
    看原始返回），留在首页会截断"粘贴路径 → 选模块"的主线。
  · 同一轮补上破坏面拦截：此前 /api/mcp/call 是裸透传，普通用户在首页
    点开就能调 vp_delete_tool_run 删掉整个运行目录。内置 AI 对话早有
    _MCP_BLOCKED，网页端口径与之对齐。

本脚本只做**不触发子进程**的检查：
  1. /mcp 返回 200 且带控制台与外部接入区块
  2. vp_delete_tool_run 被 403 拒绝，且错误文案明确指出替代入口
  3. 拦截发生在拉起 mcp_server 之前（非法工具名同样 403，不 500）

用法：python tests/_check_mcp_console.py
"""
import os
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault('VP_NO_RECOVER', '1')

FAIL = []


def check(name, cond, detail=''):
    print(('  ok   ' if cond else '  FAIL ') + name
          + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


def main():
    import app as app_module
    c = app_module.app.test_client()

    # ---- 1. 独立页面可达 ----
    r = c.get('/mcp')
    html = r.data.decode('utf-8', 'replace')
    check('/mcp 返回 200', r.status_code == 200, f'HTTP {r.status_code}')
    check('/mcp 含控制台区块（mcpCfg / mcpTool）',
          'id="mcpCfg"' in html and 'id="mcpTool"' in html)
    check('/mcp 注入本机配置路径',
          'dataset.root' in html and 'mcp_server.py' in html)
    check('/mcp 不再出现在首页',
          'id="mcpCard"' not in c.get('/').data.decode('utf-8', 'replace'))

    # ---- 1b. 打开即自动连接 + 自动定位默认工具 ----
    # 2026-10-26 用户反馈："不能自动识别需要连接的模块吗，怎么还要单独选择"。
    # MCP 只有一个服务，连接与默认工具选择都应在页面加载时完成。
    check('/mcp 在加载时自动连接（autoConnect + DOMContentLoaded）',
          'autoConnect' in html and 'DOMContentLoaded' in html
          and 'MCPC.autoConnect(false)' in html)
    check('/mcp 连接后自动定位默认工具',
          "DEFAULT_TOOL" in html and 'vp_platform_info' in html
          and 'this.pick()' in html)
    # 自动连接失败不能弹 alert（用户没主动点），必须行内提示
    check('/mcp 自动连接失败走行内提示（不 alert）',
          'autoFail' in html)

    # ---- 2. 破坏面工具被拒 ----
    r = c.post('/api/mcp/call',
               json={'name': 'vp_delete_tool_run',
                     'arguments': {'run': 'no_such_run'}})
    body = r.get_json() or {}
    msg = str(body.get('error') or '')
    check('vp_delete_tool_run 被 403 拒绝', r.status_code == 403,
          f'HTTP {r.status_code} {body}')
    check('拒绝文案指向替代入口',
          '任务中心' in msg or '运行历史' in msg, msg)

    # ---- 3. 拦截在拉起服务器之前（不需要真调工具）----
    # 若拦截缺失，这里会尝试连 MCP 服务器并因缺参/超时而返回 4xx/5xx，
    # 文案也不会是"不开放"。
    check('拦截是前置的（文案为不开放，而非工具执行错误）',
          '不开放' in msg, msg)

    print()
    if FAIL:
        print(f'✘ {len(FAIL)} 项未通过:')
        for f in FAIL:
            print('   - ' + f)
        return 1
    print('✔ MCP 控制台巡检通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
