# -*- coding: utf-8 -*-
"""手动验证用：单次调用一个 MCP 工具并打印原始结果。

用法：
    3rd\\python\\python.exe dev_tools\\mcp1.py <vp_工具名> ['{"参数JSON"}']

例：
    3rd\\python\\python.exe dev_tools\\mcp1.py vp_platform_info
    3rd\\python\\python.exe dev_tools\\mcp1.py vp_get_task "{\\"task_id\\": \\"abc\\"}"
"""
import asyncio
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'dev_tools'))
os.chdir(ROOT)

from test_mcp_e2e import PY, SERVER  # noqa: E402
from mcp import ClientSession, StdioServerParameters  # noqa: E402
from mcp.client.stdio import stdio_client  # noqa: E402


async def main():
    tool = sys.argv[1] if len(sys.argv) > 1 else 'vp_platform_info'
    args = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    params = StdioServerParameters(command=PY, args=[SERVER])
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            res = await asyncio.wait_for(s.call_tool(tool, args), 330)
            print('isError:', getattr(res, 'isError', False))
            for c in res.content:
                if getattr(c, 'type', '') == 'text':
                    print(c.text)


if __name__ == '__main__':
    asyncio.run(main())
