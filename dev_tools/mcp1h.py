# -*- coding: utf-8 -*-
"""手动验证用（HTTP 长驻版）：连 http://127.0.0.1:8915/mcp 调用一个工具。

用法：
    3rd\\python\\python.exe dev_tools\\mcp1h.py <vp_工具名> ['{"参数JSON"}']
前提：
    3rd\\python\\python.exe mcp_server.py --http 8915 --token nighttest
"""
import asyncio
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'dev_tools'))
os.chdir(ROOT)

from test_mcp_e2e import PY, SERVER  # noqa: E402
from mcp import ClientSession  # noqa: E402
from mcp.client.streamable_http import streamablehttp_client  # noqa: E402

URL = 'http://127.0.0.1:8915/mcp'
HDRS = {'Authorization': 'Bearer nighttest'}


async def main():
    tool = sys.argv[1] if len(sys.argv) > 1 else 'vp_platform_info'
    args = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    async with streamablehttp_client(URL, headers=HDRS) as (r, w, _sid):
        async with ClientSession(r, w) as s:
            await s.initialize()
            res = await asyncio.wait_for(s.call_tool(tool, args), 330)
            print('isError:', getattr(res, 'isError', False))
            for c in res.content:
                if getattr(c, 'type', '') == 'text':
                    print(c.text)


if __name__ == '__main__':
    asyncio.run(main())
