# -*- coding: utf-8 -*-
"""样本处理四工具「输出目录（可选）」的 API 层验证。

convert / fastp / prepchain 真实跑通（宿主库缺失时 prepchain 自动跳过该
阶段），hostremoval 因宿主库未建会失败——但本需求改在**创建层**（run_dir
外移 + 指针登记），失败任务同样可验证外移正确。

断言：产物落 <out_dir>/<run_name>/；tool_runs/<run>/out_link.json 登记；
/api/tool/runs 列表含该 run 且 out_dir 指向外部；产物经
/tool_runs/<run>/<file> 可下载；删除运行后外部目录与指针一并清除。
"""
import json
import os
import shutil
import sys
import time

sys.path.insert(0, os.getcwd())
os.environ['VP_NO_RECOVER'] = '1'

import app as appmod
from Virus_Platform_Core.web.state import tool_runs_root

c = appmod.app.test_client()
R1 = os.path.join(os.getcwd(), 'examples', 'example_R1.fastq.gz')
R2 = os.path.join(os.getcwd(), 'examples', 'example_R2.fastq.gz')
OUT = os.path.join(os.getcwd(), 'run', '_uitest_toolout')
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT, exist_ok=True)

fails = []


def check(label, ok, extra=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {label}' + (f'  {extra}' if extra else ''))
    if not ok:
        fails.append(label)


def run_tool(tool, params, out_dir=OUT, timeout=300):
    r = c.post('/api/tool/run', json={'tool': tool, 'params': params,
                                      'out_dir': out_dir})
    if r.status_code != 200:
        return None, None, r
    d = r.get_json()
    tid, run = d['task'], d['run']
    deadline = time.time() + timeout
    while time.time() < deadline:
        s = c.get(f'/api/task/{tid}').get_json()
        if s['status'] != 'running':
            return run, s, r
        time.sleep(1)
    return run, {'status': 'timeout'}, r


# ── 1) convert：产物外移 + 指针 + 列表可见 + 下载可达 ──────────
run, s, r = run_tool('convert', {'input': R1, 'target': 'fasta'})
check('convert 外部目录返回 out_dir', r.get_json().get('out_dir', '').startswith(OUT),
      r.get_json().get('out_dir', ''))
proddir = os.path.join(OUT, run)
check('convert 产物目录 = <out_dir>/<run_name> 且有 .fa.gz',
      os.path.isdir(proddir) and any(f.endswith('.fa.gz') for f in os.listdir(proddir)),
      os.listdir(proddir) if os.path.isdir(proddir) else 'missing')
link = os.path.join(tool_runs_root(), run, 'out_link.json')
check('tool_runs/<run>/out_link.json 登记', os.path.isfile(link))
lst = c.get('/api/tool/runs').get_json()
mine = next((x for x in lst if x['name'] == run), None)
check('历史列表含该 run 且 out_dir 指向外部目录',
      bool(mine) and mine['out_dir'].replace('/', os.sep) == proddir,
      mine and mine['out_dir'])
fa = next((f for f in os.listdir(proddir) if f.endswith('.fa.gz')), None)
dl = c.get(f'/tool_runs/{run}/{fa}')
check('产物经 /tool_runs/<run>/<file> 可下载（指针寻址）',
      dl.status_code == 200 and len(dl.data) > 0)
dl.close()   # Windows：test client 的 response 持有文件流，不关会挡住后面的删目录
lst2 = c.get('/api/tool/runs').get_json()
mine2 = next((x for x in lst if x['name'] == run), None)
check('历史列表文件清单含产物（touch 刷新签名生效）',
      bool(mine2) and any(fa in f['path'] for f in mine2['files']),
      mine2 and [f['path'] for f in mine2['files']][:4])

# ── 2) fastp：真实跑通 + 产物外移（fastp 产物在 run 目录的 00_prep/ 子层）──
run2, s2, r2 = run_tool('fastp', {'r1': R1, 'r2': R2})
check('fastp 任务完成', s2['status'] == 'done', (s2.get('error') or '')[:60])
d2 = os.path.join(OUT, run2)
walk2 = [os.path.join(dp, f) for dp, _dn, fns in os.walk(d2) for f in fns] \
    if os.path.isdir(d2) else []
check('fastp 产物在外部目录（递归含 00_prep/）',
      any(f.endswith(('.html', '.json', '.fastq.gz')) for f in walk2),
      [os.path.relpath(f, d2) for f in walk2][:5])

# ── 3) prepchain：链式运行外移（宿主库缺失自动跳过该阶段） ─────
run3, s3, r3 = run_tool('prepchain', {'input': R1, 'input2': R2}, timeout=420)
check('prepchain 有外移目录', os.path.isdir(os.path.join(OUT, run3)),
      s3['status'] if s3 else str(r3.status_code))
check('prepchain 状态非 queued（已执行或完成）',
      s3 and s3['status'] in ('done', 'failed'),
      s3 and (s3.get('error') or s3['status'])[:60])

# ── 4) hostremoval：宿主库缺失会失败，但外移与登记照常 ─────────
run4, s4, r4 = run_tool('hostremoval', {'r1': R1, 'r2': R2, 'confidence': 0.2,
                                        'db': ''}, timeout=120)
if r4.status_code == 200:
    check('hostremoval 外移目录已建', os.path.isdir(os.path.join(OUT, run4)))
    check('hostremoval 指针已登记',
          os.path.isfile(os.path.join(tool_runs_root(), run4, 'out_link.json')))
else:
    check('hostremoval 创建被拒（400，宿主库缺失）', r4.status_code == 400,
          str(r4.status_code))

# ── 5) 删除 convert 运行：外部目录与指针一并清 ─────────────────
d = c.post(f'/api/tool/runs/{run}/delete')
check('删除运行 200', d.status_code == 200)
check('外部产物目录已清', not os.path.isdir(proddir))
check('指针目录已清', not os.path.isdir(os.path.join(tool_runs_root(), run)))

shutil.rmtree(OUT, ignore_errors=True)
print('\n' + ('ALL PASS' if not fails else f'FAILED: {fails}'))
sys.exit(0 if not fails else 1)
