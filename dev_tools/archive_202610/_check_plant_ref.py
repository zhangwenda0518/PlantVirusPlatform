# -*- coding: utf-8 -*-
"""验证 plant_ref 库接入：发现、解析、索引复用、数据集卡、建库落点。"""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import app  # noqa: E402

c = app.test_client()

print('== /api/kv_index_list ==')
d = c.get('/api/kv_index_list').get_json()
for l in d['libs']:
    print(f"scope={l['scope']:<9} name={l['name']:<12} "
          f"salmon={l['salmon']} minibwa={l['minibwa']} ref={bool(l.get('reference'))}")

print('== /api/datasets 卡 ==')
d2 = c.get('/api/datasets').get_json()
for it in d2['virusref']:
    print(f"name={it['name']:<12} ref_mb={it['ref_mb']:<6} info={bool(it['ref_info'])} "
          f"kunpeng={it['kunpeng_name'] or '-'} ready={it['kunpeng_ready']}")

print('== resolve_kv_lib(plant_ref) ==')
from Virus_Platform_Core.kv_stage import (index_dir_for,  # noqa: E402
                                          resolve_kv_lib)
lib = resolve_kv_lib('plant_ref')
for k, v in lib.items():
    print(f'  {k}: {v}')
print('== 索引复用落点（salmon/minibwa，首次运行自建于库目录）==')
print('  salmon :', index_dir_for(lib['reference'], 'salmon'))
print('  minibwa:', index_dir_for(lib['reference'], 'minibwa'))

print('== kunpeng 一键建库落点预检 ==')
from Virus_Platform_Core.config import db_path  # noqa: E402
from Virus_Platform_Core.kunpeng import db_ready  # noqa: E402
kdir = os.path.join(os.path.dirname(db_path('virus', 'plant')), 'plant_ref')
print('  目标:', kdir, '| 已建:', db_ready(kdir), '（False=未建，构建页可一键）')

print('== TSWV 抽查 ==')
with io.open(lib['ref_info'], encoding='utf-8') as f:
    rows = [ln for ln in f if 'tomatomaculae' in ln.lower()]
print(f'  tomatomaculae 行: {len(rows)}')
for ln in rows:
    cols = ln.split('\t')
    print('   ', cols[0], cols[1], cols[4], cols[5])
