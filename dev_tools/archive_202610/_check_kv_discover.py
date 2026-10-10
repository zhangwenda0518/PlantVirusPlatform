# -*- coding: utf-8 -*-
"""验证鉴定库自动发现：递归、剪枝、作用域、嵌套库索引复用。"""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import app  # noqa: E402

c = app.test_client()
d = c.get('/api/kv_index_list').get_json()
print('== /api/kv_index_list ==')
for l in d['libs']:
    print(f"scope={l['scope']:<9} name={l['name']:<12} "
          f"salmon={l['salmon']} minibwa={l['minibwa']} "
          f"ref={bool(l.get('reference'))} path={l['path']}")

print('== /api/datasets 卡 ==')
d2 = c.get('/api/datasets').get_json()
for it in d2['virusref']:
    print(f"name={it['name']:<12} kunpeng={it['kunpeng_name'] or '-'} "
          f"ready={it['kunpeng_ready']} path={it['path']}")

from Virus_Platform_Core.kv_stage import index_dir_for  # noqa: E402
ref = os.path.join('databases', 'mylibs', 'nested', 'test_lib',
                   'reference.fasta')
print('嵌套库索引复用:', index_dir_for(ref, 'salmon'))
print('库内其它 fasta 不吃:', index_dir_for(
    os.path.join('databases', 'mylibs', 'nested', 'test_lib', 'other.fa'),
    'salmon'))
