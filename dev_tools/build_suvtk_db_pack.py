# -*- coding: utf-8 -*-
"""把 NCBI 提交用 BFVD 注释库（misc_db/suvtk，212M）打成可选下载件。

2026-10-25 定案：misc_db/suvtk 不再随数据库包分发——只有「NCBI 提交准备」
的 suvtk 注释阶段用它，而提交场景本就必须联网；做成 suvtk-db 可选件后，
插件管理/模块自愈一键安装，install_to 还原 databases/misc_db/suvtk 原位。

布局与坑与 ref-db 相同：文件平铺插件根（install_to 平铺语义，再包一层
子目录会嵌套）；install_to 必须指到 suvtk 本级（_swap_into_place 是整目录
替换，指到 misc_db 父级会把同级的 viroids/prob 一起换掉）。

用法:
    python dev_tools/build_suvtk_db_pack.py            # 打包 + 更新索引
    python dev_tools/build_suvtk_db_pack.py --check    # 校验索引是否最新
"""
import argparse
import hashlib
import json
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'databases', 'misc_db', 'suvtk')
DIST = os.path.join(ROOT, 'plugins-dist')
INDEX = os.path.join(ROOT, 'plugins-repo', 'index.json')
NAME = 'suvtk-db'


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def build():
    if not os.path.isdir(SRC):
        raise SystemExit(f'源目录不存在: {SRC}')
    files = sorted(f for f in os.listdir(SRC) if os.path.isfile(
        os.path.join(SRC, f)))
    if not any(f.startswith('bfvd') for f in files):
        raise SystemExit('源目录内无 bfvd* 文件，先确认库完整')
    manifest = {
        'name': NAME,
        'version': '2026-10.1',
        'kind': 'database',
        'description': 'NCBI 提交准备 BFVD 功能注释库（suvtk，~212M）——'
                       '可选下载件：安装后自动落位 databases/misc_db/suvtk，'
                       '提交准备的 suvtk 注释阶段即可用；不影响其他模块。',
        'install_to': 'databases/misc_db/suvtk',
        'databases': {'suvtk': files[0]},
    }
    os.makedirs(DIST, exist_ok=True)
    out_zip = os.path.join(DIST, f'{NAME}-plugin.zip')
    with zipfile.ZipFile(out_zip, 'w', zipfile.ZIP_STORED) as z:
        # 表格/mmseqs 库压缩率低，STORED 省时；install_to 平铺布局
        z.writestr(f'{NAME}/plugin.json',
                   json.dumps(manifest, ensure_ascii=False, indent=1))
        for f in files:
            z.write(os.path.join(SRC, f), f'{NAME}/{f}')
    size = os.path.getsize(out_zip)
    sha = sha256_file(out_zip)
    entry = {
        'name': NAME, 'version': manifest['version'], 'kind': 'database',
        'description': manifest['description'],
        'file': f'{NAME}-plugin.zip',
        'url': f'http://39.106.101.94/plugins/{NAME}-plugin.zip',
        'urls': [f'http://39.106.101.94/plugins/{NAME}-plugin.zip',
                 'https://github.com/zhangwenda0518/plant-virus-analysis-'
                 'platform/releases/download/plugin-packs/'
                 f'{NAME}-plugin.zip'],
        'size': size, 'sha256': sha,
        'install_to': manifest['install_to'],
    }
    with open(INDEX, encoding='utf-8') as f:
        idx = json.load(f)
    idx['plugins'] = [p for p in idx['plugins'] if p['name'] != NAME]
    idx['plugins'].append(entry)
    with open(INDEX, 'w', encoding='utf-8') as f:
        json.dump(idx, f, ensure_ascii=False, indent=1)
        f.write('\n')
    print(f'✔ {out_zip}  {size / 1e6:.1f}MB  sha256={sha[:16]}…'
          f'  ({len(files)} 个文件)')
    print(f'✔ 索引已登记（共 {len(idx["plugins"])} 个插件包）')
    print('  待办：zip 上传 GitHub Release (plugin-packs)；镜像站可选')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        with open(INDEX, encoding='utf-8') as f:
            idx = json.load(f)
        ent = next((p for p in idx['plugins'] if p['name'] == NAME), None)
        z = os.path.join(DIST, f'{NAME}-plugin.zip')
        ok = ent and os.path.isfile(z) and \
            os.path.getsize(z) == ent.get('size') and \
            sha256_file(z) == ent.get('sha256')
        print('索引一致' if ok else '索引过期/缺包，重跑本脚本')
        return 0 if ok else 1
    build()
    return 0


if __name__ == '__main__':
    sys.exit(main())
