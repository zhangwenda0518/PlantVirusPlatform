# -*- coding: utf-8 -*-
"""把已建成的 RefSeq 通用病毒分类库（kunpeng_db/ref）打成可选下载件。

2026-10-25 定案：数据包不再随附 kunpeng_db/ref（-680M，植物专用分发省
30%）；ref 改为「按需下载的可选件」——走插件分发链（索引 + sha256 +
三源下载），安装时经 install_to 落回 databases/kunpeng_db/ref，数据库页
扫描自动列出，工具②/④下拉即可选用。

为什么 install_to 指到 ref 本级、zip 内文件平铺：
  `_swap_into_place` 是**整目录替换**——install_to 若指到 kunpeng_db 父级，
  安装会把同级的 virus_db 一起换掉（实测代码确认）。指到 ref 本级则
  dest 不存在时干净新建、已装旧版时整目录替换，均不伤及 virus_db。
  `databases` 映射是 kind=database 的清单必填项（校验要求指向包内真实
  存在的文件），选 hash_1.k2d 仅为此；resolve_plugins 对 install_to 型
  插件的映射路径永远不在位 → 静默跳过，不会污染 cfg。

用法:
    python dev_tools/build_ref_db_pack.py            # 打包 + 更新索引
    python dev_tools/build_ref_db_pack.py --check    # 只校验索引是否最新
"""
import argparse
import hashlib
import json
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'databases', 'kunpeng_db', 'ref')
DIST = os.path.join(ROOT, 'plugins-dist')
INDEX = os.path.join(ROOT, 'plugins-repo', 'index.json')
NAME = 'ref-db'

_FILES = ('hash_1.k2d', 'taxo.k2d', 'opts.k2d', 'hash_config.k2d',
          'seqid2taxid.map', 'refvirus_acc2taxid.tsv')

MANIFEST = {
    'name': NAME,
    'version': '2026-10.1',
    'kind': 'database',
    'description': 'RefSeq 通用病毒分类库（kunpeng/k2 格式，~680M）——'
                   '可选下载件：安装后自动落位 databases/kunpeng_db/ref，'
                   '数据库页与工具②/④病毒库下拉即会列出，选中即用；'
                   '植物主库 virus_db 不受影响。',
    'install_to': 'databases/kunpeng_db/ref',
    # kind=database 的清单必填映射（校验要求指向包内真实存在的文件）；
    # install_to 型插件安装后包内不再保留文件，resolve_plugins 会跳过，
    # 不会改写 cfg——切换病毒库由用户在工具②/④下拉 / 构建页完成。
    'databases': {'ref': 'hash_1.k2d'},
}


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def build():
    missing = [f for f in _FILES if not os.path.isfile(os.path.join(SRC, f))]
    if missing:
        raise SystemExit(f'源库不完整，缺 {missing}——先建 ref 库再打包')
    os.makedirs(DIST, exist_ok=True)
    out_zip = os.path.join(DIST, f'{NAME}-plugin.zip')
    with zipfile.ZipFile(out_zip, 'w', zipfile.ZIP_STORED) as z:
        # 哈希表近乎不可压缩，STORED 省时省 CPU。
        # 布局：plugin.json + 6 个文件平铺在插件根——install_to 指向
        # databases/kunpeng_db/ref 时，sub 的顶层条目原样成为 dest 内容；
        # 若再包一层 ref/ 子目录会嵌成 kunpeng_db/ref/ref/（实测踩坑）。
        z.writestr(f'{NAME}/plugin.json',
                   json.dumps(MANIFEST, ensure_ascii=False, indent=1))
        for f in _FILES:
            z.write(os.path.join(SRC, f), f'{NAME}/{f}')
    size = os.path.getsize(out_zip)
    sha = sha256_file(out_zip)
    entry = {
        'name': NAME, 'version': MANIFEST['version'], 'kind': 'database',
        'description': MANIFEST['description'],
        'file': f'{NAME}-plugin.zip',
        # Gitee 附件上限 100MB，超限不挂；镜像自传 + GitHub Release
        'url': f'http://39.106.101.94/plugins/{NAME}-plugin.zip',
        'urls': [f'http://39.106.101.94/plugins/{NAME}-plugin.zip',
                 'https://github.com/zhangwenda0518/PlantVirusPlatform/'
                 'releases/download/plugin-packs/'
                 f'{NAME}-plugin.zip'],
        'size': size, 'sha256': sha,
        'install_to': MANIFEST['install_to'],
    }
    with open(INDEX, encoding='utf-8') as f:
        idx = json.load(f)
    idx['plugins'] = [p for p in idx['plugins'] if p['name'] != NAME]
    idx['plugins'].append(entry)
    with open(INDEX, 'w', encoding='utf-8') as f:
        json.dump(idx, f, ensure_ascii=False, indent=1)
        f.write('\n')
    print(f'✔ {out_zip}  {size / 1e6:.1f}MB  sha256={sha[:16]}…')
    print(f'✔ 索引已登记（共 {len(idx["plugins"])} 个插件包）')
    print('  待办：把 zip 上传镜像站 docs/plugins/ 与 GitHub Release '
          '(plugin-packs)')


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
