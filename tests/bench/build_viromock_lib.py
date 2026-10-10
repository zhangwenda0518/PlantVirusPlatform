# -*- coding: utf-8 -*-
"""构建 VIROMOCK 专用鉴定库：掺入株 + 11-18 全部隔离株 + 检出挑战目标参考。

产出 databases/virusref_db/viromock_kv/（reference.fasta + ref_info + 双引擎索引
由 known_virus_suite index 子命令另行构建）。
所有落盘路径均为本文件内的完整字面量常量；NCBI 请求经 _safe_fetch_ncbi 锁定
官方域名并拒绝私网解析。
"""
import ipaddress
import json
import os
import socket
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
# 平台根：默认自动推导（tests/bench/ 上溯 3 层 = 平台根），
# 可用环境变量 VP_ROOT 覆盖（换机/换目录不用改代码）。
PLATFORM = (os.environ.get('VP_ROOT')
            or os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))))
# 基准数据集（作者本机路径）：可用环境变量 VP_BENCH_DATA 覆盖
BENCH_DATA = os.environ.get('VP_BENCH_DATA', r'E:\谷歌下载\测试数据验证')
REPO_FASTA = os.path.join(BENCH_DATA, 'VIROMOCKchallenge-master',
                          'Datasets', 'fasta')
VKV = os.path.join(PLATFORM, 'databases', 'virusref_db', 'viromock_kv')
FA_P = os.path.join(VKV, 'reference.fasta')
INFO_P = os.path.join(VKV, 'reference.ref_info.tsv')
MANIFEST_P = os.path.join(VKV, 'manifest.json')
LOCAL_FASTA_PPV = os.path.join(REPO_FASTA,
                               'Plum_pox_virus_artificial_strain.fasta')
LOCAL_FASTA_PVY = os.path.join(REPO_FASTA,
                               'Potato_virus_Y_artificial_strain.fasta')

# (accession, taxid, species)
REFS = [
    # ── Dataset 11 PepMV ──
    ('DQ000985', 112229, 'Pepino mosaic virus'),
    ('AJ606359', 112229, 'Pepino mosaic virus'),
    ('MF422616', 112229, 'Pepino mosaic virus'),
    ('JQ314460', 112229, 'Pepino mosaic virus'),
    ('MK133092', 112229, 'Pepino mosaic virus'),
    ('HG313807', 112229, 'Pepino mosaic virus'),
    # ── Dataset 12 Cassava mosaic virus ──
    ('HE979770', 62079, 'African cassava mosaic virus'),
    ('HE979758', 10817, 'African cassava mosaic virus'),
    ('AJ314739', 223281, 'African cassava mosaic virus'),
    ('KR611579', 161378, 'African cassava mosaic virus'),
    # ── Dataset 13 BSV ──
    ('DQ451009', 69577, 'Banana streak virus'),
    ('KT895259', 69577, 'Banana streak virus'),
    ('DQ092436', 334778, 'Banana streak virus'),
    ('AY750155', 1411991, 'Banana streak virus'),
    ('KT895258', 69577, 'Banana streak virus'),
    ('AY493509', 328670, 'Banana streak virus'),
    # ── Dataset 14 PVY ──
    ('AB711147', 12216, 'Potato virus Y'),
    ('KC634004', 12216, 'Potato virus Y'),
    ('MF176828', 12216, 'Potato virus Y'),
    ('JQ969039', 12216, 'Potato virus Y'),
    ('FJ214726', 12216, 'Potato virus Y'),
    # ── Dataset 15 EMDV ──
    ('LN680656', 488317, 'Eggplant mottled dwarf virus'),
    ('FR751552', 488317, 'Eggplant mottled dwarf virus'),
    ('KJ082087', 488317, 'Eggplant mottled dwarf virus'),
    # ── Dataset 16 BPEV ──
    ('KX977568', 354328, 'Bell pepper endornavirus'),
    ('KR080326', 1711684, 'Bell pepper endornavirus'),
    ('JN019858', 354328, 'Bell pepper endornavirus'),
    ('JQ951943', 354328, 'Bell pepper endornavirus'),
    # ── Dataset 17 LChV1 ──
    ('MH300061', 217686, 'Little cherry virus 1'),
    ('KX192366', 217686, 'Little cherry virus 1'),
    ('EU715989', 217686, 'Little cherry virus 1'),
    ('LN794218', 217686, 'Little cherry virus 1'),
    ('MG934545', 217686, 'Little cherry virus 1'),
    # ── Dataset 18 BYDV ──
    ('EF521843', 2169986, 'Barley yellow dwarf virus'),
    ('KF523382', 224578, 'Barley yellow dwarf virus'),
    ('KY593456', 2169985, 'Barley yellow dwarf virus'),
    ('D11028', 2169984, 'Barley yellow dwarf virus'),
    ('KC559092', 2169988, 'Barley yellow dwarf virus'),
    ('EU332308', 2169986, 'Barley yellow dwarf virus'),
    # ── 半人工数据集掺入株 ──
    ('JQ911663', 10251, 'Citrus tristeza virus'),
    ('KU883267', 10251, 'Citrus tristeza virus'),
    ('MH323442', 10251, 'Citrus tristeza virus'),
    ('DQ377131', 19716, 'Grapevine yellow speckle viroid 2'),
    ('AY884983', 12216, 'Potato virus Y'),
    ('EF026076', 12216, 'Potato virus Y'),
    # ── 检出挑战目标（默认库缺失，2026-09-15 审计补充）──
    ('NC_002050', 10161, 'Tomato spotted wilt virus'),
    ('NC_002051', 10161, 'Tomato spotted wilt virus'),
    ('NC_002052', 10161, 'Tomato spotted wilt virus'),
    ('NC_040543', 1927579, 'Chenopodium quinoa mitovirus 1'),
    ('NC_005286', 43793, 'Pelargonium flower break virus'),
    ('NC_078016', 2808534, 'Pistacia emaravirus B'),
    ('NC_078017', 2808534, 'Pistacia emaravirus B'),
    ('NC_078018', 2808534, 'Pistacia emaravirus B'),
    ('NC_078019', 2808534, 'Pistacia emaravirus B'),
    ('NC_078020', 2808534, 'Pistacia emaravirus B'),
    ('NC_078021', 2808534, 'Pistacia emaravirus B'),
    ('NC_078022', 2808534, 'Pistacia emaravirus B'),
    ('NC_078023', 2808534, 'Pistacia emaravirus B'),
    ('NC_009992', 64048, 'Plum bark necrosis stem pitting-associated virus'),
]

# (本地 FASTA 常量, 登记 accession, 登记 species)
LOCAL = [
    (LOCAL_FASTA_PPV, 'PPV_MK387313_artificial',
     'Plum pox virus artificial strain'),
    (LOCAL_FASTA_PVY, 'PVY_artificial_strain',
     'Potato virus Y artificial strain'),
]

_EUTILS_HOST = 'eutils.ncbi.nlm.nih.gov'


def _safe_fetch_ncbi(path_and_query):
    """只允许 https + eutils 官方域名，解析 IP 拒绝私网/环回，防 SSRF。"""
    url = f'https://{_EUTILS_HOST}{path_and_query}'
    u = urllib.parse.urlparse(url)
    if u.scheme != 'https' or u.netloc != _EUTILS_HOST:
        raise ValueError(f'仅允许 eutils 官方域名: {url}')
    for info in socket.getaddrinfo(u.hostname, 443):
        ip = ipaddress.ip_address(info[4][0])
        if (ip.is_private or ip.is_loopback or ip.is_link_local
                or ip.is_reserved or ip.is_multicast):
            raise ValueError(f'解析到非公网地址 {ip}，拒绝请求')
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode('utf-8', 'replace')
        except Exception as e:
            print(f'  efetch 重试 {attempt+1}: {e}')
            time.sleep(5 * (attempt + 1))
    raise RuntimeError('NCBI efetch 连续失败')


def parse_fasta(text):
    recs, name, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if name:
                recs.append((name, ''.join(seq)))
            name, seq = line[1:], []
        elif line:
            seq.append(line)
    if name:
        recs.append((name, ''.join(seq)))
    return recs


def main():
    for p in (FA_P, INFO_P, MANIFEST_P):
        if not os.path.normpath(p).startswith(os.path.normpath(PLATFORM)):
            raise ValueError(f'落盘路径越界: {p}')
    os.makedirs(os.path.dirname(FA_P), exist_ok=True)
    accs = [a for a, _, _ in REFS]
    print(f'抓取 {len(accs)} 条 NCBI 参考序列 ...')
    text = _safe_fetch_ncbi(
        '/entrez/eutils/efetch.fcgi?' +
        urllib.parse.urlencode({'db': 'nucleotide', 'id': ','.join(accs),
                                'rettype': 'fasta', 'retmode': 'text'}))
    recs = parse_fasta(text)
    got = {r[0].split('.')[0] for r in recs}
    missing = [a for a in accs if a not in got]
    if missing:
        print(f'⚠ 未取到: {missing}')
    print(f'取到 {len(recs)} 条')

    meta = {a: (t, s) for a, t, s in REFS}
    cols = ['Accession', 'Taxid', 'Species_NCBI', 'Segment', 'Sequence_Type',
            'Molecule_Type2', 'Length']
    n = 0
    with open(FA_P, 'w', encoding='utf-8') as fa, \
         open(INFO_P, 'w', encoding='utf-8', newline='') as info:
        info.write('\t'.join(cols) + '\n')
        for name, seq in recs:
            acc_full = name.split()[0]
            acc = acc_full.split('.')[0]
            fa.write(f'>{acc_full} {name[len(acc_full):].strip()}\n{seq}\n')
            t, s = meta.get(acc, ('', name.split()[0]))
            info.write('\t'.join([acc_full, str(t), s, '', 'NCBI', '', str(len(seq))]) + '\n')
            n += 1
        # 本地人工株。TSWV 片段标签：NC_002052=L 8897 / NC_002050=M 4821 /
        # NC_002051=S 2916（与 NCBI 记录核对过，勿按 accession 顺序想当然）
        seg_fix = {'NC_002050': 'M', 'NC_002051': 'S', 'NC_002052': 'L'}
        for fp, acc, desc in LOCAL:
            rec = parse_fasta(open(fp, encoding='utf-8').read())[0]
            fa.write(f'>{acc} {desc}\n{rec[1]}\n')
            info.write('\t'.join([acc, '', desc, '', 'Artificial',
                                  'RNA', str(len(rec[1]))]) + '\n')
            n += 1
    # TSWV 片段标签修正
    lines = open(INFO_P, encoding='utf-8').read().splitlines()
    out = []
    for l in lines:
        p = l.split('\t')
        if p and p[0].split('.')[0] in seg_fix:
            p[3] = seg_fix[p[0].split('.')[0]]
            l = '\t'.join(p)
        out.append(l)
    open(INFO_P, 'w', encoding='utf-8', newline='').write('\n'.join(out) + '\n')
    print(f'库写入 {n} 条 -> {FA_P}')
    json.dump({'name': 'viromock_kv', 'reference': FA_P, 'ref_info': INFO_P,
               'engines': [], 'built_at': time.strftime('%Y-%m-%d %H:%M:%S')},
              open(MANIFEST_P, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)
    print('索引构建（需另行执行 known_virus_suite index：--engine salmon / minibwa 各一次，')
    print(f'  --reference "{FA_P}" --ref-info "{INFO_P}" --index-dir 库目录 --out 库目录')


if __name__ == '__main__':
    main()
