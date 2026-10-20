# -*- coding: utf-8 -*-
"""把 3rd/ 下的第三方工具打成插件 zip + 生成插件源索引（plugins-repo/index.json）。

产物：
    plugins-dist/<name>-plugin.zip     每工具一包，顶层 <name>/plugin.json + 工具本体
    plugins-repo/index.json            插件源索引（工具包 + 模块分组）

zip 本体不上传 git（体积），走 Release（tag: plugin-packs）分发，索引里记录
下载 URL 与 sha256，平台「插件管理 → 插件源」按索引自动检测缺失并一键下载解压。

模块层（MODULES）：把「工具」与「功能模块」解耦——同一个工具常被多个模块
共用（mafft 被 align/quicktree/kvchain 都用），所以**不按模块重复打包**，
而是索引里额外给一份 modules 分组；平台按模块安装时先算工具并集，再跳过
本机已有/已装的，只下载真正缺的（见 plugin_manager.plan_module_install）。

用法：
    python dev_tools/build_plugin_repo.py            # 打全部包
    python dev_tools/build_plugin_repo.py trimal     # 只打指定包（调试用）
"""
import hashlib
import json
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIST = os.path.join(ROOT, 'plugins-dist')
INDEX_DIR = os.path.join(ROOT, 'plugins-repo')

# Release 资产下载 URL 前缀（gh release create plugin-packs 上传后即生效）
RELEASE_BASE = ('https://github.com/zhangwenda0518/plant-virus-analysis-platform'
                '/releases/download/plugin-packs')
# 国内镜像站：插件 zip 放在该站 docs/plugins/ 下，由 nginx 静态服务。
# 环境变量 PVAP_PLUGIN_BASE 可覆盖（自建镜像/换域名时用）。
MIRROR_BASE = os.environ.get('PVAP_PLUGIN_BASE',
                             'http://39.106.101.94/plugins')
# Gitee Release（国内代码托管，附件单文件上限 100MB）
GITEE_BASE = ('https://gitee.com/zhang-wenda/plant-virus-analysis-platform'
              '/releases/download/plugin-packs')
GITEE_MAX_MB = 100        # 超过此大小的包不挂 Gitee（会 400 拒绝）

# 工具打包表：插件名 → (源路径(相对平台根), 版本, 描述,
#                     provides{Cfg.tools键: 包内相对路径})
# provides 键与 config.detect_tools() 的键一致；值是 zip 解出后插件目录内路径。
TOOLS = {
    'kunpeng': ('3rd/bin', '0.7.12', 'kunpeng 病毒核酸分类引擎（k-mer 精确版）',
                {'kunpeng': 'kun_peng-v0.7.12-x86_64-pc-windows-msvc.exe'}),
    'seqkit': ('3rd/bin', '2.x', 'FASTA/FASTQ 序列工具箱',
               {'seqkit': 'seqkit.exe'}),
    'crabz': ('3rd/bin', '1.x', '并行 gzip/pigz 压缩解压',
              {'crabz': 'crabz.exe'}),
    'aria2c': ('3rd/bin', '1.x', '多线程下载器（公共数据检索用）',
               {'aria2c': 'aria2c.exe'}),
    'minibwa': ('3rd/bin', '0.x', '轻量 BWA 比对（宿主去除用）',
                {'minibwa': 'minibwa.exe'}),
    'sracha': ('3rd/bin', '0.x', 'sracha 序列校正',
               {'sracha': 'sracha.exe'}),
    'blast': ('3rd/tools/Blast', '2.x', 'NCBI BLAST+ 套件（blastn/blastp/blastx）',
              {'blastn': 'bin/blastn.exe', 'makeblastdb': 'bin/makeblastdb.exe',
               'blastp': 'bin/blastp.exe', 'blastx': 'bin/blastx.exe'}),
    'spades': ('3rd/tools/SPAdes-Windows', '4.x', 'SPAdes 基因组从头组装',
               {'spades': 'bin/spades.bat'}),
    'mafft': ('3rd/tools/mafft-win', '7.x', 'MAFFT 多序列比对',
              {'mafft': 'mafft.bat'}),
    'fasttree': ('3rd/tools/FastTree', '2.1.x', 'FastTree 快速建树',
                 {'fasttree': 'FastTree.exe'}),
    'raxml-ng': ('3rd/tools/raxml-ng', '1.x', 'RAxML-NG 最大似然建树',
                 {'raxml-ng': 'raxml-ng.exe'}),
    'iqtree': ('3rd/tools/iQtree', '3.0.1', 'IQ-TREE 3 建树（v3）',
               {'iqtree3': 'iqtree-3.0.1-Windows/bin/iqtree3.exe'}),
    'trimal': ('3rd/tools/trimAl', '1.x', 'trimAl 比对清剪',
               {'trimal': 'trimal.exe'}),
    'diamond': ('3rd/tools/diamond', '2.x', 'DIAMOND 蛋白/核酸快速比对',
                {'diamond': 'diamond.exe'}),
    'mmseqs': ('3rd/tools/mmseqs', 'x', 'MMseqs2 超快序列搜索',
               {'mmseqs': 'bin/mmseqs.exe'}),
    'minimap2': ('3rd/tools/minimap2', '2.x', 'minimap2 长读长/参考比对',
                 {'minimap2': 'minimap2.exe'}),
    'samtools': ('3rd/tools/samtools', '1.x', 'SAM/BAM 处理套件',
                 {'samtools': 'bin/samtools.exe'}),
    'bcftools': ('3rd/tools/bcftools', '1.x', 'VCF/BCF 变异调用',
                 {'bcftools': 'bin/bcftools.exe'}),
    'salmon': ('3rd/tools/salmon2', '1.x', 'Salmon 转录本定量',
               {'salmon': 'salmon.exe'}),
    'fastp': ('3rd/tools/fastp', '0.x', 'fastp FASTQ 质控',
              {'fastp': 'fastp.exe'}),
    'table2asn': ('3rd/tools/table2asn', '1.x', 'NCBI table2asn 提交工具',
                  {'table2asn': 'table2asn.exe'}),
    'pandepth': ('3rd/tools/pandepth', '0.x', 'pandepth 覆盖度统计',
                 {'pandepth': 'pandepth.exe'}),
    'viral-consensus': ('3rd/tools/viral_consensus', '0.x', '病毒共识序列生成',
                        {'viral_consensus': 'viral_consensus.exe'}),
    'lsd2': ('3rd/tools/lsd2', '2.x', 'LSD2 最小进化定年',
             {'lsd2': 'lsd2.exe'}),
    'dsrnamax': ('3rd/tools/dsrnamax', '1.x', 'dsRNAmax siRNA 设计',
                 {'dsRNAmax': 'dsRNAmax_det.exe'}),
    'rdp5': ('3rd/tools/rdp5', '5.x', 'RDP5 重组检测（含 clustalw2/3seqTable）',
             {'RDP5CL': 'RDP5CL.exe'}),
    # gbdraw / orfipy 随 3rd/python 运行时分发（Scripts 入口），不单独成包：
    # 单拷 exe 依赖完整 Python 环境，拆出无法独立运行。
}

# 「硬编码路径」类工具：模块内直接拼 PLATFORM_ROOT/3rd/tools/<x>（不查
# Cfg.tools），装在 plugins/ 下永远找不到 → 用 install_to 还原到原路径。
# provides 的键仅供展示/状态判断，实际生效靠路径还原。
PATH_TOOLS = {
    'snpeff': ('3rd/tools/snpeff', '5.x', 'snpEff 变异注释（需 jre-snpeff）',
               {'snpeff': 'snpEff/snpEff.jar'}),
    'jre-snpeff': ('3rd/tools/jre-snpeff', '17', 'snpEff 专用精简 JRE',
                   {'jre-snpeff': 'bin/java.exe'}),
    'snpgenie': ('3rd/tools/snpgenie', '1.x', 'SNPGenie 选择压力分析',
                 {'snpgenie': 'snpgenie.pl'}),
    'srna-target': ('3rd/tools/srna-target', '1.x', 'siRNA 靶点预测三引擎（miranda/pita/rnahybrid）',
                    {'srna-target': 'miranda'}),
    'strawberry-perl': ('3rd/tools/strawberry-perl', '5.32', '便携 Perl（snpgenie/脚本依赖）',
                        {'strawberry-perl': 'perl/bin/perl.exe'}),
    'datasets': ('3rd/tools/datasets', '1.x', 'NCBI datasets 命令行',
                 {'datasets': 'datasets.exe'}),
}

# ----------------------------------------------------------------------
# 模块层：功能模块 → 依赖的插件（工具包）名
# 键名与 web/tools_api.py 的 TOOL_REGISTRY 对齐，另加管道阶段与专项页。
# 工具可被多个模块共用，模块安装时取并集后去重（不重复下载）。
# ----------------------------------------------------------------------
MODULES = {
    'preprocess': {
        'title': '⓪ 质控与宿主去除',
        'title_en': 'QC & host removal',
        'desc': 'fastp 质控 + seqkit 序列处理 + 宿主比对（minibwa）',
        'plugins': ['fastp', 'seqkit', 'minibwa'],
    },
    'assemble': {
        'title': '① 病毒组装',
        'title_en': 'Viral assembly',
        'desc': 'SPAdes 从头组装 + BLAST 骨架比对',
        'plugins': ['spades', 'blast'],
    },
    'identify': {
        'title': '② 病毒鉴定（kunpeng）',
        'title_en': 'Virus identification',
        'desc': 'kunpeng 核酸分类引擎 + seqkit 预处理',
        'plugins': ['kunpeng', 'seqkit'],
    },
    'contigs': {
        'title': '④ contig 分类',
        'title_en': 'Contig classification',
        'desc': 'kunpeng 分类 + 序列整理',
        'plugins': ['kunpeng', 'seqkit'],
    },
    'orf': {
        'title': '⑤ ORF 预测与注释',
        'title_en': 'ORF prediction & annotation',
        'desc': 'orfipy 预测 + BLAST 比对注释',
        'plugins': ['blast'],
    },
    'genome': {
        'title': '⑥ 基因组图谱',
        'title_en': 'Genome plots',
        'desc': 'gbdraw 出图（随 3rd/python 运行时，另装）；无则回退纯 Python 引擎',
        'plugins': [],
    },
    'align': {
        'title': '⑦ 多序列比对',
        'title_en': 'Multiple alignment',
        'desc': 'MAFFT 比对 + trimAl 清剪',
        'plugins': ['mafft', 'trimal'],
    },
    'tree': {
        'title': '⑧ 建树',
        'title_en': 'Tree building',
        'desc': 'RAxML-NG 最大似然 / FastTree 快速 / IQ-TREE 3',
        'plugins': ['raxml-ng', 'fasttree', 'iqtree', 'trimal'],
    },
    'kvsuite': {
        'title': '⑨ 已知病毒识别与定量',
        'title_en': 'Known virus suite',
        'desc': 'minimap2 比对 + samtools 处理 + snpEff 注释 + SNPGenie 选择压力',
        'plugins': ['minimap2', 'samtools', 'bcftools', 'snpeff',
                    'jre-snpeff', 'snpgenie', 'strawberry-perl', 'salmon'],
    },
    'rdp': {
        'title': 'RDP 重组分析',
        'title_en': 'Recombination (RDP5)',
        'desc': 'RDP5 重组检测（含 clustalw2 与 3seqTable）',
        'plugins': ['rdp5'],
    },
    'phylogeo': {
        'title': '系统地理分析',
        'title_en': 'Phylogeography',
        'desc': 'LSD2 定年 + 地理迁移推断',
        'plugins': ['lsd2'],
    },
    'sdt': {
        'title': 'SDT 精确分析',
        'title_en': 'SDT identity',
        'desc': 'BLAST 两两比对计算同一性',
        'plugins': ['blast'],
    },
    'verify': {
        'title': '候选序列验证',
        'title_en': 'Candidate verification',
        'desc': 'BLAST / DIAMOND / MMseqs2 三引擎交叉验证',
        'plugins': ['blast', 'diamond', 'mmseqs'],
    },
    'dsrna': {
        'title': 'dsRNA 设计',
        'title_en': 'dsRNA design',
        'desc': 'dsRNAmax siRNA 位点筛选',
        'plugins': ['dsrnamax'],
    },
    'crispr': {
        'title': 'CRISPR 检测设计',
        'title_en': 'CRISPR design',
        'desc': 'BLAST 特异性检查',
        'plugins': ['blast'],
    },
    'primer': {
        'title': '引物设计',
        'title_en': 'Primer design',
        'desc': 'BLAST 特异性检查',
        'plugins': ['blast'],
    },
    'mirna': {
        'title': 'miRNA 靶标预测',
        'title_en': 'miRNA target',
        'desc': 'miRanda / RNAhybrid / PITA 三引擎',
        'plugins': ['srna-target'],
    },
    'public_meta': {
        'title': '公共数据检索',
        'title_en': 'Public data search',
        'desc': 'NCBI datasets 命令行 + 多线程下载',
        'plugins': ['datasets', 'aria2c'],
    },
    'phylodyn': {
        'title': '进化动力学',
        'title_en': 'Evolutionary dynamics',
        'desc': 'BEAST 作业包导出与结果解析（R 由用户环境提供）',
        'plugins': ['blast'],
    },
    'misc': {
        'title': '通用工具',
        'title_en': 'General tools',
        'desc': '压缩/解压等基础设施',
        'plugins': ['crabz', 'aria2c', 'sracha', 'table2asn',
                    'viral-consensus', 'pandepth'],
    },
}


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _zip_add(z, arcname, src_path=None, data=None):
    """确定性写 zip：固定时间戳与权限位，保证同一份输入产出逐字节相同的包。

    为什么必须确定性：索引里记录每个包的 sha256，平台下载后校验。zip 默认
    把「当前时间」写进每个成员头，重打一次哈希就变，线上包与索引立刻失配，
    用户会看到"校验失败"（实测踩过）。固定成 1980-01-01（zip 纪元起点）
    后，内容不变则 sha256 不变，索引可以长期有效。
    """
    zi = zipfile.ZipInfo(arcname, date_time=(1980, 1, 1, 0, 0, 0))
    zi.compress_type = zipfile.ZIP_DEFLATED
    zi.external_attr = 0o644 << 16          # 普通文件，不带执行位
    if data is not None:
        z.writestr(zi, data)
    else:
        with open(src_path, 'rb') as f:
            z.writestr(zi, f.read())


def build(name, src_rel, version, desc, provides, install_to=None):
    src = os.path.join(ROOT, src_rel)
    if not os.path.isdir(src):
        print(f'  ✘ 源目录不存在，跳过 {name}: {src}')
        return None
    out_zip = os.path.join(DIST, f'{name}-plugin.zip')
    manifest = {'name': name, 'version': version, 'kind': 'tool',
                'description': desc, 'provides': provides}
    if install_to:
        # 内容装到平台原路径（不落 plugins/），给硬编码路径的工具用
        manifest['install_to'] = install_to
    # 3rd/bin = 自包含单 exe 清单（kunpeng/seqkit/crabz/aria2c/minibwa/sracha）：
    # 只装目标文件，不把整个 bin 目录（其他插件的目标）顺带进去。
    only_bin = os.path.normpath(src_rel) == os.path.normpath('3rd/bin')
    wanted = set(provides.values())
    members = []
    for cur, _dirs, files in os.walk(src):
        for fn in files:
            p = os.path.join(cur, fn)
            rel = os.path.relpath(p, src).replace(os.sep, '/')
            if only_bin and rel not in wanted:
                continue
            if '__pycache__' in rel or fn.endswith(('.pyc', '.log')):
                continue
            members.append((f'{name}/{rel}', p))
    members.sort()                       # 固定顺序（os.walk 顺序依赖文件系统）
    with zipfile.ZipFile(out_zip, 'w') as z:
        _zip_add(z, f'{name}/plugin.json',
                 data=json.dumps(manifest, ensure_ascii=False, indent=1))
        for arc, p in members:
            _zip_add(z, arc, src_path=p)
    size = os.path.getsize(out_zip)
    print(f'  ✔ {name}-plugin.zip  {size/1048576:.1f}MB')
    # 下载源按可用性排序：自建镜像 → Gitee（仅 <100MB 的包）→ GitHub Release
    urls = [f'{MIRROR_BASE}/{name}-plugin.zip']
    if size < GITEE_MAX_MB * 1048576:
        urls.append(f'{GITEE_BASE}/{name}-plugin.zip')
    urls.append(f'{RELEASE_BASE}/{name}-plugin.zip')
    entry = {'name': name, 'version': version, 'kind': 'tool',
             'description': desc, 'provides': provides,
             'file': f'{name}-plugin.zip',
             'url': urls[0],
             'urls': urls,
             'size': size, 'sha256': sha256_file(out_zip)}
    if install_to:
        entry['install_to'] = install_to
    return entry


def main():
    only = set(sys.argv[1:])
    os.makedirs(DIST, exist_ok=True)
    os.makedirs(INDEX_DIR, exist_ok=True)
    entries = []
    for name, spec in TOOLS.items():
        if only and name not in only:
            continue
        r = build(name, *spec)
        if r:
            entries.append(r)
    for name, (src_rel, version, desc, provides) in PATH_TOOLS.items():
        if only and name not in only:
            continue
        # install_to = 源目录相对平台根的路径（原样还原）
        r = build(name, src_rel, version, desc, provides,
                  install_to=src_rel)
        if r:
            entries.append(r)
    if only:
        return
    # 模块层：只做「模块 → 插件名」引用，不重复打包。校验引用完整性，
    # 打错的模块名会让用户装出半残环境，必须在生成索引时就报错。
    known = {e['name'] for e in entries}
    modules = {}
    for mid, m in MODULES.items():
        unknown = [p for p in m['plugins'] if p not in known]
        if unknown:
            raise SystemExit(f'[MODULES] {mid} 引用了不存在的插件: {unknown}')
        modules[mid] = {'title': m['title'], 'title_en': m['title_en'],
                        'desc': m['desc'], 'plugins': list(m['plugins'])}
    index_path = os.path.join(INDEX_DIR, 'index.json')
    with open(index_path, 'w', encoding='utf-8') as f:
        json.dump({'schema': 2, 'source': 'PVAP plugin repo',
                   'release_tag': 'plugin-packs',
                   'plugins': entries, 'modules': modules},
                  f, ensure_ascii=False, indent=1)
    total = sum(e['size'] for e in entries)
    # 共享工具统计：被多个模块引用的插件只打一个包
    refs = {}
    for m in modules.values():
        for p in m['plugins']:
            refs[p] = refs.get(p, 0) + 1
    shared = {k: v for k, v in sorted(refs.items(), key=lambda x: -x[1])
              if v > 1}
    print(f'\n索引: {index_path}（{len(entries)} 个插件, '
          f'{len(modules)} 个模块, zip 合计 {total/1048576:.0f}MB）')
    if shared:
        print('共享工具（多模块引用，只打一个包）: '
              + ', '.join(f'{k}×{v}' for k, v in shared.items()))


if __name__ == '__main__':
    main()
