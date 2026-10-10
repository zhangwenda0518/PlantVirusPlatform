# -*- coding: utf-8 -*-
"""平台级自测：所有页面 200 + 关键 API 结构 + 页面渲染含示例按钮/示例文件
+ 既有样品报告与结果中心数据。只读，不启动分析任务。"""
import json
import os
import re
import shutil
import sys
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Virus_Platform_Core.config import PLATFORM_ROOT  # noqa: E402

EX = os.path.join(PLATFORM_ROOT, 'examples')

PAGES = ['/', '/pipeline', '/samples', '/hostremoval',
         '/orf', '/annotation', '/genome', '/primer', '/build', '/results',
         '/tools', '/settings', '/meta', '/download', '/logan', '/submit',
         '/virome', '/tools?g=virus', '/tools?g=annotate', '/tools?g=compare',
         '/tools?g=detect']


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg, flush=True)
    assert cond, msg


def main():
    import app as appmod  # noqa: E402
    c = appmod.app.test_client()

    # ---------- 1. 页面 ----------
    print('--- 页面渲染 ---', flush=True)
    for p in PAGES:
        r = c.get(p)
        check(r.status_code == 200, f'GET {p} -> {r.status_code}')

    # 页面内含示例按钮（本次补充的示例入口都要真实渲染出来）
    html = c.get('/tools?g=compare').get_data(as_text=True)
    for marker, where in [('seqPrepExample()', '参考序列获取卡'),
                          ('ictvPreview()', 'ICTV 选参预览'),
                          ('id="spCascade"', 'ICTV 级联下拉容器'),
                          ('tbBuild(this)', '进化树构建·集合建树'),
                          ('id="tbSource"', '建树来源类型'),
                          ('id="alSource"', '比对来源类型'),
                          ('id="al_fa_extra"', '比对卡额外序列'),
                          ('alignRun(this)', '序列比对·运行'),
                          ('alignEditToggle()', '比对查看器·编辑模式'),
                          ('alignSave(this)', '比对查看器·保存副本'),
                          ("fillExample('al_fa', EXAMPLE_SET_FASTA)", '比对卡示例'),
                          ('id="alSample"', '比对查看器·样品/集合直选'),
                          ('id="alGroup"', '比对查看器·分组直选'),
                          ('tbSendToAlign()', '建树卡·跳转看建树比对'),
                          ('fillExample(\'qt_fa\', EXAMPLE_SET_FASTA)', 'FASTA 建树'),
                          ('fillExample(\'tv_file\', EXAMPLE_TREE_NWK)', '树查看'),
                          ('fillExample(\'sd_fa\', EXAMPLE_SET_FASTA)', 'SDT 卡'),
                          ('fillExample(\'s_files\', EXAMPLE_GB_TRIO)', '导入示例 .gb')]:
        check(marker in html, f'tools 页渲染含 {where}')
    # RDP5 重组卡输入通道（2026-09-15 补齐）：📋 粘贴 + ✨ 同种病毒示例。
    # RDP5 是「同种/同株系」比较，示例必须是同一病毒的比对（跨病毒比较在比对卡），
    # 故示例常量与比对卡的 EXAMPLE_SET_FASTA 分开——别合并，语义不同。
    for marker, where in [("pasteSeq('rdp_input', '.fasta')", 'RDP5 卡·粘贴序列'),
                          ("fillExample('rdp_input', EXAMPLE_RECOMB_FASTA)", 'RDP5 卡·同种病毒示例'),
                          ('data-i18n="tk.rdpScopeHint"', 'RDP5 卡·同种病毒范围提示'),
                          ('id="rdp_qc"', 'RDP5 卡·序列预筛开关'),
                          ('id="rdp_qc_ref"', 'RDP5 卡·预筛参考序列'),
                          ('id="rdp_qc_len"', 'RDP5 卡·预筛长度比'),
                          ('data-i18n="tk.rdpQcHint"', 'RDP5 卡·预筛判据说明')]:
        check(marker in html, f'tools 页渲染含 {where}')
    # 比较组三模块结构：参考序列获取 / 进化树构建 / SDT
    for gone in ('id="t-contigs-struct"', 'id="t-ncbi"', 'id="t-synteny-gb"',
                 'id="t-msa"', 'id="t-tree"', 'id="s_style"', 'id="s_lovis4u"',
                 'id="ms_fa"', 'id="msBox"'):
        check(gone not in html, f'旧结构已移除: {gone}')
    for present in ('id="t-seqprep"', 'id="t-treebuild"', 'id="t-sdt"'):
        check(present in html, f'新模块存在: {present}')
    # 前端按域拆为 app*.js（app.js 核心 + app-compare/runhistory/jobs/
    # browse/batch/storage/compare）：定义检查对全部文件取并集，
    # 后续再拆无需改这里
    _www = os.path.join(PLATFORM_ROOT, 'webapp', 'static')
    import glob as _glob
    js = ''.join(open(_f, encoding='utf-8').read()
                 for _f in sorted(_glob.glob(os.path.join(_www, 'app*.js'))))
    for fn in ('ictvCascadeRefetch', 'ictvPreview', 'ictvDownload', 'loadTbColls',
               'tbBuild', 'gbBuildTree', 'alignRun', 'alignLoad', 'alignRender',
               'alignEditToggle', 'alignSave', 'alignSend', 'tbSourceChanged',
               'alSourceChanged', 'alFillFromSource', 'loadAlColls'):
        check(f'function {fn}' in js or f'async function {fn}' in js,
              f'app.js/app-compare.js 定义 {fn}')
    # 旧项目归档：样品/集合列表不再出现下划线与归档项
    r = c.get('/api/samples')
    names = [x['name'] for x in r.get_json()]
    check(all(not n.startswith('_') for n in names), '样品列表无下划线项')
    # 归档样品接口：数据来自 results/_archive/ 下的历史项目。
    # 原先断言"列表必须非空"，等于把测试绑在开发机的残留数据上——清理 run/
    # 之后必然 FAIL（实测）。改为自造一个临时归档样品 → 断言接口如实反映
    # → finally 删除，测试自带 fixture、不依赖环境残留。
    from Virus_Platform_Core.config import DIRS
    arch = os.path.join(DIRS['results'], '_archive')
    probe = os.path.join(arch, '_it_platform_arch')
    os.makedirs(os.path.join(probe, '07_report'), exist_ok=True)
    with open(os.path.join(probe, '07_report', 'report.html'), 'w',
              encoding='utf-8') as f:
        f.write('<html>probe</html>')
    try:
        r = c.get('/api/samples/archived')
        items = r.get_json() if r.status_code == 200 else []
        names = [x['name'] for x in items]
        check(r.status_code == 200 and '_it_platform_arch' in names,
              f'归档样品列表可用（{len(items)} 个，含自造 fixture）')
        hit = next((x for x in items if x['name'] == '_it_platform_arch'), {})
        check(hit.get('has_report') is True, '归档项正确标识含报告')
    finally:
        shutil.rmtree(probe, ignore_errors=True)
    r = c.get('/api/gb/collections')
    check(all(not x['name'].startswith('_') for x in r.get_json()),
          '集合列表无下划线项')
    html = c.get('/annotation').get_data(as_text=True)
    check("fillExample('oa_fa')" in html, '注释页渲染含示例按钮')
    html = c.get('/genome').get_data(as_text=True)
    check('EXAMPLE_GENBANK_GB' in html, '图谱页渲染含示例 GenBank 按钮')
    # 预览缩放：linear 图原始宽 2000+px（gbdraw canvas.linear.width 固定），
    # 侧栏直接显示必横向滚很长，故默认「适应宽度」并给 100%/150% 切换。
    check('id="gp_zoom"' in html and "gpZoom('fit')" in html
          and 'data-zoom="1.5"' in html,
          '图谱页预览含缩放控制（适应宽度 / 100% / 150%）')
    html = c.get('/primer').get_data(as_text=True)
    # /primer 已重构为交互式 primer3 设计器（原「模式感知示例」随旧页面下线），
    # 断言改为新页面的示例入口：填充 pd_fa + 复用内置示例常量。
    check("fillExample('pd_fa'" in html and 'EXAMPLE_FASTA' in html,
          '引物页渲染含示例按钮')
    js = open(os.path.join(PLATFORM_ROOT, 'webapp', 'static', 'app.js'),
              encoding='utf-8').read()
    for const in ('EXAMPLE_FASTA', 'EXAMPLE_SET_FASTA', 'EXAMPLE_TREE_NWK',
                  'EXAMPLE_GENBANK_GB', 'EXAMPLE_GB_TRIO',
                  'EXAMPLE_CONSERVED_FASTA'):
        check(const in js, f'app.js 定义 {const}')

    # 示例文件齐全
    need = ['example_viral_contigs.fasta', 'example_virus_set.fasta',
            'example_conserved_set.fasta', 'example_tree.nwk',
            'example_genome.gb'] + \
        [f'example_synteny_{x}.gb' for x in 'ABC']
    for fn in need:
        check(os.path.isfile(os.path.join(EX, fn)), f'示例文件存在: {fn}')

    # ---------- 2. 关键 API ----------
    print('--- 关键 API ---', flush=True)
    for api in ('/api/tools', '/api/dbs', '/api/samples', '/api/queue',
                '/api/ictv/cascade',
                '/api/settings', '/api/tool/runs', '/api/msa/samples',
                '/api/ncbi/collections', '/api/gb/collections',
                '/api/submit/tables', '/api/submit/contig_runs',
                '/api/meta/collections', '/api/logan/jobs',
                '/api/dl/batches', '/api/tool/viral_contigs?run=x'):
        r = c.get(api)
        check(r.status_code in (200, 400), f'GET {api} -> {r.status_code}'
              + ('（空参数 400 属预期）' if r.status_code == 400 else ''))

    # ---------- 3. 既有样品报告 ----------
    print('--- 样品报告 ---', flush=True)
    _smp = os.path.join(DIRS['results'], 'it_platform_smp')
    os.makedirs(_smp, exist_ok=True)
    # 与上面归档接口同一口径：不绑开发机残留数据（清理 run/ 后必然 FAIL，
    # 2026-09-29 实测又踩中一次）。任意非下划线目录即被 /api/samples 收录，
    # 自造最小 fixture → 断言列表如实反映 → finally 删除。
    try:
        r = c.get('/api/samples')
        samples = [x['name'] for x in r.get_json()]
        check('it_platform_smp' in samples,
              f'样品列表如实反映（{len(samples)} 个，含自造 fixture）')
        for s in samples:
            r = c.get(f'/report/{s}/')
            if r.status_code == 200:
                check(True, f'样品 {s} 报告可访问')
                break
        else:
            print('  SKIP（无样品带报告）')
    finally:
        shutil.rmtree(_smp, ignore_errors=True)

    # ---------- 3a. contigs 卡：Krona 前置 / 宿主筛选与自动预测 / ID 复制 / 跳转 ----
    hv = c.get('/tools?g=virus').get_data(as_text=True)
    # Krona 旭日块已并入分类报告折叠区（summary=分类旭日图（可下钻））。
    # 报告区顺序（2026-09-09 用户确认）：分类表 → contig 明细 → 宿主预测统计 → 旭日图
    # 注：宿主预测统计在 contig 明细下方是用户明确要求（DEVELOPMENT_NOTES 五十五）；
    #     旭日图并入折叠区后排在最后，故断言为 det < host < sun。
    _i_host = hv.find('id="hostRptSec"')
    _i_det = hv.find('病毒序列分类（contig 明细）')
    _i_sun = hv.find('id="sunC"')
    check(0 < _i_det < _i_host < _i_sun,
          '报告区顺序：分类表 → contig 明细 → 宿主预测统计 → 旭日图')
    # 卡片标签现在带 data-page-node-id 等附加属性，按 id 匹配而不是整串比对
    check(len(re.findall(r'<section[^>]*\bid="t-contigs"', hv)) == 1
          and len(re.findall(r'<section[^>]*\bid="t-assemble"', hv)) == 1,
          '④ contigs / ③ assemble 卡存在且唯一')
    check(hv.count('id="hostRptSec"') == 1, '宿主预测统计块唯一（无错插副本）')
    for marker, where in [("vEx('tmv')", '示例病毒 TMV'),
                          ("vEx('pstvd')", '示例病毒 PSTVd'),
                          ("vEx('mix')", '示例病毒 Mix All')]:
        check(marker in hv, f'病毒识别组渲染含 {where}')
    for fn in ('example_tmv.fasta', 'example_pvy.fasta', 'example_cmv.fasta',
               'example_pstvd.fasta', 'example_mix.fasta'):
        check(os.path.isfile(os.path.join(EX, fn)), f'Metabuli 示例病毒存在: {fn}')
    check('id="anaRun"' not in hv and 'id="anaRunLabel"' in hv,
          '分类表「选择运行」下拉已去除（自动跟随最新运行标签）')
    for marker, where in [('id="vcHostF"', '宿主筛选下拉'),
                          ('autoHostIfMissing', '宿主预测自动运行（无宿主列时自动触发）'),
                          ('copyContigSeq(', 'ID 点击复制序列'),
                          ('anaJump(', '四件套跳转注释分析'),
                          ('processAnaJump', '跳转自动续接')]:
        check(marker in hv, f'病毒识别组渲染含 {where}')
    ha = c.get('/tools?g=annotate').get_data(as_text=True)
    check(all(f"'{a}'" in ha for a in ('blastn', 'blastx', 'primer')),
          't-hom 卡含 BLASTN / BLASTX / Primer 分析')
    for marker in ('Best E-value', 'Primer design complete', 'ncbiNuc'):
        check(marker in ha, f'metabuli 风格结果渲染含 {marker}')
    # 2026-09-29：CDD 卡已随本地 CDD 库下线，渲染里不应再有 CDD 分支
    check('ncbiCdd' not in ha and 'Conserved Domains' not in ha,
          'CDD 渲染分支已移除（cdd → t-cdd 卡不存在）')

    # 2026-09-18：引物设计由「病毒注释分析」拆到新组「检测与防治」。
    # 2026-09-29：dsRNA 设计曾整体迁出；2026-10-03 已恢复（引擎三模块 +
    # t-dsrna 工具卡 + /api/tool/dsrna_panel）——本段断言恢复后的渲染产物。
    # 断言的是**渲染出来的界面产物**（组顶栏），不是照抄 NAV_GROUPS 自证 ——
    # 这层最容易漏改的正是 pages.py 的 _PATH_TO_GROUP（漏了则 /primer 侧栏
    # 仍列注释组项、入口还指向 ?g=annotate）。
    hp = c.get('/primer').get_data(as_text=True)
    _gt = hp.split('class="group-topbar"', 1)[-1].split('</nav>', 1)[0]
    check(_gt and '检测与防治' in _gt, '/primer 组顶栏 =「检测与防治」')
    check('病毒注释分析' not in _gt, '/primer 组顶栏不再是「病毒注释分析」')
    hd = c.get('/tools?g=detect').get_data(as_text=True)
    check('id="t-dsrna"' in hd, 'dsRNA 卡已恢复（detect 组渲染 t-dsrna）')
    check('/api/tool/dsrna_panel' in hd, 'dsRNA 面板接口已随模块恢复接线')
    from Virus_Platform_Core.web.pages import NAV_GROUPS as _NG
    _g = {g['id']: g for g in _NG}
    _d_ids = {it.get('id') or it.get('href') for it in _g['detect']['items']}
    _a_ids = {it.get('id') or it.get('href') for it in _g['annotate']['items']}
    # 2026-10-06：detect 组新增 t-crispr（CRISPR/Cas12a 检测试剂设计卡），
    # 断言清单随 NAV_GROUPS 更新（旧清单缺 t-crispr 导致批量回归误报）
    check(_d_ids == {'/primer', 't-crispr', 't-dsrna', 't-mirna'},
          f'检测与防治组 = 引物设计 + CRISPR + dsRNA + miRNA（实际 {sorted(_d_ids)}）')
    check(not ({'/primer'} & _a_ids),
          '注释组不再含引物设计')

    # ---------- 3b. 模块历史运行组件（折叠/衔接/删除） ----------
    import os as _os
    # 工具工作台各组页渲染的是同一份 tools.html（卡片由前端按组显隐），
    # 故历史容器数直接从模板推导，避免硬编码数量随卡片增删而漂移。
    with open(_os.path.join(PLATFORM_ROOT, 'webapp', 'templates',
                            'tools.html'), encoding='utf-8') as _f:
        _n_rh = _f.read().count('class="rh" id="rh-')
    for pg, n in [('/tools?g=sample', _n_rh), ('/tools?g=compare', _n_rh),
                  ('/hostremoval', 1), ('/orf', 1),
                  ('/annotation', 1), ('/genome', 1), ('/primer', 1)]:
        h = c.get(pg).get_data(as_text=True)
        check(h.count('class="rh" id="rh-') == n, f'{pg} 历史容器 {n} 个')
    _os.makedirs(_os.path.join(PLATFORM_ROOT, 'run', 'tool_runs', 'it_rh_chk'),
                 exist_ok=True)
    r = c.post('/api/tool/runs/it_rh_chk/delete')
    check(r.status_code == 200
          and not _os.path.isdir(_os.path.join(PLATFORM_ROOT, 'run', 'tool_runs',
                                               'it_rh_chk')),
          '运行目录删除 API')
    r = c.post('/api/tool/runs/no_such_run/delete')
    check(r.status_code == 400, '删除不存在运行 400')

    # ---------- 3c. 导航项 ↔ 卡片一一对应 ----------
    # showModule() 只切换 <main> 的直接子元素（:scope > *）。若某个导航项
    # 对应的 id 不是 main 的直接子元素（例如被嵌在别的卡片里），点它就是
    # 空白页——2026-09-09 的 t-kvchain 正是如此（曾是 t-kvsuite 内的 div）。
    from html.parser import HTMLParser as _HP
    from Virus_Platform_Core.web.pages import NAV_GROUPS as _NAV

    class _MainKids(_HP):
        def __init__(self):
            super().__init__()
            self.depth = 0
            self.main_depth = None
            self.ids = set()

        def handle_starttag(self, tag, attrs):
            d = dict(attrs)
            if tag == 'main':
                self.main_depth = self.depth
            if (self.main_depth is not None
                    and self.depth == self.main_depth + 1 and d.get('id')):
                self.ids.add(d['id'])
            if tag not in ('br', 'img', 'input', 'meta', 'link', 'hr'):
                self.depth += 1

        def handle_endtag(self, tag):
            if tag not in ('br', 'img', 'input', 'meta', 'link', 'hr'):
                self.depth -= 1
            if tag == 'main':
                self.main_depth = None

    _mk = _MainKids()
    _mk.feed(hv)
    _want = {it['id'] for g in _NAV for it in g['items'] if it.get('id')}
    _missing = sorted(_want - _mk.ids)
    check(not _missing,
          f'导航项均有 <main> 直接子卡片（缺失: {_missing or "无"}）')

    # ---------- 3c-2. ⚡一键分析（全流程）：自有输入 + 参数 + 输出区 ----------
    # 2026-09-11 用户实测：该卡只有「运行 / 最近汇总」两个按钮——参数读的是
    # t-kvsuite 卡的 kv_* 字段（在别的模块页，看不到也改不了），也没有 .toolrun
    # 输出容器，日志与产物无处落，整页看起来「没有输入、没有输出、没有参数」。
    _kvchain = hv.split('id="t-kvchain"', 1)[-1].split('</section>', 1)[0]
    for _id, _what in (('kvc_samples', '样品输入'),
                       ('kvc_reads_r1', '直填测序数据 R1'),
                       ('kvc_reads_r2', '直填测序数据 R2'),
                       ('kvc_min_cov', '过滤阈值参数'), ('kvc_variant_qual', '变异阈值参数'),
                       ('kvc_evo', '扩展变异开关'), ('kvc_threads', '线程数'),
                       ('toolrun-kvchain', '输出区'),
                       ('kvcIdentifyTable', '鉴定结果表'),
                       ('kvcFilteredTable', '过滤结果表')):
        check(f'id="{_id}"' in _kvchain, f'一键全流程卡含{_what}（{_id}）')
    check('runKvchain(this)' in _kvchain and 'kvchainPullKvsuite()' in _kvchain,
          '一键全流程卡含运行按钮与「取该卡参数」按钮')
    # 参数必须来自本卡 kvc_* 表单：runKvchain 不得再读 kvsuite 卡的 kv_samples
    _runjs = hv.split('async function runKvchain(btn)', 1)[-1].split(
        '/* 一键全流程卡', 1)[0]
    check("val('kvc_samples')" in _runjs and "val('kv_samples')" not in _runjs,
          'runKvchain 读本卡 kvc_* 参数（不再借用 t-kvsuite 的字段）')
    # 两卡渲染互不覆盖：kvsuite 结果用 kv* 容器，一键全流程用 kvc* 容器
    check("const KV_IDS = {" in hv and "const KVC_IDS = {" in hv
          and "renderKvsuite(d, KVC_IDS)" in hv,
          '一键全流程内联复用 kvsuite 渲染（KVC_IDS 独立容器）')
    # 实时日志进卡片：任务标签表要有这两张卡的键（缺了输出区永远为空）。
    # TOOL_LABELS 定义在 app-jobs.js（2026-09-13 拆分），对全部 app*.js 取并集
    import glob as _glob
    _appjs = ''.join(
        open(_f, encoding='utf-8').read()
        for _f in sorted(_glob.glob(_os.path.join(
            PLATFORM_ROOT, 'webapp', 'static', 'app*.js'))))
    _tl = _appjs.split('const TOOL_LABELS = {', 1)[-1].split('};', 1)[0]
    for _k, _lbl in (('kvchain', '病毒定量与共识·一键'),
                     ('kvsuite', '已知病毒识别与定量')):
        check(f'{_k}:' in _tl and _lbl in _tl,
              f'TOOL_LABELS 登记 {_k}（卡片输出区可收实时日志）')

    # ---------- 3d. 基因组图谱：gbdraw 两子命令参数差异 ----------
    # gbdraw 0.14 的 circular / linear 参数集与取值集都不同，而前端 gb_opts
    # 是面向圈图语义的一套键（labels / track_type / species / label_placement…），
    # _run_gbdraw 必须按子命令翻译，否则 linear 报 unrecognized arguments
    # （退出码 2）→ 预览 / 出图 500。2026-09-09 用户实测踩到。
    try:
        import tempfile as _tf
        from Virus_Platform_Core.gbdraw_plot import _run_gbdraw as _rg
        _gopts = {'labels': 'out', 'track_type': 'tuckin', 'species': 'T',
                  'strain': 'S', 'feature_width': 20,
                  'multi_record_canvas': True, 'gc_content_width': 300,
                  'gc_skew_radius': 300, 'no_gc': True, 'no_skew': True,
                  'label_placement': 'horizontal'}
        _gp = _os.path.join(_tf.gettempdir(), 'vp_it_gbdraw', 'g')
        _made = _rg(gbk=_os.path.join(EX, 'example_genome.gb'),
                    out_prefix=_gp, mode='both', opts=_gopts)
        check(len(_made) == 2,
              f'gbdraw 圈图+线图双模式出图（前端全量参数），实际 {len(_made)} 张')
        with open(_made[0], encoding='utf-8', errors='replace') as _f:
            _svg = _f.read()
        check('hypothetical protein' in _svg,
              'gbdraw 出图含基因名（labels 翻译生效）')
    except Exception as _e:
        check(False, f'gbdraw 双模式出图失败: {_e}')

    # ---------- 3e. 内置示例体系 ----------
    # 所有模块都应能一键填示例；示例文件必须齐全，且与 UI 文案一致
    # （UI 写「CMV RNA1-3」「四病毒混合 6 条」，文件就必须是 3 段 / 6 条）。
    _exjs = open(_os.path.join(PLATFORM_ROOT, 'webapp', 'static',
                               'examples.js'), encoding='utf-8').read()
    check('CARD_EXAMPLES' in _exjs or 'var MAP' in _exjs,
          'examples.js 定义卡片示例映射')
    check('examples.js' in hv, '/tools 接入 examples.js')
    for _fn in ('example_R1.fastq.gz', 'example_R2.fastq.gz'):
        check(_os.path.isfile(_os.path.join(EX, _fn)), f'示例测序数据存在: {_fn}')
    from Virus_Platform_Core.utils import iter_fasta as _if
    _n_cmv = len(list(_if(_os.path.join(EX, 'example_cmv.fasta'))))
    _n_mix = len(list(_if(_os.path.join(EX, 'example_mix.fasta'))))
    check(_n_cmv == 3, f'CMV 示例含三分体（实际 {_n_cmv} 段）')
    check(_n_mix == 6, f'Mix All 示例含 6 条记录（实际 {_n_mix} 条）')

    # 示例结果 API（只读 examples/results/，与真实结果隔离）
    _r = c.get('/api/examples')
    _ex = _r.get_json() or []
    check(_r.status_code == 200 and len(_ex) >= 10,
          f'示例结果清单可用（{len(_ex)} 个模块）')
    _r = c.get('/api/examples/orf')
    _files = (_r.get_json() or {}).get('files', []) if _r.status_code == 200 else []
    check(len(_files) > 0, f'示例结果详情可读（orf {len(_files)} 个产物）')
    _r = c.get('/api/examples/orf/03_assembly/summary.json')
    check(_r.status_code == 200, '示例结果单文件可取')

    # 在线工具（BLASTN·BLASTX）与 LOGAN 的示例也要有输入与产物
    # （CDD 卡已于 2026-09-29 下线，示例映射同步移除）
    check(_os.path.isfile(_os.path.join(EX, 'example_contig_1.fasta')),
          '单 contig 示例文件存在（BLAST 粘贴用）')
    check("homSeqText: 'CONTIGS_TEXT'" in _exjs and
          "cddSeqText" not in _exjs,
          't-hom 示例填入序列正文（textarea），t-cdd 映射已移除')
    check("'t-hom': 'hom'" in _exjs and
          "'/logan': 'logan'" in _exjs,
          'RESULT_MAP 覆盖 hom / logan')
    _r = c.get('/api/examples/logan')
    _lgan = (_r.get_json() or {}).get('files', []) if _r.status_code == 200 else []
    check(any(f['path'].endswith('.html') and f['kind'] == 'html'
              for f in _lgan),
          f'LOGAN 示例报告可按 html 内联预览（{len(_lgan)} 个产物）')
    _r = c.get('/api/examples/sdt')
    _sdt = (_r.get_json() or {}).get('files', []) if _r.status_code == 200 else []
    check(any(f['path'].endswith('.pdf') and f['kind'] == 'pdf' for f in _sdt)
          and any(f['path'].endswith('.png') and f['kind'] == 'image'
                  for f in _sdt),
          'PDF 单列 pdf（内嵌预览）且 PNG 仍为 image')
    _man = json.load(open(_os.path.join(EX, 'results', 'manifest.json'),
                          encoding='utf-8'))
    _missing = [m for m in _man
                if not _os.path.isdir(_os.path.join(EX, 'results', m))]
    check(not _missing, f'示例结果清单与目录一致（{len(_man)} 个模块）')
    # synteny 模块已随 2026-09 重构移除，示例结果同步下线
    check(len(_man) >= 25, f'示例结果覆盖 ≥25 个模块（实际 {len(_man)}）')

    # ---------- 3f. 本地比对引擎（离线：blastn / DIAMOND / BFVD） ----------
    from Virus_Platform_Core.local_search import engine_status as _eng_status
    _st = _eng_status()
    check(bool(_st.get('ok')), '本地 blastn 就绪（病毒参考核酸库 + BLAST+）')
    check(bool(_st.get('ok_blastx')), '本地 DIAMOND blastx 就绪（viral_prot.dmnd）')
    check(bool(_st.get('ok_bfvd')), 'BFVD 病毒域库就绪（misc_db/suvtk/bfvd + mmseqs2）')
    check('hom_engine' in hv and 'value="local" selected' in hv,
          '前端可切换本地 / 在线引擎（hom_engine）')
    check('oa_lib_bfvd' in c.get('/annotation').get_data(as_text=True)
          and 'vf_methods' in hv and 'vc_methods' in hv,
          '前端含 BFVD 注释层勾选框与验证方法下拉（blastx/BFVD 双路）')
    _seq1 = ''.join(
        ln.strip() for ln in open(_os.path.join(EX, 'example_contig_1.fasta'),
                                  encoding='utf-8')
        if ln.strip() and not ln.startswith('>'))
    _r = c.get('/api/tool/analysis',
               query_string={'contig': 'example_contig_1', 'action': 'blastn',
                             'engine': 'local', 'seq': _seq1})
    _bl = _r.get_json() or {}
    check(_r.status_code == 200 and _bl.get('engine') == 'local' and
          len(_bl.get('hits') or []) > 0,
          f"本地 blastn 结果可经 API 读取（{len(_bl.get('hits') or [])} 条命中）")
    _p = _os.path.join(EX, 'results', 'hom', 'example_contig_1_blastx_local.json')
    check(_os.path.isfile(_p), f'本地示例产物存在: hom/example_contig_1_blastx_local.json')
    _d = json.load(open(_p, encoding='utf-8'))
    check(_d.get('engine') == 'local' and len(_d.get('hits') or []) > 0,
          f"hom/example_contig_1_blastx_local.json: engine=local 且有命中"
          f"（{len(_d.get('hits') or [])} 条）")

    # ---------- 4. 静态资源 ----------
    r = c.get('/static/app.js')
    check(r.status_code == 200, 'app.js 静态资源 200')
    r = c.get('/static/app-compare.js')
    check(r.status_code == 200, 'app-compare.js 静态资源 200')
    r = c.get('/static/i18n.js')
    check(r.status_code == 200, 'i18n.js 静态资源 200')

    # ---------- 5. Explorer 路由：2026-09-17 摘除入口，2026-09-27 彻底归档
    # （archive/_retire_20260927/explorer/）→ 断言 404 防回归 ----------
    # 背景：病毒浏览器在 2026-09-17 起在两个平台停用（只摘入口与蓝图注册，
    # 代码与数据留着）。这里原先逐个断言「路由存在且不是 404」，蓝图一摘就
    # 全变 404，`check()` 里的 assert 让**文件在这一节就崩**、后面所有检查跑不到
    # ——看起来像"本次改动改坏了"。所以改成反向断言：这些路由**必须 404**。
    _gone = ('/api/explorer/status', '/api/explorer/init',
             '/api/explorer/virus_options', '/api/explorer/primers',
             '/api/explorer/host', '/api/explorer/profile',
             '/api/explorer/profile_species')
    for rule in _gone:
        r = c.get(rule)
        check(r.status_code == 404, f'{rule} 已摘除（HTTP {r.status_code}，应为 404）')
    for rule in ('/api/explorer/query', '/api/explorer/charts',
                 '/api/explorer/vector', '/api/explorer/export/csv'):
        r = c.post(rule, json={})
        check(r.status_code == 404, f'{rule} 已摘除（HTTP {r.status_code}，应为 404）')
    # /reference/i18n.js 也属于 explorer 蓝图（实测 404）→ 一并纳入防回归
    r = c.get('/reference/i18n.js')
    check(r.status_code == 404, f'/reference/i18n.js 已摘除（HTTP {r.status_code}，应为 404）')

    # ---------- 5b. 进化动力学 t-phylodyn（2026-09-18 接线）+ 8 张归档卡防回归 ----------
    # 与 tools 页同一份模板，所以这里用 test_client 就能查（真浏览器在
    # _check_phylodyn_ui.py / _check_phylodyn_e2e.py）。
    r = c.get('/tools?g=phylodyn')
    html = r.data.decode('utf-8', 'replace')
    for mark in ('id="t-phylodyn"', 'id="pdFlow"', 'id="pd_input"',
                 'id="pdyn_meta"', 'id="toolrun-phylodyn"', 'id="pd_sub_on"',
                 'id="pd_coords_tsv"',
                 'function runPhylodyn(', 'function loadPhylodynResult(',
                 'const PD_STAGES', 'function drawTimeTree(',
                 # A2 迁移弧线地图（2026-09-18）：容器 + 坐标表输入 + 渲染函数。
                 # drawPgGeo 是**两卡共用**的实现，删了它系统地理卡的图也会没。
                 'id="pdGeoArc"', 'id="pd_coords"', 'function drawPgGeo(',
                 'function _pgArc('):
        check(mark in html, f'phylodyn 卡接线标记存在: {mark}')
    # ★ 撞车负控：新卡把撞车的 id 改成了 pdyn_*，而**保留的数据接入卡**仍用 pd_*
    check('id="pd_zip"' not in html and 'id="pdyn_zip"' not in html,
          'zip 来源的输入框已摘除（2026-09-18 explorer 不再导出 zip）')
    _rz = c.post('/api/phylodyn/import_export_zip', json={'zip_path': 'x.zip'})
    check(_rz.status_code == 404,
          f'/api/phylodyn/import_export_zip 已摘除（HTTP {_rz.status_code}，应为 404）')
    _gone = ('id="t-pdrand"', 'id="t-pdrrt"', 'id="t-pdtempmig"', 'id="t-pdbsp"',
             'id="t-pdrspp"', 'id="t-pdtreetime"', 'id="t-pdltt"', 'id="t-pdmjrm"',
             'function runPdBsp(', 'function runPdLtt(', 'function loadPdMjrmResult(',
             'id="t-pdrename"', 'id="t-pdgroup"', 'function runPdrename(')
    _left = [m for m in _gone if m in html]
    check(not _left, f'10 张归档卡的前端标记已消失（残留: {_left or "无"}）')
    from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY as _TR
    check('phylodyn' in _TR, 'TOOL_REGISTRY 注册了 phylodyn')
    for _t in ('pdrand', 'pdrrt', 'pdtempmig', 'pdbsp', 'pdrspp', 'pdtreetime',
               'pdltt', 'pdmjrm', 'pdrename', 'pdgroup'):
        check(_t not in _TR, f'归档工具已从注册表摘除: {_t}')
    for _t in ('pdprep',):
        check(_t in _TR, f'保留工具仍在注册表: {_t}')
    for _t in ('pdspacetime', 'pdsub'):
        check(_t not in _TR, f'已并入数据准备的卡已摘除: {_t}')
    _r2 = c.get('/tools?g=phylodyn')
    _h2 = _r2.data.decode('utf-8', 'replace')
    _left2 = [m for m in ('id="t-pdspacetime"', 'id="t-pdsub"', 'function runPdSub(',
                          'function runPdSpacetime(') if m in _h2]
    check(not _left2, f'两张已合并卡的前端标记已消失（残留: {_left2 or "无"}）')
    # 2026-09-18：摘卡留下的**非界面残影**（i18n 孤儿键 / LIGHT_TOOLS / 死函数）
    from Virus_Platform_Core.web.tools_api import LIGHT_TOOLS as _LT
    check('pdspacetime' not in _LT and 'pdsub' not in _LT,
          f'LIGHT_TOOLS 里没有已下线的工具（实际 {sorted(_LT)}）')
    _iz = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'webapp/static/i18n.js'),
        encoding='utf-8').read()
    _dead = [k for k in ('tk.pdrandBusy', 'tk.pdrrtStat', 'tk.pdbspX',
                         'tk.pdrsppStat', 'tk.pdlttY') if ("'%s'" % k) in _iz]
    check(not _dead, f'8 张归档卡的 i18n 孤儿键已清（残留: {_dead or "无"}）')
    _js = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'webapp/templates/tools.html'),
        encoding='utf-8').read()
    _df = [f for f in ('pdsubModeToggle', 'pdsubFromToggle')
           if ('function %s(' % f) in _js]
    check(not _df, f'已下线卡的死函数已清（残留: {_df or "无"}）')
    _gs = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'Virus_Platform_Core', 'geo_subsampler.py')
    check(not os.path.isfile(_gs),
          'geo_subsampler.py 已删（与 kit.subsample_fasta 同口径的重复实现，2026-09-18 收敛）')
    from Virus_Platform_Core import phylodyn_kit as _pk
    check(not hasattr(_pk, 'import_from_zip'),
          'import_from_zip 已删（zip 来源 2026-09-18 整条下线）')

    print('PLATFORM CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()
