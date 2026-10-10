# -*- coding: utf-8 -*-
"""工具箱 HTTP 边界：运行目录、任务提交、结果读取（自 app.py 拆出）。

Blueprint 名 tools，路由 /tools 与 /api/tool/*。
任务工厂本体在 Virus_Platform_Core/web/tool_jobs.py，链式在 Virus_Platform_Core/web/chains.py。
"""
import json
import os
import re
import time
import uuid

from flask import (abort, jsonify, render_template, request, send_file,
                   Response)

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT, db_path
from Virus_Platform_Core.utils import (check_path, fmt_size, safe_open)
from Virus_Platform_Core.web.state import (cfg, run_real_dir as _run_real_dir,
                                           tool_runs_root as _tool_runs_root,
                                           RUN_OUT_LINK)
from flask import Blueprint

from Virus_Platform_Core.web.tool_jobs import (
    _tool_job_convert,
    _tool_job_fastp,
    _tool_job_hostremoval,
    _tool_job_orf,
    _tool_job_orfa,
    _tool_job_genoplot,
    _tool_job_primer,
    _tool_job_identify,
    _tool_job_assemble,
    _tool_job_contigs,
    _tool_job_structcmp,
    _tool_job_verify,
    _tool_job_consensus,
    _tool_job_kvsuite,
    _tool_job_quicktree,
    _tool_job_align,
    _tool_job_sdt,
    _tool_job_identity,
    _tool_job_rdp,
    _tool_job_phylogeo,
    _tool_job_dsrna,
    _tool_job_crispr,
    _tool_job_rtt,
    _tool_job_pdprep,
    _tool_job_phylodyn,
    _tool_job_mirna,
    _tool_job_beastexport,
)
from Virus_Platform_Core.web.chains import (
    _tool_job_virchain,
    _tool_job_kvchain,
    _tool_job_prepchain,
    _tool_job_annochain,
)
from Virus_Platform_Core.web.tasks import tm

bp = Blueprint('tools', __name__)


TOOL_REGISTRY = {
    'convert':  {'title': '格式转换', 'title_en': 'Format convert',
                 'need_db': False, 'job': _tool_job_convert},
    'fastp':    {'title': '质控预处理', 'title_en': 'QC preprocess',
                 'need_db': False, 'job': _tool_job_fastp},
    'hostremoval': {'title': '宿主去除与序列提取', 'title_en': 'Host removal',
                    'need_db': False, 'job': _tool_job_hostremoval},
    # hostpredict（独立宿主预测工具）已于 2026-09-27 撤出注册表：孤儿页删除，
    # 功能保留为 contigs 分类后的自动运行（/api/tool/hostpredict_run，见
    # tool_results.py）与样品管道 ④（hostana）。
    'orf':         {'title': 'ORF 预测', 'title_en': 'ORF predict',
                    'need_db': False, 'job': _tool_job_orf},
    'orfa':        {'title': '功能注释', 'title_en': 'ORF annotate',
                    'need_db': False, 'job': _tool_job_orfa},
    'genoplot':    {'title': '基因组图谱', 'title_en': 'Genome plots',
                    'need_db': False, 'job': _tool_job_genoplot},
    'primer':      {'title': '引物设计', 'title_en': 'Primer design',
                    'need_db': False, 'job': _tool_job_primer},
    'identify': {'title': '病毒鉴定', 'title_en': 'Virus identify',
                 'need_db': True, 'job': _tool_job_identify},
    'assemble': {'title': '病毒组装', 'title_en': 'Virus assembly',
                 'need_db': False, 'job': _tool_job_assemble},
    'contigs':  {'title': 'contig分类', 'title_en': 'Contig classify',
                 'need_db': True, 'job': _tool_job_contigs},
    'sdt':      {'title': 'SDT 精确分析', 'title_en': 'SDT exact',
                 'need_db': False, 'job': _tool_job_sdt},
    'identity': {'title': 'NT+AA 同一性表', 'title_en': 'NT+AA identity',
                 'need_db': False, 'job': _tool_job_identity},
    'structcmp': {'title': '结构比较', 'title_en': 'Structure compare',
                  'need_db': False, 'job': _tool_job_structcmp},
    'verify':   {'title': '候选序列验证', 'title_en': 'Candidate verify',
                 'need_db': False, 'job': _tool_job_verify},
    'consensus': {'title': '共识序列与变异', 'title_en': 'Consensus & variants',
                  'need_db': False, 'job': _tool_job_consensus},
    'kvsuite':  {'title': '已知病毒识别与定量', 'title_en': 'Known virus suite',
                 'need_db': False, 'job': _tool_job_kvsuite},
    'align':    {'title': '序列比对', 'title_en': 'Alignment',
                 'need_db': False, 'job': _tool_job_align},
    'quicktree': {'title': '快速建树', 'title_en': 'Quick tree',
                  'need_db': False, 'job': _tool_job_quicktree},
    'virchain': {'title': '病毒识别分类·一键', 'title_en': 'Virus chain',
                 'need_db': True, 'job': _tool_job_virchain},
    'kvchain':  {'title': '病毒定量与共识·一键', 'title_en': 'KV chain',
                 'need_db': False, 'job': _tool_job_kvchain},
    'prepchain': {'title': '一键前处理', 'title_en': 'Prep chain',
                  'need_db': False, 'job': _tool_job_prepchain},
    'annochain': {'title': '注释三连·一键', 'title_en': 'Annotation chain',
                  'need_db': False, 'job': _tool_job_annochain},
    'rdp':      {'title': 'RDP 重组分析', 'title_en': 'RDP recombination',
                 'need_db': False, 'job': _tool_job_rdp},
    'phylogeo': {'title': '系统地理分析', 'title_en': 'Phylogeography',
                 'need_db': False, 'job': _tool_job_phylogeo},
    'dsrna':    {'title': 'dsRNA 设计', 'title_en': 'dsRNA design',
                 'need_db': False, 'job': _tool_job_dsrna},
    'crispr':   {'title': 'CRISPR 检测试剂设计',
                 'title_en': 'CRISPR assay design',
                 'need_db': False, 'job': _tool_job_crispr},
    'rtt':      {'title': '时间信号检验', 'title_en': 'Root-to-tip',
                 'need_db': False, 'job': _tool_job_rtt},
    # ---- 进化动力学（保留 5 张数据准备类卡）----
    'phylodyn': {'title': '⚡ 一键分析（定年 → 地理迁移 → 汇总溯源）',
                 'title_en': 'One-click analysis: dating + geographic migration',
                 'need_db': False, 'job': _tool_job_phylodyn},
    # 2026-09-18 归档 8 张卡（rand/rrt/tempmig/bsp/rspp/treetime/ltt/mjrm）：
    # 进化平台已用「时间与地理推断·本地全链」统一覆盖，保留 5 张数据准备类卡。
    # 底层 phylodyn_kit / phylodyn_trees **未删**（保留卡仍在用）；
    # 前端 section / 导航 / job 工厂已摘，详见 docs/进化动力学统一方案_20260918.md
    'pdprep':   {'title': '进化动力学数据接入', 'title_en': 'Phylodyn data prep',
                 'need_db': False, 'job': _tool_job_pdprep},
    # 2026-09-21 归档 pdrename / pdgroup 两卡（job 注册已摘）：重命名已接进
    # phylodyn 全链可选前置（t-phylodyn「🏷 重命名映射」字段，A0 三处同改）；
    # 分组属导入/收集阶段职责（online 自动跑，产物进 report）。job 函数保留于
    # tool_jobs.py，底层 phylodyn_kit 未删。
    'mirna':    {'title': 'miRNA 靶标预测', 'title_en': 'miRNA target',
                 'need_db': False, 'job': _tool_job_mirna},
    # BEAST 作业包导出（轨道 B 上半程，2026-10-08 接线）：本地只打包、不跑 MCMC
    # （引擎 beast_handoff.py，此前无任何调用方）。产物喂服务器，回传的 .log/.trees
    # 由 phylogeo 卡的「BEAST 产物导入」读入。
    'beastexport': {'title': 'BEAST 作业包导出',
                    'title_en': 'BEAST job export',
                    'need_db': False, 'job': _tool_job_beastexport},
}


# ⚠️ 这里的键必须**都是 TOOL_REGISTRY 里真实存在的工具**：
#    `pdspacetime` 已于 2026-09-18 并入数据准备卡 → 已删；
#    `pdrename` / `pdgroup` 已于 2026-09-21 归档 → 已删
LIGHT_TOOLS = {'convert', 'genoplot', 'primer', 'mirna', 'beastexport'}


# 工具②/④的病毒库内置键名。RVDB（databases/virus/rvdb）已于本会话移除，
# 保留键名只会让 /api/dbs 的“可用内置库”提示与“库未就绪”报错指向不存在的目录。
VIRUS_DB_PRESETS = {'virus': ('virus', 'plant'),
                    'refvirus': ('virus', 'ref'),
                    'k2viral': None}


@bp.route('/tools')
def page_tools():
    return render_template('tools.html')


@bp.route('/tool_runs/<run>/<path:filename>')
def page_tool_run_file(run, filename):
    """工具运行目录内文件的下载/预览（严格限制在该运行目录内）。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run) or '..' in filename:
        abort(400, '无效的运行名')
    # 外移运行：产物在自定义输出目录（out_link.json 登记），从这里寻址。
    # base 归一 + 前缀包含校验：filename 传绝对路径时 os.path.join 会整段
    # 替换 base，check_path(in_platform=True) 拦不住读平台内其它文件
    # （如 run/ai_local.json 里的 API Key）。
    base = os.path.abspath(_run_real_dir(run))
    p = check_path(os.path.join(base, filename),
                   must_exist=True, in_platform=True)
    if not p.startswith(base + os.sep):
        abort(400, '非法路径')
    resp = send_file(p)
    # .svg 必须在浏览器内联渲染：send_file 在此环境把 .svg 推断成了
    # 非标准的 image/svg（浏览器会拒绝渲染 → broken image），改回标准 MIME。
    if str(filename).lower().endswith('.svg'):
        resp.mimetype = 'image/svg+xml'
    if request.args.get('dl'):          # 显式 ?dl=1 才强制下载；
        import urllib.parse as _up      # 不带时内联渲染（report iframe 需要）
        resp.headers['Content-Disposition'] = (
            "attachment; filename*=UTF-8''"
            + _up.quote(os.path.basename(filename)))
    return resp


# 最近运行列表缓存：{signature: (ts, payload)}
# 本接口对最近 30 个运行目录做整棵 os.walk + 逐个 getsize，装了大组装/MSA
# 的运行目录动辄数千文件；而工具页与「结果恢复」都会反复拉它。签名取
# (运行名, 目录 mtime) 列表 —— 新增/删除运行或顶层变动立刻失效，嵌套变动
# 由 _RUNS_TTL 兜底。
_RUNS_CACHE = {}
_RUNS_TTL = 4.0


def _runs_signature(root):
    try:
        names = [n for n in os.listdir(root)
                 if not n.startswith('_')
                 and os.path.isdir(os.path.join(root, n))]
    except OSError:
        return None, []
    sig = []
    for n in names:
        try:
            sig.append((n, os.path.getmtime(os.path.join(root, n))))
        except OSError:
            sig.append((n, 0))
    sig.sort(key=lambda x: x[1], reverse=True)
    return tuple(sig), [n for n, _ in sig]


@bp.route('/api/tool/runs')
def api_tool_runs():
    """最近的工具运行列表（按运行目录修改时间倒序，30 条）。

    按目录名排序会让字母序靠前的前缀（contigs/genoplot…）挤掉真正
    最新的运行，这里以 mtime 为准。

    结果按「(运行名, 目录 mtime) 签名 + 短 TTL」缓存：整棵 os.walk 是本
    接口的主要开销，而它会被工具页反复轮询。
    """
    import time as _time
    root = _tool_runs_root()
    sig, names = _runs_signature(root)
    if sig is None:
        return jsonify([])
    hit = _RUNS_CACHE.get(sig[0][0] if sig else '')
    now = _time.time()
    if hit and hit[0] == sig and now - hit[1] < _RUNS_TTL:
        return jsonify(hit[2])
    out = []
    for name in names[:30]:
        d = _run_real_dir(name)
        if not os.path.isdir(d):
            continue
        files = []
        for cur, _sub, fns in os.walk(d):
            rel = os.path.relpath(cur, d)
            for fn in fns:
                p = os.path.join(cur, fn)
                try:
                    size = os.path.getsize(p)
                except OSError:
                    # MMseqs2 临时库等系统级特殊文件 stat 不到（WinError
                    # 1920），跳过大小即可，不让整个列表接口 500
                    continue
                files.append({
                    'path': fn if rel == '.' else os.path.join(rel, fn),
                    'size': fmt_size(size)})
        out.append({'name': name, 'out_dir': d, 'files': files[:40]})
    _RUNS_CACHE.clear()                     # 只保留最新签名一份，防内存增长
    _RUNS_CACHE[sig[0][0] if sig else ''] = (sig, now, out)
    return jsonify(out)


@bp.route('/api/tool/runs/<run>/delete', methods=['POST'])
def api_tool_run_delete(run):
    """删除一个工具运行目录（模块历史区的 🗑；运行名校验 + 平台内限定）。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    real = _run_real_dir(run)
    p = check_path(real, must_exist=True, in_platform=True)
    # 有任务正在这个目录里跑（任务卡 out_dir == 本运行目录）时拒绝删除：
    # 否则删掉的是运行中的结果目录，任务随后写盘失败、产物残缺且
    # 任务卡还挂在「运行中」。忙时 409（与建库接口同一口径）。
    _p_key = os.path.normcase(os.path.abspath(p))
    with tm.lock:
        _busy = [r for r in tm.tasks.values()
                 if r.get('status') in ('running', 'queued')
                 and r.get('out_dir')
                 and os.path.normcase(os.path.abspath(str(r['out_dir'])))
                 == _p_key]
    if _busy:
        abort(409, f'该运行目录上有任务正在跑（{_busy[0].get("id")} '
                   f'{_busy[0].get("name")}）；请先停止任务再删除')

    def _rmtree_retry(d):
        """Windows 边界：刚被 send_file 下载过的产物句柄可能尚未完全
        释放，首趟 rmtree(ignore_errors) 会静默留下它；回收后重试一趟。"""
        import shutil as _sh
        _sh.rmtree(d, ignore_errors=True)
        if os.path.isdir(d):
            import gc as _gc
            _gc.collect()
            time.sleep(0.3)
            _sh.rmtree(d, ignore_errors=True)

    _rmtree_retry(p)
    # 指针目录与真实目录不同（外移运行）时都清：历史区不留死登记
    base = check_path(os.path.join(_tool_runs_root(), run),
                      must_exist=False, in_platform=True)
    if os.path.abspath(base) != os.path.abspath(p):
        _rmtree_retry(base)
    _RUNS_CACHE.clear()                 # 删了运行，最近运行列表立刻失效
    return jsonify({'ok': True, 'run': run})


@bp.route('/api/tool/open', methods=['POST'])
def api_tool_open():
    """在资源管理器中打开某次工具运行的目录。"""
    body = request.get_json(force=True) or {}
    name = body.get('name') or ''
    if not name or '/' in name or '\\' in name or '..' in name:
        abort(400, '无效的运行名')
    d = check_path(_run_real_dir(name), must_exist=True, in_platform=True)
    from Virus_Platform_Core.utils import open_in_explorer
    open_in_explorer(d)
    return jsonify({'ok': True})


@bp.route('/api/tool/genoplot_preview', methods=['POST'])
def api_tool_genoplot_preview():
    """基因组图谱交互预览：接收输入 + 绘图参数，出图到临时目录，返回 SVG 内容。

    不写 tool_runs/、不建运行目录（交互模式下仅预览不落盘）。
    参数：fasta / ann（.gb/.gff）+ gb_opts（dict，透传 gbdraw CLI）。
    返回：{ok, svg, mode}（svg 为 circular SVG 文本；失败抛 400/500）。
    """
    import tempfile
    body = request.get_json(force=True) or {}
    fasta = (body.get('fasta') or '').strip()
    ann = (body.get('ann') or '').strip()
    gb_opts = body.get('gb_opts') or {}
    if not fasta and not ann:
        abort(400, '请提供 FASTA 或 GenBank 输入')
    def _pabs(v):
        v = v.strip() if v else ''
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)
    fasta_abs = _pabs(fasta) if fasta else None
    ann_abs = _pabs(ann) if ann else None
    if ann_abs and str(ann_abs).lower().endswith(
            ('.gb', '.gbk', '.gbff', '.genbank')):
        fasta_abs = None
    # 过滤 opts（去空/False/None）。预览只出 SVG，排除仅保存时生效的格式/多记录。
    _PREVIEW_SKIP = {'format', 'multi_record_canvas'}
    opts = {str(k): v for k, v in gb_opts.items()
            if v not in (None, '', False) and k not in _PREVIEW_SKIP}
    tmp = tempfile.mkdtemp(prefix='vp_genoplot_preview_')
    try:
        # 用本次请求专属的 ascii 临时目录（gbdraw 需 ascii 路径）。
        # 原先用固定路径 %TEMP%\vp_geno_preview\preview.svg 且从不清理：
        # 两个并发预览会互相覆盖，先返回的请求可能拿到后者的图。
        import shutil as _sh
        tmp_ascii = os.path.join(tmp, 'geno')
        os.makedirs(tmp_ascii, exist_ok=True)
        from Virus_Platform_Core.gbdraw_plot import _run_gbdraw
        fasta_p, gff_p, gbk_p = None, None, None
        if ann_abs and str(ann_abs).lower().endswith(
                ('.gb', '.gbk', '.gbff', '.genbank')):
            gbk_p = ann_abs
        else:
            fasta_p = fasta_abs
            # 可选 GFF：若 ann 是 .gff 且 fasta 也给了 → 配对
            if ann_abs and str(ann_abs).lower().endswith(('.gff', '.gff3')):
                gff_p = ann_abs
        made = _run_gbdraw(fasta=fasta_p, gff=gff_p, gbk=gbk_p,
                           out_prefix=os.path.join(tmp_ascii, 'preview'),
                           mode=body.get('mode') or 'circular', opts=opts)
        if not made:
            abort(400, 'gbdraw 未产出预览图')
        svg_path = made[0]
        with safe_open(svg_path) as f:
            svg = f.read()
        return jsonify({'ok': True, 'mode': 'circular', 'svg': svg})
    finally:
        _sh.rmtree(tmp, ignore_errors=True)


@bp.route('/api/tool/run', methods=['POST'])
def api_tool_run():
    """运行独立工具。body: {tool, threads, params{...}}；返回 {task, run}。"""
    from types import SimpleNamespace
    body = request.get_json(force=True) or {}
    tool = body.get('tool')
    threads = body.get('threads') or None
    p = body.get('params') or {}

    spec = TOOL_REGISTRY.get(tool)
    if not spec:
        abort(400, '未知工具')

    def _input_abs(v):
        """输入数据文件：相对路径按平台根解析；绝对路径允许电脑任意位置
        （只读）。输出产物仍限制在平台目录内。"""
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)

    def _req(key, what):
        v = (p.get(key) or '').strip()
        if not v:
            abort(400, f'请选择{what}')
        return _input_abs(v)

    def _opt(key):
        v = (p.get(key) or '').strip()
        return _input_abs(v) if v else None

    db_virus = _resolve_virus_db(p.get('db_virus'))
    if spec['need_db']:
        from Virus_Platform_Core.kunpeng import db_ready
        if not db_ready(db_virus):
            abort(400, f'所选病毒库未就绪: {db_virus}'
                       '（先到「数据库构建」页构建，或换选其他病毒库）')

    ts = time.strftime('%Y%m%d_%H%M%S')
    # 秒级时间戳会撞车：同一秒两次同工具运行（双击/多标签页）会复用同一
    # 目录并互相覆盖产物。追加 6 位随机后缀。
    run_name = f'{tool}_{ts}_{uuid.uuid4().hex[:6]}'
    # 自定义输出位置（可选）：指定目录即产物父目录，产物写 <out_dir>/<run_name>/
    # （结构与默认 tool_runs/<run_name> 一致，只换父目录，多次运行不互混）。
    # 校验同下载页：平台内或设置页已挂载输出根内，越界 400。
    # tool_runs/<run_name>/ 仍创建并写 out_link.json 指针——历史列表、产物
    # 下载、打开目录都从 tool_runs 寻址，没有登记就全部失联。
    base_dir = check_path(os.path.join(_tool_runs_root(), run_name),
                          must_exist=False, in_platform=True)
    out_root = (body.get('out_dir') or '').strip()
    if out_root:
        out_root = check_path(out_root, must_exist=False, in_platform=True)
        # 结构性目录黑名单：这些目录各有扫描器按「样品 / 库 / 任务 / 静态
        # 资源」解释其子目录。允许把工具运行写进去会造出"幽灵样品"
        # （results/<run> 出现在样品列表里、点开全 404、还删不掉），
        # databases/ 下的运行目录会被当成病毒库。前缀命中即拒绝。
        _out_key = os.path.normcase(os.path.normpath(os.path.abspath(out_root)))
        _bad_roots = (
            ('样品结果目录', DIRS['results']),
            ('数据库目录', DIRS['databases']),
            ('任务状态目录', DIRS['tasks']),
            ('前端资源目录', os.path.join(PLATFORM_ROOT, 'webapp')),
            ('第三方依赖目录', os.path.join(PLATFORM_ROOT, '3rd')),
        )
        for _label, _bad in _bad_roots:
            _bad_key = os.path.normcase(
                os.path.normpath(os.path.abspath(str(_bad))))
            if _out_key == _bad_key or _out_key.startswith(_bad_key + os.sep):
                abort(400, f'输出目录不能位于{_label}内（{_bad}）；请换一个'
                           '目录（如平台 run/ 下的自定义输出根）')
        run_dir = check_path(os.path.join(out_root, run_name),
                             must_exist=False, in_platform=True)
    else:
        run_dir = base_dir
    os.makedirs(run_dir, exist_ok=True)
    if run_dir != base_dir:
        os.makedirs(base_dir, exist_ok=True)
        try:
            with safe_open(os.path.join(base_dir, RUN_OUT_LINK),
                           'wt') as f:
                json.dump({'out_dir': run_dir, 'tool': tool}, f,
                          ensure_ascii=False, indent=2)
        except OSError:
            pass    # 指针写失败只影响历史区寻址，不阻塞运行
    # 数字参数钳制：负数/异常大值不直接传给外部工具
    if threads is not None:
        threads = max(1, min(int(threads), (os.cpu_count() or 4) * 4))

    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    try:
        job = spec['job'](ctx)
    except Exception:
        # 工厂期 abort（缺必填参数）不该留下空运行目录
        import shutil as _sh
        _sh.rmtree(run_dir, ignore_errors=True)
        if run_dir != base_dir:
            _sh.rmtree(base_dir, ignore_errors=True)
        raise

    name = cfg.tr(f"工具·{spec['title']} {ts}", f"Tool·{spec['title_en']} {ts}")

    def _job_with_run(log, prog, cancel):
        """包一层：结果 dict 注入 run 名与输出目录，供结果面板列出产物、
        明示「结果输出在哪」。"""
        res = job(log, prog, cancel)
        if isinstance(res, dict):
            res.setdefault('run', run_name)
            res['out_dir'] = str(run_dir)
        # 外移运行结束时 touch 指针目录：历史列表签名基于 tool_runs/<run>
        # 的 mtime，产物写在外部目录时没有这一下，列表里的文件清单
        # 会停在创建时刻（只剩 out_link.json）。
        try:
            os.utime(base_dir)
        except OSError:
            pass
        return res

    tid = tm.start(name, _job_with_run,
                   log_file=os.path.join(run_dir, 'run.log'),
                   weight=('light' if tool in LIGHT_TOOLS else 'heavy'),
                   out_dir=str(run_dir))
    return jsonify({'task': tid, 'run': run_name, 'out_dir': str(run_dir)})


@bp.route('/api/tool/chain_result')
def api_tool_chain_result():
    """一键流程汇总目录：返回 _chain.json 内容 + 链接清单。

    chain = virchain / kvchain；run 为汇总目录名（留空取最新一次）。
    """
    chain = (request.args.get('chain') or '').strip()
    if chain not in ('virchain', 'kvchain'):
        abort(400, '无效的链式运行类型')
    run = (request.args.get('run') or '').strip()
    root = _tool_runs_root()
    if run:
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        if not run.startswith(chain + '_'):
            abort(400, '运行名与链式类型不匹配')
    else:
        cands = []
        if os.path.isdir(root):
            for d in os.listdir(root):
                if d.startswith(chain + '_') and os.path.isfile(
                        os.path.join(root, d, '_chain.json')):
                    cands.append(d)
        if not cands:
            return jsonify({'run': None, 'out_dir': None, 'meta': None,
                            'entries': []})
        run = sorted(cands)[-1]
    cdir = check_path(os.path.join(root, run), must_exist=False,
                      in_platform=True)
    meta = None
    cj = os.path.join(cdir, '_chain.json')
    if os.path.isfile(cj):
        try:
            with safe_open(cj) as f:
                meta = json.load(f)
        except (OSError, ValueError):
            meta = None
    entries = []
    if os.path.isdir(cdir):
        for name in sorted(os.listdir(cdir)):
            if name == '_chain.json':
                continue
            full = os.path.join(cdir, name)
            entries.append({'name': name,
                            'is_dir': os.path.isdir(full),
                            'is_link': os.path.islink(full),
                            'size': (None if os.path.isdir(full)
                                     else fmt_size(os.path.getsize(full))
                                     if os.path.exists(full) else None)})
    return jsonify({'run': run, 'out_dir': cdir, 'meta': meta,
                    'entries': entries})


def _resolve_virus_db(v=None):
    """工具②/④的病毒库选择：内置键名、平台相对路径或本机绝对路径。

    缺省用平台主病毒库（cfg.databases['virus']）。返回库目录绝对路径；
    键名不识别且路径不存在时 abort(400)。
    """
    if not v:
        return cfg.databases['virus']
    v = str(v).strip()
    if v in VIRUS_DB_PRESETS:
        mapping = VIRUS_DB_PRESETS[v]
        if mapping:
            return db_path(*mapping)
        return os.path.join(DIRS['databases'], 'k2viral_db')
    if os.path.isabs(v):
        d = v
    else:
        d = os.path.join(PLATFORM_ROOT, v)
        # 兼容旧前端：'databases/<rel>' 的语义是**相对数据库根**，不是相对
        # 平台根。/api/dbs 早期就是这么输出的，而「数据库包与程序目录分离」
        # 的分发版下两者不同（DIRS['databases'] = <数据库包>\databases），
        # 按平台根解析会指向不存在的 <程序目录>\databases\...。
        rel = v.replace('\\', '/')
        if rel.startswith('databases/'):
            alt = os.path.join(DIRS['databases'], rel[len('databases/'):])
            if os.path.exists(alt):
                d = alt
    try:
        return check_path(d, must_exist=True)
    except (ValueError, FileNotFoundError):
        abort(400, f'病毒库不存在: {v}（可用内置库: '
                   f'{", ".join(VIRUS_DB_PRESETS)} 或填库目录路径）')


@bp.route('/api/tool/rdp_masked')
def api_rdp_masked():
    """导出掩蔽重组区后的干净比对（N 掩蔽，可直接送建树）。

    run=<rdp run 名>；读 rdp/events.json + summary.json 里的 input 路径。
    """
    import json as _json
    run = (request.args.get('run') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run or ''):
        abort(400, '无效的运行名')
    base = os.path.join(_run_real_dir(run), 'rdp')
    ev_file = os.path.join(base, 'events.json')
    if not os.path.isfile(ev_file):
        abort(404, '该运行无重组事件结果')
    events = _json.load(open(ev_file, encoding='utf-8')).get('events') or []
    sm_file = os.path.join(base, 'summary.json')
    inp = ''
    if os.path.isfile(sm_file):
        inp = (_json.load(open(sm_file, encoding='utf-8'))
               or {}).get('input') or ''
    if not inp or not os.path.isfile(inp):
        abort(400, '找不到原比对文件（summary.json 缺 input 记录）')
    from Virus_Platform_Core.rdp5_engine import mask_from_events
    records, n_masked, merged = mask_from_events(inp, events, mask_char='N')
    if not records:
        abort(400, '掩蔽失败：比对为空')
    nl = chr(10)
    body = nl.join(f'>{nm}{nl}{seq}' for nm, seq in records) + nl
    note = (f'; masked {n_masked} columns in {len(merged)} regions' + nl)
    out_body = body + (nl + note if note else '')
    return Response(out_body, mimetype='text/plain',
                    headers={'Content-Disposition':
                             f'attachment; filename="{run}_masked.fa"'})


@bp.route('/api/tool/rdp_masked', methods=['POST'])
def api_rdp_masked_post():
    return api_rdp_masked()


@bp.route('/api/tool/dsrna_panel')
def api_dsrna_panel():
    """dsRNA 脱靶面板清单（物种名 + 大小 + 是否配致死基因表）。

    只读 databases/dsrna/ 下的约定目录；目录缺失时返回空清单（前端提示
    自行放置面板 FASTA），不当作错误——面板是可选数据。
    """
    from Virus_Platform_Core import dsrna_offtarget as dso
    from Virus_Platform_Core.dsrna_pipeline import _norm_species
    # ⚠️ 必须走 DIRS['databases']：分发版由 config._detect_database_root() 把
    # 数据库根重定向到同级 VirusPlatform-Database，与 PLATFORM_ROOT 不等价。
    # 用 PLATFORM_ROOT 时面板清单会恒为空（本函数把"目录缺失"当正常），
    # 用户以为平台没有面板数据 —— 静默失效。2026-10-07 审查 P1-2。
    panel_dir = os.path.join(DIRS['databases'], 'dsrna', 'panel')
    lethal_dir = os.path.join(DIRS['databases'], 'dsrna',
                              'lethal_lists')
    files = []
    for f in dso.panel_files(panel_dir):
        sp = _norm_species(os.path.splitext(os.path.basename(f))[0])
        try:
            size = os.path.getsize(f)
        except OSError:
            size = 0
        files.append({'name': sp,
                      'file': os.path.basename(f),
                      'size_mb': round(size / 1e6, 1),
                      'lethal': os.path.isfile(
                          os.path.join(lethal_dir,
                                       f'{sp}_all_lethals.tsv'))})
    return jsonify({'panel_dir': panel_dir, 'files': files})


# ---- 迁移弧线动画导出 GIF（plotly + kaleido + Pillow，全离线）----
# 上限值是**请求参数**的钳制范围（防前端传超大值把渲染拖死），不是数据上限。
_GIF_MAX_FRAMES = 240
_GIF_MAX_ARCS = 60
_GIF_MAX_PX = 2400


def _phylogeo_run_dir(run):
    """校验 run 名并定位 <run>/phylogeo 结果目录。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run or ''):
        abort(400, '无效的运行名')
    d = os.path.join(_run_real_dir(run), 'phylogeo')
    if not os.path.isdir(d):
        abort(404, '该运行没有系统地理分析结果（phylogeo 目录不存在）')
    return check_path(d, must_exist=True, in_platform=True)


def _gif_int(v, lo, hi, default):
    try:
        v = int(v)
    except (TypeError, ValueError):
        return default
    return max(lo, min(v, hi))


@bp.route('/api/tool/phylogeo_gif', methods=['POST'])
def api_phylogeo_gif():
    """把 MOTP 时间窗动画导出成 GIF（plotly+kaleido+Pillow，全离线）。

    body: {run, force=0, fps=2, max_frames=80, max_arcs=24, width=960,
           height=520}
    产物固定落 <run>/phylogeo/migration.gif；summary.json 比它新（说明重跑过
    分析）或 force=1 时重渲染，否则直接返回缓存 —— 缓存时不占用任务名额。
    """
    body = request.get_json(force=True) or {}
    run = (body.get('run') or '').strip()
    gdir = _phylogeo_run_dir(run)
    sm = os.path.join(gdir, 'summary.json')
    if not os.path.isfile(sm):
        abort(404, '该运行没有 phylogeo/summary.json（结果不完整）')
    out = os.path.join(gdir, 'migration.gif')
    url = '/tool_runs/%s/phylogeo/migration.gif' % run
    force = str(body.get('force') or '').lower() not in ('', '0', 'false', 'no')
    if (not force and os.path.isfile(out)
            and os.path.getmtime(out) >= os.path.getmtime(sm)):
        st = os.stat(out)
        return jsonify({'cached': True, 'task': None, 'url': url,
                        'file': 'migration.gif', 'bytes': st.st_size,
                        'mtime': st.st_mtime})

    fps = _gif_int(body.get('fps'), 1, 10, 2)
    max_frames = _gif_int(body.get('max_frames'), 1, _GIF_MAX_FRAMES, 80)
    max_arcs = _gif_int(body.get('max_arcs'), 1, _GIF_MAX_ARCS, 24)
    width = _gif_int(body.get('width'), 320, _GIF_MAX_PX, 960)
    height = _gif_int(body.get('height'), 240, _GIF_MAX_PX, 520)

    def job(log, prog, cancel):
        from Virus_Platform_Core.phylogeo_gif import build_migration_gif
        with open(sm, encoding='utf-8') as f:
            summary = json.load(f)
        log('读取 summary.json，开始逐帧渲染（kaleido 本地出图）')
        prog('gif', 2, '准备渲染')

        def _cb(frac, msg):
            log(msg)
            prog('gif', 2 + 96 * frac, msg)

        res = build_migration_gif(summary, out, width=width, height=height,
                                  fps=fps, max_arcs=max_arcs,
                                  max_frames=max_frames, progress=_cb,
                                  cancel=cancel)
        log('GIF 落盘：%s（%d 帧 / %.0f KB / %dx%d）'
            % (res['file'], res['n_frames'], res['bytes'] / 1024.0,
               res['width'], res['height']))
        if res['frames_dropped']:
            log('窗数 %d 超过帧数上限，每 %d 窗取 1 帧（丢 %d 帧）'
                % (res['n_bins'], res['every_n_bins'], res['frames_dropped']))
        if res['n_arcs_skipped']:
            log('有 %d 条迁移两端缺坐标，未画进 GIF' % res['n_arcs_skipped'])
        res['url'] = url
        res['run'] = run
        res['out_dir'] = gdir
        return res

    tid = tm.start('导出迁移 GIF·%s' % run, job, weight='light', out_dir=gdir)
    return jsonify({'cached': False, 'task': tid, 'url': url, 'out_dir': gdir,
                    'params': {'fps': fps, 'max_frames': max_frames,
                               'max_arcs': max_arcs, 'width': width,
                               'height': height}})


@bp.route('/api/tool/pdplot_export', methods=['POST'])
def api_pdplot_export():
    """交互图一键导出 PNG/PDF（plotly → kaleido，全离线）。

    body: {name, fmt: 'png'|'pdf', data, layout, width?, height?}
    返回图片字节流（attachment 下载）。
    """
    body = request.get_json(force=True) or {}
    fmt = body.get('fmt') or 'png'
    if fmt not in ('png', 'pdf', 'svg', 'webp'):
        abort(400, f'不支持的导出格式: {fmt}')
    data = body.get('data')
    layout = body.get('layout') or {}
    if not data:
        abort(400, '缺少绘图数据')
    width = max(320, min(int(body.get('width') or 1280), 4096))
    height = max(240, min(int(body.get('height') or 800), 4096))
    import plotly.io as pio
    fig = {'data': data, 'layout': layout}
    try:
        img = pio.to_image(fig, format=fmt, width=width, height=height, scale=2)
    except Exception as e:
        abort(500, f'渲染失败（kaleido）: {e}')
    name = re.sub(r'[^\w.-]+', '_', body.get('name') or 'plot')[:60] or 'plot'
    from flask import Response
    mime = {'png': 'image/png', 'pdf': 'application/pdf',
            'svg': 'image/svg+xml', 'webp': 'image/webp'}[fmt]
    return Response(img, mimetype=mime, headers={
        'Content-Disposition': f'attachment; filename="{name}.{fmt}"'})


# ── t-mirna 候选 miRNA 库（构建与验证口径见 databases/mirna_lib/README.md）──
# 数据：服务器 /genes/gene_browser.html 下载区（plantrg + 同源种扩展，80,272 条 /
# 799 物种），清洗后 TSV；首次请求懒加载，文件 mtime 变化自动重载。
_MIRNA_LIB_CACHE = {'rows': None, 'mtime': None, 'facets': None}


def _mirna_lib_path():
    # DIRS['databases'] 而非 PLATFORM_ROOT/databases（分发版数据库外置会重定向）
    # —— 2026-10-07 审查 P1-2
    return os.path.join(DIRS['databases'], 'mirna_lib',
                        'plant_mirna_library.tsv')


def _mirna_lib_rows():
    path = _mirna_lib_path()
    if not os.path.exists(path):
        abort(500, '候选 miRNA 库缺失，请先运行 databases/mirna_lib/build_library.py --build')
    mt = os.path.getmtime(path)
    if _MIRNA_LIB_CACHE['rows'] is None or _MIRNA_LIB_CACHE['mtime'] != mt:
        rows = []
        with open(path, 'r', encoding='utf-8') as f:
            header = f.readline().rstrip('\n').split('\t')
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) != len(header):
                    continue
                r = dict(zip(header, parts))
                try:
                    r['len'] = int(r.get('len') or 0)
                except ValueError:
                    continue
                tc = r.get('target_count') or ''
                r['_tc'] = int(tc) if tc.isdigit() else -1
                rows.append(r)
        facets_s = {}
        facets_src = {}
        for r in rows:
            facets_s[r['species']] = facets_s.get(r['species'], 0) + 1
            facets_src[r['source']] = facets_src.get(r['source'], 0) + 1
        _MIRNA_LIB_CACHE['rows'] = rows
        _MIRNA_LIB_CACHE['mtime'] = mt
        _MIRNA_LIB_CACHE['facets'] = {
            'total': len(rows),
            'species': sorted(facets_s.items(), key=lambda kv: (-kv[1], kv[0])),
            'sources': sorted(facets_src.items(), key=lambda kv: -kv[1]),
        }
    return _MIRNA_LIB_CACHE['rows']


@bp.route('/api/tool/mirna_lib/facets')
def api_mirna_lib_facets():
    _mirna_lib_rows()
    return jsonify(_MIRNA_LIB_CACHE['facets'])


@bp.route('/api/tool/mirna_lib')
def api_mirna_lib():
    """检索候选 miRNA：q 子串匹配 id/物种/family；species/source 精确过滤。
    排序：plantrg 规范命名优先 → 有 target_count 注释者按注释降序 → id。"""
    rows = _mirna_lib_rows()
    q = (request.args.get('q') or '').strip().lower()
    species = (request.args.get('species') or '').strip()
    source = (request.args.get('source') or '').strip()
    try:
        limit = max(1, min(int(request.args.get('limit') or 100), 200))
    except ValueError:
        limit = 100
    out = []
    total = 0
    for r in rows:
        if species and r['species'] != species:
            continue
        if source and r['source'] != source:
            continue
        if q and (r['mirna_id'].lower().find(q) < 0
                  and r['species'].lower().find(q) < 0
                  and (r.get('family') or '').lower().find(q) < 0):
            continue
        total += 1
        if len(out) < limit:
            out.append({'mirna_id': r['mirna_id'], 'species': r['species'],
                        'seq': r['seq'], 'len': r['len'], 'source': r['source'],
                        'family': r.get('family') or '',
                        'target_count': r.get('target_count') or ''})
    out.sort(key=lambda r: (0 if r['source'] == 'plantrg' else 1,
                            -int(r['target_count']) if r['target_count'].isdigit() else 1,
                            r['mirna_id']))
    return jsonify({'total': total, 'items': out})
