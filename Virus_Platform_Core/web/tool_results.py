# -*- coding: utf-8 -*-
"""工具箱结果读取 API（自 app.py 拆出）。

kvsuite / consensus / verify / 结构比较 / MSA / 建树 / 比对编辑 /
报告数据 / 宿主预测 / contig 深度分析。"""
import json
import logging
import os
import re
import threading
import time

from flask import (Blueprint, Response, abort, jsonify, request,
                   send_file)

from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.utils import (TaskLogger, check_path, fmt_size, safe_open)
from Virus_Platform_Core.web.common import tree_overlay_for
from Virus_Platform_Core.web.state import (cfg, run_real_dir as _run_real_dir,
                                           tool_runs_root as _tool_runs_root)
from Virus_Platform_Core.web.tasks import tm

_log = logging.getLogger('Virus_Platform_Core.web.tool_results')
from Virus_Platform_Core.web.tool_jobs import _find_contig_seq

bp = Blueprint('tool_results', __name__)


def _analysis_out(run_dir, contig, action):
    """分析结果 JSON 的固定路径：run/analysis/<安全化contig>_<action>.json。

    run/contig/action 先做严格字符校验 + 规范化，杜绝 ../ 穿越。
    """
    if action not in ('blastn', 'blastx', 'primer'):
        abort(400, '无效的分析类型')
    if not re.fullmatch(r'[A-Za-z0-9_\-]{1,80}', os.path.basename(run_dir)):
        abort(400, '无效的运行名')
    adir = check_path(os.path.join(run_dir, 'analysis'), must_exist=False,
                      in_platform=True)
    os.makedirs(adir, exist_ok=True)
    safe_c = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:60]
    p = os.path.abspath(os.path.join(adir, f'{safe_c}_{action}.json'))
    if not p.startswith(os.path.abspath(adir) + os.sep):
        abort(400, '非法结果路径')
    return check_path(p, must_exist=False, in_platform=True)


# contig 分析引擎：本地（离线内置库，默认） / 在线（NCBI）
_ENGINES = ('local', 'online')


def _norm_engine(raw):
    e = (raw or 'local').strip().lower()
    return e if e in _ENGINES else 'local'


def _engine_tag(engine):
    """结果文件名后缀：在线沿用旧名（兼容历史缓存），本地加 _local。"""
    return '' if engine == 'online' else '_' + engine


@bp.route('/api/tool/viral_contigs')
def api_tool_viral_contigs():
    """某次 contigs 运行的病毒 contig 分类表（metabuli 风格列）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    tsv = check_path(os.path.join(_run_real_dir(run),
                                  'virus_classification.tsv'),
                     must_exist=False, in_platform=True)
    rows = []
    if os.path.isfile(tsv):
        import csv
        with safe_open(tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                r['near_complete'] = r.get('near_complete') == 'True'
                rows.append(r)
    # 并入宿主预测列（若该运行已做宿主预测）
    hp = check_path(os.path.join(_run_real_dir(run),
                                 '08_host_analysis', 'host_prediction.tsv'),
                    must_exist=False, in_platform=True)
    host_map = {}                      # 先初始化：未做宿主预测的运行无该文件
    if os.path.isfile(hp):
        import csv
        host_map = {}
        with safe_open(hp) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                host_map[c] = h + (f' ({conf})' if conf else '')
        for r in rows:
            r['host'] = host_map.get((r.get('contig') or '').strip(), '')
    return jsonify(rows)


@bp.route('/api/tool/verify_result')
def api_tool_verify_result():
    """候选序列验证结果：calls.tsv 行 + summary.json 摘要。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    vdir = check_path(os.path.join(_run_real_dir(run), 'verify'),
                      must_exist=False, in_platform=True)
    rows, summary = [], {}
    calls_tsv = os.path.join(vdir, 'calls.tsv')
    if os.path.isfile(calls_tsv):
        import csv
        with safe_open(calls_tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                r['len'] = int(r.get('len') or 0)
                rows.append(r)
    summ_json = os.path.join(vdir, 'summary.json')
    if os.path.isfile(summ_json):
        with safe_open(summ_json) as f:
            summary = json.load(f)
    return jsonify({'rows': rows, 'summary': summary})


@bp.route('/api/tool/verify_export', methods=['POST'])
def api_tool_verify_export():
    """候选序列验证：勾选行导出序列（FASTA 下载流）。

    body: {run: 运行名, contigs: [contig 名]}。序列从该运行的
    verify/virus_input.fasta 与 verify/viroid_input.fasta 里抽取（两条
    支路的输入 fasta 已覆盖全部参与验证的 contig，无需回溯源 fastq）。
    """
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    contigs = [str(c).strip() for c in (body.get('contigs') or [])
               if str(c).strip()]
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    if not contigs:
        abort(400, '请先勾选要导出的行')
    vdir = check_path(os.path.join(_run_real_dir(run), 'verify'),
                      must_exist=False, in_platform=True)
    if not os.path.isdir(vdir):
        abort(404, f'运行 {run} 没有验证产物（verify/ 目录不存在）')
    want = set(contigs)
    kept, missing = [], set(want)

    def _push(rec_id, seq):
        missing.discard(rec_id)
        kept.append((rec_id, seq))

    for fname in ('virus_input.fasta', 'viroid_input.fasta'):
        fp = os.path.join(vdir, fname)
        if not os.path.isfile(fp):
            continue
        rec_id, chunks = None, []
        with safe_open(fp) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if rec_id in want:
                        _push(rec_id, ''.join(chunks))
                    rec_id, chunks = line[1:].split()[0], []
                elif rec_id is not None:
                    chunks.append(line)
            if rec_id in want:
                _push(rec_id, ''.join(chunks))
    if not kept:
        abort(404, '所选 contig 在该运行的序列文件中均未找到')
    fa = ''.join(f'>{rid}\n{seq}\n' for rid, seq in kept)
    resp = Response(fa, mimetype='application/octet-stream')
    resp.headers['Content-Disposition'] = (
        f'attachment; filename="{run}_selected.fasta"')
    return resp


@bp.route('/api/tool/kvsuite_refs')
def api_tool_kvsuite_refs():
    """列出某次已知病毒运行中可用作回贴参考的病毒基因组。

    来源（按优先级）：
      <run>/kvsuite/virus-fasta/ref_<acc>/ref_<acc>.ref.fasta   过滤后的病毒基因组
      <run>/kvsuite/virus-annotations/<acc>.gb / .gtf         注释
    供 t-consensus 选参考、t-variant 查注释用。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    import csv as _csv
    base = check_path(os.path.join(_run_real_dir(run), 'kvsuite'),
                      must_exist=False, in_platform=True)

    def _read_ids(path):
        ids = []
        if not os.path.isfile(path):
            return ids
        with safe_open(path) as f:
            for line in f:
                if line.startswith('>'):
                    ids.append(line[1:].strip().split()[0])
        return ids

    fa_dir = os.path.join(base, 'virus-fasta')
    ann_dir = os.path.join(base, 'virus-annotations')
    refs = []
    if os.path.isdir(fa_dir):
        for dn in sorted(os.listdir(fa_dir)):
            if not dn.startswith('ref_'):
                continue
            sub = os.path.join(fa_dir, dn)
            if not os.path.isdir(sub):
                continue
            fasta = None
            for fn in sorted(os.listdir(sub)):
                if fn.endswith(('.fasta', '.fa', '.fna')):
                    fasta = os.path.join(sub, fn)
                    break
            if not fasta:
                continue
            acc = dn[len('ref_'):]
            try:
                size = os.path.getsize(fasta)
            except OSError:
                size = 0
            gb = os.path.join(ann_dir, f'{acc}.gb')
            gtf = os.path.join(ann_dir, f'{acc}.gtf')
            gtf_n = 0
            if os.path.isfile(gtf):
                with safe_open(gtf) as f:
                    gtf_n = sum(1 for line in f if line.strip())
            refs.append({
                'accession': acc,
                'dir': os.path.relpath(sub, base).replace('\\', '/'),
                'fasta': os.path.relpath(fasta, base).replace('\\', '/'),
                'size': fmt_size(size),
                'seq_ids': _read_ids(fasta),
                'gb': os.path.isfile(gb),
                'gtf': os.path.isfile(gtf),
                'n_cds': gtf_n,
            })

    # BAM：比对产物本体（t-variant「与上游衔接」模式用 kvsuite/bam 作输入；
    # t-consensus 固定重新比对，不再复用）
    bams = []
    bam_dir = os.path.join(base, 'bam')
    if os.path.isdir(bam_dir):
        for fn in sorted(os.listdir(bam_dir)):
            if not fn.endswith('.bam'):
                continue
            p = os.path.join(bam_dir, fn)
            try:
                nbytes = os.path.getsize(p)
            except OSError:
                continue
            bams.append({
                'name': fn,
                'path': p,
                'size': fmt_size(nbytes),
                'size_mb': round(nbytes / 1048576.0, 1),
                'indexed': os.path.isfile(p + '.bai'),
            })

    # 过滤表：给前端展示物种名（末次鉴定结果）
    filt = {}
    ft = os.path.join(base, 'filter', 'filtered.tsv')
    if os.path.isfile(ft):
        with safe_open(ft) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                acc = (r.get('Accession') or '').strip()
                if acc and acc not in filt:
                    filt[acc] = {
                        'species': (r.get('Species') or '').strip(),
                        'coverage': (r.get('Coverage(%)') or '').strip(),
                        'depth': (r.get('MeanDepth') or '').strip(),
                    }
    for r in refs:
        r['species'] = filt.get(r['accession'], {}).get('species', '')
        r['coverage'] = filt.get(r['accession'], {}).get('coverage', '')
        r['depth'] = filt.get(r['accession'], {}).get('depth', '')

    return jsonify({
        'run': run,
        'n_refs': len(refs),
        'refs': refs,
        'bams': bams,
        'has_annotations': os.path.isdir(ann_dir),
    })


_kv_fig_lock = threading.Lock()


def _kv_ensure_figures(base):
    """确保 kvsuite/figures/ 下有识别/定量可视化图（freq/coabundance，
    绘图对象 = 通过过滤的结果）。

    旧运行没有图 -> 首次访问时从 filter/filtered.tsv 现算并落盘缓存；
    已齐全则直接返回。任何失败降级为返回当前已有的图。
    """
    fig_dir = os.path.join(base, 'figures')
    _names = (os.listdir(fig_dir) if os.path.isdir(fig_dir) else [])
    listed = sorted(f for f in _names if f.lower().endswith('.png'))
    have_svg = {f.lower() for f in _names if f.lower().endswith('.svg')}
    # 旧运行只有 PNG：也要补一次 SVG（结果图整理的原件），故把
    # freq 的 SVG 一并纳入「已齐全」判据。
    need = {'freq_multi_metrics.png'}
    if need <= set(listed) and 'freq_multi_metrics.svg' in have_svg:
        return listed
    tsv = os.path.join(base, 'filter', 'filtered.tsv')
    if not os.path.isfile(tsv):
        return listed
    with _kv_fig_lock:
        try:
            from Virus_Platform_Core.known_virus_suite.kv_fig import (
                generate_figures)
            listed = generate_figures(tsv, fig_dir, logger=_log)
        except Exception as e:  # noqa: BLE001
            _log.warning('kvsuite 图表生成失败 (%s): %s', base, e)
    return listed


@bp.route('/api/tool/kvsuite_result')
def api_tool_kvsuite_result():
    """known_virus_suite 结果：鉴定 / 过滤 / 共识 QC / 变异注释 四张表 + 图文件。

    产物目录 <run>/kvsuite/{identify,filter,consensus,variant,plots}
    """
    import csv as _csv
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    base = check_path(os.path.join(_run_real_dir(run), 'kvsuite'),
                      must_exist=False, in_platform=True)

    def _read_tsv(path, limit=None):
        rows = []
        if not os.path.isfile(path):
            return rows
        with safe_open(path) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                rows.append(r)
                if limit and len(rows) >= limit:
                    break
        return rows

    ns = {'identify': [], 'filtered': [], 'discarded': [],
          'consensus_qc': [], 'variants': []}
    ns['identify'] = _read_tsv(os.path.join(base, 'identify',
                                            'all_viruses.summary.tsv'), 500)
    ns['filtered'] = _read_tsv(os.path.join(base, 'filter', 'filtered.tsv'), 500)
    ns['discarded'] = _read_tsv(os.path.join(base, 'filter', 'discarded.tsv'), 200)
    ns['consensus_qc'] = _read_tsv(
        os.path.join(base, 'consensus', 'consensus_qc.tsv'), 500)

    # 变异注释：把所有 per-genome ann.tsv 合起来
    # 注意：kv_variant 把产物直接写在 <out> 下（<out>/annotated/、
    # <out>/variant_summary.json），没有 variant/ 中间层。
    ann_dir = os.path.join(base, 'annotated')
    # QUAL/DP/AF 不在 ann.tsv 里，需从同名的 vcf 补上（键 = CHROM|POS|ALT）
    vcf_map = _kv_load_vcf_metrics(os.path.join(base, 'vcf'))
    if os.path.isdir(ann_dir):
        for fn in sorted(os.listdir(ann_dir)):
            if not fn.endswith('.ann.tsv'):
                continue
            gid = fn[: -len('.ann.tsv')]
            for r in _read_tsv(os.path.join(ann_dir, fn), 300):
                r['genome'] = gid
                k = '{}|{}|{}'.format(r.get('CHROM', ''), r.get('POS', ''),
                                     r.get('ALT', ''))
                m = vcf_map.get(k)
                if m:
                    r['QUAL'] = m['QUAL']
                    r['DP'] = m['DP']
                    r['AF'] = m['AF']
                ns['variants'].append(r)

    # 深度图（<out>/plots/<sample>/*.png 有子目录，必须递归）
    plots = []
    pdir = os.path.join(base, 'plots')
    if os.path.isdir(pdir):
        for root, _dirs, files in os.walk(pdir):
            for fn in sorted(files):
                if fn.lower().endswith(('.png', '.pdf')):
                    rel = os.path.relpath(os.path.join(root, fn), pdir)
                    plots.append(rel.replace(os.sep, '/'))
        plots.sort()

    # 变异图（<out>/variant_plots/，平铺目录）
    variant_plots = []
    vpdir = os.path.join(base, 'variant_plots')
    if os.path.isdir(vpdir):
        for fn in sorted(os.listdir(vpdir)):
            if fn.lower().endswith(('.png', '.pdf')):
                variant_plots.append(fn)
    # 只展示 png（页面内联），pdf 供下载
    variant_plots_png = [p for p in variant_plots if p.lower().endswith('.png')]

    # 扩展变异分析（分子谱 / 变异密度 / 群体遗传），与 variant_plots 同目录平铺
    evo_plots_png = [p for p in variant_plots
                     if p.lower().endswith('.png') and '_evo_' in p]
    evo_manifest = {}
    _eman = os.path.join(base, 'variant_evo', 'evo_manifest.json')
    if os.path.isfile(_eman):
        try:
            with safe_open(_eman) as f:
                evo_manifest = json.load(f)
        except (ValueError, OSError):
            evo_manifest = {}

    # 群体遗传滑窗表 + SNPGenie 汇总（只取头几行供页面预览）
    popgen_tables = {}
    _edir = os.path.join(base, 'variant_evo')
    if os.path.isdir(_edir):
        for fn in sorted(os.listdir(_edir)):
            if fn.endswith('_popgen_window.tsv'):
                acc = fn[: -len('_popgen_window.tsv')]
                popgen_tables[acc] = _read_tsv(os.path.join(_edir, fn), 200)

    snpgenie_products = {}
    if os.path.isdir(_edir):
        for fn in sorted(os.listdir(_edir)):
            if fn.endswith('.snpgenie'):
                acc = fn[: -len('.snpgenie')]
                # 优先读带 ORF 变异计数的注释版（区分无变异/无多态）
                for pname in ('product_results.annotated.txt',
                              'product_results.txt'):
                    prod = os.path.join(_edir, fn, pname)
                    if os.path.isfile(prod):
                        snpgenie_products[acc] = _read_tsv(prod, 200)
                        break

    # 识别/定量可视化图表（参考管线 freq/coabundance 两组图，matplotlib
    # 出 PNG+PDF）。绘图对象 = 通过过滤（确认存在）的结果；过滤表为运行的
    # 不可变产物，图按需生成一次即缓存；失败只记日志，不影响结果接口。
    figures = _kv_ensure_figures(base)

    summary = {}
    vsum = os.path.join(base, 'variant_summary.json')
    if os.path.isfile(vsum):
        try:
            with safe_open(vsum) as f:
                summary['variant'] = json.load(f)
            # 脱敏：results 里的 vcf/ann_vcf/ann_tsv 是绝对路径，
            # 前端渲染只用 accession/genome/n_variants/n_coding_ann，
            # 不把本机路径暴露到页面。产物文件本身保持完整（含路径）。
            var = summary.get('variant') or {}
            for r in (var.get('results') or []):
                for k in ('vcf', 'ann_vcf', 'ann_tsv'):
                    r.pop(k, None)
        except (ValueError, OSError):
            pass
    summary['n_identify'] = len(ns['identify'])
    summary['n_filtered'] = len(ns['filtered'])
    summary['n_discarded'] = len(ns['discarded'])
    summary['n_consensus'] = len(ns['consensus_qc'])
    summary['n_variants'] = len(ns['variants'])
    summary['n_variant_plots'] = len(variant_plots_png)
    summary['n_evo_plots'] = len(evo_plots_png)
    summary['n_snpgenie'] = len(evo_manifest.get('snpgenie') or [])

    return jsonify({'run': run, 'summary': summary, 'plots': plots,
                    'figures': figures,
                    'variant_plots': variant_plots_png,
                    'evo_plots': evo_plots_png,
                    'evo_manifest': evo_manifest,
                    'popgen_tables': popgen_tables,
                    'snpgenie_products': snpgenie_products,
                    'identify': ns['identify'], 'filtered': ns['filtered'],
                    'discarded': ns['discarded'],
                    'consensus_qc': ns['consensus_qc'],
                    'variants': ns['variants']})


def _kv_load_vcf_metrics(vcf_dir):
    """扫 <out>/vcf/*.vcf 取 QUAL/DP/AF，键为 CHROM|POS|ALT。

    AF 口径与 kv_variant 过滤一致：AD[1]/INFO/DP（无 AD 时退 AC/AN）。
    单行格式异常只跳过该行并 warning，不丢整个 VCF 剩余记录。
    """
    out = {}
    if not os.path.isdir(vcf_dir):
        return out
    for fn in sorted(os.listdir(vcf_dir)):
        if not fn.lower().endswith('.vcf') or fn.startswith('all.'):
            continue
        n_bad = 0
        try:
            with safe_open(os.path.join(vcf_dir, fn)) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    try:
                        f2 = line.rstrip('\n').split('\t')
                        if len(f2) < 8:
                            n_bad += 1
                            continue
                        chrom, pos, _id, _ref, alt, qual, _filt, info = f2[:8]
                        dp = af = None
                        ad = ac = an = None
                        for item in info.split(';'):
                            if item.startswith('DP='):
                                dp = item[3:]
                            elif item.startswith('AD='):
                                ad = item[3:]
                            elif item.startswith('AC='):
                                ac = item[3:]
                            elif item.startswith('AN='):
                                an = item[3:]
                        dp_i = None
                        try:
                            dp_i = int(dp) if dp is not None else None
                        except (ValueError, TypeError):
                            dp_i = None
                        alt_ad = None
                        if ad:
                            parts = ad.split(',')
                            if len(parts) >= 2:
                                try:
                                    alt_ad = int(parts[1])
                                except (ValueError, TypeError):
                                    alt_ad = None
                        if alt_ad is not None and dp_i:
                            af = round(alt_ad / dp_i, 4)
                        elif ac and an:
                            try:
                                af = round(float(ac) / float(an), 4)
                            except (ValueError, TypeError, ZeroDivisionError):
                                af = None
                        out['{}|{}|{}'.format(chrom, pos, alt)] = {
                            'QUAL': qual, 'DP': dp_i, 'AF': af}
                    except Exception as e:
                        n_bad += 1
                        if n_bad <= 3:
                            _log.warning(
                                f'VCF 行解析失败 {fn} L{line[:60]!r}: {e}')
        except OSError as e:
            _log.warning(f'VCF 读取失败 {fn}: {e}')
            continue
        if n_bad:
            _log.warning(f'VCF {fn}: {n_bad} 行解析失败已跳过')
    return out


@bp.route('/api/tool/consensus_result')
def api_tool_consensus_result():
    """共识与变异结果：coverage.tsv + variants.tsv + summary.json + 分箱覆盖曲线。

    cov 为按 bin 平均的深度序列（最多 ~800 点），避免把上万个逐位深度
    全量塞给前端。
    """
    import csv as _csv
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    cdir = check_path(os.path.join(_run_real_dir(run), 'consensus'),
                      must_exist=False, in_platform=True)
    rows, variants, summary = [], [], {}

    cov_tsv = os.path.join(cdir, 'coverage.tsv')
    if os.path.isfile(cov_tsv):
        with safe_open(cov_tsv) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                for k in ('length', 'mapped_reads', 'covered_bp',
                          'n_variants', 'n_snv', 'n_isnv'):
                    try:
                        r[k] = int(r.get(k) or 0)
                    except ValueError:
                        r[k] = 0
                for k in ('median_depth', 'min_depth', 'max_depth'):
                    try:
                        r[k] = int(float(r.get(k) or 0))
                    except ValueError:
                        r[k] = 0
                for k in ('mean_depth', 'coverage_pct',
                          'depth_ok_pct', 'called_pct'):
                    try:
                        r[k] = float(r.get(k) or 0)
                    except ValueError:
                        r[k] = 0.0
                rows.append(r)

    var_tsv = os.path.join(cdir, 'variants.tsv')
    if os.path.isfile(var_tsv):
        with safe_open(var_tsv) as f:
            for i, r in enumerate(_csv.DictReader(f, delimiter='\t')):
                if i >= 500:
                    break
                variants.append(r)

    summ_json = os.path.join(cdir, 'summary.json')
    if os.path.isfile(summ_json):
        with safe_open(summ_json) as f:
            summary = json.load(f)

    cov = {}
    for r in rows:
        name = r.get('contig') or ''
        tag = re.sub(r'[^\w\-.]+', '_', name)[:120]
        cov[name] = _bin_coverage(os.path.join(cdir, tag + '.poscounts.tsv'),
                                  int(r.get('length') or 0))
    return jsonify({'rows': rows, 'variants': variants,
                    'summary': summary, 'cov': cov})


def _bin_coverage(pc_path, length, max_bins=800):
    """poscounts.tsv → 分箱平均深度序列。"""
    if not os.path.isfile(pc_path) or length <= 0:
        return []
    bin_sz = max(1, (length + max_bins - 1) // max_bins)
    nbin = (length + bin_sz - 1) // bin_sz
    sums = [0] * nbin
    cnts = [0] * nbin
    with safe_open(pc_path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            try:
                pos = int(p[0])
                depth = int(p[1]) + int(p[2]) + int(p[3]) + int(p[4]) + int(p[5])
            except ValueError:
                continue
            b = pos // bin_sz
            if 0 <= b < nbin:
                sums[b] += depth
                cnts[b] += 1
    return [round(sums[i] / cnts[i], 1) if cnts[i] else 0 for i in range(nbin)]


@bp.route('/api/tool/structcmp_data')
def api_tool_structcmp_data():
    """结构比较结果（JSON）：identity 矩阵 + 序列名，供页内 Plotly 热图渲染。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    js = check_path(os.path.join(_run_real_dir(run),
                                 'identity_matrix.json'),
                    must_exist=True, in_platform=True)
    return send_file(js, mimetype='application/json')


@bp.route('/api/tool/msa_data')
def api_tool_msa_data():
    """structcmp 运行 → SNP-only MSA 展示数据（比较基因组 MSA 查看卡用）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    aln = check_path(os.path.join(_run_real_dir(run), 'aln.fasta'),
                     must_exist=True, in_platform=True)
    from Virus_Platform_Core.msa_view import snp_view_data
    try:
        return jsonify(snp_view_data(aln))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/tool/quicktree_data')
def api_tool_quicktree_data():
    """quicktree 运行 → 树 Newick（进化树查看器卡片直接渲染）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(_run_real_dir(run),
                         must_exist=True, in_platform=True)
    tree_file = None
    for name in ('nj.nwk', 'tree.nwk'):
        p = os.path.join(run_dir, name)
        if os.path.isfile(p):
            tree_file = p
            break
    if not tree_file:
        abort(404, '该运行没有树文件（尚未完成或建树失败，看运行日志）')
    with safe_open(tree_file) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    base = os.path.basename(tree_file)
    return jsonify({'newick': nwk, 'file': base,
                    'tool': 'NJ' if base == 'nj.nwk' else 'FastTree',
                    # 「进化树 + 基因组叠加」：与树同名的 <stem>.overlay.json，
                    # 有就带上（无则 null，前端照常渲染纯树）。
                    'overlay': tree_overlay_for(tree_file)})


_ALIGN_EXTS = ('.fasta', '.fa', '.fna', '.fas', '.aln', '.txt')
# 比对查看器单文件上限（MB）：超过即拒绝，避免整份比对进内存 + 全量 JSON 下发
_ALIGN_MAX_MB = 128


def _parse_alignment_fasta(path):
    """读比对/序列 FASTA → dict（names/seqs/长度一致性/类型）。"""
    from Virus_Platform_Core.utils import iter_fasta
    from Virus_Platform_Core.sdt_exact import detect_seqtype
    names, seqs = [], []
    for h, s in iter_fasta(path):
        names.append(re.split(r'[\s|]', (h or '').strip())[0][:60] or
                     f'seq{len(names) + 1}')
        seqs.append(s.upper())
    if not seqs:
        raise ValueError('文件中没有序列')
    lens = {len(s) for s in seqs}
    return {'names': names, 'seqs': seqs, 'n': len(seqs),
            'cols': max(len(s) for s in seqs),
            'aligned': len(lens) == 1,
            'type': detect_seqtype(seqs)}


@bp.route('/api/align/file')
def api_align_file():
    """比对查看器数据：path= 平台内或本机（只读）FASTA 比对文件。

    超过 _ALIGN_MAX_MB 的文件直接拒绝：本接口把整份比对读进内存并序列化
    成 JSON（每条序列全量下发），没有上限时一个几 GB 的 FASTA 会让工作进程
    内存爆掉、请求也永远回不来。大比对请用样品的 MSA 查看器（按页下发）。
    """
    path = (request.args.get('path') or '').strip()
    if not path.lower().endswith(_ALIGN_EXTS):
        abort(400, '请提供 FASTA/比对文件（' + '/'.join(_ALIGN_EXTS[:4]) + '…）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    try:
        sz = os.path.getsize(p)
    except OSError as e:
        abort(400, f'无法读取文件信息: {e}')
    if sz > _ALIGN_MAX_MB * 1024 * 1024:
        abort(400, f'文件 {sz / 1048576:.0f} MB 超过比对查看器上限 '
                   f'{_ALIGN_MAX_MB} MB；请改用样品的「MSA 查看」'
                   f'（按页加载）或先清剪/抽样序列')
    try:
        d = _parse_alignment_fasta(p)
    except ValueError as e:
        abort(400, str(e))
    d['path'] = os.path.relpath(p, PLATFORM_ROOT) \
        if str(p).startswith(str(PLATFORM_ROOT)) else str(p)
    return jsonify(d)


@bp.route('/api/align/save', methods=['POST'])
def api_align_save():
    """编辑保存：content 写到源文件旁 <stem>.edit.fasta（平台内）；
    源在平台外时落 tool_runs/align_edit_<ts>/edited.fa。返回新路径。"""
    body = request.get_json(force=True) or {}
    path = (body.get('path') or '').strip()
    content = str(body.get('content') or '')
    if not path or not content.strip():
        abort(400, '参数不完整（path 与 content 必填）')
    if len(content) > _ALIGN_MAX_MB * 1024 * 1024:
        abort(400, f'提交内容 {len(content) / 1048576:.0f} MB 超过比对查看器上限 '
                   f'{_ALIGN_MAX_MB} MB')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    # 校验能按 FASTA 解析（至少 2 条、只含合法字符集宽容量）
    try:
        probe = _parse_alignment_fasta_string(content)
    except ValueError as e:
        abort(400, f'内容不是合法 FASTA: {e}')
    in_plat = str(p).startswith(str(PLATFORM_ROOT) + os.sep)
    if in_plat:
        stem, _ext = os.path.splitext(p)
        dst = stem + '.edit.fasta'
    else:
        dst = check_path(os.path.join(
            _tool_runs_root(), f"align_edit_{time.strftime('%Y%m%d_%H%M%S')}",
            'edited.fasta'), must_exist=False, in_platform=True)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
    with safe_open(dst, 'wt') as f:
        f.write(content if content.endswith('\n') else content + '\n')
    return jsonify({'saved': dst, 'aligned': probe['aligned'],
                    'n': probe['n'], 'cols': probe['cols']})


def _parse_alignment_fasta_string(content):
    """与 _parse_alignment_fasta 同口径，但直接吃文本。"""
    from Virus_Platform_Core.sdt_exact import detect_seqtype
    names, seqs, cur = [], [], None
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if cur is not None:
                seqs.append(''.join(cur))
            names.append(re.split(r'[\s|]', line[1:].strip())[0][:60] or
                         f'seq{len(names) + 1}')
            cur = []
        elif cur is not None:
            cur.append(line)
        else:
            raise ValueError('第一条记录前出现序列行')
    if cur is not None:
        seqs.append(''.join(cur))
    if len(seqs) < 1 or not names:
        raise ValueError('未解析到序列')
    joined = [s.upper() for s in seqs]
    if any(not s for s in joined):
        raise ValueError('存在空序列')
    return {'names': names, 'seqs': joined, 'n': len(joined),
            'cols': max(len(s) for s in joined),
            'aligned': len({len(s) for s in joined}) == 1,
            'type': detect_seqtype(joined)}


_TREE_EXTS = ('.nwk', '.newick', '.treefile', '.contree', '.tre', '.tree')


@bp.route('/api/tree/file')
def api_tree_file():
    """树文件路径 → Newick 文本（进化树查看卡「本机树文件」输入用）。

    路径可以是平台内相对路径或本机绝对路径（只读），扩展名白名单限制。
    """
    path = (request.args.get('path') or '').strip()
    if not path.lower().endswith(_TREE_EXTS):
        abort(400, '请提供树文件（'
                    + ' / '.join(_TREE_EXTS) + '）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    with safe_open(p) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    return jsonify({'newick': nwk, 'file': os.path.basename(p),
                    'overlay': tree_overlay_for(p)})


@bp.route('/api/tool/report')
def api_tool_report():
    """生成（带缓存）并打开工具运行的结果报告（桑基/旭日/分类表）。

    ?force=1 强制重建（输入未变时默认复用既有 report.html）。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    force = request.args.get('force') in ('1', 'true', 'yes')
    run_dir = check_path(_run_real_dir(run),
                         must_exist=True, in_platform=True)
    from Virus_Platform_Core.viz import build_contig_report
    try:
        build_contig_report(run_dir, force=force)
    except RuntimeError as e:
        abort(400, str(e))
    from flask import redirect
    return redirect(f'/tool_runs/{run}/report.html')


@bp.route('/api/tool/report_data')
def api_tool_report_data():
    """分类报告数据（JSON）：供工具页内联渲染桑基/旭日/分类表/明细表。

    contigs_* 运行 → mode=contig（含 contig 明细 + 宿主列）；
    identify_* 运行 → mode=reads（含 classified sequences 表）。
    旭日图产物按输入 mtime 复用（本接口在工具页一次访问里会被调 3 次）。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    force = request.args.get('force') in ('1', 'true', 'yes')
    run_dir = check_path(_run_real_dir(run),
                         must_exist=True, in_platform=True)
    import glob as _glob
    krs = (_glob.glob(os.path.join(run_dir, '*.kreport2')) +
           _glob.glob(os.path.join(run_dir, '**', '*.kreport2'),
                      recursive=True))
    if not krs:
        abort(400, '运行目录中没有 kreport 分类报告')
    from Virus_Platform_Core.kunpeng import parse_kreport
    krows = parse_kreport(sorted(krs)[0])
    report_rows = [{'proportion': r['percent'], 'count': int(r['frags']),
                    'rank': r['rank'], 'taxid': r['taxid'],
                    'name': r['name'], 'depth': r['depth']}
                   for r in krows]

    import csv as _csv
    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    ids_tsv = os.path.join(run_dir, 'viral_ids.tsv')
    mode, vc_rows, ids_rows = ('contig', [], [])
    if os.path.isfile(tsv):
        mode = 'contig'
        with safe_open(tsv) as f:
            vc_rows = list(_csv.DictReader(f, delimiter='\t'))
    elif os.path.isfile(ids_tsv):
        mode = 'reads'
        with safe_open(ids_tsv) as f:
            ids_rows = list(_csv.DictReader(f, delimiter='\t'))

    host_map, host_cats = {}, {}
    hp = os.path.join(run_dir, '08_host_analysis', 'host_prediction.tsv')
    if os.path.isfile(hp):
        with safe_open(hp) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                if c:
                    host_map[c] = h + (f' ({conf})' if conf else '')
                if h:
                    host_cats[h] = host_cats.get(h, 0) + 1
    for r in vc_rows:
        r['host'] = host_map.get((r.get('contig') or '').strip(), '')

    n_total, n_virus = 0, 0
    if mode == 'contig':
        cfa = os.path.join(run_dir, 'contigs.filtered.fasta')
        if os.path.isfile(cfa):
            from Virus_Platform_Core.utils import count_fasta_seqs
            n_total = count_fasta_seqs(cfa)
        n_virus = len(vc_rows)
    else:
        n_total = n_virus = len(ids_rows)

    host_list = [{'cat': k, 'n': v} for k, v in
                 sorted(host_cats.items(), key=lambda x: -x[1])]

    # taxburst 旭日图（Krona 等价物，离线交互 HTML）
    from Virus_Platform_Core.viz import build_taxburst
    tb_ok, tb_unc = build_taxburst(run_dir, krows, force=force)

    dl = []
    if mode == 'contig':
        dl = [('viral_contigs.fasta?dl=1', 'Virus Sequences (FASTA)'),
              ('virus_classification.tsv?dl=1', 'Virus Classification (TSV)'),
              ('contigs.filtered.fasta?dl=1', 'All Contigs (FASTA)'),
              (os.path.basename(sorted(krs)[0]) + '?dl=1', 'Report (kreport TSV)')]
        if os.path.isfile(hp):
            dl.append(('08_host_analysis/host_prediction.tsv?dl=1',
                       'Host Prediction (TSV)'))
    else:
        import glob as _g2
        for fa in sorted(_g2.glob(os.path.join(run_dir, 'viral_sequences*.fasta'))):
            dl.append((os.path.basename(fa) + '?dl=1', 'Viral Sequences (FASTA)'))
        if os.path.isfile(ids_tsv):
            dl.append(('viral_ids.tsv?dl=1', 'Classified IDs (TSV)'))
        dl.append((os.path.basename(sorted(krs)[0]) + '?dl=1',
                   'Report (kreport TSV)'))
    downloads = [{'href': href, 'label': label} for href, label in dl
                 if os.path.isfile(os.path.join(run_dir, href.split('?')[0]))]

    return jsonify({'run': run, 'mode': mode,
                    'report': report_rows, 'vc': vc_rows, 'ids': ids_rows,
                    'host': host_list, 'host_count': len(host_map),
                    'n_total': n_total, 'n_virus': n_virus,
                    'downloads': downloads,
                    'taxburst': tb_ok})


@bp.route('/api/tool/hostpredict_run', methods=['POST'])
def api_tool_hostpredict_run():
    """对既有 contigs 运行目录原地做 ICTV 宿主预测（结果并入该运行）。"""
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(_run_real_dir(run),
                         must_exist=True, in_platform=True)
    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    if not os.path.isfile(tsv):
        abort(400, '该运行没有 virus_classification.tsv（仅 contig 分类运行支持）')
    # 可选置信度门槛（0-1 小数）；不传/空 = 不过滤（历史行为）
    min_conf = body.get('min_confidence')
    if min_conf in (None, ''):
        min_conf = None
    else:
        try:
            min_conf = float(min_conf)
        except (TypeError, ValueError):
            abort(400, 'min_confidence 需为 0-1 之间的小数')
        if not 0.0 <= min_conf <= 1.0:
            abort(400, 'min_confidence 需为 0-1 之间的小数')
    from Virus_Platform_Core.kunpeng import db_ready
    if not db_ready(cfg.databases['host']):
        abort(400, '宿主库未就绪，请先到「数据库构建」页构建宿主库')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.05, '整理输入')
        # 统一 11 列口径（blast_* 用本地 BLASTN 现算补齐）
        from Virus_Platform_Core.host_analysis import normalize_contig_table, predict_hosts
        vfa = os.path.join(run_dir, 'viral_contigs.fasta')
        normalize_contig_table(tsv, a_dir,
                               fasta=vfa if os.path.isfile(vfa) else None,
                               threads=cfg.threads, logger=logger)
        prog('predict', 0.15, 'ICTV 宿主概率级联预测')
        res = predict_hosts(run_dir, logger=logger, force=True,
                            min_confidence=min_conf)
        prog('done', 1.0, f"完成：宿主判定 {res.get('n_contigs', 0)} 条")
        logger.close()
        return res

    ts = time.strftime('%Y%m%d_%H%M%S')
    tid = tm.start(cfg.tr(f'宿主预测·{run} {ts}', f'Host prediction·{run} {ts}'),
                   job, log_file=os.path.join(run_dir, 'hostpredict.log'))
    return jsonify({'task': tid})


@bp.route('/api/tool/analysis')
def api_tool_analysis():
    """读取已完成的 contig 分析结果（JSON 缓存）。"""
    run = request.args.get('run') or ''
    contig = request.args.get('contig') or ''
    action = request.args.get('action') or ''
    seq = (request.args.get('seq') or '').strip()
    engine = _norm_engine(request.args.get('engine'))
    tag = _engine_tag(engine)
    if not contig or action not in ('blastn', 'blastx', 'primer'):
        abort(400, '参数不完整')
    # 直接用字符串序列查询（seq 输入的伪 run 缓存）
    if seq:
        import hashlib as _hl
        _key = _hl.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
        _safe = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
        run_dir = os.path.join(_tool_runs_root(), '_seq_input',
                               _safe + '_' + _key + tag)
        run_dir = check_path(run_dir, must_exist=True, in_platform=True)
        p = _analysis_out(run_dir, _safe + '_' + _key + tag, action)
    else:
        if not run or not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        run_dir = check_path(os.path.join(_tool_runs_root(), run),
                             must_exist=True, in_platform=True)
        p = _analysis_out(run_dir, contig + tag, action)
    if not os.path.isfile(p):
        abort(404, '该分析尚未运行')
    with safe_open(p) as f:
        return Response(f.read(), mimetype='application/json')


@bp.route('/api/tool/analyze', methods=['POST'])
def api_tool_analyze():
    """单 contig 深度分析三件套：blastn / blastx / primer。

    body: {run, contig, action, seq?, engine?}
      engine=local（默认）→ 本地内置库（Virus_Platform_Core.local_search：blastn 病毒参考库 /
        DIAMOND vs viral_prot），离线、秒级；
      engine=online → NCBI 在线（Virus_Platform_Core.contig_annot，域名白名单 + IP 校验）。
    （曾有的第四件 cdd（CDD 保守域搜索）已于 2026-09-29 随本地 CDD 库下线。）
    结果缓存为 run/analysis/<contig>[_local]_<action>.json。
    """
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    contig = body.get('contig') or ''
    action = body.get('action') or ''
    seq = (body.get('seq') or '').strip()
    in_seq = seq  # 供闭包安全读取（避免 UnboundLocalError）
    engine = _norm_engine(body.get('engine'))
    tag = _engine_tag(engine)
    if action not in ('blastn', 'blastx', 'primer'):
        abort(400, '无效的分析类型')
    if not run or not contig:
        abort(400, '参数不完整')
    # 直接输入序列模式（无需 run），或从运行目录取 contig 序列
    if seq:
        # 用共享伪运行缓存，避免污染具体工具运行目录
        run_dir = check_path(os.path.join(_tool_runs_root(), '_seq_input'),
                             must_exist=False, in_platform=True)
        import hashlib as _hl
        _key = _hl.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
        contig = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
        _used_for_cache = contig + '_' + _key + tag
        run_dir = os.path.join(run_dir, _used_for_cache)
        os.makedirs(run_dir, exist_ok=True)
        out_path = _analysis_out(run_dir, _used_for_cache, action)
    else:
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        run_dir = check_path(os.path.join(_tool_runs_root(), run),
                             must_exist=True, in_platform=True)
        out_path = _analysis_out(run_dir, contig + tag, action)
    if os.path.isfile(out_path):
        # 缓存命中也要探针：任务被杀/磁盘写满可能留下截断的坏 JSON，
        # 只判 isfile 就返回 cached=True，前端会永远读到坏结果且不会重算。
        try:
            with safe_open(out_path) as f:
                json.loads(f.read())
            return jsonify({'task': None, 'cached': True})
        except (OSError, ValueError):
            try:
                os.remove(out_path)          # 坏缓存删除 → 下面正常重算
            except OSError:
                pass

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import contig_annot as ca
        prog('prep', 0.1, '查找 contig 序列')
        if in_seq:
            # 清洗直接输入的序列：去掉 FASTA header 与换行，仅留核酸字母
            seq = ''.join(line.strip().upper() for line in in_seq.splitlines()
                          if line.strip() and not line.strip().startswith('>'))
            if not seq:
                raise RuntimeError('输入序列为空')
        else:
            seq = _find_contig_seq(run_dir, contig)
            if not seq:
                raise RuntimeError(f'运行目录中未找到 contig: {contig}')

        if action == 'primer':
            from Virus_Platform_Core.primer import design_primers_for_seq
            prog('primer3', 0.5, '本地 primer3 设计引物')
            pairs = design_primers_for_seq(contig, seq, num_return=5)
            result = {'action': 'primer', 'contig': contig,
                      'length': len(seq), 'primers': pairs}

        elif engine == 'local':
            from Virus_Platform_Core import local_search as ls
            qfa = os.path.join(run_dir, f'query_{action}.fasta')
            with safe_open(qfa, 'wt') as f:
                f.write(f'>{contig}\n{seq}\n')
            prog(action, 0.4, ('本地 blastn vs 病毒参考核酸库（离线）'
                               if action == 'blastn' else
                               '本地 DIAMOND blastx vs RefSeq 病毒蛋白库（离线）'))
            fn = ls.blastn_local if action == 'blastn' else ls.blastx_local
            hits = fn(qfa, run_dir, threads=cfg.threads, logger=logger)
            result = {'action': action, 'contig': contig, 'engine': 'local',
                      'hits': hits}

        else:  # blastn / blastx（NCBI URL API，virus-restricted）
            prog('submit', 0.3, f'提交 NCBI {action}（数分钟，耐心等待）')
            rid, _rtoe = ca.submit_blast(
                action, 'nt' if action == 'blastn' else 'nr',
                f'>{contig}\n{seq}\n', entrez_query='viruses[Organism]')
            prog('poll', 0.5, f'NCBI 任务 {rid} 运行中')
            hits = ca.poll_blast(rid, cancel=cancel)
            result = {'action': action, 'contig': contig, 'engine': 'online',
                      'rid': rid, 'hits': hits}

        # 原子落盘：先写 .tmp 再 os.replace，避免读到写了一半的 JSON
        # （本文件既是缓存也是 /api/tool/analysis 的直接读源，非原子写会
        # 让并发读拿到截断内容）。
        tmp_out = out_path + '.tmp'
        with safe_open(tmp_out, 'wt') as f:
            json.dump(result, f, ensure_ascii=False, indent=1)
        os.replace(tmp_out, out_path)
        prog('done', 1.0, '完成')
        logger.close()
        return {'out': out_path}

    tid = tm.start(cfg.tr(f'Contig分析·{action} {contig[:30]}',
                         f'Contig {action} {contig[:30]}'), job,
                   log_file=os.path.join(run_dir, 'run.log'))
    return jsonify({'task': tid, 'cached': False})


# ------------------------------------------------------------------
# 从序列出发的引物设计（PCR / qPCR + 探针 + 热力学评估）
# 脱靶检查默认关（specificity=true 且已配置 host-db 时才跑宿主 blastn）；
# 覆盖图数据为百分比定位，前端纯 div 渲染
# ------------------------------------------------------------------
@bp.route('/api/tool/primer_design', methods=['POST'])
def api_tool_primer_design():
    """POST 参数：

    {
      "fasta": "输入 FASTA 路径（与 fasta_text 二选一）",
      "fasta_text": "FASTA 文本内容",
      "type": "PCR" | "qPCR",
      "core": {核心层参数},
      "adv": {进阶层参数},
      "specificity": false  // 可选：宿主基因组 blastn 脱靶检查（默认关）
    }

    同步返回设计结果（单序列通常 <2s），不建后台任务。
    """
    from Virus_Platform_Core import primer_design as pd
    d = request.get_json(silent=True) or {}
    fasta = (d.get('fasta') or '').strip()
    fasta_text = d.get('fasta_text') or ''
    ptype = (d.get('type') or 'PCR').strip()
    core = d.get('core') or {}
    adv = d.get('adv') or {}
    # 默认关；打开时透传 primer.check_specificity_blast（需已配置 host-db），
    # 结果/summary 如实写 specificity_checked（2026-10-03 修复 P2-14④）
    specificity = bool(d.get('specificity'))
    if not isinstance(core, dict):
        core = {}
    if not isinstance(adv, dict):
        adv = {}
    if not fasta and not fasta_text.strip():
        return jsonify({'error': '请提供输入 FASTA（文件路径或序列内容）'}), 400
    fpath = ''
    if fasta:
        try:
            fpath = _safe_sample_fasta(fasta)
        except ValueError as e:
            return jsonify({'error': str(e)}), 400
    try:
        res = pd.run_design(fasta_text=fasta_text, fasta_path=fpath,
                            ptype=ptype, core=core, adv=adv,
                            specificity=specificity)
    except FileNotFoundError as e:
        return jsonify({'error': str(e)}), 404
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        _log.exception('引物设计失败')
        return jsonify({'error': '设计失败: %s' % e}), 500
    if res['summary']['n_pairs'] == 0:
        explain = ''
        for rec in res['records']:
            if rec.get('explain'):
                explain = rec['explain']
                break
        res['warning'] = '未设计出任何引物对'
        res['explain'] = explain
    return jsonify(res)


def _safe_sample_fasta(path):
    """校验 FASTA 路径可读，返回规范化路径。"""
    if not path:
        raise ValueError('请提供 FASTA 路径')
    p = os.path.expanduser(path.strip().strip('"'))
    if not os.path.isabs(p):
        p = os.path.join(PLATFORM_ROOT, p)
    p = os.path.abspath(p)
    if not os.path.isfile(p):
        raise ValueError('文件不存在: %s' % p)
    return p
