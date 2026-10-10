# -*- coding: utf-8 -*-
"""
从序列出发的引物设计（PCR / qPCR）

与 Virus_Platform_Core/primer.py 的区别
----------------------
Virus_Platform_Core/primer.py       : 管线阶段⑥，基于样本目录 + 多序列比对/组装产物，写 primers.tsv
Virus_Platform_Core/primer_design.py: 交互式单序列设计，输入 FASTA 内容或文件，输出覆盖图 + 表格数据
                     不写管线产物，不改动 06_primer 契约

能力
----
- PCR  ：设计引物对，无探针
- qPCR ：设计引物对 + 探针（TaqMan），探针 Tm 自动约束到引物 Tm +5~+10℃
- 热力学评估：二聚体 / 发夹 / 3' 端稳定性 / GC / Tm / 探针专项，归一化到 100
- 不评估 BLAST 脱靶（无宿主库依赖）；不评估覆盖度（单序列场景恒为满分）
"""
import io
import os
import re

# primer3 统一走 primer3_runtime（子进程探测 + 惰性导入 + bindings 校验）：
# 模块级 import 在内存不足的机器上会让 pyd 初始化段错误杀掉整个进程，
# 不可用时 _p3() 抛 RuntimeError，由调用方给出友好提示并软跳过。
from . import primer3_runtime as _p3rt


def _p3():
    try:
        return _p3rt.get_module()
    except RuntimeError as e:
        raise RuntimeError(
            'primer3 不可用：T7 引物/热力学打分在本环境未启用（' + str(e) + '）')


from . import primer_thermo as pt

# ============================================================
# 默认参数（核心层 / 进阶层）
# ============================================================
# 核心层：用户最常调的 6 项
CORE_DEFAULTS = {
    'primer_opt_size': 20,
    'primer_min_size': 18,
    'primer_max_size': 24,
    'primer_opt_tm': 60.0,
    'primer_min_tm': 57.0,
    'primer_max_tm': 63.0,
    'primer_min_gc': 40.0,
    'primer_max_gc': 60.0,
    'product_min': 100,
    'product_max': 1000,
    'num_return': 5,
    'salt_monovalent': 50.0,
    'salt_divalent': 1.5,
    'dntp_conc': 0.6,
    'dna_conc': 50.0,
}

# 进阶层：默认折叠，需要时展开
ADV_DEFAULTS = {
    'max_poly_x': 4,
    'max_self_any_th': 47.0,      # 引物自二聚体 Tm 上限
    'max_self_end_th': 47.0,      # 引物 3' 端自配对 Tm 上限
    'max_hairpin_th': 47.0,       # 发夹 Tm 上限
    'max_end_stability': 9.0,     # 3' 端稳定性上限（primer3 分值）
    'probe_opt_size': 22,
    'probe_min_size': 18,
    'probe_max_size': 27,
    'probe_opt_tm': 67.0,         # 关键：比引物 opt_tm 高约 7℃
    'probe_min_tm': 65.0,
    'probe_max_tm': 73.0,
    'probe_min_gc': 40.0,
    'probe_max_gc': 70.0,
    'region_start': 0,            # 限定设计区域（1-based，0 = 全长）
    'region_end': 0,
}

# qPCR 产物长度默认更短（扩增效率）
QPCR_PRODUCT_RANGE = (80, 300)

# 探针 Tm 差要求（step3 / varVAMP 口径）
PROBE_TM_DIFF_MIN = 5.0
PROBE_TM_DIFF_MAX = 10.0


def _f(v, default=0.0):
    """安全转 float。"""
    try:
        if v is None or v == '':
            return float(default)
        return float(v)
    except (TypeError, ValueError):
        return float(default)


def _i(v, default=0):
    try:
        if v is None or v == '':
            return int(default)
        return int(float(v))
    except (TypeError, ValueError):
        return int(default)


def build_p3_options(ptype='PCR', core=None, adv=None):
    """把前端参数组装成 primer3 选项字典。

    ptype : 'PCR' | 'qPCR'
    core  : 核心层覆盖值
    adv   : 进阶层覆盖值
    """
    c = dict(CORE_DEFAULTS)
    a = dict(ADV_DEFAULTS)
    for k, v in (core or {}).items():
        if k in c:
            c[k] = v
    for k, v in (adv or {}).items():
        if k in a:
            a[k] = v

    pmin, pmax = _i(c['product_min'], 100), _i(c['product_max'], 1000)
    if ptype == 'qPCR' and (core is None or 'product_min' not in (core or {})):
        pmin, pmax = QPCR_PRODUCT_RANGE
    if pmin < 40:
        pmin = 40
    if pmax <= pmin:
        pmax = pmin + 100

    opts = {
        'PRIMER_OPT_SIZE': _i(c['primer_opt_size'], 20),
        'PRIMER_MIN_SIZE': _i(c['primer_min_size'], 18),
        'PRIMER_MAX_SIZE': _i(c['primer_max_size'], 24),
        'PRIMER_OPT_TM': _f(c['primer_opt_tm'], 60.0),
        'PRIMER_MIN_TM': _f(c['primer_min_tm'], 57.0),
        'PRIMER_MAX_TM': _f(c['primer_max_tm'], 63.0),
        'PRIMER_MIN_GC': _f(c['primer_min_gc'], 40.0),
        'PRIMER_MAX_GC': _f(c['primer_max_gc'], 60.0),
        'PRIMER_MAX_POLY_X': _i(a['max_poly_x'], 4),
        'PRIMER_PRODUCT_SIZE_RANGE': [[pmin, pmax]],
        'PRIMER_NUM_RETURN': max(1, min(_i(c['num_return'], 5), 50)),
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 1 if ptype == 'qPCR' else 0,
        'PRIMER_SALT_MONOVALENT': _f(c['salt_monovalent'], 50.0),
        'PRIMER_SALT_DIVALENT': _f(c['salt_divalent'], 1.5),
        'PRIMER_DNTP_CONC': _f(c['dntp_conc'], 0.6),
        'PRIMER_DNA_CONC': _f(c['dna_conc'], 50.0),
    }
    # 进阶阈值：只在与默认不同的放进去（primer3 有个别取值不接受）
    for key, opt in (('max_self_any_th', 'PRIMER_MAX_SELF_ANY_TH'),
                     ('max_self_end_th', 'PRIMER_MAX_SELF_END_TH'),
                     ('max_hairpin_th', 'PRIMER_MAX_HAIRPIN_TH')):
        val = _f(a[key], 47.0)
        if val > 0:
            opts[opt] = val
    # 3' 端稳定性上限若为 0 表示不限制
    es = _f(a.get('max_end_stability'), 0.0)
    if es > 0:
        opts['PRIMER_MAX_END_STABILITY'] = es

    if ptype == 'qPCR':
        opts.update({
            'PRIMER_INTERNAL_OPT_SIZE': _i(a['probe_opt_size'], 22),
            'PRIMER_INTERNAL_MIN_SIZE': _i(a['probe_min_size'], 18),
            'PRIMER_INTERNAL_MAX_SIZE': _i(a['probe_max_size'], 27),
            'PRIMER_INTERNAL_OPT_TM': _f(a['probe_opt_tm'], 67.0),
            'PRIMER_INTERNAL_MIN_TM': _f(a['probe_min_tm'], 65.0),
            'PRIMER_INTERNAL_MAX_TM': _f(a['probe_max_tm'], 73.0),
            'PRIMER_INTERNAL_MIN_GC': _f(a['probe_min_gc'], 40.0),
            'PRIMER_INTERNAL_MAX_GC': _f(a['probe_max_gc'], 70.0),
        })
        # 自洽校验：OPT 必须落在 [MIN, MAX] 内
        po = opts['PRIMER_INTERNAL_OPT_TM']
        pn = opts['PRIMER_INTERNAL_MIN_TM']
        px = opts['PRIMER_INTERNAL_MAX_TM']
        if pn > po:
            opts['PRIMER_INTERNAL_MIN_TM'] = po
        if px < po:
            opts['PRIMER_INTERNAL_MAX_TM'] = po
    return opts


# ============================================================
# FASTA 解析
# ============================================================
def parse_fasta_text(text):
    """解析 FASTA 文本，返回 [(name, seq), ...]。支持多记录。"""
    out = []
    name, buf = None, []
    for line in (text or '').splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                out.append((name, ''.join(buf).upper()))
            name, buf = line[1:].strip() or 'seq%d' % (len(out) + 1), []
        else:
            buf.append(re.sub(r'[^A-Za-z]', '', line))
    if name is not None:
        out.append((name, ''.join(buf).upper()))
    return [(n, s) for n, s in out if s]


def _read_fasta_path(path):
    if not path:
        return []
    if not os.path.isfile(path):
        raise FileNotFoundError('FASTA 文件不存在: %s' % path)
    if path.lower().endswith('.gz'):
        import gzip
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
            return parse_fasta_text(f.read())
    with io.open(path, encoding='utf-8', errors='replace') as f:
        return parse_fasta_text(f.read())


# ============================================================
# 单条序列设计
# ============================================================
def design_for_sequence(name, seq, ptype='PCR', core=None, adv=None):
    """对一条序列设计引物对（含探针，qPCR 时）。

    返回 {'name', 'length', 'pairs': [...], 'explain': str, 'p3_opts': dict}
    """
    seq = (seq or '').upper()
    adv = dict(adv or {})
    region_s = _i(adv.get('region_start'), 0)
    region_e = _i(adv.get('region_end'), 0)
    offset = 0
    subseq = seq
    if region_s > 0 and region_e > region_s and region_e <= len(seq):
        offset = region_s - 1
        subseq = seq[offset:region_e]
    elif region_s > 0 and region_e == 0 and region_s < len(seq):
        offset = region_s - 1
        subseq = seq[offset:]

    opts = build_p3_options(ptype, core, adv)
    # 短序列时收缩产物区间，避免 primer3 直接报错
    pmin, pmax = opts['PRIMER_PRODUCT_SIZE_RANGE'][0]
    if len(subseq) < pmax:
        nmax = max(60, min(len(subseq), pmax))
        nmin = min(pmin, max(40, nmax - 20))
        opts['PRIMER_PRODUCT_SIZE_RANGE'] = [[nmin, nmax]]

    res = _p3().bindings.design_primers(
        {'SEQUENCE_ID': (name or 'seq')[:80], 'SEQUENCE_TEMPLATE': subseq}, opts)

    n = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
    pairs = []
    for i in range(n):
        f_seq = res.get('PRIMER_LEFT_%d_SEQUENCE' % i, '')
        r_seq = res.get('PRIMER_RIGHT_%d_SEQUENCE' % i, '')
        p_seq = res.get('PRIMER_INTERNAL_%d_SEQUENCE' % i, '') if ptype == 'qPCR' else ''
        f_pos = res.get('PRIMER_LEFT_%d' % i, (0, 0))
        r_pos = res.get('PRIMER_RIGHT_%d' % i, (0, 0))
        p_pos = res.get('PRIMER_INTERNAL_%d' % i, (0, 0)) if ptype == 'qPCR' else (0, 0)
        f_tm = _f(res.get('PRIMER_LEFT_%d_TM' % i), 0.0)
        r_tm = _f(res.get('PRIMER_RIGHT_%d_TM' % i), 0.0)
        p_tm = _f(res.get('PRIMER_INTERNAL_%d_TM' % i), 0.0)

        # 热力学评估
        ev = pt.evaluate_pair(
            f_seq, r_seq, p_seq, ptype, f_tm, r_tm, p_tm,
            _f(res.get('PRIMER_LEFT_%d_END_STABILITY' % i), 0.0),
            _f(res.get('PRIMER_RIGHT_%d_END_STABILITY' % i), 0.0))

        f_start = f_pos[0] + 1 + offset
        r_start = r_pos[0] + 1 + offset
        rec = {
            'idx': i + 1,
            'type': ptype,
            'f_seq': f_seq,
            'f_start': f_start,
            'f_len': f_pos[1],
            'f_tm': round(f_tm, 1),
            'f_gc': round(_f(res.get('PRIMER_LEFT_%d_GC_PERCENT' % i), 0.0), 1),
            'r_seq': r_seq,
            'r_start': r_start,
            'r_len': r_pos[1],
            'r_tm': round(r_tm, 1),
            'r_gc': round(_f(res.get('PRIMER_RIGHT_%d_GC_PERCENT' % i), 0.0), 1),
            'product': _i(res.get('PRIMER_PAIR_%d_PRODUCT_SIZE' % i), 0),
            'product_tm': round(_f(res.get('PRIMER_PAIR_%d_PRODUCT_TM' % i), 0.0), 1),
            'penalty': round(_f(res.get('PRIMER_PAIR_%d_PENALTY' % i), 0.0), 3),
            'probe_seq': p_seq,
            'probe_start': (p_pos[0] + 1 + offset) if p_seq else 0,
            'probe_len': p_pos[1] if p_seq else 0,
            'probe_tm': round(p_tm, 1) if p_seq else 0.0,
            'probe_gc': round(_f(res.get('PRIMER_INTERNAL_%d_GC_PERCENT' % i), 0.0), 1) if p_seq else 0.0,
            'amp_start': f_start,
            'amp_end': r_start,
            # primer3 的 PRIMER_RIGHT_i[0] 就是扩增子在正链上的最右碱基
            # （0-based），PRODUCT_SIZE == r_start - f_start + 1（实测：
            # left=[33,20] right=[351,20] product=319）。原实现写成
            # r_start + r_len - 1，会让 amp_end 超出真实 3' 端 r_len-1 bp、
            # span 比 PRODUCT_SIZE 大 r_len-1，前端同时显示两个互斥数字。
            'amplicon_span': r_start - f_start + 1,
            'score': ev['score'],
            'recommendation': ev['recommendation'],
            'thermo': ev,
        }
        pairs.append(rec)

    return {
        'name': name,
        'length': len(seq),
        'design_length': len(subseq),
        'region_offset': offset,
        'pairs': pairs,
        'explain': res.get('PRIMER_PAIR_EXPLAIN', '') or res.get('PRIMER_ERROR', ''),
        'p3_opts': {k: v for k, v in opts.items()
                    if k.startswith('PRIMER_') and not k.endswith('EXPLAIN_FLAG')},
    }


# ============================================================
# 覆盖图数据（百分比定位 + bp 坐标，供前端交互缩放渲染）
# ============================================================
def build_track(seq_len, pairs, seq=''):
    """生成覆盖图数据结构。

    定位算法（已对拍实测）：
      left  = 起点 / 序列长度 * 100 %
      width = 长度 / 序列长度 * 100 %
      错层  top = 2 + 14 * n（n 为第 n 条轨道，从 0 起）
      产物条 top = 竖线 top + 4

    除百分比外，每条轨道附带 1-based 闭区间 bp 坐标（f/r/p/amp_*），前端
    交互缩放视图按 bp 重新投影，不再依赖百分比；seq 为模板碱基，供放大后
    渲染正向/互补碱基层。
    """
    total = max(1, seq_len)
    tracks = []
    for k, p in enumerate(pairs):
        amp_left = (p['amp_start'] - 1) / float(total) * 100.0
        amp_width = max(p['amp_end'] - p['amp_start'] + 1, 0) / float(total) * 100.0
        tracks.append({
            'idx': p['idx'],
            'type': p['type'],
            'top': 2 + 14 * k,
            'amp_left': round(amp_left, 4),
            'amp_width': round(amp_width, 4),
            'f_left': round((p['f_start'] - 1) / float(total) * 100.0, 4),
            'f_width': 0.29,
            'r_left': round((p['r_start'] - 1) / float(total) * 100.0, 4),
            'r_width': 0.29,
            'p_left': round((p['probe_start'] - 1) / float(total) * 100.0, 4) if p['probe_seq'] else 0.0,
            'p_width': round(p['probe_len'] / float(total) * 100.0, 4) if p['probe_seq'] else 0.0,
            'has_probe': bool(p['probe_seq']),
            # bp 坐标（1-based；R 引物结合区为 [r_start, r_start+r_len-1]，
            # 其 3' 端在左侧 —— primer3 的 PRIMER_RIGHT 起点即正链最左碱基）
            'f_start': p['f_start'], 'f_len': p['f_len'],
            'r_start': p['r_start'], 'r_len': p['r_len'],
            'p_start': p['probe_start'] if p['probe_seq'] else 0,
            'p_len': p['probe_len'] if p['probe_seq'] else 0,
            'amp_start': p['amp_start'], 'amp_end': p['amp_end'],
            'product': p['product'],
            'amp_span': p.get('amplicon_span', p['product']),
            'score': p['score'],
            'recommendation': p['recommendation'],
            'f_seq': p['f_seq'],
            'r_seq': p['r_seq'],
            'probe_seq': p['probe_seq'],
        })
    return {'seq_len': seq_len, 'seq': (seq or '').upper(), 'tracks': tracks,
            'track_height': 14 * max(1, len(pairs)) + 4}


# ============================================================
# 主入口
# ============================================================
def run_design(fasta_text='', fasta_path='', ptype='PCR',
               core=None, adv=None, max_seqs=5, specificity=False, logger=None):
    """设计入口。

    fasta_text 与 fasta_path 二选一；fasta_path 优先。
    specificity=True 时对设计的扩增子跑宿主基因组 blastn 脱靶检查
    （primer.check_specificity_blast，需已配置 host-db；默认关）。
    返回 {'records': [...], 'tracks': [...], 'summary': {...}}
    """
    ptype = 'qPCR' if str(ptype).lower() == 'qpcr' else 'PCR'
    if fasta_path:
        seqs = _read_fasta_path(fasta_path)
    else:
        seqs = parse_fasta_text(fasta_text)
    if not seqs:
        raise ValueError('未解析到任何序列，请检查 FASTA 内容')

    records, tracks, n_pairs = [], [], 0
    for name, seq in seqs[:max_seqs]:
        try:
            r = design_for_sequence(name, seq, ptype, core, adv)
        except (OSError, RuntimeError) as e:
            # RuntimeError：primer3 不可用（_p3() 抛）——原只接 OSError，
            # 交互页会整体 500 而不是逐序列软跳过（2026-10-03 修复 P3-21）
            r = {'name': name, 'length': len(seq), 'pairs': [],
                 'explain': 'primer3 不可用或参数不合法: %s' % e,
                 'design_length': len(seq),
                 'region_offset': 0, 'p3_opts': {}}
        records.append({
            'name': r['name'], 'length': r['length'],
            'design_length': r.get('design_length', r['length']),
            'region_offset': r.get('region_offset', 0),
            'n_pairs': len(r['pairs']),
            'explain': r.get('explain', ''),
            'pairs': r['pairs'],
        })
        tr = build_track(r['length'], r['pairs'], seq)
        tr['name'] = r['name']
        tracks.append(tr)
        n_pairs += len(r['pairs'])

    # 特异性检查（可选，默认关；2026-10-03 修复 P2-14④ 接入交互路径）
    spec_checked = False
    spec_note = '未做脱靶特异性检查：未勾选特异性检查'
    if specificity:
        amp_fa = None
        try:
            import tempfile
            from .config import current_host_genome
            from .primer import check_specificity_blast
            host = current_host_genome()
            if not host or not os.path.isfile(host):
                # 宿主库缺失时 check_specificity_blast 会"跳过并返回空"，
                # 与"检查了但无命中"不可区分 → 先自行判定，如实标注未检查
                spec_note = ('未做脱靶特异性检查：未配置宿主基因组（host-db）'
                             '或未构建 BLAST 库')
                if logger:
                    logger.log(spec_note, 'WARN')
            else:
                fd, amp_fa = tempfile.mkstemp(prefix='vp_pd_amp_', suffix='.fa')
                os.close(fd)
                n_amp = 0
                amp_names = {}      # 扩增子 ID → 引物对 dict（回填 host_hit）
                with io.open(amp_fa, 'w', encoding='utf-8', newline='\n') as f:
                    seqmap = {nm: (sq or '').upper() for nm, sq in seqs[:max_seqs]}
                    for rec in records:
                        up = seqmap.get(rec['name'], '')
                        for p in rec['pairs']:
                            a0 = int(p.get('amp_start') or 0)
                            a1 = int(p.get('amp_end') or 0)
                            amp = up[a0 - 1:a1] if 0 < a0 <= a1 <= len(up) else ''
                            if not amp:
                                continue
                            pid = '%s_P%s' % (rec['name'],
                                              p.get('idx') or (n_amp + 1))
                            f.write('>%s\n%s\n' % (pid, amp))
                            amp_names[pid] = p
                            n_amp += 1
                if n_amp:
                    hits = check_specificity_blast(amp_fa, logger=logger) or {}
                    spec_checked = True
                    spec_note = (f'宿主基因组 blastn 脱靶检查：{len(hits)}/{n_amp} '
                                 f'对扩增子命中宿主')
                    for pid, p in amp_names.items():
                        p['host_hit'] = 'YES' if pid in hits else 'NO'
                else:
                    spec_note = '未做脱靶特异性检查：没有可用的扩增子序列'
        except Exception as e:
            spec_note = f'脱靶特异性检查失败（{type(e).__name__}）：未检查'
            if logger:
                logger.log(spec_note + f': {e}', 'WARN')
        finally:
            if amp_fa and os.path.isfile(amp_fa):
                try:
                    os.remove(amp_fa)
                except OSError:
                    pass

    scores = [p['score'] for rec in records for p in rec['pairs']]
    summary = {
        'type': ptype,
        'n_sequences': len(records),
        'n_pairs': n_pairs,
        'best_score': round(max(scores), 1) if scores else 0.0,
        'avg_score': round(sum(scores) / len(scores), 1) if scores else 0.0,
        'n_recommended': sum(1 for s in scores if s >= pt.RECOMMENDED_MIN),
        'n_usable': sum(1 for s in scores if pt.USABLE_MIN <= s < pt.RECOMMENDED_MIN),
        'base_score': pt.QPCR_BASE_SCORE if ptype == 'qPCR' else pt.PCR_BASE_SCORE,
        'probe_present': any(p['probe_seq'] for rec in records for p in rec['pairs']),
        # 恒写特异性声明（true=检查真的跑成；false=未检查，不代表无脱靶）
        'specificity_requested': bool(specificity),
        'specificity_checked': bool(spec_checked),
        'specificity_note': spec_note,
    }
    return {'records': records, 'tracks': tracks, 'summary': summary}
