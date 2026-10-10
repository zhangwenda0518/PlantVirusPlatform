# -*- coding: utf-8 -*-
"""
引物热力学评估与评分（独立模块）

设计原则
--------
1. 不评估 BLAST 脱靶（平台无宿主库依赖，且耗时不可控）
2. 不评估覆盖度（单序列输入时无多序列参考，恒为满分属空转项）
3. 评分归一化到 100：score = (基准分 - 实际扣分) / 基准分 * 100
   - PCR  基准 45 = 二聚体 30 + GC/Tm 15
   - qPCR 基准 75 = 二聚体 30 + GC/Tm 15 + 探针专项 30
4. 无探针的 qPCR 按"探针不合格"扣满 30 分，不降级到 PCR 基准
   （qPCR 无探针属方法学缺陷，降级评分会把缺陷藏起来）

扣分档位沿用 step3_validate_primers.py 口径，便于与既有结果对拍。
"""

# primer3 统一走 primer3_runtime 的子进程探测 + 惰性导入：
# 模块级 import 在内存不足的机器上会让 pyd 初始化段错误杀掉整个进程。
from . import primer3_runtime as _p3rt

try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.SeqUtils import gc_fraction
    BIO_OK = True
except ImportError:
    BIO_OK = False


# ---------------- 基准分（归一化分母） ----------------
PCR_BASE_SCORE = 45.0     # 二聚体 30 + GC/Tm 15
QPCR_BASE_SCORE = 75.0    # 二聚体 30 + GC/Tm 15 + 探针专项 30

# ---------------- 扣分上限（用于文档与校验） ----------------
MAX_DIMER_PENALTY = 30.0
MAX_GCTM_PENALTY = 15.0
MAX_PROBE_PENALTY = 30.0

# ---------------- 推荐等级阈值（归一化后） ----------------
RECOMMENDED_MIN = 80.0
USABLE_MIN = 60.0


def _safe_call(func_name, *args, **kwargs):
    """primer3 个别序列可能抛异常（如超长同聚物），统一兜底为 None。

    func_name 是 primer3 模块上的函数名：primer3 不可用（探测失败）
    时返回 None 降级，不抛异常、更不段错误。
    """
    try:
        mod = _p3rt.get_module()
        return getattr(mod, func_name)(*args, **kwargs)
    except RuntimeError:
        return None
    except Exception:
        return None


# ================================================================
# 1. 二聚体分析
# ================================================================
def dimer_analysis(fwd_seq, rev_seq):
    """自二聚体 + 交叉二聚体。

    返回 ΔG (kcal/mol, 负值越强越差) 与 Tm (℃)。
    primer3 的 dg 单位是 cal/mol，除以 1000 转 kcal/mol。
    """
    out = {'self_fwd_dg': 0.0, 'self_rev_dg': 0.0, 'cross_dg': 0.0,
           'self_fwd_tm': 0.0, 'self_rev_tm': 0.0, 'cross_tm': 0.0}

    def _one(func, *args):
        r = _safe_call(func, *args)
        if r is None:
            return 0.0, 0.0
        return round(r.dg / 1000.0, 2), round(r.tm, 1)

    if fwd_seq:
        out['self_fwd_dg'], out['self_fwd_tm'] = _one('calc_homodimer', fwd_seq)
    if rev_seq:
        out['self_rev_dg'], out['self_rev_tm'] = _one('calc_homodimer', rev_seq)
    if fwd_seq and rev_seq:
        out['cross_dg'], out['cross_tm'] = _one('calc_heterodimer', fwd_seq, rev_seq)
    return out


# ================================================================
# 2. 发夹与 3' 端稳定性
# ================================================================
def hairpin_analysis(fwd_seq, rev_seq):
    """发夹结构 Tm（来自 primer3 设计输出，此处用于独立复核）。"""
    out = {'self_fwd_hairpin_tm': 0.0, 'self_rev_hairpin_tm': 0.0}
    if fwd_seq:
        r = _safe_call('calc_hairpin', fwd_seq)
        if r is not None:
            out['self_fwd_hairpin_tm'] = round(r.tm, 1)
    if rev_seq:
        r = _safe_call('calc_hairpin', rev_seq)
        if r is not None:
            out['self_rev_hairpin_tm'] = round(r.tm, 1)
    return out


def end_stability_analysis(fwd_seq, rev_seq):
    """3' 端相邻稳定性（kcal/mol）。

    用 calc_end_stability(seq1, seq2) 计算「seq1 的 3' 端退火到 seq2」的 ΔG。
    这里算的是 **引物对之间** 的 3' 端互配（Fwd 3' 端 vs Rev 原序列）：
    calc_end_stability 把 seq2 按自身 5'→3' 读入，F 的 3' 端与之反平行互补
    即真实互配。原先传 revcomp(rev) 等于把 R 摆成模板链方向，真实 3' 互配
    被算成不结合（2026-10-03 修复 P2-10；实测 F=ACGTA / R=TACGT：
    传 R 得 -2749 cal/mol，传 revcomp(R) 得 +1700 cal/mol）。
    单条引物的 3' 端稳定性请取 primer3 设计输出的
    PRIMER_LEFT_i_END_STABILITY 字段（设计阶段已给出，无需重算）。
    """
    out = {'fwd_rev_end_dg': 0.0}
    if not (fwd_seq and rev_seq):
        return out
    r = _safe_call('calc_end_stability', fwd_seq, rev_seq)
    if r is not None:
        out['fwd_rev_end_dg'] = round(r.dg / 1000.0, 2)
    return out


# ================================================================
# 3. GC / Tm 复核（Biopython Nearest-Neighbor）
# ================================================================
def gc_tm_check(fwd_seq, rev_seq, probe_seq=''):
    """用 Biopython Tm_NN 独立复核 Tm，GC 用 gc_fraction 计算。

    口径旁注（2026-10-03 修复 P3-19，取值与字段名保持不变）：
      * tm_note：Tm_NN 用 Biopython 默认盐条件（Na 50 mM / Mg 1.5 mM /
        dNTP 0.6 mM / DNA 50 nM），与 primer3 内部条件不同，仅供排序与
        交叉核对，不可与 primer3 的 Tm 混用比较。
      * gc_note：gc_fraction 默认 ambiguous='remove'——含 N 等简并碱基的
        引物只统计非 N 子集（分母同时剔除），序列偏短时 GC% 波动更大。
    """
    out = {'gc_fwd': 0.0, 'gc_rev': 0.0, 'gc_probe': 0.0,
           'tm_fwd_nn': 0.0, 'tm_rev_nn': 0.0, 'tm_probe_nn': 0.0,
           'tm_note': 'Biopython Tm_NN 默认盐条件（50 mM Na, 1.5 mM Mg, '
                      '0.6 mM dNTP, 50 nM DNA），与 primer3 内部条件不同，仅供排序',
           'gc_note': "gc_fraction(ambiguous='remove')：含 N 等简并碱基时只统计非 N 子集"}
    if not BIO_OK:
        return out

    def _gc(s):
        try:
            return round(gc_fraction(Seq(s)) * 100, 1)
        except Exception:
            return 0.0

    def _tm(s):
        try:
            return round(mt.Tm_NN(Seq(s)), 1)
        except Exception:
            return 0.0

    if fwd_seq:
        out['gc_fwd'], out['tm_fwd_nn'] = _gc(fwd_seq), _tm(fwd_seq)
    if rev_seq:
        out['gc_rev'], out['tm_rev_nn'] = _gc(rev_seq), _tm(rev_seq)
    if probe_seq:
        out['gc_probe'], out['tm_probe_nn'] = _gc(probe_seq), _tm(probe_seq)
    return out


# ================================================================
# 4. 探针专项验证
# ================================================================
def probe_validation(probe_seq, fwd_seq, rev_seq,
                     probe_tm=0.0, fwd_tm=0.0, rev_tm=0.0):
    """探针四项质量检查（varVAMP / qprimer_designer 口径）。

    1. 发夹 Tm ≤ 47 ℃
    2. 自二聚体 ΔG ≥ -6 kcal/mol
    3. 与 Fwd / Rev 交叉二聚体 ΔG ≥ -6 kcal/mol
    4. 探针 Tm 比引物 Tm 高 5-10 ℃
    """
    out = {'probe_hairpin_tm': 0.0, 'probe_self_dg': 0.0,
           'probe_fwd_dg': 0.0, 'probe_rev_dg': 0.0,
           'probe_tm_diff_ok': True, 'probe_warnings': []}
    if not probe_seq or len(probe_seq) < 14:
        return out

    warns = []
    r = _safe_call('calc_hairpin', probe_seq)
    if r is not None:
        out['probe_hairpin_tm'] = round(r.tm, 1)
        if r.tm > 47:
            warns.append('probe hairpin Tm %.0f C (>47)' % r.tm)

    r = _safe_call('calc_homodimer', probe_seq)
    if r is not None:
        out['probe_self_dg'] = round(r.dg / 1000.0, 2)
        if r.dg / 1000.0 < -6.0:
            warns.append('probe self-dimer dG %.1f' % (r.dg / 1000.0))

    if fwd_seq:
        r = _safe_call('calc_heterodimer', probe_seq, fwd_seq)
        if r is not None:
            out['probe_fwd_dg'] = round(r.dg / 1000.0, 2)
            if r.dg / 1000.0 < -6.0:
                warns.append('probe-Fwd dimer dG %.1f' % (r.dg / 1000.0))
    if rev_seq:
        r = _safe_call('calc_heterodimer', probe_seq, rev_seq)
        if r is not None:
            out['probe_rev_dg'] = round(r.dg / 1000.0, 2)
            if r.dg / 1000.0 < -6.0:
                warns.append('probe-Rev dimer dG %.1f' % (r.dg / 1000.0))

    if probe_tm > 0 and fwd_tm > 0 and rev_tm > 0:
        ref = max(fwd_tm, rev_tm)
        if not (ref + 5 <= probe_tm <= ref + 10):
            out['probe_tm_diff_ok'] = False
            warns.append('probe Tm %.1f C not 5-10 C above primers' % probe_tm)

    out['probe_warnings'] = warns
    return out


# ================================================================
# 5. 评分（归一化到 100）
# ================================================================
def compute_score(ptype, dimer, gctm, probe=None, has_probe=False):
    """按扣分档位计分，再归一化到 100。

    ptype    : 'PCR' | 'qPCR'
    dimer    : dimer_analysis() 返回值
    gctm     : {gc_fwd, gc_rev, tm_fwd_nn, tm_rev_nn}（用 NN 复核值）
    probe    : probe_validation() 返回值
    has_probe: 该引物对是否实际带探针

    返回 (score_0_100, breakdown_dict)
    """
    base = QPCR_BASE_SCORE if ptype == 'qPCR' else PCR_BASE_SCORE
    deduct = 0.0
    detail = []

    # --- 二聚体（最多 -30）---
    d_dim = 0.0
    for key, label in (('self_fwd_dg', 'self_fwd'),
                       ('self_rev_dg', 'self_rev'),
                       ('cross_dg', 'cross')):
        dg = dimer.get(key, 0.0) or 0.0
        if dg < -9.0:
            d_dim += 12
            detail.append('%s_dg %.1f (-12)' % (label, dg))
        elif dg < -6.0:
            d_dim += 6
            detail.append('%s_dg %.1f (-6)' % (label, dg))
        elif dg < -3.0:
            d_dim += 2
            detail.append('%s_dg %.1f (-2)' % (label, dg))
    d_dim = min(d_dim, MAX_DIMER_PENALTY)
    deduct += d_dim

    # --- GC 偏差 + Tm 差（最多 -15）---
    d_gc = 0.0
    for key in ('gc_fwd', 'gc_rev'):
        gc = gctm.get(key)
        # 用 is None 判定"算不出"，不能用 `or 50.0`：GC=0%（poly-A/T 引物）
        # 会被当成缺失值 → 直接落进无偏差区间，扣分永不生效。
        if gc is None:
            continue
        if gc < 35 or gc > 65:
            d_gc += 5
        elif gc < 40 or gc > 60:
            d_gc += 2
    tf = gctm.get('tm_fwd_nn', 0.0) or 0.0
    tr = gctm.get('tm_rev_nn', 0.0) or 0.0
    if tf > 0 and tr > 0 and abs(tf - tr) > 3:
        d_gc += 3
        detail.append('Tm diff %.1f (-3)' % abs(tf - tr))
    d_gc = min(d_gc, MAX_GCTM_PENALTY)
    deduct += d_gc

    # --- 探针专项（仅 qPCR，最多 -30）---
    d_probe = 0.0
    if ptype == 'qPCR':
        if not has_probe:
            d_probe = MAX_PROBE_PENALTY
            detail.append('no probe (-30)')
        elif probe:
            if not probe.get('probe_tm_diff_ok', True):
                d_probe += 20
                detail.append('probe Tm off (-20)')
            if probe.get('probe_hairpin_tm', 0) > 47:
                d_probe += 8
                detail.append('probe hairpin (-8)')
            if probe.get('probe_self_dg', 0) < -6.0:
                d_probe += 8
                detail.append('probe self-dimer (-8)')
            for k, lbl in (('probe_fwd_dg', 'probe-Fwd'), ('probe_rev_dg', 'probe-Rev')):
                if probe.get(k, 0) < -6.0:
                    d_probe += 6
                    detail.append('%s dimer (-6)' % lbl)
        d_probe = min(d_probe, MAX_PROBE_PENALTY)
    deduct += d_probe

    raw = max(0.0, base - deduct)
    score = round(raw / base * 100.0, 1) if base > 0 else 0.0
    return score, {
        'base': base,
        'deduct_total': round(deduct, 1),
        'deduct_dimer': round(d_dim, 1),
        'deduct_gctm': round(d_gc, 1),
        'deduct_probe': round(d_probe, 1),
        'raw_before_norm': round(raw, 1),
        'detail': detail,
    }


def recommendation(score):
    if score >= RECOMMENDED_MIN:
        return 'RECOMMENDED'
    if score >= USABLE_MIN:
        return 'USABLE'
    return 'NOT_RECOMMENDED'


# ================================================================
# 6. 一次性评估一个引物对
# ================================================================
def evaluate_pair(fwd_seq, rev_seq, probe_seq='',
                  ptype='PCR', p3_fwd_tm=0.0, p3_rev_tm=0.0, p3_probe_tm=0.0,
                  end_fwd=0.0, end_rev=0.0):
    """整合所有热力学项并给出归一化评分。

    p3_*_tm 为 primer3 设计时给出的 Tm，用于探针 Tm 差的判断；
    郭缺时退化为 Biopython NN 的 Tm。
    end_fwd/end_rev 为 primer3 设计输出的 3' 端稳定性（可选）。
    """
    dimer = dimer_analysis(fwd_seq, rev_seq)
    hairpin = hairpin_analysis(fwd_seq, rev_seq)
    ends = end_stability_analysis(fwd_seq, rev_seq)
    gctm = gc_tm_check(fwd_seq, rev_seq, probe_seq)

    ends['fwd_end_stability'] = round(float(end_fwd or 0.0), 2)
    ends['rev_end_stability'] = round(float(end_rev or 0.0), 2)

    eff_probe_tm = p3_probe_tm or gctm.get('tm_probe_nn', 0.0)
    eff_fwd_tm = p3_fwd_tm or gctm.get('tm_fwd_nn', 0.0)
    eff_rev_tm = p3_rev_tm or gctm.get('tm_rev_nn', 0.0)

    probe = probe_validation(probe_seq, fwd_seq, rev_seq,
                             eff_probe_tm, eff_fwd_tm, eff_rev_tm)
    has_probe = bool(probe_seq and len(probe_seq) >= 14)

    score, breakdown = compute_score(ptype, dimer, gctm, probe, has_probe)

    out = {}
    out.update(dimer)
    out.update(hairpin)
    out.update(ends)
    out.update(gctm)
    out.update(probe)
    out['probe_present'] = has_probe
    out['score'] = score
    out['recommendation'] = recommendation(score)
    out['score_breakdown'] = breakdown
    return out
