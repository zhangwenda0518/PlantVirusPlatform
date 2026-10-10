# -*- coding: utf-8 -*-
"""dsRNA 臂/siRNA 效价与安全评分（dsRIP 论文公式的平台统一实现）。

合并了 dsRNAmax 管线中 dsrip_scoring.py 与 safety_scoring.py 两份复制实现
（Cedden et al., BMC Biology 2025, 23:114 的忠实移植），作为平台唯一评分口径：
- 效价：asymScore（lookup.db 四核苷酸热力学不对称查表，权重 5）
        + self_fold_energy（ViennaRNA MFE，权重 3）+ 各 bonus
- 安全：逐 siRNA 脱靶加权计数（致死基因 ×safety_essential_coeff，默认 20）
        臂内相对归一化 + 与效价按 eff/safety 优先级（默认 50/50）合成，
        滑窗（250–350 siRNA）取 combined 最优窗口
- 绝对红线：致死基因 0 错配命中一票否决（由调用方依据扫描结果判定）

与 C-host 管线版本的既知差异（均为有意修正，勿"对齐"回去）：
1. find_best_window / pick_best_window 的 nt_end = 末位 siRNA 起点 + 20
   （C-host 侧 dsrip_scoring.py 历史版本短了 20nt，safety_scoring.py 早已修；
   dsrip_scoring.py 也已于 2026-09-17 补齐，并经 2026-09-18 数值对拍确认
   与本模块逐 siRNA total_score 一致、w933/w3806/w11863 最优窗坐标完全相同）
2. edge_asym 用修正规则（antisense[0]∈{A,U} 且 antisense[18]∈{G,C}），
   上游 dsRIP 第二子句写错成 antisense[0]=='C'，管线文档已留档
3. 脱靶来源解耦：本模块不做任何 bowtie/扫描，只消费调用方传入的
   {siRNA序号: 加权命中数}（平台用 dsrna_offtarget 的 uint64 扫描器产出）

依赖：viennarna（import RNA）、numpy；数据文件 databases/dsrna/{lookup.db,
siRNA_parameters.txt}。
"""
import sqlite3
from pathlib import Path

import numpy as np

try:
    import RNA
except Exception as _rna_err:      # noqa: BLE001
    # RNA/__init__py 里 `except Exception: raise ImportError(...)` 把底层
    # 错误（通常是冻结分发下缺某个 DLL）掩盖成一句固定文案，无从排查。
    # 这里绕开包 __init__ 直接加载底层扩展，把**真实原因**带进报错。
    import glob as _glob
    import importlib.machinery as _mach
    import importlib.util as _util
    import os as _os

    # find_spec('RNA') 只定位不执行，安全；包目录取其 submodule_search_locations
    _spec = _util.find_spec('RNA')
    _pkg_dir = ''
    if _spec is not None and _spec.submodule_search_locations:
        _pkg_dir = list(_spec.submodule_search_locations)[0]
    _detail = f'{type(_rna_err).__name__}: {_rna_err}'
    if _pkg_dir:
        for _pyd in sorted(_glob.glob(_os.path.join(_pkg_dir, '_RNA*.pyd'))):
            try:
                _loader = _mach.ExtensionFileLoader('_RNA_diag', _pyd)
                _mod = _util.module_from_spec(_util.spec_from_loader(
                    '_RNA_diag', _loader))
                _loader.exec_module(_mod)
            except Exception as _e2:            # noqa: BLE001
                _detail = (f'底层扩展 {_os.path.basename(_pyd)} 加载失败 — '
                           f'{type(_e2).__name__}: {_e2}')
                break
    raise ImportError(f'ViennaRNA(RNA) 加载失败：{_detail}') from _rna_err

from .config import PLATFORM_ROOT

DSRNA_DATA_DIR = Path(PLATFORM_ROOT) / 'databases' / 'dsrna'
LOOKUP_DB = DSRNA_DATA_DIR / 'lookup.db'
PARAMS_PATH = DSRNA_DATA_DIR / 'siRNA_parameters.txt'

# 21nt siRNA（Dicer 产物）；窗口以 siRNA 条数计
SIRNA_LEN = 21

DEFAULT_PARAMS = {
    'asymScore_weight': 5,
    'self_fold_energy_weight': 3,
    'edge_asyym_bonus': 10,
    'anti_GC_bonus': 0,
    'anti_GC_upper_bound': 65,
    'anti_10th_A_bonus': 10,
    'anti_9_14_GC_bonus': 5,
    'anti_9_14_GC_lower_bound': 35,
    'GGGG_CCCC_penalty': 0,
    'safety_essential_coeff': 20,
    'min_window': 250,
    'max_window': 350,
}

# T7 class III 启动子（+ 转录起始 G）；正反引物都带 → PCR 产物两端各一个
# 启动子，IVT 双向转录自发退火成 dsRNA
T7_PROMOTER = 'TAATACGACTCACTATAGGG'


class ScoringDataError(RuntimeError):
    """lookup.db / 参数文件缺失等数据面错误。"""


def load_params(path=None):
    """读 siRNA_parameters.txt（存在且键合法的 int 覆盖默认值）。"""
    p = dict(DEFAULT_PARAMS)
    path = Path(path) if path else PARAMS_PATH
    if path.exists():
        for line in path.read_text(encoding='utf-8', errors='ignore').splitlines():
            line = line.strip()
            if not line or '=' not in line:
                continue
            k, v = line.split('=', 1)
            if k.strip() in p:
                try:
                    p[k.strip()] = int(v)
                except ValueError:
                    pass
    return p


# ---------------------------------------------------------------- 序列工具

_COMP = str.maketrans('ACGTUacgtu', 'TGCAAtgcaa')


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def comp_rna(seq):
    """RNA 互补，保持 5'→3'（对齐 dsRIP Comp_RNA）。"""
    return seq.translate(str.maketrans('ACGU', 'UGCA'))


def read_fasta_text(text):
    """FASTA 文本 -> [(name, seq)]；序列行剔除非字母并转大写。"""
    seqs, name, buf = [], None, []
    import re
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                seqs.append((name, ''.join(buf).upper()))
            name = line[1:].strip().split()[0] or f'seq{len(seqs) + 1}'
            buf = []
        else:
            buf.append(re.sub(r'[^A-Za-z]', '', line).upper())
    if name is not None:
        seqs.append((name, ''.join(buf).upper()))
    return [(n, s) for n, s in seqs if s]


def read_fasta(path):
    return read_fasta_text(Path(path).read_text(encoding='utf-8', errors='ignore'))


# ---------------------------------------------------------------- dsRIP 核心

def generate_siRNA(sequence, sirna_len=SIRNA_LEN):
    """dsRIP generate_siRNA 忠实移植。

    反义链 5'→3'，两端补虚拟 'XX' 模拟 Dicer 2nt 突出端邻居上下文，再切
    (k+2)-mer。返回 [(antisense_5_3, sense_3_5)]（基因组顺序），
    窗口数 = L - k + 1（全部 21-mer）。
    """
    whole_antisense_5_3 = revcomp(sequence).replace('T', 'U')
    whole_antisense_XX = 'XX' + whole_antisense_5_3
    whole_sense_XX = 'XX' + comp_rna(whole_antisense_5_3)
    k = sirna_len + 2
    pairs = []
    for i in range(len(whole_antisense_XX) - k + 1):
        a = whole_antisense_XX[i:i + k][2:]
        s = whole_sense_XX[i:i + k][:-2]
        pairs.append((a, s))
    return pairs


class ThermoLookup:
    """dsRIP lookup.db：256 条 RNA 四核苷酸最近邻 dG 值。"""

    def __init__(self, db_path=None):
        db_path = Path(db_path) if db_path else LOOKUP_DB
        if not db_path.exists():
            raise ScoringDataError(
                f'缺少热力学查找表 {db_path}（应随发行包提供 databases/dsrna/lookup.db）')
        conn = sqlite3.connect(str(db_path))
        try:
            self.table = dict(conn.execute('SELECT key, value FROM lookup'))
        finally:
            conn.close()

    def get(self, key):
        return self.table.get(key)


def thermo_asym(siRNA_tuple, lookup):
    """asymScore = dG(反义链 5' 端 4-mer) - dG(sense 3' 端 4-mer 反向)。

    正值越大 = 反义链 5' 端越弱配对 → RISC 越倾向选它当 guide。
    """
    antisense, sense = siRNA_tuple
    anti_val = lookup.get(antisense[:4])
    sense_val = lookup.get(sense[-4:][::-1])
    if anti_val is None or sense_val is None:
        raise ScoringDataError(
            f'lookup.db 缺少四核苷酸: {antisense[:4]} / {sense[-4:][::-1]}')
    return anti_val - sense_val


def _find_substring(s):
    for sub in ('GGGG', 'CCCC'):
        if sub in s:
            return sub
    return None


def siRNA_feature_prediction(siRNAs, sirna_len=SIRNA_LEN, temp=25, lookup=None):
    """逐 siRNA 双链特征提取（不含 mRNA 可及性/ORF——病毒场景无意义，
    上游权重本就是 0）。"""
    lookup = lookup or ThermoLookup()
    md = RNA.md()
    md.temperature = temp
    feats = {}
    for i, (antisense, sense) in enumerate(siRNAs):
        anti_9_14 = antisense[8:14]
        gc_9_14 = round((anti_9_14.count('C') + anti_9_14.count('G')) / 6 * 100, 1)
        anti_gc = round((antisense.count('C') + antisense.count('G'))
                        / len(antisense) * 100, 1)
        asym = round(thermo_asym((antisense, sense), lookup), 3)
        fc = RNA.fold_compound(antisense, md)
        structure, mfe = fc.mfe()
        # 修正版边缘不对称（上游 dsRIP 第二子句笔误，见模块 docstring）
        edge_asym = (1 if antisense[0] in ('A', 'U') else 0) \
            + (1 if antisense[18] in ('G', 'C') else 0)
        feats[f'sirna{i + 1}_{i + sirna_len}'] = {
            'antisense_5_3': antisense,
            'sense_3_5': sense,
            'anti_GC': anti_gc,
            'anti_9_14_GC': gc_9_14,
            'asymScore': asym,
            'self_fold_energy': round(mfe, 4),
            'anti_10th_A': 1 if antisense[9] == 'A' else 0,
            'G_C_repeat': _find_substring(antisense),
            'edge_asym': edge_asym,
        }
    return feats


def normalize(value, lo, hi, tlo=0, thi=100):
    """dsRIP normalize_score：把 [lo, hi] 线性映射到 [tlo, thi]。

    调用方**故意传反向区间**——`(7.5, -7.5)` 与 `(0, -15)`——让"弱 5' 端 /
    弱折叠"得高分，所以反向区间给出反向输出是设计意图，不是 bug。

    2026-09-18 与本仓 `dsrip_scoring.normalize` 逐位对齐：改写前本模块用
    手写等价式 `(v + 7.5) / 15.0 * 100.0`，数学恒等但 IEEE-754 运算顺序不同，
    840 条 siRNA 里有 14 条 `total_score` 差 0.1 分（末位差落在 round(,1) 的
    .x5 边界）。改用同一写法后两侧**逐位相同**（见 `_crosscheck_platform_20260918.py`）。
    """
    return (value - lo) / (hi - lo) * (thi - tlo) + tlo


def score_efficiency(features, params):
    """效价总分：归一化 asym（±7.5→0/100，弱 5' 端高分）与 MFE（0/-15→100/0，
    弱折叠高分）加权平均，再加 bonus。"""
    a_w = params['asymScore_weight']
    f_w = params['self_fold_energy_weight']
    for feats in features.values():
        asym_n = normalize(feats['asymScore'], 7.5, -7.5, 100, 0)
        fold_n = normalize(feats['self_fold_energy'], 0, -15, 100, 0)
        base = (asym_n * a_w + fold_n * f_w) / (a_w + f_w)
        extra = 0
        if feats['G_C_repeat']:
            extra -= params['GGGG_CCCC_penalty']
        if feats['edge_asym']:
            extra += params['edge_asyym_bonus']
        if feats['anti_GC'] <= params['anti_GC_upper_bound']:
            extra += params['anti_GC_bonus']
        if feats['anti_9_14_GC'] >= params['anti_9_14_GC_lower_bound']:
            extra += params['anti_9_14_GC_bonus']
        if feats['anti_10th_A']:
            extra += params['anti_10th_A_bonus']
        feats['total_score'] = round(base + extra, 1)
    return features


def _sorted_names(features):
    return sorted(features, key=lambda n: int(n.split('_')[0][5:]))


def assign_safety_scores(features, off_counts):
    """把逐 siRNA 加权脱靶计数挂到特征上并做臂内相对归一化。

    off_counts: {siRNA序号(1-based): 加权命中数}（未出现在 dict 里的 = 0）。
    加权（致死 ×coeff）由调用方完成——本函数只消费已加权计数。
    归一化：臂内最安全 = 100，最差 = 0（dsRIP 口径）。注意这是**相对**分，
    绝对风险由调用方另行判定（红线/阈值），本模块在 report 里保留 raw 计数。
    """
    names = _sorted_names(features)
    raw = {n: int(off_counts.get(i, 0)) for i, n in enumerate(names, 1)}
    vals = list(raw.values())
    vmin, vmax = (min(vals), max(vals)) if vals else (0, 0)
    rng = vmax - vmin
    for n in names:
        feats = features[n]
        feats['raw_off_target_count'] = raw[n]
        feats['norm_off_target'] = round((1 - (raw[n] - vmin) / rng) * 100, 1) \
            if rng > 0 else 100.0
    return features


def pick_best_window(features, efficiency_priority=50, safety_priority=50,
                     min_window=250, max_window=350):
    """combined = (eff×效价 + safety×归一化安全)/100，滑窗取均值最大窗口。

    返回窗口的 nt_start..nt_end 为**基因组 1-based 闭区间**（nt_end =
    末位 siRNA 起点 + 20，即末位 siRNA 的 3' 端）。臂短于 min_window 时
    回退整臂打分并在 note 标注（不返回非有限值）。
    """
    names = _sorted_names(features)
    if not names:
        return {'mean_combined': None, 'start_siRNA': None, 'end_siRNA': None,
                'nt_start': None, 'nt_end': None, 'note': 'no siRNA'}
    total = np.array([features[n]['total_score'] for n in names], dtype=float)
    safety = np.array([features[n].get('norm_off_target', 0.0) for n in names],
                      dtype=float)
    combined = (efficiency_priority * total + safety_priority * safety) / 100.0
    n = len(combined)

    def _window(best, scores):
        if n < min_window:
            return {'mean': float(scores.mean()), 's': 0, 'e': n,
                    'note': f'arm shorter than min_window ({n} < {min_window}); '
                            f'whole arm scored'}
        best_v, best_s, best_e = -np.inf, 0, 0
        csum = np.concatenate([[0.0], np.cumsum(scores)])
        for w in range(min_window, min(max_window, n) + 1):
            means = (csum[w:] - csum[:-w]) / w
            i = int(np.argmax(means))
            if means[i] > best_v:
                best_v, best_s, best_e = float(means[i]), i, i + w
        return {'mean': best_v, 's': best_s, 'e': best_e, 'note': ''}

    weff = _window(None, total)
    wsafe = _window(None, combined)
    win = wsafe                                   # 选窗以 combined 为准
    s, e = win['s'], win['e']
    # 臂（sense）帧坐标（2026-10-03 增加）：_sorted_names 沿**反义链（guide）**
    # 5'→3' 排列（序号 1 = guide 5' 端窗口 = 臂 3' 末端窗口），故 nt_start/
    # nt_end 是 guide 帧坐标（dsRIP/C-host safety_scoring 同口径，逐位一致，
    # 勿改）。窗口覆盖 guide 位置 [s, e+19]（0-based），换算到臂帧即
    # [L-1-(e+19), L-1-s]（L = 臂长 = n + SIRNA_LEN - 1）→ 1-based 闭区间。
    # 引物切片/覆盖度判定必须用 sense_nt_*：直接 arm_seq[nt_start-1:nt_end]
    # 取到的是镜像区间（臂长 > 窗口+20 时两者错位，见 2026-10-03 审计）。
    _L = len(names) + SIRNA_LEN - 1
    return {
        'mean_combined': win['mean'],
        'start_siRNA': names[s],
        'end_siRNA': names[e - 1],
        'nt_start': s + 1,
        'nt_end': (e - 1) + SIRNA_LEN,
        'sense_nt_start': _L - ((e - 1) + SIRNA_LEN) + 1,
        'sense_nt_end': _L - (s + 1) + 1,
        'frame': 'nt_start/nt_end = guide(反义链) 帧；sense_nt_* = 臂(正义链) 帧',
        'window_efficiency_mean': float(total[s:e].mean()),
        'window_safety_mean': float(safety[s:e].mean()),
        'window_raw_offtarget_sum': int(sum(features[n].get('raw_off_target_count', 0)
                                            for n in names[s:e])),
        'note': win['note'],
        'efficacy_best_window': {                    # 纯效价窗（参考信息）
            'nt_start': weff['s'] + 1, 'nt_end': weff['e'] - 1 + SIRNA_LEN,
            'sense_nt_start': _L - (weff['e'] - 1 + SIRNA_LEN) + 1,
            'sense_nt_end': _L - (weff['s'] + 1) + 1,
            'mean_score': weff['mean'],
        },
    }


# ---------------------------------------------------------------- 驱动入口

def score_sequence(name, seq, params=None, off_counts=None,
                   efficiency_priority=50, safety_priority=50, top=5,
                   lookup=None):
    """单条臂完整评分：效价 + （可选）安全 → combined 最优窗口。

    off_counts: {siRNA序号: 加权命中数}；None/空 = 无脱靶数据（安全维度
    全部按 100 处理，窗口退化为纯效价选窗）。
    返回 JSON 可序列化 dict。
    """
    params = params or load_params()
    seq = (seq or '').upper()
    siRNAs = generate_siRNA(seq)
    feats = siRNA_feature_prediction(siRNAs, lookup=lookup)
    feats = score_efficiency(feats, params)
    # 无脱靶数据也走归一化（全 0 → 全 100），combined 量纲与有数据时一致
    feats = assign_safety_scores(feats, off_counts or {})
    best = pick_best_window(feats, efficiency_priority, safety_priority,
                            params['min_window'], params['max_window'])
    scores = [f['total_score'] for f in feats.values()]
    top_sirnas = sorted(feats.items(), key=lambda kv: -kv[1]['total_score'])[:top]
    return {
        'name': name,
        'length': len(seq),
        'n_siRNAs': len(feats),
        'median_score': round(float(np.median(scores)), 2) if scores else None,
        'mean_score': round(float(np.mean(scores)), 2) if scores else None,
        'best_window': best,
        'top_siRNAs': [
            {'name': n, 'score': f['total_score'],
             'antisense': f['antisense_5_3'], 'asymScore': f['asymScore'],
             'self_fold': f['self_fold_energy'], '9_14_GC': f['anti_9_14_GC'],
             '10th_A': f['anti_10th_A'], 'raw_off_target':
                 f.get('raw_off_target_count', 0)}
            for n, f in top_sirnas],
        'efficacy_only': not bool(off_counts),
    }


def score_fasta(path_or_seqs, params=None, top=5):
    """批口径：仅效价评分（无脱靶维度），兼容旧 dsrip_scoring.score_fasta 用法。"""
    params = params or load_params()
    seqs = (path_or_seqs if isinstance(path_or_seqs, list)
            else read_fasta(path_or_seqs))
    out = {}
    for name, seq in seqs:
        r = score_sequence(name, seq, params=params, off_counts=None, top=top)
        r.pop('efficacy_only', None)
        out[name] = r
    return out


def amplicon_qc(name, seq, params=None, lookup=None):
    """成品扩增子质检口径：整段打效价（扩增子通常 < min_window，整臂打分）。"""
    return score_sequence(name, seq, params=params, off_counts=None,
                          efficiency_priority=100, safety_priority=0,
                          top=5, lookup=lookup)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        rep = score_fasta(sys.argv[1])
        for n, r in rep.items():
            bw = r['best_window']
            print(f"=== {n} ({r['length']} nt, {r['n_siRNAs']} siRNAs) "
                  f"median {r['median_score']}")
            print(f"  best window nt {bw['nt_start']}..{bw['nt_end']} "
                  f"mean {bw['mean_combined'] and round(bw['mean_combined'], 1)} "
                  f"{bw.get('note', '')}")
    else:
        print('usage: python dsrna_scoring.py arms.fa')
