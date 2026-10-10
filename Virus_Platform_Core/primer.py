# -*- coding: utf-8 -*-
"""
引物设计（primer3）——管线/工具侧编排层。

- conserved 模式：基于 05_phylo 多比对找保守区，在一致序列上设计（跨近缘参考通用）
- plain 模式：对病毒 contigs 直接设计覆盖全长的引物对
- 可选：blastn 扩增物 vs 宿主基因组特异性检查

设计引擎统一在 `Virus_Platform_Core/primer_design.py`（新版 primer3 接口 + `primer_thermo`
热力学评分 + TaqPCR/qPCR 探针支持）。本模块只保留管线特有逻辑：保守区
判定、全长分窗、`primers.tsv` 落盘、特异性 BLAST。

历史：本模块原先自持一套 primer3 参数（`P3_GLOBAL`）并用已弃用的
`primer3.bindings.designPrimers`，不产出任何质量评分——同一条序列经
「阶段⑧/工具页」与「/primer 交互页」得到不可比结论。合并后仅保留管线
特有的产物区间 `PIPELINE_PRODUCT_RANGE`，其余参数一律取自
`primer_design.CORE_DEFAULTS`；实测在产物区间与 num_return 相同时，
两条路径给出的引物对完全一致。
"""
import os
import json
from collections import Counter

from .config import get_config, current_host_genome
from .utils import (check_path, safe_open, run_cmd, iter_fasta, is_step_done, mark_step_done)

# 管线口径的产物区间（常规 RT-PCR 扩增子长度）。
# 其余参数（引物长度/Tm/GC/盐浓度/自二聚体阈值）与 /primer 交互页共用
# primer_design.CORE_DEFAULTS 与 ADV_DEFAULTS，不再各自维护。
PIPELINE_PRODUCT_RANGE = (300, 1500)


def design_primers_for_seq(name, seq, num_return=3, product_range=None):
    """对单条序列设计引物对。返回引物 dict 列表。

    字段名与历史版本一致（F_seq/F_pos/F_len/F_tm/F_gc、R_*、product、penalty、
    self_any_max、hairpin_max），并新增 score / recommendation / probe_* 与
    扩增子坐标。序列短于默认产物区间时自动收缩区间，避免 primer3 直接报错。
    另有 `amp_seq`：从模板截取的真实扩增子序列（模板缺失/坐标越界时为空串）。
    """
    from . import primer_design as pd
    pr = product_range or ((60, max(len(seq), 80)) if len(seq) < 300
                           else PIPELINE_PRODUCT_RANGE)
    core = {'product_min': int(pr[0]), 'product_max': int(pr[1]),
            'num_return': int(num_return)}
    rec = pd.design_for_sequence(name, seq, ptype='PCR', core=core, adv={})
    tpl = (seq or '').upper()
    pairs = []
    for i, p in enumerate(rec.get('pairs') or []):
        th = p.get('thermo') or {}
        # 真实扩增子序列（供宿主特异性 BLAST 用）：amp_start/amp_end 为
        # 1-based、均落在模板上（amp_end 即 R 引物 5' 端所在的正链碱基）
        _a0, _a1 = int(p.get('amp_start', 0) or 0), int(p.get('amp_end', 0) or 0)
        amp_seq = tpl[_a0 - 1:_a1] if 0 < _a0 <= _a1 <= len(tpl) else ''
        pairs.append({
            'name': f"{name}_P{i + 1}",
            'F_seq': p.get('f_seq', ''),
            'F_pos': p.get('f_start', 0),
            'F_len': p.get('f_len', 0),
            'F_tm': round(p.get('f_tm', 0.0) or 0.0, 1),
            'F_gc': round(p.get('f_gc', 0.0) or 0.0, 1),
            'R_seq': p.get('r_seq', ''),
            'R_pos': p.get('r_start', 0),
            'R_len': p.get('r_len', 0),
            'R_tm': round(p.get('r_tm', 0.0) or 0.0, 1),
            'R_gc': round(p.get('r_gc', 0.0) or 0.0, 1),
            'product': p.get('product', 0),
            'penalty': round(p.get('penalty', 0.0) or 0.0, 3),
            # primer3 对"无自二聚体/无发夹"报 0.0（负 Tm 钳到 0）；沿用该
            # 口径，使 history 表格语义不变（primer_thermo 给的是未钳制值）
            'self_any_max': round(max(0.0, th.get('self_fwd_tm') or 0.0,
                                      th.get('self_rev_tm') or 0.0), 1),
            'hairpin_max': round(max(0.0, th.get('self_fwd_hairpin_tm') or 0.0,
                                     th.get('self_rev_hairpin_tm') or 0.0), 1),
            'score': p.get('score'),
            'recommendation': p.get('recommendation', ''),
            'probe_seq': p.get('probe_seq', ''),
            'probe_tm': p.get('probe_tm', 0.0),
            'probe_gc': p.get('probe_gc', 0.0),
            'amp_start': p.get('amp_start', 0),
            'amp_end': p.get('amp_end', 0),
            'amp_seq': amp_seq,
        })
    return pairs


def conserved_regions_and_consensus(aln_fasta, min_len=400, min_ident=0.9,
                                    min_cov=0.7):
    """从比对找保守区段。返回 [(consensus_seq, aln_start, aln_end)]。"""
    names, seqs = [], []
    for h, s in iter_fasta(aln_fasta):
        names.append(h.split()[0])
        seqs.append(s.upper())
    if not seqs:
        return []
    n = len(seqs)
    L = len(seqs[0])
    # 长度不齐 = 输入不是多序列比对。zip(*seqs) 只产出最短列，而下面按
    # L（首条长度）索引 ok_cols → IndexError；工具页 conserved 模式直接
    # 复制用户 FASTA 当 aln.fasta，没有等长校验，必须在这里明确报错。
    if any(len(s) != L for s in seqs):
        raise ValueError(
            'conserved 模式需要等长的多序列比对（当前序列长度不一致，'
            '请先用 MAFFT 等工具比对后再上传）')
    cons_chars, ok_cols = [], []
    for col in zip(*seqs):
        non_gap = [c for c in col if c in 'ACGTU']
        if len(non_gap) >= n * min_cov:
            most, cnt = Counter(non_gap).most_common(1)[0]
            frac = cnt / len(non_gap)
            ok_cols.append(frac >= min_ident)
            cons_chars.append(most if frac >= min_ident else 'N')
        else:
            ok_cols.append(False)
            cons_chars.append('-')
    # 连续 ok 区段
    regions = []
    start = None
    for i in range(L):
        if ok_cols[i] and start is None:
            start = i
        elif not ok_cols[i] and start is not None:
            regions.append((start, i))
            start = None
    if start is not None:
        regions.append((start, L))
    out = []
    for s, e in regions:
        if e - s >= min_len:
            sub = ''.join(cons_chars[s:e]).replace('-', '').replace('U', 'T')
            if len(sub) >= min_len:
                out.append((sub, s, e))
    return out


def check_specificity_blast(amplicons_fasta, host_genome=None, logger=None):
    """可选：扩增子 blastn vs 宿主基因组，返回命中 {引物名: [宿主 contig, ...]}。

    命中判据按查询长度自适应：比对长度 >= min(100, 查询长度)。真实扩增子
    （数百 bp）沿用 ≥100bp 的实质命中口径；占位序列（两引物串联，约 40bp）
    要求整段命中，避免"永远不可能命中"的死阈值把 host_hit 固定成 NO。
    """
    from .assembly import _ascii_work_base
    cfg = get_config()
    blastn = cfg.tool('blastn')
    makeblastdb = cfg.tool('makeblastdb')
    host = host_genome or current_host_genome()
    # current_host_genome() 无宿主库时返回 None，os.path.isfile(None) 会抛
    # TypeError 并被上层吞成"特异性检查失败: ...not NoneType"。
    if not host or not os.path.isfile(host):
        if logger:
            logger.log('未配置宿主基因组（host-db），跳过引物特异性检查', 'WARN')
        return {}
    qlen = {h.split()[0]: len(s) for h, s in iter_fasta(amplicons_fasta)}
    db_dir = _ascii_work_base('vp_blast')
    prefix = os.path.join(db_dir, 'host')
    import glob as _g
    if not _g.glob(prefix + '.n??'):
        if logger:
            logger.log("构建宿主 BLAST 库（仅首次，1.8GB 基因组约需 10-30 分钟）...", "WARN")
        run_cmd([makeblastdb, '-in', host, '-dbtype', 'nucl', '-out', prefix],
                logger=logger)
    out_tsv = os.path.join(db_dir, 'amp_hits.tsv')
    # -max_target_seqs 1 → 10：宿主基因组高度重复时 1 条命中可能不是最强
    # 的那条，只留 1 条会低估脱靶（2026-10-03 修复 P2-14③）
    run_cmd([blastn, '-query', amplicons_fasta, '-db', prefix,
             '-outfmt', '6 qseqid sseqid pident length', '-evalue', '1e-5',
             '-max_target_seqs', '10', '-out', out_tsv], logger=logger)
    hits = {}
    with safe_open(out_tsv) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            aln_len = float(parts[3])
            thr = min(100.0, float(qlen.get(parts[0], 0) or 100))
            if aln_len >= thr:
                hits.setdefault(parts[0], []).append(parts[1])
    return hits


def design_primers(sample_dir, mode='conserved', num_return=3,
                   logger=None, force=False, do_specificity=False,
                   assembly_dir=None, phylo_dir=None):
    """阶段⑥ 主入口。mode: conserved | plain"""
    step = 'primer'
    out_dir = check_path(os.path.join(sample_dir, '06_primer'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段⑥引物设计已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    all_pairs = []
    groups_used = []

    if mode == 'conserved':
        p_dir = check_path(phylo_dir or os.path.join(sample_dir, '05_phylo'),
                           must_exist=True)
        p_summary_p = os.path.join(p_dir, 'summary.json')
        with safe_open(p_summary_p) as f:
            p_summary = json.load(f)
        for g in p_summary.get('groups', []):
            if g.get('skipped'):
                continue
            aln = check_path(os.path.join(p_dir, g['dir'], 'aln.fasta'),
                             must_exist=True)
            regions = conserved_regions_and_consensus(aln)
            if logger:
                logger.log(f"组 {g['group']}: 保守区 {len(regions)} 个")
            for k, (cons, s, e) in enumerate(regions[:5]):
                pairs = design_primers_for_seq(f"{g['group']}_C{k + 1}", cons,
                                               num_return=num_return)
                for p in pairs:
                    p['target'] = f"{g['group']} 保守区{k + 1}(比对列 {s + 1}-{e})"
                all_pairs.extend(pairs)
            if regions:
                groups_used.append(g['group'])
    else:
        a_dir = check_path(assembly_dir or os.path.join(sample_dir, '03_assembly'),
                           must_exist=True)
        vfa = os.path.join(a_dir, 'viral_contigs.fasta')
        if not os.path.isfile(vfa):
            from .orf import extract_viral_contigs
            vfa, _ = extract_viral_contigs(a_dir)
        for h, s in iter_fasta(vfa):
            cid = h.split()[0]
            # 分窗设计覆盖全长
            win, step_len = 3000, 2500
            idx = 1
            for start in range(0, max(len(s) - 400, 1), step_len):
                sub = s[start:start + win]
                if len(sub) < 600:
                    # 丢弃的末端窗口必须留痕（2026-10-03 修复 P3-20）：
                    # 原先静默 break，用户以为全长已覆盖
                    if logger:
                        logger.log(f"{cid}: 末端窗口 {start + 1}-{len(s)}bp "
                                   f"仅 {len(sub)}bp（< 最短设计长度 600bp），"
                                   f"该窗口已丢弃", 'WARN')
                    break
                pairs = design_primers_for_seq(f"{cid}_W{idx}", sub,
                                               num_return=max(1, num_return - 1))
                for p in pairs:
                    p['target'] = f"{cid} {start + 1}-{min(start + win, len(s))}bp"
                all_pairs.extend(pairs)
                idx += 1

    if not all_pairs:
        msg = "未设计出引物（保守区不足或序列问题）"
        if logger:
            logger.log(msg, "WARN")
        summary = {'stage': step, 'n_primers': 0, 'reason': msg,
                   # 恒写特异性声明（2026-10-03 修复 P2-14①）：无产物也要
                   # 如实标注，避免下游把「没跑」当成「已检查且通过」
                   'specificity_checked': bool(do_specificity),
                   'specificity_amplicon': 'skipped'}
        with safe_open(summary_file, 'wt') as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        mark_step_done(out_dir, step)
        return summary

    # 特异性检查（可选）
    spec_hits = {}
    amp_mode = 'skipped'
    spec_ok = False        # 检查是否真正跑成（失败时产物须标 FAILED 而非 NO）
    if not do_specificity and logger:
        # 默认关时必须留痕（2026-10-03 修复 P2-14②）：产物 host_hit=NA /
        # summary specificity_checked=false 之外，日志再 WARN 一次，
        # 防止把「未检查」误读成「已检查且无脱靶」
        logger.log("未做脱靶特异性检查（do_specificity=False）：primers.tsv 的 "
                   "host_hit 列全为 NA，不代表无脱靶", "WARN")
    if do_specificity:
        from .utils import write_fasta_record
        amp_fa = os.path.join(out_dir, 'amplicons.fa')
        n_real, n_ph = 0, 0
        with safe_open(amp_fa, 'wt') as f:
            for p in all_pairs:
                seq = p.get('amp_seq') or ''
                if seq:
                    n_real += 1
                else:
                    # 模板缺失/坐标越界时的降级：左臂 + 右臂反向互补串联。
                    # 仅覆盖引物区（约 40bp），不代表真实扩增子。
                    from Bio.Seq import Seq as BSeq
                    seq = p['F_seq'] + str(BSeq(p['R_seq']).reverse_complement())
                    n_ph += 1
                write_fasta_record(f, p['name'], seq)
        amp_mode = 'real' if n_real and not n_ph else ('placeholder' if not n_real else 'mixed')
        if n_ph and logger:
            logger.log(f"有 {n_ph}/{len(all_pairs)} 对引物未能取到扩增子模板序列，"
                       f"其 amplicons.fa 记录为两引物串联占位（约 40bp），"
                       f"BLAST 结果只反映引物区结合性、不代表完整扩增子", 'WARN')
        try:
            spec_hits = check_specificity_blast(amp_fa, logger=logger)
            spec_ok = True
        except Exception as e:
            if logger:
                logger.log(f"特异性检查失败（{type(e).__name__}）: {e}", "WARN")
        else:
            if logger:
                logger.log(f"特异性检查（{amp_mode} 扩增子）：{len(spec_hits)}/{len(all_pairs)} "
                           f"对引物在宿主基因组中命中（占比 {len(spec_hits) / max(1, len(all_pairs)):.1%}）")

    tsv = os.path.join(out_dir, 'primers.tsv')
    with safe_open(tsv, 'wt') as f:
        f.write("pair\ttarget\tF_primer(5'-3')\tF_pos\tF_Tm\tF_GC(%)\t"
                "R_primer(5'-3')\tR_pos\tR_Tm\tR_GC(%)\tproduct(bp)\t"
                "penalty\tself_dimer_Th\tself_hairpin_Th\t"
                "score\trecommendation\thost_hit\n")
        for p in all_pairs:
            host_hit = 'NA'          # 未检查（不代表无脱靶）
            if do_specificity:
                # 检查失败 → FAILED（如实），不得写成 NO 冒充"已检查无脱靶"
                host_hit = (('YES' if p['name'] in spec_hits else 'NO')
                            if spec_ok else 'FAILED')
            f.write(f"{p['name']}\t{p.get('target', '')}\t{p['F_seq']}\t{p['F_pos']}\t"
                    f"{p['F_tm']}\t{p['F_gc']}\t{p['R_seq']}\t{p['R_pos']}\t"
                    f"{p['R_tm']}\t{p['R_gc']}\t{p['product']}\t{p['penalty']}\t"
                    f"{p['self_any_max']}\t{p['hairpin_max']}\t"
                    f"{p.get('score', '')}\t{p.get('recommendation', '')}\t"
                    f"{host_hit}\n")

    summary = {'stage': step, 'mode': mode, 'n_primers': len(all_pairs),
               'tsv': 'primers.tsv',
               'specificity_checked': bool(do_specificity and spec_ok),
               'specificity_amplicon': amp_mode,
               'groups': groups_used}
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段⑥ 完成: {len(all_pairs)} 对引物 -> primers.tsv")
    return summary
