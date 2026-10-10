# -*- coding: utf-8 -*-
"""dsRNA 全链路集成测试（需先启动平台？不需要——直接驱动引擎层）。

覆盖：
1. run_chain 端到端：合成目标 → 选窗 → 效价 → Coccinella 面板脱靶 → 安全 →
   T7 引物 → 成品 QC，校验产物文件与结构
2. 确定性：同 seed 两次选窗结果逐字节一致（det 构建）
3. 内存守卫：离谱 iterations 触发 MemoryBudgetExceeded
4. primer 限定最优窗（scope=best_window 时扩增子落在窗内）

环境前置（2026-10-03 起自足化，不再静默失败或依赖仓库外数据）：
- databases/dsrna/{lookup.db, siRNA_parameters.txt}（小文件，随发行包/仓库旁
  数据提供；缺失时本测试直接失败并说明）；
- 脱靶面板：优先用 databases/dsrna/panel/ 下真实的
  Coccinella_septempunctata* 面板；缺失时在平台 run/ 下合成一份小面板并
  临时把引擎的 PANEL_DIR 指过去（面板扫描代码路径照跑，命中数不参与断言）。
- 工作目录放平台 run/ 下（gitignored）：引擎产物写入走 safe_open，要求落在
  平台写根内 —— 与 web 层真实调用口径一致（此前写 %TEMP%，真实运行不存在
  这种配置）。

运行：python tests/_it_dsrna.py
"""
import json
import os
import random
import shutil
import sys
import time
from pathlib import Path

import numpy as np

PLATFORM = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PLATFORM))

from Virus_Platform_Core import dsrna_pipeline as dp          # noqa: E402
from Virus_Platform_Core import dsrna_scoring as dscore       # noqa: E402
from Virus_Platform_Core import dsrna_offtarget as dso        # noqa: E402

# 数据面前置：lookup.db 是 dsRIP 效价的必需查找表（256 条 RNA 四核苷酸 dG），
# 不带它会以 ScoringDataError 中断——这里先给一句可执行的说明再退出。
if not Path(dscore.LOOKUP_DB).exists():
    print(f'✘ 缺少热力学查找表 {dscore.LOOKUP_DB}\n'
          f'  → 把 lookup.db（与 siRNA_parameters.txt）放入 databases/dsrna/ 后重跑')
    sys.exit(2)
if not Path(dscore.PARAMS_PATH).exists():
    print(f'△ 参数文件 {dscore.PARAMS_PATH} 不存在：本次按引擎内置默认参数跑'
          f'（load_params 的既定行为）')


def _random_seq(n, seed):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def _ensure_panel(dp, scratch):
    """脱靶面板：真实面板优先，缺失则合成小面板并临时重指引擎 PANEL_DIR。

    返回 (是否用了合成面板, 说明文字)。
    """
    real = sorted(Path(dp.PANEL_DIR).glob('Coccinella_septempunctata*'))
    if real:
        return False, f'真实面板 {real[0].name}'
    pan = Path(scratch) / 'panel'
    pan.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20260916)
    with open(pan / 'Coccinella_septempunctata_synthetic_test.fa', 'w',
              encoding='utf-8') as fh:
        for i, n in enumerate((4000, 3500)):
            s = ''.join(rng.choice('ACGT') for _ in range(n))
            fh.write(f'>XM_TEST{i + 1}\n{s}\n')
    dp.PANEL_DIR = pan
    dp.CACHE_DIR = pan / '_cache'
    return True, f'合成小面板（真实 panel 未装）{pan}'


def main():
    ok = True
    # 工作目录在平台 run/ 下（safe_open 写盘要求平台写根内；见模块 docstring）
    tmp = PLATFORM / 'run' / f'_it_dsrna_{os.getpid()}_{int(time.time())}'
    if tmp.exists():
        shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    print(f'临时工作目录: {tmp}')
    targets = tmp / 'targets.fa'
    t1 = _random_seq(900, 1)
    t2 = _random_seq(900, 2)
    targets.write_text(f'>t1\n{t1}\n>t2\n{t2}\n', encoding='utf-8')

    synthetic, panel_note = _ensure_panel(dp, tmp)
    print(f'脱靶面板: {panel_note}'
          + ('（仅测结构，命中数不参与断言）' if synthetic else ''))

    def log(m):
        print('  [log]', m)

    def prog(stage, frac, msg):
        print(f'  [prog {stage} {frac:.2f}] {msg}')

    params = {
        'construct_len': 300, 'kmer_len': 21, 'iterations': 30, 'seed': 42,
        'candidates': 2, 'n_windows': 1,
        'panel': 'Coccinella_septempunctata', 'max_mm': 1,
        'primer_scope': 'best_window',
        'primer_product_min': 200, 'primer_product_max': 300,
    }
    run_dir = tmp / 'run1'
    run_dir.mkdir()
    print('== 1) run_chain 端到端 ==')
    summary, report = dp.run_chain(str(targets), params, run_dir, log, prog, None)
    out_dir = Path(summary['out_dir'])
    for f in ('report.json', 'primers.tsv', 'arms.fa', 'amplicons.fa',
              'offtarget_lethal.tsv', 'result.csv'):
        p = out_dir / f
        print(f'  artifact {f}:', 'OK' if p.exists() else 'MISSING')
        ok &= p.exists()
    assert report['arms'], '没有产出臂'
    assert report['primers'] and report['primers'][0]['pairs'], '没有产出引物'
    assert report['offtarget_panel'], '面板扫描未运行'
    assert 'summary' in report['offtarget_panel']
    # QC：scope=best_window 时扩增子应在最优窗内
    for q in report['amplicon_qc']:
        if q['ok']:
            assert q['window_overlap'] is not None, 'QC 缺窗内占比'
            assert q['window_overlap'] > 0.99, \
                f"扩增子不在最优窗内（占比 {q['window_overlap']}）"
    print('  QC 窗内占比断言通过')
    # report.json 可 JSON 化（无 numpy 类型泄漏）
    json.dumps(report, ensure_ascii=False)[:100]
    print('  report.json 可序列化')

    print('== 2) 同 seed 确定性 ==')
    arms2, _logs2 = dp.select_windows(
        [('t1', t1), ('t2', t2)], dp.probe_exe(), tmp, params, log, prog, None)
    assert arms2 and arms2[0]['seq'] == report['arms'][0]['seq'], \
        '同 seed 两次选窗结果不一致'
    print('  两次臂序列逐字节一致')

    print('== 3) 内存守卫 ==')
    try:
        dp.guard_memory([t1, t2], 21, 10**7)
        print('  FAIL: 未触发 MemoryBudgetExceeded')
        ok = False
    except dp.MemoryBudgetExceeded as e:
        print('  正确触发:', str(e)[:60], '…')

    # 平台自带夹具目录（先于 step 4 定义：step 4 的 GCVA 臂已内置于此）
    fx = PLATFORM / 'tests' / 'fixtures' / 'dsrna'

    print('== 4) 效价模块口径（对照 C-host 管线 v4 结果）==')
    # 夹具已内置到平台（2026-09-18 前此处是跨仓绝对路径
    # r'D:/桌面/C-host_classify/…/gcva_cocktail_arms.fa'，迁移任一侧就会断链；
    # 现改为平台自带，且在 .exists() 缺失时直接失败而不是静默跳过 —— 静默跳过
    # 会让这个回归护栏在迁移后无声失效）
    gcva = fx / 'gcva_cocktail_arms.fa'
    assert gcva.exists(), f'缺少夹具 {gcva}'
    rep = dscore.score_fasta(gcva)
    w933 = rep.get('GCVA_dsRNA_w933', {})
    bw = w933.get('best_window', {})
    assert bw.get('nt_start') == 7 and bw.get('nt_end') == 296, \
        f"w933 最优窗漂移: {bw.get('nt_start')}..{bw.get('nt_end')}"
    print('  w933 best window 7..296 与管线 v4 一致')

    # ================================================================
    # 论文已发表序列做夹具（Fletcher et al. 2025, NAR Genom Bioinform,
    # lqaf064 Supplement Table 1/2 + Figure 3）。这些臂有三种基础生测数据
    # 背书（Suppl. Tables 3-5）。可机器核验的保真锚点：
    #   - 长度恰为 300 nt（PDF 文本层抽取核对）
    #   - TEF-17 前 168 nt 与 Ce-EF1a 比对同源度 = 85.1%（Figure 3 图注），
    #     TEF-21 与其 300nt 比对段 = 79.3%（实测 85.12 / 79.00，容差 0.7pp；
    #     任何一处抄写错误都会打破这两个百分比）
    #   注意：Figure 1/4 截图里的运行产物与 Table 1 最终发表序列不是同一
    #   分子（同前缀、不同迭代，GC 45.9/45.7 vs Table1 的 45.00/42.67），
    #   图注 GC 不能作为 Table 1 序列的断言依据（论文内部口径差异，勿改）。
    # ================================================================

    print('== 5) 已发表 TEF 臂：转录保真 + dsRIP 效价 + T7 约定 ==')
    arms_fa = fx / 'tef17_tef21_arms.fa'
    pub = dict(dscore.read_fasta(arms_fa))
    tef17 = pub['TEF-17']
    tef21 = pub['TEF-21']

    def _gc(s):
        return 100.0 * sum(c in 'GC' for c in s) / len(s)

    assert len(tef17) == 300 and len(tef21) == 300, \
        f'TEF 臂长度异常: {len(tef17)}/{len(tef21)}（PDF Table 1 = 300nt）'
    print(f'  长度: {len(tef17)}/{len(tef21)} · GC: {_gc(tef17):.2f}% / '
          f'{_gc(tef21):.2f}%（图内 GC 45.9/45.7 属 Figure 运行产物，非 Table 1）')
    # dsRIP 效价在真实已发表分子上可算、量纲正常（lookup.db 覆盖全部 4-mer）
    for nm, s in (('TEF-17', tef17), ('TEF-21', tef21)):
        r = dscore.score_sequence(nm, s)
        assert 0 <= r['median_score'] <= 130 and r['n_siRNAs'] == len(s) - 20, nm
        print(f'  dsRIP {nm}: median {r["median_score"]} · '
              f'best window nt {r["best_window"]["nt_start"]}'
              f'-{r["best_window"]["nt_end"]}')
    # T7 启动子约定与论文 Table 2 完全一致（6 条已发表引物全部以此为 5' 尾）
    paper_primers = [
        'TAATACGACTCACTATAGGGAAGTACTACATCACCATCATCGATGC',   # CB-c.el-TEF-T7-F1
        'TAATACGACTCACTATAGGGGTATCCGATCTTCTTGATGAATCCAG',   # CB-c.el-TEF-T7-R1
        'TAATACGACTCACTATAGGGAAATACTATGTCACAATTATCG',       # CB-mu-21nt-TEF-F
        'TAATACGACTCACTATAGGGGTAACCAATCTTCTTAATAAAGC',      # CB-mu-21nt-TEF-R
        'TAATACGACTCACTATAGGGTCAAGAACATGATTACTGGTAC',       # CB-mu-17nt-TEF-F
        'TAATACGACTCACTATAGGGTAGCATTTCCATCCTTCCTTTC',       # CB-mu-17nt-TEF-R
    ]
    assert all(p.startswith(dscore.T7_PROMOTER) for p in paper_primers), \
        '论文 T7 尾与平台 T7_PROMOTER 常量不一致'
    assert all(18 <= len(p) - len(dscore.T7_PROMOTER) <= 28 for p in paper_primers)
    print(f'  T7 约定: 平台 {dscore.T7_PROMOTER} == 论文 Table 2 全部 6 条引物的 5\' 尾')

    print('== 6) Figure 3 同源度对拍 + 扫描器错配分层（真实数据）==')
    aln = dict(dscore.read_fasta(fx / 'ce_ef1a_alnpartners.fa'))
    ce21 = aln['ce_ef1a_aln_partner_for_TEF-21_Fig3_300nt']
    ce17 = aln['ce_ef1a_aln_partner_for_TEF-17_Fig3_168nt']
    # 等长无缝比对下的同源度对拍（±0.7pp 容差；双向校验两段抄写）
    ident21 = 100.0 * sum(a == b for a, b in zip(tef21, ce21)) / len(ce21)
    ident17 = 100.0 * sum(a == b for a, b in zip(tef17[:len(ce17)], ce17)) / len(ce17)
    assert abs(ident21 - 79.3) <= 0.7, f'TEF-21 vs Ce-EF1a {ident21:.2f}% != 79.3%'
    assert abs(ident17 - 85.1) <= 0.7, f'TEF-17 vs Ce-EF1a {ident17:.2f}% != 85.1%'
    print(f'  同源度对拍: TEF-21 {ident21:.1f}% (图注 79.3) · '
          f'TEF-17 {ident17:.1f}% (图注 85.1)')
    # 发表设计的语义 = 与 C. elegans 无连续 ≥17nt 精确匹配 → 21-mer 口径下
    # mm0 必为 0；~20% 错配密度下 1mm 最多零星出现
    frag = ce21 + ce17
    frag_rc = dso.revcomp(frag)
    both = list(dso.encode_kmers(frag.encode())) + \
        list(dso.encode_kmers(frag_rc.encode()))
    idx = {'CeEF1a': {'arr': np.unique(np.array(both, dtype=np.uint64)),
                      'lethal_arr': None}}
    res = dso.scan_arms([('TEF-17', tef17), ('TEF-21', tef21)], idx, max_mm=2)
    for nm in ('TEF-17', 'TEF-21'):
        d = res['arms'][nm]
        c = d['per_species']['CeEF1a']
        assert c['mm0'] == 0, f'{nm} 与 Ce-EF1a 出现 21-mer 精确命中（违背发表设计语义）'
        print(f'  {nm}: mm0=0 · mm1={c["mm1"]} · mm2={c["mm2"]} '
              f'（与 Figure 3 同源度口径自洽）')

    print('== 7) 已发表分子跑全链路（v1.1.15 constructLen==L 边界 + QC）==')
    # TEF-17 作为目标喂回链路：det 构建 v1.1.15 修复了 off-by-one，
    # constructLen == 目标长（300 ≤ 303）应正常出臂且为目标的子串
    pub_fa = tmp / 'tef17_target.fa'
    pub_fa.write_text(f'>TEF-17\n{tef17}\n', encoding='utf-8')
    params_pub = dict(params, construct_len=300, candidates=1, n_windows=1,
                      seed=42, panel='Coccinella_septempunctata', max_mm=1)
    run_dir2 = tmp / 'run_tef'
    run_dir2.mkdir()
    summary2, report2 = dp.run_chain(str(pub_fa), params_pub, run_dir2,
                                     log, prog, None)
    arm_pub = report2['arms'][0]['seq']
    assert arm_pub in tef17, '单目标链路的臂应是目标的精确子串'
    assert report2['primers'][0]['pairs'], '已发表分子未产出引物'
    q0 = report2['amplicon_qc'][0]
    assert q0['ok'] and q0['window_overlap'] and q0['window_overlap'] > 0.99
    print(f'  臂=目标子串 OK（{len(arm_pub)} nt）· 引物 '
          f'{summary2["n_primer_pairs"]} 对 · 扩增子 QC 窗内占比 '
          f'{q0["window_overlap"]}')

    # ================================================================
    # 审查修复（2026-09-17 第二轮）回归：以下分支在首轮集成测试中漏测
    # ================================================================

    print('== 8) isolates 覆盖度分支（修复：日志行 TypeError）==')
    # 分离株池 = 两个目标自身 + 一条 5% 突变的 t1 衍生株
    r = random.Random(9)
    mut = list(t1)
    for i in range(len(t1)):
        if r.random() < 0.05:
            mut[i] = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut[i]]
    iso_fa = tmp / 'isolates.fa'
    iso_fa.write_text(f'>t1\n{t1}\n>t2\n{t2}\n>t1_mut5pct\n{"".join(mut)}\n',
                      encoding='utf-8')
    run_dir3 = tmp / 'run_iso'
    run_dir3.mkdir()
    params_iso = dict(params, candidates=1, seed=43, isolates=str(iso_fa),
                      panel='')                      # 顺带覆盖空面板跳过
    summary3, report3 = dp.run_chain(str(targets), params_iso, run_dir3,
                                     log, prog, None)
    iso_cov = report3['isolate_coverage']
    assert iso_cov, 'isolates 参数被忽略（isolate_coverage 为空）'
    for aname, m in iso_cov.items():
        assert set(m) == {'overall', 'per_isolate'}, f'{aname} 结构异常: {set(m)}'
        assert 0.0 <= m['overall'] <= 1.0
        # 臂是某个目标的子串 → 对分离株并集覆盖应为 1，且至少对源目标完美覆盖
        assert m['overall'] == 1.0, f'{aname} overall={m["overall"]} != 1'
        assert max(m['per_isolate'].values()) == 1.0, \
            f'{aname} 应对至少一个分离株完美覆盖'
    # t1 的 5% 突变衍生株：若臂源自 t1，则对该衍生株覆盖必须 <1（逃逸信号）
    t1_cov = iso_cov[report3['arms'][0]['name']]['per_isolate']
    if t1_cov.get('t1') == 1.0:
        assert t1_cov['t1_mut5pct'] < 1.0, \
            '5% 突变衍生株的覆盖应 <1（否则覆盖度计算失真）'
    assert report3['offtarget_panel'] is None, \
        "panel='' 应显式跳过脱靶扫描（此前会错误地扫全部面板）"
    src = [k for k, v in t1_cov.items() if v == 1.0]
    print(f'  isolates 覆盖度结构 OK · 臂源分离株 {src} · overall=1 · '
          f"panel='' 跳过扫描 OK")

    print('== 9) result.csv 中选臂 + max_mm=0 语义 ==')
    # result.csv（exe 原生格式，末行为 sense arm）必须是中选臂本身——
    # 修复前固定拷 0 号候选的 CSV，maximin 重排后可能张冠李戴
    csv_txt = (tmp / 'run1' / 'dsrnamax' / 'result.csv').read_text(
        encoding='utf-8', errors='replace')
    assert report['arms'][0]['seq'] in csv_txt, \
        'result.csv 不含中选臂序列（拷贝了落选候选的 CSV）'
    print('  result.csv 含中选臂序列 OK')
    # max_mm=0：1mm 窗口必须整体无解（None），不许漏进 mm1 档
    host_arr = np.unique(np.concatenate(
        [dso.encode_kmers(t1), dso.encode_kmers(dso.revcomp(t1))]))
    mut_q = list(t1)
    mut_q[10] = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut_q[10]]
    res0 = dso.scan_arms([('mm1q', ''.join(mut_q))],
                         {'t1': {'arr': host_arr, 'lethal_arr': None}},
                         max_mm=0)
    tiers0 = res0['arms']['mm1q']['sirna_tiers']['t1']['tiers']
    # 突变在第 10 位：覆盖该位点的 11 个窗口（0..10）对 t1 无 ≤0 错配解，
    # 其余 880-11=869 个窗口仍是精确匹配（tier 0）；任何窗口都不许落 mm1
    assert tiers0.count(None) == 11 and tiers0.count(0) == 869, \
        f'max_mm=0 分层异常: None={tiers0.count(None)}, 0={tiers0.count(0)}'
    assert res0['arms']['mm1q']['per_species']['t1'] == {'mm0': 869, 'mm1': 0,
                                                         'mm2': 0}
    print('  max_mm=0：11 个跨突变窗口=None · 869 个精确窗=0 · mm1 计数 0 OK')


    if ok:
        print('DSRNA_IT_OK')
    else:
        print('DSRNA_IT_FAIL')
    shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
