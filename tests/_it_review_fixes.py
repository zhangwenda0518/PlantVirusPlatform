# -*- coding: utf-8 -*-
"""审查报告（2026-10-10）全量修复的回归测试。

覆盖：P0-1（零命中旧阳性污染）、P1-2（随机子采样）、P1-6（ani_status）、
P1-1/P1-3/P1-5/P1-7/P1-8/P1-9/P1-10 与 P2 的接线断言（源级）+ 行为断言
（可离线验证的部分）。

用法: python tests/_it_review_fixes.py
退出码: 0 = 全部通过；1 = 存在失败。
"""
import gzip
import io
import json
import os
import sys
import tempfile
from types import SimpleNamespace

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault('VP_NO_RECOVER', '1')

FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name
          + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


class _L:
    def info(self, m):
        pass

    def warning(self, m):
        pass


def _platform_tmp(tag):
    """平台内临时目录（safe_open 白名单要求）。"""
    d = os.path.join(ROOT, 'run', 'logs', f'_rvw_{tag}')
    os.makedirs(d, exist_ok=True)
    return d


def test_p0_1():
    print('== P0-1 零命中旧阳性污染 ==')
    from Virus_Platform_Core.known_virus_suite.kv_identify import IdentifyStage
    d = _platform_tmp('p01')
    ident = os.path.join(d, 'identify')
    os.makedirs(ident, exist_ok=True)
    ref = os.path.join(d, 'ref.fa')
    open(ref, 'w').write('>NC_1\n' + 'ACGT' * 1500 + '\n')
    hdr = ('Sample\tAccession\tLength\tSpecies\tEM_Reads\tAvg_Read_Len\t'
           'Coverage(%)\tAvg_Read_ANI\n')
    old = os.path.join(ident, 'all_viruses.best.summary.tsv')
    open(old, 'w', encoding='utf-8').write(
        hdr + 's\tNC_1\t6000\tPiVB\t100\t150\t80.0\t95.0\n')
    open(os.path.join(ident, 'all_viruses.summary.tsv'), 'w',
         encoding='utf-8').write(hdr + 'x\n')
    st = IdentifyStage(SimpleNamespace(reference=ref, ref_info=None),
                       None, {}, _L(), d)
    s, c, n = st.summarize([])
    check('零命中：旧 best 表被清除', not os.path.isfile(old))
    body = open(os.path.join(ident, 'all_viruses.summary.tsv'),
                encoding='utf-8').read()
    check('零命中：写带表头空表', body.startswith('Sample\t')
          and body.count('\n') == 1)
    check('新契约列就位（ani_status/Cov_Depth_Basis）',
          'ani_status' in body and 'Cov_Depth_Basis' in body)

    from Virus_Platform_Core.known_virus_suite.kv_filter import FilterStage
    fdir = os.path.join(d, 'filter')
    os.makedirs(fdir, exist_ok=True)
    # 上一轮残留的旧通过表
    open(os.path.join(fdir, 'filtered.tsv'), 'w').write(hdr + 'stale\n')
    # 2026-10-20：polars → pandas（FilterStage 已改用 pandas）
    import pandas as pd
    fs = FilterStage(SimpleNamespace(min_reads=0, min_depth=0, min_cov=0,
                                     ratio=0, ani_thresh=95, min_tpm=0),
                     _L(), d)
    fs.apply(pd.DataFrame({'Sample': pd.Series(dtype='object')}))
    check('filter：空输入覆盖旧 filtered.tsv',
          open(os.path.join(fdir, 'filtered.tsv'), encoding='utf-8').read()
          .count('\n') <= 1)


def test_p1_6_ani_status():
    print('== P1-6 ani_status 三态 ==')
    from Virus_Platform_Core.known_virus_suite.kv_identify import IdentifyStage
    d = _platform_tmp('p16')
    ref = os.path.join(d, 'ref.fa')
    open(ref, 'w').write('>A\n' + 'ACGT' * 1500 + '\n>B\n' + 'ACGT' * 1500
                         + '\n>C\n' + 'ACGT' * 1500 + '\n')
    base = {'Length': 6000, 'taxid': 1, 'Segment': '',
            'Molecule_Type2': 'RNA', 'Host': '', 'Geo_Location': '',
            'GenBank_Title': '', 'Coverage(%)': 90.0, 'MeanDepth': 20.0,
            'Cov_Depth_Basis': 'est_reads_x_len', 'Uniq_Reads': 100,
            'Sample_Total_Mapped': 1000, 'Avg_Read_Len': 150.0}
    rows = [
        dict(base, Sample='s', Accession='A', Species='spA', EM_Reads=100.0,
             Avg_Read_ANI=None),
        dict(base, Sample='s', Accession='B', Species='spB', EM_Reads=90.0,
             Avg_Read_ANI=97.0),
        dict(base, Sample='s', Accession='C', Species='spC', EM_Reads=80.0,
             Avg_Read_ANI=90.0),
    ]
    st = IdentifyStage(
        SimpleNamespace(reference=ref, ref_info=None, min_reads=0,
                        min_depth=0, min_cov=0, ratio=0.0, ani_thresh=95,
                        genes_cov=None),
        None, {}, _L(), d)
    summ, conf, novel = st.summarize(rows)
    st_map = dict(zip(conf['Accession'].to_list(),
                      conf['ani_status'].to_list()))
    check('null ANI → ani_unavailable（在 best 表内但自带标注）',
          st_map.get('A') == 'ani_unavailable')
    check('97% → ani_confirmed', st_map.get('B') == 'ani_confirmed')
    n_map = dict(zip(novel['Accession'].to_list(),
                     novel['ani_status'].to_list()))
    check('90% → ani_below_threshold（novel 表）',
          n_map.get('C') == 'ani_below_threshold')


def test_p1_2_random():
    print('== P1-2 随机子采样 ==')
    from Virus_Platform_Core.pipeline import _subsample_random
    d = _platform_tmp('p12')
    r1 = os.path.join(d, 'a_R1.fq.gz')
    r2 = os.path.join(d, 'a_R2.fq.gz')

    def wq(p, n):
        with gzip.open(p, 'wt') as f:
            for i in range(n):
                f.write(f'@a{i}\nACGTACGTAC\n+\nIIIIIIIIII\n')
    wq(r1, 100)
    wq(r2, 100)
    d1, d2 = os.path.join(d, 'o1'), os.path.join(d, 'o2')
    os.makedirs(d1)
    os.makedirs(d2)
    s1, s2 = _subsample_random(r1, r2, d1, 10, None, None, seed=7)
    ids = [l[1:].strip() for l in gzip.open(s1, 'rt') if l.startswith('@')]
    ids2 = [l[1:].strip() for l in gzip.open(s2, 'rt') if l.startswith('@')]
    check('数量正确', len(ids) == 10)
    check('R1/R2 配对同步', ids == ids2)
    s1b, _ = _subsample_random(r1, r2, d2, 10, None, None, seed=7)
    ids_b = [l[1:].strip() for l in gzip.open(s1b, 'rt') if l.startswith('@')]
    check('同 seed 完全复现', ids == ids_b)
    check('覆盖后段（非前缀）', max(int(i[1:]) for i in ids) > 10)



def test_p1_1_host_reuse():
    """P1-1 host 复用三路：命中 / 相对路径同库（口径修复）/ 参数变重跑。"""
    import gzip as gzm
    from Virus_Platform_Core.utils import mark_step_done
    from Virus_Platform_Core.host_removal import remove_host, _finfo
    d = _platform_tmp('p11h')
    sd = os.path.join(d, 's')
    h1 = os.path.join(sd, '01_host_removal')
    os.makedirs(h1, exist_ok=True)
    fq = os.path.join(h1, 'in_R1.fastq.gz')
    with gzm.open(fq, 'wt') as f:
        f.write('@a\nACGT\n+\nIIII\n')
    db_abs = os.path.abspath(os.path.join(d, 'hostdb'))
    stats = {'stage': 'host_removal', 'marker': 'MINE',
             'classify_inputs': [_finfo(fq)],
             'classify_params': {'confidence': 0.0, 'threads': None,
                                 'min_hit_groups': 2},
             'host_db': {'path': db_abs, 'name': 'hostdb'}}

    def arm():
        json.dump(stats, open(os.path.join(h1, 'stats.json'), 'w',
                              encoding='utf-8'), ensure_ascii=False,
                  default=str)
        mark_step_done(h1, 'host_removal')
    arm()
    r = remove_host(sd, fq, None, db_abs, confidence=0.0)
    check('host 复用：全一致命中', r.get('marker') == 'MINE')
    arm()
    try:
        r2 = remove_host(sd, fq, None,
                         os.path.relpath(db_abs, os.getcwd()),
                         confidence=0.0)
        check('host 复用：相对路径同库仍命中（abspath 口径）',
              r2.get('marker') == 'MINE')
    except FileNotFoundError:
        check('host 复用：相对路径同库仍命中（abspath 口径）', False)
    arm()
    try:
        remove_host(sd, fq, None, db_abs, confidence=0.5)
        check('host 复用：confidence 变 → 重跑', False)
    except FileNotFoundError:
        check('host 复用：confidence 变 → 重跑', True)


def test_wiring():
    print('== 接线断言（P1-1/3/5/7/8/9/10 + P2）==')
    def has(path, needle):
        try:
            return needle in io.open(os.path.join(ROOT, path),
                                     encoding='utf-8').read()
        except OSError:
            return False

    check('P1-1 subsample 指纹含 mode/seed/inputs',
          has('Virus_Platform_Core/pipeline.py', "prev.get('mode') == C.subsample_mode"))
    check('proportion 子采样三路接线（CLI/Web/AI）',
          has('main.py', 'subsample_mode=args.subsample_mode')
          and has('Virus_Platform_Core/web/samples.py', "subsample_mode=_val('subsample_mode'")
          and has('Virus_Platform_Core/ai_agent/executors.py', "'subsample_mode'"))
    check('P1-1 fastp 指纹含 params',
          has('Virus_Platform_Core/preprocess.py', "'params': {'dedup': bool(dedup), 'min_len': min_len}"))
    check('P1-1 host 复用校验输入/库/参数',
          has('Virus_Platform_Core/host_removal.py', '_db_ok = (_prev_db in (None'))
    check('P1-1 assembly .done 校验 inputs_fp',
          has('Virus_Platform_Core/assembly.py', "_prev.get('inputs_fp') == _cur_fp"))
    check('P1-1 kvsuite 外层指纹含 engine/thresholds',
          has('Virus_Platform_Core/kv_stage.py', "'thresholds': _default_thresholds(),"))
    check('P1-3 未分类候选池（assembly 产出）',
          has('Virus_Platform_Core/assembly.py', 'unclassified_contigs.fasta'))
    check('P1-3 未分类独立验证（pipeline 通道）',
          has('Virus_Platform_Core/pipeline.py', "'03b_verify_unclassified'"))
    check('P1-5 检索源状态落盘 + 全失败退出',
          has('Virus_Platform_Core/public_meta/search_engine.py', 'source_status.json')
          and has('Virus_Platform_Core/public_meta/search_engine.py', 'sys.exit(2)'))
    check('P1-7 run_one 非零退出作废',
          has('docs/paper/benchmark/run_one.py', "res['error'] = (f'identify 引擎退出码"))
    check('P1-7 run_bench 按 spec_sha1 续跑',
          has('docs/paper/benchmark/run_bench.py',
              "done.add(_e)"))
    tj = json.load(io.open(os.path.join(
        ROOT, 'docs/paper/benchmark/truth_viromock_d110_groups.json'),
        encoding='utf-8'))
    check('P1-8 D9 真值标 ambiguous', bool(tj['9'][0].get('ambiguous')))
    check('P1-8 评分隔离 ambiguous（敏感性单独成行）',
          has('docs/paper/benchmark/compute_and_plot.py',
              'VIROMOCK-D110(敏感性)'))
    check('P1-9 窗口表 N_EFF_BASIS + 探索性声明',
          has('Virus_Platform_Core/known_virus_suite/kv_variant_evo.py',
              'EXPLORATORY Pool-Seq statistic')
          and has('Virus_Platform_Core/known_virus_suite/kv_variant_evo.py',
                  'N_EFF_BASIS'))
    check('P1-10 BEAST 等距抽样',
          has('Virus_Platform_Core/beast_io.py', 'linspace'))
    check('P1-10 无 log 的树导入带探索性警告',
          has('Virus_Platform_Core/phylogeo.py', '未做 ESS/收敛诊断'))
    from Virus_Platform_Core.rtt_drt import DEFAULT_N_PERM
    check('P2 DRT 默认置换 999', DEFAULT_N_PERM == 999)
    ms = io.open(os.path.join(ROOT, 'docs/paper/manuscript_body.md'),
                 encoding='utf-8').read()
    check('P2 手稿 min contig 200bp', 'default 200 bp' in ms)
    check('P2 手稿 BFVD（非 CDD）', 'BFVD virus-specific domain library' in ms
          and 'ani_unavailable' in ms)
    check('P2 步数口径统一（README 15 阶段）',
          has('README.md', '15 个阶段') and '13 步级联' not in
          io.open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read())


def main():
    test_p0_1()
    test_p1_6_ani_status()
    test_p1_2_random()
    test_prop_subsample()
    test_p1_1_host_reuse()
    test_wiring()
    import shutil
    shutil.rmtree(os.path.join(ROOT, 'run', 'logs', '_rvw_p01'), ignore_errors=True)
    shutil.rmtree(os.path.join(ROOT, 'run', 'logs', '_rvw_p16'), ignore_errors=True)
    shutil.rmtree(os.path.join(ROOT, 'run', 'logs', '_rvw_p12'), ignore_errors=True)
    shutil.rmtree(os.path.join(ROOT, 'run', 'logs', '_rvw_psub'), ignore_errors=True)
    shutil.rmtree(os.path.join(ROOT, 'run', 'logs', '_rvw_p11h'), ignore_errors=True)
    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 审查修复回归全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
