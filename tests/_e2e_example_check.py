# -*- coding: utf-8 -*-
"""E2E：用内置示例 reads 真实跑一遍 ⓪~⑩ 可用阶段，验证审查修复在真实
数据流下的行为。跑完自动清理样品目录。

用法: python tests/_e2e_example_check.py
"""
import io
import json
import os
import sys
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault('VP_NO_RECOVER', '1')

SAMPLE = 'example_e2e'
FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name
          + (f' — {detail}' if detail and not cond else ''), flush=True)
    if not cond:
        FAIL.append(name)


def main():
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.utils import check_path, safe_sample_name
    from Virus_Platform_Core.pipeline import save_sample_input
    r1 = os.path.join(ROOT, 'examples', 'example_R1.fastq.gz')
    r2 = os.path.join(ROOT, 'examples', 'example_R2.fastq.gz')
    sample = safe_sample_name(SAMPLE)
    sd = check_path(os.path.join(DIRS['results'], sample),
                    must_exist=False, in_platform=True)
    import shutil
    shutil.rmtree(sd, ignore_errors=True)      # 干净起点
    os.makedirs(sd, exist_ok=True)
    save_sample_input(sd, r1, r2, sample)

    from Virus_Platform_Core.pipeline import (DEFAULT_ANALYZE_STAGES,
                                              run_analysis)
    t0 = time.time()
    log_path = os.path.join(sd, 'e2e_run.log')
    logf = io.open(log_path, 'w', encoding='utf-8')

    class Log:
        def log(self, m, lvl=None):
            logf.write(m + '\n')
            logf.flush()

        def warning(self, m):
            self.log('WARN ' + m)

        def close(self):
            logf.close()
    lg = Log()
    from Virus_Platform_Core.config import get_config
    _cfg = get_config()
    try:
        run_analysis(sample=sample, r1=r1, r2=r2,
                     stages=list(DEFAULT_ANALYZE_STAGES),
                     db_host=_cfg.databases['host'],
                     db_virus=_cfg.databases['virus'],
                     threads=8, logger=Log(), progress=None)
    except Exception as e:
        print(f'✗ 流水线异常: {type(e).__name__}: {e}')
        logf.close()
        return 1
    finally:
        logf.close()
    wall = time.time() - t0
    print(f'流水线完成: {wall:.0f}s（日志 {log_path}）')

    # ---------- 修复行为核对（真实产物） ----------
    def rd(rel):
        p = os.path.join(sd, rel)
        if not os.path.isfile(p):
            return None
        try:
            return json.load(open(p, encoding='utf-8'))
        except (OSError, ValueError):
            return None

    fp_json = rd('00_prep/fastp.json')
    # fastp 不在 DEFAULT_ANALYZE_STAGES（可选阶段）：跑了才验 params
    check('P1-1 fastp 摘要含 params（dedup/min_len）',
          (not fp_json) or ('params' in fp_json))
    hs = rd('01_host_removal/stats.json')
    check('P1-1 host 摘要含 classify_inputs/host_db（复用判据数据源）',
          bool(hs) and 'classify_inputs' in hs and 'host_db' in hs)
    kv = rd('02b_kvsuite/summary.json')
    ifp = (kv or {}).get('inputs_fp') or {}
    check('P1-1 kvsuite 指纹含 engine/ref_info/thresholds',
          'engine' in ifp and 'ref_info' in ifp and 'thresholds' in ifp)
    _full = os.path.join(sd, (kv or {}).get('identify_tsv', '')
                         or '02b_kvsuite/identify/all_viruses.summary.tsv')
    _full_txt = (io.open(_full, encoding='utf-8', errors='replace')
                 .readline() if os.path.isfile(_full) else '')
    check('P1-6 全量表含 ani_screened/Cov_Depth_Basis 列',
          'ani_screened' in _full_txt and 'Cov_Depth_Basis' in _full_txt,
          _full_txt[:120])
    _best = os.path.join(sd, '02b_kvsuite', 'identify',
                         'all_viruses.best.summary.tsv')
    if os.path.isfile(_best):
        _bt = io.open(_best, encoding='utf-8', errors='replace').readline()
        check('P1-6 best 表含 ani_status 列', 'ani_status' in _bt,
              _bt[:120])
    else:
        check('P1-6 本轮无 best 表（全被剔除/未出候选，属合法结果）', True)
    asm = rd('03_assembly/summary.json')
    check('P1-1 assembly 摘要含 inputs_fp', bool((asm or {}).get('inputs_fp')))
    check('P1-3 unclassified 候选池计数就位',
          'unclassified_pool_n' in (asm or {}))
    ufa = os.path.join(sd, '03_assembly', 'unclassified_contigs.fasta')
    check('P1-3 unclassified_contigs.fasta 已产出',
          os.path.isfile(ufa))
    vb = rd('03b_verify/summary.json')
    check('P1-3 主验证摘要就绪', bool(vb))
    # P0-1：identify 空结果契约——本轮有命中则为真实表；零命中则为空表
    itsv = os.path.join(sd, kv.get('identify_tsv', '')
                        if kv else '02b_kvsuite/identify/all_viruses.summary.tsv')
    check('P0-1 identify 表本轮落盘（有命中=数据表；零命中=带表头空表）',
          os.path.isfile(itsv))
    rep = os.path.join(sd, '07_report', 'report.html')
    check('⑩ 报告已产出', os.path.isfile(rep))

    # ---------- 第二遍：新格式指纹的复用命中（不重跑） ----------
    from Virus_Platform_Core.pipeline import run_analysis as _ra
    t1 = time.time()
    logf2 = io.open(os.path.join(sd, 'e2e_run2.log'), 'w',
                    encoding='utf-8')

    class Log2:
        def log(self, m, lvl=None):
            logf2.write(m + '\n')
            logf2.flush()

        def warning(self, m):
            self.log('WARN ' + m)

        def close(self):
            logf2.close()
    run_analysis(sample=sample, r1=r1, r2=r2,
                 stages=list(DEFAULT_ANALYZE_STAGES),
                 db_host=_cfg.databases['host'],
                 db_virus=_cfg.databases['virus'],
                 threads=8, logger=Log2(), progress=None)
    wall2 = time.time() - t1
    check(f'第二遍全部复用命中（{wall2:.0f}s，远快于首遍 {wall:.0f}s）',
          wall2 < max(wall / 3, 10))
    # 各阶段 .done/摘要未被重写（复用即不重算）：比对 mtime 前后不变由
    # 耗时差佐证，另抽 kvsuite summary 的 elapsed_sec 不变
    kv2 = rd('02b_kvsuite/summary.json')
    check('第二遍 kvsuite 摘要未被重算',
          (kv2 or {}).get('elapsed_sec') == (kv or {}).get('elapsed_sec'))

    # ---------- 清理 ----------
    shutil.rmtree(sd, ignore_errors=True)
    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 示例数据 E2E 全部通过（样品目录已清理）')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
