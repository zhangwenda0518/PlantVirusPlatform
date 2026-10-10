# -*- coding: utf-8 -*-
"""单实验执行器：由 run_bench.py 以子进程方式拉起（便于 psutil 监控整棵进程树）。

用法: python run_one.py <spec.json>
spec: {"type": "kv"|"asm"|"ctrl_salmon_k31"|"ctrl_bwa", ...}
结束前向 stdout 打印一行 JSON 结果（RESULT: {...}）。
"""
import csv
import glob
import io
import json
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..', '..', '..')))

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
# 基准鉴定库（评测快照 = 2026-09 的 kv_index，2026-09-29 起被 plant_ref 取代）。
# 路径优先级：spec['db_dir'] > 环境变量 VP_BENCH_DB > 默认旧路径
# （databases/virusref_db/kv_index/，现已不存在）。缺库时 main() 报清晰中文
# 错误（见 resolve_db）。BLAST 库默认 databases/virusref_db/blast/virus，
# 可用环境变量 VP_BENCH_BLAST 覆盖。
KV_DIR_DEFAULT = os.path.join(REPO, 'databases', 'virusref_db', 'kv_index')
KV_DIR = KV_DIR_DEFAULT
KV_REF = os.path.join(KV_DIR, 'reference.fasta')
KV_INFO = os.path.join(KV_DIR, 'reference.ref_info.tsv')
BLAST_DB_DEFAULT = os.path.join(REPO, 'databases', 'virusref_db', 'blast', 'virus')
BLAST_DB = os.environ.get('VP_BENCH_BLAST') or BLAST_DB_DEFAULT


def resolve_db(spec):
    """按 spec['db_dir'] > VP_BENCH_DB > 默认旧路径解析基准库并校验存在。

    返回 (kv_ref, kv_info, kv_dir)；缺失时抛 FileNotFoundError 并附中文修复
    指引（旧路径 kv_index 已被 plant_ref 取代，复现需重建等价快照库）。
    """
    db_dir = (spec.get('db_dir') or os.environ.get('VP_BENCH_DB')
              or KV_DIR_DEFAULT)
    ref = os.path.join(db_dir, 'reference.fasta')
    info = os.path.join(db_dir, 'reference.ref_info.tsv')
    if not (os.path.isfile(ref) and os.path.isfile(info)):
        raise FileNotFoundError(
            f'[基准库缺失] 找不到参考库：{ref}\n'
            f'说明：本评测快照绑定 2026-09 的 kv_index 库（含 TSWV 注入记录，'
            f'manifest 2026-09-19），2026-09-29 起平台默认库已改为 '
            f'virusref_db/plant_ref/，旧 kv_index 目录已移除。\n'
            f'复现请先重建等价库；或显式指定库目录：环境变量 '
            f'VP_BENCH_DB=<库目录>，或 spec 里传 "db_dir"'
            f'（目录内需有 reference.fasta / reference.ref_info.tsv）。')
    return ref, info, db_dir


def load_ref_info(kv_info=None):
    out = {}
    with open(kv_info or KV_INFO, encoding='utf-8') as f:
        for row in csv.reader(f, delimiter='\t'):
            if row and not row[0].startswith('#') and len(row) > 2:
                out[row[0]] = row[2]
    return out


def run_kv(spec, out_dir):
    """平台鉴定（known_virus_suite identify，salmon/minibwa 引擎，平台默认库）。

    库路径见 resolve_db（spec['db_dir'] > VP_BENCH_DB > 旧默认 kv_index）。
    """
    name = spec['name']
    kv_ref, kv_info, kv_dir = resolve_db(spec)
    os.makedirs(out_dir, exist_ok=True)
    r1, r2 = spec['r1'], spec.get('r2')
    # 深度稀释实验：随机抽样 N 对/条 reads（seed 默认可复现，记录进输出元数据）
    if spec.get('subsample'):
        n = int(spec['subsample'])
        seed = int(spec.get('subsample_seed', 42))
        sub1 = os.path.join(out_dir, f'sub_R1.fastq')
        _sample_fastq(r1, sub1, n, seed)
        r1 = sub1
        if r2:
            sub2 = os.path.join(out_dir, 'sub_R2.fastq')
            _sample_fastq(r2, sub2, n, seed)
            r2 = sub2
    sheet = os.path.join(out_dir, 'sheet.tsv')
    with open(sheet, 'w', encoding='utf-8') as f:
        f.write('name\tr1\n')
        f.write(f'{name}\t{r1}\n')
    engine = spec.get('engine', 'salmon')
    cmd = [sys.executable, '-m',
           'Virus_Platform_Core.known_virus_suite.known_virus_suite',
           'identify', '--out', out_dir, '--reference', kv_ref,
           '--ref-info', kv_info, '--engine', engine,
           '--align-threads', '4', '--index-dir', kv_dir,
           '--sample-sheet', sheet]
    if r2:
        # 批量表支持 r2 列
        with open(sheet, 'w', encoding='utf-8') as f:
            f.write('name\tr1\tr2\n')
            f.write(f'{name}\t{r1}\t{r2}\n')
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True,
                       encoding='utf-8', errors='replace', timeout=1800,
                       cwd=REPO)
    wall = time.time() - t0
    log = (r.stdout or '') + (r.stderr or '')
    open(os.path.join(out_dir, 'identify_run.log'), 'w',
         encoding='utf-8').write(log)

    res = {'wall_s': round(wall, 1), 'returncode': r.returncode}
    # spec 指纹（P1-7）：续跑判定据此比对——同名实验但输入/参数/脚本变了
    # 不能沿用旧 RESULT。
    import hashlib
    res['spec_sha1'] = hashlib.sha1(
        json.dumps(spec, ensure_ascii=False, sort_keys=True)
        .encode('utf-8')).hexdigest()
    if r.returncode != 0:
        # 非零退出 = 本轮无效（审查 P1-7）：不解析任何盘上产物表，直接
        # 以 error 收场（RESULT 由文件尾的统一出口打印，此处只返回，
        # 否则 stdout 会出现两行 RESULT）。此前旧 TSV 残留会被读成本轮结果。
        res['error'] = (f'identify 引擎退出码 {r.returncode}（本轮作废，'
                        '详见 identify_run.log 尾部）')
        tail = log[-400:]
        res['error_tail'] = tail
        return res
    # 采样元数据（仅深度稀释实验；记录 seed 与抽样量以便复现）
    if spec.get('subsample'):
        res['subsample'] = {'n_reads': int(spec['subsample']),
                            'seed': int(spec.get('subsample_seed', 42)),
                            'mode': 'random'}
    # 索引构建提示行（k15/k31 自动切换是否触发）
    m = re.search(r'检测到超短 reads[^\n]*', log)
    res['srna_adapted'] = bool(m)
    if m:
        res['srna_log'] = m.group(0).strip()
    bm = re.search(r'\[salmon\] 索引完成 ([\d.]+)s', log)
    if bm:
        res['index_build_s'] = float(bm.group(1))
    hm = re.search(r'\[' + name + r'\][^\n]*?mapped=(\d+)/(\d+)', log)
    if hm:
        res['mapped'], res['total'] = int(hm.group(1)), int(hm.group(2))
    # 确诊表（平台双轨过滤后的最终判定）
    conf = os.path.join(out_dir, 'identify', 'all_viruses.best.summary.tsv')
    species = []
    if os.path.isfile(conf):
        with open(conf, encoding='utf-8') as f:
            rd = csv.DictReader(f, delimiter='\t')
            for row in rd:
                species.append({'species': row.get('Species', ''),
                                'acc': row.get('Accession', ''),
                                'cov': row.get('Coverage(%)', ''),
                                'depth': row.get('MeanDepth', ''),
                                'em_reads': row.get('EM_Reads', '')})
    res['confirmed'] = species
    # 物种级定量汇总：quant.sf 的 EM reads / refstats 的 mapped 记录数，
    # 按 ref_info 物种名 join（含全库所有参考，不限确诊条目）
    info = load_ref_info()
    em, mapped = {}, {}
    qf = os.path.join(out_dir, 'quant', name, 'quant.sf')
    if os.path.isfile(qf):
        with open(qf, encoding='utf-8') as f:
            f.readline()
            for line in f:
                p_ = line.rstrip('\n').split('\t')
                if len(p_) >= 5 and float(p_[4]) > 0:
                    sp = info.get(p_[0])
                    if sp:
                        em[sp] = em.get(sp, 0.0) + float(p_[4])
    rf = os.path.join(out_dir, 'stats', f'{name}.refstats.tsv')
    if os.path.isfile(rf):
        with open(rf, encoding='utf-8') as f:
            rd = csv.DictReader(f, delimiter='\t')
            for row in rd:
                sp = info.get(row['Accession'])
                if sp:
                    mapped[sp] = mapped.get(sp, 0) + int(row['Mapped_Reads'] or 0)
    res['species_em'] = em
    res['species_mapped'] = mapped
    return res


def run_asm(spec, out_dir):
    """小RNA组装（平台 srna 模式）+ contig BLAST 物种判定（≥60bp，论文口径）。"""
    import shutil
    from Virus_Platform_Core.assembly import run_spades
    # SPAdes 输出必须在平台根内（check_path in_platform 约束）：
    # 用平台内 scratch 目录，跑完取数后清理
    scratch = os.path.join(REPO, '_zz_bench_asm', os.path.basename(out_dir))
    os.makedirs(scratch, exist_ok=True)
    try:
        t0 = time.time()
        sp = run_spades(spec['r1'], spec.get('r2'), scratch, mode='srna',
                        threads=4, memory_gb=8)
        wall = time.time() - t0
        contigs = []
        name, seq = None, []
        for ln in open(sp, encoding='utf-8'):
            if ln.startswith('>'):
                if name:
                    contigs.append((name, ''.join(seq)))
                name, seq = ln[1:].strip(), []
            else:
                seq.append(ln.strip())
        if name:
            contigs.append((name, ''.join(seq)))
        n60 = sum(1 for _, s in contigs if len(s) >= 60)
        n200 = sum(1 for _, s in contigs if len(s) >= 200)
        res = {'wall_s': round(wall, 1), 'n_contigs': len(contigs),
               'n_ge60': n60, 'n_ge200': n200,
               'total_bp': sum(len(s) for _, s in contigs)}
        # ≥60bp contigs BLAST 平台库（word_size 11 默认，e=1e-3）
        fa60 = [c for c in contigs if len(c[1]) >= 60]
        species = {}
        if fa60:
            import tempfile
            with tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False,
                                             encoding='utf-8') as tf:
                for nm, s in fa60:
                    tf.write(f'>{nm}\n')
                    for i in range(0, len(s), 70):
                        tf.write(s[i:i + 70] + '\n')
                tmpfa = tf.name
            blastn = glob.glob(os.path.join(REPO, '3rd', 'tools', 'Blast', 'bin',
                                            'blastn.exe'))[0]
            r = subprocess.run([blastn, '-query', tmpfa, '-db', BLAST_DB,
                                '-outfmt', '6 qseqid sseqid pident length',
                                '-max_target_seqs', '1', '-evalue', '1e-3',
                                '-num_threads', '4'], capture_output=True,
                               text=True, encoding='utf-8', errors='replace')
            os.unlink(tmpfa)
            info = load_ref_info()
            for ln in r.stdout.splitlines():
                p = ln.split('\t')
                if len(p) >= 4 and int(p[3]) >= 60:
                    sp_name = info.get(p[1], p[1])
                    species.setdefault(sp_name, 0)
                    species[sp_name] += 1
        res['species_contigs'] = species
        return res
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _sample_fastq(src, dst, n_records, seed=42):
    """从 src 随机抽取 n_records 条 reads（4 行/条），支持 gz 输入，输出明文 fastq。

    2026-10-03 方法变更：此前为"顺序截取前 N 条"（无种子、受 fastq 排序影响）；
    现改为 random.Random(seed).sample 全量索引后随机抽样，seed 默认 42 并写入
    输出元数据（res['subsample']）。注意：已发表的 50K 深度稀释数字基于旧的
    顺序截取实现（见 BENCHMARK_zh.md §2.3 方法变更记录），结论不受影响但抽样
    实现不同。
    """
    import gzip
    import random
    opener = gzip.open if str(src).lower().endswith('.gz') else open
    # 先流式计数记录数，再抽索引，第二遍按索引写出（不做全量驻留）
    n_total = 0
    with opener(src, 'rt', encoding='utf-8', errors='replace') as f:
        for _ in f:
            n_total += 1
    n_total //= 4
    idx = set(random.Random(seed).sample(range(n_total),
                                         min(n_records, n_total)))
    with opener(src, 'rt', encoding='utf-8', errors='replace') as f, \
            open(dst, 'w', encoding='utf-8') as w:
        for i, ln in enumerate(f):
            if i // 4 in idx:
                w.write(ln)


def run_chain(spec, out_dir):
    """sRNA 全链路实测：①fastp 质控（平台预处理，min_len=15 适配小RNA）
    → ②salmon k15 鉴定 → ③srna 组装+BLAST。分阶段计时。"""
    import shutil
    from Virus_Platform_Core.preprocess import run_fastp
    from Virus_Platform_Core.assembly import run_spades
    stages = {}
    scratch = os.path.join(REPO, '_zz_bench_chain', os.path.basename(out_dir))
    os.makedirs(scratch, exist_ok=True)
    try:
        # ① fastp
        t0 = time.time()
        prep_dir = os.path.join(scratch, '00_prep')
        run_fastp(scratch, spec['r1'], spec.get('r2'), threads=4)
        stages['fastp_s'] = round(time.time() - t0, 1)
        qc1 = os.path.join(prep_dir, 'fastp_R1.fastq.gz')
        qc2 = os.path.join(prep_dir, 'fastp_R2.fastq.gz')
        if not os.path.isfile(qc1):
            raise FileNotFoundError(f'fastp 未产出 {qc1}')

        # ② 鉴定（复用 run_kv，输入为质控后 reads）
        t0 = time.time()
        kv_res = run_kv({**spec, 'r1': qc1, 'r2': qc2 if os.path.isfile(qc2) else None,
                         'subsample': None}, os.path.join(out_dir, 'identify'))
        stages['identify_s'] = round(time.time() - t0, 1)

        # ③ 组装 + BLAST
        t0 = time.time()
        asm_res = run_asm({**spec, 'r1': qc1, 'r2': None},
                          os.path.join(out_dir, 'asm'))
        stages['assemble_s'] = round(time.time() - t0, 1)

        res = {'stages': stages,
               'confirmed': kv_res.get('confirmed', []),
               'species_em': kv_res.get('species_em', {}),
               'mapped': kv_res.get('mapped'), 'total': kv_res.get('total'),
               'n_contigs': asm_res.get('n_contigs'),
               'species_contigs': asm_res.get('species_contigs', {})}
        # 质控后 reads 存活数：直接数 fastp 输出 fastq 记录数
        # （fastp 的 JSON 报告在 Windows 下含未转义反斜杠路径，json 解析会炸）
        import gzip as _g
        n_lines = 0
        with _g.open(qc1, 'rt', encoding='utf-8', errors='replace') as f:
            for i, _ln in enumerate(f):
                n_lines = i + 1
        res['qc_survived_reads'] = n_lines // 4
        return res
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def run_ctrl_salmon_k31(spec, out_dir):
    """对照组：原始 k31 索引对小RNA数据的定量（预期 0 映射）。"""
    t0 = time.time()
    idx = spec['index']
    outq = os.path.join(out_dir, 'q')
    salmon = glob.glob(os.path.join(REPO, '3rd', 'tools', 'salmon2',
                                    'salmon.exe'))[0]
    subprocess.run([salmon, 'quant', '-i', idx, '-l', 'A', '-r', spec['r1'],
                    '--validateMappings', '-o', outq, '-p', '4'],
                   capture_output=True, text=True, timeout=1200)
    wall = time.time() - t0
    meta = json.load(open(os.path.join(outq, 'aux_info', 'meta_info.json'),
                          encoding='utf-8'))
    return {'wall_s': round(wall, 1), 'percent_mapped': meta.get('percent_mapped'),
            'num_mapped': meta.get('num_mapped')}


def run_ctrl_bwa(spec, out_dir):
    """对照组：minibwa 引擎遇小RNA应明确拒绝（预期报错，不产出全零定量）。"""
    from Virus_Platform_Core.known_virus_suite.kv_engines import BwaEngine
    t0 = time.time()
    be = BwaEngine(glob.glob(os.path.join(REPO, '3rd', 'bin', 'minibwa.exe'))[0],
                   samtools='x', threads=1)
    try:
        be.align('unused-index', {'name': spec['name'], 'r1': spec['r1']},
                 out_dir)
        return {'wall_s': round(time.time() - t0, 1), 'rejected': False}
    except ValueError as e:
        return {'wall_s': round(time.time() - t0, 1), 'rejected': True,
                'message': str(e)[:120]}


def main():
    spec = json.load(open(sys.argv[1], encoding='utf-8'))
    out_dir = spec['out']
    os.makedirs(out_dir, exist_ok=True)
    typ = spec['type']
    fn = {'kv': run_kv, 'asm': run_asm, 'chain': run_chain,
          'ctrl_salmon_k31': run_ctrl_salmon_k31, 'ctrl_bwa': run_ctrl_bwa}[typ]
    res = fn(spec, out_dir)
    res['type'] = typ
    res['name'] = spec['name']
    print('RESULT:' + json.dumps(res, ensure_ascii=False), flush=True)
    sys.stdout.flush()
    os._exit(0)          # 硬退出：不等待任何残留非守护线程/句柄（防挂起）


if __name__ == '__main__':
    main()
