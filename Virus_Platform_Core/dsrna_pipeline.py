# -*- coding: utf-8 -*-
"""dsRNA 设计全链路编排（工具任务 t-dsrna 的引擎侧）。

链路：dsRNAmax 选窗（maximin 多候选）→ siRNA 效价（dsrna_scoring）→
脱靶扫描（dsrna_offtarget，uint64 排序数组）→ 安全合成评分 → Primer3 引物
（限定在安全×效价最优窗内）→ 成品扩增子 QC（对真实交付分子重打分）。

设计要点（均来自 2026-09 对 dsRNAmax 管线的审计结论）：
- 打分对象必须 = 交付分子：引物设计限定在 best window 内，引物定稿后对
  **真实扩增子**重跑效价 + 脱靶扫描，杜绝"300nt 臂上打分、215bp 产物交付"
- median 目标会遮蔽离群株：多候选按 (最差目标, 中位, 总和) maximin 重排，
  命中 0 / 覆盖 <25% 最优目标的目标显式告警
- 内存 ≈ 160B × 去重 k-mer × iterations（iterations=100 实测 13-16KB/kmer），
  超预算 fail-fast 而不是等 OOM
- exe 兼容性：-seed 只有 det 构建认识、-kmerLen<21 必须显式 -otKmerLen，
  probe_exe() 一次性探测并缓存
"""
import json
import os
import re
import shutil
import subprocess
import threading
import time
from pathlib import Path

from .config import db_asset
from . import dsrna_scoring as dscore
from .utils import decode_output, safe_open, task_register_proc

# 数据库根优先、程序内置兜底（分发版数据库外置时 PLATFORM_ROOT 下没有面板数据）
PANEL_DIR = Path(db_asset('dsrna', 'panel'))
LETHAL_DIR = Path(db_asset('dsrna', 'lethal_lists'))
CACHE_DIR = PANEL_DIR / '_cache'

EXE_CANDIDATES = ('3rd/tools/dsrnamax/dsRNAmax_det.exe',
                  '3rd/tools/dsrnamax/dsRNAmax.exe',
                  '3rd/bin/dsrnamax.exe')

RSS_PER_KMER_PER_ITER = 160      # B；iterations=100 实测 13-16KB/kmer 反推
DEFAULT_RSS_CAP_GB = 32.0
# 单次 exe 调用的挂死上限（秒）：可用 VP_DSRNA_TIMEOUT_SEC 覆盖
EXEC_TIMEOUT_S = float(os.environ.get('VP_DSRNA_TIMEOUT_SEC') or 0) or 3600.0

_EXE_PROBE = None


class DsrnaExeError(RuntimeError):
    pass


class MemoryBudgetExceeded(ValueError):
    """预算内跑不下，fail-fast 而不是让 OOM killer 决定。"""


# ---------------------------------------------------------------- exe 探测

def _run_exe(cmd, timeout=EXEC_TIMEOUT_S, cancel=None, what='dsRNAmax'):
    """跑 exe 并回收输出：进程登记 + 超时/取消看门狗 + 平台解码口径。

    - CREATE_NO_WINDOW：GUI（pywebview）模式下不弹控制台黑框；
    - task_register_proc：任务停止（TaskManager.cancel）时可硬杀本进程，
      否则一次选窗最长挂到 timeout，用户点了「停止」也只能干等；
    - 输出按字节收集后 decode_output（utf-8→GBK）：dsRNAmax 报错会回显
      含中文的输入/输出路径，固定 utf-8+replace 会打成 U+FFFD 不可读。
    返回 (returncode, 输出文本, 耗时秒)。取消抛 RuntimeError，超时抛 DsrnaExeError。
    """
    cmd = [str(c) for c in cmd]
    creationflags = 0
    if os.name == 'nt':
        creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                creationflags=creationflags)
    except OSError as e:
        raise DsrnaExeError(f'{what} 无法执行: {e}') from e
    try:
        task_register_proc(proc)         # 任务取消 → 硬杀（非任务线程内为空操作）
    except Exception:                    # noqa: BLE001 —— 登记失败不阻塞运行
        pass

    chunks = []

    def _reader():
        try:
            while True:
                blk = proc.stdout.read(65536)
                if not blk:
                    break
                chunks.append(blk)
        except Exception:                # noqa: BLE001 —— 读流异常交给 returncode 判定
            pass

    th = threading.Thread(target=_reader, daemon=True)
    th.start()

    def _kill():
        try:
            proc.kill()
        except Exception:                # noqa: BLE001 —— 进程可能已退出
            pass
        try:
            proc.wait(timeout=10)
        except Exception:                # noqa: BLE001
            pass

    t0 = time.time()
    timed_out = cancelled = False
    while proc.poll() is None:
        if cancel is not None and cancel.is_set():
            cancelled = True
            _kill()
            break
        if timeout and time.time() - t0 > timeout:
            timed_out = True
            _kill()
            break
        time.sleep(0.15)
    th.join(timeout=5)
    out = decode_output(b''.join(chunks))
    try:
        proc.stdout.close()
    except Exception:                    # noqa: BLE001 —— 管道已随进程回收
        pass
    dt = time.time() - t0
    if cancelled:
        raise RuntimeError('用户取消')
    if timed_out:
        raise DsrnaExeError(f'{what} 超时（>{timeout:g}s 未返回）')
    return proc.returncode, out, dt


def probe_exe(force=False):
    """定位 dsRNAmax 并探测版本号与 -seed 支持（结果缓存）。"""
    global _EXE_PROBE
    if _EXE_PROBE and not force:
        return _EXE_PROBE
    from .config import get_config
    cands = []
    try:
        p = get_config().tool('dsRNAmax')
        if p:
            cands.append(p)
    except Exception:                    # noqa: BLE001 —— 配置未就绪时退回候选路径
        pass
    cands += [os.path.join(PLATFORM_ROOT, *c.split('/')) for c in EXE_CANDIDATES]
    exe = next((c for c in cands if c and os.path.isfile(c)), None)
    if not exe:
        raise DsrnaExeError(
            '未找到 dsRNAmax 可执行文件（3rd/tools/dsrnamax/），'
            '请确认发行包完整或在设置页登记路径')
    # -h 走同一条子进程通道（不弹黑框、可用平台解码口径读回中文路径回显）
    _rc, text, _dt = _run_exe([exe, '-h'], timeout=30, what='dsRNAmax -h')
    m = re.search(r'Version:\s*([\d.]+)', text)
    _EXE_PROBE = {
        'path': exe,
        'version': m.group(1) if m else 'unknown',
        'has_seed': '-seed' in text,
    }
    return _EXE_PROBE


# ---------------------------------------------------------------- 内存守卫

def distinct_target_kmers(seqs, k=21):
    """目标集去重正向 k-mer 数（dsRNAmax 只索引正向链，口径一致）。

    seqs 接受 [seq] 或 [(name, seq)] 两种形态。
    """
    seen = set()
    acgt = set('ACGT')
    for item in seqs:
        s = (item[1] if isinstance(item, (tuple, list)) else item).upper()
        for i in range(max(0, len(s) - k + 1)):
            w = s[i:i + k]
            if set(w) <= acgt:
                seen.add(w)
    return len(seen)


def guard_memory(seqs, kmer_len, iterations, cap_gb=None):
    """按 160B × k-mer × iterations 估峰值 RSS，超预算抛 MemoryBudgetExceeded。"""
    n_kmer = distinct_target_kmers(seqs, kmer_len)
    need_gb = n_kmer * RSS_PER_KMER_PER_ITER * max(1, iterations) / 1e9 * 1.15 + 0.5
    cap = float(os.environ.get('DSRNA_RSS_CAP_GB', 0) or 0)
    avail = _mem_available_gb()
    if avail is not None:
        cap = min(cap, max(0.5, 0.7 * avail)) if cap else max(0.5, 0.7 * avail)
    if not cap:
        cap = cap_gb or DEFAULT_RSS_CAP_GB
    if need_gb > cap:
        raise MemoryBudgetExceeded(
            f'目标集含 {n_kmer:,} 个去重 {kmer_len}-mer，{iterations} 次迭代预计峰值 '
            f'内存 ≈ {need_gb:.1f} GB，超出预算 {cap:.1f} GB。请降低迭代次数、'
            f'减少/拉近目标序列，或在环境变量 DSRNA_RSS_CAP_GB 显式放宽。')
    return n_kmer, need_gb, cap


def _mem_available_gb():
    """系统可用物理内存（GB）；取不到返回 None。

    Linux 读 /proc/meminfo；Windows 走 GlobalMemoryStatusEx（ctypes，无新依赖）。
    修复 2026-10-03：原实现只读 /proc/meminfo，Windows 恒 None —— 内存守卫
    「按可用内存 70% 收紧预算」的分支在主平台（Windows）上是死代码，
    只剩固定 32GB 上限（大目标集时该收紧反而是主要的 OOM 前置保护）。
    """
    try:
        with open('/proc/meminfo') as fh:
            for line in fh:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) / 1024 / 1024
        return None                       # 能读文件但没有该字段（非 Linux）
    except OSError:
        pass
    if os.name == 'nt':
        try:
            import ctypes

            class _MEMSTATUSEX(ctypes.Structure):
                _fields_ = [('dwLength', ctypes.c_ulong),
                            ('dwMemoryLoad', ctypes.c_ulong),
                            ('ullTotalPhys', ctypes.c_ulonglong),
                            ('ullAvailPhys', ctypes.c_ulonglong),
                            ('ullTotalPageFile', ctypes.c_ulonglong),
                            ('ullAvailPageFile', ctypes.c_ulonglong),
                            ('ullTotalVirtual', ctypes.c_ulonglong),
                            ('ullAvailVirtual', ctypes.c_ulonglong),
                            ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]

            st = _MEMSTATUSEX()
            st.dwLength = ctypes.sizeof(_MEMSTATUSEX)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
                return st.ullAvailPhys / 1e9
        except Exception:                # noqa: BLE001 —— 取不到就当未知
            return None
    return None


# ---------------------------------------------------------------- FASTA 工具

def write_fasta(path, seqs):
    """写 FASTA（平台 safe_open：UTF-8 + LF 行尾 + 平台内写入校验）。"""
    with safe_open(path, 'wt') as f:
        for n, s in seqs:
            f.write(f'>{n}\n')
            for i in range(0, len(s), 70):
                f.write(s[i:i + 70] + '\n')
    return path


# ---------------------------------------------------------------- 选窗（Stage 1）

def parse_dsrnamax_table(out):
    """解析 dsRNAmax stdout 表格行 -> [cells]（tablewriter 竖线格式）。"""
    rows = []
    for line in out.splitlines():
        if line.startswith('|') and 'TARGET' not in line.upper() \
                and '---' not in line:
            cells = [c.strip() for c in line.strip('|').split('|')]
            if len(cells) >= 7:
                rows.append(cells)
    return rows


def candidate_key(rows):
    """maximin 重排键：(最差目标命中, 中位, 总和)——median 会遮蔽离群株。"""
    hits = []
    for r in rows:
        try:
            hits.append(float(r[1]))
        except (ValueError, IndexError):
            continue
    if not hits:
        return (0.0, 0.0, 0.0)
    hits.sort()
    n = len(hits)
    med = hits[n // 2] if n % 2 else (hits[n // 2 - 1] + hits[n // 2]) / 2
    return (hits[0], med, sum(hits))


def run_dsrnamax_once(exe_info, fa, construct_len, kmer_len, iterations, seed,
                      off_targets='', csv_out=None, timeout=EXEC_TIMEOUT_S,
                      cancel=None):
    """跑一次 exe，返回 (sense_arm, stdout, 耗时)。失败抛 DsrnaExeError。

    注意：exe 的 -csv 写失败时退出码仍为 0 且 stdout 不再输出序列
    （上游 output.go early-return），此处靠「解析不到臂序列」兜底报错。

    exe 与中间文件路径都可能是中文（平台根含中文，2026-10-03 实测
    dsRNAmax_det.exe 对中文 targets/-csv 路径读写正常，无需 ASCII 中转）。
    子进程通道（无黑框 / 可取消 / 输出解码）统一走 _run_exe。
    """
    exe = exe_info['path']
    cmd = [exe, '-targets', str(fa),
           '-constructLen', str(construct_len),
           '-kmerLen', str(kmer_len),
           '-otKmerLen', str(kmer_len),          # kmerLen<21 时不显式传会 fatal
           '-iterations', str(iterations)]
    if exe_info.get('has_seed'):
        cmd += ['-seed', str(seed)]              # stock exe 不认识 -seed（exit 2）
    if off_targets:
        cmd += ['-offTargets', str(off_targets)]
    if csv_out:
        cmd += ['-csv', str(csv_out)]
    p_rc, out, dt = _run_exe(cmd, timeout=timeout, cancel=cancel)
    if p_rc != 0:
        raise DsrnaExeError(f'dsRNAmax 退出码 {p_rc}: {out[-800:]}')
    m = re.search(r'dsRNA sense-arm sequence.*?\n([ACGTNacgtn]+)', out, re.S)
    if not m:
        raise DsrnaExeError('未能从 dsRNAmax 输出解析出臂序列:\n' + out[-800:])
    return m.group(1).strip().upper(), out, dt


def _unmasked_runs(seqs):
    """每条目标的最长连续 ACGT 段。"""
    out = []
    for _n, s in seqs:
        run = best = 0
        for ch in s:
            run = run + 1 if ch in 'ACGTacgt' else 0
            best = max(best, run)
        out.append(best)
    return out


def mask_targets(current, arm, kmer_len, construct_len, log=None):
    """屏蔽已选臂对应区段，供下一窗不重复。

    优先精确子串定位（单目标时共识必为子串）；嵌合共识臂不是任何目标的
    子串（s.find 失效，控制台版会静默不 mask 导致重复窗），此处回退为
    「臂 21-mer 在目标中的全部精确命中点，向后延臂长屏蔽」。
    """
    masked, exact_hit = [], False
    arm_kmers = {arm[i:i + kmer_len] for i in range(len(arm) - kmer_len + 1)}
    for n, s in current:
        i = s.find(arm)
        if i >= 0:
            exact_hit = True
            s = s[:i] + 'N' * len(arm) + s[i + len(arm):]
        else:
            hits = [j for j in range(max(0, len(s) - kmer_len + 1))
                    if s[j:j + kmer_len] in arm_kmers]
            for j in hits:
                span = min(len(s), j + max(construct_len, kmer_len))
                s = s[:j] + 'N' * (span - j) + s[span:]
        masked.append((n, s))
    if not exact_hit and hits_any(current, arm, kmer_len) and log:
        log('下一窗屏蔽：臂为嵌合共识、非目标精确子串，改用 k-mer 命中回退屏蔽'
            '（命中点向后延臂长）')
    return masked


def hits_any(current, arm, kmer_len):
    arm_kmers = {arm[i:i + kmer_len] for i in range(len(arm) - kmer_len + 1)}
    for _n, s in current:
        for j in range(max(0, len(s) - kmer_len + 1)):
            if s[j:j + kmer_len] in arm_kmers:
                return True
    return False


def select_windows(seqs, exe_info, work_dir, params, log, prog, cancel):
    """多窗 × 多候选选窗。返回 (arms, logs)。"""
    n_windows = max(1, int(params.get('n_windows') or 1))
    candidates = max(1, int(params.get('candidates') or 3))
    construct_len = int(params.get('construct_len') or 300)
    kmer_len = int(params.get('kmer_len') or 21)
    iterations = int(params.get('iterations') or 100)
    base_seed = int(params.get('seed') or 42)
    off = (params.get('off_targets') or '').strip()

    arms, logs = [], []
    current = list(seqs)
    for w in range(n_windows):
        if cancel is not None and cancel.is_set():
            raise RuntimeError('用户取消')
        runs = _unmasked_runs(current)
        avail, cap = max(runs), min(runs)
        if avail < kmer_len + 10:
            logs.append(f'第 {w + 1} 窗前停止：目标未屏蔽序列只剩 {avail} nt')
            break
        clen = construct_len
        if clen > cap:
            clen = max(30, cap)
            logs.append(f'第 {w + 1} 窗：constructLen {construct_len} → {clen}'
                        f'（最短目标仅剩 {cap} nt 未屏蔽序列）')
        fa = write_fasta(Path(work_dir) / f'targets_w{w}.fa', current)
        n_cand = candidates if len(current) > 1 else 1
        cands = []
        for c in range(n_cand):
            if cancel is not None and cancel.is_set():
                raise RuntimeError('用户取消')
            prog('select', (w * n_cand + c + 1) / (n_windows * n_cand) * 0.35,
                 f'选窗 {w + 1}/{n_windows} · 候选 {c + 1}/{n_cand}')
            seed_c = base_seed + 1000 * c        # 候选间种子域不重叠
            a, o, dt = run_dsrnamax_once(
                exe_info, fa, clen, kmer_len, iterations, seed_c, off,
                csv_out=Path(work_dir) / f'cand_w{w}_c{c}.csv', cancel=cancel)
            rows = parse_dsrnamax_table(o)
            # orig 记录候选原始编号：result.csv 兼容产物必须拷中选者自己的
            # CSV（cand_w{w}_c{orig}.csv），maximin 重排后 winner 未必是 c0
            cands.append({'orig': c, 'arm': a, 'rows': rows,
                          'key': candidate_key(rows), 'dt': dt, 'seed': seed_c})
        cands.sort(key=lambda c: c['key'], reverse=True)
        best = cands[0]
        arm = best['arm']
        # 中选候选的原生 CSV 固化为 result.csv（旧 UI/下游兼容）——
        # 按候选原始编号取，maximin 重排后 winner 未必是 c0
        best_csv = Path(work_dir) / f'cand_w{w}_c{best["orig"]}.csv'
        if best_csv.exists():
            shutil.copyfile(best_csv, Path(work_dir) / 'result.csv')
        if n_cand > 1:
            logs.append(
                f'第 {w + 1} 窗：{n_cand} 个候选按最差目标重排，保留最优 '
                + '；'.join(f"min {c['key'][0]:.0f}/med {c['key'][1]:.0f}"
                            for c in cands))
        # 盲目标 / 弱覆盖告警（median 口径的已知盲区）
        hits = []
        for r in best['rows']:
            try:
                hits.append((r[0], float(r[1])))
            except ValueError:
                continue
        blind = [n for n, h in hits if h == 0]
        if blind:
            logs.append('⚠ 警告：' + '、'.join(blind)
                        + ' 的 21-mer 命中为 0 —— 该臂无法沉默这些目标。单臂只能'
                          '覆盖共享 21-mer 的同源目标；异源目标请各自出臂后鸡尾酒混合')
        elif hits:
            best_hit = max(h for _n, h in hits)
            weak = [n for n, h in hits if 0 < h < 0.25 * best_hit]
            if weak:
                logs.append('⚠ 警告：' + '、'.join(weak)
                            + ' 覆盖不足最优目标的 25%（median 看不出来），'
                              '建议为这些目标单独补臂')
        arms.append({
            'name': f'arm_{w + 1}',
            'seq': arm,
            'length': len(arm),
            'gc': round(100 * sum(arm.count(c) for c in 'GC') / max(1, len(arm)), 1),
            'seed': best['seed'],
            'runtime_s': round(best['dt'], 1),
            'candidates_tried': n_cand,
            'target_stats': [{'target': r[0], 'kmer_matches': r[1], 'swg': r[2],
                              'gc': r[3]} for r in best['rows']],
            'candidate_keys': [{'min': c['key'][0], 'median': c['key'][1],
                                'sum': c['key'][2]} for c in cands],
        })
        logs.append(f'第 {w + 1} 窗：{len(arm)} nt · {best["dt"]:.1f}s · '
                    f'目标数 {len(best["rows"])}')
        current = mask_targets(current, arm, kmer_len, clen, log)
    return arms, logs


# ---------------------------------------------------------------- 脱靶（Stage 3/4）

def _norm_species(stem):
    """文件 stem → 物种名（剥 _refseq_rna 等常见后缀）。"""
    for suf in ('_refseq_rna', '_mrna', '_transcripts', '_cds'):
        if stem.endswith(suf):
            return stem[:-len(suf)]
    return stem


def scan_panel(arm_seqs, panel_names=None, max_mm=1, log=None, prog=None,
               cancel=None):
    """建/取面板索引并扫描臂。panel_names=None 时用面板目录全部文件。

    panel_names 支持裸种名 / 文件名 / 文件 stem，宽松匹配（双向前缀）。
    返回 (scanres, reason)；面板为空时返回 (None, reason)。
    """
    from . import dsrna_offtarget as dso
    files = dso.panel_files(str(PANEL_DIR))
    if panel_names:
        want = [str(w) for w in panel_names]

        def _match(f):
            st = Path(f).stem
            sp = _norm_species(st)
            return (f in want or Path(f).name in want or st in want
                    or sp in want
                    or any(sp.startswith(w) or w.startswith(sp) for w in want))
        files = [f for f in files if _match(f)]
    if not files:
        return None, ('面板目录为空或所选面板不存在：' + str(PANEL_DIR)
                      + '（可将 RefSeq mRNA FASTA 放入该目录后重跑）')
    indexes = {}
    for i, f in enumerate(files):
        if cancel is not None and cancel.is_set():
            raise RuntimeError('用户取消')
        stem = Path(f).stem
        sp = _norm_species(stem)
        if prog:
            prog('offtarget', 0.45 + 0.15 * i / max(1, len(files)),
                 f'索引面板 {sp}（{i + 1}/{len(files)}）')
        # 索引面板是整条链里最重的一步（真实面板 GB 级），索引内也响应取消
        info = dso.build_species_index(f, cache_dir=str(CACHE_DIR), log=log,
                                       cancel=cancel)
        # 致死表按裸种名命名（Apis_mellifera_all_lethals.tsv）
        lethal = dso.load_lethal_ids(str(LETHAL_DIR), sp)
        larr = (dso.build_lethal_index(f, lethal, cancel=cancel)
                if lethal and info['arr'].size else None)
        indexes[stem] = {'arr': info['arr'], 'lethal_arr': larr,
                         'bp': info['bp'], 'n_records': info['n_records']}
        if log:
            log(f'面板索引 {sp}: {info["n_records"]} 条转录本 / '
                f'{info["bp"]:,} bp / 去重 21-mer {info["arr"].size:,}'
                + (f' / 致死基因转录本 {len(lethal)}' if lethal else ''))
    if prog:
        prog('offtarget', 0.62, f'扫描 {len(arm_seqs)} 条臂 × {len(indexes)} 物种'
                                f'（≤{max_mm} 错配）')
    scanres = dso.scan_arms(arm_seqs, indexes, max_mm=max_mm, log=log,
                            cancel=cancel)
    return scanres, None


def safety_counts(scanres, coeff=20):
    """扫描结果 → {arm: {siRNA序号: 加权命中数}} + 汇总/红线。

    加权口径同 dsRIP：加权数 = (普通命中) + 致死命中 × coeff。
    """
    out = {}
    if not scanres:
        return out
    # mm 阈值必须显式取：max_mm=0（界面「0 错配」档）时 `or 1` 会静默放宽成
    # ≤1 错配，把 1mm 命中也算进脱靶加权（2026-10-03 修复）。
    mm_cap = scanres.get('max_mm')
    mm_cap = 1 if mm_cap is None else int(mm_cap)
    for arm, d in scanres['arms'].items():
        per = {}
        tiers_all = d.get('sirna_tiers') or {}
        n = d['n_sirna']
        all_hit = [0] * n
        lethal_hit = [0] * n
        for _sp, td in tiers_all.items():
            tiers = td.get('tiers') if isinstance(td, dict) else None
            if tiers:
                for i, t in enumerate(tiers):
                    if i < n and t is not None and t <= mm_cap:
                        all_hit[i] += 1
            lt = td.get('lethal_tiers') if isinstance(td, dict) else None
            if lt:
                for i, t in enumerate(lt):
                    if i < n and t is not None and t <= mm_cap:
                        lethal_hit[i] += 1
        for i in range(n):
            per[i + 1] = (all_hit[i] - lethal_hit[i]) + lethal_hit[i] * coeff
        out[arm] = per
    return out


def scan_summary(scanres, arm_seqs, max_mm):
    """给报告用的臂级汇总（含绝对计数与红线）。"""
    summary = {}
    for arm_name, _seq in arm_seqs:
        d = (scanres or {}).get('arms', {}).get(arm_name)
        if not d:
            continue
        per_sp = {}
        for sp, c in (d.get('per_species') or {}).items():
            per_sp[sp] = {'mm0': c.get('mm0', 0), 'mm1': c.get('mm1', 0),
                          'mm2': c.get('mm2', 0)}
        lethal = {}
        for sp, c in (d.get('lethal') or {}).items():
            lethal[sp] = {'mm0': c.get('mm0', 0), 'mm1': c.get('mm1', 0)}
        summary[arm_name] = {
            'n_sirna': d.get('n_sirna'),
            'per_species': per_sp,
            'lethal': lethal,
            'red_line': bool(d.get('red_line')),
            'lethal_rows': (d.get('lethal_rows') or [])[:50],
        }
    return summary


# ---------------------------------------------------------------- 引物（Stage 5）

def design_primers(arms, params, log=None):
    """在 best window（或全臂）内设计 T7 引物，返回每臂引物与真实扩增子。

    复用平台统一引物引擎 primer_design.design_for_sequence（含热力学评分），
    不再维护第三套 primer3 参数。

    切片坐标系（2026-10-03 修复）：best_window 的 nt_start/nt_end 是
    **guide（反义链）帧**（dsRIP 口径），臂上与之等价的区间是
    sense_nt_start..sense_nt_end（= L-nt_end+1 .. L-nt_start+1）。此前直接
    arm_seq[nt_start-1:nt_end] 切的是**镜像区间**——引物落在偏离最优窗的位置，
    且 QC 的「窗内占比」用同一组镜像坐标自比恒 ≈ 1，掩盖了错位。旧报告
    （无 sense_nt_*）按上式就地换算，行为一致。
    """
    from . import primer_design as pd
    product_min = int(params.get('primer_product_min') or 200)
    product_max = int(params.get('primer_product_max') or 300)
    num_return = int(params.get('primer_num_return') or 3)
    scope = params.get('primer_scope') or 'best_window'
    out = []
    for a in arms:
        arm_seq = a['seq'].upper()
        bw = a.get('best_window') or {}
        nt_s = bw.get('sense_nt_start')
        nt_e = bw.get('sense_nt_end')
        if (nt_s is None or nt_e is None) and bw.get('nt_start') \
                and bw.get('nt_end'):
            # guide 帧 → 臂帧（臂长 L：sense = L-nt_end+1 .. L-nt_start+1）
            _L = len(arm_seq)
            nt_s, nt_e = _L - bw['nt_end'] + 1, _L - bw['nt_start'] + 1
        if scope == 'best_window' and nt_s and nt_e and nt_e > nt_s:
            sub = arm_seq[nt_s - 1:nt_e]
            off0 = nt_s - 1
            scope_used = f'best_window(臂 nt {nt_s}-{nt_e})'
        else:
            sub = arm_seq
            off0 = 0
            scope_used = 'full_arm'
        try:
            rec = pd.design_for_sequence(
                a['name'], sub, ptype='PCR',
                core={'product_min': product_min, 'product_max': product_max,
                      'num_return': num_return}, adv={})
        except RuntimeError as e:
            # primer3 不可用（分发版未携带/环境缺失）：该臂无引物，但不让
            # 整条 dsRNA 链失败 —— 选窗/效价产物照常输出。
            entry = {
                'arm': a['name'], 'scope': scope_used, 'scope_offset': off0,
                'template_len': len(sub),
                'explain': f'T7 引物设计跳过：{e}',
                'pairs': [],
            }
            out.append(entry)
            continue
        pairs = rec.get('pairs') or []
        entry = {
            'arm': a['name'], 'scope': scope_used, 'scope_offset': off0,
            'template_len': len(sub),
            'explain': rec.get('explain', ''),
            'pairs': [],
        }
        for p in pairs[:num_return]:
            amp_s, amp_e = p.get('amp_start'), p.get('amp_end')
            amp_seq = sub[max(0, amp_s - 1):amp_e] if amp_s and amp_e else ''
            entry['pairs'].append({
                'idx': p.get('idx'),
                'f_T7': dscore.T7_PROMOTER + p['f_seq'],
                'r_T7': dscore.T7_PROMOTER + p['r_seq'],
                'f_seq': p['f_seq'], 'r_seq': p['r_seq'],
                'f_tm': p.get('f_tm'), 'r_tm': p.get('r_tm'),
                'f_gc': p.get('f_gc'), 'r_gc': p.get('r_gc'),
                'product_bp': p.get('product'),
                'amp_start_in_arm': (amp_s + off0) if amp_s else None,
                'amp_end_in_arm': (amp_e + off0) if amp_e else None,
                'penalty': p.get('penalty'),
                'score': p.get('score'),
                'recommendation': p.get('recommendation'),
                'amplicon': amp_seq,
            })
        if not pairs:
            if log:
                log(f'⚠ {a["name"]}：在 {scope_used} 内未找到满足 Tm/GC 约束的'
                    f'引物对（可放宽产物区间或改用 full_arm 范围）')
        out.append(entry)
    return out


def window_overlap(amp_s, amp_e, win_s, win_e):
    """扩增子与 best window 的重叠占比（对扩增子长度）。"""
    if not (amp_s and amp_e and win_s and win_e):
        return None
    inter = max(0, min(amp_e, win_e) - max(amp_s, win_s) + 1)
    return round(inter / max(1, amp_e - amp_s + 1), 3)


# ---------------------------------------------------------------- 隔离株覆盖（可选）

def isolate_coverage(isolates_fa, arm_seqs, k=21):
    """臂/扩增子 21-mer 对分离株池的完美覆盖（双向），0-1。

    双侧都并入反向互补集合：臂 21-mer 只要与分离株正/负链任一方向
    完全一致即算覆盖。
    """
    from . import dsrna_offtarget as dso
    if not isolates_fa or not os.path.isfile(isolates_fa):
        return None

    def _both_strand_set(seq):
        enc = dso.encode_kmers(seq.upper(), k)
        if enc.size == 0:
            return set()
        rc = dso.encode_kmers(dso.revcomp(seq.upper()), k)
        return set(int(x) for x in enc) | set(int(x) for x in rc)

    iso_sets = []
    for name, seq in dscore.read_fasta(isolates_fa):
        s = _both_strand_set(seq)
        if s:
            iso_sets.append((name, s))
    if not iso_sets:
        return None
    iso_any = set().union(*[s for _n, s in iso_sets])
    out = {}
    for name, seq in arm_seqs:
        q = set(int(x) for x in dso.encode_kmers(seq.upper(), k))
        if not q:
            continue
        covered = len(q & iso_any) / len(q)
        out[name] = {'overall': round(covered, 3),
                     'per_isolate': {n: round(len(q & s) / len(q), 3)
                                     for n, s in iso_sets}}
    return out


# ---------------------------------------------------------------- 主链

def run_chain(targets_path, params, run_dir, log, prog, cancel):
    """全链路入口。targets_path 为已就位的目标 FASTA 绝对路径。

    返回 (summary, report) 两个 dict：summary 给 UI 摘要（含旧字段兼容），
    report 为全量 JSON 报告（落 report.json）。
    """
    t0 = time.time()
    run_dir = Path(run_dir)
    out_dir = run_dir / 'dsrnamax'
    out_dir.mkdir(parents=True, exist_ok=True)

    exe_info = probe_exe()
    log(f'dsRNAmax {exe_info["version"]} · {exe_info["path"]}'
        + (' · 支持 -seed（可复现）' if exe_info['has_seed'] else
           ' · 不支持 -seed（结果不可复现，建议换用 det 构建）'))

    seqs = dscore.read_fasta(targets_path)
    if not seqs:
        raise RuntimeError('目标 FASTA 中没有序列')
    log(f'目标 {len(seqs)} 条，总 {sum(len(s) for _n, s in seqs):,} nt')

    iterations = int(params.get('iterations') or 100)
    kmer_len = int(params.get('kmer_len') or 21)
    # max_mm 显式解析（0 = 界面「0 错配」档是合法值）：`or 1` 会把 0 吞成 1，
    # 静默把 1mm 命中也计入脱靶/红线判定（2026-10-03 修复）
    _mm_raw = params.get('max_mm')
    max_mm = 1 if _mm_raw is None or str(_mm_raw).strip() == '' \
        else max(0, min(int(_mm_raw), 2))
    n_kmer, need_gb, cap = guard_memory(seqs, kmer_len, iterations)
    log(f'内存守卫：去重 {kmer_len}-mer {n_kmer:,} · 预计峰值 ≈{need_gb:.1f} GB'
        f' / 预算 {cap:.1f} GB')

    prog('select', 0.05, '开始选窗')
    arms, logs = select_windows(seqs, exe_info, out_dir, params, log, prog, cancel)
    if not arms:
        raise RuntimeError('未能选出任何臂（目标过短或构造长度超过最短目标）')
    write_fasta(out_dir / 'arms.fa', [(a['name'], a['seq']) for a in arms])

    # ---- 效价 + 安全合成评分
    prog('score', 0.62, 'siRNA 效价评分')
    coeff = int((params.get('safety_essential_coeff')
                 or dscore.DEFAULT_PARAMS['safety_essential_coeff']))
    # 面板三态：None（API 未给）= 扫全部已装面板；''/'none' = 显式跳过
    # （UI 全不勾选，提示文案承诺了"仅跳过脱靶扫描与安全评分"）；列表 = 所选
    panel_raw = params.get('panel')
    if panel_raw is None:
        panel_sel = None
    else:
        txt = str(panel_raw).strip()
        if txt.lower() in ('', 'none', 'off'):
            panel_sel = []
        else:
            panel_sel = [x.strip() for x in txt.split(',') if x.strip()] or None
    if panel_sel == []:
        scanres, why = None, '未勾选任何脱靶面板——按选择跳过脱靶扫描与安全评分'
    else:
        scanres, why = scan_panel(
            [(a['name'], a['seq']) for a in arms],
            panel_names=panel_sel, max_mm=max_mm,
            log=log, prog=prog, cancel=cancel)
    off_counts = safety_counts(scanres, coeff=coeff) if scanres else {}
    if scanres:
        log(f'脱靶扫描完成：{scanres["k"]}-mer，≤{scanres["max_mm"]} 错配，'
            f'{len(scanres["species"])} 物种')
    elif panel_sel:
        # 用户勾了面板却没扫成（名称对不上/目录空/索引坏）：必须显式告警，
        # 否则 report 里 offtarget_panel=None 会被误读成"本来就没启用"
        logs.append('⚠ 所选脱靶面板未能匹配任何面板文件，脱靶扫描与安全评分已跳过：'
                    + why)
        log(f'未启用面板脱靶扫描：{why}')
    else:
        log(f'未启用面板脱靶扫描：{why}')

    eff_p = int(params.get('efficiency_priority') or 50)
    saf_p = int(params.get('safety_priority') or 50)
    if eff_p + saf_p != 100:
        s = max(1, eff_p + saf_p)
        eff_p, saf_p = eff_p * 100 // s, 100 - eff_p * 100 // s
    scoring = {}
    for a in arms:
        r = dscore.score_sequence(a['name'], a['seq'],
                                  off_counts=off_counts.get(a['name']),
                                  efficiency_priority=eff_p,
                                  safety_priority=saf_p)
        a['best_window'] = r['best_window']
        scoring[a['name']] = r
        bw = r['best_window']
        log(f'{a["name"]}: best window nt {bw["nt_start"]}-{bw["nt_end"]} · '
            f'combined {bw["mean_combined"]:.1f} · 效价 {r["median_score"]}'
            + (f' · 原始脱靶 {bw["window_raw_offtarget_sum"]}'
               if scanres else '（无脱靶数据，纯效价选窗）'))

    # ---- 红线（绝对安全阈，不随臂内相对归一化掩盖）
    red_lines = []
    if scanres:
        for arm_name, d in scanres['arms'].items():
            if d.get('red_line'):
                rows = (d.get('lethal_rows') or [])[:5]
                where = '、'.join(f"{r['species']}(siRNA#{r['sirna_idx']} 0mm)"
                                  for r in rows if r.get('mm') == 0) or '见 lethal_rows'
                red_lines.append(f'⛔ {arm_name}：致死基因 0 错配命中（{where}）'
                                 f'——建议弃用该臂或改窗')
                log(red_lines[-1])

    # ---- 引物 + 成品 QC
    prog('primer', 0.75, 'Primer3 设计 T7 引物（限定最优窗）')
    primer_entries = design_primers(arms, params, log=log)

    prog('qc', 0.88, '成品扩增子 QC（对真实交付分子重打分）')
    # 扩增子脱靶复扫（复用已建索引，代价低）
    amp_pairs = []
    for pe in primer_entries:
        for p in pe['pairs']:
            if p.get('amplicon'):
                amp_pairs.append((f"{pe['arm']}_p{p['idx']}", p['amplicon']))
    amp_scan = None
    if scanres and amp_pairs:
        amp_scan, _ = scan_panel(amp_pairs, panel_names=panel_sel,
                                 max_mm=max_mm,
                                 log=None, prog=None, cancel=cancel)

    isolates = (params.get('isolates') or '').strip()
    iso_cov = None
    if isolates:
        prog('qc', 0.92, '分离株覆盖度')
        iso_cov = isolate_coverage(isolates, [(a['name'], a['seq']) for a in arms]
                                   + amp_pairs)
        if iso_cov:
            for aname, m in iso_cov.items():
                worst = min(m['per_isolate'].values()) if m.get('per_isolate') \
                    else m.get('overall', 1.0)
                if worst < 0.9:
                    log(f'⚠ {aname}：对分离株池最低 21-mer 覆盖仅 '
                        f'{worst * 100:.0f}%，存在逃逸风险')

    # ---- QC 汇总
    qc = []
    for pe in primer_entries:
        p0 = (pe['pairs'] or [None])[0]
        if not p0:
            qc.append({'arm': pe['arm'], 'ok': False,
                       'reason': '无可引物对'})
            continue
        amp = p0['amplicon']
        q = dscore.amplicon_qc(f'{pe["arm"]}_amplicon', amp) if amp else None
        ad = (amp_scan or {}).get('arms', {}).get(f'{pe["arm"]}_p{p0["idx"]}')
        amp_off = None
        if ad:
            tot = sum(c.get('mm0', 0) for c in (ad.get('per_species') or {}).values())
            lrow = ad.get('lethal_rows') or []
            amp_off = {'mm0_total': tot, 'red_line': bool(ad.get('red_line')),
                       'lethal_rows': lrow[:10]}
        bw = next((a['best_window'] for a in arms if a['name'] == pe['arm']), {})
        # 窗内占比按**臂帧**比较：引物按 sense_nt_* 切片（见 design_primers），
        # 这里也必须用同一帧，否则又变回"拿镜像坐标自比"。旧报告无
        # sense_nt_* 时用臂长换算（同 design_primers 的口径）。
        w_lo, w_hi = bw.get('sense_nt_start'), bw.get('sense_nt_end')
        if (w_lo is None or w_hi is None) and bw.get('nt_start') and bw.get('nt_end'):
            _L = len(next(a['seq'] for a in arms if a['name'] == pe['arm']))
            w_lo, w_hi = _L - bw['nt_end'] + 1, _L - bw['nt_start'] + 1
        ov = window_overlap(p0['amp_start_in_arm'], p0['amp_end_in_arm'],
                            w_lo, w_hi)
        qc.append({
            'arm': pe['arm'], 'ok': True,
            'pair': {kk: p0[kk] for kk in ('f_T7', 'r_T7', 'product_bp',
                                           'f_tm', 'r_tm')},
            'amplicon_len': len(amp),
            'amplicon_score_median': q and q['median_score'],
            'amplicon_best_window': q and q['best_window'],
            'amplicon_offtarget': amp_off,
            'window_overlap': ov,
        })
        log(f'QC {pe["arm"]}: 扩增子 {len(amp)} bp · 效价中位 '
            f'{q and q["median_score"]} · 窗内占比 {ov}'
            + (f' · 面板 0mm 命中 {amp_off["mm0_total"]}' if amp_off else ''))

    # ---- 落盘
    report = {
        'tool': 'dsrna',
        'exe': {'path': exe_info['path'], 'version': exe_info['version'],
                'has_seed': exe_info['has_seed']},
        'params': {k: params.get(k) for k in
                   ('construct_len', 'kmer_len', 'iterations', 'seed',
                    'candidates', 'n_windows', 'off_targets', 'panel',
                    'max_mm', 'efficiency_priority', 'safety_priority',
                    'safety_essential_coeff', 'primer_scope',
                    'primer_product_min', 'primer_product_max',
                    'primer_num_return', 'isolates')},
        'memory': {'distinct_kmers': n_kmer, 'predicted_gb': round(need_gb, 2),
                   'cap_gb': cap},
        'logs': logs,
        'arms': arms,
        'scoring': scoring,
        'offtarget_panel': None if not scanres else {
            'species': {k: {'bp': v['bp'], 'n_records': v['n_records']}
                        for k, v in scanres['species'].items()},
            'max_mm': scanres['max_mm'],
            'summary': scan_summary(scanres, [(a['name'], a['seq']) for a in arms],
                                    scanres['max_mm']),
        },
        'red_lines': red_lines,
        'primers': primer_entries,
        'amplicon_qc': qc,
        'isolate_coverage': iso_cov,
        'elapsed_s': round(time.time() - t0, 1),
    }
    # 落盘统一走平台 safe_open（编码 utf-8 + 平台根内校验），与其余工具产物一致
    with safe_open(out_dir / 'report.json', 'wt') as f:
        json.dump(report, f, ensure_ascii=False, indent=1)

    # TSV 产物（兼容管线 TSV 契约：arm/sirna_idx/species/category/gene）
    if scanres:
        lines = ['arm\tsirna_idx\tspecies\tcategory\tgene']
        for arm_name, d in scanres['arms'].items():
            for r in (d.get('lethal_rows') or []):
                lines.append(f'{arm_name}\t{r["sirna_idx"]}\t{r["species"]}\t'
                             f'{r["mm"]}mm_lethal\tlethal')
        with safe_open(out_dir / 'offtarget_lethal.tsv', 'wt') as f:
            f.write('\n'.join(lines) + '\n')

    plines = ['arm\tidx\tforward_T7\treverse_T7\tproduct_bp\tf_tm\tr_tm\t'
              'amp_start\tamp_end\tamplicon']
    for pe in primer_entries:
        for p in pe['pairs']:
            plines.append('\t'.join(str(x) for x in (
                pe['arm'], p['idx'], p['f_T7'], p['r_T7'], p['product_bp'],
                p['f_tm'], p['r_tm'], p['amp_start_in_arm'],
                p['amp_end_in_arm'], p['amplicon'])))
    with safe_open(out_dir / 'primers.tsv', 'wt') as f:
        f.write('\n'.join(plines) + '\n')
    write_fasta(out_dir / 'amplicons.fa',
                [(f"{pe['arm']}_p{p['idx']}", p['amplicon'])
                 for pe in primer_entries for p in pe['pairs'] if p['amplicon']])

    summary = {
        'tool': 'dsrna',
        'csv': str(out_dir / 'result.csv'),
        'n_targets': len(seqs),
        'construct_len': int(params.get('construct_len') or 300),
        'kmer_len': kmer_len,
        'iterations': iterations,
        'seed': int(params.get('seed') or 42),
        'off_targets': (params.get('off_targets') or '').strip() or '(无)',
        'out_dir': str(out_dir),
        # 新链路字段
        'exe_version': exe_info['version'],
        'n_arms': len(arms),
        'panel_species': (scanres and len(scanres['species'])) or 0,
        'red_lines': red_lines,
        'n_primer_pairs': sum(len(pe['pairs']) for pe in primer_entries),
        'report': 'dsrnamax/report.json',
        'primers_tsv': 'dsrnamax/primers.tsv',
        'arms_fa': 'dsrnamax/arms.fa',
        'elapsed_s': report['elapsed_s'],
    }
    # exe 原生 CSV 若存在则保留路径（候选循环以最后一次为准，仅作原始存档）
    return summary, report
