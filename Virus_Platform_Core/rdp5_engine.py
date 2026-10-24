# -*- coding: utf-8 -*-
"""RDP5 引擎集成——原生 Windows 命令行版（RDP5CL.exe）。

RDP5 本体是 Windows 程序（服务器上经 Wine 运行 /root/RDP5/RDP5CL.exe）；
本平台即 Windows，直接原生调用，无需 Wine。发行版已拉取到
3rd/tools/rdp5/（含 3seqTable、clustalw2、FastTree 等伴生工具）。

CLI 口径（对齐服务器 /opt/rdp5_run_39.sh）：
    RDP5CL.exe -f <input.fasta> -ofp <prefix> -ds -rbdp
输入文件须与 RDP5CL 同目录运行（服务器即此做法），产物落在当前目录：
    <prefix>.csv                                  主结果表
    <prefix> breakpoint distribution.csv          断点分布
    <prefix>*BreakpointPositions*.csv             断点位置
事件解析器移植自 virome_phylo_pipeline/utils/rdp5_validate.py:parse_rdp5_csv
（RDP5 主表列序：event,markers,bp_start,bp_end,…,recombinant,minor,major,
  RDP,GENECONV,Bootscan,Maxchi,Chimaera,SiScan,PhylPro,LARD,3Seq）。

mask（吸收 virome_phylo_pipeline/mask_recombination.py）：
检出事件的重叠区间合并后，比对中对应列掩蔽为 N（可进 BEAST）——
产出「干净比对」直接送建树，即参考管线「重组必须在建树前」的语义。

⚠️ 已知怪癖（2026-10-12 实测定案）：-rbdp 的断点分布分析只统计
**双亲本都已鉴定** 的事件——主表事件只要带 Unknown 亲本（真实病毒
数据极常见），`<prefix> breakpoint distribution.csv` 会静默输出全 0、
`…BreakpointPositions.csv` 为空或缺失，前端照画就是一条 0 平线。
`breakpoint_distribution_degenerate()` + `synthesize_breakpoint_distribution()`
即为该情况兜底：用事件断点重建 200nt 滑窗分布，格式与 RDP5 原产 CSV 兼容。
"""

import os
import queue
import re
import shutil
import subprocess
import threading
import time
import uuid

from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.utils import decode_output, to_utf8_file

PLATFORM_TOOLS = ('3rd/tools/rdp5',)

# ⚠️ 模块级锁：RDP5CL 固定把产物写到**自己所在的目录**（bin_dir），而平台可以
# 并发发起多个 RDP5 任务 → 两个调用会互相删除/覆盖对方的输入与产物
# （2026-10-03 修复，双保险之一；另一个是每次调用的独立前缀）。进程内串行化
# 代价可接受：RDP5CL 本身是单机重任务，并发跑两个只会互相拖慢。
_RDP5_LOCK = threading.Lock()

# RDP5 主表九方法 p 值列（紧跟 recombinant/minor/major 之后）
RDP5_METHODS = ['RDP', 'GENECONV', 'Bootscan', 'Maxchi', 'Chimaera',
                'SiScan', 'PhylPro', 'LARD', '3Seq']


def find_rdp5cl():
    """定位 RDP5CL.exe（3rd/tools/rdp5/ 或 platform.json tools 登记）。"""
    from Virus_Platform_Core.config import get_config
    try:
        p = get_config().tool('RDP5CL')
        if p and os.path.isfile(p):
            return p
    except Exception:
        pass
    for cand in PLATFORM_TOOLS:
        p = os.path.join(PLATFORM_ROOT, *cand.split('/'), 'RDP5CL.exe')
        if os.path.isfile(p):
            return p
    return None


def _cleanup_rdp5_products(bin_dir, run_prefix):
    """删除本次调用在 RDP5 目录留下的文件（临时输入 + 半成品产物）。

    只匹配本次调用的完整前缀（含 uuid），不会碰到别的任务的文件。
    """
    try:
        names = os.listdir(bin_dir)
    except OSError:
        return
    for fn in names:
        if fn.startswith(run_prefix):
            try:
                os.remove(os.path.join(bin_dir, fn))
            except OSError:
                pass


def run_rdp5cl(fasta, out_dir, prefix='rdp5', timeout=3600, logger=None,
               cancel=None):
    """运行 RDP5CL。返回主结果 CSV 路径（失败抛 RuntimeError）。

    做法对齐服务器封装：把输入复制进 RDP5 目录、cd 到该目录运行
    （RDP5CL 把产物写在当前目录），结束后把 <prefix>* 产物拷回 out_dir。

    ⚠️ 并发安全（2026-10-03 修复）：
      · 进程内用模块级 `_RDP5_LOCK` 串行化（RDP5CL 的产物落在共享目录，
        并发调用会互删/互覆）；
      · 每次调用用 **`rdp5_<uuid4 前 8 位>_` 独立前缀**（输入与产物名全带它），
        回拷只搬本前缀的文件 —— 即使锁被绕过（多进程并发）也不会串件；
      · 拷回 out_dir 时去掉 uuid 段，产物名仍是调用方给的 `<prefix>*`
        （如 `rdp5.csv`），下游命名契约不变。
    """
    with _RDP5_LOCK:
        return _run_rdp5cl_impl(fasta, out_dir, prefix=prefix, timeout=timeout,
                                logger=logger, cancel=cancel)


def _run_rdp5cl_impl(fasta, out_dir, prefix='rdp5', timeout=3600,
                     logger=None, cancel=None):
    exe = find_rdp5cl()
    if not exe:
        raise RuntimeError('未找到 RDP5CL.exe（3rd/tools/rdp5/）')
    bin_dir = os.path.dirname(exe)
    os.makedirs(out_dir, exist_ok=True)
    # 每次调用独立前缀（uuid4 前 8 位）：并发/多进程下输入与产物互不覆盖
    run_prefix = f'{prefix}_{uuid.uuid4().hex[:8]}_'
    input_name = f'{run_prefix}input.fasta'
    local_fa = os.path.join(bin_dir, input_name)
    shutil.copyfile(fasta, local_fa)

    cmd = [exe, '-f', input_name, '-ofp', run_prefix, '-ds', '-rbdp']
    if logger:
        logger.log('$ ' + ' '.join(cmd) + f'   (cwd={bin_dir})')
    # 按字节迭代再 decode_output（utf-8→GBK 兜底）：RDP5CL 把含中文路径的
    # 提示行按系统 ANSI 打印，text=True+utf-8 会在捕获时打成 U+FFFD（不可逆）
    proc = subprocess.Popen(cmd, cwd=bin_dir, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    deadline = time.time() + timeout
    # ⚠️ cancel/超时不能只在"读到新输出行"时检查（2026-10-03 修复）：
    #    RDP5CL 可能长时间**静默**（大比对 / 捆绑工具计算中），旧实现下
    #    cancel 与 timeout 双双失效。改为：stdout 交给独立读取线程（Windows
    #    下 select 不能用于管道），主线程按 0.5s 轮询排队，静默期间也照常
    #    检查 cancel/超时 —— 命中立刻 kill 并报错。
    #    （不依赖 kill 关闭管道：RDP5CL 派生的孙进程可能仍攥着 stdout 句柄。）
    lines = queue.Queue()
    _EOF = object()

    def _drain():
        try:
            for ln in proc.stdout:
                lines.put(ln)
        except Exception:
            pass
        finally:
            lines.put(_EOF)

    reader = threading.Thread(target=_drain, daemon=True, name='rdp5-stdout')
    reader.start()
    while True:
        try:
            line = lines.get(timeout=0.5)
        except queue.Empty:
            line = None
        if line is _EOF:
            break
        if line is not None and logger:
            logger.log(decode_output(line).rstrip('\n'))
        if cancel is not None and cancel.is_set():
            proc.kill()
            _cleanup_rdp5_products(bin_dir, run_prefix)
            raise RuntimeError('用户取消')
        if time.time() > deadline:
            proc.kill()
            _cleanup_rdp5_products(bin_dir, run_prefix)
            raise RuntimeError(f'RDP5 超时（>{timeout}s）')
        if line is None and proc.poll() is not None:
            # 进程已退出但管道没关（孙进程还攥着句柄）：把余下行读净后收尾，
            # 不再空等 EOF（否则会把正常退出的任务拖到超时）。
            time.sleep(0.1)
            while True:
                try:
                    rest = lines.get_nowait()
                except queue.Empty:
                    break
                if rest is _EOF:
                    break
                if logger:
                    logger.log(decode_output(rest).rstrip('\n'))
            break
    rc = proc.wait()
    reader.join(timeout=2)
    # 产物拷回（文本产物顺带就地规范为 UTF-8：RDP5 按系统 ANSI 自写，
    # 与平台其余产物统一；二进制/纯 ASCII 文件 to_utf8_file 原样跳过）。
    # 只搬**本次前缀**的文件；拷回时去掉 uuid 段，保持 `<prefix>*` 命名契约。
    for fn in os.listdir(bin_dir):
        if fn.startswith(run_prefix) and fn != input_name:
            src = os.path.join(bin_dir, fn)
            if os.path.isfile(src):
                tgt_name = prefix + fn[len(run_prefix):]
                tgt = os.path.join(out_dir, tgt_name)
                shutil.copyfile(src, tgt)
                if tgt_name.lower().endswith(('.csv', '.txt')):
                    to_utf8_file(tgt)
                try:
                    os.remove(src)
                except OSError:
                    pass
    try:
        os.remove(local_fa)
    except OSError:
        pass
    if rc != 0:
        raise RuntimeError(f'RDP5CL 退出码 {rc}')
    csv_path = os.path.join(out_dir, f'{prefix}.csv')
    if not os.path.isfile(csv_path):
        raise RuntimeError('RDP5CL 未产出结果 CSV（可能无事件或输入格式问题）')
    return csv_path


def parse_rdp5_csv(csv_path):
    """解析 RDP5 主结果表 → 事件列表（移植 parse_rdp5_csv）。

    行格式（逗号分隔）：event#,markers,bp_start,bp_end,…,
    recombinant,minor_parent,major_parent,九方法 p 值（NS=不显著）。
    """
    events = []
    with open(csv_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('Table key', '~', '*', '$', '^',
                                            'Minor', 'Major', 'Unknown', 'NS')):
                continue
            m = re.match(r'\s*(\d+)\s*,\s*(\d+)([~\^\*\$]*)\s*,\s*(\d+)\s*,\s*(\d+)',
                         line)
            if not m:
                continue
            fields = [x.strip() for x in line.split(',')]
            if len(fields) <= 10:
                continue
            rec_id = fields[8].split('\n')[0].strip().rstrip(',').lstrip('^')

            def get_parent(field):
                if field.startswith('Unknown'):
                    m2 = re.search(r'\(([^)]+)\)', field)
                    return m2.group(1) if m2 else ''
                return field.split('\n')[0].strip().rstrip(',')

            minor_id = get_parent(fields[9])
            major_id = get_parent(fields[10])
            pvals = {}
            for i, mn in enumerate(RDP5_METHODS):
                j = 11 + i
                if j < len(fields) and fields[j] and fields[j] != 'NS':
                    try:
                        pvals[mn] = float(fields[j])
                    except ValueError:
                        pass
            events.append({
                'num': int(m.group(1)), 'markers': m.group(3),
                'bp_start': int(m.group(4)), 'bp_end': int(m.group(5)),
                'recombinant': rec_id, 'minor_parent': minor_id,
                'major_parent': major_id,
                'methods': pvals, 'n_methods': len(pvals),
            })
    return events


def merge_intervals(intervals, aln_len):
    """合并区间（1-based，含跨末端环状拆分）。移植 mask_recombination。"""
    segs = []
    for b, e in intervals:
        b, e = max(1, b), min(aln_len, e)
        if b <= e:
            segs.append((b, e))
        else:
            segs.append((b, aln_len))
            if e >= 1:
                segs.append((1, e))
    segs.sort()
    merged = []
    for b, e in segs:
        if merged and b <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((b, e))
    return merged


def mask_alignment(fasta_path, intervals, mask_char='N'):
    """把比对中给定区间（1-based）掩蔽为 mask_char。

    返回 (records [(name, seq)], n_masked_cols)。
    移植 mask_recombination.mask_alignment（不依赖 dnasp）。
    """
    names, seqs = [], []
    name, buf = None, []
    with open(fasta_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    seqs.append(''.join(buf))
                    names.append(name)
                name = line[1:].split()[0] if line[1:].split() else line[1:]
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        seqs.append(''.join(buf))
        names.append(name)
    if not seqs:
        return [], 0
    rows = [list(s) for s in seqs]
    n_masked = 0
    for b, e in intervals:
        for i in range(b - 1, min(e, len(rows[0]))):
            for r in rows:
                if r[i] not in ('-', 'N'):
                    r[i] = mask_char
            n_masked += 1
    records = list(zip(names, (''.join(r) for r in rows)))
    return records, n_masked


def _names_from(fasta_path):
    out = []
    with open(fasta_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                h = line[1:].strip()
                out.append(h.split()[0] if h.split() else h)
    return out


def mask_from_events(fasta_path, events, mask_char='N'):
    """按事件表掩蔽：合并全部 (bp_start, bp_end) → 掩蔽比对。

    返回 (records [(name, seq)], n_masked_cols, merged_intervals)。
    """
    # 比对长度 = **单条**序列长度（读第一条序列，遇到下一个 '>' 即停）。
    # 早先写成"累加所有序列的长度"，会把 3×100bp 的比对算成 300——
    # 环状区间（bp_start > bp_end）走 merge_intervals 的跨末端拆分时，
    # 就会给出 (90, 300) 这种越界区间。
    aln_len, seen_header = 0, False
    with open(fasta_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if seen_header:
                    break
                seen_header = True
                continue
            if seen_header:
                aln_len += len(line.strip())
    # 兼容两种事件键：RDP5 的 bp_start/bp_end 与本地引擎的 start/end
    intervals = []
    for e in events:
        b = e.get('bp_start', e.get('start'))
        en = e.get('bp_end', e.get('end'))
        if b is None or en is None:
            continue
        intervals.append((int(b), int(en)))
    merged = merge_intervals(intervals, aln_len)
    records, n_masked = mask_alignment(fasta_path, merged, mask_char)
    return records, n_masked, merged


# ── 断点分布兜底（2026-10-12）────────────────────────────────────
# 背景：RDP5CL 的 -rbdp 只统计双亲本都已鉴定的事件（见模块 docstring）。
# 实测复现：10 条 CMV RNA3 真实序列 → 2 个事件（各带一个 Unknown 亲本）
# → 分布 CSV 1146 行全 0；示例比对删到只剩 3 条序列（迫使次亲本 Unknown）
# 同样复现；而双亲本齐全的运行分布正常——因果闭环。
BDP_CSV = 'rdp5 breakpoint distribution.csv'
BDP_HEADER = ('Position in alignment, '
              'Recombination breakpoint number (200nt win)')


def breakpoint_distribution_degenerate(csv_path):
    """判定 RDP5 断点分布 CSV 是否为"静默全 0"退化产物。

    True = 文件缺失、无数据行、或 RBN 列（第 2 列）全 0。
    RDP5 正常产物在没有断点的位置 RBN 也接近 0，但 CI 列 >0；此处只看
    RBN 列，且只有"事件存在"时上层才会调它，误判面足够小。
    """
    if not os.path.isfile(csv_path):
        return True
    saw_row = False
    with open(csv_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit():
                continue        # 表头 / 空行
            parts = line.split(',')
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            saw_row = True
            if v != 0:
                return False
    return True


def synthesize_breakpoint_distribution(csv_path, events, aln_len=None,
                                       window=200):
    """用检出事件的断点重建 200nt 滑窗断点分布 CSV（RDP5 退化兜底）。

    每个 200nt 窗（起点 p，覆盖 [p, p+200)）计落入的事件断点数；环状
    事件（start > end）的两个断点直接按位置计，不跨末端绕回。位置网格
    步长 2 对齐 RDP5 原产 CSV。aln_len 缺省时按最靠后的断点补足。
    返回写入的数据行数。
    """
    bps = []
    for e in events:
        for k in ('bp_start', 'start', 'bp_end', 'end'):
            v = e.get(k)
            if v is not None:
                bps.append(int(v))
    if not bps:
        return 0
    if not aln_len:
        aln_len = max(bps) + window
    with open(csv_path, 'w', encoding='utf-8', newline='') as f:
        f.write(BDP_HEADER + '\n')
        n_rows = 0
        for p in range(1, aln_len + 1, 2):
            cnt = sum(1 for b in bps if p <= b < p + window)
            f.write(f'{p},{cnt}\n')
            n_rows += 1
    return n_rows


def orf_csv_to_aln_map(orf_csv, aln_fasta, ref_id=None):
    """基因坐标 CSV（Rgrapher ORFCoords 风格）→ 比对坐标基因图谱。

    CSV 表头按 Rgrapher 惯例（Gene.symbol / Start / Stop，大小写不敏感）；
    坐标是**参考序列**的核苷酸坐标（1-based，GenBank 口径）。参考行取
    ref_id 前缀匹配的序列，缺省用比对第一条。参考行里的 gap 列自动跨过
    （与 RDP5 断点坐标同一坐标系）。越界/畸形行跳过并计数。
    返回 {'genes': [{'name','start','end'}, …], 'skipped': int, 'ref': str}。
    """
    import csv as _csv
    name_col, start_col, stop_col = None, None, None
    genes, skipped = [], 0
    with open(orf_csv, encoding='utf-8-sig', errors='replace') as f:
        rows = list(_csv.reader(f))
    if not rows:
        raise ValueError('基因坐标 CSV 为空')
    hdr = [h.strip().lower() for h in rows[0]]
    for i, h in enumerate(hdr):
        if h in ('gene.symbol', 'gene_symbol', 'gene', 'symbol', 'name'):
            name_col = i
        elif h in ('start', 'begin', 'from'):
            start_col = i
        elif h in ('stop', 'end', 'to'):
            stop_col = i
    for row in (rows[1:] if name_col is not None else rows):
        try:
            if name_col is not None:
                gname = str(row[name_col]).strip()
                start, stop = int(row[start_col]), int(row[stop_col])
            else:
                nums = [c for c in row if c.strip() and
                        c.strip().lstrip('.').isdigit()]
                texts = [c for c in row if c.strip() and
                         not c.strip().lstrip('.').isdigit()]
                if len(nums) < 2 or not texts:
                    raise ValueError
                gname = texts[0].strip()
                start, stop = int(nums[0]), int(nums[1])
        except (ValueError, IndexError):
            skipped += 1
            continue
        if not gname or start < 1 or stop < start:
            skipped += 1
            continue
        genes.append({'name': gname, 'ref_start': start, 'ref_stop': stop})
    if not genes:
        if skipped:
            raise ValueError(f'基因坐标 CSV 无有效行（{skipped} 行畸形已跳过）')
        raise ValueError('基因坐标 CSV 无有效行')

    # 参考行 → nt→col 映射
    cur_name, buf = None, []
    rows_fa = []
    with open(aln_fasta, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_name is not None:
                    rows_fa.append((cur_name, ''.join(buf)))
                cur_name = line[1:].split()[0] if line[1:].split() else line[1:]
                buf = []
            elif line:
                buf.append(line)
    if cur_name is not None:
        rows_fa.append((cur_name, ''.join(buf)))
    if not rows_fa:
        raise ValueError('比对 FASTA 为空，无法映射基因坐标')
    ref_name, ref_seq = rows_fa[0]
    for nm, sq in rows_fa:
        if ref_id and nm.startswith(ref_id):
            ref_name, ref_seq = nm, sq
            break
    nt2col, n = {}, 0
    for i, ch in enumerate(ref_seq, 1):
        if ch != '-':
            n += 1
            nt2col[n] = i
    mapped, out = 0, []
    for g in genes:
        s, e = nt2col.get(g['ref_start']), nt2col.get(g['ref_stop'])
        if s is None or e is None:
            skipped += 1
            continue
        out.append({'name': g['name'], 'start': min(s, e), 'end': max(s, e)})
        mapped += 1
    if not out:
        raise ValueError(f'基因坐标全部越界（参考序列 {ref_name} 共 {n} nt，'
                         f'{skipped} 行无法映射）——坐标应按参考序列 nt 给')
    return {'genes': out, 'skipped': skipped, 'ref': ref_name,
            'ref_nt_count': n}


