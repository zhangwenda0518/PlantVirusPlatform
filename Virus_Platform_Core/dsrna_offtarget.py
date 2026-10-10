# -*- coding: utf-8 -*-
"""
dsRNA 臂 siRNA 宿主脱靶扫描（k-mer 排序 uint64 数组 + 二分，纯 numpy 实现）。

移植自已验证脚本 C-host_classify/plant_virus_db_pipeline/dsRNAmax/probe_out/
host_cds_offtarget.py，语义与其保持一致：
  - 2-bit 编码（A=0 C=1 G=2 T=3），滚动编码，跳过含非 ACGT 的窗口；
  - 索引双链入集（seq 与 revcomp(seq) 都 encode），查询只走正向；
  - 每条 siRNA 按 mm0 -> mm1 -> mm2 分层顺序查（best_tier 语义）；
  - np.searchsorted 二分查找（纯 Python set 存 ~1e8 级 k-mer 要 ~15GB 指针，
    排序数组 + 二分是参考实现在 246 服务器上验证过的方案）。

红线（均来自真实 bug，改动前先读对应注释）：
  1. mm 分层必须 0->1->2 顺序查：先查 1mm 会把精确命中记成 1mm
     （blastn 曾抓到一条 21/21 精确匹配被本扫描记进 1mm 档）；
  2. revcomp 必须 complement AND reverse：只 complement 的 k-mer 集合不是
     负链 k-mer 集合，会漏掉所有负链命中（同样是 blastn 抓出来的）；
  3. 索引必须双链入集，查询单链即可（负链命中无需二次查询）；
  4. searchsorted 返回的 idx 必须 clip 到 [0, size-1] 再比较
     （查询值大于数组最大值时 idx == size，直接下标会越界/误判）；
  5. 2-bit 编码顺序固定 A=0 C=1 G=2 T=3（与参考实现 CODE 表一致）。

仅依赖 numpy + 标准库（无 biopython）。Python 3.12。
"""
import json
import random
import sys
from pathlib import Path

import numpy as np

K_DEFAULT = 21

# lethal_rows 逐条记录的截断上限（防止极端臂把结果撑爆）
_LETHAL_ROW_CAP = 500

# 2-bit 编码表：A=0 C=1 G=2 T=3（红线 5），大小写都认，其余字节=255（坏窗口）
_CODE = np.full(256, 255, dtype=np.uint8)
for _b, _c in zip(b"ACGT", range(4), strict=True):
    _CODE[_b] = _c
    _CODE[_b + 32] = _c

# 反向互补表：bytes 用 _RC_BYTES，str 用 _RC_STR（都同时处理大小写）
_RC_BYTES = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
_RC_STR = str.maketrans("ACGTacgt", "TGCAtgca")


class PanelError(Exception):
    """panel/索引层面的配置或数据错误（文件缺失、FASTA 无记录、索引结构非法）。"""


def revcomp(seq: str) -> str:
    """反向互补：complement AND reverse（红线 2）。

    seq.translate(RC) 单独使用只是 complement，其 k-mer 集合并不等于负链的
    k-mer 集合——这个错误曾让所有真实负链命中全部漏检（blastn 抓到一条
    21/21 匹配是本函数漏掉的）。兼容 bytes 输入（内部 FASTA 解析用），bytes
    进 bytes 出。注意字母表是 DNA（ACGTacgt）：含 U 的 RNA 串不会被转换，
    本模块所有调用方（臂/扩增子/面板 FASTA）约定为 DNA 编码。
    """
    if isinstance(seq, bytes):
        return seq.translate(_RC_BYTES)[::-1]
    return seq.translate(_RC_STR)[::-1]


def encode_kmers(seq, k: int = K_DEFAULT) -> np.ndarray:
    """滚动 2-bit 编码序列的全部 k-mer，跳过含非 ACGT 的窗口。

    输入 bytes 或 str（str 里非 ASCII 字符按非 ACGT 处理，不会错位）。
    返回 uint64 一维数组（基因组顺序，未排序）；n < k 时返回空数组。
    """
    if not 1 <= k <= 32:
        raise ValueError(f"k 必须在 1..32（uint64 2-bit 编码上限），收到 k={k}")
    if isinstance(seq, str):
        # errors='replace' 保长度：非 ASCII 字符替换成 '?'，窗口整体作废
        seq = seq.upper().encode("ascii", "replace")
    c = _CODE[np.frombuffer(seq, dtype=np.uint8)]
    n = c.size
    if n < k:
        return np.empty(0, dtype=np.uint64)
    m = n - k + 1
    ok = np.ones(m, dtype=bool)
    acc = np.zeros(m, dtype=np.uint64)
    for j in range(k):
        w = c[j:j + m]
        ok &= w != 255                      # 窗口内含非 ACGT -> 整窗丢弃
        acc = (acc << np.uint64(2)) | w.astype(np.uint64)
    return acc[ok]


# ---------------- 错配变体表（与参考实现 build_variant_table 同思路） ----------------

_VARIANT_CACHE: dict = {}


def _variant_tables(k: int):
    """k-mer 的 1mm / 2mm XOR 掩码表（按 k 缓存）。

    每个位置 12 个掩码（a^b, a!=b，共 4*3），对任一具体 k-mer 每位置恰有 3 个
    掩码非平凡 -> 252 个掩码恰好覆盖全部 63 个单替换邻居；双掩码同理覆盖全部
    C(k,2)*9 = 1890 个双替换邻居。np.unique 去重后与逐一枚举 63/1890 个变体
    完全等价，但一次位运算生成、无需按碱基分支。
    """
    tables = _VARIANT_CACHE.get(k)
    if tables is not None:
        return tables
    single, double = [], []
    for i in range(k):
        sh = 2 * (k - 1 - i)
        for a in range(4):
            for b in range(4):
                if a != b:
                    single.append((a ^ b) << sh)
    for i in range(k):
        for j in range(i + 1, k):
            shi, shj = 2 * (k - 1 - i), 2 * (k - 1 - j)
            for a in range(4):
                for b in range(4):
                    if a == b:
                        continue
                    for c in range(4):
                        for d in range(4):
                            if c == d:
                                continue
                            double.append(((a ^ b) << shi) | ((c ^ d) << shj))
    tables = (np.array(single, dtype=np.uint64), np.array(double, dtype=np.uint64))
    _VARIANT_CACHE[k] = tables
    return tables


def _has_hit(kmer: int, host: np.ndarray, masks: np.ndarray) -> bool:
    """kmer 的全部掩码变体在排序数组 host 里是否有命中（参考实现 has_hit）。

    红线 4：searchsorted 之后必须 clip 到数组范围内再比较——查询值大于数组
    最大值时 idx == host.size，直接 host[idx] 会越界。
    """
    if host.size == 0:
        return False
    q = np.unique(np.uint64(kmer) ^ masks)
    idx = np.searchsorted(host, q)
    np.clip(idx, 0, host.size - 1, out=idx)
    return bool(np.any(host[idx] == q))


def _best_tier(kmer: int, host: np.ndarray, k: int, max_mm: int):
    """kmer 对宿主库的最佳命中档：0/1/2 = 最小错配数，None = max_mm 内无命中。

    红线 1：必须 0 -> 1 -> 2 顺序查。先查 1mm 会把精确命中静默记成 1mm
    （参考实现实测踩过：blastn 报 21/21 精确匹配，扫描却归到 1mm 档）。
    """
    # 精确命中用 searchsorted 而不是 `in`：`in` 在上亿元素数组上是逐个线性扫描
    idx = np.searchsorted(host, np.uint64(kmer))
    if idx < host.size and host[idx] == np.uint64(kmer):
        return 0
    single, double = _variant_tables(k)
    if max_mm >= 1 and _has_hit(kmer, host, single):
        return 1
    if max_mm >= 2 and _has_hit(kmer, host, double):
        return 2
    return None


# ---------------- 内部小工具 ----------------

def _fasta_iter(path):
    """极简 FASTA 迭代器，yield (header 首词, 序列 bytes)。只依赖标准库。"""
    name, buf = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(buf)
                parts = line[1:].split()
                name = parts[0].decode("ascii", "replace") if parts else ""
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, b"".join(buf)


def _strip_version(name: str) -> str:
    """去掉 RefSeq 风格版本号后缀：XM_012345.1 -> XM_012345。"""
    i = name.rfind(".")
    if i > 0 and name[i + 1:].isdigit():
        return name[:i]
    return name


def _log(log, msg: str) -> None:
    """兼容平台 logger 对象（.log(msg[, level])）与普通可调用；失败静默。"""
    if log is None:
        return
    try:
        if callable(log):
            log(msg)
        else:
            try:
                log.log(msg)
            except TypeError:
                log.log(msg, "INFO")
    except Exception:
        pass


# ---------------- 对外 API ----------------

_PANEL_EXTS = {".fa", ".fna", ".fasta"}


def panel_files(panel_dir) -> list:
    """列出 panel_dir 下的 *.fa / *.fna / *.fasta（按路径名排序）。

    目录不存在返回 []（调用方自行决定是否报错）。
    """
    d = Path(panel_dir)
    if not d.is_dir():
        return []
    return sorted(str(p) for p in d.iterdir()
                  if p.is_file() and p.suffix.lower() in _PANEL_EXTS)


def build_species_index(fa_path, k: int = K_DEFAULT, cache_dir=None, log=None,
                        cancel=None) -> dict:
    """对单个物种 FASTA 建脱靶索引：双链全部 k-mer -> 排序去重 uint64 数组。

    双链入集是红线 3：seq 与 revcomp(seq) 都 encode 后合并去重，之后查询单链
    正向即可覆盖正负两条链的命中。

    cache_dir 给定时做磁盘缓存：
      <stem>.k{k}.npy            排序去重后的 uint64 数组
      <stem>.k{k}.json           元数据（源文件 mtime_ns / size、k、数组长度
                                 n、总 bp、记录数）——任一不匹配即重建

    cancel 给定时逐记录检查（threading.Event 风格），置位抛 RuntimeError：
    真实面板是 GB 级，重建索引可能要几分钟，必须可中断（2026-10-03）。

    返回 {'arr': np.ndarray, 'bp': int, 'n_records': int, 'cached': bool}。
    源文件不存在、或 FASTA 里没有任何记录时抛 PanelError。
    """
    fa = Path(fa_path)
    if not fa.is_file():
        raise PanelError(f"FASTA 不存在: {fa}")

    npy = jsn = None
    if cache_dir is not None:
        cdir = Path(cache_dir)
        cdir.mkdir(parents=True, exist_ok=True)
        npy = cdir / f"{fa.stem}.k{k}.npy"
        jsn = cdir / f"{fa.stem}.k{k}.json"

    # --- 尝试读缓存：元数据（mtime_ns/size/k/n/bp）与实际数组逐一核对 ---
    if npy is not None and npy.is_file() and jsn.is_file():
        try:
            meta = json.loads(jsn.read_text(encoding="utf-8"))
            st = fa.stat()
            if (meta.get("k") == k
                    and meta.get("mtime_ns") == st.st_mtime_ns
                    and meta.get("size") == st.st_size
                    and isinstance(meta.get("bp"), int)
                    and isinstance(meta.get("n_records"), int)):
                arr = np.load(npy, allow_pickle=False)
                if arr.dtype == np.uint64 and arr.size == meta.get("n"):
                    _log(log, f"[dsrna_offtarget] 缓存命中: {fa.name} "
                              f"k={k} n={arr.size:,}")
                    return {"arr": arr,
                            "bp": int(meta["bp"]),
                            "n_records": int(meta["n_records"]),
                            "cached": True}
        except (OSError, ValueError, json.JSONDecodeError):
            pass                            # 缓存损坏/元数据不匹配 -> 走重建

    # --- 重建：双链全部 k-mer 合并去重（红线 3） ---
    arrs, bp, n_records = [], 0, 0
    for _name, seq in _fasta_iter(fa):
        if cancel is not None and cancel.is_set():
            raise RuntimeError('用户取消')
        n_records += 1
        bp += len(seq)
        for strand in (seq, revcomp(seq)):  # 红线 3：双链入集
            a = encode_kmers(strand, k)
            if a.size:
                arrs.append(a)
    if n_records == 0:
        raise PanelError(f"FASTA 中没有任何记录: {fa}")
    arr = (np.unique(np.concatenate(arrs)) if arrs
           else np.empty(0, dtype=np.uint64))

    if npy is not None:
        np.save(npy, arr)
        jsn.write_text(json.dumps({
            "source": fa.name, "k": int(k),
            "mtime_ns": fa.stat().st_mtime_ns, "size": fa.stat().st_size,
            "n": int(arr.size), "bp": int(bp), "n_records": int(n_records),
        }, ensure_ascii=False, indent=1), encoding="utf-8")

    _log(log, f"[dsrna_offtarget] 建索引: {fa.name} {n_records} 条记录 / "
              f"{bp:,} bp -> {arr.size:,} 个去重 {k}-mer（双链）")
    return {"arr": arr, "bp": int(bp), "n_records": int(n_records),
            "cached": False}


def load_lethal_ids(lethal_dir, species_stem: str) -> set:
    """读 <species_stem>_all_lethals.tsv 的 gene_id 集合（首行表头，第一列）。

    文件不存在返回空 set。
    """
    p = Path(lethal_dir) / f"{species_stem}_all_lethals.tsv"
    if not p.is_file():
        return set()
    ids = set()
    with open(p, encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\r\n")
            if not line or i == 0:          # 首行是表头
                continue
            gene_id = line.split("\t")[0].strip()
            if gene_id:
                ids.add(gene_id)
    return ids


# 致死子索引进程内缓存：键含源文件 mtime 与 lethal_ids 全集，同任务内
# scan_panel（臂扫描 + 扩增子复扫）调用两次时免整文件重读（Apis 88MB）。
_LETHAL_INDEX_CACHE: dict = {}


def build_lethal_index(fa_path, lethal_ids, k: int = K_DEFAULT,
                       cancel=None) -> np.ndarray:
    """只对 lethal_ids 命中的转录本建同款排序数组（双链，红线 3）。

    命中判定：record 名（header 第一词）本身在 lethal_ids 中，或其去版本号
    形式在 lethal_ids 中。RefSeq 记录名形如 XM_012345.1，而 lethal 列表里的
    id 可能带或不带版本号——record 名与去版本号名字各查一次，两种来源都能
    匹配上。结果按 (文件, mtime_ns, k, lethal_ids) 进程内缓存。
    """
    if not lethal_ids:
        return np.empty(0, dtype=np.uint64)
    fa = Path(fa_path)
    try:
        mtime = fa.stat().st_mtime_ns
    except OSError:
        mtime = None
    key = (str(fa), mtime, k, tuple(sorted(lethal_ids)))
    hit = _LETHAL_INDEX_CACHE.get(key)
    if hit is not None:
        return hit
    arrs = []
    for name, seq in _fasta_iter(fa):
        if cancel is not None and cancel.is_set():
            raise RuntimeError('用户取消')
        if name in lethal_ids or _strip_version(name) in lethal_ids:
            for strand in (seq, revcomp(seq)):
                a = encode_kmers(strand, k)
                if a.size:
                    arrs.append(a)
    arr = (np.unique(np.concatenate(arrs)) if arrs
           else np.empty(0, dtype=np.uint64))
    if mtime is not None:
        if len(_LETHAL_INDEX_CACHE) >= 16:      # 简单上限，防长期驻留进程膨胀
            _LETHAL_INDEX_CACHE.pop(next(iter(_LETHAL_INDEX_CACHE)))
        _LETHAL_INDEX_CACHE[key] = arr
    return arr


def scan_arms(arm_seqs, species_indexes, k: int = K_DEFAULT, max_mm: int = 1,
              log=None, cancel=None) -> dict:
    """对每条 dsRNA 臂逐 siRNA（k-mer）扫描各物种索引的脱靶情况。

    参数
    ----
    arm_seqs : [(arm_name, seq_str)]          候选臂（按基因组顺序）
    species_indexes : {物种名: {'arr': ndarray,          build_species_index 产物
                                'lethal_arr': ndarray|None,  致死子索引（可无）
                                'bp': int, 'n_records': int}}
    cancel : 可选，threading.Event 风格（.is_set()）；置位时抛 RuntimeError
             （大面板 ≤2mm 扫描耗时较长，需可中断）

    语义（与参考实现一致）：
      - 查询只走正向——索引双链入集（红线 3），负链命中无需二次查询；
      - 每条 siRNA 先 mm0 再 mm1（max_mm>=2 再 mm2），按最小错配档计一次，
        即 best_tier 语义（红线 1）；
      - lethal（致死基因）只判定到 mm1 档；
      - 逐 siRNA 批量：变体集合 np.unique 后一次 searchsorted（参考实现
        has_hit 思路），按 siRNA 分批以保留逐 siRNA 归属。

    返回（全部可 JSON 序列化，无 numpy 数组/标量）：
    {
      'k': k, 'max_mm': max_mm,
      'species': {name: {'bp', 'n_records', 'distinct_kmers'}},
      'arms': {arm_name: {
          'n_sirna': int,
          'per_species': {sp: {'mm0','mm1','mm2'}},   # 命中 siRNA 条数，按最小错配档计一次
          'lethal':    {sp: {'mm0','mm1'}},           # 只对有 lethal_arr 的物种出键
          'lethal_rows': [{'sirna_idx','species','mm'}],  # 逐条致死命中，上限 500 截断
          'truncated': bool,                          # lethal_rows 是否被截断
          'red_line': bool,                           # 任一 lethal mm0 命中
          'sirna_tiers': {sp: {'tiers': [0/1/2/None, ...],   # 长度=n_sirna，best_tier
                               'lethal_tiers': [...] | None}},  # 无 lethal_arr 时为 None
      }},
    }
    """
    if not species_indexes:
        raise PanelError("species_indexes 为空：至少需要一个物种索引")
    for sp, idx in species_indexes.items():
        if not isinstance(idx, dict) or "arr" not in idx:
            raise PanelError(f"物种 {sp} 的索引缺少 'arr' 键（须为 "
                             f"build_species_index 的产物）")

    species_info = {}
    for sp, idx in species_indexes.items():
        species_info[sp] = {
            "bp": int(idx.get("bp", 0)),
            "n_records": int(idx.get("n_records", 0)),
            "distinct_kmers": int(np.asarray(idx["arr"]).size),
        }

    arms = {}
    for arm_name, seq in arm_seqs:
        # 臂自身 k-mer 按基因组顺序即 siRNA 顺序；查询单链正向（红线 3）
        fwd = encode_kmers(seq, k)
        n_si = int(fwd.size)
        per_species = {sp: {"mm0": 0, "mm1": 0, "mm2": 0}
                       for sp in species_indexes}
        lethal = {}
        lethal_rows = []
        truncated = False
        red_line = False
        sirna_tiers = {}

        for sp, idx in species_indexes.items():
            host = np.asarray(idx["arr"])
            larr = idx.get("lethal_arr")
            larr = None if larr is None else np.asarray(larr)
            tiers = [None] * n_si
            lethal_tiers = None if larr is None else [None] * n_si

            for i in range(n_si):
                if cancel is not None and cancel.is_set():
                    raise RuntimeError('用户取消')
                qi = int(fwd[i])
                # 红线 1：_best_tier 内部固定 0 -> 1 -> 2 顺序
                t = _best_tier(qi, host, k, max_mm)
                tiers[i] = t
                if t is not None:
                    per_species[sp][f"mm{t}"] += 1
                if larr is not None and larr.size:
                    # 致死子索引只判到 mm1 档
                    lt = _best_tier(qi, larr, k, 1)
                    lethal_tiers[i] = lt
                    if lt is not None:
                        lethal[sp] = lethal.get(sp, {"mm0": 0, "mm1": 0})
                        lethal[sp][f"mm{lt}"] += 1
                        if lt == 0:
                            red_line = True
                        if len(lethal_rows) < _LETHAL_ROW_CAP:
                            # sirna_idx 用 1 基（与 TSV 契约/server_species_scan
                            # 及 UI 展示口径一致；i 是 0 基内部下标）
                            lethal_rows.append({"sirna_idx": i + 1,
                                                "species": sp, "mm": lt})
                        else:
                            truncated = True
            sirna_tiers[sp] = {"tiers": tiers, "lethal_tiers": lethal_tiers}

        arms[str(arm_name)] = {
            "n_sirna": n_si,
            "per_species": per_species,
            "lethal": lethal,
            "lethal_rows": lethal_rows,
            "truncated": truncated,
            "red_line": bool(red_line),
            "sirna_tiers": sirna_tiers,
        }
        _log(log, f"[dsrna_offtarget] 臂 {arm_name}: {n_si} 个 siRNA, "
                  f"per_species={ {sp: c for sp, c in per_species.items()} }, "
                  f"red_line={red_line}")

    return {"k": int(k), "max_mm": int(max_mm),
            "species": species_info, "arms": arms}


# ---------------- 自测 ----------------

def selftest() -> bool:
    """移植参考实现 selftest 的全部断言，另补两条本模块特有检查。

    全部通过返回 True（供集成方以 exit code / bool 消费）。
    """
    k = K_DEFAULT
    # --- 参考实现断言 1/2：臂自身作为宿主库（双链入集），正反向查询全 tier 0 ---
    arm = "ACTTTCAGTTAGATGATTCTCTAGATGAACATTCACACTGATATAGTAGGCAAAGGTCTTAA"
    host = np.unique(np.concatenate(
        [encode_kmers(arm, k), encode_kmers(revcomp(arm), k)]))
    ok = all(_best_tier(int(x), host, k, 2) == 0 for x in encode_kmers(arm, k))
    ok_rc = all(_best_tier(int(x), host, k, 2) == 0
                for x in encode_kmers(revcomp(arm), k))

    # --- 参考实现断言 3：负链-only 命中必须能被正向查询找到 ---
    # 宿主正链只含 revcomp(q)（q 本身不在正链里），双链入集后正向查 q 应 tier 0
    q = arm[:k]
    host_plus = "AAAA" + revcomp(q) + "AAAA"
    both = np.unique(np.concatenate(
        [encode_kmers(host_plus, k), encode_kmers(revcomp(host_plus), k)]))
    ok_minus = _best_tier(int(encode_kmers(q, k)[0]), both, k, 2) == 0

    # --- 参考实现断言 4（反例）：只 complement 不 reverse 建库 -> 查不到 ---
    complement_only = np.unique(encode_kmers(host_plus.translate(_RC_STR), k))
    old_bug = _best_tier(int(encode_kmers(q, k)[0]), complement_only, k, 2)

    # --- 参考实现断言 5/6：1mm 邻居不得记成精确；2mm 邻居落 mm2 档 ---
    # （分层顺序是本 selftest 最初抓到的 bug：先查 1mm 会吞掉精确命中）
    k0 = int(encode_kmers(arm, k)[0])
    host_fwd = np.unique(encode_kmers(arm, k))
    t_one = _best_tier(k0 ^ 1, host_fwd, k, 2)
    t_two = _best_tier(k0 ^ 3 ^ (3 << (2 * (k - 1))), host_fwd, k, 2)

    # --- 补充断言 1：编码顺序（红线 5）、N 断窗、k 参数化 ---
    # "ACGT" 4-mer 手算值 0b00_01_10_11 = 27，锁死 A=0 C=1 G=2 T=3 顺序
    ok_code = int(encode_kmers("ACGT", k=4)[0]) == 0b00011011
    s44 = "ACGT" * 11
    s_n = s44[:22] + "N" + s44[23:]
    # 含 N 序列的 k-mer 集合 == 两侧干净片段 k-mer 集合的并（跨 N 窗口全丢）
    ok_n = np.array_equal(np.unique(encode_kmers(s_n, k)),
                          np.union1d(encode_kmers(s44[:22], k),
                                     encode_kmers(s44[23:], k)))
    ok_kpar = encode_kmers(s44, k=15).size == len(s44) - 15 + 1

    # --- 补充断言 2：scan_arms 端到端（2000bp 随机基因组嵌入 60nt 外源片段）---
    rng = random.Random(20260916)
    genome = "".join(rng.choice("ACGT") for _ in range(2000))
    frag = "".join(rng.choice("ACGT") for _ in range(60))
    emb = genome[:500] + frag + genome[560:]
    idx_arr = np.unique(np.concatenate(
        [encode_kmers(emb, k), encode_kmers(revcomp(emb), k)]))
    larr = np.unique(np.concatenate(
        [encode_kmers(frag, k), encode_kmers(revcomp(frag), k)]))

    def _mut(s, pos):
        return s[:pos] + {"A": "C", "C": "A", "G": "T", "T": "G"}[s[pos]] + s[pos + 1:]

    arm_1mm = _mut(frag, 10)      # 位点 10 落在窗口 0..10 -> 11 个 mm1、29 个 mm0
    arm_2mm = _mut(_mut(frag, 5), 10)  # 窗口 0..5 双突变 -> 6 个 mm2；6..10 mm1
    res = scan_arms(
        [("exact", frag), ("one_mm", arm_1mm), ("two_mm", arm_2mm)],
        {"syn": {"arr": idx_arr, "lethal_arr": larr,
                 "bp": len(emb), "n_records": 1}},
        k=k, max_mm=2)
    res1 = scan_arms(
        [("two_mm", arm_2mm)],
        {"syn": {"arr": idx_arr, "lethal_arr": None,
                 "bp": len(emb), "n_records": 1}},
        k=k, max_mm=1)
    ps = {n: res["arms"][n]["per_species"]["syn"]
          for n in ("exact", "one_mm", "two_mm")}
    # max_mm=1 时 two_mm 的 6 个双突变窗口必须整体无解（None），不许被误记进
    # mm1 档；真正只差 1 个错配的窗口 6..10 仍应落在 mm1
    ok_e2e = (
        res["arms"]["exact"]["n_sirna"] == 40
        and ps["exact"] == {"mm0": 40, "mm1": 0, "mm2": 0}
        and ps["one_mm"] == {"mm0": 29, "mm1": 11, "mm2": 0}
        and ps["two_mm"] == {"mm0": 29, "mm1": 5, "mm2": 6}
        and res1["arms"]["two_mm"]["per_species"]["syn"]
        == {"mm0": 29, "mm1": 5, "mm2": 0}
        and res1["arms"]["two_mm"]["sirna_tiers"]["syn"]["tiers"][:6]
        == [None] * 6
        and res["arms"]["exact"]["red_line"] is True
        and res["arms"]["exact"]["lethal"]["syn"] == {"mm0": 40, "mm1": 0}
        and res["arms"]["two_mm"]["lethal"]["syn"] == {"mm0": 29, "mm1": 5}
        # 嵌入片段的 40 个 siRNA 逐条 tier 全 0（宿主与 lethal 皆是）
        and res["arms"]["exact"]["sirna_tiers"]["syn"]["tiers"] == [0] * 40
        and res["arms"]["exact"]["sirna_tiers"]["syn"]["lethal_tiers"] == [0] * 40
        and isinstance(json.dumps(res), str)   # 返回值可 JSON 序列化
    )

    print(f"exact-vs-host (fwd queries)   -> {ok}")
    print(f"exact-vs-host (rc queries)    -> {ok_rc}")
    print(f"minus-strand-only hit found   -> {ok_minus}")
    print(f"  complement-only build misses it: tier={old_bug} (expect != 0)")
    print(f"1 mm neighbour -> tier {t_one} (expect 1)")
    print(f"2 mm neighbour -> tier {t_two} (expect 2)")
    print(f"code order / N-split / k-param-> "
          f"{ok_code and ok_n and ok_kpar}")
    print(f"scan_arms e2e (mm0/1mm/2mm tiers, lethal, JSON) -> {ok_e2e}")
    good = (ok and ok_rc and ok_minus and old_bug != 0
            and t_one == 1 and t_two == 2
            and ok_code and ok_n and ok_kpar and ok_e2e)
    print("SELFTEST_OK" if good else "SELFTEST_FAIL")
    return good


if __name__ == "__main__":
    try:                                    # Windows 控制台中文输出兜底
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, OSError):
        pass
    sys.exit(0 if selftest() else 1)
