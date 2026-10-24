#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
prescreen_engine.py — SRA / GSA 病毒样本预筛引擎

移植自 MMPV-RNA `public_metadata_pipeline/stat_prescreen.py`（2026-10-01）。
用途不变：给定一批 SRA / GSA Run accession，取回上游**已算好的**物种分类
组成，**不下载任何测序数据**即可判断哪些样本含病毒，用于把上万条候选 run
先筛一遍，只把真正含病毒的送进下载与下游分析。

    1) stat_run_summary.tsv   —— 每 run 一行：测序量、病毒占比、主要物种
    2) stat_viral_taxa.tsv    —— 每 run 的病毒类 taxon 明细（长表）
    3) prescreen_pass.tsv     —— （可选）通过阈值的 run 清单
    4) unresolved.tsv / .list —— 上游没有分类数据的 run（需本地比对）

数据来源
--------
    NCBI (SRR/ERR/DRR)
      GET https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy
            ?acc={RUN}&cluster_name=public
      ⚠ 必须带 cluster_name=public，否则 HTTP 200 + 空 body（隐蔽坑）。

    CNCB GSA (CRR)
      GET https://ngdc.cncb.ac.cn/gsa/taxonomyAnalysis/treePath/{CRA}/{CRR}/{CRR}.json
      URL 需要 (CRA, CRR) 两个号；先解析 CRA 再验证。

口径说明（来自 NCBI 官方答复 ngs-tools#43）
-------------------------------------------
    analysed     : 分母。被赋予 taxid 的 spot 总数。
    total_count  : 该节点**含全部后代**的累计计数。
    self_count   : 该节点**自身**命中数（不含后代）。
    identified   : 至少命中 1 个 k-mer 的 spot 数。
    占比 pct     = total_count / analysed * 100

    ⚠ GSA 的 percent 分母是**总 reads**（identified + unidentified），
      STAT 的分母是 analysed（已鉴定 reads）——两者不可直接比较，
      跨来源比较一律用归一化列 `virus_pct_ident`（病毒命中 / 已鉴定命中）。

覆盖范围（上游实测 2026-10，470 条真实样本）
--------------------------------------------
    SRR (NCBI) / DRR (DDBJ) / ERR (ENA)  ✅ 基本都有
    CRR (CNCB/GSA)                       ⚠ 只有一部分有（实测 148 条中 60 条）
    缺失主要但不完全按项目聚集 —— 必须逐 run 判，不能按项目判。

本平台版本与上游的差异（移植时按平台约定做的改动，逐条可查）
------------------------------------------------------------
1. **病毒判定复用平台 taxonomy**：上游自带一个 nodes.dmp 解析类，需要用户
   准备 `--taxdump`。平台 `databases/tax_db` 三件套随分发，改用
   `taxonomy.TaxonomyGraph`（且它额外处理 merged.dmp）。`--taxdump` 从
   必需变为可选覆盖项，新增 `--no-taxdump` 退回名称启发式。
   同时**修掉上游一个死参数**：`--taxids` 上游只用于打印提示、判定时从未
   生效；本版让它真正加入目标子树集合。
   ⚠ 注意：子树判定与名称启发式是**互补**关系，最终判定取并集 ——
   实测单用任一边都会漏（详见 `is_viral` 的注释，有 470 run 的实测数字）。
2. **出站请求过安全校验**：上游只有超时/重试，没有任何地址边界校验。
   本版每次请求前走 `public_data.check_url_host`（主机白名单 + 解析 IP
   阻断私网/环回），并复用平台的 `_OPENER` 使**重定向逐跳复校**。
3. **去掉 `pipeline_utils.UI` 依赖**（平台无此模块），内置等价 UI。
4. **默认输出目录**改到平台 `run/meta_search/` 下，避免落到 CWD。

保留的上游关键纪律（改动时不得回退）
------------------------------------
  * `cluster_name=public` 不能省。
  * GSA 的 `CRA000112` 是搜索页模板噪声，必须排除。
  * GSA 只有**确认态**（NO_ANALYSIS / NOT_FOUND）才写缓存；429/网络类
    `FAILED` 绝不缓存 —— 上游注释指出，早期把临时故障写成"无数据"是
    一个真 bug，会永久污染结论。缓存带 `_v` 版本号，旧格式一律重抓。
  * 「最特异的检出物种」用 self_count 选（用 total_count 会永远选中
    Viruses 这类高层节点）。
  * `top_organism` 必须排除 root(1) 与 cellular organisms(131567)。
  * 有病毒信号但没细分到具名节点时，填 "Viruses (未细分)"，不留空。

用法
----
    # 从文本列表
    python prescreen_engine.py --runs sra.list -o <outdir>
    # 从平台「公共数据检索」产出的元数据 CSV 取 Run 列
    python prescreen_engine.py --csv Global_Unified_Metadata_Core14.csv -o <outdir>
    # 只看病毒占比 > 0.5%，并单独导出通过清单
    python prescreen_engine.py --runs sra.list -o <outdir> --virus-pct 0.5
    # 断点续跑（默认开启，已下载的 JSON 直接复用）
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# ---- 平台模块导入：独立脚本（Web 子进程）与包内 import 两种方式都要成立 ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))          # public_meta → 平台根
try:                                                     # 独立脚本 / --run-engine
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from Virus_Platform_Core.config import DIRS                      # type: ignore
    from Virus_Platform_Core.public_data import (                    # type: ignore
        _OPENER, check_url_host)
    from Virus_Platform_Core.taxonomy import (                       # type: ignore
        TaxonomyGraph, taxonomy_ready)
except ImportError:                                      # 作为包被 import 时
    from ..config import DIRS                            # type: ignore
    from ..public_data import _OPENER, check_url_host     # type: ignore
    from ..taxonomy import TaxonomyGraph, taxonomy_ready  # type: ignore


# ============================================================================
# UI（内置；平台无 pipeline_utils，且 Web 侧以子进程捕获 stdout 进任务日志）
# ============================================================================
class UI:
    C = {"cyan": "\033[96m", "green": "\033[92m", "yellow": "\033[93m",
         "red": "\033[91m", "gray": "\033[90m", "purple": "\033[95m",
         "bold": "\033[1m", "reset": "\033[0m"}

    @staticmethod
    def banner(t): print(f"\n{UI.C['purple']}{UI.C['bold']}  {t}{UI.C['reset']}\n")

    @classmethod
    def stage(cls, n, s="start"): print(f"{cls.C['cyan']}[>>>] {n}{cls.C['reset']}")

    @classmethod
    def ok(cls, m): print(f"  {cls.C['green']}[OK]{cls.C['reset']} {m}")

    @classmethod
    def warn(cls, m): print(f"  {cls.C['yellow']}[WARN]{cls.C['reset']} {m}")

    @classmethod
    def err(cls, m): print(f"  {cls.C['red']}[ERR]{cls.C['reset']} {m}")

    @classmethod
    def info(cls, m): print(f"  {cls.C['gray']}> {m}{cls.C['reset']}")


# ============================================================================
# 常量
# ============================================================================
STAT_URL = ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy"
            "?acc={acc}&cluster_name=public")

TAXID_VIRUSES = 10239          # NCBI: Viruses
TAXID_ROOT = 1
TAXID_CELLULAR = 131567        # cellular organisms

# 名字兜底启发式：没有 taxonomy 时用它从 tax_table 里挑病毒行。
# 「够用但不严谨」——报告里必须标注为名称启发式。
VIRAL_NAME_RE = re.compile(
    r"(viridae|virales|virinae|virus|viruses|viroid|virophage|"
    r"\bphage\b|bacteriophage|retroelement|Caulimoviridae)",
    re.IGNORECASE,
)

USER_AGENT = "virus-platform-prescreen/1.0 (research; local pipeline)"

# ---- GSA (CNCB / NGDC) ----
GSA_BASE = "https://ngdc.cncb.ac.cn"
GSA_SEARCH = GSA_BASE + "/gsa/search?searchTerm={acc}"
GSA_TREE = GSA_BASE + "/gsa/taxonomyAnalysis/treePath/{cra}/{acc}/{acc}.json"
# 搜索页模板里恒定出现的噪声 CRA，必须排除，否则会把它当成候选
GSA_CRA_NOISE = {"CRA000112"}
# GSA 缓存格式版本。改动缓存结构时 +1，旧缓存会被自动作废并重抓。
GSA_CACHE_V = 2

# GSA 树里的病毒节点名（GSA 有顶层 "Viruses" 节点，名字判定比 STAT 更可靠）
GSA_VIRAL_RE = re.compile(
    r"(^Viruses$|viridae|virales|virinae|virus|viruses|viroid|virophage|phage)",
    re.IGNORECASE,
)


# ============================================================================
# 出站请求安全（复用平台校验：主机白名单 + 解析 IP + 重定向逐跳复校）
# ============================================================================
_VERIFIED_HOSTS: Set[str] = set()
_VERIFY_LOCK = threading.Lock()


def _guard(url: str) -> None:
    """请求前安全校验。同一 host 校验一次即可 —— 重定向的逐跳校验由
    平台 `_OPENER` 的 `_SafeRedirectHandler` 负责，不受此缓存影响。"""
    host = (urllib.parse.urlsplit(url).hostname or "").lower()
    with _VERIFY_LOCK:
        if host in _VERIFIED_HOSTS:
            return
    ok, why = check_url_host(url)
    if not ok:
        raise ValueError(f"出站请求被拒绝: {why}")
    with _VERIFY_LOCK:
        _VERIFIED_HOSTS.add(host)


def _open_url(url: str, timeout: int):
    """校验 + 打开（走平台 opener，重定向逐跳复校）。"""
    _guard(url)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return _OPENER.open(req, timeout=timeout)


# ============================================================================
# 病毒子树判定（复用平台 TaxonomyGraph）
# ============================================================================
def make_viral_judge(tax_dir: Optional[str] = None,
                     extra_taxids: Tuple[int, ...] = ()
                     ) -> Optional[Callable[[int], bool]]:
    """构造「taxid 是否落在目标子树内」的判定函数。

    默认目标是 NCBI Viruses(10239)；`extra_taxids` 可追加自定义目标
    （对应 CLI 的 `--taxids`，上游该参数是死参数，本版让它真正生效）。

    taxonomy 不可用时返回 None，调用方退回名称启发式。
    判定结果按 taxid 缓存（上游等价物是 Taxonomy._cache）。
    """
    if not taxonomy_ready(tax_dir):
        return None
    tg = TaxonomyGraph(tax_dir)
    targets = {TAXID_VIRUSES} | {int(t) for t in extra_taxids}
    cache: Dict[int, bool] = {}

    def judge(taxid: Any) -> bool:
        t = int(taxid)
        if t in cache:
            return cache[t]
        hit = any(x in targets for x, _r, _n in tg.lineage(t))
        cache[t] = hit
        return hit

    return judge


# ============================================================================
# 抓取：NCBI STAT
# ============================================================================
class Fetcher:
    """带重试/退避/限速的 STAT 抓取器，落盘缓存以支持断点续跑。"""

    def __init__(self, raw_dir: str, retries: int = 4, timeout: int = 90,
                 min_interval: float = 0.0) -> None:
        self.raw_dir = raw_dir
        self.retries = retries
        self.timeout = timeout
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._last = 0.0
        os.makedirs(raw_dir, exist_ok=True)

    def _throttle(self) -> None:
        if self.min_interval <= 0:
            return
        with self._lock:
            gap = time.time() - self._last
            if gap < self.min_interval:
                time.sleep(self.min_interval - gap)
            self._last = time.time()

    def cache_path(self, acc: str) -> str:
        return os.path.join(self.raw_dir, f"{acc}.json")

    def fetch(self, acc: str, use_cache: bool = True
              ) -> Tuple[str, Optional[Dict[str, Any]], str, float]:
        """返回 (acc, payload|None, status, seconds)；status ∈ {OK, NO_STAT, FAILED}。"""
        cp = self.cache_path(acc)
        if use_cache and os.path.exists(cp):
            try:
                with open(cp, "r", encoding="utf-8") as f:
                    txt = f.read()
                if txt.strip() == "":          # 缓存了「无数据」这一事实
                    return acc, None, "NO_STAT", 0.0
                return acc, json.loads(txt)[0], "OK", 0.0
            except Exception:
                pass                            # 缓存损坏 → 重新抓

        url = STAT_URL.format(acc=urllib.parse.quote(acc))
        last_err = ""
        t0 = time.time()
        for attempt in range(self.retries):
            self._throttle()
            try:
                with _open_url(url, self.timeout) as resp:
                    raw = resp.read()
                if not raw:
                    # 空 body = 该 accession 无 STAT 结果（CRR / 新 run / 不存在）
                    with open(cp, "w", encoding="utf-8") as f:
                        f.write("")
                    return acc, None, "NO_STAT", time.time() - t0
                payload = json.loads(raw)[0]
                with open(cp, "w", encoding="utf-8") as f:
                    json.dump(payload, f, ensure_ascii=False)
                return acc, payload, "OK", time.time() - t0
            except urllib.error.HTTPError as e:
                last_err = f"HTTP {e.code}"
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(min(2 ** attempt + random.random(), 30))
                    continue
                break
            except ValueError as e:              # 安全校验拒绝 —— 不该重试
                return acc, None, f"FAILED (blocked: {e})", time.time() - t0
            except Exception as e:               # 超时 / 连接重置
                last_err = f"{type(e).__name__}: {e}"
                time.sleep(min(2 ** attempt + random.random(), 30))
        return acc, None, f"FAILED ({last_err})", time.time() - t0


# ============================================================================
# 抓取：CNCB GSA（CRR）
# ============================================================================
class GsaFetcher:
    """抓取 NGDC/GSA 的 CRR 分类树。

    与 NCBI 不同，GSA 的 URL 需要 (CRA, CRR) 两个 accession，而输入通常
    只有 CRR，因此要先解析 CRA：
      1. `--cra-map` 提供的映射（最快，推荐大批量）
      2. 抓 gsa/search?searchTerm={CRR}，正则出候选 CRA（排除模板噪声），
         再逐个用 treePath 验证 —— 验证过才算数，不靠猜。
    """

    def __init__(self, raw_dir: str, retries: int = 4, timeout: int = 60,
                 min_interval: float = 0.4,
                 cra_map: Optional[Dict[str, str]] = None) -> None:
        self.raw_dir = raw_dir
        self.retries = retries
        self.timeout = timeout
        self.min_interval = min_interval
        self.cra_map = cra_map or {}
        self._lock = threading.Lock()
        self._last = 0.0
        self._cra_cache: Dict[str, List[str]] = {}
        os.makedirs(raw_dir, exist_ok=True)

    def _throttle(self) -> None:
        with self._lock:
            gap = time.time() - self._last
            if gap < self.min_interval:
                time.sleep(self.min_interval - gap)
            self._last = time.time()

    def _get(self, url: str) -> str:
        last_err: Optional[Exception] = None
        for attempt in range(self.retries):
            self._throttle()
            try:
                with _open_url(url, self.timeout) as r:
                    return r.read().decode("utf-8", errors="replace")
            except urllib.error.HTTPError as e:
                last_err = e
                if e.code in (429, 500, 502, 503, 504):
                    # NGDC 对高频访问会 429，退避要够长
                    time.sleep(min(3 * (attempt + 1) + random.random(), 45))
                    continue
                raise
            except ValueError:
                raise                                # 安全校验拒绝，不重试
            except Exception as e:
                last_err = e
                time.sleep(min(2 ** attempt + random.random(), 30))
        raise last_err if last_err else RuntimeError("unknown")

    def cra_candidates(self, acc: str) -> List[str]:
        if acc in self._cra_cache:
            return self._cra_cache[acc]
        if acc in self.cra_map:
            cands = [self.cra_map[acc]]
        else:
            html = self._get(GSA_SEARCH.format(acc=acc))
            cands = [c for c in sorted(set(re.findall(r"CRA\d{6}", html)))
                     if c not in GSA_CRA_NOISE]
        self._cra_cache[acc] = cands
        return cands

    def cache_path(self, acc: str) -> str:
        return os.path.join(self.raw_dir, f"{acc}.json")

    def _write_cache(self, cp: str, status: str, cra: str = "",
                     detail: str = "", tree: Optional[Dict[str, Any]] = None) -> None:
        """统一缓存包装。

        带 `_v` 版本号是有意的：旧格式（空文件 / 裸 tree）一律视为不可信
        并重抓，否则旧版把「临时故障」误记成「无数据」的结果会被一直沿用。
        """
        payload: Dict[str, Any] = {"_v": GSA_CACHE_V, "_status": status,
                                   "_cra": cra, "_detail": detail}
        if tree is not None:
            payload["_tree"] = tree
        with open(cp, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False)

    def fetch(self, acc: str, use_cache: bool = True, recheck_nodata: bool = False
              ) -> Tuple[str, Optional[Dict[str, Any]], str, float, str, str]:
        """返回 (acc, tree|None, status, seconds, cra, detail)。

        status：
          OK           拿到分类树
          NO_ANALYSIS  CRA 解析成功但该 run 无分类分析（treePath 404）—— 确认态
          NOT_FOUND    搜索页找不到 CRA 候选，accession 不在 GSA —— 确认态
          FAILED       429 / 网络 / 超时等**临时**故障 —— 绝不缓存

        只有「确认态」才写缓存。临时故障一旦被写成 NO_DATA，就会永久污染
        结果（上游注释指出这正是早期版本的一个 bug）。
        """
        cp = self.cache_path(acc)

        if use_cache and os.path.exists(cp):
            try:
                obj = json.loads(open(cp, "r", encoding="utf-8").read() or "null")
            except Exception:
                obj = None
            # 只认当前版本的缓存包装；旧格式（空文件 / 裸 tree）一律重抓
            if isinstance(obj, dict) and obj.get("_v") == GSA_CACHE_V:
                st = obj.get("_status", "")
                if st == "OK":
                    tree = obj.get("_tree")
                    if isinstance(tree, dict):
                        return acc, tree, "OK", 0.0, obj.get("_cra", ""), ""
                elif not recheck_nodata:
                    # --recheck-nodata：忽略「无数据」标记，重新问一次
                    # （GSA 会陆续补做分析，昨天的无数据今天可能就有了）
                    return (acc, None, st, 0.0,
                            obj.get("_cra", ""), obj.get("_detail", ""))

        t0 = time.time()

        # ---- 1) 解析 CRA ----
        try:
            cands = self.cra_candidates(acc)
        except Exception as e:
            # 解析失败多半是限流/网络 —— 临时故障，不缓存
            return acc, None, f"FAILED (CRA解析 {type(e).__name__})", time.time() - t0, "", ""

        if not cands:
            d = "搜索页无 CRA 候选（accession 不在 GSA）"
            self._write_cache(cp, "NOT_FOUND", detail=d)
            return acc, None, "NOT_FOUND", time.time() - t0, "", d

        # ---- 2) 逐个候选 CRA 取树 ----
        confirmed_missing = 0        # 明确 404 的候选数
        transient: List[str] = []    # 临时故障的候选
        for cra in cands:
            try:
                raw = self._get(GSA_TREE.format(cra=cra, acc=acc))
                if not raw.strip():
                    confirmed_missing += 1
                    continue
                tree = json.loads(raw)
                if not isinstance(tree, dict):
                    confirmed_missing += 1
                    continue
                self._write_cache(cp, "OK", cra=cra, tree=tree)
                return acc, tree, "OK", time.time() - t0, cra, ""
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    confirmed_missing += 1
                    continue
                transient.append(f"{cra}:HTTP{e.code}")
            except json.JSONDecodeError:
                confirmed_missing += 1
            except Exception as e:
                transient.append(f"{cra}:{type(e).__name__}")

        # 只要有一个候选是**临时**故障，就不能断定「没有数据」
        if transient:
            return (acc, None, f"FAILED ({';'.join(transient)})", time.time() - t0, "",
                    "候选 CRA 全部请求失败（临时故障，未缓存）")

        cra_str = "/".join(cands)
        detail = (f"候选 CRA {cra_str} 均返回 404，"
                  f"GSA 未对该 run 做分类分析（可日后用 --recheck-nodata 复查）")
        self._write_cache(cp, "NO_ANALYSIS", cra=cra_str, detail=detail)
        return acc, None, "NO_ANALYSIS", time.time() - t0, cra_str, detail


# ============================================================================
# 解析：NCBI STAT
# ============================================================================
def parse_payload(payload: Dict[str, Any],
                  judge: Optional[Callable[[Any], bool]]
                  ) -> Dict[str, Any]:
    """把 STAT payload 提炼成一行汇总 + 病毒明细。"""
    tot = payload.get("tax_totals", {}) or {}
    table = payload.get("tax_table", []) or []

    analysed = int(tot.get("analysed") or tot.get("total") or 0)
    identified = int(tot.get("identified") or 0)
    denom = analysed if analysed > 0 else 0

    def pct(c: int) -> float:
        return round(100.0 * c / denom, 6) if denom else 0.0

    def is_viral(tid: Any, org: str) -> bool:
        """病毒判定 = taxonomy 子树 **∪** 名称启发式（互补，不是二选一）。

        为什么必须取并集 —— 实测依据（2026-10-01，470 条真实 run，
        上游 all470 基线的 viral_taxa 明细逐条对比）：

          * **taxonomy 独有** 4777 条明细 / 91 个 taxon 名：全是
            Kitrinoviricota、Orthornavirae、Caudoviricetes、Alsuviricetes、
            Betasatellite 这类**新版 ICTV 高阶分类单元**（realm/kingdom/
            phylum/class）与卫星，名字里**没有** "virus" 字样 ——
            名称启发式 100% 漏判（实测「仅 B 有」集合里名字匹配数为 0）。
          * **名称启发式独有** 89 条明细 / 4 个名：Viroids、
            Cytorhabdovirus、Viruses incertae sedis、unclassified
            Kanagawavirus —— 这些 taxid 在平台 tax_db 的 nodes.dmp /
            names.dmp / merged.dmp 里**都查不到**（STAT 的 tax_table 出自
            另一版 taxonomy），纯子树判定会漏（实测 lineage 为空）。

        单用任何一边都会漏，并集两边都覆盖。注意 `virus_pct` / `virus_count`
        不受本判定影响（直接读 Viruses(10239) 节点的 total_count），
        本判定只决定 `stat_viral_taxa.tsv` 明细与 `top_virus` 的选取。
        """
        if judge is not None and judge(tid):
            return True
        return bool(VIRAL_NAME_RE.search(org))

    row: Dict[str, Any] = {
        "source": "NCBI_STAT",
        "analysed": analysed,
        "identified": identified,
        "identified_pct": pct(identified),
        "n_taxa": len(table),
        "virus_count": 0,
        "virus_pct": 0.0,
        "virus_pct_ident": 0.0,
        "top_organism": "",
        "top_organism_pct": 0.0,
    }

    viral_rows: List[Dict[str, Any]] = []
    virus_node = next((r for r in table if r.get("tax_id") == TAXID_VIRUSES), None)
    if virus_node:
        row["virus_count"] = int(virus_node.get("total_count", 0))
        row["virus_pct"] = pct(row["virus_count"])

    for r in table:
        tid = r.get("tax_id")
        if tid is None or tid in (TAXID_ROOT, TAXID_CELLULAR, TAXID_VIRUSES):
            continue
        if is_viral(tid, r.get("org", "")):
            viral_rows.append({
                "tax_id": int(tid),
                "org": r.get("org", ""),
                "path": "",
                "self_count": int(r.get("self_count", 0)),
                "total_count": int(r.get("total_count", 0)),
                "pct": pct(int(r.get("total_count", 0))),
            })

    # 病毒明细按 total_count 降序（父节点在前，便于看层级）
    viral_rows.sort(key=lambda x: -x["total_count"])
    row["_viral_rows"] = viral_rows

    # 若没有显式的 Viruses(10239) 行，用检出的病毒节点里 total_count 最大者兜底
    if not virus_node and viral_rows:
        row["virus_count"] = viral_rows[0]["total_count"]
        row["virus_pct"] = viral_rows[0]["pct"]

    # 归一化口径：病毒命中 / 已鉴定命中。用于和 GSA 跨来源比较
    # （GSA 的 percent 分母是总 reads，STAT 的是 analysed，直接比会低估 GSA）
    if identified > 0:
        row["virus_pct_ident"] = round(100.0 * row["virus_count"] / identified, 6)

    # 「最特异的病毒」：优先从**有具体名字**的行里选（名字匹配病毒正则），
    # 再按 self_count 取最大（直接命中该节点的读数，不含后代继承）。
    #
    # 为什么不能直接在全集上取 max(self_count)：
    #   taxonomy 子树含 Riboviria / Orthornavirae / Caudoviricetes 这类
    #   realm / kingdom / class 级节点，STAT 会把**无法细分到低阶**的 reads
    #   直接记在这些节点上 —— 于是它们的 self_count 不为 0，会把具体物种名
    #   挤掉，选出来毫无信息量（实测 26 处退化，如
    #   `Cauliflower mosaic virus` → `Riboviria`）。
    #   上游纯名称启发式恰好天然排除了这类节点（名字不含 "virus"），
    #   这里显式保留该效果：有具体名可选就优先具体名；一个具体名都没有
    #   （上游会填 "Viruses (未细分)"）才退回高阶节点 —— 那反而是增益，
    #   例如上游 `Viruses (未细分)` → 本版 `Caudoviricetes sp.`。
    if viral_rows:
        named = [x for x in viral_rows if VIRAL_NAME_RE.search(x["org"])]
        pool = named or viral_rows
        best_v = max(pool, key=lambda x: (x["self_count"], x["total_count"]))
        row["top_virus"] = best_v["org"]
        row["top_virus_pct"] = best_v["pct"]
    elif row["virus_count"] > 0:
        # 有病毒信号但没有细分到具名节点 —— 别留空，否则下游会以为没检出
        row["top_virus"] = "Viruses (未细分)"
        row["top_virus_pct"] = row["virus_pct"]
    else:
        row["top_virus"] = ""
        row["top_virus_pct"] = 0.0

    # 非病毒的主要物种：同样排除病毒子树，用 self_count 选最特异的那个
    non_viral = [
        r for r in table
        if r.get("tax_id") not in (TAXID_ROOT, TAXID_CELLULAR, TAXID_VIRUSES)
        and not is_viral(r.get("tax_id"), r.get("org", ""))
    ]
    if non_viral:
        best_o = max(non_viral, key=lambda x: (int(x.get("self_count", 0)),
                                               int(x.get("total_count", 0))))
        row["top_organism"] = best_o.get("org", "")
        row["top_organism_pct"] = pct(int(best_o.get("total_count", 0)))
    return row


# ============================================================================
# 解析：CNCB GSA
# ============================================================================
def _gsa_flat(node: Dict[str, Any], path: Tuple[str, ...] = ()) -> List[Dict[str, Any]]:
    """把 GSA 的嵌套树拍平成列表，保留完整路径（GSA 没有 tax_id，路径即层级证据）。"""
    p = path + (node.get("name", ""),)
    out = [{"name": node.get("name", ""), "hits": node.get("hits"),
            "percent": node.get("percent"), "path": " > ".join(p)}]
    for c in node.get("children") or []:
        out.extend(_gsa_flat(c, p))
    return out


def _gsa_pct(v: Any) -> float:
    """'0.04%' -> 0.04；无效返回 0。"""
    if not isinstance(v, str):
        return 0.0
    try:
        return float(v.rstrip("%").strip())
    except ValueError:
        return 0.0


def _gsa_hits(v: Any) -> int:
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return 0


def parse_gsa_tree(tree: Dict[str, Any]) -> Dict[str, Any]:
    """把 GSA 分类树提炼成与 STAT 路径同构的一行汇总。

    口径差异：GSA 的 percent 分母是**总 reads**，STAT 是 analysed。
    故额外算 virus_pct_ident = 病毒命中 / 已鉴定命中，用于跨来源比较。
    """
    rows = _gsa_flat(tree)
    root = tree
    ident_pct = _gsa_pct(root.get("identified"))
    unident_pct = _gsa_pct(root.get("unidentified"))

    # 用命中最多的节点反推总 reads： total = hits / (percent/100)
    total_reads = 0.0
    for r in rows:
        h, p = _gsa_hits(r["hits"]), _gsa_pct(r["percent"])
        if h > 0 and p > 0:
            total_reads = max(total_reads, h / (p / 100.0))
    ident_hits = int(round(total_reads * ident_pct / 100.0)) if total_reads else 0

    row: Dict[str, Any] = {
        "source": "GSA",
        # analysed 统一语义 = virus_pct 的分母。GSA 的分母是总 reads。
        "analysed": int(round(total_reads)) if total_reads else 0,
        "identified": ident_hits,
        "identified_pct": round(ident_pct, 4),
        # 已解析的未鉴定占比也要落进汇总行（此前只解析不入行，GSA 行缺这列）
        "unidentified_pct": round(unident_pct, 4),
        "n_taxa": len(rows) - 1,
        "virus_count": 0,
        "virus_pct": 0.0,
        "virus_pct_ident": 0.0,
        "top_organism": "",
        "top_organism_pct": 0.0,
    }

    viral: List[Dict[str, Any]] = []
    virus_root = next((r for r in rows
                       if (r["name"] or "").strip().lower() == "viruses"), None)
    for r in rows:
        nm = r["name"] or ""
        if nm.strip().lower() == "viruses":
            continue
        if not GSA_VIRAL_RE.search(nm):
            continue
        viral.append({
            "tax_id": "",
            "org": nm,
            "path": r["path"],
            "self_count": "",
            "total_count": _gsa_hits(r["hits"]),
            "pct": round(_gsa_pct(r["percent"]), 6),
        })

    if virus_root is not None:
        row["virus_count"] = _gsa_hits(virus_root["hits"])
        row["virus_pct"] = round(_gsa_pct(virus_root["percent"]), 6)

    # 没有顶层 Viruses 节点时，用病毒节点里的最大值兜底
    if not row["virus_count"] and viral:
        row["virus_count"] = max(v["total_count"] for v in viral)
        row["virus_pct"] = max(v["pct"] for v in viral)

    if ident_hits > 0:
        row["virus_pct_ident"] = round(100.0 * row["virus_count"] / ident_hits, 6)

    # 最特异的病毒 = 「没有病毒后代」的那些节点里命中最多的
    viral_names = {v["org"] for v in viral}
    leaves = [v for v in viral
              if not any(o["path"].startswith(v["path"] + " > ") and o["org"] in viral_names
                         for o in viral)]
    pool = leaves or viral
    if pool:
        best_v = max(pool, key=lambda x: x["total_count"])
        row["top_virus"] = best_v["org"]
        row["top_virus_pct"] = best_v["pct"]
    elif row["virus_count"] > 0:
        # 有病毒信号但没有细分到具名节点 —— 别留空，否则下游会以为没检出
        row["top_virus"] = "Viruses (未细分)"
        row["top_virus_pct"] = row["virus_pct"]

    viral.sort(key=lambda x: -x["total_count"])
    row["_viral_rows"] = viral

    # 主要非病毒物种：最深层且命中最多的非病毒、非根节点
    non_viral = [r for r in rows
                 if (r["name"] or "").strip().lower() not in ("", "viruses")
                 and not GSA_VIRAL_RE.search(r["name"] or "")
                 and r["name"] not in ("cellular organisms",)]
    if non_viral:
        maxdepth = max(r["path"].count(" > ") for r in non_viral)
        deep = [r for r in non_viral if r["path"].count(" > ") == maxdepth]
        best_o = max(deep, key=lambda x: _gsa_hits(x["hits"]))
        row["top_organism"] = best_o["name"]
        row["top_organism_pct"] = round(_gsa_pct(best_o["percent"]), 6)
    return row


# ============================================================================
# 输入读取
# ============================================================================
def fmt_cell(v: Any) -> Any:
    """TSV 里把浮点写成定点小数，避免 4e-06 这种科学计数法影响可读性。"""
    if isinstance(v, float):
        s = f"{v:.8f}".rstrip("0").rstrip(".")
        return s if s not in ("", "-0") else "0"
    return v


def read_accessions(args: argparse.Namespace) -> List[str]:
    accs: List[str] = []
    if args.runs:
        with open(args.runs, "r", encoding="utf-8-sig", errors="replace") as f:
            for line in f:
                s = line.strip().split()[0] if line.strip() else ""
                if s and not s.startswith("#"):
                    accs.append(s)
    if args.csv:
        with open(args.csv, "r", encoding="utf-8-sig", errors="replace") as f:
            for rec in csv.DictReader(f):
                v = (rec.get(args.run_col) or "").strip()
                if v:
                    accs.append(v)
    if args.acc:
        accs.extend(args.acc)

    seen: Set[str] = set()
    out: List[str] = []
    for a in accs:
        a = a.strip()
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return out


# ============================================================================
# 主流程
# ============================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="SRA / GSA 病毒样本预筛（基于 NCBI STAT 与 CNCB GSA 预计算分类结果）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
双通路自动路由
--------------
按 accession 前缀自动分流，无需手工指定：

  SRR / ERR / DRR  ->  NCBI STAT
      扁平表，带 tax_id；分母 analysed（已鉴定 reads）

  CRR              ->  CNCB GSA
      嵌套树，**无 tax_id**（只能按名字判病毒）；分母是**总 reads**

  => 两者分母不同，直接比较会低估 GSA。故输出多一列 virus_pct_ident
     （病毒命中 / 已鉴定命中），跨来源比较请用这一列。

口径与边界（务必转达使用者）
----------------------------
* NCBI / GSA 都是 **k-mer 快速分类**，不是严谨比对。低丰度病毒可能漏检，
  也可能把宿主内源病毒序列（EVM/ERV）算进去。定位是**便宜的预筛**，
  结论必须用下游比对 / 组装复核。
* CRR 不在 SRA 内，永远不会有 STAT 结果；GSA 也只有一部分 CRR 做过分析。
  拿不到分类结果的 run 会被单列（unresolved.tsv / .list），
  这批只能下载 reads 走本地比对（kraken2 / 平台病毒识别模块）。
* 想要「反向查询」（全库里有谁含某病毒），本接口做不到 ——
  那需要 BigQuery / Athena 上的 sra_tax_analysis_tool 表做 SQL 过滤。
""",
    )
    p.add_argument("--runs", help="accession 列表文件（每行一个）")
    p.add_argument("--csv", help="含 Run 列的 CSV/TSV（如平台检索产出的元数据表）")
    p.add_argument("--run-col", default="Run", help="CSV 中的 Run 列名（默认 Run）")
    p.add_argument("--acc", nargs="*", help="直接给 accession")
    p.add_argument("-o", "--outdir",
                   default=os.path.join(DIRS['meta_search'], 'prescreen_out'),
                   help="输出目录（默认在平台 run/meta_search/ 下）")
    p.add_argument("-j", "--jobs", type=int, default=6, help="并发数（默认 6，勿过高）")
    p.add_argument("--gsa-jobs", type=int, default=3,
                   help="GSA 并发数（NGDC 有反爬，默认 3，勿调高）")
    p.add_argument("--min-interval", type=float, default=0.0,
                   help="NCBI 最小请求间隔秒（0=不限）")
    p.add_argument("--gsa-min-interval", type=float, default=0.4,
                   help="GSA 最小请求间隔秒（默认 0.4，太低会 429）")
    p.add_argument("--retries", type=int, default=4, help="失败重试次数")
    p.add_argument("--timeout", type=int, default=90, help="单请求超时秒")
    p.add_argument("--virus-pct", type=float, default=None,
                   help="病毒占比阈值(%%); 达到则标记 pass 并导出 prescreen_pass.tsv")
    p.add_argument("--taxids", help="自定义目标 taxid 列表文件（每行一个），追加到 Viruses(10239)")
    p.add_argument("--taxdump", help="taxdump 目录（含 nodes.dmp）；默认用平台 databases/tax_db")
    p.add_argument("--no-taxdump", action="store_true",
                   help="不用 taxonomy，退回名称启发式判定（近似，不推荐）")
    p.add_argument("--cra-map", help="CRR->CRA 映射文件（TSV 两列），可跳过 GSA 网页解析")
    p.add_argument("--no-gsa", action="store_true", help="禁用 GSA 通路，只跑 NCBI")
    p.add_argument("--recheck-nodata", action="store_true",
                   help="忽略「无数据」标记，重新问一次 GSA（GSA 会陆续补做分析，可定期复查）")
    p.add_argument("--no-cache", action="store_true", help="忽略已有 JSON 缓存，强制重抓")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    accs = read_accessions(args)
    if not accs:
        UI.err("没有读到任何 accession。用 --runs / --csv / --acc 指定输入。")
        return 2

    os.makedirs(args.outdir, exist_ok=True)
    raw_dir = os.path.join(args.outdir, "stat_raw")
    gsa_raw_dir = os.path.join(args.outdir, "gsa_raw")

    # ---- 目标 taxid（--taxids 为追加项）----
    extra_taxids: Tuple[int, ...] = ()
    if args.taxids:
        extra: List[int] = []
        with open(args.taxids, "r", encoding="utf-8-sig") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    try:
                        extra.append(int(line.split()[0]))
                    except ValueError:
                        pass
        extra_taxids = tuple(extra)
        if extra_taxids:
            UI.info(f"追加目标 taxid: {extra_taxids}")

    # ---- 病毒判定器（默认复用平台 taxonomy）----
    judge: Optional[Callable[[Any], bool]] = None
    if not args.no_taxdump:
        tax_dir = args.taxdump or DIRS['taxonomy']
        try:
            judge = make_viral_judge(tax_dir, extra_taxids)
        except Exception as e:
            UI.warn(f"taxonomy 加载失败({e})，回退名称启发式")
            judge = None
        if judge is not None:
            UI.ok(f"已加载 taxonomy: {tax_dir}")
            UI.info("判定口径: taxonomy 子树 ∪ 名称启发式（互补；单用任一边都会漏）")
        else:
            UI.warn(f"taxonomy 不可用（{tax_dir}），回退名称启发式（近似）")
    else:
        UI.warn("--no-taxdump：病毒明细用名称启发式（近似；会漏掉新版 ICTV 高阶单元）")

    # ---- 按前缀路由 ----
    gsa_accs = [a for a in accs if a.upper().startswith("CRR")]
    stat_accs = [a for a in accs if a not in set(gsa_accs)]
    if args.no_gsa:
        stat_accs, gsa_accs = accs, []

    UI.banner("SRA / GSA 病毒样本预筛")
    UI.info(f"待处理 run 数 : {len(accs):,}  (NCBI STAT {len(stat_accs)} / GSA {len(gsa_accs)})")
    UI.info(f"并发          : NCBI {args.jobs} 线程 (间隔 {args.min_interval}s) | "
            f"GSA {args.gsa_jobs} 线程 (间隔 {args.gsa_min_interval}s)")
    UI.info(f"输出目录      : {os.path.abspath(args.outdir)}")
    if gsa_accs:
        UI.info("GSA 通路：无 tax_id，病毒按名字判定；分母是总 reads（见 virus_pct_ident）")

    cra_map: Dict[str, str] = {}
    if args.cra_map and os.path.exists(args.cra_map):
        with open(args.cra_map, "r", encoding="utf-8-sig", errors="replace") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    cra_map[parts[0]] = parts[1]
        UI.ok(f"已加载 CRR->CRA 映射 {len(cra_map)} 条")

    fetcher = Fetcher(raw_dir, retries=args.retries, timeout=args.timeout,
                      min_interval=args.min_interval)
    gsa = GsaFetcher(gsa_raw_dir, retries=args.retries, timeout=args.timeout,
                     min_interval=args.gsa_min_interval, cra_map=cra_map)
    gsa_sem = threading.Semaphore(max(1, args.gsa_jobs))

    summary_rows: List[Dict[str, Any]] = []
    viral_rows_all: List[Dict[str, Any]] = []
    status_count: Dict[str, int] = {}
    t_start = time.time()
    done = 0
    gsa_set = set(gsa_accs)

    def work(acc: str) -> Tuple[str, str, Optional[Dict[str, Any]], str, float, str, str]:
        if acc in gsa_set:
            with gsa_sem:
                a, payload, status, dt, cra, detail = gsa.fetch(
                    acc, use_cache=not args.no_cache, recheck_nodata=args.recheck_nodata)
            return a, "GSA", payload, status, dt, detail, cra
        a, payload, status, dt = fetcher.fetch(acc, use_cache=not args.no_cache)
        return a, "NCBI_STAT", payload, status, dt, "", ""

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(work, a): a for a in accs}
        for fut in as_completed(futs):
            acc, source, payload, status, dt, detail, cra = fut.result()
            done += 1
            key = f"{source}:{status.split(' ')[0]}"
            status_count[key] = status_count.get(key, 0) + 1

            if payload is None:
                summary_rows.append({
                    "acc": acc, "source": source, "status": status, "detail": detail,
                    "cra": cra,
                    "analysed": 0, "identified": 0, "identified_pct": 0.0,
                    "unidentified_pct": 0.0,
                    "n_taxa": 0, "virus_count": 0, "virus_pct": 0.0, "virus_pct_ident": 0.0,
                    "top_virus": "", "top_virus_pct": 0.0,
                    "top_organism": "", "top_organism_pct": 0.0,
                    "fetch_s": round(dt, 2), "pass": "",
                })
            else:
                r = (parse_gsa_tree(payload) if source == "GSA"
                     else parse_payload(payload, judge))
                vr = r.pop("_viral_rows", [])
                for x in vr:
                    viral_rows_all.append({"acc": acc, "source": source, **x})
                r["acc"] = acc
                r["source"] = source
                r["status"] = status
                r["detail"] = detail
                r["cra"] = cra
                r["fetch_s"] = round(dt, 2)
                if args.virus_pct is not None:
                    r["pass"] = "PASS" if r["virus_pct"] >= args.virus_pct else ""
                else:
                    r["pass"] = ""
                summary_rows.append(r)

            if done % 25 == 0 or done == len(accs):
                el = time.time() - t_start
                rate = done / el if el > 0 else 0
                eta = (len(accs) - done) / rate if rate > 0 else 0
                print(f"  {UI.C['gray']}[{done}/{len(accs)}] {rate:.2f} run/s "
                      f"ETA {eta/60:.1f} min{UI.C['reset']}")

    summary_rows.sort(key=lambda x: x["acc"])

    # ---------------- 输出 ----------------
    # virus_pct      : 各来源自身口径（STAT=/analysed；GSA=/总 reads）
    # virus_pct_ident: 归一化口径（病毒命中/已鉴定命中），**跨来源比较用这一列**
    cols = ["acc", "source", "cra", "status", "detail", "virus_pct", "virus_pct_ident",
            "virus_count", "top_virus", "top_virus_pct",
            "analysed", "identified", "identified_pct", "unidentified_pct", "n_taxa",
            "top_organism", "top_organism_pct", "pass", "fetch_s"]
    sum_path = os.path.join(args.outdir, "stat_run_summary.tsv")
    with open(sum_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in summary_rows:
            w.writerow({k: fmt_cell(v) for k, v in r.items()})

    vir_path = os.path.join(args.outdir, "stat_viral_taxa.tsv")
    # ⚠ 必须显式排序：viral_rows_all 是按**并发完成顺序**追加的，不排就
    # 会「同一输入两次跑，文件行序不同」——上游遗留的不确定性。
    # sorted 是稳定排序，acc 内部仍保持 parse 时的 -total_count 降序。
    viral_rows_sorted = sorted(viral_rows_all, key=lambda x: (x["acc"], x["source"]))
    with open(vir_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["acc", "source", "tax_id", "org", "path",
                                          "self_count", "total_count", "pct"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in viral_rows_sorted:
            w.writerow({k: fmt_cell(v) for k, v in r.items()})

    pass_path = None
    if args.virus_pct is not None:
        passed = [r for r in summary_rows if r.get("pass") == "PASS"]
        pass_path = os.path.join(args.outdir, "prescreen_pass.tsv")
        with open(pass_path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in sorted(passed, key=lambda x: -x["virus_pct"]):
                w.writerow({k: fmt_cell(v) for k, v in r.items()})

    # ---- 未决清单：没有分类数据的 run，需要走本地比对路线 ----
    unresolved = [r for r in summary_rows if r["status"] != "OK"]
    unrel_path = os.path.join(args.outdir, "unresolved.tsv")
    with open(unrel_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["acc", "source", "cra", "status", "detail"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(unresolved, key=lambda x: (x["source"], x["acc"])):
            w.writerow({k: fmt_cell(v) for k, v in r.items()})
    # 纯 accession 列表，直接可喂给下游下载/比对
    unrel_list = os.path.join(args.outdir, "unresolved.list")
    with open(unrel_list, "w", encoding="utf-8") as f:
        for r in sorted(unresolved, key=lambda x: x["acc"]):
            f.write(r["acc"] + "\n")

    el = time.time() - t_start
    ok_rows = [r for r in summary_rows if r["status"] == "OK"]
    with_virus = [r for r in ok_rows if r["virus_pct"] > 0]

    UI.stage("完成", "done")
    UI.ok(f"总耗时 {el/60:.1f} min，平均 {len(accs)/el:.2f} run/s" if el > 0 else "完成")
    UI.info(f"状态分布: {status_count}")
    UI.ok(f"有分类数据: {len(ok_rows)}/{len(accs)}")
    for src in ("NCBI_STAT", "GSA"):
        sub = [r for r in summary_rows if r["source"] == src]
        if sub:
            okn = sum(1 for r in sub if r["status"] == "OK")
            UI.info(f"  {src}: {okn}/{len(sub)} 有数据")
    UI.ok(f"检出病毒(>0%): {len(with_virus)}")
    if args.virus_pct is not None:
        UI.ok(f"通过阈值 {args.virus_pct}%: "
              f"{len([r for r in ok_rows if r['virus_pct'] >= args.virus_pct])}")
    UI.info(f"汇总   -> {sum_path}")
    UI.info(f"病毒明细 -> {vir_path}")
    if pass_path:
        UI.info(f"通过清单 -> {pass_path}")

    if unresolved:
        UI.warn(f"{len(unresolved)} 个 run 没有分类数据，已单列出来 -> {unrel_path}")
        UI.info(f"  （纯 accession 列表: {unrel_list}，可直接喂给下游本地比对）")
        by_status: Dict[str, int] = {}
        for r in unresolved:
            k = f"{r['source']}:{r['status'].split(' ')[0]}"
            by_status[k] = by_status.get(k, 0) + 1
        for k in sorted(by_status):
            UI.info(f"    {k}: {by_status[k]}")
        if any(r["status"].startswith("FAILED") for r in unresolved):
            UI.warn("其中有 FAILED —— 属于临时故障，**未写入缓存**，直接重跑即可补上")

    # 摘要报告
    rep = os.path.join(args.outdir, "STAT_PRESCREEN_REPORT.md")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("# SRA / GSA 病毒样本预筛报告\n\n")
        f.write(f"- 生成时间: {datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"- 输入 run 数: {len(accs)}\n")
        f.write(f"- 有分类数据: {len(ok_rows)}\n")
        f.write(f"- 检出病毒(>0%): {len(with_virus)}\n")
        f.write(f"- 状态分布: {status_count}\n")
        f.write(f"- 无分类数据: {len(unresolved)}（见 unresolved.tsv）\n")
        f.write(f"- NCBI 病毒判定: "
                f"{'taxonomy 子树 ∪ 名称启发式（互补）' if judge is not None else '名称启发式(近似，未用 taxonomy)'}\n")
        f.write("- GSA 病毒判定: 名称判定（GSA 分类树不含 tax_id）\n")
        if judge is not None:
            f.write("- ⚠ 为什么不是「纯 taxonomy」：STAT 的 tax_table 出自另一版\n")
            f.write("  taxonomy，其中部分 taxid（Viroids / Cytorhabdovirus 等）在平台\n")
            f.write("  tax_db 里已不存在；而新版 ICTV 高阶单元（Kitrinoviricota /\n")
            f.write("  Orthornavirae 等）名字不含 \"virus\"，名称启发式又会漏。\n")
            f.write("  两边单用都会漏，故取并集。\n")
        f.write("\n")
        if unresolved:
            f.write("## 无分类数据的 run（需本地比对）\n\n")
            f.write("这些 run 拿不到预计算分类结果，**不是分析失败**，"
                    "而是上游没有这份数据。要判断物种只能下载 reads 本地比对。\n\n")
            f.write(f"完整清单见 `unresolved.tsv`（含原因）与 `unresolved.list`"
                    f"（纯 accession，可直接喂下游）。共 {len(unresolved)} 条。\n\n")

            # 按项目聚合 —— 缺失几乎总是整项目缺失，聚合后一眼能看出该等哪个项目
            by_proj: Dict[str, List[str]] = {}
            for r in unresolved:
                by_proj.setdefault(r.get("cra") or "(未解析)", []).append(r["acc"])
            f.write("### 按项目（CRA）聚合\n\n")
            f.write("| CRA | 无数据 run 数 | 状态 |\n|---|---|---|\n")
            for k in sorted(by_proj, key=lambda x: -len(by_proj[x])):
                st = next((r["status"].split(" ")[0] for r in unresolved
                           if (r.get("cra") or "(未解析)") == k), "")
                f.write(f"| {k} | {len(by_proj[k])} | {st} |\n")
            f.write("\n")

            f.write("### 明细（前 40 条）\n\n")
            f.write("| run | 来源 | CRA | 状态 | 原因 |\n|---|---|---|---|---|\n")
            for r in sorted(unresolved, key=lambda x: (x["source"], x["acc"]))[:40]:
                f.write(f"| {r['acc']} | {r['source']} | {r.get('cra','')} | "
                        f"{r['status'].split(' ')[0]} | {r['detail']} |\n")
            if len(unresolved) > 40:
                f.write(f"\n_（共 {len(unresolved)} 条，此处仅列前 40）_\n")
            f.write("\n")
        f.write("## 口径警告\n\n")
        f.write("- `virus_pct` 是各来源**自身口径**：NCBI STAT 分母为 analysed（已鉴定 reads），\n")
        f.write("  GSA 分母为**总 reads**。直接比较会低估 GSA。\n")
        f.write("- 跨来源比较请用 `virus_pct_ident`（病毒命中 / 已鉴定命中）。\n")
        f.write("- 两者均为 **k-mer 快速分类**，不是严谨比对：低丰度病毒会漏检，\n")
        f.write("  宿主内源病毒序列会误报。**这是便宜的预筛，结论须用下游比对/组装复核。**\n\n")
        f.write("## 病毒占比 Top 30\n\n")
        f.write("| run | 来源 | 病毒% | 病毒%(已鉴定口径) | 主要病毒 | 测序量 | 主要非病毒物种 |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in sorted(with_virus, key=lambda x: -x["virus_pct_ident"])[:30]:
            f.write(f"| {r['acc']} | {r['source']} | {r['virus_pct']:.4f} | "
                    f"{r['virus_pct_ident']:.4f} | {r['top_virus']} | "
                    f"{r['analysed']:,} | {r['top_organism']} |\n")
    UI.info(f"报告   -> {rep}")

    # ---- 自动出图（失败只提示，不阻断预筛结果本身）----
    try:
        from .prescreen_plot import plot_prescreen      # 包内
    except ImportError:
        try:
            from prescreen_plot import plot_prescreen   # 独立脚本直跑
        except Exception:
            plot_prescreen = None
    if plot_prescreen is not None:
        try:
            figs = plot_prescreen(args.outdir, threshold=args.virus_pct)
            UI.ok(f"预筛图组 -> {os.path.join(args.outdir, 'figures')}（{len(figs)} 个文件）")
        except Exception as _e:                          # noqa: BLE001
            UI.warn(f"预筛图组生成失败（不影响结果表）: {_e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
