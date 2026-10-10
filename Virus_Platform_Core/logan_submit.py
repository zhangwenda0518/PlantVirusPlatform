#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logan 序列查询提交自动化（kmviz 网页 DOM 直操作版）— 邮箱池轮换版
================================================================
- --email 支持逗号分隔多个邮箱, 每提交 N 条轮换下一个(分散限额风险)
- 每 RESTART_EVERY 条提交自动重启浏览器, 防 msedgedriver session 失效雪崩
- 轮询持续 HTTP 400 超 dead_after 秒判 session 死亡, 重启浏览器重提一次
其余逻辑同原版。用法:
    python logan_submit.py 输入.txt -o 结果目录 --email a@qq.com,b@qq.com

本文件已并入植物病毒分析平台（vp 包）：单文件自包含，仅依赖
标准库 + selenium（pip install selenium）+ 本机 Edge/Chrome，
可独立拷出平台任意目录使用；平台「LOGAN 溯源」批量模式由
Virus_Platform_Core.logan_trace.batch_submit 以子进程方式调用本脚本。
(并入时修复: result.csv 的 vis_url 引用未定义变量 sid 的 NameError;
 结果文件名统一 safe_stem 安全化 + _in_base 路径包含校验, 防穿越;
 关键节点输出机器可读进度行 "PROGRESS {json}" 供平台实时抓取)
"""
import argparse
import csv
import io
import logging
import os
import random
import re
import shutil
import sys
import time
import json
from pathlib import Path

# ---------------------------------------------------------------- 输入解析

_ACC_RE = re.compile(r"[A-Za-z]{1,2}_?[A-Za-z]{0,2}\d{4,9}(?:\.\d{1,2})?")
_SEQ_RE = re.compile(r"^[ACGTUNRYKMSWBDHVacgtunrykmswbdhv\-\s]+$")
_SYSRAND = random.SystemRandom()      # 邮箱轮换用系统级随机源


def safe_stem(acc: str) -> str:
    """结果文件名安全化：仅保留字母数字 . _ -，折叠连续点，防空防穿越。"""
    stem = re.sub(r"[^A-Za-z0-9._\-]", "_", str(acc))
    stem = re.sub(r"\.{2,}", ".", stem).strip(".")
    return stem or "query"


def _in_base(base: str, name: str) -> str:
    """规范化拼接到 base 下的路径；逃逸 base 即拒绝（防 ../ 穿越）。"""
    p = os.path.normpath(os.path.join(base, name))
    if os.path.isabs(name) or not p.startswith(base + os.sep):
        raise ValueError(f"非法输出文件名: {name}")
    return p


def emit_progress(**kw):
    """输出机器可读进度行 "PROGRESS {json}"（flush 保证管道实时可读；
    平台批量模式解析此行驱动进度条，人类阅读日志不受影响）。"""
    kw.setdefault('t', round(time.time(), 1))
    print('PROGRESS ' + json.dumps(kw, ensure_ascii=False), flush=True)


# curl：Windows 10+ 自带；缺失则下载链必挂，启动时探测一次快速失败（P3-17c）
_CURL = None


def probe_curl(log=None):
    """探测 curl（curl.exe 优先，其次 PATH 里的 curl）。返回可执行名或 None。"""
    global _CURL
    _CURL = shutil.which('curl.exe') or shutil.which('curl')
    if log:
        log.info("curl: %s", _CURL or "未找到")
    return _CURL


def _ascii_space_free(path) -> bool:
    """路径是否 ASCII 且不含空格（curl -o 落非 ASCII 路径时的中转判定）。"""
    s = str(path)
    return ' ' not in s and s.isascii()


def sanitize_seq(seq: str) -> str:
    """kmviz 实测只接受 ACGTN(其余 IUPAC 码触发 Alphabet is not 'dna',
    N 在比对区会被通配, 故统一把非 ACGTN 字母替换为 N)。"""
    import re as _re
    return _re.sub(r"[^ACGTNacgtn]", "N", seq).upper()


def parse_input(text: str):
    """解析输入文本 -> [(acc, seq_or_None), ...]"""
    entries = []
    lines = text.replace("\r\n", "\n").split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            acc = line[1:].split()[0] if len(line) > 1 else ""
            seq_parts = []
            while i < len(lines):
                s = lines[i].strip()
                if s.startswith(">"):
                    break
                if s and not _SEQ_RE.match(s):
                    break
                if s:
                    seq_parts.append(s)
                i += 1
            seq = "".join(seq_parts) if seq_parts else None
            if acc:
                entries.append((acc, seq))
            continue
        if "|" in line:
            left, right = line.split("|", 1)
            left, right = left.strip(), right.strip()
            if right and _SEQ_RE.match(right) and len(right) > 30:
                entries.append((left, right))
            elif left:
                entries.append((left, None))
            continue
        if "\t" in line:
            left, right = line.split("\t", 1)
            right = right.strip()
            if right and _SEQ_RE.match(right) and len(right) > 30:
                entries.append((left.strip(), right))
            else:
                entries.append((left.strip(), None))
            continue
        for tok in re.split(r"[,;\s]+", line):
            tok = tok.strip()
            if tok and _ACC_RE.fullmatch(tok):
                entries.append((tok, None))
    return entries


# ---------------------------------------------------------------- 浏览器自动化

def make_driver(headless: bool):
    from selenium import webdriver
    from selenium.webdriver.edge.options import Options as EdgeOptions
    from selenium.webdriver.chrome.options import Options as ChromeOptions

    errs = []
    for opt_cls, name in ((EdgeOptions, "edge"), (ChromeOptions, "chrome")):
        try:
            opts = opt_cls()
            if headless:
                opts.add_argument("--headless=new")
            # 瘦身: 禁图片/动画/扩展, 提速降耗
            opts.add_argument("--disable-gpu")
            opts.add_argument("--blink-settings=imagesEnabled=false")
            opts.add_argument("--disable-extensions")
            opts.add_argument("--disable-features=AnimateResizeWindow")
            opts.add_argument("--window-size=1920,1080")
            opts.add_argument("--log-level=3")
            opts.add_argument("--no-first-run")
            if name == "edge":
                return webdriver.Edge(options=opts), name
            return webdriver.Chrome(options=opts), name
        except Exception as e:
            errs.append(f"{name}: {e}")
    raise RuntimeError("无法启动浏览器(Edge/Chrome 都失败):\n" + "\n".join(errs))


KMVIZ_URL = "https://logan-search.org/dashboard"
DEFAULT_GROUP = "Fast_No_human"

SEL_TEXTAREA = '[id="{\\"index\\":\\"text\\",\\"type\\":\\"kmviz-input\\"}"]'
SEL_EMAIL = '[id="{\\"index\\":\\"notif-Email\\",\\"type\\":\\"kmviz-config-notif\\"}"]'
SEL_MODAL = "[role='dialog']"
SEL_GROUP_WRAP = ".kmviz-dmc-select-input-root .mantine-MultiSelect-wrapper"


def submit_one(driver, wait, acc: str, seq, retries: int, log, group=DEFAULT_GROUP, email=None):
    """返回 (status, message, session_url, session_time)。"""
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC

    sub_clicked = False
    for attempt in range(1, retries + 1):
        try:
            driver.get(KMVIZ_URL)
            modal = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, SEL_MODAL)))
            try:
                raw = driver.execute_script("return localStorage.getItem('user-sessions')")
                known_sessions = set(json.loads(raw)) if raw else set()
            except Exception:
                known_sessions = set()

            ta = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, SEL_TEXTAREA)))
            ta.clear()
            # seq 可能为 None（纯 accession 行 / 头部无序列）：sanitize_seq
            # 内部 re.sub(seq) 会抛 TypeError，被外层 except 吞成"第 N 次尝试
            # 失败"——先归一成空串，让 Logan-Search 按 accession 自己取序列。
            ta.send_keys(f">{acc}\n{sanitize_seq(seq or '')}")

            for b in modal.find_elements(By.TAG_NAME, "button"):
                if b.text.strip() == "Load":
                    b.click(); break
            time.sleep(2)

            if email:
                modal.find_element(By.CSS_SELECTOR, SEL_EMAIL).send_keys(email)

            wrap = modal.find_element(By.CSS_SELECTOR, SEL_GROUP_WRAP)
            wrap.click(); time.sleep(1.5)
            active = driver.switch_to.active_element
            active.send_keys(group); time.sleep(1.5)
            clicked = False
            for o in driver.find_elements(By.CSS_SELECTOR, "[role='option']"):
                if o.text.strip() == group:
                    try:
                        o.click()
                    except Exception:
                        driver.execute_script("arguments[0].click();", o)
                    clicked = True; break
            if not clicked:
                return "ERROR", f"Groups 选项 {group} 未找到", None, time.time()
            time.sleep(1)
            try:
                modal.find_element(By.TAG_NAME, "h1").click()
            except Exception:
                pass
            time.sleep(1)

            sub_btn = None
            dl = time.time() + 60
            while time.time() < dl:
                subs = [b for b in modal.find_elements(By.TAG_NAME, "button") if b.text.strip() == "Submit"]
                if subs and subs[0].get_attribute("disabled") is None:
                    sub_btn = subs[0]
                    break
                time.sleep(2)
            if not sub_btn:
                return "ERROR", "Submit 按钮持续 disabled", None, time.time()
            time.sleep(2)
            if sub_btn.get_attribute("disabled") is not None:
                return "ERROR", "Submit 按钮仍 disabled", None, time.time()
            try:
                sub_btn.click()
            except Exception:
                driver.execute_script("arguments[0].click();", sub_btn)
            sub_clicked = True

            deadline = time.time() + 600
            while time.time() < deadline:
                try:
                    raw = driver.execute_script("return localStorage.getItem('user-sessions')")
                    sids = json.loads(raw) if raw else []
                    new = [s for s in sids if s not in known_sessions]
                    if new:
                        known_sessions.update(new)
                        sid = new[-1]
                        return ("SUBMITTED", f"session={sid}",
                                f"https://logan-search.org/api/download/{sid}", time.time())
                except Exception as e:
                    log.debug("session 轮询异常(下轮重试): %s: %s", type(e).__name__, str(e)[:120])
                time.sleep(3)
            return "TIMEOUT", "600s 未见新 session (Submit 已点击, 可能已在服务器排队)", None, time.time()

        except Exception as e:
            msg = f"第 {attempt} 次尝试失败: {type(e).__name__}: {str(e)[:200]}"
            log.warning(msg)
            if attempt == retries:
                return "ERROR", msg, None, time.time()
            if sub_clicked:
                return "ERROR", msg + " | Submit 已点击, 不重试避免重复提交", None, time.time()
            time.sleep(3 * attempt)


def dom_status_probe(driver, sid: str, log):
    """读取浏览器当前页面(即刚提交那条的实时结果视图)的文本,
    判断服务端是否已渲染 FAILED / No match found。
    返回 "FAILED" / None(无法判断)。不主动填任何 session id,
    只被动读当前页面——页面停在哪个 session 就看哪个。"""
    try:
        if "logan-search" not in (driver.current_url or ""):
            return None
        txt = driver.execute_script("return document.body.innerText") or ""
        up = txt.upper()
        if "NO MATCH FOUND" in up:
            return "FAILED"
        return None
    except Exception as e:
        log.debug("DOM 结果页探测失败(忽略): %s: %s", type(e).__name__, str(e)[:100])
        return None


def wait_and_download(session_url_or_id, out_dir: Path, acc: str, log,
                      first_wait: int = 300, poll: int = 60, max_wait: int = 1800,
                      session_time: float = None, max_age: int = 10800,
                      dead_after: int = 900, driver_ref=None,
                      progress_cb=None) -> str:
    """等待结果并下载。返回: tsv 路径 / ""(未就绪) / "DL_FAIL" / "DEAD"(持续400判死)

    progress_cb(dict)：每次状态探测后回调（stage=waiting/downloading/
    downloaded），供平台实时抓取进度。预热期(first_wait)不再整段盲睡，
    改为每 30s 一次轻量 HEAD 探测——结果提前就绪即提前下载，总等待
    窗口(first_wait+max_wait)与原逻辑一致。"""
    import subprocess
    import zipfile

    def _cb(**extra):
        if progress_cb:
            try:
                progress_cb(extra)
            except Exception:
                pass

    sid = session_url_or_id.rstrip("/").split("/")[-1]
    url = f"https://logan-search.org/api/download/{sid}"
    stem = safe_stem(acc)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = os.path.normpath(os.path.abspath(str(out_dir)))
    zip_p = _in_base(base, f"{stem}.zip")
    tsv_p = _in_base(base, f"{stem}.tsv")
    # curl 落盘点：输出目录含非 ASCII（如平台根在中文路径下）时先写
    # %TEMP% 的 ASCII 临时文件，下载完再拷回，避免 curl 编码问题（P3-17d）
    dl_p = zip_p
    _tmp_zip = None
    if not _ascii_space_free(base):
        import tempfile
        fd, _tmp_zip = tempfile.mkstemp(prefix='vp_logan_', suffix='.zip')
        os.close(fd)
        dl_p = _tmp_zip

    def _ret(v):
        """返回前清理 %TEMP% 中转 zip（非 ASCII 落点专用，P3-17d）。"""
        if _tmp_zip and os.path.isfile(_tmp_zip):
            try:
                os.remove(_tmp_zip)
            except OSError:
                pass
        return v
    start_t = time.time()
    if driver_ref is None:
        driver_ref = [None]

    age_deadline = (session_time or time.time()) + max_age
    warm_until = start_t + first_wait          # 预热期：30s 间隔提前探测
    deadline = min(start_t + first_wait + max_wait, age_deadline)
    last_probe_t = 0.0  # DOM 结果页探测节流: 每 5 分钟看一眼当前提交的实时结果
    if first_wait:
        log.info("轮询 %s（前 %ds 每 30s 探测一次，就绪即下载）", sid, first_wait)
    while time.time() < deadline:
        try:
            r = subprocess.run(
                [_CURL or "curl.exe", "-s", "-o", os.devnull,
                 "-w", "%{http_code}", "-I", "-L", url],
                capture_output=True, text=True, timeout=60)
            code = r.stdout.strip()
        except Exception as e:
            log.debug("HEAD 探测异常: %s", e)
            code = "ERR"
        if code == "200":
            _cb(stage='downloading', sid=sid)
            log.info("%s 就绪 (HTTP 200), 开始下载", sid)
            try:
                r = subprocess.run([_CURL or "curl.exe", "-L", "--fail",
                                    "-o", dl_p, url],
                                   capture_output=True, timeout=600)
                if r.returncode != 0:
                    log.warning("%s 下载失败 rc=%d", sid, r.returncode)
                    return _ret("DL_FAIL")
            except Exception as e:
                log.warning("%s 下载异常: %s", sid, e)
                return _ret("DL_FAIL")
            if dl_p != zip_p and os.path.isfile(dl_p):
                # 非 ASCII 落点 → 拷回真实产物路径（P3-17d）
                shutil.copyfile(dl_p, zip_p)
            if not os.path.isfile(zip_p) or os.path.getsize(zip_p) == 0:
                log.warning("%s 下载文件为空/缺失", sid)
                return _ret("DL_FAIL")
            tsv_names = []
            try:
                with zipfile.ZipFile(zip_p) as z:
                    tsv_names = [nm for nm in z.namelist() if nm.endswith(".tsv")]
                    if not tsv_names:
                        log.warning("%s zip 内无 tsv", sid)
                        return _ret("DL_FAIL")
                    Path(tsv_p).write_bytes(z.read(tsv_names[0]))
                    for j, nm in enumerate(tsv_names[1:], 1):
                        alt_p = _in_base(base, f"{stem}_{j}.tsv")
                        Path(alt_p).write_bytes(z.read(nm))
            except zipfile.BadZipFile:
                log.warning("%s zip 损坏", sid)
                return _ret("DL_FAIL")
            _cb(stage='downloaded', sid=sid, tsv=os.path.basename(tsv_p))
            return _ret(tsv_p)
        if time.time() >= age_deadline:
            break
        if code == "400" and time.time() - start_t >= dead_after:
            log.warning("%s 持续 HTTP 400 超 %ds, 判定 session 失效", sid, dead_after)
            return _ret("DEAD")
        # 每 5 分钟用浏览器看一眼「刚提交的这条」的实时结果页, 渲染出 FAILED 即提前判死
        if driver_ref[0] is not None and time.time() - last_probe_t >= 300:
            last_probe_t = time.time()
            st = dom_status_probe((driver_ref[0]() if callable(driver_ref[0]) else driver_ref[0]), sid, log)
            if st == "FAILED":
                log.warning("%s 结果页渲染 FAILED/No match found, 提前判定为阴性", sid)
                return _ret("DEAD")
        _cb(stage='waiting', sid=sid, waited=int(time.time() - start_t),
            http=code, span=first_wait + max_wait)
        wait_step = 30 if time.time() < warm_until else poll
        log.info("%s 未完成 (HTTP %s), %ds 后重试", sid, code, wait_step)
        time.sleep(wait_step)
    if time.time() >= age_deadline:
        log.warning("%s 距提交已超 3 小时未完成, 判定失败", sid)
    else:
        log.info("%s 本轮等满 %ds 未完成, 先继续下一条 (稍后补漏)", sid, max_wait)
    return _ret("")


# ---------------------------------------------------------------- 主流程

def main():
    ap = argparse.ArgumentParser(description="Logan 序列查询批量提交 (kmviz 自动化, 邮箱池轮换)")
    ap.add_argument("input", help="输入 txt 文件")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--headless", action="store_true")
    ap.add_argument("--delay", type=float, default=1.0)
    ap.add_argument("--group", default=DEFAULT_GROUP,
                    choices=["All", "All_No_viral_human", "Fast", "Fast_No_human", "Fast_No_RefSeq",
                             "Transcriptomic", "Metatranscriptomic", "Metagenomic", "GenBank_RefSeq"])
    ap.add_argument("--email", default="",
                    help="通知邮箱, 逗号分隔多个则**每次提交随机轮换**"
                         "（避免连续重复；必填；不再内置任何个人默认邮箱）")
    ap.add_argument("--email-rotate", type=int, default=5,
                    help="[已弃用/不生效] 历史参数：现实现为每条随机轮换邮箱，"
                         "不做「每 N 条顺序轮换」；保留仅为兼容旧命令行")
    ap.add_argument("--first-wait", type=int, default=300)
    ap.add_argument("--max-wait", type=int, default=1800)
    args = ap.parse_args()

    emails = [e.strip() for e in args.email.split(",") if e.strip()]
    if not emails:
        print("错误: --email 不能为空", file=sys.stderr)
        sys.exit(1)

    in_path = Path(args.input)
    if not in_path.exists():
        print(f"错误: 输入文件不存在: {in_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = (Path(args.outdir) if args.outdir
               else in_path.parent / "logan_submit_out").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base = os.path.normpath(os.path.abspath(str(out_dir)))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(base, "run.log"), mode="w",
                                encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    log = logging.getLogger("logan_submit")
    if len(emails) > 1:
        # 2026-10-03 修复 P3-17a：实际是「每条随机轮换」，不再谎称每 N 条
        log.info("邮箱池: %d 个, 每次提交随机轮换（避免连续重复）", len(emails))
    if probe_curl(log) is None:
        log.error("未找到 curl（curl.exe / PATH 中的 curl）：Logan 结果下载"
                  "依赖 curl，请先启用 Windows 自带 curl 或安装 curl 后重试")
        sys.exit(3)

    entries = parse_input(in_path.read_text(encoding="utf-8-sig", errors="replace"))
    total = len(entries)
    entries = [(a, s) for a, s in entries
               if not os.path.isfile(_in_base(base, f"{safe_stem(a)}.tsv"))]
    log.info("解析到 %d 条待查询记录, 断点续跑跳过 %d 条已有结果, 实际待跑 %d 条",
             total, total - len(entries), len(entries))
    if not entries:
        log.error("没有可跑的记录, 退出")
        sys.exit(2)

    driver, browser = make_driver(args.headless)
    from selenium.webdriver.support.ui import WebDriverWait
    wait = WebDriverWait(driver, 20)
    log.info("浏览器启动: %s, headless=%s", browser, args.headless)
    submit_count = 0
    RESTART_EVERY = 10
    total_submitted = 0  # 邮箱轮换用跨重启累计计数, 不随浏览器重启清零

    def restart_driver(reason=""):
        nonlocal driver, wait, submit_count
        log.info("重启浏览器 (%s)...", reason or "定期保养")
        try:
            driver.quit()
        except Exception:
            pass
        time.sleep(5)
        driver, browser = make_driver(args.headless)
        wait = WebDriverWait(driver, 20)
        submit_count = 0

    def pick_email():
        if len(emails) == 1:
            return emails[0]
        # 随机轮换邮箱, 分散单邮箱频率; 避免连续重复同一邮箱
        for _ in range(5):
            c = _SYSRAND.choice(emails)
            if c != pick_email.last:
                pick_email.last = c
                return c
        return c
    pick_email.last = None

    results = []
    unconfirmed = []          # 提交后未见 session（UNKNOWN）：纳入兜底回收与汇总
    try:
        pending = []
        for n, (acc, seq) in enumerate(entries, 1):
            if submit_count >= RESTART_EVERY:
                restart_driver("已提交 %d 条" % submit_count)
            email = pick_email()
            log.info("[%d/%d] 提交 %s (email=%s)", n, len(entries), acc, email)
            emit_progress(stage='submit', acc=acc, done=n - 1,
                          total=len(entries), msg=f'正在提交 {acc}')
            t0 = time.time()

            def _pcb(d, _acc=acc, _n=n, _total=len(entries)):
                d.update({'acc': _acc, 'done': _n - 1, 'total': _total,
                          'span': args.first_wait + args.max_wait})
                emit_progress(**d)

            status, msg, url, sess_t = submit_one(driver, wait, acc, seq, args.retries, log,
                                                   group=args.group, email=email)
            submit_count += 1
            total_submitted += 1
            sid = (url or "").rstrip("/").split("/")[-1]
            emit_progress(stage='submitted' if status == 'SUBMITTED' else 'submit_fail',
                          acc=acc, sid=sid, done=n - 1, total=len(entries),
                          msg=msg)
            rec = {"accession": acc, "status": status, "message": msg,
                   "session_url": url or "",
                   "vis_url": (f"https://logan-search.org/dashboard/{sid}"
                               if sid else ""),
                   "tsv": "", "submit_time": t0,
                   "session_time": sess_t}
            if status == "SUBMITTED" and url:
                tsv = wait_and_download(url, out_dir, acc, log,
                                        first_wait=args.first_wait, max_wait=args.max_wait,
                                        session_time=sess_t, driver_ref=[lambda: driver],
                                        progress_cb=_pcb)
                if tsv == "DEAD":
                    restart_driver("%s session 失效" % acc)
                    email = pick_email()
                    status, msg, url, sess_t = submit_one(driver, wait, acc, seq, args.retries, log,
                                                           group=args.group, email=email)
                    submit_count += 1
                    total_submitted += 1
                    sid = (url or "").rstrip("/").split("/")[-1]
                    rec.update({"status": status, "message": msg + " | 原session失效后重提",
                                "session_url": url or "",
                                "vis_url": (f"https://logan-search.org/dashboard/{sid}"
                                            if sid else ""),
                                "session_time": sess_t})
                    tsv = wait_and_download(url, out_dir, acc, log,
                                            first_wait=args.first_wait, max_wait=args.max_wait,
                                            session_time=sess_t, driver_ref=[lambda: driver],
                                            progress_cb=_pcb) if url else ""
                if tsv == "DEAD":
                    tsv = ""
                if tsv and tsv != "DL_FAIL":
                    rec["tsv"] = tsv
                    log.info("[%d/%d] %s 结果已下载", n, len(entries), acc)
                    emit_progress(stage='segment_done', acc=acc, done=n,
                                  total=len(entries), msg=f'{acc} 结果已下载')
                else:
                    if tsv == "DL_FAIL":
                        rec["message"] += " | 就绪但下载失败, 稍后补漏"
                    log.info("[%d/%d] %s 挂起, 继续下一条", n, len(entries), acc)
                    pending.append(rec)
            elif status == "TIMEOUT":
                rec["status"] = "UNKNOWN"
                log.warning("[%d/%d] %s 提交后未确认 session, 不重提", n, len(entries), acc)
                # 2026-10-03 修复 P3-17b：Submit 已点击但 600s 内未见新 session
                # 不等于失败（可能在服务器排队）——纳入兜底回收清单并在汇总输出，
                # 不再静默丢在结果表里。
                unconfirmed.append(rec)
            results.append(rec)
            time.sleep(args.delay)

        # 会话兜底回收：末尾再扫一次 localStorage，把未归属的新 session 作为
        # 候选下载链接列出（供人工核对），并写入对应条目的 message
        if unconfirmed:
            cand = []
            try:
                raw = driver.execute_script(
                    "return localStorage.getItem('user-sessions')")
                sids = json.loads(raw) if raw else []
                known = {r["session_url"].rstrip("/").split("/")[-1]
                         for r in results if r.get("session_url")}
                cand = [s for s in sids if s not in known]
            except Exception as e:
                log.debug("session 兜底扫描失败(忽略): %s: %s",
                          type(e).__name__, str(e)[:120])
            log.warning("提交后未确认 session 的条目 %d 条（可能已在服务器排队，"
                        "本轮未重提）：%s", len(unconfirmed),
                        ", ".join(r["accession"] for r in unconfirmed))
            if cand:
                cand_urls = ["https://logan-search.org/api/download/%s" % s
                             for s in cand]
                log.warning("兜底扫描发现 %d 个未归属 session（候选下载链接，"
                            "需人工核对归属）：%s", len(cand), ", ".join(cand_urls))
                for rec in unconfirmed:
                    rec["message"] += (" | 候选 session（未确认归属，需人工核对）: "
                                       + ", ".join(cand_urls))

        round_no = 1
        while pending:
            round_no += 1
            expired = [r for r in pending if time.time() - r["session_time"] > 10800]
            for rec in expired:
                rec["status"] = "FAILED"
                rec["message"] += " | 超3小时未完成, 判定失败"
                log.warning("%s 超 3 小时未完成, 判定失败", rec["accession"])
            pending = [r for r in pending if r["status"] != "FAILED"]
            if not pending:
                break
            log.info("===== 补漏轮 %d: %d 条挂起 =====", round_no, len(pending))
            still = []
            for rec in pending:
                log.info("补查 %s", rec["accession"])
                tsv = wait_and_download(rec["session_url"], out_dir, rec["accession"], log,
                                        first_wait=0, max_wait=args.max_wait,
                                        session_time=rec["session_time"])
                if tsv and tsv not in ("DL_FAIL", "DEAD"):
                    rec["tsv"] = tsv
                    rec["status"] = "SUBMITTED"
                    log.info("%s 结果已下载", rec["accession"])
                elif tsv == "DL_FAIL":
                    time.sleep(30)
                    tsv2 = wait_and_download(rec["session_url"], out_dir, rec["accession"], log,
                                             first_wait=0, poll=30, max_wait=120,
                                             session_time=rec["session_time"])
                    if tsv2 and tsv2 not in ("DL_FAIL", "DEAD"):
                        rec["tsv"] = tsv2
                        rec["status"] = "SUBMITTED"
                        log.info("%s 重试下载成功", rec["accession"])
                    else:
                        still.append(rec)
                elif tsv == "DEAD":
                    rec["status"] = "FAILED"
                    rec["message"] += " | 补漏轮 session 失效(400), 判定失败"
                    log.warning("%s 补漏轮 session 失效, 判定失败", rec["accession"])
                else:
                    still.append(rec)
            if len(still) == len(pending):
                time.sleep(300)
            pending = still
            if pending:
                log.warning("补漏轮后仍有 %d 条未落定, 直接判 FAILED 收尾(防死循环)", len(pending))
                for rec in pending:
                    rec["status"] = "FAILED"
                    rec["message"] += " | 补漏轮穷尽后仍无结果, 判定失败(需人工复查)"
                break
    finally:
        try:
            driver.quit()
        except Exception:
            pass

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=["accession", "status", "message",
                                        "session_url", "vis_url", "tsv",
                                        "submit_time", "session_time"])
    w.writeheader()
    w.writerows(results)
    Path(_in_base(base, "result.csv")).write_text(buf.getvalue(), encoding="utf-8-sig")

    ok = sum(1 for r in results if r["tsv"])
    unk = [r for r in results if r["status"] == "UNKNOWN"]
    if unk:
        log.warning("未确认 session（UNKNOWN）%d 条——已计入 result.csv，"
                    "请人工到 logan-search.org 核对或重跑批量补漏: %s",
                    len(unk), ", ".join(r["accession"] for r in unk))
    log.info("完成: %d/%d 成功, 结果: %s", ok, len(results),
             os.path.join(base, "result.csv"))


if __name__ == "__main__":
    main()
