# -*- coding: utf-8 -*-
"""用 246 服务器上 MMPV-paper 的 sra.info.csv（真实采集元数据）跑 GCVA-32 地理迁移。"""
import csv
import io
import json
import os
import re
import shutil
import sys
import types
import threading

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\桌面\植物病毒分析平台")
os.chdir(r"D:\桌面\植物病毒分析平台")

import paramiko                                               # noqa: E402
HOST, USER, PWD = "202.119.189.246", "zhangwenda", "zwd123456"
REMOTE = ("/home/zhangwenda/MMPV-paper/goji-virome/01_data/sra_rna.data4/"
          "8.Virus_full.result/Cytorhabdovirus_sp_lycii_OR489165.1/sra.info.csv")

SRC = r"D:\桌面\重组分析\gcva_32.fasta"
WORK = os.path.abspath("analysis/gcva_phylogeo")
os.makedirs(WORK, exist_ok=True)
FASTA = os.path.join(WORK, "gcva_32.fasta")
shutil.copyfile(SRC, FASTA)

# ---- 1. 服务器取 sra.info.csv ----
c = paramiko.SSHClient()
c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
c.connect(HOST, username=USER, password=PWD, timeout=20,
          allow_agent=False, look_for_keys=False)
sftp = c.open_sftp()
LOCAL_INFO = os.path.join(WORK, "sra.info.csv")
sftp.get(REMOTE, LOCAL_INFO)
sftp.close()
c.close()
print("downloaded sra.info.csv:", os.path.getsize(LOCAL_INFO), "bytes")

info = {}
with open(LOCAL_INFO, encoding="utf-8-sig", newline="") as f:
    for r in csv.DictReader(f):
        run = (r.get("Run") or "").strip()
        if run:
            info[run] = r

# ---- 2. 建 meta.csv（region=省份段，year=CollectionDate 年份） ----
names = [l[1:].strip() for l in open(FASTA, encoding="utf-8") if l.startswith(">")]
rows, miss = [], []
for n in names:
    acc = n.split("_")[0]
    r = info.get(acc)
    if not r:
        miss.append(acc)
        rows.append((n, "", ""))
        continue
    loc = (r.get("Location") or "").strip()
    parts = [x.strip() for x in loc.split(",") if x.strip()]
    region = parts[1] if len(parts) >= 2 else (parts[0] if parts else "")
    m_y = re.search(r"(19|20)\d\d", r.get("CollectionDate") or "")
    year = m_y.group(0) if m_y else ""
    rows.append((n, region, year))
    if not region or not year:
        miss.append(acc)

with open(os.path.join(WORK, "gcva_32.meta.csv"), "w", encoding="utf-8",
          newline="") as f:
    f.write("seq_id,region,year\n")
    for n, r_, y in rows:
        # 空行不写：参考号 OR489165.1 是所有 run 名的子串，空行会被
        # 引擎的「包含匹配」抢先命中，把整表拖成 Unknown（实测踩中）
        if r_ and y:
            f.write(f"{n},{r_},{y}\n")
from collections import Counter
cnt = Counter((r_, y) for _, r_, y in rows if r_ and y)
print("meta：region+year 齐全", sum(1 for _, r_, y in rows if r_ and y), "/ 32")
print("分布:", dict(cnt))
print("缺:", miss)

# ---- 3. 系统地理（MOTP）+ GIF + 补图 ----
from Virus_Platform_Core.web.tool_jobs import _tool_job_phylogeo   # noqa: E402
ctx = types.SimpleNamespace(run_dir=WORK,
                            p={"input": FASTA,
                               "meta": os.path.join(WORK, "gcva_32.meta.csv"),
                               "motp_bin": 1})
ctx.req = lambda k, w=None: ctx.p[k]
ctx.opt = lambda k, d=None: ctx.p.get(k, d)
_tool_job_phylogeo(ctx)(lambda m: None, lambda *a: None, threading.Event())

from Virus_Platform_Core.toolbox_plots import add_plots            # noqa: E402
from Virus_Platform_Core import phylogeo_gif as g                  # noqa: E402
pd = os.path.join(WORK, "phylogeo")
sd = json.load(open(os.path.join(pd, "summary.json"), encoding="utf-8"))
try:
    print("GIF:", g.build_migration_gif(sd, os.path.join(pd, "migration.gif")))
except Exception as e:
    print("GIF 未出：", repr(e)[:140])
print("FIGS:", [p.split("gcva_phylogeo")[-1] for p in add_plots(WORK)])
print("dir:", sorted(os.listdir(pd)))
