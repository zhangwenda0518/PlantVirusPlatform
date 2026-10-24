# -*- coding: utf-8 -*-
"""sRNA 宿主去除（blastn-short 路径）。

背景（2026-10-18 基准实证）：平台现有两条宿主去除路径对 21–24nt 小RNA
读长均无效——kunpeng k=35 固定大于读长、bwa-mem 种子机制不覆盖；
blastn-short（BLAST ≤30nt 专用任务）可正常比对。实测 Bari S3/S9 分别
去除 46.8% / 66.1% 宿主 reads，病毒检出不受影响（宿主 reads 本就不映射
病毒库），收益在组装轨降噪与少量 FP 消除。

用法:
  python scripts/srna_host_deplete.py --reads in.fastq[.gz] \
      --host-genome host.fa[.gz] --out-dir <dir> [--threads 8] [--evalue 1e-3]
      [--paired R2.fastq[.gz]]

产物: <out-dir>/kept_R1.fastq(.gz)（及 kept_R2）+ deplete_stats.json
"""
import argparse
import gzip
import itertools
import json
import os
import subprocess
import sys
import time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BLASTN = os.path.join(REPO, "3rd", "tools", "Blast", "bin", "blastn.exe")
MKDB = os.path.join(REPO, "3rd", "tools", "Blast", "bin", "makeblastdb.exe")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    if r.returncode != 0:
        raise RuntimeError(f"命令失败 rc={r.returncode}: {' '.join(cmd)}\n"
                           f"{(r.stderr or '')[-500:]}")


def ensure_plain(src, dst):
    if src.endswith(".gz"):
        if os.path.isfile(dst) and os.path.getsize(dst) > 0:
            return dst
        log(f"解压 {src}")
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            while True:
                b = fi.read(1 << 24)
                if not b:
                    break
                fo.write(b)
        return dst
    return src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", required=True, help="R1 fastq(.gz)")
    ap.add_argument("--paired", help="R2 fastq(.gz)（可选）")
    ap.add_argument("--host-genome", required=True, help="宿主基因组 fasta(.gz)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", default="1e-3")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    db_tag = "".join(ch for ch in os.path.basename(a.host_genome)
                     if ch.isalnum())[:24] or "hostdb"
    db = os.path.join(a.out_dir, f"db_{db_tag}")
    genome = ensure_plain(a.host_genome, os.path.join(a.out_dir, "host.fa"))
    if not os.path.isfile(db + ".nsq"):
        log("makeblastdb")
        run([MKDB, "-in", genome, "-dbtype", "nucl", "-out", db])

    reads = ensure_plain(a.reads, os.path.join(a.out_dir, "reads_R1.fastq"))
    fq = {os.path.basename(a.reads): reads}
    if a.paired:
        fq[os.path.basename(a.paired)] = ensure_plain(
            a.paired, os.path.join(a.out_dir, "reads_R2.fastq"))

    host_names = set()
    stats = {"db": db, "files": {}}
    for tag, rp in fq.items():
        hits = os.path.join(a.out_dir, tag + ".hits.txt")
        log(f"blastn-short: {tag}")
        run([BLASTN, "-task", "blastn-short", "-query", rp, "-db", db,
             "-outfmt", "6 qseqid", "-evalue", a.evalue,
             "-num_threads", str(a.threads), "-dust", "no", "-out", hits])
        with open(hits) as f:
            host_names |= {l.strip() for l in f if l.strip()}

    kept = {t: os.path.join(a.out_dir, "kept_" + os.path.basename(p))
            for t, p in fq.items()}
    n_in = n_host = 0
    writers = {t: open(p, "w", encoding="utf-8") for t, p in kept.items()}
    srcs = list(fq.items())
    for name, rp in srcs:
        opener = gzip.open if rp.endswith(".gz") else open
        mode = "rt" if rp.endswith(".gz") else "r"
        with opener(rp, mode) as f:
            for h, s, plus, q in itertools.zip_longest(*[iter(f)] * 4):
                if h is None:
                    break
                rid = h[1:].strip().split()[0]
                if name == srcs[0][0]:
                    n_in += 1
                if rid in host_names:
                    if name == srcs[0][0]:
                        n_host += 1
                    continue
                writers[name].write(f"@{rid}\n{s.strip()}\n+\n{q.strip()}\n")
    for w in writers.values():
        w.close()

    stats.update({"input": n_in, "host_removed": n_host,
                  "kept": n_in - n_host,
                  "host_pct": round(100.0 * n_host / max(n_in, 1), 2),
                  "kept_files": kept})
    with open(os.path.join(a.out_dir, "deplete_stats.json"), "w",
              encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=1)
    log(f"输入 {n_in} | 去除宿主 {n_host}（{stats['host_pct']}%）| "
        f"kept {n_in-n_host} -> {a.out_dir}")


if __name__ == "__main__":
    main()
