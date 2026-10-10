#!/usr/bin/env python3
"""管线 E1/E2/E4/E5 探针 (2026-09-16)

在临时 work_dir 上直接调 virome_phylo_pipeline 的三个函数, 断言:

E1  report_builder: PDF 产物 → <embed class="pdf-embed">; 超限/不可渲染 → 显式 ⚠ 而非静默丢
E2  report_builder: 声明了但磁盘上没有的 summary → missing_declared_summaries 列出;
                    stage 未完成时不误报
E4  spread3_viz: 离线地图零外链 (external_refs == []) + strict 下 CDN 不存在; Leaflet 仅在 offline=False
E5  spread3_viz: 未知地点 → latitude=None + coord_missing; build_phylogeo_report 不崩

运行: PYTHONIOENCODING=utf-8 C:/Python312/python.exe tests/_probe_pipeline_e145.py
"""
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

# 外部参照仓（延伸基因组/MMPV-RNA 工作区）：默认仍是维护者本机路径，
# 可用环境变量 VP_EXT_ROOT 覆盖（指向 MMPV-RNA 工作区根）。
EXT_ROOT = os.environ.get('VP_EXT_ROOT', r'D:\桌面\延伸基因组\MMPV-RNA')
PIPE = Path(EXT_ROOT) / 'virome_phylo_pipeline'
if not PIPE.is_dir():
    print(f'[SKIP] 外部参照仓不存在: {PIPE}'
          '（如已迁移，请设 VP_EXT_ROOT 指向 MMPV-RNA 工作区根）')
    sys.exit(0)
sys.path.insert(0, str(PIPE))
sys.path.insert(0, str(PIPE / "utils"))

FAIL = []
OK = []


def check(name, cond, detail=""):
    (OK if cond else FAIL).append(name)
    print(("  PASS  " if cond else "  FAIL  ") + name + (f"  | {detail}" if detail else ""))


import report_builder as rb            # noqa: E402
import spread3_viz as s3               # noqa: E402
from utils.self_contained_html import external_refs   # noqa: E402

tmp = Path(tempfile.mkdtemp(prefix="probe_e145_"))
try:
    # ═══════════════ E1: PDF / 不可渲染产物 ═══════════════
    print("\n[E1] PDF 入报告 + 未内联显式化")
    vdir = tmp / "DemoVirus"
    (vdir / "geography" / "geo_analysis").mkdir(parents=True)
    (vdir / "qc").mkdir(parents=True)

    pdf_ok = vdir / "geography" / "geo_analysis" / "spatial_map.pdf"
    pdf_ok.write_bytes(b"%PDF-1.4\n" + b"x" * 2000)
    # 唯一 stem 的 PDF: 不会被同名 png 顶掉, 必须真正以 <embed> 内联
    pdf_only = vdir / "geography" / "geo_analysis" / "route_timeline.pdf"
    pdf_only.write_bytes(b"%PDF-1.4\n" + b"x" * 3000)
    pdf_big = vdir / "geography" / "geo_analysis" / "huge_map.pdf"
    pdf_big.write_bytes(b"%PDF-1.4\n" + b"x" * (rb.MAX_PDF_KB * 1024 + 5000))
    png_dup = vdir / "geography" / "geo_analysis" / "spatial_map.png"
    png_dup.write_bytes(b"\x89PNG\r\n\x1a\n" + b"y" * 500)
    bad_ext = vdir / "qc" / "notes.txt"
    bad_ext.write_text("not a plot", encoding="utf-8")

    skips = []
    pats = ["geography/geo_analysis/*.pdf", "geography/geo_analysis/*.png", "qc/*.txt"]
    got = rb.collect_plots(str(vdir), pats, skipped=skips)
    names = sorted(Path(p).name for p in got)
    check("E1 collect_plots 收录 pdf 且同 stem 取 png",
          names == ["huge_map.pdf", "route_timeline.pdf", "spatial_map.png"], f"got={names}")
    check("E1 低优先级同名单 + 非白名单扩展名 都进 skipped",
          any("spatial_map.pdf" in p for p, _ in skips) and
          any("notes.txt" in p for p, _ in skips),
          f"skips={[Path(p).name for p, _ in skips]}")

    card_pdf, why_pdf = rb.asset_card(str(pdf_ok), str(tmp))
    check("E1 pdf → <embed class=\"pdf-embed\">", '<embed class="pdf-embed"' in card_pdf and not why_pdf,
          why_pdf or "ok")
    card_big, why_big = rb.asset_card(str(pdf_big), str(tmp))
    check("E1 超限 pdf → 不内联但显式给原因 + 下载链接",
          bool(why_big) and "⚠" in card_big and "download" in card_big, why_big)
    card_txt, why_txt = rb.asset_card(str(bad_ext), str(tmp))
    check("E1 非白名单扩展名 → 显式原因", bool(why_txt), why_txt)
    card_html_fake = vdir / "qc" / "interactive.html"
    card_html_fake.write_text("<html>ok</html>", encoding="utf-8")
    card_h, why_h = rb.asset_card(str(card_html_fake), str(tmp))
    check("E1 .html 产物 → 链接卡 (不假装内联)",
          "在新标签打开" in card_h and not why_h, why_h or "ok")

    # ═══════════════ E2: 声明产物缺失护栏 ═══════════════
    print("\n[E2] 声明但未找到产物 → 显式列出")
    cnt_csv = vdir / "geography" / "geo_analysis" / "spatiotemporal_counts.csv"
    cnt_csv.write_text("region,year,count\nBJ,2020,3\n", encoding="utf-8")
    fake_mod = {"dir": "geography/geo_analysis", "stages": ["virspacetime"],
                "summaries": ["spatiotemporal_counts.csv", "missing_one.csv"]}
    miss = rb.missing_declared_summaries(str(vdir), fake_mod, ["virspacetime"])
    check("E2 存在的不报 / 缺失的报出", miss == ["geography/geo_analysis/missing_one.csv"], f"{miss}")
    check("E2 stage 未完成 → 不误报",
          rb.missing_declared_summaries(str(vdir), fake_mod, []) == [], "")
    # 用真 MODULES 声明做一次契约核对: 报告 glob 与产出端路径必须一致
    geo_mod = [m for m in rb.MODULES if m["dir"] == "geography"][0]
    miss_real = rb.missing_declared_summaries(str(vdir), geo_mod, ["virspacetime", "geo"])
    check("E2 真声明能命中 geo_analysis/spatiotemporal_counts.csv (路径对齐)",
          not any("spatiotemporal_counts" in x for x in miss_real), f"missing={miss_real}")
    check("E2 真声明下其余未产出项被列出", len(miss_real) >= 1, f"missing={miss_real}")
    (vdir / "geography" / "geo_analysis" / "missing_one.csv").write_text("a\n1\n", encoding="utf-8")
    miss2 = rb.missing_declared_summaries(str(vdir), fake_mod, ["virspacetime"])
    check("E2 全部补齐 → 空 (不刷屏)", miss2 == [], f"{miss2}")

    # ═══════════════ E4/E5: spread3_viz ═══════════════
    print("\n[E4/E5] 坐标 fail-loud + 离线地图零外链")
    locs = ["Beijing", "Ningxia", "AtlantisDeep"]   # 第三个不在内置表 → 必须缺坐标
    n_pairs = len(locs) * (len(locs) - 1)
    cols = ["location.indicators.%s.%s" % (a, b)
            for a in locs for b in locs if a != b]
    assert len(cols) == n_pairs
    log_path = tmp / "phylogeo.log"
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("state\t" + "\t".join(cols) + "\n")
        for i in range(10):
            # Beijing→Ningxia 后验 0.95 (显著), Ningxia→Beijing 0.2 (不显著), 其余 0.05
            vals = []
            for c in cols:
                if c.endswith("Beijing.Ningxia"):
                    vals.append("0.95")
                elif c.endswith("Ningxia.Beijing"):
                    vals.append("0.2")
                else:
                    vals.append("0.05")
            f.write(f"{i}\t" + "\t".join(vals) + "\n")

    r = s3.generate_spread3_json(str(log_path), None, locs, str(tmp / "out"))
    check("E5 解析出 6 条路线", r.get("n_routes") == 6, f"n_routes={r.get('n_routes')} err={r.get('error')}")
    check("E5 未知地点进 coord_missing", r.get("coord_missing") == ["AtlantisDeep"],
          f"{r.get('coord_missing')}")

    jp = r.get("json_path")
    check("E5 JSON 已生成", bool(jp) and os.path.exists(jp or ""), str(r.get("error")))
    data = json.load(open(jp, encoding="utf-8"))
    by_id = {n["id"]: n for n in data["nodes"]}
    check("E5 未知地点 lat/lon = None (不伪造)",
          by_id["AtlantisDeep"]["latitude"] is None and by_id["AtlantisDeep"]["longitude"] is None,
          str(by_id["AtlantisDeep"]))
    check("E5 已知地点坐标来自内置表",
          abs(by_id["Beijing"]["latitude"] - 39.9) < 0.2, str(by_id["Beijing"]))
    check("E5 metadata 带 coord_missing",
          data["metadata"].get("coord_missing") == ["AtlantisDeep"],
          str(data["metadata"].get("coord_missing")))
    check("E5 links 带 source/target 整数索引 (P2-11 回归护栏)",
          all(isinstance(l.get("source"), int) and isinstance(l.get("target"), int)
              for l in data["links"]), str(data["links"][:1]))

    html_off = s3.generate_interactive_map_html(jp, str(tmp / "map_offline.html"), offline=True)
    htxt = Path(html_off).read_text(encoding="utf-8")
    refs = external_refs(htxt)
    check("E4 离线地图零外链", refs == [], f"refs={refs[:4]}")
    check("E4 离线地图无 leaflet/unpkg script",
          "unpkg.com/leaflet" not in htxt and "cdn" not in htxt.lower().split("</style>")[-1][:400].lower(),
          "")
    check("E4 离线地图自绘 SVG", "<svg" in htxt and "自绘 SVG" in htxt, "")
    check("E5 离线地图显式标注缺坐标地点",
          "AtlantisDeep" in htxt and "没有坐标" in htxt and "不编造坐标" in htxt, "")
    check("E4 离线地图含显著路线与不显著路线区分",
          "路线后验 BF" in htxt or "显著" in htxt, "")

    html_on = s3.generate_interactive_map_html(jp, str(tmp / "map_online.html"), offline=False)
    otxt = Path(html_on).read_text(encoding="utf-8")
    check("E4 online 分支保留 Leaflet (opt-in)",
          "leaflet" in otxt.lower() and len(external_refs(otxt)) > 0,
          f"refs={len(external_refs(otxt))}")

    # ═══════════════ E5 联动: 报告不崩 ═══════════════
    rep = s3.build_phylogeo_report(jp, str(cnt_csv), str(tmp / "repdir"))
    rtxt = Path(rep).read_text(encoding="utf-8")
    check("E5 报告生成且 None 坐标渲染为 —", "—" in rtxt and "AtlantisDeep" in rtxt, "")
    check("E5 报告显式提示缺坐标", "解析不到经纬度" in rtxt and "未画点" in rtxt, "")
    check("E5 报告未把 None 格式化成 '?' 或崩", "None" not in rtxt and "nan" not in rtxt.lower(), "")

    # ═══════════════ E1 端到端: 真跑 build_html_report ═══════════════
    print("\n[E1-e2e] build_html_report 整份报告")
    # 模块声明核对: geography 段 (virspacetime) 声明 spatiotemporal_counts.csv 已存在;
    # 另外塞一张 PDF (必须 <embed>) 与一张超限 PDF (必须显式 ⚠)
    results = [{
        "virus": "DemoVirus", "phase": "test", "success": True,
        "virus_dir": str(vdir),
        "outputs": {}, "metrics": {},
        "stages_done": ["virspacetime", "spread3"],
        "stages_skipped": [], "stages_failed": [],
    }]
    outdir = tmp / "report_out"
    outdir.mkdir(parents=True, exist_ok=True)
    rpath = rb.build_html_report(results, str(outdir))
    full = Path(rpath).read_text(encoding="utf-8")
    check("E1-e2e 报告已生成", os.path.exists(rpath) and len(full) > 2000, rpath)
    check("E1-e2e PDF 以 <embed class=\"pdf-embed\"> 入报告",
          'class="pdf-embed"' in full, "")
    check("E1-e2e 超限 PDF 显式列在「有产物未内联」",
          "有产物未内联" in full and "huge_map.pdf" in full, "")
    check("E1-e2e 同名单跳过项以「未嵌入的候选」列出",
          "未嵌入的候选" in full and "spatial_map.pdf" in full, "")
    check("E1-e2e png 走 <img> 内联",
          'src="data:image/png' in full, "")
    decl_mod = [m for m in rb.MODULES if m["dir"] == "geography/geo_analysis"]
    print(f"  (geography 模块声明产物: {decl_mod[0]['summaries'] if decl_mod else 'N/A'})")

finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"\n{'=' * 60}\nPASS {len(OK)} / FAIL {len(FAIL)}")
for f in FAIL:
    print(f"  ✗ {f}")
sys.exit(1 if FAIL else 0)
