# -*- coding: utf-8 -*-
"""
VirusPlatform 功能树 —— 对称双向思维导图 + 底部结果图网格
用 matplotlib 绘制（避免 R 的长度限制和转义陷阱）
用法: python plot_functional_tree_sides.py <EN|CN> <FULL|SIMPLE>
输出: figures/functional_tree_sides_<lang>_<dens>.pdf/.png
"""
import sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
import numpy as np
from PIL import Image

lang = (sys.argv[1] if len(sys.argv) > 1 else 'CN').upper()
DENS = (sys.argv[2] if len(sys.argv) > 2 else 'FULL').upper()
if lang not in ('EN', 'CN'): lang = 'CN'
if DENS not in ('FULL', 'SIMPLE'): DENS = 'FULL'

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ASSETS = os.path.join(ROOT_DIR, 'figures', 'assets', 'grid')
OUTDIR = os.path.join(ROOT_DIR, 'figures')

# ---- 字体 --------------------------------------------------------------
if lang == 'CN':
    fp = fm.FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')
    fname_b = r'C:\Windows\Fonts\msyhbd.ttc'
    FB = fm.FontProperties(fname=fname_b)
else:
    fp = fm.FontProperties(fname=r'C:\Windows\Fonts\arial.ttf')
    fname_b = r'C:\Windows\Fonts\arialbd.ttf'
    FB = fm.FontProperties(fname=fname_b)

# ---- 树结构（域, 模块, [末端...]）---------------------------------------
SPEC = {
'CN': {
'FULL': [
 ('数据资源与获取', '数据资源', [
   '公共数据检索（SRA + GSA）', '公共数据批量下载（aria2c + sracha）',
   '基因组/宿主基因组下载', 'Open-Virome 公共病毒组',
   '病毒样本预筛', '本地数据库构建']),
 ('病毒检测与鉴定', '样本处理', [
   '样品创建与批量导入', '一键前处理', '测序质控（fastp）',
   '格式转换（fq2fa）', '宿主去除（kunpeng）']),
 ('病毒检测与鉴定', '定量与共识', [
   '已知病毒识别与定量（salmon）', '覆盖深度与 BAM',
   '共识序列', '变异分析']),
 ('病毒检测与鉴定', '识别与组装', [
   '一键分析（②→④）', '宏病毒组组装（SPAdes）',
   '候选序列验证（BLASTx / BFVD）', '组装结果再鉴定',
   '宿主预测（ICTV 级联）']),
 ('注释与比较基因组', '基因组注释', [
   'ORF 预测（pyrodigal）', '功能注释（DIAMOND / Pfam / BFVD）',
   '同源比对（BLASTN / BLASTX）', '基因组图谱（gbdraw / DFV）']),
 ('注释与比较基因组', '比较基因组', [
   '参考获取与 CDS/PEP 提取', '多序列比对（MAFFT + trimAl）',
   '进化树构建（RAxML-NG / IQ-TREE）', 'SDT 同一性分析']),
 ('进化与动力学分析', '进化动力学', [
   '元数据治理与比对 QC', '时间信号与定年', 'ML 分子钟（TreeTime）',
   '天际线与 BSP', 'RTT / DRT 回归', '系统地理 + 迁移动画',
   '重组分析（RDP5）', 'BEAST 交接（Markov jump）',
   '一键定年→地理→溯源']),
 ('检测与防治设计', '检测与防治试剂设计', [
   '引物设计（primer3）', 'CRISPR 检测试剂设计',
   'dsRNA 设计全链路', 'miRNA 靶标预测（5 引擎）']),
 ('溯源与提交', '溯源与提交', [
   'LOGAN 溯源（2340 万 SRA）', 'NCBI 提交准备',
   'AI 助手（结果解读 · 检索执行）']),
],
'SIMPLE': [
 ('数据资源与获取', '', [
   '公共数据检索（SRA + GSA）', '参考与宿主基因组下载',
   '病毒样本预筛', '数据库构建']),
 ('病毒检测与鉴定', '', [
   '宿主去除', '已知病毒识别与定量', '宏基因组组装', '宿主预测']),
 ('注释与比较基因组', '', [
   'ORF 预测与功能注释', '基因组图谱',
   '进化树构建', 'SDT 同一性分析']),
 ('进化与动力学分析', '', [
   '分子钟定年', '系统地理分析', '重组分析（RDP5）', '天际线与 BSP']),
 ('检测与防治设计', '', [
   '引物设计', 'CRISPR 检测试剂设计', 'dsRNA 设计', 'miRNA 靶标预测']),
 ('溯源与提交', '', [
   'LOGAN 溯源', 'NCBI 提交准备', 'AI 助手']),
]},
'EN': {
'FULL': [
 ('Data acquisition', 'Data resources', [
   'Public data search (SRA + GSA)', 'Bulk download (aria2c + sracha)',
   'Genome / host genome download', 'Open-Virome public viromes',
   'Sample virome prescreening', 'Local database construction']),
 ('Virus Detection & Identification', 'Sample processing', [
   'Sample creation & batch import', 'One-click preprocessing',
   'Read QC (fastp)', 'Format conversion (fq2fa)',
   'Host read removal (kunpeng)']),
 ('Virus Detection & Identification', 'Quantification & consensus', [
   'Known virus ID & quantification', 'Coverage depth & BAM',
   'Consensus sequence calling', 'Variant analysis']),
 ('Virus Detection & Identification', 'Identification & assembly', [
   'One-click pipeline (stage 2-4)', 'Metaviral assembly (SPAdes)',
   'Candidate validation', 'Assembly re-identification',
   'Host prediction (ICTV cascade)']),
 ('Annotation & genomics', 'Genome annotation', [
   'ORF prediction (pyrodigal)', 'Functional annotation (DIAMOND / Pfam)',
   'Homology search (BLASTN / BLASTX)', 'Genome maps (gbdraw / DFV)']),
 ('Annotation & genomics', 'Comparative genomics', [
   'Reference acquisition & CDS/PEP', 'Multiple alignment (MAFFT)',
   'Tree inference (RAxML-NG / IQ-TREE)', 'SDT pairwise identity']),
 ('Evolution & Phylodynamics', 'Phylodynamics', [
   'Metadata governance & alignment QC', 'Temporal signal & dating',
   'ML molecular clock (TreeTime)', 'Skyline & BSP plot',
   'RTT / DRT regression', 'Phylogeography + migration',
   'Recombination (RDP5)', 'BEAST interface (Markov jump)',
   'One-click dating → geography → tracing']),
 ('Assay design', 'Assay & control design', [
   'Primer design (primer3)', 'CRISPR detection assay',
   'dsRNA design pipeline', 'miRNA target prediction']),
 ('Tracing & submission', 'Tracing & submission', [
   'LOGAN tracing (23.4M SRA runs)', 'NCBI submission',
   'AI assistant (interpret / retrieval)']),
],
'SIMPLE': [
 ('Data acquisition', '', [
   'Public data search (SRA + GSA)', 'Reference & host genome dl',
   'Sample virome prescreening', 'Database construction']),
 ('Virus Detection & ID', '', [
   'Host read removal', 'Known virus ID & quant',
   'Metaviral assembly', 'Host prediction']),
 ('Annotation & genomics', '', [
   'ORF prediction & annot', 'Genome maps',
   'Tree inference', 'SDT identity analysis']),
 ('Evolution & Phylodyn', '', [
   'Molecular clock dating', 'Phylogeography',
   'Recombination (RDP5)', 'Skyline & BSP']),
 ('Assay design', '', [
   'Primer design', 'CRISPR assay design',
   'dsRNA design', 'miRNA target prediction']),
 ('Tracing & submission', '', [
   'LOGAN tracing', 'NCBI submission', 'AI assistant']),
]}}

TREE = SPEC[lang][DENS]

# ---- 布局参数 -----------------------------------------------------------
ROW_H = 0.42           # 每末端行距（英寸）
DOM_GAP = 0.85         # 功能域之间间隙（英寸）
MOD_GAP = 0.45         # 模块之间间隙（英寸）
ROOT_H = 0.9           # 根区域高
TRUNK_X = 0.0          # 主干 x（中心）
DOM_BUS_X = 1.5        # 功能域总线偏移
MOD_BUS_X = 2.8        # 模块总线偏移
TIP_X = 3.8            # 末端 x
LABEL_PAD = 0.12       # 末端标签偏移

# ---- 行分配（平铺所有末端）---------------------------------------------
tips_flat = []  # (domain, module, tip, global_row)
row = 0
for dm, mod, tips in TREE:
    for t in tips:
        tips_flat.append((dm, mod, t, row))
        row += 1
ntips_total = row
total_height = row * ROW_H + ROOT_H + MOD_GAP * 2

# 域中心
dom_centers = {}
for dm, mod, tips in TREE:
    ds = [r for d2, m2, t2, r in tips_flat if d2 == dm]
    if ds: dom_centers[dm] = (min(ds) + max(ds)) / 2
# 模块中心
mod_centers = {}
for dm, mod, tips in TREE:
    ms = [r for d2, m2, t2, r in tips_flat if d2 == dm and m2 == mod]
    if ms: mod_centers[(dm, mod)] = (min(ms) + max(ms)) / 2

# ---- 域配色 -------------------------------------------------------------
DOMAIN_COLORS = {
 '数据资源与获取': ('#0072B2', '#EAF2F9'),
 '病毒检测与鉴定': ('#009E73', '#E7F4EE'),
 '注释与比较基因组': ('#E69F00', '#FBF2DF'),
 '进化与动力学分析': ('#CC79A7', '#F8EDF4'),
 '检测与防治设计': ('#D55E00', '#FAECE3'),
 '溯源与提交': ('#4D4D4D', '#EEEEEE'),
 'Data acquisition': ('#0072B2', '#EAF2F9'),
 'Virus Detection & Identification': ('#009E73', '#E7F4EE'),
 'Annotation & genomics': ('#E69F00', '#FBF2DF'),
 'Evolution & Phylodynamics': ('#CC79A7', '#F8EDF4'),
 'Assay design': ('#D55E00', '#FAECE3'),
 'Tracing & submission': ('#4D4D4D', '#EEEEEE'),
}
DOMAIN_ORDER = [d for d, m, ts in TREE]

# ---- 图像网格数据 -------------------------------------------------------
ASSET_DIR = os.path.join(ROOT_DIR, 'figures', 'assets', 'grid')
IMAGES = {
'CN': [
 ('ui/homepage.png', 1, '平台首页 · 总览'),
 ('meta_temporal.png', 1, '检索元数据 · 时间分布'),
 ('meta_regions.png', 1, '检索元数据 · 地理分布'),
 ('meta_organizations.png', 1, '检索元数据 · 机构分布'),
 ('meta_tissues.png', 1, '检索元数据 · 组织分布'),
 ('meta_stages.png', 1, '检索元数据 · 发育阶段'),
 ('kvsuite_depth.png', 2, '覆盖深度图（kvsuite）'),
 ('viromock_quant.png', 2, '检出与定量基准（VIROMOCK）'),
 ('variant_map.png', 2, '变异位点图（snpEff）'),
 ('variant_effect.png', 2, '变异效应图（snpEff）'),
 ('variant_overview.png', 2, '变异总览'),
 ('stage_funnel.png', 2, '分析流程漏斗'),
 ('gcva_genome_circular.png', 3, 'GCVA 基因组圈图'),
 ('gcva_sdt_heatmap.png', 3, 'GCVA SDT 同一性热图'),
 ('gcva_identity_composite.png', 3, 'GCVA 同一性矩阵'),
 ('gcva_nt_heatmap.png', 3, 'GCVA NT 同一性热图'),
 ('gcva_tree.png', 3, 'GCVA 进化树'),
 ('gcva_rdp_breakpoints.png', 4, 'GCVA 重组断点分布'),
],
'EN': [
 ('ui/homepage.png', 1, 'Platform home'),
 ('meta_temporal.png', 1, 'Metadata: temporal'),
 ('meta_regions.png', 1, 'Metadata: geographic'),
 ('meta_organizations.png', 1, 'Metadata: organizations'),
 ('meta_tissues.png', 1, 'Metadata: tissues'),
 ('meta_stages.png', 1, 'Metadata: stages'),
 ('kvsuite_depth.png', 2, 'Depth plot (kvsuite)'),
 ('variant_map.png', 2, 'Variant landscape'),
 ('variant_effect.png', 2, 'Variant effects'),
 ('variant_overview.png', 2, 'Variant overview'),
 ('gcva_genome_circular.png', 3, 'GCVA genome map'),
 ('gcva_sdt_heatmap.png', 3, 'GCVA SDT heatmap'),
 ('gcva_identity_composite.png', 3, 'GCVA identity matrix'),
 ('gcva_nt_heatmap.png', 3, 'GCVA NT heatmap'),
 ('gcva_tree.png', 3, 'GCVA phylogeny'),
 ('gcva_rdp_breakpoints.png', 4, 'GCVA recombination'),
 ('gcva_sdt_distribution.png', 4, 'GCVA identity distribution'),
]}

# ---- 绘图 ---------------------------------------------------------------
W_IN = 14.5
DPI = 300
total_h = ROOT_H + ntips_total * ROW_H + DOM_GAP * len(TREE) + MOD_GAP * len(mod_centers) + 2.5
fig, ax = plt.subplots(1, 1, figsize=(W_IN, total_h))
ax.set_xlim(-0.5, W_IN - 0.5)
ax.set_ylim(total_h, -0.5)
ax.set_aspect('equal')
ax.axis('off')

# y 坐标：向下增大（y = 行号 × ROW_H）

# ---- 1) 根标签（顶部居中）----------------------------------------------
root_y = 0.5
ax.text(W_IN / 2, root_y,
        '植物病毒分析平台' if lang == 'CN' else 'VirusPlatform',
        ha='center', va='center', fontsize=14, fontweight='bold',
        fontfamily=fp if lang == 'CN' else 'Arial', color='black')

# ---- 2) 功能域色带 + 域标签 + 模块标签 ----------------------------------
trunk_x = W_IN / 2 - 2.5
dom_bus_x = trunk_x + 1.2
mod_bus_x = trunk_x + 2.6
tip_x = trunk_x + 3.8

# 功能域中心 y
def dom_ymid(dm):
    ds = [r for d2, m2, t2, r in tips_flat if d2 == dm]
    return (min(ds) + max(ds)) / 2 if ds else 0

def mod_ymid(dm, mod):
    ms = [r for d2, m2, t2, r in tips_flat if d2 == dm and m2 == mod]
    return (min(ms) + max(ms)) / 2 if ms else 0

# ---- 3) 树边 ------------------------------------------------------------
for g_idx, (dm, mod, tips) in enumerate(TREE):
    c = DOMAIN_COLORS[dm][0]
    yd = dom_ymid(dm)
    # 根→域分支
    ax.plot([trunk_x, dom_bus_x], [ROOT_TRUNK_END, yd],
            color=c, linewidth=1.2, solid_capstyle='round')
    # 域→模块
    for mi, (m2, ts) in enumerate(zip([x[1] for x in TIP_GROUPS[dm]],
                                       [x[1] for x in TIP_GROUPS[dm]])):
        pass  # 简化：不画模块连线

# ---- 简化绘制（直接从 tips_flat + TREE 布局）----------------------------
# 域总线（垂直线）
for dm, mod, tips in TREE:
    c = DOMAIN_COLORS[dm][0]
    yd = dom_ymid(dm)
    ds = sorted([r for d2, m2, t2, r in tips_flat if d2 == dm])
    if not ds: continue
    y_top = ds[0] - 0.5; y_bot = ds[-1] + 0.5
    ax.plot([dom_bus_x, dom_bus_x], [y_top, y_bot], color=c, linewidth=1.5)
    # 根→域分支（水平线）
    ax.plot([trunk_x, dom_bus_x], [ROOT_TRUNK_END, yd],
            color=c, linewidth=1.2, solid_capstyle='round')
    # 模块分支
    for d2, m2, t2, r in tips_flat:
        if d2 != dm: continue
        ax.plot([dom_bus_x, mod_bus_x], [r, r], color=c, linewidth=0.6)
        ax.plot([mod_bus_x, mod_bus_x], [y0_mod, y1_mod], color=c, linewidth=0.6)
        ax.plot([mod_bus_x, tip_x], [r, r], color='grey35', linewidth=0.35)

# ---- 4) 域标签（色带上方）----------------------------------------------
for dm, mod, tips in TREE:
    c = DOMAIN_COLORS[dm][0]
    yd = dom_ymid(dm)
    ax.text(dom_bus_x - 0.3, yd, dm, ha='right', va='center',
            fontsize=8, fontweight='bold', fontfamily=fp, color=c)

# ---- 5) 模块标签（色带内侧）--------------------------------------------
for (dm, mod), tips in MOD_GROUPS.items():
    ym = mod_ymid(dm, mod)
    ax.text(mod_bus_x + 0.15, ym, mod, ha='left', va='center',
            fontsize=6.5, fontweight='bold', fontfamily=fp, color='#2B2B2B')

# ---- 6) 末端标签 --------------------------------------------------------
for dm, mod, tips in TREE:
    for t in tips:
        r = row_of[t]
        y = r * ROW_H + ROOT_H
        ax.text(tip_x + 0.3, y, t, ha='left', va='center',
                fontsize=6, fontfamily=fp, color='#1A1A1A')

# ---- 7) 底部结果图网格 ---------------------------------------------------
grid_top = total_height + 0.5
grid_imgs = IMAGES[lang]
ncol_g = 6
cell_w_g = W_IN / ncol_g
cell_h_g = 1.55

for k, (fname, dom_id, cap) in enumerate(grid_imgs):
    kr = k
    col = kr % ncol_g; rr = kr // ncol_g
    x0 = col * cell_w_g + 0.1
    y0 = grid_top + rr * cell_h_g
    fpath = os.path.join(ASSET_DIR, fname)
    if not os.path.exists(fpath):
        continue
    img = Image.open(fpath)
    iw, ih = img.size
    # 等比缩放
    draw_w = cell_w_g - 0.2
    draw_h = draw_w * ih / iw
    if draw_h > cell_h_g - 0.55:
        draw_h = cell_h_g - 0.55
        draw_w = draw_h * iw / ih
    img_reshaped = np.array(img.resize((int(draw_w * 100), int(draw_h * 100))))
    ax_img = fig.add_axes((x0 / W_IN, 1 - (y0 + draw_h + 0.45) / total_h,
                           draw_w / W_IN, draw_h / total_h))
    ax_img.imshow(img_reshaped)
    ax_img.axis('off')
    # 灰边框
    ax_img.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False,
                     edgecolor='grey55', linewidth=0.5))
    # 图注
    fig.text((x0 + 0.1) / W_IN, 1 - (y0 + draw_h + 0.55) / total_h,
             cap, ha='left', va='top', fontsize=5.5, style='italic',
             color='#444444', fontfamily=fp)

fig.savefig(os.path.join(OUTDIR, f'functional_tree_sides_{lang}_{DENS.lower()}.pdf'),
            format='pdf', bbox_inches='tight', pad_inches=0.3, facecolor='white')
fig.savefig(os.path.join(OUTDIR, f'functional_tree_sides_{lang}_{DENS.lower()}.png'),
            dpi=300, bbox_inches='tight', pad_inches=0.3, facecolor='white')
plt.close(fig)
print(f'saved: figures/functional_tree_sides_{lang}_{DENS.lower()}')
