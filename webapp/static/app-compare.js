// ---------------- 比较基因组：⑤ 参考序列下载 / ⑥ GenBank 集合与共线性 ----------------
// 2026-09-13 自 app.js 拆出（原 3545-4349 行）。普通 script 全局共享作用域：
// 本文件在 app.js 之后加载，可直接使用 $/esc/t/toast/jfetch 等核心工具；
// 页面内联脚本调用本文件函数同样按全局解析，无需 export。
// （本块只服务「比较基因组」页；后续拆分沿用同一模式。）
// ---------------- 比较基因组：⑤ 参考序列下载 / ⑥ GenBank 集合与共线性 ----------------
/* 原 tools.html 内联脚本迁移至此；输入取值用 _v（避免与页面内联 val 重名）。 */
function _v(id) { return ($(id)?.value || '').trim(); }

function _watchTask(tid, onDone, btn, busyText) {
  const unlock = lockBtn(btn, busyText || '⏳ 运行中…');
  taskLogOpen.add(tid);
  startPolling();
  const timer = setInterval(async () => {
    try {
      const snap = await (await fetch('/api/task/' + tid)).json();
      if (['done', 'failed', 'cancelled'].includes(snap.status)) {
        clearInterval(timer);
        unlock();
        if (onDone) onDone(snap);
      }
    } catch (e) { /* 轮询继续 */ }
  }, 5000);
}

async function loadNcbiCollections() {
  const box = $('ncbiCollections');
  if (!box) return;
  try {
    const cols = await (await fetch('/api/ncbi/collections')).json();
    box.innerHTML = cols.length
      ? '<table class="table" style="width:100%;border-collapse:collapse"><tr><th>集合</th><th>序列数</th><th>库</th><th>检索式</th><th>日期</th></tr>' +
        cols.map(c => `<tr><td><b>${esc(c.name)}</b></td><td>${c.n_seqs}</td><td>${esc(c.db)}</td>` +
          `<td title="${esc(c.query)}">${esc((c.query || '').slice(0, 50))}</td><td>${esc(c.date)}</td></tr>`).join('') + '</table>' +
        '<p class="hint">分析管道 → 「NCBI 参考集合」填集合名即可把这些参考追加进进化树比对。</p>'
      : '<p class="hint">（暂无集合，先搜索并下载）</p>';
  } catch (e) { box.innerHTML = '<p class="hint">加载失败</p>'; }
}

async function gbDownload(btn) {
  const name = collName(), term = _v('s_term'), accs = _v('s_acc');
  if ((!term && !accs) || (term && accs)) {
    alert(t('ui.cmp.needQueryOrAcc')); return;
  }
  try {
    const r = await fetch('/api/gb/download', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, term, accessions: accs,
                             max: +_v('s_max') || 50 })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 下载中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

async function gbImport(btn) {
  const files = _v('s_files');
  if (!files) { alert(t('ui.cmp.needGbPath')); return; }
  const name = collName();
  try {
    const r = await fetch('/api/gb/import', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name,
                             files: files.split(',').map(s => s.trim()).filter(Boolean) })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 导入中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

/* ================= CDS/PEP 提取（PhyloSuite 布局） ================= */
async function gbExtract(name, btn) {
  try {
    const r = await fetch('/api/gb/extract', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => {
      loadGbCollections();
    }, btn, '⏳ 提取中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

/* 显式巡检：重解析全部 .gb，重建清单与警告（列表页只读已生成的清单） */
async function gbInspect(name, btn) {
  try {
    const r = await fetch('/api/gb/inspect', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 巡检中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

/* 集合建树：三种来源——
   genome  全基因组树（集合 extract/genome.fa → MAFFT → 建树，molecule=genome）
   selected 挑选的序列集（extract/selected/CDS.fa|PEP.fa → MAFFT → 建树）
   aligned  已比对 FASTA（跳过 MAFFT 直接建树） */
async function gbBuildTree(name, btn) {
  const treeTool = ($('s_tree_tool')?.value) || 'fasttree';
  const src = ($('tbSource') && $('tbSource').value) || 'genome';
  try {
    if (src === 'genome') {
      const r = await fetch('/api/gb/phylo', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name, tree_tool: treeTool })});
      if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
      const d = await r.json();
      _watchTask(d.task, loadGbCollections, btn, '⏳ 建树中…');
      return;
    }
    const mol = ($('tbSelectedMol')?.value) || 'CDS';
    const seqs = src === 'aligned'
      ? (_v('tbAlignedPath') || '')
      : `run/gb_collections/${name}/extract/selected/${mol}.fa`;
    if (!seqs) { alert(t('ui.cmp.needAlnPath')); return; }
    const method = treeTool === 'nj' ? 'nj' : 'fasttree';
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tool: 'quicktree',
        params: { seqs, method,
                  aligned: src === 'aligned' ? '1' : '0' } })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    if (typeof taskLogOpen !== 'undefined') taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, null, btn, '⏳ 建树中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

/* 建树来源切换：按来源显示/隐藏对应控件 */
function tbSourceChanged() {
  const src = ($('tbSource') && $('tbSource').value) || 'genome';
  const selItem = $('tbSelectedItem'), alItem = $('tbAlignedItem');
  const useAl = $('tbUseAlPath'), coll = $('tbColl');
  if (selItem) selItem.style.display = src === 'selected' ? '' : 'none';
  if (alItem) alItem.style.display = src === 'aligned' ? '' : 'none';
  if (useAl) useAl.style.display = src === 'aligned' ? '' : 'none';
  /* 全基因组树必须选集合；另两种来源集合仍用于定位 selected/ 目录 */
  if (coll) coll.style.opacity = (src === 'aligned') ? '0.5' : '';
}

/* 把比对卡里当前打开的比对路径带入建树卡 */
function tbUseAlPath() {
  const p = (alData && alData.path) || _v('alPath');
  if (!p) { alert(t('ui.cmp.noAlnOpen')); return; }
  const inp = $('tbAlignedPath');
  if (inp) inp.value = p;
}

/* 建树卡 → 比对卡：带着集合/分组跳到「序列比对」的比对查看器看建树比对 */
function tbSendToAlign() {
  const coll = _v('tbColl');
  if (!coll) { alert(t('ui.cmp.needGbSetOrAln')); return; }
  const sample = 'gb:' + coll;
  location.hash = '#t-align';
  const sel = $('alSample');
  if (!sel) return;
  const opt = Array.from(sel.options).find(o => o.value === sample);
  if (!opt) { alert(t('ui.cmp.setNoAln')); return; }
  sel.value = sample;
  alPickGroups();
}

/* 比对卡来源切换：全长序列 / 挑选序列集 / 手填 */
function alSourceChanged() {
  const src = ($('alSource') && $('alSource').value) || 'manual';
  const item = $('alSelItem');
  if (item) item.style.display = src === 'selected' ? '' : 'none';
}

/* 按来源填入比对输入框 */
async function alFillFromSource() {
  const src = ($('alSource') && $('alSource').value) || 'manual';
  if (src === 'manual') { alert(t('ui.cmp.manualNeedInput')); return; }
  const coll = _v('alColl');
  if (!coll) { alert(t('ui.cmp.needGbSet')); return; }
  const inp = $('al_fa');
  if (!inp) return;
  if (src === 'genome') {
    inp.value = `run/gb_collections/${coll}/extract/genome.fa`;
  } else {
    const mol = ($('alSelMol')?.value) || 'CDS';
    inp.value = `run/gb_collections/${coll}/extract/selected/${mol}.fa`;
  }
}

/* 比对卡集合下拉（来源为集合时用） */
async function loadAlColls() {
  const sel = $('alColl');
  if (!sel) return;
  try {
    const cols = await (await fetch('/api/gb/collections')).json();
    sel.innerHTML = cols.length
      ? cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}（${c.n_records} 条记录）</option>`).join('')
      : '<option value="">（暂无集合——先到「参考序列获取」下载）</option>';
  } catch (e) { sel.innerHTML = '<option value="">（无法连接）</option>'; }
}

async function loadGbCollections() {
  /* 集合列表双渲染：
     #gbCollectionsMgmt（GenBank 集合管理卡）= 🧬 提取 / 🔍 巡检；
     #gbWarnings = 巡检警告汇总。 */
  const mgmtBox = $('gbCollectionsMgmt');
  if (!mgmtBox) return;
  try {
    const cols = await (await fetch('/api/gb/collections')).json();
    const gbRow = c => {
      const cds = (c.records || []).reduce((s, r) => s + (+r.cds || 0), 0);
      const src = c.source === 'query' ? esc((c.term || '').slice(0, 36))
        : (c.source === 'accessions' ? esc(c.accessions) + ' 个 accession' : '本机导入');
      const results = [
        c.has_phylo ? `<a href="/results#msa" title="结果中心查看比对/树/SDT">🌳 MSA·树</a>` : '',
      ].filter(Boolean).join(' ') || '—';
      const warnN = (c.warnings || []).length;
      const warnTip = warnN ? esc(c.warnings.join(' | ')) : '';
      return `<tr><td><b>${esc(c.name)}</b>${warnN ? ` <span class="hint" style="color:#b45309" title="${warnTip}">⚠${warnN}</span>` : ''}</td>` +
        `<td>${c.n_records}</td><td>${cds}</td>` +
        `<td>${src}</td><td>${esc(c.date)}</td><td>${results}</td>` +
        `<td style="white-space:nowrap">[[ACTIONS]]</td></tr>`;
    };
    const mkTable = rows => rows.length
      ? '<table class="table" style="width:100%;border-collapse:collapse">' +
        '<tr><th>集合</th><th>记录数</th><th>CDS 数</th><th>来源</th><th>日期</th><th>结果</th><th></th></tr>' +
        rows.join('') + '</table>'
      : '<p class="hint">（暂无集合）</p>';
    if (mgmtBox) {
      mgmtBox.innerHTML = mkTable(cols.map(c => gbRow(c).replace('[[ACTIONS]]',
        `<button class="btn small primary" data-name="${esc(c.name)}" onclick="gbExtract(this.dataset.name, this)" title="提取 genome / CDS / PEP（分类分目录 + 按基因拆分）">🧬 提取</button> ` +
        `<button class="btn small" data-name="${esc(c.name)}" onclick="gbInspect(this.dataset.name, this)" title="重解析全部 .gb，重建清单与警告">🔍 巡检</button>`)));
    }
    const wbox = $('gbWarnings');
    if (wbox) {
      const warned = cols.filter(c => (c.warnings || []).length);
      wbox.innerHTML = warned.map(c =>
        `<details><summary class="hint" style="color:#b45309;cursor:pointer">⚠ ${esc(c.name)}：${c.warnings.length} 条巡检警告</summary>` +
        '<ul class="hint" style="margin:4px 0 8px">' +
        c.warnings.map(w => `<li>${esc(w)}</li>`).join('') + '</ul></details>').join('');
    }
  } catch (e) {
    if (mgmtBox) mgmtBox.innerHTML = '<p class="hint">加载失败</p>';
  }
}

loadNcbiCollections();
loadGbCollections();
ictvCascadeRefetch();
loadTbColls();
loadAlColls();
tbSourceChanged();
alSourceChanged();

/* ================= 比较基因组独立工作区：⑦ MSA / ⑧ 进化树 / ⑨ SDT ================= */
/* 三张卡片各自的输入/参数/输出都独立于结果中心；渲染器为共享实现。 */

/* ---- 样品/集合列表（MSA / 树 / SDT 卡与 results 页共用，缓存一次） ---- */
let sampleListCache = null;
async function fetchSampleList() {
  if (!sampleListCache) {
    try {
      sampleListCache = await (await fetch('/api/msa/samples')).json();
    } catch (e) { sampleListCache = []; }
  }
  return sampleListCache;
}
function sampleLabel(x) { return x.label || x.sample; }

/* ---- 共享渲染：SNP-only 热图 ---- */
const SNP_COLORS = { A: '#4daf4a', C: '#377eb8', G: '#ffb200',
                     T: '#e41a1c', U: '#e41a1c', '-': '#e8edf2' };

function snpTableHtml(d, cols, pageNo, showCons, showDiv) {
  /* d = /api/msa/data 或 /api/tool/msa_data 的返回。返回 {html, info, meta}。 */
  const nChunks = Math.max(1, Math.ceil(d.n_snp / cols));
  const start = pageNo * cols, end = Math.min(d.n_snp, start + cols);
  const pos = d.positions.slice(start, end);
  const tips = pos.map((p, ci) => {
    const counts = {};
    for (const r of d.rows) counts[r[ci]] = (counts[r[ci]] || 0) + 1;
    return `${t('c.snpPos')} ${p} · ` + Object.entries(counts)
      .sort((a, b) => b[1] - a[1]).map(([c, n]) => `${c}:${n}`).join(' ');
  });
  let html = '<table style="border-collapse:collapse;font-family:Consolas,monospace;font-size:11px">';
  if (showDiv) {
    html += '<tr><td style="position:sticky;left:0;background:#fff;font-size:9px;color:#94a3b8;padding-right:6px;text-align:right">变异度</td>';
    for (let ci = start; ci < end; ci++) {
      const h = Math.round((d.diversity[ci] || 0) * 26);
      html += `<td style="width:16px;height:28px;vertical-align:bottom;padding:0">` +
        `<div style="height:${h}px;background:#dc2626;opacity:.75" title="${esc(tips[ci - start])}"></div></td>`;
    }
    html += '</tr>';
  }
  if (showCons) {
    html += '<tr><td style="position:sticky;left:0;background:#fff;font-size:9px;color:#94a3b8;padding-right:6px;text-align:right">' + t('c.consensus') + '</td>';
    for (let ci = start; ci < end; ci++) {
      const c = d.consensus[ci];
      html += `<td style="width:16px;text-align:center;background:#f1f5f9;font-weight:bold" ` +
        `title="${esc(tips[ci - start])}">${esc(c)}</td>`;
    }
    html += '</tr>';
  }
  for (let ri = 0; ri < d.names.length; ri++) {
    const nm = d.names[ri];
    const label = nm.length > 26 ? nm.slice(0, 25) + '…' : nm;
    html += `<tr><td style="position:sticky;left:0;background:#fff;font-size:10px;max-width:180px;` +
      `overflow:hidden;text-overflow:ellipsis;white-space:nowrap;padding-right:6px" title="${esc(nm)}">${esc(label)}</td>`;
    const row = d.rows[ri];
    for (let ci = start; ci < end; ci++) {
      const c = row[ci];
      const bg = SNP_COLORS[c] || '#984ea3';
      const fg = c === '-' ? '#cbd5e1' : '#fff';
      html += `<td style="width:16px;text-align:center;background:${bg};color:${fg}" ` +
        `title="${esc(nm)} · ${esc(tips[ci - start])}">${c === '-' ? '·' : esc(c)}</td>`;
    }
    html += '</tr>';
  }
  html += '</table>';
  return { html,
    info: d.n_snp ? ' ' + t('c.pageOf').replace('{p}', pageNo + 1).replace('{n}', nChunks).replace('{a}', pos[0]).replace('{b}', pos[pos.length - 1]) + ' ' : '',
    meta: `序列 ${d.n_seq} 条${d.seqs_truncated ? '（超出上限仅显示前 80）' : ''}` +
      ` · ${t('c.vsFull').replace('{l}', d.aln_len).replace('{n}', d.n_snp)}` };
}

/* ---- 共享渲染：identity 热图（plotly）与矩阵表 ---- */
function identityHeatmap(divId, names, matrix, title, zminManual) {
  const n = names.length;
  const lowest = Math.min(...matrix.flat());
  const zmin = Math.max(0, zminManual != null && zminManual !== '' ?
                        +zminManual : lowest - 5);
  const box = $(divId);
  Plotly.newPlot(divId, [{
    type: 'heatmap', z: matrix, x: names, y: names,
    zmin, zmax: 100,
    colorscale: [[0, '#b2182b'], [0.5, '#f7f7f7'], [1, '#1a9850']],
    colorbar: { title: { text: 'identity %' }, thickness: 12 },
    hovertemplate: '%{y} vs %{x}<br>%{z:.1f}%<extra></extra>',
  }], {
    title, width: Math.min(Math.max(box.clientWidth, 560), 1100),
    height: Math.min(Math.max(n * 24 + 160, 480), 900),
    margin: { l: 150, b: 120, t: 60, r: 20 },
    xaxis: { tickangle: -45, automargin: true },
    yaxis: { automargin: true, autorange: 'reversed' },
  }, { responsive: false });
  return zmin;
}

function identityMatrixTable(names, matrix) {
  const head = '<tr><th></th>' + names.map(n => `<th title="${esc(n)}">${esc(n.slice(0, 14))}</th>`).join('') + '</tr>';
  return '<tbody>' + head + names.map((n, i) =>
    `<tr><th style="text-align:left;white-space:nowrap" title="${esc(n)}">${esc(n.slice(0, 22))}</th>` +
    matrix[i].map(v => `<td>${v == null ? '' : (+v).toFixed(1)}</td>`).join('') + '</tr>').join('') + '</tbody>';
}

/* ---- 共享渲染：Archaeopteryx.js 树查看器 ----
   面板自带：矩形/环形/无根布局、系统发育图/聚类图切换、支持值与支持圆点、
   梯化排序、中点重根、字号/节点/枝宽滑块、双搜索框、子树上钻、
   Download 菜单（PNG/SVG/PDF/phyloXML/Newick/Nexus）。 */

let _tvViewer = null;

function _treeInternalLabelsAllNumeric(newick) {
  /* FastTree tree.nwk 与 RAxML-NG raxml.raxml.support（旧 IQ-TREE treefile/contree
     同款）把支持值写成内部节点名（如 )0.753: 或 )95.2:）；全部为数字时按支持值
     解析，否则保留为内部名。 */
  const labels = [];
  const re = /\)[ \t]*([^,():;\[\]\s]*)[ \t]*:/g;
  let m;
  while ((m = re.exec(newick)) !== null) {
    const lab = (m[1] || '').trim();
    if (lab) labels.push(lab);
  }
  return labels.length > 0 && labels.every(l => /^[0-9.]+$/.test(l));
}

function renderTreeTo(box, newick, opts) {
  /* opts: {layout, genomeOverlay}：'rectangular' | 'circular' | 'unrooted'。返回 tips 数。
     genomeOverlay（可选）=「基因组分带 / 每叶基因轨道 / 同源连线」叠加数据，
     只在矩形布局绘制；由 /api/tree/* 按「与树同名的 .overlay.json」给出。
     有 track 与叶名对不上时该叶不画轨道，不影响树本身。 */
  if (_tvViewer) { try { _tvViewer.destroy(); } catch (e) { /* 重建即可 */ } _tvViewer = null; }
  box.innerHTML = '';
  const asConf = _treeInternalLabelsAllNumeric(newick);
  /* 同一容器重新 launch 即为官方的换树方式；launch 会整体重建面板。 */
  const launchOpts = {
    layout: opts.layout || 'rectangular',
    nhConfidenceValuesAsInternalNames: asConf,
    nhConfidenceValuesInBrackets: !asConf,
    enableDownloads: true,
    enableAccessToDatabases: false,
    pngExportScale: 4,
  };
  if (opts.genomeOverlay) launchOpts.genomeOverlay = opts.genomeOverlay;
  _tvViewer = archaeopteryx.launchArchaeopteryx(box, 'tree.nwk', newick, launchOpts);
  /* 面板中支持值复选框默认不勾选；树带支持值时自动勾上。 */
  setTimeout(() => {
    const cb = box.querySelector('#conf_cb');
    if (asConf && cb && !cb.checked) cb.click();
  }, 60);
  let tips = 0;
  try {
    tips = forester.getAllExternalNodes(archaeopteryx.parseTree('tree.nwk', newick)).length;
  } catch (e) { tips = 0; }
  return tips;
}

/* ---- ② 比对查看器：从样品 / 集合直接选比对 ---- */
async function alPickGroups() {
  const s = _v('alSample'), sel = $('alGroup');
  if (!s) { sel.innerHTML = '<option value="">（先选样品）</option>'; return; }
  const list = await fetchSampleList();
  const item = list.find(x => x.sample === s);
  sel.innerHTML = (item?.groups || []).map(g =>
    `<option value="${esc(g.group)}">${esc(g.group)}${g.aln === 'trim' ? '（清剪后）' : ''}</option>`).join('')
    || '<option value="">（该样品没有比对）</option>';
  if (sel.value) alPickView();
}

async function alPickView() {
  const s = _v('alSample'), g = _v('alGroup');
  if (!s || !g) return;
  const box = $('alBox');
  box.innerHTML = '<p class="hint">加载中…</p>';
  try {
    const r = await fetch(`/api/msa/path?sample=${encodeURIComponent(s)}&group=${encodeURIComponent(g)}`);
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    const d = await r.json();
    $('alPath').value = d.path || '';
    await alignLoad();
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

async function alPickInit() {
  const sel = $('alSample');
  if (!sel) return;
  const list = await fetchSampleList();
  sel.innerHTML = list.length
    ? list.map(x => `<option value="${esc(x.sample)}">${esc(sampleLabel(x))}</option>`).join('')
    : '<option value="">（暂无比对结果）</option>';
  if (sel.value) alPickGroups();
}

/* ---- ⑧ 进化树查看器卡片 ---- */
let tvData = null, tvSource = '';

async function treeToolLoadFile() {
  const p = _v('tv_file');
  if (!p) { alert(t('ui.cmp.needNewick')); return; }
  try {
    const r = await fetch(`/api/tree/file?path=${encodeURIComponent(p)}`);
    if (!r.ok) { $('tvMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    const d = await r.json();
    tvData = { newick: d.newick, meta: { file: d.file, tool: '' },
               overlay: d.overlay || null };
    tvSource = 'file';
    treeToolRender();
  } catch (e) { $('tvMeta').textContent = '无法连接: ' + esc(e); }
}

async function treeToolLoadGroups() {
  const s = _v('tvSample'), sel = $('tvGroup');
  if (!s) { sel.innerHTML = '<option value="">（先选样品）</option>'; return; }
  const list = await fetchSampleList();
  const item = list.find(x => x.sample === s);
  const groups = (item?.groups || []).filter(g => g.trees && g.trees.length);
  sel.innerHTML = groups.map(g => `<option value="${esc(g.group)}">${esc(g.group)}</option>`).join('')
    || '<option value="">（该样品没有树文件）</option>';
  if (sel.value) treeToolLoadFiles();
}

function treeToolLoadFiles() {
  const sel = $('tvFileSel');
  fetchSampleList().then(list => {
    const item = list.find(x => x.sample === _v('tvSample'));
    const g = (item?.groups || []).find(x => x.group === _v('tvGroup'));
    sel.innerHTML = ((g && g.trees) || []).map(tr =>
      `<option value="${esc(tr.file)}">${esc(tr.label)}（${esc(tr.file)}）</option>`).join('')
      || '<option value="">（无树文件）</option>';
    if (sel.value) treeToolView();
  });
}

async function treeToolView() {
  const s = _v('tvSample'), g = _v('tvGroup'), f = _v('tvFileSel');
  if (!s || !g || !f) return;
  $('tvMeta').textContent = '加载中…';
  try {
    const r = await fetch(`/api/tree/data?sample=${encodeURIComponent(s)}&group=${encodeURIComponent(g)}&file=${encodeURIComponent(f)}`);
    if (!r.ok) { tvData = null; $('tvMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    const d = await r.json();
    tvData = { newick: d.newick, meta: d, overlay: d.overlay || null };
    tvSource = 'sample';
    treeToolRender();
  } catch (e) { $('tvMeta').textContent = '无法连接: ' + esc(e); }
}

function treeToolRender() {
  const box = $('tvBox');
  if (!tvData) { box.innerHTML = `<p class="hint">${t('c.noData2')}</p>`; $('tvMeta').textContent = ''; return; }
  let tips = 0;
  try {
    tips = renderTreeTo(box, tvData.newick,
                        { layout: 'rectangular', genomeOverlay: tvData.overlay || null });
  } catch (e) {
    box.innerHTML = '<p class="hint" style="color:#b91c1c">树解析失败: ' + esc(e) + '</p>';
    return;
  }
  const m = tvData.meta || {};
  $('tvMeta').textContent = [
    `${tips} 个序列`, m.tool, m.model ? '模型 ' + m.model : '',
    m.logl ? 'logL ' + m.logl : '', m.file, _tvOverlayNote(tvData.overlay),
    '查看器面板：布局切换 / 显示与支持值 / 缩放 / 搜索 / Download 导出',
  ].filter(Boolean).join(' · ');
}

function _tvOverlayNote(ov) {
  /* 叠加数据的一行摘要（无叠加返回 ''，不占位）。 */
  if (!ov || !(ov.tracks || []).length) return '';
  const nGene = (ov.tracks || []).reduce((s, t) => s + ((t.genes || []).length), 0);
  const nLink = (ov.links || []).length;
  return `基因组叠加：${ov.tracks.length} 条轨道 / ${nGene} 个基因`
    + (nLink ? ` / ${nLink} 条同源连线（按蛋白同一性着色）` : '')
    + '（仅矩形布局绘制）';
}

async function treeToolInit() {
  const sel = $('tvSample');
  if (!sel) return;
  const list = await fetchSampleList();
  const withTrees = list.map(x => ({
    sample: x.sample, label: sampleLabel(x),
    groups: (x.groups || []).filter(g => g.trees && g.trees.length)
  })).filter(x => x.groups.length);
  sel.innerHTML = withTrees.length
    ? withTrees.map(x => `<option value="${esc(x.sample)}">${esc(x.label)}</option>`).join('')
    : '<option value="">（暂无树文件，先建树）</option>';
  if (sel.value) treeToolLoadGroups();
}

/* ---- ⑧① 快速建树（FASTA → MAFFT → NJ / FastTree） ---- */
let qtRun = null;

async function quickTreeRun(btn) {
  const seqs = _v('qt_fa');
  if (!seqs) { alert(t('ui.common.needFasta2')); return; }
  const method = _v('qtMethod') || 'nj';
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tool: 'quicktree',
                             params: { seqs, method,
                                       max_n: +_v('qt_maxn') || 100 } })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    qtRun = d.run;
    $('tvMeta').textContent =
      `⏳ MAFFT 比对 + ${method === 'fasttree' ? 'FastTree' : 'NJ'} 建树中…`;
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => quickTreeLoad(), btn, '⏳ 建树中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

async function quickTreeLoad() {
  try {
    const r = await fetch(`/api/tool/quicktree_data?run=${encodeURIComponent(qtRun)}`);
    if (!r.ok) { $('tvMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    const d = await r.json();
    tvData = { newick: d.newick, meta: { file: d.file, tool: d.tool },
               overlay: d.overlay || null };
    tvSource = 'quicktree';
    treeToolRender();
    const dl = $('qtDl');
    dl.style.display = '';
    dl.href = `/tool_runs/${encodeURIComponent(qtRun)}/${encodeURIComponent(d.file)}`;
    dl.download = d.file;
  } catch (e) { $('tvMeta').textContent = '无法连接: ' + esc(e); }
}

/* ---- ⑨ SDT 卡片（精确引擎：逐对 MAFFT + Get_Similarity，替代 SDT exe） ---- */
let sdData = null, sdRun = null;

async function sdtRunMatrix(btn) {
  const seqs = _v('sd_fa');
  if (!seqs) { alert(t('ui.common.needFasta2')); return; }
  const mode = _v('sd_mode') || 'nt';          // nt | aa | ntaa
  const ntaa = mode === 'ntaa';
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(ntaa ? {
        /* NT+AA 同一性表模式（BioAider 口径） */
        tool: 'identity',
        threads: +_v('sd_threads') || null,
        params: { nt_seqs: seqs, aa_seqs: _v('sd_aa'),
                  max_n: +_v('sd_maxn') || 30,
                  aligned: !!$('sd_aligned').checked,
                  palette: _v('sd_palette') || 'sdt' } } : {
        /* NT / AA · SDT 精确矩阵模式（序列类型自动判别兜底） */
        tool: 'sdt',
        threads: +_v('sd_threads') || null,
        params: { seqs, max_n: +_v('sd_maxn') || 30,
                  seqtype: mode,
                  orient: !!$('sd_orient').checked,
                  aligned: !!$('sd_aligned').checked,
                  palette: _v('sd_palette') || 'sdt' } })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    sdRun = d.run;
    sdData = null;
    $('idtResult').style.display = 'none';
    $('sdExact').innerHTML = '';
    $('sdMeta').textContent = ntaa
      ? '⏳ NT+AA 同一性计算中（NT 与 AA 各做逐对比对）…'
      : (mode === 'aa'
        ? '⏳ AA 蛋白同一性逐对比对中…'
        : '⏳ SDT 逐对精确比对中（耗时与序列对数成正比，可用线程越多越快）…');
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => ntaa ? identityLoad(sdRun) : sdtExactLoad(),
               btn, '⏳ 矩阵分析中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

async function identityLoad(run) {
  try {
    const r = await fetch(`/tool_runs/${encodeURIComponent(run)}/identity.json`);
    if (!r.ok) { alert(t('ui.common.loadFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    const base = `/tool_runs/${encodeURIComponent(run)}`;
    $('idtHeat').innerHTML = d.aa
      ? `<a href="${base}/identity_composite.png" target="_blank"><img src="${base}/identity_composite.png"
           style="max-width:100%;max-height:720px;border:1px solid #dfe5ec;border-radius:8px"></a>`
      : '<p class="hint">有效 AA 序列不足 2 条，未出复合热图（仅 NT 矩阵）</p>';
    const rows = (d.pairs || []).map((p, i) =>
      `<tr><td>${i + 1}</td><td style="font-family:monospace" title="${esc(p.a)}">${esc(p.a)}</td>` +
      `<td style="font-family:monospace" title="${esc(p.b)}">${esc(p.b)}</td>` +
      `<td>${p.nt == null ? '—' : p.nt.toFixed(2)}</td>` +
      `<td>${p.aa == null ? '—' : p.aa.toFixed(2)}</td></tr>`).join('');
    $('idtTable').innerHTML =
      '<tbody><tr><th>#</th><th>序列 A</th><th>序列 B</th><th>NT Identity (%)</th><th>AA Identity (%)</th></tr>' +
      rows + '</tbody>';
    $('idtDl').innerHTML =
      `<a class="btn small" href="${base}/identity_table.csv?dl=1">⬇ 逐对同一性表 CSV</a>` +
      ` <a class="btn small" href="${base}/nt_matrix.csv?dl=1">⬇ NT 矩阵 CSV</a>` +
      (d.aa ? ` <a class="btn small" href="${base}/aa_matrix.csv?dl=1">⬇ AA 矩阵 CSV</a>` +
        ` <a class="btn small" href="${base}/identity_composite.png?dl=1">⬇ 复合热图 PNG</a>` +
        ` <a class="btn small" href="${base}/identity_composite.pdf?dl=1">⬇ 复合热图 PDF（矢量）</a>` : '') +
      ` <a class="btn small" href="${base}/nt_heatmap.png?dl=1">⬇ NT 热图 PNG</a>`;
    $('idtResult').style.display = '';
    $('idtResult').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) { alert(t('ui.common.connShort') + ': ' + e); }
}

async function sdtExactLoad() {
  try {
    const r = await fetch(`/tool_runs/${encodeURIComponent(sdRun)}/sdt_matrix.json`);
    if (!r.ok) { $('sdMeta').textContent = '加载失败: ' + (await r.json()).error; return; }
    sdData = await r.json();
    const base = `/tool_runs/${encodeURIComponent(sdRun)}`;
    const names = sdData.names, m = sdData.matrix;
    const stype = (sdData.seqtype || 'nt').toUpperCase();
    /* 交互热图：可切配色、悬浮看值、缩放（复用 examples.js 的 heatmap 渲染）。
       只有 heatmap 渲染函数与配色表都就绪才走交互；否则回退静态 PNG。 */
    const csSrc = (window.VPExamples && window.VPExamples.colorscales) || [];
    /* SDT 经典色阶：与 sdt_exact.SDT_COLORS 同步（后端 10 档红→黄硬分界）。
       Plotly 没有同名内置色标——按同一色值序列写成插值 stops（结构与其它
       matplotlib 系色标一致），并在配色表头部补一条 'SDT' 作为默认选中项。 */
    const SDT_CS = [[0, '#7d050d'], [0.1111, '#a81c0d'], [0.2222, '#cc3d0f'],
                    [0.3333, '#e66b17'], [0.4444, '#f5992b'], [0.5556, '#fcbf52'],
                    [0.6667, '#ffdb82'], [0.7778, '#ffedb2'], [0.8889, '#fff7db'],
                    [1, '#fffff0']];
    const csList = csSrc.length ? [['SDT', SDT_CS]].concat(csSrc) : [];
    const hasHM = !!(window.VPExamples && window.VPExamples.heatmap) &&
                  csList.length > 0;
    const wantLbl = { sdt: 'SDT', cividis: 'Cividis', viridis: 'Viridis',
                      'RdYlBu': 'RdYlBu', Spectral: 'Spectral', YlGnBu: 'YlGnBu',
                      coolwarm: 'coolwarm', magma: 'Magma' }[sdData.palette] || 'Viridis';
    const want = csList.some(c => c[0] === wantLbl) ? wantLbl : (csList[0]?.[0] || '');
    const distHtml = `
      <h4 style="font-size:14px;color:#1a5276;margin:14px 0 4px">Identity 分布</h4>
      <a href="${base}/sdt_distribution.png" target="_blank" style="display:inline-block;max-width:520px">
        <img src="${base}/sdt_distribution.png" style="max-width:520px;width:100%;border:1px solid #dfe5ec;border-radius:8px"></a>`;
    const staticHtml = `
      <h4 style="font-size:14px;color:#1a5276;margin:8px 0 8px">SDT 热图（聚类排序 · 三角）</h4>
      <a href="${base}/sdt_heatmap.png" target="_blank"><img src="${base}/sdt_heatmap.png"
         style="max-width:100%;max-height:720px;border:1px solid #dfe5ec;border-radius:8px"></a>`;
    if (hasHM) {
      /* 一次性写完整块：之后**不能**再 `innerHTML +=`——那会把 <select> 与
         Plotly 图整段重建，挂在旧节点上的 change 监听随之失效，
         表现为「改配色没反应」（曾出现的 bug）。 */
      $('sdExact').innerHTML = `
        <h4 style="font-size:14px;color:#1a5276;margin:8px 0 6px">SDT 热图（交互 · 聚类排序 · 三角）</h4>
        <div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:0 0 6px">
          <span class="hint">配色</span>
          <select id="sdReadyCs">${csList.map(c =>
            `<option value="${c[0]}" ${c[0] === want ? 'selected' : ''}>${c[0]}</option>`).join('')}</select>
          <span class="hint">悬浮看值 · 可缩放 · 下载见下方</span>
        </div>
        <div id="sdReadyHeat" style="width:100%;overflow:auto;background:#fff;border-radius:8px"></div>` +
        distHtml;
      const md = { names: names, nt: m, aa: sdData.aa || null, seqtype: stype.toLowerCase(), n: names.length };
      const sel = $('sdReadyCs');
      const sdDraw = () => {
        const k = sel.value;
        const c = csList.filter(x => x[0] === k)[0] || csList[0] || ['', 'YlGnBu'];
        window.VPExamples.heatmap($('sdReadyHeat'), md, c[1], stype.toLowerCase());
      };
      if (sel) sel.addEventListener('change', sdDraw);
      try {
        sdDraw();
      } catch (e) {
        /* 兜底：交互渲染失败时退回静态 PNG，绝不把错误抛给用户 */
        $('sdExact').innerHTML = staticHtml + distHtml;
      }
    } else {
      $('sdExact').innerHTML = staticHtml + distHtml;
    }
    $('sdDl').innerHTML =
      `<a class="btn small" href="${base}/sdt_matrix.csv?dl=1">⬇ 矩阵 CSV</a>` +
      ` <a class="btn small" href="${base}/sdt_heatmap.png?dl=1">⬇ 热图 PNG</a>` +
      ` <a class="btn small" href="${base}/sdt_heatmap.pdf?dl=1">⬇ 热图 PDF（矢量）</a>` +
      ` <a class="btn small" href="${base}/sdt_distribution.png?dl=1">⬇ 分布图 PNG</a>` +
      ` <a class="btn small" href="${base}/sdt_distribution.pdf?dl=1">⬇ 分布图 PDF（矢量）</a>`;
    $('sdTable').innerHTML = identityMatrixTable(names, m);
    const st = (sdData.seqtype || 'nt').toUpperCase();
    $('sdMeta').textContent =
      `${names.length} 条序列（${st}） · ${sdData.pairs} 对独立比对 · SDT v1.3 口径` +
      (sdData.aligned ? ' · 已比对直算' : '') +
      (sdData.palette === 'sdt' ? ' · SDT 经典色阶' : ` · ${sdData.palette} 色阶`);
    $('sdExact').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) { $('sdMeta').textContent = '无法连接: ' + esc(e); }
}

async function sdtToolLoadGroups() {
  const s = _v('sdSample'), sel = $('sdGroup');
  if (!s) { sel.innerHTML = '<option value="">（先选样品）</option>'; return; }
  const list = await fetchSampleList();
  const item = list.find(x => x.sample === s);
  sel.innerHTML = (item?.groups || []).map(g =>
    `<option value="${esc(g.group)}">${esc(g.group)}</option>`).join('')
    || '<option value="">（该样品没有矩阵/比对）</option>';
  if (sel.value) sdtToolView();
}

async function sdtToolView() {
  const s = _v('sdSample'), g = _v('sdGroup');
  if (!s || !g) return;
  $('sdMeta').textContent = '计算/加载矩阵中…';
  try {
    const r = await fetch(`/api/sdt/data?sample=${encodeURIComponent(s)}&group=${encodeURIComponent(g)}`);
    if (!r.ok) { sdData = null; $('sdMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    sdData = await r.json();
    sdRun = null;
    $('sdDl').innerHTML = '';
    sdtToolRender();
  } catch (e) { $('sdMeta').textContent = '无法连接: ' + esc(e); }
}

function sdtToolRender() {
  const box = $('sdHeat'), meta = $('sdMeta');
  if (!sdData || !sdData.names?.length) {
    box.style.display = 'none'; meta.textContent = ''; return;
  }
  const zmin = identityHeatmap('sdHeat', sdData.names, sdData.matrix,
    `SDT 全长成对 identity 矩阵（${esc(sdRun || (_v('sdSample') + ' · ' + _v('sdGroup')))}）`,
    _v('sdZmin'));
  box.style.display = '';
  $('sdTable').innerHTML = identityMatrixTable(sdData.names, sdData.matrix);
  meta.textContent =
    `${sdData.names.length} 条序列 · 色标 ${zmin.toFixed(0)}–100%` +
    (sdData.source ? ` · 数据来源 ${sdData.source}` : '');
}

async function sdtToolInit() {
  const sel = $('sdSample');
  if (!sel) return;
  const list = await fetchSampleList();
  sel.innerHTML = list.length
    ? list.map(x => `<option value="${esc(x.sample)}">${esc(sampleLabel(x))}</option>`).join('')
    : '<option value="">（暂无矩阵/比对，先跑流程或集合建树）</option>';
  if (sel.value) sdtToolLoadGroups();
}

alPickInit();
treeToolInit();
sdtToolInit();

