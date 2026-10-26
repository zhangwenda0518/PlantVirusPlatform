/* ==========================================================================
 * CDS / PEP 提取产物（独立模块）
 * 页面：/cds-export   接口：/api/cds/table · /api/cds/selection/save · /api/cds/export
 * 依赖 app.js 提供的全局：$ / esc / t
 * ========================================================================== */
let _cpData = null;                 // /api/cds/table 返回值
let _cpSel = new Map();             // rid -> 基因名（人工确认后）

/* ---------- 集合下拉（本模块自拉，不依赖工具页的批量加载） ---------- */
async function cdsLoadCollections() {
  const sel = $('cpColl');
  if (!sel) return;
  try {
    const r = await fetch('/api/gb/collections');
    if (!r.ok) return;
    const cols = await r.json();
    const cur = sel.value;
    sel.innerHTML = `<option value="">${t('cd.selPlaceholder', '（选择集合）')}</option>` +
      cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}` +
        `${t('cd.nRecPre', '（')}${c.n_records}${t('cd.nRecPost', ' 条）')}</option>`).join('');
    if (cur && cols.some(c => c.name === cur)) sel.value = cur;
    else if (cols.length === 1) sel.value = cols[0].name;
  } catch (e) { /* 集合列表拉取失败不阻塞页面 */ }
}

/* ---------- 解析并挑选 ---------- */
async function cdsPickOpen() {
  const name = ($('cpColl') && $('cpColl').value) || '';
  const msg = $('cpMsg');
  if (!name) { if (msg) msg.textContent = t('cd.needColl', '先选择集合'); return; }
  const panel = $('cdsPickPanel');
  const tbl = $('cpTable');
  if (panel) panel.style.display = '';
  if (tbl) tbl.innerHTML = `<p class="hint">${t('c.parsing', '解析中…')}</p>`;
  try {
    const r = await fetch('/api/cds/table?name=' + encodeURIComponent(name));
    if (!r.ok) {
      const e = await r.json().catch(() => ({}));
      if (tbl) tbl.innerHTML = '<p class="hint">' + esc(e.error || t('cd.parseFail', '解析失败')) + '</p>';
      return;
    }
    _cpData = await r.json();
    _cpSel = new Map();
    const saved = (_cpData.selection && _cpData.selection.items) || [];
    saved.forEach(x => { if (x && x.rid) _cpSel.set(x.rid, x.gene || 'gene'); });
    const fill = (id, vals) => {
      const el = $(id); if (!el) return;
      const cur = el.value;
      el.innerHTML = `<option value="">${t('cd.all', '（全部）')}</option>` +
        vals.map(v => `<option value="${esc(v)}">${esc(v)}</option>`).join('');
      el.value = cur;
    };
    fill('cpFam', _cpData.families || []);
    fill('cpGen', _cpData.genera || []);
    cdsPickRender();
    if (msg) msg.textContent = t('cd.parsedPre', '已解析 ') + _cpData.n_cds +
      t('cd.parsedCds', ' 条 CDS / ') + _cpData.n_virus + t('cd.parsedVirus', ' 个病毒') +
      (saved.length ? t('cd.parsedRestore', '，恢复上次选择 ') + saved.length + t('cd.parsedRestoreEnd', ' 条') : '');
  } catch (e) {
    if (tbl) tbl.innerHTML = `<p class="hint">${t('c.connFail', '无法连接: ')}` + esc(e) + '</p>';
  }
}

/* ---------- 过滤 ---------- */
function isHypothetical(r) {
  const s = ((r.product || '') + ' ' + (r.gene || '')).toLowerCase();
  return s.includes('hypothetical') || s.includes('unknown') ||
         s.includes('putative uncharacterized') ||
         /^(hp|orf\d+)$/.test((r.gene || '').toLowerCase());
}

function cdsPickRows() {
  if (!_cpData) return [];
  const fam = ($('cpFam') && $('cpFam').value) || '';
  const gen = ($('cpGen') && $('cpGen').value) || '';
  const showHyp = $('cpShowHyp') && $('cpShowHyp').checked;
  return _cpData.rows.filter(r => {
    if (fam && r.family !== fam) return false;
    if (gen && r.genus !== gen) return false;
    if (!showHyp && isHypothetical(r)) return false;
    return true;
  });
}

/* ---------- 渲染主表（按病毒分组平铺） ---------- */
function cdsPickRender() {
  const box = $('cpTable');
  if (!box || !_cpData) return;
  const rows = cdsPickRows();
  const cnt = $('cpCount');
  if (cnt) cnt.textContent = t('cd.showPre', '显示 ') + rows.length + ' / ' + _cpData.n_cds +
    t('cd.showMid', ' 条，已选 ') + _cpSel.size + t('cd.showTail', ' 条');
  if (!rows.length) { box.innerHTML = `<p class="hint">${t('c.noMatch', '无匹配记录')}</p>`; cdsPickSelRender(); return; }

  const groups = new Map();
  rows.forEach(r => {
    if (!groups.has(r.virus)) groups.set(r.virus, []);
    groups.get(r.virus).push(r);
  });
  let html = '<table style="width:100%;border-collapse:collapse;font-size:13px">' +
    '<thead><tr><th style="width:34px"></th>' +
    `<th style="text-align:left">${t('cd.thGene', '基因名')}</th>` +
    `<th style="text-align:left">${t('cd.thProduct', '产物')}</th>` +
    `<th style="width:78px">${t('cd.thStart', '起始')}</th>` +
    `<th style="width:78px">${t('cd.thEnd', '终止')}</th>` +
    `<th style="width:74px">${t('cd.thLen', '长度')}</th>` +
    `<th style="width:44px">${t('cd.thStrand', '链')}</th></tr></thead><tbody>`;
  for (const [virus, list] of groups) {
    const nsel = list.filter(r => _cpSel.has(r.rid)).length;
    const g0 = list[0];
    html += `<tr><td colspan="7" style="padding:6px 8px;background:var(--line-100,#f2f2f2);` +
      `font-weight:600">🦠 ${esc(virus)} <span class="hint" style="font-weight:400">` +
      `${esc(g0.family || '')}${g0.genus ? ' / ' + esc(g0.genus) : ''} · ${list.length}` +
      `${t('cd.nRecPlain', ' 条')}` +
      `${nsel ? ` · ${t('cd.picked', '已选')} ${nsel}` : ''}</span></td></tr>`;
    list.forEach(r => {
      const on = _cpSel.has(r.rid);
      const showGene = on ? _cpSel.get(r.rid) : r.gene;
      html += '<tr style="border-bottom:1px solid var(--line-100,#eee)">' +
        `<td style="text-align:center"><input type="checkbox" ${on ? 'checked' : ''} ` +
        `onchange="cdsPickToggle('${esc(r.rid)}', this.checked)"></td>` +
        `<td class="cp-gene" data-rid="${esc(r.rid)}" ondblclick="cdsPickEdit(this)" ` +
        `title="${t('cd.dblRename', '双击可改名')}" style="cursor:text">${esc(showGene)}</td>` +
        `<td class="hint">${esc(r.product || '')}</td>` +
        `<td style="text-align:right">${r.start}</td>` +
        `<td style="text-align:right">${r.end}</td>` +
        `<td style="text-align:right">${r.length}</td>` +
        `<td style="text-align:center">${esc(r.strand)}</td></tr>`;
    });
  }
  html += '</tbody></table>';
  box.innerHTML = html;
  cdsPickSelRender();
}

function cdsPickToggle(rid, on) {
  if (on) {
    const row = (_cpData.rows || []).find(r => r.rid === rid);
    _cpSel.set(rid, (row && row.gene) || 'gene');
  } else {
    _cpSel.delete(rid);
  }
  cdsPickRender();
}

/* 双击改基因名；改名即纳入选择（否则改动只留在 DOM 上，导出取不到） */
function cdsPickEdit(td) {
  const rid = td.dataset.rid;
  const cur = _cpSel.has(rid) ? _cpSel.get(rid) : td.textContent.trim();
  const inp = document.createElement('input');
  inp.type = 'text';
  inp.value = cur;
  inp.style.cssText = 'width:96%;font:inherit;padding:1px 3px';
  let done = false;
  const commit = (ok) => {
    if (done) return;
    done = true;
    const v = (inp.value || '').trim();
    if (ok && v) _cpSel.set(rid, v);
    cdsPickRender();
  };
  inp.onblur = () => commit(true);
  inp.onkeydown = (e) => {
    if (e.key === 'Enter') { e.preventDefault(); commit(true); }
    if (e.key === 'Escape') { e.preventDefault(); commit(false); }
  };
  td.innerHTML = '';
  td.appendChild(inp);
  inp.focus();
  inp.select();
}

/* ---------- 已选面板（按基因名归并） ---------- */
function cdsPickSelRender() {
  const box = $('cpSel');
  if (!box) return;
  if (!_cpSel.size) { box.innerHTML = `<p class="hint">${t('cd.none', '（未选择）')}</p>`; return; }
  const byGene = new Map();
  for (const [rid, gene] of _cpSel) {
    const r = (_cpData.rows || []).find(x => x.rid === rid) || {};
    if (!byGene.has(gene)) byGene.set(gene, []);
    byGene.get(gene).push({ rid, r });
  }
  let html = '<table style="width:100%;border-collapse:collapse;font-size:12px">';
  for (const [gene, list] of [...byGene.entries()].sort()) {
    html += `<tr><td colspan="3" style="padding:4px 6px;background:var(--line-100,#f2f2f2);` +
      `font-weight:600">📕 ${esc(gene)} <span class="hint" style="font-weight:400">` +
      `${list.length}${t('cd.nRecPlain', ' 条')}</span></td></tr>`;
    list.forEach(({ rid, r }) => {
      html += '<tr style="border-bottom:1px solid var(--line-100,#eee)">' +
        `<td style="padding:3px 6px" title="${esc(r.virus || '')}">${esc((r.virus || '').slice(0, 22))}</td>` +
        `<td style="text-align:right;white-space:nowrap">${r.length || 0}nt</td>` +
        `<td style="width:22px"><button class="btn small" style="padding:0 5px" ` +
        `onclick="cdsPickToggle('${esc(rid)}', false)">×</button></td></tr>`;
    });
  }
  html += '</table>';
  box.innerHTML = html;
}

function cdsPickClear() { _cpSel.clear(); cdsPickRender(); }

function cdsPickItems() {
  return [..._cpSel.entries()].map(([rid, gene]) => ({ rid, gene }));
}

/* ---------- 保存 / 导出 ---------- */
async function cdsPickSave() {
  const name = ($('cpColl') && $('cpColl').value) || '';
  const msg = $('cpMsg');
  if (!name) return;
  try {
    const r = await fetch('/api/cds/selection/save', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, items: cdsPickItems() })
    });
    const j = await r.json();
    if (msg) msg.textContent = r.ok
      ? t('cd.savedPre', '✔ 已保存 ') + j.n + t('cd.savedEnd', ' 条')
      : (t('cd.saveFail', '保存失败: ') + (j.error || ''));
  } catch (e) { if (msg) msg.textContent = t('cd.saveFail', '保存失败: ') + e; }
}

async function cdsPickExport() {
  const name = ($('cpColl') && $('cpColl').value) || '';
  const msg = $('cpMsg');
  if (!name) return;
  if (!_cpSel.size) { if (msg) msg.textContent = t('cd.needPick', '请先勾选 CDS'); return; }
  if (msg) msg.textContent = t('cd.exporting', '导出中…');
  try {
    // async：导出要重新全量解析集合每个 .gb 并现场翻译，科级集合分钟级。
    // 改后台任务 → 任务中心可见进度/日志/取消；完成后读结果预览里的字段。
    const r = await fetch('/api/cds/export', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, items: cdsPickItems(), async: true })
    });
    const j0 = await r.json();
    if (!r.ok) { if (msg) msg.textContent = t('cd.exportFail', '导出失败: ') + (j0.error || ''); return; }
    let snap = null;
    for (;;) {
      await new Promise(res => setTimeout(res, 1500));
      const sr = await fetch('/api/task/' + encodeURIComponent(j0.task) + '?log_lines=10');
      snap = await sr.json();
      if (msg) msg.textContent = `${t('cd.exporting', '导出中…')} ${trTaskMsg(snap.stage || '')} `
        + (snap.pct != null ? Math.round((snap.pct || 0) * 100) + '%' : '');
      if (snap.status !== 'running') break;
    }
    if (snap.status !== 'done') {
      if (msg) msg.textContent = t('cd.exportFail', '导出失败: ') + (snap.error || snap.status);
      return;
    }
    const j = (snap.result && snap.result.stats) || {};
    const nGenes = (j.n_genes != null) ? j.n_genes : (j.genes || []).length;
    const nMiss = (j.n_missing != null) ? j.n_missing : (j.missing || []).length;
    if (msg) msg.innerHTML = t('cd.okPre', '✔ 导出 ') + j.n_cds + t('cd.okCds', ' 条 CDS / ') +
      j.n_pep + t('cd.okPep', ' 条蛋白，') + nGenes + t('cd.okGenes', ' 个基因 → ') +
      `<span title="${esc(j.dir || '')}">extract/selected/</span>` +
      (nMiss ? t('cd.okMissPre', '（') + nMiss + t('cd.okMissEnd', ' 条未找到）') : '');
  } catch (e) { if (msg) msg.textContent = t('cd.exportFail', '导出失败: ') + e; }
}

/* ---------- 初始化 ---------- */
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', cdsLoadCollections);
} else {
  cdsLoadCollections();
}
