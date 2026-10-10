/* miRNA 靶标预测（t-mirna）前端逻辑 —— PAmiRDB 展示层移植。
   独立文件原因：tools.html 尾部内联脚本曾被并行会话回滚两次，卡片 HTML 保留
   但整块 JS 被删；放静态文件后模板只依赖一行 <script src>，改动互不影响。
   依赖（须先于本文件加载）：i18n.js(t) → app.js($, esc) → app-jobs.js(_submitTool)。 */
'use strict';

const mirS = { results: [], filtered: [], page: 1, pageSize: 20,
               sortCol: 'best_energy', sortAsc: true, run: null };
const MIR_CONS_ORDER = {NONE:0, LOW:1, 'LOW+':2, MEDIUM:3, 'MEDIUM+':4,
                        HIGH:5, 'HIGH+':6, 'VERY HIGH':7};
/* 引擎勾选（PAmiRDB 风格）：miRanda 必选锁死；其余按勾选参与共识投票 */
const MIR_ENGINE_KEYS = ['psrna', 'rnahybrid', 'rna22', 'tapir', 'psrobot'];

function mirSelectedAlgorithms() {
  return ['miranda', ...MIR_ENGINE_KEYS.filter(k => {
    const el = $('mir-alg-' + k);
    return el && el.checked;
  })];
}

async function runMirna(btn) {
  const params = {
    mirna_text: $('mirna_text').value.trim(),
    virus_text: $('virus_text').value.trim(),
    mirna_fasta: $('mirna_fasta').value.trim(),
    virus_fasta: $('virus_fasta').value.trim(),
    algorithms: mirSelectedAlgorithms().join(',')
  };
  if (!params.mirna_text && !params.mirna_fasta) {
    alert(t('ui.mi.needMirna')); return;
  }
  if (!params.virus_text && !params.virus_fasta) {
    alert(t('ui.mi.needVirus')); return;
  }
  _submitTool('mirna', params, null, btn, '预测中…', loadMirnaResult);
}

async function loadMirnaResult(run) {
  try {
    const s = await fetch(`/tool_runs/${run}/mirna_target/summary.json`).then(r => r.json());
    const rows = await fetch(`/tool_runs/${run}/mirna_target/results.json`).then(r => r.json());
    mirS.run = run;
    mirS.results = Array.isArray(rows) ? rows : [];
    const nAlg = (s.algorithms || []).length;
    $('mirStats').textContent =
      `miRNA ${s.n_mirna} × 病毒 ${s.n_virus} = ${s.n_pairs} 对 · 命中 ${s.n_hits} · ` +
      `${nAlg || (s.full_mode ? 6 : 3)} 引擎 · ${s.elapsed_sec}s`;
    $('mirResult').style.display = '';
    $('mir-detail').style.display = 'none';
    mirResetFilters();
  } catch (e) {
    alert(t('ui.mi.readResultFail') + ': ' + e);
  }
}

function mirConsClass(c) {
  if (!c || c === 'NONE') return 'mir-c-none';
  if (c.indexOf('VERY') >= 0) return 'mir-c-vhigh';
  if (c.indexOf('HIGH') >= 0) return 'mir-c-high';
  if (c.indexOf('MEDIUM') >= 0) return 'mir-c-medium';
  return 'mir-c-low';
}

function mirApplyFilters() {
  const fTier = $('mir-f-tier').value;
  const fSeed = $('mir-f-seed').value;
  const fE = parseFloat($('mir-f-energy').value);
  const fId = parseFloat($('mir-f-identity').value);
  const fMm = parseInt($('mir-f-mm').value, 10);
  const fSearch = $('mir-f-search').value.toLowerCase().trim();
  const fQ2 = $('mir-f-query2') && $('mir-f-query2').checked;
  const MIN_TIER = { VHIGH: 7, HIGH: 5, MEDIUM: 3, LOW: 1 };
  mirS.filtered = mirS.results.filter(r => {
    if (fTier && (MIR_CONS_ORDER[r.consensus_level] || 0) < (MIN_TIER[fTier] || 0)) return false;
    if (fSeed === '1' && r.seed_perfect_wc !== 7) return false;
    if (fQ2 && !mirIsQuery2(r)) return false;
    if (!isNaN(fE) && (r.best_energy ?? 0) > fE) return false;
    if (!isNaN(fId) && (r.best_identity || 0) < fId) return false;
    if (!isNaN(fMm) && (r.num_mismatches ?? 0) > fMm) return false;
    if (fSearch && (r.mirna_id || '').toLowerCase().indexOf(fSearch) < 0 &&
        (r.target_virus_id || '').toLowerCase().indexOf(fSearch) < 0) return false;
    return true;
  });
  mirS.page = 1;
  mirSortResults();
}

function mirResetFilters() {
  ['mir-f-energy', 'mir-f-identity', 'mir-f-mm', 'mir-f-search'].forEach(id => { $(id).value = ''; });
  $('mir-f-tier').value = ''; $('mir-f-seed').value = '';
  if ($('mir-f-query2')) $('mir-f-query2').checked = false;
  mirApplyFilters();
}

function mirSortBy(col) {
  if (mirS.sortCol === col) mirS.sortAsc = !mirS.sortAsc;
  else { mirS.sortCol = col; mirS.sortAsc = col !== 'best_energy'; }
  mirS.page = 1;
  mirSortResults();
}

function mirSortResults() {
  const col = mirS.sortCol, asc = mirS.sortAsc;
  mirS.filtered.sort((a, b) => {
    let av = a[col], bv = b[col];
    if (col === 'consensus_level') { av = MIR_CONS_ORDER[av] || 0; bv = MIR_CONS_ORDER[bv] || 0; }
    if (av == null) av = asc ? Infinity : -Infinity;
    if (bv == null) bv = asc ? Infinity : -Infinity;
    return asc ? (av > bv ? 1 : av < bv ? -1 : 0) : (av < bv ? 1 : av > bv ? -1 : 0);
  });
  mirRenderTable();
}

function mirRenderTable() {
  const total = mirS.filtered.length;
  const totalPages = Math.max(Math.ceil(total / mirS.pageSize), 1);
  if (mirS.page > totalPages) mirS.page = totalPages;
  if (mirS.page < 1) mirS.page = 1;
  const start = (mirS.page - 1) * mirS.pageSize;
  const items = mirS.filtered.slice(start, start + mirS.pageSize);
  $('mir-f-count').textContent = `${total} / ${mirS.results.length} 条（第 ${mirS.page}/${totalPages} 页）`;
  const tbody = $('mir-tbody');
  if (!total) {
    tbody.innerHTML = '<tr><td colspan="9" style="text-align:center;padding:20px;color:#8892a6">当前过滤条件下没有结果</td></tr>';
    $('mir-pager').style.display = 'none';
    return;
  }
  tbody.innerHTML = items.map((r, i) => {
    const idx = start + i;
    const cLevel = r.consensus_level || 'NONE';
    const algoStr = r.consensus_detail || 'miRanda';
    const seedCls = r.seed_perfect_wc === 7 ? 'mir-seed-ok' : 'mir-seed-no';
    return '<tr onclick="mirShowDetail(' + idx + ')" style="cursor:pointer">'
      + '<td style="font-weight:600;color:#2c6e49">' + esc(r.mirna_id) + '</td>'
      + '<td style="font-size:10px">' + esc(r.target_virus_name || r.target_virus_id || '') + '</td>'
      + '<td><span class="mir-badge ' + mirConsClass(cLevel) + '">' + esc(cLevel) + '</span></td>'
      + '<td>' + (r.best_energy != null ? r.best_energy.toFixed(1) : '-') + '</td>'
      + '<td>' + (r.best_identity != null ? r.best_identity.toFixed(0) + '%' : '-') + '</td>'
      + '<td>' + (r.num_mismatches != null ? r.num_mismatches : '-') + '</td>'
      + '<td class="' + seedCls + '">' + (r.seed_perfect_wc != null ? r.seed_perfect_wc + '/7' : '-') + '</td>'
      + '<td style="font-size:10px;color:#67748e">' + esc(algoStr) + '</td>'
      + '<td><button onclick="event.stopPropagation();mirShowDetail(' + idx + ')" '
      + 'style="padding:3px 8px;font-size:10px;border:1px solid #2c6e49;background:none;'
      + 'color:#2c6e49;border-radius:4px;cursor:pointer">View</button></td>'
      + '</tr>';
  }).join('');
  const pager = $('mir-pager');
  if (totalPages > 1) {
    pager.style.display = 'flex';
    $('mir-page-info').textContent = `第 ${mirS.page} / ${totalPages} 页`;
    const btns = pager.querySelectorAll('button');
    btns[0].disabled = btns[1].disabled = (mirS.page <= 1);
    btns[2].disabled = btns[3].disabled = (mirS.page >= totalPages);
  } else pager.style.display = 'none';
}

function mirGoPage(p) {
  const totalPages = Math.max(Math.ceil(mirS.filtered.length / mirS.pageSize), 1);
  mirS.page = Math.max(1, Math.min(totalPages, p));
  mirRenderTable();
  $('mir-table').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function mirScoreItem(v, l) {
  return '<div class="mir-s-item"><div class="v">' + v + '</div><div class="l">' + l + '</div></div>';
}

function mirShowDetail(idx) {
  const r = mirS.filtered[idx];
  if (!r) return;
  document.querySelectorAll('#mir-tbody tr').forEach((tr, i) => {
    tr.className = i === idx - (mirS.page - 1) * mirS.pageSize ? 'mir-active' : '';
  });
  $('mir-detail').style.display = '';
  $('mir-det-title').textContent = r.mirna_id + ' → ' + (r.target_virus_name || r.target_virus_id || 'target');
  const cLevel = r.consensus_level || 'NONE';
  $('mir-det-scores').innerHTML =
    '<div class="mir-s-item mir-span-all"><div class="v"><span class="mir-badge ' + mirConsClass(cLevel) + '">'
    + esc(cLevel) + '</span>&nbsp; ' + (r.num_algorithms || 0) + ' 个引擎通过</div>'
    + '<div class="l">' + esc(r.consensus_detail || '') + '</div></div>'
    + mirScoreItem(r.best_energy != null ? r.best_energy.toFixed(1) + ' kcal/mol' : '-', 'Energy (MFE)')
    + mirScoreItem(r.seed_perfect_wc != null ? r.seed_perfect_wc + '/7 WC' : '-', 'Seed 质量')
    + mirScoreItem(r.best_identity != null ? r.best_identity.toFixed(1) + '%' : '-', 'Identity')
    + mirScoreItem(r.num_mismatches != null ? r.num_mismatches : '-', '错配')
    + mirScoreItem(r.miranda_score != null ? r.miranda_score.toFixed(0) : '-', 'miRanda Score')
    + mirScoreItem(r.target_start != null ? (r.target_start + 1) + '–' + r.target_end : '-', '位点 (nt)')
    + mirScoreItem(r.genome_length != null ? r.genome_length + ' nt' : '-', '基因组长度');
  $('mir-tab-aln').innerHTML = r.alignment
    ? mirDuplexSVG(r.alignment) + '<div class="mir-aln-wrap">' + mirAlnHTML(r.alignment) + '</div>'
    : '<p style="color:#8892a6;padding:20px;text-align:center">无比对数据</p>';
  $('mir-tab-algos').innerHTML = mirAlgoTable(r);
  $('mir-tab-map').innerHTML = mirPosMap(r);
  $('mir-tab-struct').innerHTML = mirStructSVG(r.mirna_seq || '', r.mirna_id || '');
  mirShowTab('aln');
  $('mir-detail').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

function mirShowTab(tp) {
  ['aln', 'algos', 'map', 'struct'].forEach(x => {
    $('mir-tab-' + x).className = x === tp ? 'mir-tab-pane on' : 'mir-tab-pane';
  });
  document.querySelectorAll('#mir-detail .mir-tab-btn').forEach((b, i) => {
    b.className = (['aln', 'algos', 'map', 'struct'][i] === tp) ? 'mir-tab-btn on' : 'mir-tab-btn';
  });
}

/* —— 双链 SVG（移植 PAmiRDB renderDuplexSVG）—— */
function mirDuplexSVG(aln) {
  const lines = aln.split('\n');
  if (lines.length < 3) return '';
  let mCore = lines[0].replace(/^miRNA\s+5'\s*/, '').replace(/\s*3'$/, '');
  let pCore = lines[1].replace(/^\s+/, '');
  let tCore = lines[2].replace(/^target\s+3'\s*/, '').replace(/\s*5'$/, '');
  const L = Math.max(mCore.length, pCore.length, tCore.length);
  while (mCore.length < L) mCore += ' ';
  while (pCore.length < L) pCore += ' ';
  while (tCore.length < L) tCore += ' ';
  const CELL = 20, PAD_L = 52, B = 8, W = PAD_L + L * CELL + 24, H = 98;
  const baseFill = { A: '#e74c3c', U: '#3498db', G: '#27ae60', C: '#e67e22', T: '#3498db' };
  let mPos = 0, seedX1 = -1, seedX2 = -1;
  for (let i = 0; i < L; i++) {
    const c = mCore[i];
    if (c && c !== '-' && c !== ' ') {
      mPos++;
      if (mPos === 2) seedX1 = PAD_L + i * CELL;
      if (mPos === 8) seedX2 = PAD_L + (i + 1) * CELL;
    }
  }
  let bg = '', circles = '', connectors = '', posNums = '';
  if (seedX1 >= 0 && seedX2 > seedX1) {
    bg = '<rect x="' + seedX1 + '" y="4" width="' + (seedX2 - seedX1) + '" height="72" fill="rgba(241,196,15,0.08)" rx="3"/>'
       + '<rect x="' + seedX1 + '" y="4" width="' + (seedX2 - seedX1) + '" height="72" fill="none" stroke="rgba(241,196,15,0.35)" stroke-width="1" rx="3" stroke-dasharray="3,2"/>';
  }
  mPos = 0;
  for (let i = 0; i < L; i++) {
    const cx = PAD_L + i * CELL + CELL / 2;
    const mc = (mCore[i] || ' ').toUpperCase().replace('T', 'U');
    const tc = (tCore[i] || ' ').toUpperCase().replace('T', 'U');
    const pc = pCore[i] || ' ';
    if (mc !== '-' && mc !== ' ') {
      circles += '<circle cx="' + cx + '" cy="14" r="' + B + '" fill="' + (baseFill[mc] || '#888') + '"/>'
               + '<text x="' + cx + '" y="18" text-anchor="middle" font-size="9" fill="#fff" font-weight="700">' + mc + '</text>';
      mPos++;
    } else if (mc === '-') {
      circles += '<text x="' + cx + '" y="18" text-anchor="middle" font-size="11" fill="#555">–</text>';
    }
    if (tc !== '-' && tc !== ' ') {
      circles += '<circle cx="' + cx + '" cy="64" r="' + B + '" fill="' + (baseFill[tc] || '#888') + '"/>'
               + '<text x="' + cx + '" y="68" text-anchor="middle" font-size="9" fill="#fff" font-weight="700">' + tc + '</text>';
    } else if (tc === '-') {
      circles += '<text x="' + cx + '" y="68" text-anchor="middle" font-size="11" fill="#555">–</text>';
    }
    if (pc === '|') {
      connectors += '<line x1="' + cx + '" y1="23" x2="' + cx + '" y2="55" stroke="#27ae60" stroke-width="1.5"/>';
    } else if (pc === ':') {
      connectors += '<line x1="' + cx + '" y1="23" x2="' + cx + '" y2="55" stroke="#e67e22" stroke-width="1.5" stroke-dasharray="2.5,2"/>';
    } else if (pc === '.') {
      connectors += '<line x1="' + (cx - 3) + '" y1="30" x2="' + (cx + 3) + '" y2="48" stroke="#e74c3c" stroke-width="1.5"/>'
                  + '<line x1="' + (cx + 3) + '" y1="30" x2="' + (cx - 3) + '" y2="48" stroke="#e74c3c" stroke-width="1.5"/>';
    }
    if (mc !== '-' && mc !== ' ' && mPos % 5 === 0) {
      posNums += '<text x="' + cx + '" y="86" text-anchor="middle" font-size="8" fill="#555">' + mPos + '</text>';
    }
  }
  const labels = '<text x="2" y="18" font-size="9" fill="#8b949e" font-family="monospace">5\'</text>'
    + '<text x="' + (W - 10) + '" y="18" text-anchor="end" font-size="9" fill="#8b949e" font-family="monospace">3\'</text>'
    + '<text x="2" y="38" font-size="8" fill="#6e7681">miRNA</text>'
    + '<text x="2" y="52" font-size="8" fill="#6e7681">target</text>'
    + '<text x="2" y="68" font-size="9" fill="#8b949e" font-family="monospace">3\'</text>'
    + '<text x="' + (W - 10) + '" y="68" text-anchor="end" font-size="9" fill="#8b949e" font-family="monospace">5\'</text>';
  const svg = '<svg width="' + W + '" height="' + H + '" xmlns="http://www.w3.org/2000/svg" style="max-width:100%;display:block">'
    + '<rect width="' + W + '" height="' + H + '" fill="#0d1117"/>' + bg + connectors + circles + labels + posNums + '</svg>';
  return '<div class="mir-duplex-wrap"><h4>miRNA : Target 双链结构</h4>' + svg
    + '<div class="mir-duplex-legend">'
    + '<span><svg width="12" height="10"><line x1="6" y1="0" x2="6" y2="10" stroke="#27ae60" stroke-width="2"/></svg>Watson-Crick</span>'
    + '<span><svg width="12" height="10"><line x1="6" y1="0" x2="6" y2="10" stroke="#e67e22" stroke-width="2" stroke-dasharray="2,1.5"/></svg>Wobble (G:U)</span>'
    + '<span><svg width="12" height="10"><line x1="2" y1="0" x2="10" y2="10" stroke="#e74c3c" stroke-width="2"/><line x1="10" y1="0" x2="2" y2="10" stroke="#e74c3c" stroke-width="2"/></svg>错配</span>'
    + '<span><svg width="12" height="10"><rect x="2" y="0" width="8" height="10" fill="rgba(241,196,15,0.25)" stroke="rgba(241,196,15,0.6)" stroke-width="1"/></svg>Seed 区 (2–8)</span>'
    + '</div></div>';
}

/* —— 彩色比对照排（移植 PAmiRDB renderAlnHTML）—— */
function mirAlnHTML(aln) {
  const lines = aln.split('\n');
  if (lines.length < 3) return esc(aln);
  // 2026-10-03 适配新帧：比对串改为 miRanda 官方书写约定（上行 miRNA 3'→5'，
  // 与官方 miranda.exe 输出一致）。兼容旧格式（上行 5'→3'）：按上行标签定帧，
  // 种子位（miRNA 5' 第 2-8 位）在新帧里自右向左数。
  const mLabel = lines[0].match(/^miRNA\s+([53])'/);
  const mir3first = mLabel ? mLabel[1] === '3' : false;
  const mTop = mir3first ? "3'" : "5'", mBot = mir3first ? "5'" : "3'";
  const mCore = lines[0].replace(/^miRNA\s+[53]'\s*/, '').replace(/\s*[53]'$/, '');
  const pCore = lines[1].replace(/^\s+/, '');
  const tCore = lines[2].replace(/^target\s+[53]'\s*/, '').replace(/\s*[53]'$/, '');
  const nMir = mCore.replace(/-/g, '').length;
  function colorLine(seq, pipes, isMirna) {
    let pos = 0, html = '';
    for (let i = 0; i < seq.length; i++) {
      const c = seq[i];
      if (c === '-') { html += '<span class="mir-a-gp">' + c + '</span>'; continue; }
      pos++;
      // 5' 起算位次：新帧上行左端是 3' 端，需换算
      const pos5 = (isMirna && mir3first) ? (nMir - pos + 1) : pos;
      const p = pipes[i] || ' ';
      const cls = p === '|' ? 'mir-a-wc' : p === ':' ? 'mir-a-wo' : p === '.' ? 'mir-a-mm' : '';
      const seedCls = (isMirna && pos5 >= 2 && pos5 <= 8) ? ' mir-a-seed' : '';
      html += '<span class="' + cls + seedCls + '">' + c + '</span>';
    }
    return html;
  }
  let pHtml = '';
  for (let k = 0; k < pCore.length; k++) {
    pHtml += pCore[k] === '|' ? '<span class="mir-a-wc">|</span>'
      : pCore[k] === ':' ? '<span class="mir-a-wo">:</span>' : pCore[k];
  }
  return "miRNA  " + mTop + " " + colorLine(mCore, pCore, true) + " " + mBot + "\n"
    + "       " + pHtml + "\ntarget " + mBot + " " + colorLine(tCore, pCore, false) + " " + mTop;
}

/* —— 各引擎状态表（移植 PAmiRDB renderAlgoTable）—— */
function mirAlgoTable(r) {
  const rows = [
    ['miRanda', r.miranda_score,
      'Score: ' + (r.miranda_score != null ? r.miranda_score.toFixed(0) : 'N/A')
      + ', Energy: ' + (r.miranda_energy != null ? r.miranda_energy.toFixed(1) : 'N/A')
      + ', Identity: ' + (r.best_identity != null ? r.best_identity.toFixed(1) + '%' : 'N/A')],
    ['psRNATarget', r.psrna_score,
      'Score: ' + (r.psrna_score != null ? r.psrna_score.toFixed(2) : 'N/A')
      + ', Expectation: ' + (r.psrna_expectation != null ? r.psrna_expectation.toFixed(2) : 'N/A')],
    ['RNAhybrid', r.rnahybrid_mfe,
      'MFE: ' + (r.rnahybrid_mfe != null ? r.rnahybrid_mfe.toFixed(1) : 'N/A')
      + ', Seed WC: ' + (r.rnahybrid_seed != null ? (r.rnahybrid_seed ? 'OK' : 'NO') : 'N/A')],
    ['RNA22', r.rna22_mfe,
      'MFE: ' + (r.rna22_mfe != null ? r.rna22_mfe.toFixed(1) : 'N/A')
      + ', Pattern: ' + (r.rna22_pattern_score != null ? r.rna22_pattern_score.toFixed(1) : 'N/A')],
    ['TAPIR', r.tapir_score,
      'Score: ' + (r.tapir_score != null ? r.tapir_score : 'N/A')
      + ', 可及性: ' + (r.tapir_accessibility != null ? r.tapir_accessibility : 'N/A')],
    ['psRobot', r.psrobot_score,
      'Score: ' + (r.psrobot_score != null ? r.psrobot_score : 'N/A')]
  ];
  return '<table class="mir-algo-tbl"><thead><tr><th>Algorithm</th><th>Status</th><th>Key Metrics</th></tr></thead><tbody>'
    + rows.map(row => {
      const cls = row[1] == null ? 'mir-fail' : (row[1] ? 'mir-pass' : 'mir-fail');
      const status = row[1] == null ? 'SKIP' : (row[1] ? 'PASS' : 'FAIL');
      return '<tr><td>' + row[0] + '</td><td class="' + cls + '">' + status + '</td>'
        + '<td style="color:#67748e">' + row[2] + '</td></tr>';
    }).join('') + '</tbody></table>';
}

/* —— miRNA 二级结构弧图（移植 PAmiRDB Nussinov MBP + rendermiRNAStructSVG）—— */
function mirCanPairRNA(a, b) {
  a = a.toUpperCase().replace('T', 'U'); b = b.toUpperCase().replace('T', 'U');
  if ((a === 'G' && b === 'U') || (a === 'U' && b === 'G')) return true;
  return ({ A: 'U', U: 'A', G: 'C', C: 'G' })[a] === b;
}

function mirNussinovFold(seq) {
  seq = seq.toUpperCase().replace(/T/g, 'U');
  const n = seq.length;
  const dp = [];
  for (let i = 0; i < n; i++) { dp[i] = []; for (let j = 0; j < n; j++) dp[i][j] = 0; }
  for (let d = 4; d < n; d++) {
    for (let i = 0; i + d < n; i++) {
      const j = i + d;
      dp[i][j] = Math.max(dp[i + 1][j], dp[i][j - 1]);
      if (mirCanPairRNA(seq[i], seq[j])) {
        const inner = (i + 1 <= j - 1) ? dp[i + 1][j - 1] : 0;
        dp[i][j] = Math.max(dp[i][j], inner + 1);
      }
      for (let k = i + 1; k < j; k++) {
        dp[i][j] = Math.max(dp[i][j], dp[i][k] + dp[k + 1][j]);
      }
    }
  }
  const pairs = {};
  function tb(i, j) {
    if (i >= j || !dp[i] || dp[i][j] === 0) return;
    if (i + 1 <= j && dp[i + 1][j] === dp[i][j]) { tb(i + 1, j); return; }
    if (dp[i][j - 1] === dp[i][j]) { tb(i, j - 1); return; }
    if (mirCanPairRNA(seq[i], seq[j])) {
      const inner = (i + 1 <= j - 1) ? dp[i + 1][j - 1] : 0;
      if (inner + 1 === dp[i][j]) { pairs[i] = j; pairs[j] = i; tb(i + 1, j - 1); return; }
    }
    for (let k = i + 1; k < j; k++) {
      if (dp[i][k] + dp[k + 1][j] === dp[i][j]) { tb(i, k); tb(k + 1, j); return; }
    }
  }
  if (n > 1) tb(0, n - 1);
  return pairs;
}

function mirStructSVG(seq, mirnaId) {
  if (!seq || seq.length < 8) {
    return '<p style="color:#8892a6;padding:20px;text-align:center">无 miRNA 序列</p>';
  }
  seq = seq.toUpperCase().replace(/T/g, 'U');
  const n = seq.length;
  const pairs = mirNussinovFold(seq);
  const numPairs = Object.keys(pairs).length / 2 | 0;
  const CELL = 22, PAD_L = 28, PAD_R = 20, seqY = 100, arcMaxH = 78;
  const W = PAD_L + n * CELL + PAD_R, H = 140;
  const bFill = { A: '#e74c3c', U: '#3498db', G: '#27ae60', C: '#e67e22' };
  let arcs = '';
  for (let i = 0; i < n; i++) {
    if (pairs[i] !== undefined && i < pairs[i]) {
      const j = pairs[i];
      const xi = PAD_L + i * CELL + CELL / 2, xj = PAD_L + j * CELL + CELL / 2;
      const arcH = Math.min((j - i) * CELL * 0.3, arcMaxH);
      const mx = (xi + xj) / 2, cy = seqY - arcH;
      arcs += '<path d="M ' + xi + ' ' + (seqY - 10) + ' Q ' + mx.toFixed(1) + ' ' + cy.toFixed(1)
        + ' ' + xj + ' ' + (seqY - 10) + '" fill="none" stroke="#2c6e49" stroke-width="1.8" stroke-opacity="0.65" stroke-linecap="round"/>';
    }
  }
  let circles = '', posNums = '';
  for (let i = 0; i < n; i++) {
    const cx = PAD_L + i * CELL + CELL / 2;
    const base = seq[i];
    const isPaired = pairs[i] !== undefined;
    circles += '<circle cx="' + cx + '" cy="' + seqY + '" r="9" fill="' + (bFill[base] || '#888') + '" opacity="' + (isPaired ? '1' : '0.55') + '"/>';
    circles += '<text x="' + cx + '" y="' + (seqY + 4) + '" text-anchor="middle" font-size="9" fill="#fff" font-weight="700">' + base + '</text>';
    if (i === 0 || (i + 1) % 5 === 0 || i === n - 1) {
      posNums += '<text x="' + cx + '" y="' + (seqY + 22) + '" text-anchor="middle" font-size="8" fill="#aaa">' + (i + 1) + '</text>';
    }
  }
  const ends = '<text x="' + (PAD_L - 6) + '" y="' + (seqY + 4) + '" text-anchor="end" font-size="9" fill="#aaa" font-style="italic">5\'</text>'
    + '<text x="' + (PAD_L + n * CELL + 6) + '" y="' + (seqY + 4) + '" text-anchor="start" font-size="9" fill="#aaa" font-style="italic">3\'</text>';
  const gcN = (seq.match(/[GC]/g) || []).length;
  const info = '<text x="' + PAD_L + '" y="14" font-size="9" fill="#8b949e">'
    + n + ' nt · GC: ' + (gcN / n * 100).toFixed(0) + '% · ' + numPairs + ' base pairs (Nussinov MBP)</text>';
  const svg = '<svg width="' + W + '" height="' + H + '" xmlns="http://www.w3.org/2000/svg" style="font-family:sans-serif;max-width:100%;display:block">'
    + info + arcs + circles + posNums + ends + '</svg>';
  let db = '';
  for (let i = 0; i < n; i++) db += (pairs[i] === undefined ? '.' : (i < pairs[i] ? '(' : ')'));
  return '<div style="background:#0d1117;border-radius:8px;padding:14px 16px;overflow-x:auto">'
    + '<div style="font-size:10px;font-weight:700;color:#8b949e;text-transform:uppercase;letter-spacing:.5px;margin-bottom:10px">'
    + esc(mirnaId) + ' — 二级结构弧图</div>'
    + '<div style="overflow-x:auto">' + svg + '</div>'
    + '<div style="margin-top:10px;font-family:Consolas,monospace;font-size:11px;white-space:nowrap;overflow-x:auto">'
    + '<span style="color:#e67e22">5\'</span><span style="color:#58a6ff"> ' + seq + ' </span><span style="color:#e67e22">3\'</span><br>'
    + '<span style="color:#555">   </span><span style="color:#7ee787"> ' + db + '</span></div>'
    + '<div style="margin-top:10px;display:flex;gap:12px;flex-wrap:wrap;align-items:center">'
    + '<span style="font-size:9px;color:#8b949e"><span style="display:inline-block;width:9px;height:9px;border-radius:50%;background:#e74c3c;vertical-align:middle;margin-right:3px"></span>A</span>'
    + '<span style="font-size:9px;color:#8b949e"><span style="display:inline-block;width:9px;height:9px;border-radius:50%;background:#3498db;vertical-align:middle;margin-right:3px"></span>U</span>'
    + '<span style="font-size:9px;color:#8b949e"><span style="display:inline-block;width:9px;height:9px;border-radius:50%;background:#27ae60;vertical-align:middle;margin-right:3px"></span>G</span>'
    + '<span style="font-size:9px;color:#8b949e"><span style="display:inline-block;width:9px;height:9px;border-radius:50%;background:#e67e22;vertical-align:middle;margin-right:3px"></span>C</span>'
    + '<span style="font-size:9px;color:#555;margin-left:8px">亮=已配对 · 淡=未配对 · 弧=预测碱基对</span></div></div>';
}

/* —— 位点位置图（移植 PAmiRDB renderPosMap；genome_length 由任务写入每行）—— */
function mirPosMap(r) {
  const tLen = r.genome_length || 0;
  const tStart = r.target_start || 0;
  const tEnd = r.target_end || 0;
  if (!tLen) return '<p style="color:#8892a6;padding:20px;text-align:center">位置数据不可用</p>';
  const W = 580, PAD = 44, H = 90, genomeY = 42, innerW = W - PAD * 2;
  const x1 = Math.max(PAD + Math.round(Math.max(tStart - 1, 0) / tLen * innerW), PAD);
  let x2 = PAD + Math.round(tEnd / tLen * innerW);
  x2 = Math.min(x2, PAD + innerW);
  const bw = Math.max(x2 - x1, 10);
  const svg = '<svg width="' + W + '" height="' + H + '" xmlns="http://www.w3.org/2000/svg" style="max-width:100%;font-family:sans-serif">'
    + '<rect x="' + PAD + '" y="' + (genomeY + 1) + '" width="' + innerW + '" height="4" fill="#d0d0d0" rx="2"/>'
    + '<text x="' + (PAD - 4) + '" y="' + (genomeY + 5) + '" text-anchor="end" font-size="9" fill="#888">genome</text>'
    + '<rect x="' + x1 + '" y="' + (genomeY - 8) + '" width="' + bw + '" height="20" fill="rgba(44,110,73,0.3)" stroke="#2c6e49" stroke-width="1.5" rx="2"/>'
    + '<text x="' + (x1 + bw / 2) + '" y="' + (genomeY - 11) + '" text-anchor="middle" font-size="9" fill="#2c6e49" font-weight="bold">▼</text>'
    + '<line x1="' + (x1 + bw / 2) + '" y1="' + (genomeY + 5) + '" x2="' + (x1 + bw / 2) + '" y2="' + (genomeY + 16) + '" stroke="#2c6e49" stroke-width="1" stroke-dasharray="2,2"/>'
    + '<text x="' + (x1 + bw / 2) + '" y="' + (genomeY + 26) + '" text-anchor="middle" font-size="8.5" fill="#2c6e49" font-weight="600">' + (tStart + 1) + '–' + tEnd + ' nt</text>'
    + '<line x1="' + PAD + '" y1="' + (genomeY + 5) + '" x2="' + PAD + '" y2="' + (genomeY + 9) + '" stroke="#ccc" stroke-width="1"/>'
    + '<line x1="' + (PAD + innerW / 2) + '" y1="' + (genomeY + 5) + '" x2="' + (PAD + innerW / 2) + '" y2="' + (genomeY + 9) + '" stroke="#ccc" stroke-width="1"/>'
    + '<line x1="' + (PAD + innerW) + '" y1="' + (genomeY + 5) + '" x2="' + (PAD + innerW) + '" y2="' + (genomeY + 9) + '" stroke="#ccc" stroke-width="1"/>'
    + '<text x="' + PAD + '" y="' + (H - 10) + '" font-size="8" fill="#aaa">1</text>'
    + '<text x="' + (PAD + innerW / 2) + '" y="' + (H - 10) + '" text-anchor="middle" font-size="8" fill="#aaa">' + Math.round(tLen / 2) + '</text>'
    + '<text x="' + (PAD + innerW) + '" y="' + (H - 10) + '" text-anchor="end" font-size="8" fill="#aaa">' + tLen + ' bp</text>'
    + '</svg>';
  return '<div class="mir-posmap"><h4>结合位点在病毒基因组上的位置</h4>' + svg + '</div>';
}

function mirExportCSV() {
  const rows = [['miRNA_ID', 'Virus_ID', 'Virus_name', 'Consensus', 'Energy', 'Identity%',
    'Mismatches', 'Seed_WC', 'Algos_pass', 'Target_start', 'Target_end', 'miRanda_score']];
  mirS.filtered.forEach(r => {
    rows.push([r.mirna_id, r.target_virus_id || '', r.target_virus_name || '',
      r.consensus_level || 'NONE', r.best_energy ?? '', r.best_identity ?? '',
      r.num_mismatches ?? '', r.seed_perfect_wc ?? '', r.num_algorithms ?? '',
      r.target_start ?? '', r.target_end ?? '', r.miranda_score ?? '']);
  });
  const csv = '\uFEFF' + rows.map(r => r.join(',')).join('\n');
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' }));
  a.download = 'mirna_target_prediction.csv';
  a.click();
}

function mirExportJSON() {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([JSON.stringify(mirS.filtered, null, 2)], { type: 'application/json' }));
  a.download = 'mirna_target_prediction.json';
  a.click();
}

/* ── Query2 过滤（PAmiRDB 口径）：miRNA 第 10、11 位与靶标 A-U / U-A
   Watson-Crick 配对（slicing 切割位点约束）。解析行内 alignment（miRanda
   三行格式），miRNA 链跳过 gap 计位，与后端 compute_seed_quality 同口径。 */
function mirIsQuery2(r) {
  const lines = (r.alignment || '').split('\n');
  if (lines.length < 3) return false;
  // 2026-10-03 适配新帧：上行改为 miRNA 3'→5'（5' 端在右），位次需自右向左数；
  // 配对判定改用碱基互补（本模块比对串的 '|' 同时覆盖 WC 与 G:U wobble，
  // 只看杠线会把 wobble 误判为 WC）。兼容旧格式（上行 5'→3'）。
  const mLabel = lines[0].match(/^miRNA\s+([53])'/);
  const mir3first = mLabel ? mLabel[1] === '3' : false;
  const mirnaLine = lines[0].replace(/^miRNA\s+[53]'\s*/, '').replace(/\s*[53]'$/, '');
  const pipeLine = lines[1].trim();
  const targetLine = lines[2].replace(/^target\s+[53]'\s*/, '').replace(/\s*[53]'$/, '');
  const nMir = mirnaLine.replace(/-/g, '').length;
  const ok = { 10: false, 11: false };
  let pos = 0;
  for (let i = 0; i < mirnaLine.length && i < pipeLine.length; i++) {
    if (mirnaLine[i] === '-') continue;
    pos++;
    const pos5 = mir3first ? (nMir - pos + 1) : pos;
    if (pos5 > 11) continue;
    if (pos5 === 10 || pos5 === 11) {
      const b = mirnaLine[i].toUpperCase();
      const t = (targetLine[i] || '').toUpperCase();
      const partner = { 'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G' }[b];
      ok[pos5] = (b === 'A' || b === 'U') && !!partner && t === partner;
    }
  }
  return ok[10] && ok[11];
}

/* ── 候选 miRNA 库选择器（数据与验证见 databases/mirna_lib/README.md）── */
const mirLibS = { loaded: false, sel: new Map(), lastItems: [], debT: null };

function mirOpenLib() {
  $('mirLibMask').style.display = 'flex';
  if (!mirLibS.loaded) {
    mirLibS.loaded = true;
    fetch('/api/tool/mirna_lib/facets').then(r => {
      if (!r.ok) throw new Error('HTTP ' + r.status);
      return r.json();
    }).then(f => {
      const fill = (id, anyKey, pairs) => {
        $(id).innerHTML = '<option value="">' + t(anyKey) + '</option>'
          + pairs.map(([name, n]) => '<option value="' + esc(name) + '">'
            + esc(name) + ' (' + n + ')</option>').join('');
      };
      fill('mir-lib-species', 'tk.mirLibAnySpecies', f.species);
      fill('mir-lib-source', 'tk.mirLibAnySrc', f.sources);
      mirLibSearch();
    }).catch(e => {
      mirLibS.loaded = false;
      $('mirLibMask').style.display = 'none';
      alert(t('ui.mi.candLibFail') + ': ' + e);
    });
  }
}

function mirCloseLib() { $('mirLibMask').style.display = 'none'; }

function mirLibSearchDebounced() {
  clearTimeout(mirLibS.debT);
  mirLibS.debT = setTimeout(mirLibSearch, 250);
}

async function mirLibSearch() {
  const url = '/api/tool/mirna_lib?q=' + encodeURIComponent($('mir-lib-q').value.trim())
    + '&species=' + encodeURIComponent($('mir-lib-species').value)
    + '&source=' + encodeURIComponent($('mir-lib-source').value) + '&limit=100';
  try {
    const res = await fetch(url).then(r => r.json());
    mirLibS.lastItems = res.items || [];
    $('mir-lib-count').textContent = '命中 ' + res.total + ' 条'
      + (res.total > mirLibS.lastItems.length ? '，显示前 ' + mirLibS.lastItems.length + ' 条' : '');
    const tb = $('mir-lib-tbody');
    if (!mirLibS.lastItems.length) {
      tb.innerHTML = '<tr><td colspan="5" style="text-align:center;padding:14px;color:#8892a6">'
        + t('tk.mirLibEmpty') + '</td></tr>';
      return;
    }
    tb.innerHTML = mirLibS.lastItems.map((r, i) => {
      const on = mirLibS.sel.has(r.mirna_id + '|' + r.species);
      const note = r.source === 'plantrg'
        ? (r.target_count ? 'targets ' + esc(r.target_count) : esc(r.source))
        : esc(r.source);
      return '<tr style="border-bottom:1px solid #eef1f5;' + (on ? 'background:#eef7f1' : '') + '">'
        + '<td style="padding:4px 6px"><input type="checkbox" ' + (on ? 'checked' : '')
        + ' onchange="mirLibToggle(' + i + ', this.checked)"></td>'
        + '<td style="padding:4px 6px;font-weight:600;color:#2c6e49">' + esc(r.mirna_id) + '</td>'
        + '<td style="padding:4px 6px">' + esc(r.species) + '</td>'
        + '<td style="padding:4px 6px;font-family:Consolas,monospace">' + esc(r.seq)
        + ' <span style="color:#8892a6">(' + r.len + ')</span></td>'
        + '<td style="padding:4px 6px;font-size:10px;color:#67748e">' + note + '</td></tr>';
    }).join('');
  } catch (e) {
    alert(t('ui.mi.searchFail') + ': ' + e);
  }
}

function mirLibToggle(idx, on) {
  const r = mirLibS.lastItems[idx];
  if (!r) return;
  const k = r.mirna_id + '|' + r.species;
  if (on) mirLibS.sel.set(k, r); else mirLibS.sel.delete(k);
  $('mir-lib-selinfo').textContent = '已选 ' + mirLibS.sel.size + ' 条';
}

function mirLibAdd() {
  if (!mirLibS.sel.size) { alert(t('ui.mi.needCand')); return; }
  const ta = $('mirna_text');
  const lines = [];
  mirLibS.sel.forEach(r => {
    // plantrg 行 mirna_id 是短名（miR156a-5p），拼物种前缀避免跨物种撞名
    const id = r.source === 'plantrg'
      ? r.species.replace(/ /g, '_') + '-' + r.mirna_id
      : r.mirna_id;
    lines.push('>' + id + ' ' + r.species);
    lines.push(r.seq);
  });
  const block = lines.join('\n');
  ta.value = ta.value.trim() ? ta.value.replace(/\s+$/, '') + '\n' + block : block;
  mirCloseLib();
}
