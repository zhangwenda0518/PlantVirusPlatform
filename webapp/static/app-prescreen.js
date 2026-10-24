/* 病毒样本预筛页（/prescreen）—— 数据资源组
 *
 * 依赖 app.js 的 jfetch / val / esc / watchTaskBtn 与 i18n.js 的 t()。
 * 口径红线（不要在前端改写这些措辞）：
 *   - 两条通路分母不同：virus_pct 不可跨来源比，跨来源看 virus_pct_ident
 *   - 这是 k-mer 快速分类，不是严谨比对；只出清单，不干预下载范围
 */
'use strict';

const PS = { data: null, chart: null };

/* ── 输入来源切换 ── */
function psSrcChange() {
  const v = (document.querySelector('input[name="psSrc"]:checked') || {}).value || 'paste';
  [['paste', 'psPanePaste'], ['project', 'psPaneProject'], ['file', 'psPaneFile']]
    .forEach(function (p) {
      const el = document.getElementById(p[1]);
      if (el) el.classList.toggle('on', p[0] === v);
    });
  if (v === 'project') psLoadProjects();
}

async function psLoadProjects() {
  const sel = document.getElementById('psProject');
  if (!sel || sel.dataset.loaded === '1') return;
  try {
    const d = await jfetch('/api/prescreen/projects');
    const ps = d.projects || [];
    sel.innerHTML = ps.length
      ? ps.map(function (p) {
          const tb = p.tables.length ? p.tables.join('/') : t('ps.noTable', '无元数据表');
          const mark = p.has_result ? ' · ' + t('ps.hasResult', '已有预筛结果') : '';
          return '<option value="' + esc(p.name) + '">' + esc(p.name)
            + ' — ' + esc(tb) + esc(mark) + '</option>';
        }).join('')
      : '<option value="">' + t('ps.noProj', '（暂无检索项目，请先去「公共数据检索」跑一次）') + '</option>';
    sel.dataset.loaded = '1';
  } catch (e) {
    sel.innerHTML = '<option value="">' + esc(t('ps.loadFail', '加载失败') + ': ' + e.message) + '</option>';
  }
}

/* ── 启动预筛 ── */
async function psRun(btn) {
  const err = document.getElementById('psErr');
  err.textContent = '';
  const name = val('psName').trim();
  if (!name) { err.textContent = t('ps.needName', '请填写项目名（结果按此归档）'); return; }
  if (!/^[A-Za-z0-9_\-]+$/.test(name)) {
    err.textContent = t('ps.badName', '项目名只能含字母/数字/下划线/连字符');
    return;
  }
  const source = (document.querySelector('input[name="psSrc"]:checked') || {}).value || 'paste';
  const payload = { name: name, source: source };
  if (source === 'paste') {
    payload.runs = val('psRuns');
    if (!payload.runs.trim()) { err.textContent = t('ps.needRuns', '请粘贴至少一条 Run 编号'); return; }
  } else if (source === 'project') {
    payload.project = val('psProject');
    if (!payload.project) { err.textContent = t('ps.needProj', '请选择一个检索项目'); return; }
  } else {
    payload.path = val('psFile');
    if (!payload.path) { err.textContent = t('ps.needFile', '请选择文件'); return; }
  }
  const pct = val('psPct');
  if (pct !== '') payload.virus_pct = Number(pct);
  payload.jobs = Number(val('psJobs')) || 6;
  payload.gsa_jobs = Number(val('psGsaJobs')) || 3;
  payload.no_taxdump = document.getElementById('psNoTaxdump').checked;
  payload.no_gsa = document.getElementById('psNoGsa').checked;
  payload.recheck_nodata = document.getElementById('psRecheck').checked;

  try {
    const d = await jfetch('/api/prescreen/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload)
    });
    psSaveName();
    psAttachTask(d.task);
  } catch (e) { err.textContent = e.message; }
}

/* 项目名持久化：跳页回来后输入框不再为空，「读取已有结果」直接可用 */
function psSaveName() {
  const name = val('psName').trim();
  if (!name) return;
  try { localStorage.setItem('ps.lastName', name); } catch (e) { /* 隐私模式等 */ }
}

/* 看护一个预筛任务：日志跟随 + 按钮锁定 + 终态后自动读结果。
   首次提交与切页回来恢复共用这一条路径。 */
function psAttachTask(tid) {
  psPollLog(tid);
  watchTaskBtn(tid, document.getElementById('psRunBtn'), '⏳ 预筛中…',
    psTaskEnd, { cancelable: true });
}

function psTaskEnd(snap) {
  const log = document.getElementById('psLog');
  const hint = document.getElementById('psHint');
  if (hint) hint.textContent = '';
  if (snap && snap.status === 'done') {
    log.textContent += '\n' + t('ps.done', '✓ 预筛完成') + '\n';
    psRefresh(true);
  } else if (snap) {
    log.textContent += '\n✘ ' + (snap.error || snap.status) + '\n';
  }
}

/* ── 切页回来重挂运行中任务 ──
 * 预筛在 TaskManager 线程里跑，与本页生命周期无关；但日志轮询与按钮看护
 * 都是页面内存态，跳走即失（这是「点了运行去别页、回来找不到任务」的原因）。
 * 页面加载时先问后端本页是否还有任务在跑：
 *   有 → 回填项目名、重挂日志与按钮看护，完成后自动读取结果；
 *   没有 → 回填上次的项目名，已有结果则静默加载（回来即见结果）。
 */
async function psResume() {
  const nameInput = document.getElementById('psName');
  let lastName = '';
  try { lastName = localStorage.getItem('ps.lastName') || ''; } catch (e) { /* ignore */ }
  try {
    const d = await jfetch('/api/prescreen/running');
    const t0 = (d.tasks || [])[0];          // 最新一条（列表新→旧序）
    if (t0) {
      nameInput.value = t0.name;
      psSaveName();
      const hint = document.getElementById('psHint');
      if (hint) hint.textContent = t('ps.resumed',
        '已恢复运行中的预筛任务，完成后自动读取结果；离开本页不影响运行，也可去任务中心查看');
      psAttachTask(t0.task);
      return;
    }
  } catch (e) { /* 恢复失败不影响正常使用 */ }
  if (lastName) {
    nameInput.value = lastName;
    psRefresh(true);
  }
}

function psPollLog(tid) {
  const box = document.getElementById('psLog');
  box.style.display = '';
  box.textContent = '';
  let seen = 0, fails = 0;
  const timer = setInterval(async function () {
    let s = null;
    try {
      const r = await fetch('/api/task/' + tid);
      if (!r.ok) throw new Error('HTTP ' + r.status);
      s = await r.json();
      fails = 0;
    } catch (e) {
      if (++fails >= 15) { clearInterval(timer); box.textContent += '\n' + t('ps.pollLost', '日志跟随中断'); }
      return;
    }
    if (s.log && s.log.length) {
      const fresh = s.log.slice(Math.min(seen, s.log.length));
      seen = s.log.length;
      if (fresh.length) {
        box.textContent += fresh.join('\n') + '\n';
        box.scrollTop = box.scrollHeight;
      }
    }
    if (['done', 'failed', 'cancelled'].includes(s.status)) clearInterval(timer);
  }, 2000);
}

/* ── 读结果 ──
   quiet=true（切页回来自动恢复时）不弹「还没有结果」的提示，只隐藏结果卡 */
async function psRefresh(quiet) {
  const name = val('psName').trim();
  if (!name) { alert(t('ps.needName2', '请先填写项目名')); return; }
  psSaveName();
  const err = document.getElementById('psErr');
  err.textContent = '';
  try {
    const d = await jfetch('/api/prescreen/result?name=' + encodeURIComponent(name));
    if (!d.exists) {
      document.getElementById('psResultCard').style.display = 'none';
      if (!quiet) err.textContent = t('ps.noResult', '该项目还没有预筛结果，请先运行一次');
      return;
    }
    PS.data = d;
    // 换项目要清掉上一个项目的病毒谱，否则会串数据
    if (PSV.name !== d.name) {
      PSV.name = d.name; PSV.data = null; PSV.promise = null; PSV.err = '';
      PSV.level = 'specific'; PSV.search = '';
      const s = document.getElementById('psVirSearch');
      if (s) s.value = '';
      const sel = document.getElementById('psVirLevel');
      if (sel) { sel.innerHTML = ''; sel.dataset.n = ''; }
      const vt = document.getElementById('psVirTable');
      if (vt) vt.innerHTML = '';
      psVirNote('');
    }
    document.getElementById('psResultCard').style.display = '';
    psRenderStats();
    psRenderFiles();
    psRender();
    psChart();
    psVirPrefetch();                 // 静默把病毒谱算出来（首次慢，之后缓存）
    if (document.getElementById('psPaneVirus').style.display !== 'none') psVirEnsure();
  } catch (e) { err.textContent = e.message; }
}

function psRenderStats() {
  const s = PS.data.stats || {};
  const items = [
    [s.total, t('ps.sTotal', '输入 run 数')],
    [s.ok, t('ps.sOk', '有分类数据')],
    [s.with_virus, t('ps.sVirus', '检出病毒 (>0%)')],
    [s.pass, t('ps.sPass', '通过阈值')],
    [s.unresolved, t('ps.sUnres', '无分类数据')],
    ['NCBI ' + s.ncbi_ok + '/' + s.ncbi_total, t('ps.sNcbi', 'NCBI 通路')],
    ['GSA ' + s.gsa_ok + '/' + s.gsa_total, t('ps.sGsa', 'GSA 通路')]
  ];
  document.getElementById('psStats').innerHTML = items.map(function (x) {
    return '<div class="ps-stat"><b>' + esc(String(x[0])) + '</b><span>' + esc(x[1]) + '</span></div>';
  }).join('');
}

function psRenderFiles() {
  const box = document.getElementById('psFiles');
  const name = encodeURIComponent(PS.data.name);
  const files = PS.data.files || [];
  box.innerHTML = files.length ? files.map(function (f) {
    const url = '/prescreen_files/' + name + '/' + encodeURIComponent(f.name);
    return '<a class="btn small" href="' + url + '?dl=1">⬇ ' + esc(f.name)
      + ' <span class="hint">(' + esc(psBytes(f.size)) + ')</span></a>';
  }).join('') : '<span class="hint">' + esc(t('ps.noFiles', '暂无产物')) + '</span>';
}

function psBytes(n) {
  if (!n) return '—';
  if (n >= 1e6) return (n / 1e6).toFixed(1) + ' MB';
  if (n >= 1e3) return (n / 1e3).toFixed(0) + ' KB';
  return n + ' B';
}

/* ── 表格 ── */
function psNum(v, d) {
  if (v === '' || v === null || v === undefined) return '—';
  const n = Number(v);
  return isNaN(n) ? String(v) : n.toFixed(d === undefined ? 4 : d);
}

function psRows() {
  const d = PS.data;
  const v = (document.querySelector('input[name="psView"]:checked') || {}).value || 'all';
  if (v === 'pass') return d.passed || [];
  if (v === 'unres') return d.unresolved || [];
  if (v === 'virus') return (d.summary || []).filter(function (r) {
    return Number(r.virus_pct || 0) > 0;
  });
  return d.summary || [];
}

function psRender() {
  const rows = psRows();
  const el = document.getElementById('psTable');
  if (!rows.length) {
    el.innerHTML = '<p class="hint">' + esc(t('ps.empty', '该视图没有记录')) + '</p>';
    return;
  }
  const head = ['run', t('ps.cSource', '来源'), t('ps.cStatus', '状态'),
    t('ps.cPct', '病毒%'), t('ps.cPctI', '病毒%(已鉴定口径)'),
    t('ps.cTopV', '主要病毒'), t('ps.cReads', '测序量'),
    t('ps.cTopO', '主要非病毒物种'), t('ps.cPass', '通过')];
  let html = '<table class="tbl"><thead><tr>'
    + head.map(function (h) { return '<th>' + esc(h) + '</th>'; }).join('')
    + '</tr></thead><tbody>';
  html += rows.map(function (r) {
    const isV = Number(r.virus_pct || 0) > 0;
    const okst = r.status === 'OK';
    return '<tr style="cursor:pointer" onclick="psDetail(\'' + esc(r.acc) + '\')">'
      + '<td class="mono">' + esc(r.acc) + '</td>'
      + '<td>' + esc(r.source || '') + '</td>'
      + '<td><span class="ps-tag ' + (okst ? 'ok' : 'no') + '">'
      + esc(okst ? 'OK' : String(r.status || '').split(' ')[0]) + '</span></td>'
      + '<td class="num' + (isV ? ' ps-hit' : '') + '">' + esc(psNum(r.virus_pct)) + '</td>'
      + '<td class="num' + (isV ? ' ps-hit' : '') + '">' + esc(psNum(r.virus_pct_ident)) + '</td>'
      + '<td>' + esc(r.top_virus || '—') + '</td>'
      + '<td class="num">' + esc(r.analysed ? Number(r.analysed).toLocaleString() : '—') + '</td>'
      + '<td>' + esc(r.top_organism || '—') + '</td>'
      + '<td>' + (r.pass === 'PASS' ? '<span class="ps-tag ok">PASS</span>' : '') + '</td>'
      + '</tr>';
  }).join('');
  html += '</tbody></table>';
  el.innerHTML = html;
}

async function psDetail(acc) {
  const box = document.getElementById('psDetail');
  box.innerHTML = '<p class="hint">' + esc(t('ps.loadingDetail', '加载明细…')) + '</p>';
  try {
    const d = await jfetch('/api/prescreen/viral?name='
      + encodeURIComponent(PS.data.name) + '&acc=' + encodeURIComponent(acc));
    const rows = d.rows || [];
    if (!rows.length) {
      box.innerHTML = '<p class="hint">' + esc(acc + ' — ' + t('ps.noViral', '没有病毒类 taxon 明细（可能未检出病毒，或该 run 无分类数据）')) + '</p>';
      return;
    }
    box.innerHTML = '<h3 style="font-size:14px;font-weight:500">'
      + esc(acc) + ' — ' + esc(t('ps.detailTitle', '病毒类 taxon 明细'))
      + ' <span class="hint">(' + rows.length + ')</span></h3>'
      + '<div class="ps-tbl" style="max-height:260px"><table class="tbl"><thead><tr>'
      + ['tax_id', 'org', 'path', 'self_count', 'total_count', 'pct']
        .map(function (h) { return '<th>' + esc(h) + '</th>'; }).join('')
      + '</tr></thead><tbody>'
      + rows.map(function (r) {
          return '<tr><td class="mono">' + esc(r.tax_id || '—') + '</td>'
            + '<td>' + esc(r.org || '') + '</td>'
            + '<td class="hint">' + esc(r.path || '—') + '</td>'
            + '<td class="num">' + esc(r.self_count === '' ? '—' : r.self_count) + '</td>'
            + '<td class="num">' + esc(r.total_count || '') + '</td>'
            + '<td class="num">' + esc(psNum(r.pct, 6)) + '</td></tr>';
        }).join('')
      + '</tbody></table></div>';
  } catch (e) {
    box.innerHTML = '<p class="err">' + esc(e.message) + '</p>';
  }
}

/* ── 图：病毒占比 Top 20（用归一化口径，跨来源可比）── */
function psChart() {
  const el = document.getElementById('psChart');
  if (!el || typeof Plotly === 'undefined') return;
  const rows = (PS.data.summary || [])
    .filter(function (r) { return Number(r.virus_pct_ident || 0) > 0; })
    .sort(function (a, b) { return b.virus_pct_ident - a.virus_pct_ident; })
    .slice(0, 20)
    .reverse();
  if (!rows.length) {
    el.innerHTML = '<p class="hint">' + esc(t('ps.noChart', '没有检出病毒的 run，无图可绘')) + '</p>';
    return;
  }
  const dark = document.documentElement.getAttribute('data-theme') === 'dark';
  const fg = dark ? '#e8e8e8' : '#2c2c2a';
  const grid = dark ? 'rgba(255,255,255,.12)' : 'rgba(0,0,0,.08)';
  const colors = rows.map(function (r) {
    return r.source === 'GSA' ? '#BA7517' : '#185FA5';
  });
  Plotly.newPlot(el, [{
    type: 'bar', orientation: 'h',
    x: rows.map(function (r) { return r.virus_pct_ident; }),
    y: rows.map(function (r) { return r.acc; }),
    text: rows.map(function (r) { return r.top_virus || ''; }),
    textposition: 'none',
    marker: { color: colors },
    customdata: rows.map(function (r) { return [r.source, r.top_virus || '—', r.virus_pct]; }),
    hovertemplate: '%{y}<br>病毒%(已鉴定口径): %{x:.4f}%<br>来源: %{customdata[0]}'
      + '<br>主要病毒: %{customdata[1]}<br>本来源口径: %{customdata[2]}%<extra></extra>'
  }], {
    margin: { l: 110, r: 20, t: 36, b: 44 },
    title: { text: t('ps.chartTitle', '病毒占比 Top 20（已鉴定口径，跨来源可比）'), font: { size: 13, color: fg } },
    xaxis: { title: { text: t('ps.chartX', '病毒命中 / 已鉴定命中 %'), font: { size: 11, color: fg } },
             tickfont: { size: 11, color: fg }, gridcolor: grid, zeroline: false },
    yaxis: { tickfont: { size: 11, color: fg }, automargin: true },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    showlegend: false, font: { color: fg }
  }, { displayModeBar: false, responsive: true });
}

/* ══════════════════════════════════════════════════════════════
   病毒谱（病毒视角）
   ──────────────────────────────────────────────────────────────
   把同一批样本翻过来看：检出了哪些病毒、什么分类、占比多少、
   在多少个样本里出现过。
   后端一次算完所有层级（贵的是读分类库，与层级无关），
   故这里的层级切换是纯前端，不再请求。
   ══════════════════════════════════════════════════════════════ */
const PSV = {
  data: null, level: 'specific', promise: null, err: '',
  sort: { col: 'n_samples', dir: -1 },
  search: ''
};

const PSV_RANK_LABEL = {
  specific: '最特异检出', species: '种 species', genus: '属 genus',
  family: '科 family', order: '目 order', class: '纲 class',
  phylum: '门 phylum', kingdom: '界 kingdom', realm: '域 realm'
};

/* 按量级换精度 —— 病毒占比普遍是痕量级，定点 3 位会全抹成 0 */
function psPct(v) {
  const x = Number(v) || 0;
  if (x === 0) return '0';
  if (x >= 1) return x.toFixed(3);
  if (x >= 0.001) return x.toFixed(4);
  return x.toExponential(2);
}

function psMainTab() {
  const el = document.querySelector('input[name=psMain]:checked');
  const v = el ? el.value : 'sample';
  document.getElementById('psPaneSample').style.display = v === 'sample' ? '' : 'none';
  document.getElementById('psPaneVirus').style.display = v === 'virus' ? '' : 'none';
  if (v === 'virus') psVirEnsure();
}

/* 共享同一个请求 Promise —— 预取与用户点击不会重复触发，也不会互相踩 */
function psVirLoad() {
  if (PSV.data) return Promise.resolve(PSV.data);
  if (PSV.promise) return PSV.promise;
  const name = PS.data && PS.data.name;
  if (!name) return Promise.resolve(null);
  PSV.promise = jfetch('/api/prescreen/virome?name=' + encodeURIComponent(name))
    .then(function (d) {
      PSV.promise = null;
      if (d && d.specific) { PSV.data = d; PSV.err = ''; }
      return PSV.data;
    })
    .catch(function (e) {
      PSV.promise = null; PSV.err = e.message; return null;
    });
  return PSV.promise;
}

/* 静默预取：结果一加载就把病毒谱算出来（首次要读分类库，十几秒），
   等用户点开「病毒谱」通常已就绪 */
function psVirPrefetch() { psVirLoad(); }

function psVirNote(msg) {
  const box = document.getElementById('psVirStats');
  if (box) box.innerHTML = msg ? '<div class="ps-loading">' + esc(msg) + '</div>' : '';
}

async function psVirEnsure() {
  if (PSV.data) { psVirRender(); return; }
  psVirNote(t('ps.virLoading', '正在计算病毒谱… 首次需要读取平台分类库（约 15~60 秒），之后会缓存，再次打开秒开。'));
  const d = await psVirLoad();
  if (!d) {
    psVirNote('');
    const e = document.getElementById('psVirErr');
    if (e) e.textContent = (PSV.err || t('ps.virFail', '病毒谱计算失败'));
    return;
  }
  psVirRender();
}

function psVirLevels() {
  const base = PSV.data && (PSV.data.specific || PSV.data[Object.keys(PSV.data)[0]]);
  return (base && base.levels_available) || ['specific'];
}

/* 层级切换是纯前端的：后端一次算完所有层级，切换不重新请求 */
function psVirLevelChange(v) {
  PSV.level = v;
  psVirRender();
}

function psVirRender() {
  if (!PSV.data) return;
  const levels = psVirLevels();
  if (levels.indexOf(PSV.level) < 0) PSV.level = 'specific';
  const sel = document.getElementById('psVirLevel');
  if (sel && sel.dataset.n !== String(levels.length)) {
    sel.innerHTML = levels.map(function (lv) {
      return '<option value="' + esc(lv) + '">' + esc(PSV_RANK_LABEL[lv] || lv) + '</option>';
    }).join('');
    sel.dataset.n = String(levels.length);
  }
  if (sel) sel.value = PSV.level;

  const lv = PSV.data[PSV.level];
  if (!lv) return;
  const err = document.getElementById('psVirErr');
  if (err) err.textContent = '';
  psVirStats(lv);
  psVirTable(lv);
  psVirChart(lv);
  psVirHeat(lv);
  psVirFam(lv);
}

function psVirStats(lv) {
  const rich = (lv.richness || []).reduce(function (m, r) {
    return Math.max(m, Number(r.n) || 0);
  }, 0);
  const items = [
    [lv.n_total, t('ps.vTotal', '样本总数（分母）')],
    [lv.n_with_data, t('ps.vData', '有分类数据')],
    [lv.n_with_virus, t('ps.vHit', '检出病毒样本')],
    [lv.n_distinct, t('ps.vDistinct', '去重病毒数')],
    [rich, t('ps.vRich', '单样本最多病毒数')]
  ];
  document.getElementById('psVirStats').innerHTML = items.map(function (x) {
    return '<div class="ps-stat"><b>' + esc(String(x[0])) + '</b><span>' + esc(x[1]) + '</span></div>';
  }).join('');
  const d = document.getElementById('psVirDenom');
  if (d) {
    d.textContent = t('ps.vDenomTip', '发生率分母 = 全部样本 {n}（含 STAT/GSA 本来没做分类的 {m} 个，「无数据」按未检出计）')
      .replace('{n}', lv.n_total).replace('{m}', lv.n_total - lv.n_with_data);
  }
}

/* ── 清单表 ── */
function psVirRows(lv) {
  const q = PSV.search.trim().toLowerCase();
  let rows = lv.viruses || [];
  if (q) {
    rows = rows.filter(function (v) {
      return (v.name || '').toLowerCase().indexOf(q) >= 0
        || (v.family || '').toLowerCase().indexOf(q) >= 0
        || (v.genus || '').toLowerCase().indexOf(q) >= 0
        || (v.realm || '').toLowerCase().indexOf(q) >= 0;
    });
  }
  const c = PSV.sort.col, dir = PSV.sort.dir;
  return rows.slice().sort(function (a, b) {
    let x = a[c], y = b[c];
    if (typeof x === 'string' || typeof y === 'string') {
      x = String(x || ''); y = String(y || '');
      return dir * x.localeCompare(y);
    }
    return dir * ((Number(x) || 0) - (Number(y) || 0));
  });
}

function psVirSearchInput() {
  const el = document.getElementById('psVirSearch');
  PSV.search = el ? el.value : '';
  const lv = PSV.data && PSV.data[PSV.level];
  if (lv) psVirTable(lv);          // 只重画表，不动图
}

function psVirSortBy(col) {  if (PSV.sort.col === col) PSV.sort.dir = -PSV.sort.dir;
  else PSV.sort = { col: col, dir: (col === 'name' || col === 'family') ? 1 : -1 };
  const lv = PSV.data && PSV.data[PSV.level];
  if (lv) psVirTable(lv);
}

function psVirTable(lv) {
  const rows = psVirRows(lv);
  const box = document.getElementById('psVirTable');
  const cnt = document.getElementById('psVirCount');
  if (cnt) cnt.textContent = t('ps.vShown', '显示 {a} / {b} 条')
    .replace('{a}', rows.length).replace('{b}', (lv.viruses || []).length);
  if (!rows.length) {
    box.innerHTML = '<p class="hint">' + esc(t('ps.vEmpty', '没有匹配的病毒')) + '</p>';
    return;
  }
  const cols = [
    ['name', t('ps.cName', '病毒'), ''],
    ['rank', t('ps.cRank', '层级'), ''],
    ['family', t('ps.cFamily', '科'), ''],
    ['genus', t('ps.cGenus', '属'), ''],
    ['n_samples', t('ps.cN', '检出样本'), 'num'],
    ['occ', t('ps.cOcc', '发生率 %'), 'num'],
    ['pct_mean', t('ps.cMean', '均值 %*'), 'num'],
    ['pct_max', t('ps.cMax', '最大 %'), 'num']
  ];
  const arrow = function (c) {
    return PSV.sort.col === c ? (PSV.sort.dir > 0 ? ' ▲' : ' ▼') : '';
  };
  let h = '<table class="tbl"><thead><tr>' + cols.map(function (c) {
    return '<th class="sortable" style="cursor:pointer" onclick="psVirSortBy(\'' + c[0] + '\')">'
      + esc(c[1]) + arrow(c[0]) + '</th>';
  }).join('') + '</tr></thead><tbody>';
  h += rows.map(function (v) {
    const src = Object.keys(v.sources || {}).join('/');
    return '<tr>'
      + '<td title="' + esc(src) + '">' + esc(v.name) + '</td>'
      + '<td class="ps-rank">' + esc(v.rank) + '</td>'
      + '<td class="ps-fam">' + esc(v.family || '—') + '</td>'
      + '<td class="ps-fam">' + esc(v.genus || '—') + '</td>'
      + '<td class="num">' + v.n_samples + '</td>'
      + '<td class="num">' + (Number(v.occ) || 0).toFixed(2) + '</td>'
      + '<td class="num">' + psPct(v.pct_mean) + '</td>'
      + '<td class="num">' + psPct(v.pct_max) + '</td>'
      + '</tr>';
  }).join('');
  box.innerHTML = h + '</tbody></table>'
    + '<p class="hint" style="margin:6px 0 0">'
    + esc(t('ps.cMeanTip', '* 均值 = 在该病毒**检出到的样本**内求平均；占比由 total_count / analysed 精确重算（产物里的 pct 只有 2 位小数，痕量级会被抹成 0）。'))
    + '</p>';
}

function psVirCsv() {
  const lv = PSV.data && PSV.data[PSV.level];
  if (!lv) return;
  const rows = psVirRows(lv);
  const head = ['name', 'rank', 'taxid', 'realm', 'family', 'genus',
    'n_samples', 'occ_pct', 'pct_mean_detected', 'pct_max',
    'pct_mean_all_samples', 'sources'];
  const esc2 = function (s) {
    s = (s === null || s === undefined) ? '' : String(s);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  };
  const lines = [head.join(',')];
  rows.forEach(function (v) {
    lines.push([v.name, v.rank, v.taxid, v.realm, v.family, v.genus,
      v.n_samples, v.occ, v.pct_mean, v.pct_max, v.pct_mean_all,
      Object.keys(v.sources || {}).join('|')].map(esc2).join(','));
  });
  const blob = new Blob(['\ufeff' + lines.join('\r\n')], { type: 'text/csv;charset=utf-8' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = (PS.data.name || 'prescreen') + '_virome_' + PSV.level + '.csv';
  document.body.appendChild(a); a.click();
  setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}

/* ── 图 ── */
function psPlotTheme() {
  const dark = document.documentElement.getAttribute('data-theme') === 'dark';
  return {
    fg: dark ? '#e8e8e8' : '#2c2c2a',
    grid: dark ? 'rgba(255,255,255,.12)' : 'rgba(0,0,0,.08)'
  };
}

const PSV_RANK_COLOR = {
  species: '#185FA5', genus: '#1F7A5C', family: '#BA7517', order: '#7B4FA8',
  class: '#B34A6B', phylum: '#3E7CA6', kingdom: '#5B8C3A', realm: '#8A6D3B'
};

function psVirChart(lv) {
  const el = document.getElementById('psVirChart');
  if (!el || typeof Plotly === 'undefined') return;
  const rows = (lv.viruses || []).slice(0, 25).slice().reverse();
  if (!rows.length) {
    el.innerHTML = '<p class="hint">' + esc(t('ps.noChart', '没有检出病毒的 run，无图可绘')) + '</p>';
    return;
  }
  const th = psPlotTheme();
  Plotly.newPlot(el, [{
    type: 'bar', orientation: 'h',
    x: rows.map(function (v) { return v.occ; }),
    y: rows.map(function (v) { return v.name; }),
    marker: {
      color: rows.map(function (v) { return PSV_RANK_COLOR[v.rank] || '#8A8A8A'; })
    },
    customdata: rows.map(function (v) {
      return [v.rank, v.n_samples, psPct(v.pct_mean), psPct(v.pct_max),
        v.family || '—', Object.keys(v.sources || {}).join('/')];
    }),
    hovertemplate: '%{y}<br>发生率: %{x:.2f}% (%{customdata[1]}/'
      + lv.n_total + ' 样本)<br>层级: %{customdata[0]}'
      + '<br>检出内均值: %{customdata[2]}% · 最大: %{customdata[3]}%'
      + '<br>科: %{customdata[4]}<br>来源: %{customdata[5]}<extra></extra>'
  }], {
    margin: { l: 240, r: 30, t: 36, b: 44 },
    title: {
      text: t('ps.virChartTitle', '检出率 Top 25（发生率 = 出现该病毒的样本数 / 全部样本数）'),
      font: { size: 13, color: th.fg }
    },
    xaxis: {
      title: { text: t('ps.virChartX', '发生率 %'), font: { size: 11, color: th.fg } },
      tickfont: { size: 11, color: th.fg }, gridcolor: th.grid, zeroline: false
    },
    yaxis: { tickfont: { size: 11, color: th.fg }, automargin: true },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    showlegend: false, font: { color: th.fg }
  }, { displayModeBar: false, responsive: true });
}

function psVirHeat(lv) {
  const el = document.getElementById('psVirHeat');
  if (!el || typeof Plotly === 'undefined') return;
  const m = lv.matrix || {};
  if (!m.viruses || !m.viruses.length) {
    el.innerHTML = '<p class="hint">' + esc(t('ps.virHeatNone', '没有可画的检出矩阵')) + '</p>';
    return;
  }
  const th = psPlotTheme();
  // log10 变换：占比跨 1e-6 ~ 1 六个数量级，线性色阶会全黑。
  // 真实占比放在 hover 里，不隐藏原始量级。
  const z = m.z.map(function (row) {
    return row.map(function (v) { return v > 0 ? Math.log10(v) : null; });
  });
  const text = m.z.map(function (row) {
    return row.map(function (v) { return v > 0 ? psPct(v) + '%' : '未检出'; });
  });
  Plotly.newPlot(el, [{
    type: 'heatmap', z: z, x: m.samples, y: m.viruses, text: text,
    hovertemplate: '样本 %{x}<br>%{y}<br>占比 %{text}<extra></extra>',
    colorscale: 'Blues', zmin: -6, zmax: 0,
    colorbar: {
      title: { text: t('ps.virHeatCb', 'log10(占比%)'), font: { size: 11, color: th.fg } },
      tickfont: { size: 10, color: th.fg },
      tickvals: [-6, -5, -4, -3, -2, -1, 0],
      ticktext: ['1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1']
    }
  }], {
    margin: { l: 240, r: 10, t: 20, b: 60 },
    xaxis: {
      tickfont: { size: 9, color: th.fg }, showticklabels: false,
      title: { text: t('ps.virHeatX', '样本（按检出丰富度降序）'), font: { size: 11, color: th.fg } }
    },
    yaxis: { tickfont: { size: 11, color: th.fg }, automargin: true },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: th.fg }
  }, { displayModeBar: false, responsive: true });
}

function psVirFam(lv) {
  const el = document.getElementById('psVirFam');
  if (!el || typeof Plotly === 'undefined') return;
  const rows = (lv.by_family || []).slice(0, 15).slice().reverse();
  if (!rows.length) {
    el.innerHTML = '<p class="hint">' + esc(t('ps.virFamNone', '没有可归类的病毒')) + '</p>';
    return;
  }
  const th = psPlotTheme();
  Plotly.newPlot(el, [{
    type: 'bar', orientation: 'h',
    x: rows.map(function (f) { return f.n_samples; }),
    y: rows.map(function (f) { return f.family; }),
    marker: { color: '#185FA5' },
    customdata: rows.map(function (f) { return [f.n_hits, f.viruses]; }),
    hovertemplate: '%{y}<br>涉及样本: %{x}<br>样本×病毒对: %{customdata[0]}'
      + '<br>该科下病毒数: %{customdata[1]}<extra></extra>'
  }], {
    margin: { l: 200, r: 30, t: 36, b: 44 },
    title: {
      text: t('ps.virFamTitle2', '各科涉及样本数 Top 15（横轴为去重样本数）'),
      font: { size: 13, color: th.fg }
    },
    xaxis: { tickfont: { size: 11, color: th.fg }, gridcolor: th.grid, zeroline: false },
    yaxis: { tickfont: { size: 11, color: th.fg }, automargin: true },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    showlegend: false, font: { color: th.fg }
  }, { displayModeBar: false, responsive: true });
}

/* ── 初始化 ── */
(function () {
  psSrcChange();
  psLoadProjects();
  psResume();          // 切页回来：重挂运行中任务 / 回填上次项目名并读结果
})();
