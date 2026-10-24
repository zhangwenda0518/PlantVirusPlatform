/* 植物病毒分析平台 GUI 逻辑（原生 JS，无依赖；多语言文案来自 i18n.js） */
'use strict';

const $ = (id) => document.getElementById(id);

function esc(s) {
  return String(s ?? '').replace(/[&<>\"']/g,
    c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

/* 结果图整理：内置 svguest 的 ?open= 只接受同源 URL（见 pages.py /svguest/）。
   图形在浏览器内编辑/导出，不上传；href 用 encodeURIComponent 整体编码。 */
function svgEditHref(href) {
  return '/svguest/?open=' + encodeURIComponent(href);
}
function svgEditLink(href, label, cls) {
  return `<a class="${cls || 'btn small'}" target="_blank" rel="noopener"`
    + ` href="${esc(svgEditHref(href))}"`
    + ` title="在本地 svguest 中整理（编辑与导出都在浏览器完成，不上传）">`
    + `🎨 ${esc(label || '整理')}</a>`;
}

// 输入框取值（trim）；meta/download/submit 等页的历史页面依赖此全局函数
function val(id) {
  const el = $(id);
  return el && el.value != null ? String(el.value).trim() : '';
}

// 样品创建统一入口（P2-12）：pipeline 页与 samples 页共用，避免重复 fetch 样板。
// payload 字段与 /api/pipeline/create 契约一致：sample/r1/r2/project/subsample
async function apiCreateSample(payload) {
  try {
    const r = await fetch('/api/pipeline/create', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(payload)});
    if (!r.ok) {
      let msg = t('c.fail', '创建失败');
      try { msg = (await r.json()).error || msg; } catch (e) {}
      return { ok: false, error: msg };
    }
    return { ok: true, data: await r.json() };
  } catch (e) {
    return { ok: false, error: t('c.connFail', '无法连接平台服务') + ': ' + e, conn: true };
  }
}

// ---------------- 统一 API fetch（全站共用） ----------------
/* 非 2xx 抛 Error(d.error || 'HTTP xxx')；响应非 JSON 按 {} 处理。
   meta/submit 页的 jfetch、tasks 页的 jf 均委托到这里，行为口径一致。 */
async function jfetch(url, opts) {
  const r = await fetch(url, opts);
  const d = await r.json().catch(() => ({}));
  if (!r.ok) throw new Error(d.error || ('HTTP ' + r.status));
  return d;
}

// ---------------- 按钮防重复点击 / 任务看护 ----------------
// lockBtn(btn, '⏳ 检索中…') → 返回解锁函数；按钮置灰防二次点击。
function lockBtn(btn, busyText) {
  if (!btn) return () => {};
  if (btn.dataset.origHtml == null) btn.dataset.origHtml = btn.innerHTML;
  btn.disabled = true;
  btn.classList.add('btn-busy');
  btn.innerHTML = busyText || '⏳ 运行中…';
  let done = false;
  return function unlock() {
    if (done) return;
    done = true;
    btn.disabled = false;
    btn.classList.remove('btn-busy');
    btn.innerHTML = btn.dataset.origHtml;
  };
}

// 瞬时请求包装：请求期间锁定按钮（fetch 顺序操作也适用）
async function withBtn(btn, fn, busyText) {
  if (btn && btn.disabled) return undefined;   // 已在执行，忽略本次点击
  const unlock = lockBtn(btn, busyText);
  try {
    return await fn(unlock);
  } finally {
    unlock();
  }
}

// 任务看护：锁定按钮直到后台任务进入终态（done/failed/cancelled）。
// 返回终态任务快照；onEnd(t) 在解锁前回调（刷新列表等）。
// opts.cancelable: 锁定期间在按钮旁显示「⏹ 停止」，终态后自动移除。
async function watchTaskBtn(tid, btn, busyText, onEnd, opts) {
  const unlock = lockBtn(btn, busyText);
  let stopBtn = null;
  if (opts && opts.cancelable) {
    stopBtn = document.createElement('button');
    stopBtn.className = 'btn small danger';
    stopBtn.textContent = '⏹ 停止';
    stopBtn.onclick = async () => {
      if (!confirm(t('tk.confirmStop', '停止该任务？'))) return;
      stopBtn.disabled = true; stopBtn.textContent = '停止中…';
      try { await fetch('/api/task/' + tid + '/cancel', {method: 'POST'}); } catch (e) {}
      stopBtn.disabled = false; stopBtn.textContent = '⏹ 停止';
    };
    btn.insertAdjacentElement('afterend', stopBtn);
  }
  const finish = (snap) => {
    if (stopBtn) stopBtn.remove();
    unlock();
    if (onEnd) onEnd(snap);
    return snap;
  };
  let misses = 0;              // 任务快照不可达计数（404/网络）
  for (;;) {
    await new Promise(r => setTimeout(r, 3000));
    let snap = null, httpOk = true;
    try {
      const r = await fetch('/api/task/' + tid);
      httpOk = r.ok;
      snap = await r.json();
    } catch (e) { httpOk = false; }
    if (!httpOk) {
      // 服务重启后内存任务丢失（404）或网络故障：最多重试 10 次（30s）后
      // 解锁按钮，避免永久卡在「运行中」。
      if (++misses >= 10) return finish(null);
      continue;
    }
    misses = 0;
    if (snap && ['done', 'failed', 'cancelled'].includes(snap.status)) {
      return finish(snap);
    }
  }
}

function fmtSize(n) {
  if (n == null) return '';
  if (n > 1e9) return (n / 1e9).toFixed(2) + ' GB';
  if (n > 1e6) return (n / 1e6).toFixed(1) + ' MB';
  if (n > 1e3) return (n / 1e3).toFixed(0) + ' KB';
  return n + ' B';
}

// 秒 → 人读时长
function fmtDur(sec) {
  if (sec == null || sec < 0) return '';
  const s = Math.round(sec);
  if (s >= 3600) return `${Math.floor(s / 3600)} h ${String(Math.floor((s % 3600) / 60)).padStart(2, '0')} min`;
  if (s >= 60) return `${Math.floor(s / 60)} min ${String(s % 60).padStart(2, '0')} s`;
  return `${s} s`;
}

// ---------------- Toast 通知 ----------------
function toast(title, body, opts) {
  let box = $('toastBox');
  if (!box) {
    box = document.createElement('div');
    box.id = 'toastBox';
    document.body.appendChild(box);
  }
  const el = document.createElement('div');
  el.className = 'toast' + ((opts && opts.kind) ? ' ' + opts.kind : '');
  el.innerHTML = `<div class="t-title">${esc(title)}</div>` +
    (body ? `<div class="t-body">${esc(body)}</div>` : '') +
    ((opts && opts.actions)
      ? `<div class="t-actions">${opts.actions.map((a, i) =>
          `<button class="btn small ${a.primary ? 'primary' : ''}" data-ta="${i}">${esc(a.label)}</button>`).join('')}</div>`
      : '');
  if (opts && opts.actions) {
    el.querySelectorAll('[data-ta]').forEach(b => {
      b.onclick = () => { opts.actions[+b.dataset.ta].onClick && opts.actions[+b.dataset.ta].onClick(); el.remove(); };
    });
  }
  box.appendChild(el);
  const ttl = (opts && opts.ttl) != null ? opts.ttl : 8000;
  if (ttl > 0) setTimeout(() => el.remove(), ttl);
  return el;
}

/* 后台完成时发系统通知（页面不可见才打扰；需用户授权） */
function sysNotify(title, body) {
  try {
    if (document.visibilityState === 'visible' || !window.Notification) return;
    if (Notification.permission === 'granted') new Notification(title, {body});
    else if (Notification.permission === 'default') Notification.requestPermission();
  } catch (e) { /* 通知不可用忽略 */ }
}

// ---------------- 服务断连检测 ----------------
let connFails = 0;
function setConnBanner(down) {
  let b = $('conn-banner');
  if (down) {
    if (!b) {
      b = document.createElement('div');
      b.id = 'conn-banner';
      b.style.cssText = 'position:sticky;top:0;z-index:99;background:#c62828;color:#fff;'
        + 'padding:10px 16px;font-size:14px;line-height:1.6;';
      document.body.prepend(b);
    }
    b.innerHTML = '<b>⚠ ' + t('conn.lost', '与平台服务的连接已断开') + '</b><br>'
      + t('conn.lostHint',
          '1. 黑色控制台窗口是否仍在运行（关闭该窗口 = 退出平台）；<br>'
        + '2. 若已关闭：重新双击 VirusPlatform.exe（或 启动平台.bat）；<br>'
        + '3. 大任务占满内存/CPU 时服务可能暂时无响应，任务结束后本提示自动消失。');
    startConnProbe();
  } else if (b) {
    b.remove();
  }
}

// 仅真网络失败（fetch 抛 TypeError）才算断连：服务端 500/JSON 解析失败
// 属局部接口问题，不应谎报"与平台服务的连接已断开"
function netFail(e) { return e instanceof TypeError; }

// 断连横幅自愈探针：每 5s 探测一次，服务恢复即自动摘掉横幅。
// 探测必须带超时：fetch 默认永不超时，第一轮若挂死在半断开的
// keep-alive 连接上，await 永不返回 → 横幅永远无法自愈。
let _connProbe = null;
function startConnProbe() {
  if (_connProbe) return;
  _connProbe = setInterval(async () => {
    const ac = new AbortController();
    const timer = setTimeout(() => ac.abort(), 4000);
    try {
      const r = await fetch('/api/tasks', { signal: ac.signal, cache: 'no-store' });
      if (r.ok) { clearInterval(_connProbe); _connProbe = null; setConnBanner(false); }
    } catch (e) { /* 仍未恢复（或本次探测超时），继续下一轮 */ }
    finally { clearTimeout(timer); }
  }, 5000);
}

// 日志框跟随尾部（运行中任务日志持续追加时自动滚到底）
function autoscrollLogs() {
  document.querySelectorAll('pre[data-autoscroll="1"]').forEach(el => {
    el.scrollTop = el.scrollHeight;
  });
}

// ---------------- 粘贴序列输入 ----------------
/* pasteSeq(inputId, ext)：弹出粘贴对话框，确认后写入 uploads/ 并回填路径。
   ext: '.fasta' | '.fastq' | '.tsv' | '.txt'
   宽屏大文本框（约 1040px 宽 / 22 行高，可拖拽放大、窗口矮时自动收缩）：长序列/多记录 FASTA
   粘贴时不再挤压在窄框里；实时显示字符数与记录数，Ctrl+Enter 确认、Esc 取消。 */
/* 输入框旁的 📋：在**文件输入行下方**展开/收起一个独立的大文本框
   （尺寸对齐 BLAST 卡：rows=10 / min-height 220px），与文件选择彻底分开。
   粘贴内容自动写入临时文件（防抖 700ms + 失焦）并回填路径，各卡原有运行逻辑
   无需改动；实时显示字符数/记录数，Ctrl+Enter 立即写入并收起，Esc 收起，
   再点一次 📋 也收起。 */
var _pasteTimers = {};
/* 与后端 /api/paste_input 的 20MB 上限一致：超限文本不改走粘贴，提示用文件 */
var _PASTE_MAX_CHARS = 20 * 1024 * 1024;
function _pasteWrite(inputId, ext, ta, onDone) {
  var text = (ta.value || '').trim();
  if (!text) return;
  if (text.length > _PASTE_MAX_CHARS) {
    alert(t('ui.app.pasteTooBig'));
    return;
  }
  fetch('/api/paste_input', {
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ text: text, ext: ext })
  }).then(function (r) { return r.json(); }).then(function (d) {
    if (!d.path) { alert(d.error || t('ui.app.writeFail')); return; }
    var inp = $(inputId);
    if (inp) {
      inp.value = d.path;
      /* 通知页面「值已变」：如 /genome 的实时预览需要立刻重绘 */
      inp.dispatchEvent(new Event('input', { bubbles: true }));
      inp.dispatchEvent(new Event('change', { bubbles: true }));
    }
    if (typeof onDone === 'function') onDone(d.path);
  }).catch(function (e) { alert(t('ui.common.connShort') + ': ' + e); });
}
function pasteSeq(inputId, ext) {
  ext = ext || '.fasta';
  var inp = $(inputId);
  if (!inp) return;
  var boxId = 'pastebox-' + inputId;
  var old = document.getElementById(boxId);
  if (old) { old.remove(); return; }               // 再点一次 = 收起
  var isFq = ext === '.fastq' || ext === '.fq';
  var isFa = !isFq && ['.fasta', '.fa', '.fna', '.fas', '.faa'].indexOf(ext) >= 0;
  var fmtHint = isFq ? 'FASTQ（@ 开头）' : (isFa ? 'FASTA（> 开头）' : '文本');
  var box = document.createElement('div');
  box.id = boxId;
  box.style.margin = '6px 0 2px';
  box.innerHTML =
    '<textarea rows="10" spellcheck="false" style="width:100%;min-height:220px;' +
    'font-family:var(--mono);font-size:12.5px;line-height:1.55;padding:10px 12px;' +
    'resize:vertical;border-radius:6px" placeholder="直接 Ctrl+V 粘贴 ' + fmtHint +
    ' 序列…"></textarea>' +
    '<div class="hint" style="margin:4px 0 0;display:flex;justify-content:space-between;gap:12px">' +
    '<span class="paste-stat">0 字符</span>' +
    '<span>粘贴后自动写入临时文件 · Ctrl+Enter 立即写入并收起 · Esc 收起</span></div>';
  var row = inp.closest('.filerow') || inp.parentElement;
  row.insertAdjacentElement('afterend', box);
  var ta = box.querySelector('textarea');
  var stat = box.querySelector('.paste-stat');
  function updStat() {
    var v = ta.value;
    var n = v.replace(/\s/g, '').length;
    var recs = (v.match(/^[>@]/gm) || []).length;
    stat.textContent = n.toLocaleString() + ' 字符' + (recs ? ' · ' + recs + ' 条记录' : '');
  }
  function schedule() {
    updStat();
    clearTimeout(_pasteTimers[inputId]);
    _pasteTimers[inputId] = setTimeout(function () { _pasteWrite(inputId, ext, ta); }, 700);
  }
  ta.addEventListener('input', schedule);
  ta.addEventListener('blur', function () { _pasteWrite(inputId, ext, ta); });
  ta.addEventListener('keydown', function (e) {
    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
      e.preventDefault();
      clearTimeout(_pasteTimers[inputId]);
      _pasteWrite(inputId, ext, ta, function () {
        if (typeof toast === 'function') toast('已写入临时文件', '', { ttl: 2200 });
        box.remove();
      });
    } else if (e.key === 'Escape') { e.preventDefault(); box.remove(); }
  });
  ta.focus();
}

// ---------------- 内置示例数据 ----------------
// 平台自带示例（两条植物病毒完整基因组，~6-7kb），各工具「✨ 示例」按钮共用
const EXAMPLE_FASTA = 'examples/example_viral_contigs.fasta';
// 6 条同属近缘基因组（3 参考 + 3 受控突变衍生株）：MSA / 结构比较 / SDT 等多序列示例
const EXAMPLE_SET_FASTA = 'examples/example_virus_set.fasta';
// 10 条**同种病毒**（黄瓜花叶病毒 CMV）RNA3 完整基因组流系，已 MAFFT 比对（2292 列，gap 4.6%）：
// 专供 RDP5 重组检测（同种比较）。实测 13.7s 报 2 个事件（OQ514051.1 印度株，RDP/GENECONV/
// Maxchi/3Seq 四法支持 p=4.7e-3；PP942736.1 中国株，SiScan p=5.3e-26）。
// ⚠ 别拿它当「序列比对」卡的示例——那张卡是跨病毒比较，两者用途不同。
const EXAMPLE_RECOMB_FASTA = 'examples/example_recomb_set.fasta';
// 5 条同种近缘序列（≥90% 一致）：保守区引物设计示例（需全表保守区段）
const EXAMPLE_CONSERVED_FASTA = 'examples/example_conserved_set.fasta';
// 示例树（上集建树产物）：进化树查看器「✨ 示例」
const EXAMPLE_TREE_NWK = 'examples/example_tree.nwk';
// 示例树 + 基因组叠加：3 条共线性基因组的 NJ 树 + 基因轨道 / 同源连线
// （叠加来自与树同名的 example_synteny_tree.overlay.json，由后端按名自动带上）
const EXAMPLE_SYNTENY_TREE_NWK = 'examples/example_synteny_tree.nwk';
// 示例 GenBank（含 CDS 注释）：基因组图谱 GenBank 模式「✨ 示例」
const EXAMPLE_GENBANK_GB = 'examples/example_genome.gb';
// 共线性比较离线示例（3 条同属小基因组 .gb，逗号分隔供 s_files 导入）
const EXAMPLE_GB_TRIO = 'examples/example_synteny_A.gb,examples/example_synteny_B.gb,examples/example_synteny_C.gb';
// 时间与地理推断·本地全链（t-phylodyn）的三个示例：树 / 比对 / 元数据
const EXAMPLE_DATING_TREE = 'examples/example_treedater.nwk';
const EXAMPLE_PHYLOGEO_FASTA = 'examples/example_phylogeo.fasta';
const EXAMPLE_PHYLOGEO_META = 'examples/example_phylogeo.meta.csv';
// 示例数据的**实际绝对路径**表（服务端给出）：示例目录可能位于程序目录
// 之外（程序/数据库/示例三分离打包），此时下面这些相对路径不成立，
// fillExample 会用 /api/example_paths 换成真实绝对路径。
let _EX_PATHS = null;
async function _exampleAbs(path) {
  const raw = String(path || '');
  if (!raw) return raw;
  try {
    if (!_EX_PATHS) {
      const r = await fetch('/api/example_paths');
      const d = r.ok ? await r.json() : {};
      _EX_PATHS = d.files || {};
    }
  } catch (e) { _EX_PATHS = _EX_PATHS || {}; }
  return raw.split(',').map(function (one) {
    const s = one.trim();
    if (!s) return s;
    const base = s.split(/[\\/]/).pop();
    return (_EX_PATHS && _EX_PATHS[base]) || s;
  }).join(',');
}
async function fillExample(inputId, path) {
  var inp = $(inputId);
  if (!inp) return;
  inp.value = await _exampleAbs(path || EXAMPLE_FASTA);
  if (typeof toast === 'function') {
    toast('已填入内置示例', inp.value, { ttl: 3000 });
  }
}

/* 「参考序列获取」卡「✨ 示例」：预填递进选择（界→科→属）+ 集合名 */
function seqPrepExample() {
  var name = $('s_name');
  if (name) name.value = 'dianthovirus_ictv';
  ['Realm', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'].forEach(c => {
    delete document.body.dataset['ictv_' + c];
  });
  document.body.dataset['ictv_Realm'] = 'Riboviria';
  document.body.dataset['ictv_Family'] = 'Tombusviridae';
  document.body.dataset['ictv_Genus'] = 'Dianthovirus';
  ictvCascadeRefetch().then(function () {
    if (typeof toast === 'function')
      toast('已填入示例：界→科→属 递进选择', 'Dianthovirus（6 条），点「🔍 预览」离线查看', { ttl: 3500 });
  });
}

const ICTV_RANK_MAIN = ['Realm', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'];
let ictvRankMeta = [];      // [{col, zh}]（含 Subfamily 等附加级）

/* 参考库切换（plant = 植物口径 / ictv = 全病毒界，库本身决定宿主范围）
   切换后清空已选深层级并重拉选项——两库可见的分类范围不同，
   旧选择可能在新库下不存在。 */
function ictvDbChanged() {
  ICTV_RANK_MAIN.forEach(c => { delete document.body.dataset['ictv_' + c]; });
  const pv = $('spPreview');
  if (pv) pv.innerHTML = '';
  ictvCascadeRefetch();
}

/* 当前选中的参考库（缺省 plant） */
function ictvDb() {
  return _v('sp_db') || 'plant';
}

/* 级联选参：任意一级变更 → 清空更深等级 → 重新拉取各级选项 */
async function ictvCascadeRefetch() {
  const box = $('spCascade');
  if (!box) return;
  const qs = ICTV_RANK_MAIN.filter(c => _v('sp_c_' + c))
    .map(c => `${c}=${encodeURIComponent(_v('sp_c_' + c))}`).join('&');
  const qsAll = 'db=' + encodeURIComponent(ictvDb()) + (qs ? '&' + qs : '');
  try {
    const r = await fetch('/api/ictv/cascade?' + qsAll);
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    const d = await r.json();
    ictvRankMeta = d.ranks.map(x => ({ col: x.col, zh: x.zh }));
    const deepestSel = ICTV_RANK_MAIN.filter(c => _v('sp_c_' + c)).pop();
    $('spScope').textContent = deepestSel
      ? `当前最深选择：${deepestSel} = ${_v('sp_c_' + deepestSel)}`
      : '（尚未选择，先从「界」开始逐级下钻）';
    box.innerHTML = ICTV_RANK_MAIN.map(col => {
      const meta = d.ranks.find(x => x.col === col);
      const opts = (meta && meta.options) || [];
      return `<div class="gp-item"><label>${esc(col === 'Realm' ? '界 Realm' : meta ? meta.zh + ' ' + col : col)}</label>` +
        `<select id="sp_c_${col}" onchange="ictvCascadeChanged('${col}')">` +
        `<option value="">（全部）</option>` +
        opts.slice(0, 600).map(o => `<option value="${esc(o.name)}">${esc(o.name)}（${o.n}）</option>`).join('') +
        `</select></div>`;
    }).join('');
    // 回填已选值（后端按已选过滤返回，值仍在选项里）
    ICTV_RANK_MAIN.forEach(c => {
      const sel = $('sp_c_' + c);
      const v = document.body.dataset['ictv_' + c] || '';
      if (sel && v && [...sel.options].some(o => o.value === v)) sel.value = v;
    });
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

function ictvCascadeChanged(col) {
  // 该级变更 → 更深等级全部清空（dataset 记录当前选择）
  const idx = ICTV_RANK_MAIN.indexOf(col);
  ICTV_RANK_MAIN.slice(idx).forEach(c => { delete document.body.dataset['ictv_' + c]; });
  const v = _v('sp_c_' + col);
  if (v) document.body.dataset['ictv_' + col] = v;
  ictvCascadeRefetch();
}

function ictvSelectedLevels() {
  const out = {};
  ICTV_RANK_MAIN.forEach(c => {
    const v = _v('sp_c_' + c) || document.body.dataset['ictv_' + c] || '';
    if (v) out[c] = v;
  });
  return out;
}

async function ictvPreview() {
  const levels = ictvSelectedLevels();
  if (!Object.keys(levels).length) { alert(t('ui.app.needRank')); return; }
  const box = $('spPreview');
  box.innerHTML = '<p class="hint">选参中…</p>';
  try {
    const r = await fetch('/api/ictv/preview', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(Object.assign({}, levels, {
        genome: _v('sp_genome') || 'complete',
        db: ictvDb(),
        limit: +_v('sp_limit') || 200,
        per_genus: +(_v('sp_per_genus') || 0),
        per_species: +(_v('sp_per_species') || 0) }))});
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">选参失败: ' + esc((await r.json()).error || '') + '</p>'; return; }
    const d = await r.json();
    box.innerHTML = `<p class="hint"><b>${esc(d.scope)}</b> 命中 <b>${d.total}</b> 条 · 抽样 <b>${esc(d.sampling)}</b>`
      + `（${d.total > d.n_shown ? '按本地优先显示前 ' + d.n_shown : '全部显示'}）</p>`
      + '<table class="table" style="width:100%;border-collapse:collapse">'
      + '<tr><th>Accession</th><th>来源</th><th>完整度</th><th>物种</th><th>基因组</th><th>宿主</th></tr>'
      + d.rows.map(x => `<tr><td>${esc(x.acc)}</td>`
        + `<td>${x.source === 'ncbi' ? '🌐 ncbi（需下载）' : '💾 ' + esc(x.source)}</td>`
        + `<td>${esc(x.coverage || '—')}</td><td>${esc(x.species)}</td>`
        + `<td>${esc(x.genome || '—')}</td><td>${esc(x.host || '—')}</td></tr>`).join('')
      + '</table>';
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

/* 集合名：留空时自动生成。
   前缀取最深已选分类级（属 > 科 > 目 …）的清洗后小写名，
   后缀 YYYYMMDD_HHMM 保证多次下载不撞名。 */
function defaultCollName() {
  let base = '';
  for (let i = ICTV_RANK_MAIN.length - 1; i >= 0; i--) {
    const v = _v('sp_c_' + ICTV_RANK_MAIN[i]);
    if (v) { base = v; break; }
  }
  if (!base) base = (_v('s_term') || '').trim().split(/\s+/)[0] || 'collection';
  base = base.toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '')
    || 'collection';
  const d = new Date();
  const p = n => String(n).padStart(2, '0');
  return base + '_' + d.getFullYear() + p(d.getMonth() + 1) + p(d.getDate())
    + '_' + p(d.getHours()) + p(d.getMinutes());
}

/* 取集合名：用户填了就用，留空自动生成并回填输入框（可见可改） */
function collName() {
  const el = $('s_name');
  let v = el ? String(el.value || '').trim() : '';
  if (!v) {
    v = defaultCollName();
    if (el) el.value = v;
  }
  return v;
}

async function ictvDownload(btn) {
  const levels = ictvSelectedLevels();
  const coll = collName();
  if (!Object.keys(levels).length) { alert(t('ui.app.needRank')); return; }
  try {
    const r = await fetch('/api/ictv/download', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(Object.assign({}, levels, {
        coll, genome: _v('sp_genome') || 'complete',
        db: ictvDb(),
        limit: +_v('sp_limit') || 200,
        per_genus: +(_v('sp_per_genus') || 0),
        per_species: +(_v('sp_per_species') || 0) }))});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => { loadGbCollections(); loadNcbiCollections(); loadTbColls(); },
               btn, '⏳ ICTV 序列下载中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}

/* ================= 序列比对模块：运行 / 彩色查看器 / 拖拽修剪 / 编辑保存 ================= */
let alData = null, alEditing = false;
let alSelA = -1, alSelB = -1;   // 拖拽选中的列范围（0-based 含端点）
let alUndoStack = [];           // 修剪前快照（撤销用，最多 20 步）
let _alDragA = null;            // 拖拽进行中的锚点列

const AL_NT_COLORS = { A: '#2e7d32', T: '#c62828', U: '#c62828', G: '#e65100',
                       C: '#1565c0', N: '#78909c', '-': '#eceff1', '?': '#eceff1' };
const AL_AA_COLORS = (() => {
  const groups = {
    hydrophobic: 'AVLIMFWC', positive: 'KRH', negative: 'DE',
    polar: 'STNQ', special: 'PG', stop: '*' };
  const colors = { hydrophobic: '#90caf9', positive: '#ef9a9a',
                   negative: '#a5d6a7', polar: '#ffe082', special: '#ce93d8',
                   stop: '#b0bec5' };
  const map = { '-': '#eceff1', X: '#cfd8dc', U: '#cfd8dc', B: '#cfd8dc',
                Z: '#cfd8dc', '*': '#b0bec5' };
  for (const g in groups)
    for (const ch of groups[g]) map[ch] = colors[g];
  return map;
})();
function alColor(ch, type) {
  const m = type === 'aa' ? AL_AA_COLORS : AL_NT_COLORS;
  return (m[ch] || '#cfd8dc');
}

async function alignRun(btn) {
  const seqs = _v('al_fa');
  const extra = _v('al_fa_extra');
  if (!seqs && !extra) {
    alert(t('ui.common.needFasta2')); return;
  }
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tool: 'align',
                             params: { seqs: seqs || extra,
                                       seqs_extra: extra || '',
                                       strategy: _v('al_strategy') || 'auto',
                                       trimal: _v('al_trimal') || 'automated1',
                                       max_n: +_v('al_maxn') || 100 } })});
    if (!r.ok) { alert(t('ui.common.startFail') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => {
      // 完成后自动打开清剪结果（无清剪则原比对）
      const res = (TMlastResult(d.task) || {});
      const p = res.aln_used === 'aln.trim.fasta'
        ? `tool_runs/${d.run}/aln.trim.fasta` : `tool_runs/${d.run}/aln.fasta`;
      $('alPath').value = p;
      alignLoad();
    }, btn, '⏳ 比对中…');
  } catch (e) { alert(t('ui.common.connFail') + ': ' + e); }
}
function TMlastResult(tid) {
  // 任务面板轮询的数据不在 JS 侧持久化：退化为直接请求比对文件即可
  return null;
}

async function alignLoad() {
  const p = _v('alPath');
  if (!p) { alert(t('ui.app.needAlnPath')); return; }
  const box = $('alBox');
  box.innerHTML = '<p class="hint">加载中…</p>';
  try {
    const r = await fetch('/api/align/file?path=' + encodeURIComponent(p));
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    alData = await r.json();
    alSelA = alSelB = -1; alUndoStack = [];
    alignRender();
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

function _alWrapW() { return +($('alCols') && $('alCols').value) || 100; }

function _alLineHtml(seq, b) {
  // [b*W, b*W+W) 区间渲染：相邻同色且同选中态的字符合并为一个 span（RLE）。
  // 长 gap 段 / 同色连续段被大幅压缩，万列比对也能整条渲染不卡。
  const W = _alWrapW();
  const start = b * W, end = Math.min(seq.length, start + W);
  const a = Math.min(alSelA, alSelB), z = Math.max(alSelA, alSelB);
  let html = '', i = start;
  while (i < end) {
    const ch = seq[i];
    const sel = alSelA >= 0 && i >= a && i <= z;
    let j = i + 1;
    while (j < end && seq[j] === ch && (alSelA >= 0 && j >= a && j <= z) === sel) j++;
    html += `<span style="color:${alColor(ch, alData.type)};font-weight:700;${sel ? 'background:#fbbf24;color:#111;' : ''}">${esc(seq.slice(i, j))}</span>`;
    i = j;
  }
  return html;
}

function _alBlockHtml(b) {
  // 一个"块" = 全部序列的第 [b*W, b*W+W) 列；整条比对 = 块序列横向滚动。
  const W = _alWrapW(), cols = alData.cols, nameW = alData._nameW;
  const start = b * W, end = Math.min(cols, start + W);
  let h = `<div class="al-block" data-b="${b}" style="margin-bottom:8px">`;
  h += `<div style="color:#94a3b8;white-space:pre">${' '.repeat(nameW)} `
    + `${String(start + 1).padStart(6)} ──── ${String(end).padEnd(6)}</div>`;
  alData.names.forEach((name, i) => {
    const s = (alData.seqs[i] + '-'.repeat(cols)).slice(0, cols);
    h += `<div style="white-space:pre">`
       + `<span style="color:#334155;font-weight:600;display:inline-block;width:${nameW}ch;vertical-align:top;overflow:hidden">${esc(name)}</span> `
       + `<span class="al-seq" data-b="${b}">${_alLineHtml(s, b)}</span></div>`;
  });
  return h + '</div>';
}

function alignRender() {
  const box = $('alBox'), nav = $('alNav');
  if (!alData) { box.innerHTML = '<p class="hint">暂无数据</p>'; nav.style.display = 'none'; return; }
  if (!alData.aligned)
    $('alMeta').textContent = '⚠ 序列长度不一致（未比对或 FASTA）：彩色视图按最长序列展示，建议先运行比对';
  else
    $('alMeta').textContent = '';
  const W = _alWrapW();
  alData._nameW = Math.max(8, ...alData.names.map(n => n.length));
  const nB = Math.max(1, Math.ceil(alData.cols / W));
  let html = '<div style="font-family:Consolas,monospace;font-size:12.5px;user-select:none">';
  for (let b = 0; b < nB; b++) html += _alBlockHtml(b);
  box.innerHTML = html + '</div>';
  _alBindDrag(box);
  nav.style.display = 'flex';
  $('alSaveBtn').style.display = '';
  $('alTrimDel').style.display = alSelA >= 0 ? '' : 'none';
  $('alTrimKeep').style.display = alSelA >= 0 ? '' : 'none';
  $('alUndo').style.display = alUndoStack.length ? '' : 'none';
  _alSelInfo();
  // 图例
  const key = alData.type === 'aa'
    ? ['AVLIMFWC|疏水', 'KRH|正电', 'DE|负电', 'STNQ|极性', 'PG|特殊', '-|gap']
    : ['A', 'T/U', 'G', 'C', 'N', '-'].map(x => x + '|' + ({ A: 'A', 'T/U': 'T/U', G: 'G', C: 'C', N: 'N', '-': 'gap' }[x]));
  $('alLegend').innerHTML = key.map(kv => {
    const [chars, label] = kv.split('|');
    const c = alColor(chars[0], alData.type);
    return `<span style="margin-right:12px;font-size:12px"><span style="display:inline-block;width:12px;height:12px;background:${c};vertical-align:-2px;margin-right:3px;border-radius:2px"></span>${esc(label)}</span>`;
  }).join('');
  const dl = $('alDl');
  dl.style.display = '';
  dl.href = '/tool_runs/' + alData.path.replace(/^.*tool_runs[\/\\]/, '');
  dl.download = alData.path.split(/[\\/]/).pop();
  dl.textContent = t('tk.dlRef', '⬇ 下载') + ' ' + alData.path.split(/[\\/]/).pop();
}

/* ---- 拖拽选列与修剪（✂删选中列 / 🎯只留选中列 / ↩撤销）---- */
function _alColFromEvent(e) {
  const el = e.target && e.target.closest ? e.target.closest('.al-seq') : null;
  if (!el || !alData) return null;
  const b = +el.dataset.b, W = _alWrapW();
  const rect = el.getBoundingClientRect();
  const len = Math.min(W, alData.cols - b * W);
  if (rect.width <= 0 || len <= 0) return null;
  const cw = rect.width / len;   // 等宽字体：整段宽 / 列数 = 列宽
  return b * W + Math.min(len - 1, Math.max(0, Math.floor((e.clientX - rect.left) / cw)));
}

function _alUpdateBlock(b) {
  const el = document.querySelector(`#alBox .al-block[data-b="${b}"]`);
  if (el) el.outerHTML = _alBlockHtml(b);
}

function _alSelInfo() {
  const el = $('alSelInfo');
  if (!el) return;
  if (alSelA < 0 || !alData) { el.textContent = ''; return; }
  const a = Math.min(alSelA, alSelB), b = Math.max(alSelA, alSelB);
  el.textContent = t('tk.alSelRange', '已选列 {a}–{b}（{n} 列）')
    .replace('{a}', a + 1).replace('{b}', b + 1).replace('{n}', String(b - a + 1));
}

function _alBindDrag(box) {
  if (box._alBound) return;   // 事件委托只挂一次（render 会重建内部 DOM）
  box._alBound = true;
  box.addEventListener('mousedown', e => {
    if (!alData || alEditing) return;
    const c = _alColFromEvent(e);
    if (c === null) return;
    e.preventDefault();
    _alDragA = c; alSelA = alSelB = c;
    _alUpdateBlock(+e.target.closest('.al-seq').dataset.b);
    _alSelInfo();
  });
  box.addEventListener('mousemove', e => {
    if (_alDragA === null || !alData) return;
    const c = _alColFromEvent(e);
    if (c === null || c === alSelB) return;
    alSelB = c;
    // 拖拽中只重画锚点/光标两个块；mouseup 后一次性铺满全部块
    _alUpdateBlock(Math.floor(c / _alWrapW()));
    _alUpdateBlock(Math.floor(_alDragA / _alWrapW()));
    _alSelInfo();
  });
  document.addEventListener('mouseup', () => {
    if (_alDragA === null) return;
    _alDragA = null;
    if (alData && alSelA >= 0) alignRender();
  });
}

function _alApply(op) {
  if (!alData || alSelA < 0) return;
  const a = Math.min(alSelA, alSelB), b = Math.max(alSelA, alSelB);
  alUndoStack.push({ names: alData.names.slice(), seqs: alData.seqs.map(s => s.slice()) });
  if (alUndoStack.length > 20) alUndoStack.shift();
  alData.seqs = alData.seqs.map(s => op === 'del'
    ? s.slice(0, a) + s.slice(b + 1)
    : s.slice(a, b + 1));
  alData.cols = alData.seqs.reduce((m, s) => Math.max(m, s.length), 0);
  alData.aligned = alData.seqs.every(s => s.length === alData.cols);
  alSelA = alSelB = -1;
  alignRender();
}

function alTrimDel() { _alApply('del'); }
function alTrimKeep() { _alApply('keep'); }

function alUndo() {
  const st = alUndoStack.pop();
  if (!st || !alData) { return; }
  alData.names = st.names;
  alData.seqs = st.seqs;
  alData.cols = st.seqs.reduce((m, s) => Math.max(m, s.length), 0);
  alData.aligned = st.seqs.every(s => s.length === alData.cols);
  alSelA = alSelB = -1;
  alignRender();
}

function alignEditToggle() {
  alEditing = !alEditing;
  const wrap = $('alEditWrap');
  wrap.style.display = alEditing ? '' : 'none';
  $('alEditBtn').textContent = alEditing ? '👁 退出编辑' : '✏️ 编辑模式';
  if (alEditing && alData) {
    // FASTA 化（70 列换行）
    let out = '';
    alData.names.forEach((n, i) => {
      out += '>' + n + '\n';
      const s = alData.seqs[i];
      for (let j = 0; j < s.length; j += 70) out += s.slice(j, j + 70) + '\n';
    });
    $('alEdit').value = out;
    $('alEditMsg').textContent = '';
  }
}

async function alignSave(btn) {
  // 两个入口：✏️ 编辑模式存 textarea 内容；查看器模式直接把当前
  // （可能已拖拽修剪过的）状态序列化为 FASTA 存副本
  let content = '';
  if (alEditing) {
    content = $('alEdit').value;
    if (!content.trim()) { alert(t('ui.app.editorEmpty')); return; }
  } else {
    if (!alData) { alert(t('ui.app.needAlnOpen')); return; }
    let out = '';
    alData.names.forEach((n, i) => {
      out += '>' + n + '\n';
      const s = alData.seqs[i];
      for (let j = 0; j < s.length; j += 70) out += s.slice(j, j + 70) + '\n';
    });
    content = out;
  }
  const msg = alEditing ? $('alEditMsg') : $('alSaveMsg');
  btn.disabled = true;
  try {
    const r = await fetch('/api/align/save', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ path: alData.path, content }) });
    const d = await r.json();
    if (!r.ok) { msg.textContent = '❌ ' + (d.error || '保存失败'); return; }
    msg.textContent = `✅ 已保存 ${d.n} 条${d.aligned ? '（等长比对）' : '（⚠ 长度不一致——未比对状态）'}：${d.saved}`;
    // 打开保存后的副本继续查看（后续可直接送建树 / SDT）
    $('alPath').value = d.saved;
    alignLoad();
  } catch (e) { msg.textContent = '无法连接: ' + e; }
  btn.disabled = false;
}

function alignSend(target) {
  if (!alData) { alert(t('ui.app.needAlnOpen2')); return; }
  const p = alData.path;
  if (target === 'treebuild') {
    $('qt_fa').value = p;
    location.hash = '#t-treebuild';
  } else {
    $('sd_fa').value = p;
    $('sd_aligned').checked = alData.aligned;
    location.hash = '#t-sdt';
  }
  if (typeof renderModuleTree === 'function') renderModuleTree();
}

/* ================= 进化树构建：集合下拉 + 建树 ================= */
async function loadTbColls() {
  const sel = $('tbColl');
  if (!sel) return;
  try {
    const cols = await (await fetch('/api/gb/collections')).json();
    sel.innerHTML = cols.length
      ? cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}（${c.n_records} 条记录）</option>`).join('')
      : '<option value="">（暂无集合——先到「参考序列获取」下载）</option>';
  } catch (e) { sel.innerHTML = '<option value="">（无法连接）</option>'; }
}

async function tbBuild(btn) {
  const coll = _v('tbColl');
  if (!coll) { alert(t('ui.app.needGbSet')); return; }
  await gbBuildTree(coll, btn);
  loadGbCollections();
}

/* 通用文本复制（引物序列等小段文本） */
async function copyText(text, tag) {
  try {
    await navigator.clipboard.writeText(text);
  } catch (e) {
    const ta = document.createElement('textarea');
    ta.value = text; document.body.appendChild(ta);
    ta.select(); document.execCommand('copy'); ta.remove();
  }
  if (typeof toast === 'function') toast('已复制', tag || '', { ttl: 2500 });
}

/* Metabuli 同款示例病毒（TMV / PVY / CMV / PSTVd / Mix）：
   来源 <DEMO-IP>/metabuli examples，已固化到 examples/ */
const EXAMPLE_METABULI = {
  tmv: 'examples/example_tmv.fasta',
  pvy: 'examples/example_pvy.fasta',
  cmv: 'examples/example_cmv.fasta',
  pstvd: 'examples/example_pstvd.fasta',
  mix: 'examples/example_mix.fasta'
};
function vEx(key, inputId) {
  const path = EXAMPLE_METABULI[key];
  if (!path) return;
  const inp = $(inputId || 'c_fa');
  if (inp) inp.value = path;
  if (key === 'pstvd') {                 // 359nt 类病毒：低于默认 500bp 过滤线
    const ml = $('c_minlen');
    if (ml && +ml.value > 300) { ml.value = 200; }
  }
  const names = { tmv: 'TMV 烟草花叶病毒', pvy: 'PVY 马铃薯Y病毒',
                  cmv: 'CMV 黄瓜花叶病毒(RNA1-3)', pstvd: 'PSTVd 马铃薯纺锤块茎类病毒',
                  mix: 'Mix 四病毒混合(6条)' };
  if (typeof toast === 'function')
    toast('已填入示例病毒', names[key] + ' · ' + path, { ttl: 3500 });
}

// ---------------- 表格分页（每页 20 行） ----------------
const PER_PAGE = 20;
// fnExpr 里可用 {p} 占位目标页码，如 "metaGo" 或 "dlFileGo('abc',{p})"
function pagerHtml(page, totalPages, fnExpr) {
  totalPages = Math.max(1, Math.floor(totalPages || 1));
  if (totalPages <= 1) return '';
  page = Math.min(Math.max(1, page || 1), totalPages);
  const call = p => fnExpr.includes('{p}')
    ? fnExpr.replace('{p}', p) : `${fnExpr}(${p})`;
  const nums = [];
  for (let n = 1; n <= totalPages; n++) {
    if (n === 1 || n === totalPages || Math.abs(n - page) <= 2) nums.push(n);
  }
  let html = `<div class="pager" style="display:flex;gap:4px;align-items:center;flex-wrap:wrap;margin:8px 0">`
    + `<button class="btn small" ${page <= 1 ? 'disabled' : ''} onclick="${call(page - 1)}">◀</button>`;
  let prev = 0;
  for (const n of nums) {
    if (n - prev > 1) html += '<span class="hint">…</span>';
    html += `<button class="btn small${n === page ? ' primary' : ''}" onclick="${call(n)}">${n}</button>`;
    prev = n;
  }
  html += `<button class="btn small" ${page >= totalPages ? 'disabled' : ''} onclick="${call(page + 1)}">▶</button>`
    + `<span class="hint" style="margin-left:4px">第 ${page} / ${totalPages} 页</span></div>`;
  return html;
}

// ---------------- 跨模块文件交接（下载页 → 其它模块预填） ----------------
// 把文件参数暂存到 sessionStorage，再跳到目标模块页；目标页用 takePrefill 取出并预填输入框。
function sendTo(page, params) {
  try {
    sessionStorage.setItem('vp_prefill',
      JSON.stringify(Object.assign({ page: page, ts: Date.now() }, params)));
  } catch (e) {}
  if (page === 'fastp' || page === 'identify') {
    location.href = '/tools?g=sample#t-' + page;
  } else {
    location.href = '/' + page;
  }
}

// 取出本页的预填参数（page 不匹配则不动；取出即销毁，避免二次误填）
function takePrefill(page) {
  try {
    const d = JSON.parse(sessionStorage.getItem('vp_prefill') || 'null');
    if (!d || d.page !== page) return null;
    sessionStorage.removeItem('vp_prefill');
    return d;
  } catch (e) { return null; }
}

// ---------------- 全局上下文（当前样品 / 当前项目） + 全局清理 ----------------
/* 为什么需要：一级导航是普通 <a href>，每次切模块都是整页刷新，页面级 JS 变量
   （curSample / curStages / curBranch …）必然归零。以前"上次在分析哪个样品"
   只写进 localStorage['vp_last_sample'] 却**从来没有读取点**，所以切回来就失忆。
   这里把「当前样品 / 当前项目」提升为跨模块共享的全局上下文：
     - 真源 = localStorage（同源同浏览器，所有模块页可读）
     - 顶部导航常驻两个可点胶囊，任何模块页都能看到/切换
     - 变更通过 document 上的 'vp-ctx' 事件广播，各页按需重挂自己
   并配一个全局清理总闸，把"还在跑的 + 排队的 + 互相交接的"一次归零
   （只清运行态与关联，不删任何结果文件）。 */
const VP_CTX_SAMPLE_KEY = 'vp_ctx_sample';
const VP_CTX_PROJECT_KEY = 'vp_ctx_project';
// 全局清理要清的浏览器态（跨模块"相互关联"）；
// 主题 / 收藏 / 语言、文件浏览器上次目录属于个人偏好，不在清理范围。
const VP_CTX_RESET_KEYS = ['vp_ctx_sample', 'vp_ctx_project', 'vp_last_sample',
                           'vp_prefill', 'vp_runfill', 'vp_logclosed',
                           'vp_anajump', 'vp_verifyjump'];

function lsGet(k) {
  try { return localStorage.getItem(k) || ''; } catch (e) { return ''; }
}
function lsSet(k, v) {
  try { if (v) localStorage.setItem(k, v); else localStorage.removeItem(k); }
  catch (e) { /* 隐私模式等场景忽略 */ }
}

const VP_CTX = {
  sample: '',
  project: '',
  load() {
    this.sample = lsGet(VP_CTX_SAMPLE_KEY);
    this.project = lsGet(VP_CTX_PROJECT_KEY);
    return this;
  },
  setSample(name) {
    name = String(name || '');
    if (name === this.sample) return false;
    this.sample = name;
    lsSet(VP_CTX_SAMPLE_KEY, name);
    this._emit();
    return true;
  },
  setProject(name) {
    name = String(name || '');
    if (name === this.project) return false;
    this.project = name;
    lsSet(VP_CTX_PROJECT_KEY, name);
    this._emit();
    return true;
  },
  clear() {
    if (!this.sample && !this.project) return false;
    this.sample = ''; this.project = '';
    lsSet(VP_CTX_SAMPLE_KEY, ''); lsSet(VP_CTX_PROJECT_KEY, '');
    this._emit();
    return true;
  },
  _emit() {
    renderNavCtx();
    try {
      document.dispatchEvent(new CustomEvent('vp-ctx',
        {detail: {sample: this.sample, project: this.project}}));
    } catch (e) { /* 老浏览器兜底：忽略广播 */ }
  },
};
VP_CTX.load();
window.VP_CTX = VP_CTX;      // 各页内联脚本与 tools.html 选样弹窗取用

/* 项目筛选下拉：以全局当前项目为默认值；用户在本页手动选过则尊重本页选择。
   返回生效的筛选值。所有页面共用一份实现，避免各页各写一遍导致口径不一。 */
function gxSyncProjSelect(sel, projs) {
  if (!sel) return '';
  let cur = sel.dataset.touched ? sel.value : (VP_CTX.project || '');
  if (cur && projs.indexOf(cur) < 0) cur = '';
  sel.innerHTML = `<option value="">${esc(t('gx.allProjects', '全部项目'))}</option>`
    + projs.map(p => `<option value="${esc(p)}">${esc(p)}</option>`).join('');
  sel.value = cur;
  return cur;
}

function injectGlobalCtx() {
  // 2026-10-21：导航改版（左侧栏 + 顶部工具条）。上下文胶囊属于「全局工具」，
  // 挂到顶栏的 .navctx-host；旧写法挂进 .navlinks，改版后会掉进竖排侧栏。
  const host = document.querySelector('.navctx-host') ||
               document.querySelector('.navlinks');
  if (!host || document.getElementById('navCtx')) return;
  const wrap = document.createElement('span');
  wrap.id = 'navCtx';
  wrap.className = 'navctx';
  wrap.innerHTML =
    `<span class="navctx-chip" id="navCtxSample" onclick="gxPickSample()"></span>`
    + `<span class="navctx-chip" id="navCtxProject" onclick="gxPickProject()"></span>`
    + `<button class="navctx-reset" id="navCtxReset" type="button" `
    + `onclick="gxGlobalReset(this)"></button>`;
  host.appendChild(wrap);
  // 旧版这里还要把 #langBtn / 📂 收进 .navtail 防换行；改版后它们本来就在
  // 顶栏的 .navtail 里，无需再搬。
  renderNavCtx();
}

function renderNavCtx() {
  const s = $('navCtxSample'), p = $('navCtxProject'), r = $('navCtxReset');
  if (!s || !p) return;
  const tip = t('gx.ctxTip', '');
  const none = t('gx.none', '未选择');
  const all = t('gx.allProjects', '全部项目');
  // 胶囊宽度有限（max-width 116px，窄窗口还会再收），长样品名会被省略号截断，
  // 所以 title 里必须给出**完整**名字，否则用户没法确认选中的是哪一个。
  s.className = 'navctx-chip' + (VP_CTX.sample ? '' : ' is-empty');
  s.title = `${t('gx.sample', '样品')}：${VP_CTX.sample || none}`
    + (tip ? '\n' + t('gx.pickSample', '') + ' — ' + tip : '');
  s.innerHTML = `<span class="k">${esc(t('gx.sample', '样品'))}</span>`
    + `<span class="v">${esc(VP_CTX.sample || none)}</span>`;
  p.className = 'navctx-chip' + (VP_CTX.project ? '' : ' is-empty');
  p.title = `${t('gx.project', '项目')}：${VP_CTX.project || all}`
    + (tip ? '\n' + t('gx.pickProject', '') + ' — ' + tip : '');
  p.innerHTML = `<span class="k">${esc(t('gx.project', '项目'))}</span>`
    + `<span class="v">${esc(VP_CTX.project || all)}</span>`;
  if (r) {
    r.textContent = t('gx.reset', '🧹 全局清理');
    r.title = t('gx.resetTitle', '');
  }
}

/* 通用模态框（全局上下文选择 / 全局清理确认共用） */
function gxModal(titleText, bodyHtml, onMount) {
  const box = document.createElement('div');
  box.style.cssText = 'position:fixed;inset:0;z-index:9999;background:rgba(0,0,0,.45);'
    + 'display:flex;align-items:center;justify-content:center';
  const panel = document.createElement('div');
  panel.style.cssText = 'background:var(--paper,#fff);color:var(--ink-700,#33475b);'
    + 'border-radius:10px;padding:16px 18px;min-width:460px;max-width:92vw;'
    + 'max-height:88vh;overflow:auto;box-shadow:0 10px 40px rgba(0,0,0,.25)';
  panel.innerHTML = `<div style="display:flex;align-items:center;gap:10px;margin-bottom:10px">`
    + `<b style="font-size:14px">${esc(titleText)}</b>`
    + `<button class="btn small" id="gxModalX" style="margin-left:auto">✕</button></div>`
    + `<div id="gxModalBody">${bodyHtml}</div>`;
  box.appendChild(panel);
  document.body.appendChild(box);
  const onKey = e => { if (e.key === 'Escape') close(); };
  function close() {
    document.removeEventListener('keydown', onKey);
    box.remove();
  }
  document.addEventListener('keydown', onKey);
  panel.querySelector('#gxModalX').onclick = close;
  box.addEventListener('click', e => { if (e.target === box) close(); });
  if (onMount) onMount(panel.querySelector('#gxModalBody'), close);
  return close;
}

async function gxFetchSamples() {
  try { return await (await fetch('/api/samples')).json(); }
  catch (e) { return null; }
}

async function gxPickSample() {
  const list = await gxFetchSamples();
  if (list === null) {
    toast(t('gx.resetFail', '加载失败'), '', {kind: 'failed', ttl: 8000});
    return;
  }
  if (!list.length) { alert(t('gx.noSamples', '暂无样品')); return; }
  gxModal(t('gx.pickSample', '切换当前样品'),
    `<p class="hint" style="margin:0 0 8px">${esc(t('gx.pickHint', ''))}</p>`
    + `<div class="gx-pick-row"><input type="text" id="gxQ" `
    + `placeholder="${esc(t('gx.searchPh', ''))}">`
    + `<button class="btn small" id="gxNone">${esc(t('gx.clearSample', '清空当前样品'))}</button></div>`
    + `<div class="gx-pick-list" id="gxList"></div>`,
    (body, close) => {
      const box = body.querySelector('#gxList');
      const q = body.querySelector('#gxQ');
      const draw = () => {
        const s = q.value.trim().toLowerCase();
        const shown = list.filter(x => !s || String(x.name).toLowerCase().indexOf(s) >= 0);
        box.innerHTML = shown.length ? shown.map(x => `
          <div class="gx-item" data-name="${esc(x.name)}">
            <span class="nm">${esc(x.name)}${x.name === VP_CTX.sample ? ' ✔' : ''}</span>
            ${x.project ? `<span style="color:var(--ink-400,#8b98a5)">🏷 ${esc(x.project)}</span>` : ''}
            <span class="meta">${x.done || 0}/${x.total || 7}</span>
          </div>`).join('') : `<div class="gx-empty">—</div>`;
        box.querySelectorAll('.gx-item').forEach(el => {
          el.onclick = () => { VP_CTX.setSample(el.dataset.name); close(); };
        });
      };
      q.oninput = draw;
      body.querySelector('#gxNone').onclick = () => {
        VP_CTX.setSample('');
        close();
      };
      q.focus();
      draw();
    });
}

async function gxPickProject() {
  const list = (await gxFetchSamples()) || [];
  const projs = [...new Set(list.map(x => x.project).filter(Boolean))].sort();
  gxModal(t('gx.pickProject', '切换当前项目'),
    `<p class="hint" style="margin:0 0 8px">${esc(t('gx.projHint', ''))}</p>`
    + `<div class="gx-pick-list" id="gxList"></div>`,
    (body, close) => {
      const box = body.querySelector('#gxList');
      const opts = [{v: '', label: t('gx.allProjects', '全部项目'), n: list.length}]
        .concat(projs.map(p => ({
          v: p, label: p,
          n: list.filter(x => x.project === p).length})));
      box.innerHTML = opts.map(o => `
        <div class="gx-item" data-v="${esc(o.v)}">
          <span class="nm">${esc(o.label)}${o.v === VP_CTX.project ? ' ✔' : ''}</span>
          <span class="meta">${o.n}</span>
        </div>`).join('');
      box.querySelectorAll('.gx-item').forEach(el => {
        el.onclick = () => { VP_CTX.setProject(el.dataset.v); close(); };
      });
    });
}

/* 全局清理：非破坏性 —— 停任务 + 清队列 + 清任务记录 + 清浏览器关联，
   结果文件（results/ tool_runs/ downloads/ submissions/ …）一律不动。 */
function gxGlobalReset(btn) {
  if (btn && btn.disabled) return;
  const li = k => `<li>${esc(t(k, ''))}</li>`;
  gxModal(t('gx.resetTitle', '全局清理'),
    `<p style="margin:0 0 8px">${esc(t('gx.resetIntro', ''))}</p>`
    + `<ul style="margin:0 0 10px 18px;padding:0;font-size:12.5px;line-height:1.9">`
    + li('gx.resetItem1') + li('gx.resetItem2') + li('gx.resetItem3')
    + li('gx.resetItem4') + `</ul>`
    + `<p style="margin:0 0 10px;padding:8px 10px;border-radius:8px;`
    + `background:var(--green-50,#eef4f0);color:var(--green-900,#14532d);font-size:12.5px">`
    + `✅ ${esc(t('gx.resetKeep', ''))}</p>`
    + `<label class="chk" style="display:flex;gap:6px;align-items:flex-start;font-size:12.5px">`
    + `<input type="checkbox" id="gxArch" checked> `
    + `<span>${esc(t('gx.resetArchive', ''))}</span></label>`,
    (body, close) => {
      const row = document.createElement('div');
      row.style.cssText = 'display:flex;gap:8px;justify-content:flex-end;margin-top:16px';
      row.innerHTML = `<button class="btn small" id="gxNo">`
        + `${esc(t('gx.resetCancel', '取消'))}</button>`
        + `<button class="btn primary" id="gxYes">`
        + `${esc(t('gx.resetConfirm', '确认清理'))}</button>`;
      body.appendChild(row);
      body.querySelector('#gxNo').onclick = close;
      const yes = body.querySelector('#gxYes');
      yes.onclick = async () => {
        yes.disabled = true;
        yes.textContent = t('gx.resetting', '⏳ 清理中…');
        let d = null;
        try {
          const r = await fetch('/api/global/reset', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            // confirm 是"防误触"令牌（非鉴权）：本接口在空请求体下就会做
            // 破坏性清理，任何对路由发空 POST 的巡检/脚本都会误伤。
            body: JSON.stringify(
              {confirm: 'reset',
               include_archive: !!body.querySelector('#gxArch').checked})});
          d = await r.json();
          if (!r.ok) throw new Error((d && d.error) || ('HTTP ' + r.status));
        } catch (e) {
          yes.disabled = false;
          yes.textContent = t('gx.resetConfirm', '确认清理');
          const s = String(e);
          // 最可能的两种失败：服务未重启（404）或端口断开。
          // 前者给一句能直接照做的提示，比裸 HTTP 404 有用。
          const msg = s.indexOf('404') >= 0
            ? t('gx.resetNeedRestart', s) : s;
          toast(t('gx.resetFail', '全局清理失败'), msg,
                {kind: 'failed', ttl: 14000});
          return;
        }
        gxClearBrowserCtx();
        // 结果经 sessionStorage 带到 reload 之后展示：整页刷新是复位全部
        // 页面级 DOM 态最可靠的方式，但会吃掉 toast，所以先存后 reload。
        try { sessionStorage.setItem('vp_reset_result', JSON.stringify(d)); }
        catch (e) { /* 存不下就直接刷新，结果以计数为准 */ }
        close();
        location.reload();
      };
    });
}

function gxClearBrowserCtx() {
  for (const k of VP_CTX_RESET_KEYS) {
    try { localStorage.removeItem(k); } catch (e) { /* ignore */ }
    try { sessionStorage.removeItem(k); } catch (e) { /* ignore */ }
  }
  VP_CTX.sample = ''; VP_CTX.project = '';
}

function gxResetDetail(d) {
  return t('gx.resetDetail', '')
    .replace('{c}', (d && d.cancelled_tasks) || 0)
    .replace('{q}', (d && d.cleared_queue) || 0)
    .replace('{r}', (d && d.removed_task_records) || 0);
}

/* 项目筛选下拉的变更入口：模板 onchange 直接调它。
   为什么不用 addEventListener：属性处理器在解析期就注册，早于
   DOMContentLoaded 里挂的监听器；只在 DOMContentLoaded 挂监听的话，
   页面自身的 onchange（会重刷下拉）先跑，会把用户刚选的旧值刷回去。
   这里先同步全局项目，再由 'vp-ctx' 广播驱动各页刷新列表。 */
function gxProjChanged(sel) {
  if (!sel) return;
  sel.dataset.touched = '1';
  VP_CTX.setProject(sel.value || '');
}

/* 结果中心：按全局项目筛选样品行（归档表不参与，它没有项目维度） */
function gxApplyResultsFilter() {
  const tbl = $('sampleTbl');
  const sel = $('resProjFilter');
  if (!tbl || !sel) return;
  const rows = [...tbl.querySelectorAll('tr[data-project]')];
  if (!rows.length) {
    sel.style.display = 'none';
    return;
  }
  const projs = [...new Set(rows.map(r => r.dataset.project).filter(Boolean))].sort();
  const cur = gxSyncProjSelect(sel, projs);
  rows.forEach(r => {
    const ok = !cur || r.dataset.project === cur;
    r.style.display = ok ? '' : 'none';
    const nxt = r.nextElementSibling;
    if (nxt && nxt.classList.contains('filerow-tr') && !ok) nxt.style.display = 'none';
  });
}

/* 全局上下文变更广播的订阅者：各页按需重挂自己。
   管道页重刷样品列表并在样品变化时切换选中样品；结果中心重筛。 */
document.addEventListener('vp-ctx', () => {
  const pf = $('projFilter');
  if (pf) {
    if ($('samples')) loadSamples();
    else if (typeof scLoadSamples === 'function') scLoadSamples();
  }
  if ($('samples') && typeof curSample !== 'undefined'
      && VP_CTX.sample !== (curSample || '')) {
    if (VP_CTX.sample) {
      selectSample(VP_CTX.sample);
    } else {
      curSample = null;
      loadSamples();
      const p = $('pipe');
      if (p) p.innerHTML = `<div class="card"><p class="hint" `
        + `style="margin:4px 0">${esc(t('pp.empty', ''))}</p></div>`;
    }
  }
  gxApplyResultsFilter();
});

/* 整页加载后统一接线：注入导航上下文条、回显上一轮清理结果、接上全局样品 */
function gxApplyPage() {
  injectGlobalCtx();
  try {
    const raw = sessionStorage.getItem('vp_reset_result');
    if (raw) {
      sessionStorage.removeItem('vp_reset_result');
      const d = JSON.parse(raw);
      let msg = gxResetDetail(d);
      if (d && d.active_left) {
        msg += ' · ' + t('gx.resetLeft', '').replace('{n}', d.active_left);
      }
      toast(t('gx.resetDone', '全局清理完成'), msg, {ttl: 9000});
    }
  } catch (e) { /* 结果回显失败不影响页面 */ }

  const pf = $('projFilter');
  // 管道页 / 样品页：项目下拉变更走模板 onchange="gxProjChanged(this)"；
  // 这里只兜住未来新增、漏写 onchange 的同类下拉。
  if (pf && !pf.getAttribute('onchange')) {
    pf.addEventListener('change', () => gxProjChanged(pf));
  }
  // 管道页：导航栏的 /pipeline 不带 ?sample=，以前就因此永远回不到上次的样品；
  // 现在用全局当前样品自动接上（?sample= 显式指定时仍以 URL 为准）。
  // 第二参 'auto'：恢复失败（上次样品已不在磁盘）走静默清理 + toast，
  // 不弹阻塞式 alert（见 app-batch.js selectSample）。
  if ($('samples') && VP_CTX.sample) {
    const q = new URLSearchParams(location.search);
    if (!q.get('sample') && typeof selectSample === 'function') {
      selectSample(VP_CTX.sample, 'auto');
    }
  }
  gxApplyResultsFilter();
}

// ---------------- 库状态 ----------------
async function loadDbs() {
  try {
    const r = await fetch('/api/dbs');
    const d = await r.json();
    const badge = (ok, label) =>
      `<span class="badge ${ok ? 'ok' : 'no'}">${label} ${ok ? '✔' : '✘'}</span>`;
    const html = badge(d.taxonomy.ready, 'NCBI Taxonomy')
               + badge(d.host.ready, t('hm.db.host'))
               + badge(d.virus.ready, t('hm.db.virus'));
    const put = (id, text) => { const el = $(id); if (el) el.innerHTML = text; };
    put('dbStatus', html);
    put('taxStat', badge(d.taxonomy.ready, 'Taxonomy'));
    put('hostStat', badge(d.host.ready, t('hm.db.host')));
    put('virusStat', badge(d.virus.ready, t('hm.db.virus')));
  } catch (e) {}
}

async function loadTools() {
  const box = $('toolList');
  if (!box) return;
  try {
    const r = await fetch('/api/tools');
    const d = await r.json();
    box.innerHTML = Object.entries(d.tools).map(([k, v]) =>
      `<div class="trow"><span>${esc(k)}</span>
       <span class="${v ? 'ok' : 'no'}">${v ? '✔ ' + esc(v.split('\\').pop()) : '✘ ' + t('bd.notFound', '未找到')}</span></div>`
    ).join('');
  } catch (e) {}
}

// ---------------- 管道（模块化逐级运行） ----------------
let curSample = null;
let curStages = [];          // 当前样品的管道状态（runUpTo 计算用）
let SERVER_DEFAULTS = {};    // 设置页的分析默认参数（渲染初值用）

async function loadDefaults() {
  try {
    const d = await (await fetch('/api/settings')).json();
    SERVER_DEFAULTS = d.defaults || {};
  } catch (e) { SERVER_DEFAULTS = {}; }
}

function defVal(key, fallback) {
  const v = SERVER_DEFAULTS[key];
  return (v === undefined || v === null || v === '') ? fallback : v;
}

async function loadSamples() {
  try {
    const r = await fetch('/api/samples');
    const list = await r.json();
    const box = $('samples');
    if (!box) return;
    if (!list.length) {
      box.innerHTML = `<p class="hint">${t('pp.noSamples')}</p>`;
      return;
    }
    // 项目筛选下拉（去重）：默认值取全局当前项目，用户在本页选过则用本页的
    const sel = $('projFilter');
    let filt = '';
    if (sel) {
      const projs = [...new Set(list.map(s => s.project).filter(Boolean))];
      filt = gxSyncProjSelect(sel, projs);
    }
    const shown = filt ? list.filter(s => s.project === filt) : list;
    box.innerHTML = shown.length ? shown.map(s => {
      const pct = Math.round(s.done / (s.total || 7) * 100);
      // 样品名来自磁盘目录名（可含引号等字符）：走 data-* + dataset 传参，
      // 内联 onclick="'${esc(x)}'" 会被 HTML 实体解码还原成引号闭合字符串
      return `<div class="sample-item ${s.name === curSample ? 'active' : ''}"
            data-sample="${esc(s.name)}" onclick="selectSample(this.dataset.sample)">
        <button class="btn small si-del" title="${t('rs.confirmClear', '清除结果（保留样品与输入信息）')}"
                data-sample="${esc(s.name)}" onclick="event.stopPropagation();clearSample(this.dataset.sample, this)">🧹</button>
        <button class="btn small danger si-del" title="${t('rs.del', '删除样品及其全部文件')}"
                data-sample="${esc(s.name)}" onclick="event.stopPropagation();deleteSample(this.dataset.sample, this)">🗑</button>
        <div class="nm"><span>${esc(s.name)}</span><span class="pct">${pct}%</span></div>
        ${s.project ? `<div class="hint" style="margin:2px 0 0">🏷 ${esc(s.project)}</div>` : ''}
        <div class="pr-bar"><div style="width:${pct}%"></div></div>
        <div class="pr">${t('pp.progress')} ${s.done}/${s.total} ${t('pp.steps')}</div></div>`;
    }).join('') : `<p class="hint">${t('pp.noSamples')}</p>`;
  } catch (e) { if (netFail(e)) setConnBanner(true); else console.error('加载失败:', e); }
}

// ---------------- 建库 ----------------
async function buildTaxonomy(btn) {
  let r;
  try {
    r = await fetch('/api/build_taxonomy', {method: 'POST'});
  } catch (e) { alert(t('c.startFail', '启动失败') + ': ' + e); return; }
  if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json().catch(() => ({}))).error || r.status)); return; }
  const d = await r.json();
  taskLogOpen.add(d.task);
  startPolling();
  watchTaskBtn(d.task, btn, '⏳ 构建中…', () => loadDbs());
}

async function buildHostDb(btn) {
  const taxid = +$('hostTaxid').value;
  if (!taxid) { alert(t('bd.needTaxid', '请填写宿主物种的 NCBI TaxID')); return; }
  const body = {
    genome: $('hostGenome').value.trim(),
    taxid,
    hash_capacity: $('hostCap').value || '256M',
    threads: +$('hostThreads').value || undefined,
    rebuild: $('hostRebuild').checked,
    clean_mid: $('hostCleanMid').checked,
  };
  /* 输出目录留空 = 后端按物种自动命名 host-db/<TaxID>_<源目录名>_host_db */
  const outDirEl = $('hostOutDir');
  if (outDirEl && outDirEl.value.trim()) body.out_dir = outDirEl.value.trim();
  let r;
  try {
    r = await fetch('/api/build_host_db', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(body)});
  } catch (e) { alert(t('c.startFail', '启动失败') + ': ' + e); return; }
  if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json().catch(() => ({}))).error || r.status)); return; }
  const d = await r.json();
  taskLogOpen.add(d.task);
  startPolling();
  if (d.out_dir && typeof toast === 'function') {
    toast(t('bd.hostOutTo', '输出目录: ') + d.out_dir, '', {ttl: 5000});
  }
  /* 建库成功后后端会把新库设为「当前宿主库」→ 刷新列表与徽章 */
  watchTaskBtn(d.task, btn, '⏳ 建库中…', () => { loadDbs(); loadHostDbs(); });
}

/* 在线下载宿主基因组（gget·Ensembl / NCBI taxon / NCBI accession）→
   完成后回填建库表单，上下游在数据库构建页一条龙衔接 */
/* ---------- 🔍 先查后下：datasets summary 预览（基因组查询下载页共用） ---------- */
async function dsSummaryQuery(out, body) {
  out.innerHTML = '⏳ ' + t('rd.qryRunning', '查询中…');
  let r;
  try {
    r = await fetch('/api/datasets_summary', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(body)});
  } catch (e) { out.innerHTML = '❌ ' + t('c.startFail', '启动失败') + ': ' + e; return null; }
  const d = await r.json().catch(() => ({}));
  if (!r.ok) { out.innerHTML = '❌ ' + t('rd.qryFail', '查询失败') + ': ' + esc((d && d.error) || r.status); return null; }
  if (!d.rows || !d.rows.length) {
    out.innerHTML = '⚠ ' + t('rd.qryNone', '没有匹配记录（换个 taxon/条件再试）');
    return null;
  }
  return d;
}

function dsQryNote(d) {
  return d.truncated
    ? t('rd.qryTrunc', '匹配数超出预览上限，显示前 {n} 条（下载会取完整集合）').replace('{n}', d.returned)
    : t('rd.qryMatched', '共 {n} 条匹配').replace('{n}', d.total_count || d.returned);
}

function dsQryTable(head, rows) {
  return '<div style="max-height:380px;overflow:auto;margin-top:4px">'
    + '<table class="tbl"><thead><tr>'
    + head.map(h => `<th>${h}</th>`).join('')
    + '</tr></thead><tbody>' + rows.join('') + '</tbody></table></div>';
}

function fmtLen(n) { return n ? (+n).toLocaleString() : '-'; }

function renderHostAsmTable(d, out) {
  const head = [t('rd.qryThAcc', 'Accession'), t('rd.qryThName', '组装名称'),
    t('rd.qryThOrg', '物种'), t('rd.qryThLen', '长度 (bp)'),
    t('rd.qryThLevel', '组装级别'), t('rd.qryThSrc', '来源'),
    t('rd.qryThRel', '发布日期'), ''];
  const rows = d.rows.map(x => '<tr>'
    + `<td><code>${esc(x.accession)}</code></td>`
    + `<td>${esc(x.name)}</td>`
    + `<td>${esc(x.organism)}${x.common ? ' (' + esc(x.common) + ')' : ''}</td>`
    + `<td class="num">${fmtLen(x.length)}</td>`
    + `<td>${esc(x.level)}</td><td>${esc(x.source)}</td><td>${esc(x.released)}</td>`
    + `<td><button class="btn" style="padding:2px 8px" onclick="pickHostAcc('${x.accession}')">${t('rd.qryPick', '选用')}</button></td></tr>`);
  out.innerHTML = '✅ ' + dsQryNote(d) + dsQryTable(head, rows);
}

/* ① 宿主基因组：查 NCBI 组装（taxon / TaxID / accession 三入口），行可选用 */
async function qryHostAsm(btn) {
  const acc = $('hostDlAcc').value.trim();
  const sp = $('hostDlSpecies').value.trim();
  const tid = $('hostDlTaxid').value.trim();
  const out = $('hostQryResult');
  if (!acc && !sp && !tid) {
    alert(t('bd.hostDlNeedSp', '请填写宿主物种拉丁名（如 Oryza sativa）或 Assembly Accession'));
    return;
  }
  const body = { kind: 'genome', limit: 200, accession: acc || undefined,
                 taxon: acc ? undefined : (sp || tid) };
  const d = await dsSummaryQuery(out, body);
  if (!d) return;
  renderHostAsmTable(d, out);
}

function pickHostAcc(acc) {
  $('hostDlAcc').value = acc;
  if (typeof toast === 'function') {
    toast(t('rd.qryPicked', '已填入 Assembly Accession') + ': ' + acc, '', {ttl: 4000});
  }
}

async function dlHostGenome(btn) {
  const species = $('hostDlSpecies').value.trim();
  let accession = ($('hostDlAcc') && $('hostDlAcc').value.trim()) || '';
  if (!species && !accession) {
    alert(t('bd.hostDlNeedSp', '请填写宿主物种拉丁名（如 Oryza sativa）或 Assembly Accession'));
    return;
  }
  const source = $('hostDlSource').value || 'gget';
  if (source === 'ncbi' && !accession) {
    /* 先查后下闸门：同一物种常有多个组装版本，盲下会把多个版本合并进
       同一个宿主库（互相污染）。单版本自动选用后继续；多版本渲染选用表
       并拦下本次下载，用户选定 accession 后再点下载。 */
    const out = $('hostQryResult') || $('hostDlResult');
    const d = await dsSummaryQuery(out, { kind: 'genome', taxon: species,
                                          limit: 200 });
    if (!d) return;                        // 查询失败/无结果：out 里已有原因
    if (d.returned > 1) {
      renderHostAsmTable(d, out);
      out.insertAdjacentHTML('afterbegin', '⚠ ' +
        t('rd.hostDlGate', '该物种匹配到 {n} 个组装/版本，请从下方表格选用一个 accession 再下载').replace('{n}', d.returned) + '<br>');
      return;
    }
    accession = d.rows[0].accession;
    if ($('hostDlAcc')) $('hostDlAcc').value = accession;
    if (typeof toast === 'function') {
      toast(t('rd.hostDlAuto', '该物种仅 1 个组装，已自动选用') + ': ' + accession,
            '', {ttl: 5000});
    }
  }
  const body = {
    species,
    source,
  };
  if (accession) body.accession = accession;
  if ($('hostDlChr') && $('hostDlChr').value.trim()) body.chromosomes = $('hostDlChr').value.trim();
  if ($('hostDlRelease').value) body.release = +$('hostDlRelease').value;
  if ($('hostDlTaxid').value) body.taxid = +$('hostDlTaxid').value;
  let r;
  try {
    r = await fetch('/api/host_genome_download', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(body)});
  } catch (e) { alert(t('c.startFail', '启动失败') + ': ' + e); return; }
  if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json().catch(() => ({}))).error || r.status)); return; }
  const d = await r.json();
  taskLogOpen.add(d.task);
  startPolling();
  /* 任务 result 带回 genome_fa/taxid：建库页同屏则回填表单（旧布局），
     「基因组查询下载」独立页则写入结果行并给出建库跳转 */
  watchTaskBtn(d.task, btn, '⏳ 下载中…', (snap) => {
    loadHostDbs();
    const res = snap && snap.result;
    if (res && res.genome_fa) {
      if ($('hostGenome')) {
        $('hostGenome').value = res.genome_fa;
        if (res.taxid) $('hostTaxid').value = res.taxid;
      }
      const box = $('hostDlResult');
      if (box) {
        box.innerHTML = `✅ ${esc(res.genome_fa)}` +
          (res.taxid ? ` · TaxID ${res.taxid}` : '') +
          ` <a href="/build#m-host">${t('rd.goHostBuild', '去构建宿主库 →')}</a>`;
      }
      if (typeof toast === 'function') {
        toast(t('bd.hostDlDone', '宿主基因组已就绪，可直接构建宿主库'), '', {ttl: 6000});
      }
    }
  });
}

/* 宿主库列表：当前生效的 + 已有的（按物种目录），可一键「设为当前」 */
async function loadHostDbs() {
  const box = $('hostDbList');
  if (!box) return;
  let d;
  try { d = await (await fetch('/api/host_dbs')).json(); } catch (e) { return; }
  const items = d.items || [];
  if (!items.length) {
    box.innerHTML = `<p class="hint">${t('bd.hostNone', '还没有宿主库。填好基因组与 TaxID 后点「构建宿主库」。')}</p>`;
    return;
  }
  box.innerHTML = items.map(it => {
    const who = it.taxid
      ? `${it.species || '?'} <span class="hint">(taxid=${it.taxid})</span>`
      : `<span class="hint">${t('bd.hostAnon', '身份未记录（无 host_db.json 清单）')}</span>`;
    const tags = [];
    if (it.active) tags.push(`<span class="badge ok">${t('bd.hostActive', '当前')}</span>`);
    if (it.legacy) tags.push(`<span class="badge">${t('bd.hostLegacy', '旧布局')}</span>`);
    if (it.conflicted) tags.push(`<span class="badge no">${t('bd.hostConflict', '元数据冲突')}</span>`);
    if (!it.ready) tags.push(`<span class="badge no">${t('bd.hostNotReady', '不完整')}</span>`);
    const btn = it.active ? ''
      : ` <button class="btn small" onclick="setHostDb('${(it.selector || it.name).replace(/'/g, "\\'")}', this)">${t('bd.hostUse', '设为当前')}</button>`;
    const sub = it.conflicted
      ? `<div class="hint">⚠ seqid2taxid.map 含多个 taxid: ${(it.taxids_in_map || []).join(', ')}</div>`
      : '';
    return `<div style="margin:5px 0"><b>${it.name}</b> ${tags.join(' ')} — ${who}${btn}${sub}</div>`;
  }).join('');
}

async function setHostDb(name, btn) {
  const old = btn ? btn.disabled : false;
  if (btn) btn.disabled = true;
  try {
    const r = await fetch('/api/set_host_db', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({name})});
    const d = await r.json().catch(() => ({}));
    if (!r.ok) { alert(t('bd.opFail', '操作失败') + ': ' + (d.error || r.status)); return; }
    if (typeof toast === 'function') {
      toast(t('bd.hostSwitched', '✔ 当前宿主库已切换: ') + name, '', {ttl: 5000});
    }
    loadDbs(); loadHostDbs();
  } catch (e) {
    alert(t('bd.connFail', '无法连接平台服务') + ': ' + e);
  } finally { if (btn) btn.disabled = old; }
}

// ---------------- 结果页 ----------------
async function openDir(sample) {
  await fetch('/api/open_report_dir/' + encodeURIComponent(sample));
}

async function deleteSample(sample, btn) {
  const msg = t('pv.confirmDel', '删除样品「{name}」将连同其全部结果文件与输入档案一并移除，不可恢复。')
    .replace('{name}', sample);
  if (!confirm(msg)) return;
  btn.disabled = true;
  try {
    const r = await fetch(`/api/samples/${encodeURIComponent(sample)}/delete`,
                          {method: 'POST'});
    if (!r.ok) {
      alert(t('rs.delFail', '删除失败') + ': ' + ((await r.json()).error || ''));
      btn.disabled = false;
      return;
    }
    toast(t('rs.del', '🗑 删除') + ' · ' + sample, '', {ttl: 3000});
    btn.closest('tr')?.remove();
    btn.closest('.sample-item')?.remove();
    /* 流程页：删除的是当前选中样品时清除选中态并刷新列表 */
    if (typeof curSample !== 'undefined' && curSample === sample) {
      curSample = null;
      VP_CTX.setSample('');
    }
    loadSamples();   /* 无样品列表容器时内部自动跳过 */
  } catch (e) { btn.disabled = false; if (netFail(e)) setConnBanner(true); else console.error('样品操作失败:', e); }
}

/* 清除结果 ≠ 删除样品：只删分析产物，保留样品与输入档案（00_prep/input.json），
   之后重跑管道即从第一步重新分析 */
async function clearSample(sample, btn) {
  const msg = t('rs.confirmClear',
    '清除样品「{name}」的全部分析结果文件？（保留样品与输入信息，可重新分析）')
    .replace('{name}', sample);
  if (!confirm(msg)) return;
  btn.disabled = true;
  try {
    const r = await fetch(`/api/samples/${encodeURIComponent(sample)}/clear`,
                          {method: 'POST'});
    if (!r.ok) {
      alert(t('rs.clearFail', '清除失败') + ': ' + ((await r.json()).error || ''));
      btn.disabled = false;
      return;
    }
    toast(t('rs.cleared', '🧹 结果已清除') + ' · ' + sample, '', {ttl: 3000});
    btn.disabled = false;
    /* 流程页：刷新列表并重载当前管道视图；结果页：整页刷新完成度 */
    if ($('pipe')) {
      loadSamples();
      if (typeof curSample !== 'undefined' && curSample === sample) selectSample(sample);
    } else {
      location.reload();
    }
  } catch (e) { btn.disabled = false; if (netFail(e)) setConnBanner(true); else console.error('样品操作失败:', e); }
}

async function toggleFiles(sample, btn) {
  const row = btn.closest('tr').nextElementSibling;
  const box = row.querySelector('.filelist');
  if (row.style.display === 'none') {
    row.style.display = '';
    if (!box.dataset.loaded) {
      const r = await fetch('/api/sample_files/' + encodeURIComponent(sample));
      const d = await r.json();
      box.innerHTML = (d.files || []).map(f => {
        const segs = String(f.path).split('/').map(encodeURIComponent).join('/');
        const base = `/api/samples/${encodeURIComponent(sample)}/${segs}`;
        const edit = /\.(svg|png)$/i.test(f.path)
          ? ' ' + svgEditLink(base, '整理', 'mono') : '';
        return `<a href="${base}?dl=1" download>${esc(f.path)} (${esc(f.size)})</a>${edit}`;
      }).join('') || `<span class="hint">${t('rs.noFiles')}</span>`;
      box.dataset.loaded = '1';
    }
  } else {
    row.style.display = 'none';
  }
}

// ---------------- 主页总览 ----------------
async function loadHome() {
  try {
    const [r, rkv] = await Promise.all([
      fetch('/api/dbs'),
      fetch('/api/kv_index_list').catch(() => null),
    ]);
    const dbs = await r.json();
    // 病毒鉴定库（salmon 比对索引）不在 /api/dbs 体系内，单独探
    let kvLibs = [];
    try { kvLibs = ((await rkv.json()).libs) || []; } catch (e) {}
    // 四张库卡均指向 /build（库的查看与构建入口都在那里，与分析流程页无关）
    const defs = [
      ['0️⃣', t('hm.db.tax'), 'taxonomy', t('hm.db.tax.d'), '/build'],
      ['1️⃣', t('hm.db.host'), 'host', t('hm.db.host.d'), '/build'],
      ['🦠', t('hm.db.virus'), 'virus', t('hm.db.virus.d'), '/build'],
    ];
    let html = defs.map(([icon, name, key, desc, href]) => {
      const st = (dbs[key] || {}).ready;
      const label = st ? t('hm.ready')
        : (key === 'virus' ? t('hm.virusNotReady') : t('hm.notReady'));
      return `<a class="module-card" href="${href}">
        <span class="mc-state ${st ? 'st-ok' : 'st-no'}">${label}</span>
        <div class="mc-icon">${icon}</div>
        <div class="mc-name">${name}</div>
        <div class="mc-desc">${desc}</div></a>`;
    }).join('');
    // 第 4 张：病毒鉴定库（有库才称就绪）
    const kvSt = kvLibs.length > 0;
    const kvEng = (kvLibs[0] && kvLibs[0].salmon) ? 'salmon' : '';
    const kvDesc = t('hm.db.kvidx.d') + (kvSt && kvEng ? ` · ${kvEng}` : '');
    html += `<a class="module-card" href="/build">
      <span class="mc-state ${kvSt ? 'st-ok' : 'st-no'}">
        ${kvSt ? t('hm.db.kvidxN') + ' (' + kvLibs.length + ')' : t('hm.db.kvidxNone')}</span>
      <div class="mc-icon">🔍</div>
      <div class="mc-name">${t('hm.db.kvidx')}</div>
      <div class="mc-desc">${kvDesc}</div></a>`;
    $('dbCards').innerHTML = html;
  } catch (e) { if (netFail(e)) setConnBanner(true); else console.error('加载失败:', e); }
  try {
    const r = await fetch('/api/samples');
    const list = (await r.json()).slice(0, 6);
    $('homeSamples').innerHTML = list.length ? list.map(s => {
      const pct = Math.round(s.done / (s.total || 7) * 100);
      /* 2026-10-11：/pipeline 归档后样品卡改跳结果中心（按样品看产物） */
      return `<div class="module-card home-sample"
        onclick="location.href='/results'">
        <div class="mc-name">🧪 ${esc(s.name)}</div>
        <div class="pr-bar" style="height:6px;background:var(--line-100);border-radius:4px;overflow:hidden;margin:8px 0 4px">
          <div style="height:100%;width:${pct}%;background:linear-gradient(90deg,var(--green-600),#43b47c)"></div></div>
        <div class="mc-desc">${t('hm.progress')} ${s.done}/${s.total} ${t('pp.steps')}（${pct}%）</div></div>`;
    }).join('') : `<p class="hint">${t('hm.noSamples')}</p>`;
  } catch (e) { if (netFail(e)) setConnBanner(true); else console.error('加载失败:', e); }
}

// ---------------- 阶段结果预览（表格/JSON 就地弹窗，HTML 新窗口打开） ----------------
function previewStageFile(sample, stageDir, filename, stageKey) {
  // view 路径可含子目录（如 02b 的 filter/filtered.tsv），逐段编码
  const enc = filename.split('/').map(encodeURIComponent).join('/');
  // HTML 报告类（fastp/桑基/旭日/主报告）直接整页打开，体验更好
  if (filename.endsWith('.html')) {
    window.open(`/api/samples/${encodeURIComponent(sample)}/${encodeURIComponent(stageDir)}/${enc}`, '_blank');
    return;
  }
  const mask = document.createElement('div');
  mask.className = 'dlgmask';
  mask.style.zIndex = 120;
  mask.innerHTML = `<div class="dlg" style="width:900px">
    <div class="dlghead"><b>📄 ${esc(filename)}</b>
      <span>
        <button class="btn small" onclick="window.open('/api/samples/${encodeURIComponent(sample)}/${encodeURIComponent(stageDir)}/${enc}?dl=1','_blank')">⬇ 下载</button>
        <button class="btn small" onclick="this.closest('.dlgmask').remove()">✕</button>
      </span></div>
    <div class="pv-body" style="padding:6px 18px 16px;overflow:auto;max-height:70vh">
      <p class="hint">${t('c.loading')}</p></div></div>`;
  document.body.appendChild(mask);
  mask.addEventListener('click', e => { if (e.target === mask) mask.remove(); });
  const body = mask.querySelector('.pv-body');
  fetch(`/api/samples/${encodeURIComponent(sample)}/${encodeURIComponent(stageDir)}/${enc}`)
    .then(r => { if (!r.ok) throw new Error('HTTP ' + r.status); return r.text(); })
    .then(text => {
      if (filename.endsWith('.json')) {
        let obj;
        try { obj = JSON.parse(text); } catch (e) { obj = null; }
        body.innerHTML = obj
          ? `<pre style="background:#0d1f16;color:#b8e6c9;padding:12px;border-radius:8px;font-size:11.5px;white-space:pre-wrap">${esc(JSON.stringify(obj, null, 2))}</pre>`
          : `<pre>${esc(text.slice(0, 200000))}</pre>`;
        return;
      }
      // TSV/CSV → 表格（每页 20 行分页，全量行数保留）
      const rows = text.split(/\r?\n/).filter(l => l.trim());
      if (!rows.length) { body.innerHTML = `<p class="hint">${t('c.noData')}</p>`; return; }
      const sep = filename.endsWith('.csv') ? ',' : '\t';
      pvCells = rows.map(l => l.split(sep));
      pvPage = 1;
      renderPvTable(body);
    })
    .catch(e => { body.innerHTML = `<p class="err">${t('c.loadFail', '读取失败')}: ${esc(e)}</p>`; });
}

// 阶段文件预览表格分页 + 排序
let pvCells = null, pvPage = 1, pvSortI = -1, pvSortDir = 1;
function pvGo(p) {
  pvPage = p;
  const body = document.querySelector('.pv-body');
  if (body) renderPvTable(body);
}
function pvSortBy(i) {
  if (pvSortI === i) pvSortDir = -pvSortDir;
  else { pvSortI = i; pvSortDir = 1; }
  const body = document.querySelector('.pv-body');
  if (body) renderPvTable(body);
}
function renderPvTable(body) {
  if (!pvCells || !pvCells.length) return;
  const ncol = Math.max(...pvCells.map(c => c.length));
  let dataRows = pvCells.slice(1);
  if (pvSortI >= 0) {
    dataRows.sort((a, b) => {
      const va = a[pvSortI] || '', vb = b[pvSortI] || '';
      const na = parseFloat(va), nb = parseFloat(vb);
      const c = (!isNaN(na) && !isNaN(nb)) ? na - nb
        : String(va).localeCompare(String(vb), 'zh');
      return c * pvSortDir;
    });
  }
  const pages = Math.max(1, Math.ceil(dataRows.length / PER_PAGE));
  pvPage = Math.min(Math.max(1, pvPage), pages);
  const rows = dataRows.slice((pvPage - 1) * PER_PAGE, pvPage * PER_PAGE);
  let html = '<table class="tbl"><thead><tr>';
  for (let i = 0; i < ncol; i++) {
    const arrow = pvSortI === i ? (pvSortDir === 1 ? ' ▲' : ' ▼') : '';
    html += `<th class="sortable" onclick="pvSortBy(${i})">${esc((pvCells[0][i] || '').slice(0, 40))}${arrow}</th>`;
  }
  html += '</tr></thead><tbody>';
  for (const cells of rows) {
    html += '<tr>';
    for (let i = 0; i < ncol; i++) {
      const v = (cells[i] || '');
      html += `<td title="${esc(v.slice(0, 300))}">${esc(v.length > 60 ? v.slice(0, 57) + '…' : v)}</td>`;
    }
    html += '</tr>';
  }
  html += '</tbody></table>';
  html += `<p class="hint">共 ${dataRows.toLocaleString()} 行（不含表头）· 每页 ${PER_PAGE} 行</p>`
    + pagerHtml(pvPage, pages, 'pvGo');
  body.innerHTML = html;
}

// ---------------- 工具注册表 + 常用收藏（跨页共享） ----------------
const TOOL_REGISTRY = [
  { id: 't-convert',  href: '/tools#t-convert',  ic: '🔄', zh: '⓪ 格式转换',         en: '⓪ Format convert' },
  { id: 't-fastp',    href: '/tools#t-fastp',    ic: '🧹', zh: '① 序列质控',         en: '① Sequence QC' },
  { id: 't-identify', href: '/tools#t-identify', ic: '🦠', zh: '② 病毒识别和分类',   en: '② Virus classify' },
  { id: 't-assemble', href: '/tools#t-assemble', ic: '🧩', zh: '③ 病毒组装',         en: '③ Assembly' },
  { id: 't-contigs',  href: '/tools#t-contigs',  ic: '🔎', zh: '④ 病毒 contig 深度分析', en: '④ Contig deep-dive' },
  { id: 't-seqprep',  href: '/tools#t-seqprep',  ic: '⬇', zh: '⑤ 参考序列下载',    en: '⑤ NCBI references' },
  { id: 't-align',    href: '/tools#t-align',    ic: '🔤', zh: '⑦ 多序列比对/MSA 查看', en: '⑦ MSA viewer' },
  { id: 't-treebuild', href: '/tools#t-treebuild', ic: '🌳', zh: '⑧ 进化树查看器',     en: '⑧ Tree viewer' },
  { id: 't-sdt',      href: '/tools#t-sdt',      ic: '📐', zh: '⑨ SDT 分析和绘制',   en: '⑨ SDT matrix' },
  { id: 't-seqview',  href: '/results#seqview',  ic: '📄', zh: '序列查看器',         en: 'Sequence viewer' },
  { id: 'dl',         href: '/download',          ic: '📥', zh: '公共数据下载',       en: 'Public data downloads' },
];

function getFavs() {
  try { return JSON.parse(localStorage.getItem('vp_fav_tools') || '[]'); }
  catch (e) { return []; }
}

function toggleFav(id) {
  const favs = getFavs();
  const i = favs.indexOf(id);
  if (i >= 0) favs.splice(i, 1); else favs.push(id);
  try { localStorage.setItem('vp_fav_tools', JSON.stringify(favs)); } catch (e) {}
  document.querySelectorAll(`.starbtn[data-tool="${id}"]`)
    .forEach(b => b.classList.toggle('on', favs.includes(id)));
  if (typeof renderFavStrip === 'function') renderFavStrip();
  return favs.includes(id);
}

function toolName(tool) {
  return VP_LANG === 'en' ? tool.en : tool.zh;
}

function renderFavStrip() {
  const box = $('favStrip');
  if (!box) return;
  const favs = getFavs();
  const section = $('favSection');
  if (!favs.length) {
    box.innerHTML = '';
    if (section) section.style.display = 'none';
    return;
  }
  if (section) section.style.display = '';
  box.innerHTML = favs.map(id => {
    const tool = TOOL_REGISTRY.find(x => x.id === id);
    return tool ? `<a class="fav-chip" href="${tool.href}"><span class="ic">${tool.ic}</span>${esc(toolName(tool))}</a>` : '';
  }).join('');
}

// ---------------- 深色主题 ----------------
function applyTheme(theme) {
  if (theme === 'dark') document.documentElement.setAttribute('data-theme', 'dark');
  else document.documentElement.removeAttribute('data-theme');
  try { localStorage.setItem('vp_theme', theme); } catch (e) {}
}

// ---------------- 拖拽文件 → 上传 → 填路径 ----------------
/* 路径输入框统一接线：粘贴清洗 + 「最近使用」历史下拉。
   粘贴时把「复制文件地址」的引号 / file:// 前缀裁掉，失焦（change）再兜底
   洗一遍手输内容；聚焦/点击时弹出该输入框的历史路径，点选即回填。 */
function wirePathInput(input) {
  if (input.dataset.pathWired) return;
  input.dataset.pathWired = '1';
  input.addEventListener('paste', e => {
    const cd = e.clipboardData || window.clipboardData;
    const txt = cd ? cd.getData('text') : null;
    if (txt == null) return;
    const clean = cleanPathText(txt);
    if (clean && clean !== txt) {
      e.preventDefault();
      input.value = clean;
      input.dispatchEvent(new Event('input', { bubbles: true }));
      input.dispatchEvent(new Event('change', { bubbles: true }));
    }
  });
  input.addEventListener('change', () => {
    const clean = cleanPathText(input.value);
    if (clean && clean !== input.value) input.value = clean;
    pushPathHist(input.id, input.value);
  });
  input.addEventListener('focus', () => showPathHist(input));
  input.addEventListener('click', () => showPathHist(input));
  input.addEventListener('keydown', e => {
    if (e.key === 'Escape') closePathHist();
  });
}


// 找成对输入框：id 含 _r1/_r2、_R1/_R2 等成对标记时返回另一半个
function _pairSibling(input) {
  const id = input.id || '';
  const m = id.match(/^(.*_)([rR])(1|2)$/);
  if (!m) return null;
  const other = m[1] + m[2] + (m[3] === '1' ? '2' : '1');
  const el = $(other);
  return (el && el.tagName === 'INPUT') ? el : null;
}

// ── 桌面拖拽代理（scripts/drag_proxy.py 悬浮窗）轮询 ──
// 代理悬浮窗收到资源管理器拖入的真实路径后 POST 到平台；此轮询取走
// 并填入「最后点过的输入框」。仅在本页有输入行且窗口聚焦时轮询。
let _dropPollTs = 0;
function fillInputWithPaths(input, paths) {
  const sib = (typeof _pairSibling === 'function') ? _pairSibling(input) : null;
  if (paths.length >= 2 && sib) {
    input.value = paths[0]; sib.value = paths[1];
    pushPathHist(input.id, paths[0]); pushPathHist(sib.id, paths[1]);
  } else {
    input.value = paths.join(',');
    pushPathHist(input.id, input.value);
  }
  toast(t('dz.ok', '已导入上传文件'),
    '经拖拽代理填入: ' + paths.map(p => p.split(/[\/]/).pop()).join(', '),
    { ttl: 4000 });
}

setInterval(async () => {
  if (!document.hasFocus() || !document.querySelector('.filerow')) return;
  if (!window._lastFilerowInput) return;
  try {
    const d = await (await fetch('/api/dropped_paths/last?since=' + _dropPollTs)).json();
    if (d.ts) _dropPollTs = d.ts;
    if (d.fresh && d.paths && d.paths.length && window._lastFilerowInput) {
      fillInputWithPaths(window._lastFilerowInput, d.paths);
    }
  } catch (e) {}
}, 3000);

function wireDropzones() {
  document.querySelectorAll('.filerow').forEach(row => {
    if (row.dataset.dropWired) return;
    const input = row.querySelector('input[type=text]');
    if (!input) return;
    row.classList.add('drop-able');
    row.dataset.dropWired = '1';
    wirePathInput(input);
    input.addEventListener('focus', () => { window._lastFilerowInput = input; });
    row.addEventListener('dragover', e => {
      if ([...(e.dataTransfer.types || [])].includes('Files')) {
        e.preventDefault();
        row.classList.add('dragover');
      }
    });
    row.addEventListener('dragleave', () => row.classList.remove('dragover'));
    row.addEventListener('drop', async e => {
      row.classList.remove('dragover');
      e.preventDefault();
      // 文件夹拖放：浏览器拿不到本机文件夹路径，指引到 📂 按钮
      const entries = [...(e.dataTransfer.items || [])]
        .map(it => { try { return it.webkitGetAsEntry && it.webkitGetAsEntry(); }
                     catch (err) { return null; } })
        .filter(Boolean);
      if (entries.some(en => en.isDirectory)) {
        toast(t('dz.isDir', '拖入的是文件夹'),
          t('dz.isDirHint', '浏览器无法获取文件夹本机路径——请用行尾 📂 按钮选择目录'),
          { ttl: 6000 });
        return;
      }
      const files = [...(e.dataTransfer.files || [])];
      if (!files.length) return;               // 非文件拖放（如文本）不接管
      // R1/R2 成对输入框：拖 2 个文件时自动分装
      const sib = _pairSibling(input);
      if (files.length >= 2 && sib) {
        const up = async f => {
          const fd = new FormData(); fd.append('file', f);
          const r = await fetch('/api/upload', { method: 'POST', body: fd });
          const d = await r.json();
          if (!r.ok) throw new Error(d.error || 'upload failed');
          return d.path;
        };
        const old = input.placeholder;
        try {
          input.value = await up(files[0]);
          sib.value = await up(files[1]);
          pushPathHist(input.id, input.value);
          pushPathHist(sib.id, sib.value);
          input.placeholder = old;
          toast(t('dz.pair', 'R1/R2 已分别填入'), `${files[0].name} + ${files[1].name}`,
            { ttl: 4000 });
        } catch (err) { input.placeholder = old; toast(t('dz.fail', '上传失败'), String(err), {kind: 'failed', ttl: 6000}); }
        return;
      }
      // 常规：多文件全部上传，逗号拼接（供吃逗号列表的字段）
      const old = input.placeholder;
      input.placeholder = t('dz.uploading', '⬆ 上传中…') + ' ' + files[0].name;
      try {
        const paths = [];
        for (const f of files) {
          const fd = new FormData();
          fd.append('file', f);
          const r = await fetch('/api/upload', { method: 'POST', body: fd });
          const d = await r.json();
          if (!r.ok) throw new Error(d.error || 'upload failed');
          paths.push(d.path);
        }
        input.value = paths.join(',');
        pushPathHist(input.id, input.value);
        input.placeholder = old;
        toast(t('dz.ok', '已导入上传文件'), paths.map(p => p.split(/[\/]/).pop()).join(', '), {ttl: 4000});
      } catch (err) {
        input.placeholder = old;
        toast(t('dz.fail', '上传失败'), String(err), {kind: 'failed', ttl: 6000});
      }
    });
  });
  // 设置页的目录输入框（.s-ctrl，多数未写 type=text）不拖文件，
  // 但同样支持粘贴路径自动清洗
  document.querySelectorAll('.s-ctrl input').forEach(wirePathInput);
}

// ---------------- 序列查看器（专项分析页） ----------------
let svPath = '', svPage = 0;
async function seqviewLoad(page) {
  const p = ($('svPath')?.value || '').trim();
  const box = $('svBox');
  if (!p) { box.innerHTML = `<p class="hint">${t('sv.needPath', '请先选择或拖入 FASTA 文件')}</p>`; return; }
  if (p !== svPath) svPage = 0;
  svPage = Math.max(0, page || 0);
  box.innerHTML = `<p class="hint">${t('c.loading')}</p>`;
  try {
    const r = await fetch(`/api/seqview?path=${encodeURIComponent(p)}&page=${svPage}`);
    const d = await r.json();
    if (!r.ok) { box.innerHTML = `<p class="hint" style="color:#b91c1c">${esc(d.error || '读取失败')}</p>`; return; }
    svPath = p;
    $('svMeta').textContent =
      `${d.total} ${t('sv.records', '条')} · ${fmtSize(d.total_bp)}` +
      (d.truncated ? ` · ${t('sv.capped', '仅统计前 20000 条')}` : '');
    if (!d.rows || !d.rows.length) {
      box.innerHTML = `<p class="hint">${t('c.noData')}</p>`;
      $('svPager').style.display = 'none';
      return;
    }
    box.innerHTML = '<table class="tbl"><tr>' +
      `<th>ID</th><th>${t('sv.len', '长度')}</th><th>GC%</th><th>${t('sv.deg', '简并%')}</th><th>${t('sv.preview', '序列预览')} (300bp)</th></tr>` +
      d.rows.map(x => `<tr><td class="mono" title="${esc(x.id)}">${esc(x.id.length > 42 ? x.id.slice(0, 40) + '…' : x.id)}</td>` +
        `<td class="mono">${x.len.toLocaleString()}</td><td class="mono">${x.gc}</td><td class="mono">${x.deg}</td>` +
        `<td><span class="seqview-preview" title="${esc(x.preview || '')}">${esc(x.preview || '')}</span></td></tr>`).join('') +
      '</table>';
    $('svPager').style.display = 'flex';
    $('svPageInfo').textContent = `${t('pp.progress').split(' ')[0] || ''} ${svPage + 1} / ${d.pages}`;
  } catch (e) { box.innerHTML = `<p class="hint" style="color:#b91c1c">${t('c.loadFail', '读取失败')}: ${esc(e)}</p>`; }
}

// ---------------- 初始化 ----------------
window.addEventListener('DOMContentLoaded', () => {
  const q = new URLSearchParams(location.search);
  // 建库页
  if ($('hostTaxid') !== null) loadDbs();
  // 主页
  if ($('dbCards')) loadHome();
  // 分析管道页
  if ($('samples')) {
    loadSamples();
    loadQueue();
    setInterval(loadQueue, 4000);
    loadDefaults().then(() => {
      // 全局参数初值应用设置页默认
      const cf = $('confidence');
      if (cf && !cf.value && SERVER_DEFAULTS.confidence) cf.value = SERVER_DEFAULTS.confidence;
      if (q.get('sample')) selectSample(q.get('sample'));
    });
    if (q.get('focus')) setTimeout(() => {
      const el = $('sample'); if (el) { el.focus(); el.scrollIntoView({ behavior: 'smooth' }); }
    }, 200);
    const sub = $('subsample_on');
    if (sub) sub.addEventListener('change', () => {
      $('sub_row').style.display = sub.checked ? '' : 'none';
    });
    if (location.hash.startsWith('#grp-')) {
      setTimeout(() => {
        const g = $(location.hash.slice(1));
        if (g) g.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }, 800);
    }
  }
  if ($('tasks')) startPolling();
  // 全页通用：拖拽接线 / 收藏星 / 常用条
  wireDropzones();
  wireStars();
  renderFavStrip();
  // 全局上下文（当前样品/项目）接线 + 上一轮全局清理结果回显。
  // 放最后：各页自身的 init（loadSamples / 下拉填充）已完成，接线不会互相踩。
  gxApplyPage();
});

// ---------------- 收藏星渲染（专项分析各工具卡标题） ----------------
function wireStars() {
  // id 挂在 <section class="card" id="t-..."> 上（不是 h2）；旧选择器
  // 'section.card > h2[id]' 永远选不中 → 收藏星从不渲染、"常用工具"整块失效。
  document.querySelectorAll('section.card[id^="t-"]').forEach(card => {
    const id = card.id;
    if (!id.startsWith('t-') || card.querySelector('.starbtn')) return;
    const host = card.querySelector('h2') || card;
    const btn = document.createElement('button');
    btn.className = 'starbtn';
    btn.type = 'button';
    btn.dataset.tool = id;
    btn.title = t('fav.toggle', '加入/移出常用工具');
    btn.textContent = '☆';
    btn.classList.toggle('on', getFavs().includes(id));
    btn.onclick = () => {
      const on = toggleFav(id);
      btn.textContent = on ? '★' : '☆';
      toast(t('fav.strip', '常用工具'),
            on ? t('fav.added', '★ 已加入常用（总览页可见）')
               : t('fav.removed', '已移出常用'), {ttl: 2000});
    };
    host.appendChild(btn);
  });
}

// ---------------- 顶部导航任务徽章 ----------------
// 「任务中心」链接旁显示运行中任务数；每 15s 静默刷新，
// 仅当页面渲染了徽章（topnav）时启用。
(function navTasksBadgeLoop() {
  const badge = document.getElementById('navTasksBadge');
  if (!badge) return;
  const tick = async () => {
    try {
      const list = await (await fetch('/api/tasks?logs=0')).json();
      const arr = Array.isArray(list) ? list : (list.active || []);
      const n = arr.filter(x => x.status === 'running').length;
      badge.textContent = n > 99 ? '99+' : (n || '');
      badge.style.display = n ? 'inline-block' : 'none';
      badge.title = n ? ('运行中任务 ' + n + ' 个') : '';
    } catch (e) { /* 服务不可达时静默 */ }
    setTimeout(tick, 15000);
  };
  tick();
})();

/* ── 样品选择对话框 + 按样品批量运行（tools / host_removal 等页共用）── */

async function kvPickSamples(targetId) {
  const target = $(targetId) ? targetId : 'kv_samples';
  // 拉样品列表 → 弹窗多选（复选框）
  let list = [];
  try {
    const r = await fetch('/api/samples');
    list = await r.json();
  } catch (e) { alert(t('ui.app.sampleListFail') + ': ' + e); return; }
  if (!list.length) { alert(t('ui.app.noSamplesYet')); return; }
  const cur = new Set((val(target) || '').split(',').map(s => s.trim()).filter(Boolean));
  // 全局上下文（顶部导航）：当前样品默认勾上、当前项目默认作为筛选条件，
  // 省掉"切到本模块后重新找一遍样品"的重复操作。
  const gxSample = (window.VP_CTX && VP_CTX.sample) || '';
  const gxProject = (window.VP_CTX && VP_CTX.project) || '';
  if (gxSample && !cur.size && list.some(x => x.name === gxSample)) cur.add(gxSample);
  const html = `
    <div style="padding:4px 2px">
      <p class="hint" style="margin:0 0 8px">${t('tk.pickSamplesHintPre', '勾选要分析的样品（已选 ')}${cur.size}）</p>
      <div class="filerow" style="margin:0 0 8px">
        <select id="kvPickProj" style="flex:1">
          <option value="">${t('gx.allProjects', '全部项目')}</option>
          ${[...new Set(list.map(x => x.project).filter(Boolean))].sort()
              .map(p => `<option value="${esc(p)}">${esc(p)}</option>`).join('')}
        </select>
      </div>
      <div class="scroll-md" style="border:1px solid var(--line);border-radius:6px">
        <table class="tbl" style="margin:0">
          <thead><tr><th style="width:40px">${t('tk.thPick', '选')}</th><th>${t('tk.thSampleName', '样品名')}</th><th>${t('tk.thProject', '项目')}</th><th>${t('c.thStatus', '状态')}</th></tr></thead>
          <tbody>${list.map((x, i) => `
            <tr data-proj="${esc(x.project || '')}">
              <td><input type="checkbox" class="kv-pick" value="${esc(x.name)}" ${cur.has(x.name) ? 'checked' : ''}></td>
              <td class="mono">${esc(x.name)}</td>
              <td>${esc(x.project || '')}</td>
              <td>${esc((x.last_status || '') + ' ' + (x.done || 0) + '/' + (x.total || 0))}</td>
            </tr>`).join('')}
          </tbody>
        </table>
      </div>
      <div style="margin-top:10px;display:flex;gap:8px">
        <button class="btn primary" id="kvPickOk">${t('c.ok', '确定')}</button>
        <button class="btn" id="kvPickCancel">${t('c.cancel', '取消')}</button>
        <button class="btn small" id="kvPickAll">${t('c.selectAll', '全选')}</button>
        <button class="btn small" id="kvPickNone">${t('c.clearAll', '清空')}</button>
      </div>
    </div>`;
  const box = document.createElement('div');
  box.style.cssText = 'position:fixed;inset:0;z-index:9999;background:rgba(0,0,0,.45);display:flex;align-items:center;justify-content:center';
  const panel = document.createElement('div');
  panel.style.cssText = 'background:var(--bg-card,#fff);border-radius:10px;padding:16px;min-width:520px;max-width:92vw;box-shadow:0 10px 40px rgba(0,0,0,.25)';
  panel.innerHTML = html;
  box.appendChild(panel);
  document.body.appendChild(box);
  const close = () => box.remove();
  panel.querySelector('#kvPickCancel').onclick = close;
  // 项目筛选：默认落到全局当前项目；全选/清空只作用于当前可见（已筛选）的行
  const projSel = panel.querySelector('#kvPickProj');
  if (projSel && gxProject
      && [...projSel.options].some(o => o.value === gxProject)) {
    projSel.value = gxProject;
  }
  const applyProj = () => {
    const p = projSel ? projSel.value : '';
    panel.querySelectorAll('tbody tr[data-proj]').forEach(tr => {
      tr.style.display = (!p || tr.dataset.proj === p) ? '' : 'none';
    });
  };
  if (projSel) projSel.onchange = applyProj;
  applyProj();
  panel.querySelector('#kvPickAll').onclick = () =>
    panel.querySelectorAll('tbody tr[data-proj]').forEach(tr => {
      if (tr.style.display !== 'none') tr.querySelector('.kv-pick').checked = true;
    });
  panel.querySelector('#kvPickNone').onclick = () =>
    panel.querySelectorAll('tbody tr[data-proj]').forEach(tr => {
      if (tr.style.display !== 'none') tr.querySelector('.kv-pick').checked = false;
    });
  panel.querySelector('#kvPickOk').onclick = () => {
    const picked = [...panel.querySelectorAll('.kv-pick:checked')].map(c => c.value);
    $(target).value = picked.join(',');
    close();
  };
  box.onclick = ev => { if (ev.target === box) close(); };
}

/* ── 按样品批量运行（样品级模块通用）──────────────────────────
   样品名列表（逗号分隔，📋 多选对话框带项目筛选）→ 逐个读样品登记的
   R1/R2 → 逐个提交工具任务。任务服务器自带排队（heavy 顺序执行），
   任务中心每样品一张卡，输入/输出/进度天然按样品衔接。
   makeParams(name, info) 返回该样品的任务参数；返回 null 或抛错则跳过。 */
async function batchRunForSamples(btn, tool, field, makeParams) {
  const names = (val(field) || '').split(',').map(s => s.trim()).filter(Boolean);
  if (!names.length) { alert(t('tk.batchNeed', '请先用 📋 选择样品（可按项目全选）')); return; }
  return withBtn(btn, async () => {
    let ok = 0;
    const skip = [];
    for (const nm of names) {
      let d;
      try {
        const r = await fetch('/api/pipeline/' + encodeURIComponent(nm));
        if (!r.ok) { skip.push(nm + '(' + t('c.loadFail', '加载失败') + ')'); continue; }
        d = await r.json();
      } catch (e) { skip.push(nm + '(' + t('sc.connFail', '连接失败') + ')'); continue; }
      let p = null;
      try { p = makeParams(nm, d); } catch (e) { skip.push(nm + '(' + e.message + ')'); continue; }
      if (!p) continue;
      try {
        const r = await fetch('/api/tool/run', {
          method: 'POST', headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ tool, params: p }) });
        if (!r.ok) {
          const e = await r.json().catch(() => ({}));
          skip.push(nm + '(' + (e.error || t('tk.batchSubmitFail', '提交失败')) + ')');
          continue;
        }
        ok++;
      } catch (e) { skip.push(nm + '(' + t('sc.connFail', '连接失败') + ')'); }
    }
    if (typeof startPolling === 'function') startPolling();
    toast(t('tk.batchDone', '批量提交'),
      t('tk.batchSummary', '已提交 {n} 个任务（任务中心看进度）')
        .replace('{n}', ok)
      + (skip.length ? '；' + t('tk.batchSkipped', '跳过') + ' ' + skip.length
        + '：' + skip.join('、') : ''),
      { ttl: 8000 });
  }, t('tk.batchSubmitting', '⏳ 提交中…'));
}

/* 按样品登记的输入构造 r1/r2 参数（双端 {r1,r2}；单端 {r1}） */
function _peParams(nm, d) {
  if (!d.r1) throw new Error(t('tk.batchNoInput', '样品无输入序列'));
  return d.r2 ? { r1: d.r1, r2: d.r2 } : { r1: d.r1 };
}

