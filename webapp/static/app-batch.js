// ---------------- 管道卡片 / 批处理队列 ----------------
// 2026-09-13 自 app.js 拆出。普通 script 共享全局作用域：本文件在
// app.js 之后加载，可直接使用 $/esc/t/toast/jfetch 等核心工具；页面
// 内联脚本与其它 app-*.js 里的同名调用按全局解析（加载顺序见模板）。
// 注意：本块内引用的其它 app-*.js 函数均为调用期解析（事件/回调触发），
// 不存在加载期交叉依赖。
// （原「批量导入（TSV）」入口已于同日移除：与「从文件夹批量导入」重复，
//  其独有能力——分散路径/自定义命名/创建后入队——已并入文件夹扫描流程。）

// 入队后不会自动跑（2026-09-20 用户要求手动启动）：
//   queued = 待启动（入队初始态，点「▶ 启动」才开跑）
//   ready  = 排队中（已点启动，等前面的跑完 / 等恢复）
//   running/done/failed/cancelled 同旧语义
const Q_ST = { queued: t('stt.queued'), ready: t('stt.ready'),
               running: t('stt.running'), done: t('stt.done'),
               failed: t('stt.failed'), cancelled: t('stt.cancelled') };
async function loadQueue() {
  const box = $('queueBox');
  if (!box) return;
  try {
    const q = await (await fetch('/api/queue')).json();
    syncQueuePauseUI(!!q.paused);
    if (!q.items || !q.items.length) {
      box.innerHTML = `<p class="hint">${t('pp.queueEmpty')}</p>`;
      return;
    }
    const stLabel = it => {
      // 暂停中的就绪条目明确标注「已暂停」，否则用户分不清是排队还是卡死
      if (it.status === 'ready' && q.paused) return t('pp.qPausedItem', '已暂停');
      return Q_ST[it.status] || it.status;
    };
    box.innerHTML = '<table class="tbl dl-mini"><tr>' +
      `<th>${t('pp.qSample')}</th><th>${t('pp.qStatus')}</th><th>${t('pp.qOps')}</th></tr>` +
      q.items.map(it => {
        const st = stLabel(it);
        const statCls = it.status === 'done' ? 'done'
          : it.status === 'failed' ? 'failed'
          : it.status === 'running' ? 'running' : 'queued';
        const ops = it.status === 'queued'
          ? `<button class="btn small primary" onclick="queueStart('${esc(it.id)}')" title="${esc(t('pp.qStartTip', '开始运行该样品（入队后不会自动开跑）'))}">▶ ${esc(t('pp.qStart', '启动'))}</button> ` +
            `<button class="btn small danger" onclick="queueRemove('${esc(it.id)}')">✖</button>`
          : it.status === 'ready'
            ? `<button class="btn small danger" onclick="queueRemove('${esc(it.id)}')" title="${esc(t('pp.qRemoveReadyTip', '移出排队（未开跑可移除）'))}">✖</button>`
            : '';
        return `<tr><td><b>${esc(it.sample)}</b>` +
          (it.project ? ` <span class="hint" style="margin:0">🏷 ${esc(it.project)}</span>` : '') +
          `</td><td><span class="tstat ${statCls}">${st}</span>` +
          (it.error ? `<div class="hint" style="color:#b03a2e;margin:2px 0 0">${esc(it.error)}</div>` : '') +
          `</td><td>${ops}</td></tr>`;
      }).join('') + '</table>';
  } catch (e) { /* 静默，下一次轮询重试 */ }
}

async function queueStart(id) {
  try {
    const r = await fetch(`/api/queue/${encodeURIComponent(id)}/start`, {method: 'POST'});
    if (!r.ok) {
      const d = await r.json().catch(() => ({}));
      toast(t('pp.qStartFail', '启动失败'), d.error || (r.status + ''), {kind: 'failed', ttl: 12000});
      loadQueue();
      return;
    }
    if (QUEUE_PAUSED) {
      toast(t('pp.qReadyPaused', '已就绪'),
            t('pp.qReadyPausedHint', '队列处于暂停状态，点「▶ 恢复队列」后才会开跑'), {ttl: 6000});
    }
  } catch (e) {
    toast(t('pp.qStartFail', '启动失败'), String(e), {kind: 'failed', ttl: 12000});
  }
  loadQueue();
}

async function queueStartAll() {
  try {
    const r = await fetch('/api/queue/start_all', {method: 'POST'});
    if (!r.ok) {
      const d = await r.json().catch(() => ({}));
      toast(t('pp.qStartFail', '启动失败'), d.error || (r.status + ''), {kind: 'failed', ttl: 12000});
      loadQueue();
      return;
    }
    const d = await r.json();
    if (d.started > 0) {
      toast(t('pp.qStartedAll', '已全部启动'),
            t('pp.qStartedN', '{n} 个样品进入队列开始依序运行').replace('{n}', d.started),
            {ttl: 5000});
    }
  } catch (e) {
    toast(t('pp.qStartFail', '启动失败'), String(e), {kind: 'failed', ttl: 12000});
  }
  loadQueue();
}

let QUEUE_PAUSED = false;   // 最近一次 /api/queue 快照的暂停态（切换按钮用）
function syncQueuePauseUI(paused) {
  QUEUE_PAUSED = !!paused;
  const btn = $('queuePauseBtn');
  if (!btn) return;
  btn.textContent = paused ? t('pp.qResume', '▶ 恢复队列') : t('pp.qPause', '⏸ 暂停队列');
  const hint = $('queuePauseHint');
  if (hint) {
    hint.textContent = paused
      ? t('pp.qPausedHint', '已暂停：新入队样品不会自动开跑（正在运行的样品不受影响）')
      : '';
  }
}

async function queueTogglePause() {
  try {
    await fetch(QUEUE_PAUSED ? '/api/queue/resume' : '/api/queue/pause', {method: 'POST'});
  } catch (e) {
    toast(t('pp.qPauseFail', '操作失败'), String(e), {kind: 'failed', ttl: 12000});
  }
  loadQueue();
}

async function queueRemove(id) {
  try {
    const r = await fetch(`/api/queue/${encodeURIComponent(id)}/remove`, {method: 'POST'});
    if (!r.ok) {
      const d = await r.json().catch(() => ({}));
      toast(t('pp.qRemoveFail', '移除条目失败'), d.error || (r.status + ''), {kind: 'failed', ttl: 12000});
    }
  } catch (e) {
    toast(t('pp.qRemoveFail', '移除条目失败'), String(e), {kind: 'failed', ttl: 12000});
  }
  loadQueue();
}

async function queueClear() {
  try {
    const r = await fetch('/api/queue/clear_finished', {method: 'POST'});
    if (!r.ok) {
      const d = await r.json().catch(() => ({}));
      toast(t('pp.qClearFail', '清除已完成条目失败'), d.error || (r.status + ''), {kind: 'failed', ttl: 12000});
    }
  } catch (e) {
    toast(t('pp.qClearFail', '清除已完成条目失败'), String(e), {kind: 'failed', ttl: 12000});
  }
  loadQueue();
}

/* source='auto'：页面加载时自动恢复上次选中的样品（app.js 初始化调用）。
   这条路径失败几乎都是 localStorage 残留的样品名已不在磁盘上（被删除/
   归档/换过数据根），属过期状态：清回全局上下文（vp-ctx 既有逻辑会
   curSample=null、重载样品列表、恢复占位文案）+ toast 提示，不能弹
   阻塞式 alert——alert 会卡住事件循环，样品列表就停在「加载中…」。
   手动点击/URL 显式指定（source 缺省）失败时才弹 alert，且带上样品名
   与后端返回的错误详情。 */
async function selectSample(name, source) {
  curSample = name;
  // 写进全局上下文（跨模块共享）；取代原来只写不读的 vp_last_sample
  VP_CTX.setSample(name || '');
  loadSamples();
  try {
    const r = await fetch('/api/pipeline/' + encodeURIComponent(name) +
                          '?lang=' + encodeURIComponent(VP_LANG));
    if (!r.ok) {
      const err = (await r.json().catch(() => ({}))).error || ('HTTP ' + r.status);
      if (source === 'auto') {
        VP_CTX.setSample('');
        toast(t('pp.staleSample', '上次选中的样品已不可用'),
              t('pp.staleSampleD',
                '「{s}」已不在当前数据目录（可能已删除、归档或切换过数据目录），已自动取消选中')
                .replace('{s}', name) + '：' + err,
              { kind: 'failed', ttl: 10000 });
      } else {
        alert(t('c.loadFail', '读取失败') + '（' + name + '）：' + err);
      }
      return;
    }
    renderPipe(await r.json());
  } catch (e) { if (netFail(e)) setConnBanner(true); else console.error('pipeline 加载失败:', e); }
}

/* 管道 DAG 视图：③组装后的下游分支互不依赖（依赖关系见 Virus_Platform_Core/pipeline.py STAGE_REGISTRY），
   并排展示为紧凑节点，点击节点展开完整卡片（参数/产物/日志）。
   hostana 排在 verify 前：与「组装·kunpeng鉴定分类·宿主预测·再结构域鉴定」段口径一致
   （2026-09-20 STAGE_ORDER 已同步调整）。 */
const PIPE_FANOUT = ['hostana', 'verify', 'consensus', 'orf', 'orfa', 'phylo', 'primer', 'gbdraw'];

/* 分段一键：按业务段串联运行（stages 传给 runStages，段内顺序由后端
   topo_sort 依 STAGE_ORDER 定）。段与段之间靠 .done 断点衔接，
   某段失败只需重跑该段，不必从头再来。 */
const PIPE_TEMPLATES = [
  { id: 'prep', label: '⚡ 前处理',
    tip: '子采样(可选) → fastp 质控 → FASTQ→FASTA → 宿主去除',
    stages: ['subsample', 'fastp', 'fq2fa', 'host'] },
  { id: 'kv',   label: '⚡ 定量与共识',
    tip: '已知病毒识别与定量 → 共识序列与变异',
    stages: ['kvsuite', 'consensus'] },
  { id: 'asm',  label: '⚡ 组装·kunpeng鉴定分类·宿主预测·再结构域鉴定',
    tip: '③ 组装（kunpeng 鉴定分类·提取 viral_contigs）→ ④ ICTV 宿主预测 → ③b blastx/BFVD 再结构域鉴定',
    stages: ['assembly', 'hostana', 'verify'] },
  { id: 'anno', label: '⚡ 注释三连',
    tip: '⑥ ORF 预测 → ⑥b 功能注释 → ⑨ 基因组图谱；完成后可到「数据提交」页准备 NCBI 提交包',
    stages: ['orf', 'orfa', 'gbdraw'] },
];

function runTemplate(sample, id) {
  const tpl = PIPE_TEMPLATES.find(x => x.id === id);
  if (!tpl) return;
  let usable = curStages.filter(s => s.status !== 'unavailable').map(s => s.stage);
  if (tpl.exclude) usable = usable.filter(k => !tpl.exclude.includes(k));
  const stages = tpl.stages ? tpl.stages.filter(k => usable.includes(k)) : usable;
  if (!stages.length) { alert(t('pp.tplNone', '该段所需阶段当前不可用')); return; }
  runStages(sample, stages);
}

/* 注释产物 → 提交项目一键衔接（/api/submit/adopt_sample）：
   建/复用同名提交项目 + 导入 ③组装分类表 + 关联 orf 注释，
   成功后跳提交页并自动展开该项目（?open=）。 */
async function adoptSubmit(sample) {
  try {
    const r = await fetch('/api/submit/adopt_sample', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sample }) });
    const d = await r.json().catch(() => ({}));
    if (!r.ok) { alert((d.error || t('c.startFail', '操作失败'))); return; }
    const parts = [];
    if (d.imported) parts.push(t('pp.adoptImp', '导入 contig {n} 行').replace('{n}', d.imported.added));
    if (d.linked) parts.push(t('pp.adoptLink', '关联 CDS {n} 条').replace('{n}', d.linked.n_cds));
    if ((d.missing || []).length) parts.push(t('pp.adoptMiss', '缺：{m}').replace('{m}', d.missing.join('；')));
    toast(t('pp.submitGo', '一键生成提交项目（NCBI）'),
      (d.created ? t('pp.adoptNew', '已建档 ') : t('pp.adoptReuse', '已更新 '))
        + d.name + '：' + parts.join('，'), { ttl: 8000 });
    location.href = '/submit?open=' + encodeURIComponent(d.name);
  } catch (e) {
    if (netFail(e)) setConnBanner(true);
    alert(t('c.connFail', '无法连接平台服务') + ': ' + e);
  }
}

let curBranch = '';   // 当前展开的分支节点（fan-out 详情卡）
let cardIdx = 0;      // 阶段卡片序号（renderPipe 每次渲染前重置）

function toggleBranch(stage) {
  curBranch = curBranch === stage ? '' : stage;
  document.querySelectorAll('#pipe .pipe-branch-details > .stage-card').forEach(el => {
    el.style.display = el.dataset.stage === curBranch ? '' : 'none';
  });
  document.querySelectorAll('#pipe .pipe-node').forEach(el =>
    el.classList.toggle('on', el.dataset.stage === curBranch));
}

function renderPipe(d) {
  const stMap = {done: t('pp.done'), ready: t('pp.ready'),
                 blocked: t('pp.blocked'), unavailable: t('pp.unavailable'),
                 running: t('pp.running'), skipped: t('pp.skipped'),
                 failed: t('pp.failed')};
  curStages = d.stages;
  cardIdx = 0;                     // 每次渲染重置序号
  const byKey = {}; d.stages.forEach(s => byKey[s.stage] = s);
  const saved = collectParams();   // 重渲染前保存参数值
  const gmap = {};
  (d.groups || []).forEach(([g, ks]) => ks.forEach(k => gmap[k] = g));

  const card = (s, hidden) => {
    const i = cardIdx++;
    const canRun = s.status === 'ready' || s.status === 'done'
                || s.status === 'skipped' || s.status === 'failed';
    const isFresh = s.status === 'ready' || s.status === 'skipped'
                 || s.status === 'failed';
    const btns = canRun ? `
      <button class="btn small ${isFresh ? 'primary' : ''}"
              data-sample="${esc(d.sample)}"
              onclick="runStages(this.dataset.sample, ['${s.stage}'])">
        ${s.status === 'done' ? t('pp.rerun') : t('pp.runThis')}</button>
      <button class="btn small" data-sample="${esc(d.sample)}"
              onclick="runUpTo(this.dataset.sample, '${s.stage}')">
        ${t('pp.runTo')}</button>` :
      (s.status === 'unavailable' ?
        `<span class="hint" style="margin:0">${t('pp.installHint')}</span>` : '');
    const viewBtn = s.view
      ? ` <button class="btn small" data-sample="${esc(d.sample)}" data-dir="${esc(s.dir)}" data-view="${esc(s.view)}" data-stage="${esc(s.stage)}"
              onclick="previewStageFile(this.dataset.sample,this.dataset.dir,this.dataset.view,this.dataset.stage)">📊 ${t('c.view')}</button>` : '';
    const files = (s.outputs || []).map(o =>
      `<span>📄 ${esc(o.name.split('/').pop())} · ${fmtSize(o.size)}</span>`).join('');
    const params = renderStageParams(s.stage, saved);
    return `<div class="stage-card ${s.status}" data-stage="${s.stage}"${hidden ? ' style="display:none"' : ''}>
      <div class="stage-head"><span class="no">${i + 1}</span>
        <span class="nm">${esc(s.name)}</span>
        <span class="sbadge ${s.status}">${stMap[s.status] || s.status}</span></div>
      <div class="stage-sum">${esc(s.summary || t('pp.notRun'))}</div>
      ${files ? `<div class="stage-files">${files}</div>` : ''}
      ${params}
      <div class="stage-btns">${btns}${viewBtn}</div>
      <div class="stage-log" id="log-${s.stage}"></div></div>`;
  };
  const arrow = () => `<div class="pipe-arrow">↓</div>`;

  // 依赖分层：预处理链 → ②b 已知病毒识别与定量 → ③组装 → 下游并行分支 → ⑩报告；
  // 未来新增的未知阶段兜底追加在报告前。
  // 注意：② 病毒筛查（kraken2 分类）已于 2026-09-10 退役，阶段键由 'virus'
  // 换成 'kvsuite'；这里曾漏改，导致 ②b 卡片既不匹配 byKey.virus、
  // 也不在 known 集合里，被当成 extras 追加到 ⑩报告**之后**（顺序错乱）。
  const pre = ((d.groups || [])[0] || [null, []])[1].filter(k => byKey[k]).map(k => byKey[k]);
  const fanout = PIPE_FANOUT.map(k => byKey[k]).filter(Boolean);
  const known = new Set([...pre.map(s => s.stage), 'kvsuite', 'assembly', 'report', ...PIPE_FANOUT]);
  const extras = d.stages.filter(s => !known.has(s.stage));

  let html = `
    <div style="display:flex;justify-content:space-between;align-items:center;
                flex-wrap:wrap;gap:8px;margin-bottom:8px">
      <div><b style="font-size:17px;color:var(--green-900)">🧪 ${esc(d.sample)}</b>
        <div class="hint" style="margin:4px 0 0">${t('pp.input')}${esc(d.r1 || t('pp.none'))}
          ${d.r2 ? ' ＋ ' + esc(d.r2) : ` ${t('pp.single')}`}</div></div>
      <button class="btn" title="${t('pp.runAllTip', '依次运行全部剩余阶段（推荐改用下方分段一键，失败只重跑本段）')}"
              data-sample="${esc(d.sample)}"
              onclick="runStages(this.dataset.sample, null)">
        ${t('pp.runAll')}</button>
    </div>
    <div class="pipe-templates">
      <span class="hint" style="margin:0">${t('pp.tplHint', '⚡ 分段一键：')}</span>
      ${PIPE_TEMPLATES.map(tpl =>
        `<button class="btn small" title="${esc(tpl.tip)}"
                 data-sample="${esc(d.sample)}" onclick="runTemplate(this.dataset.sample, '${tpl.id}')">${esc(tpl.label)}</button>`).join('')}
    </div>` +
    // 功能注释完成后给「提交准备」一键衔接：建/复用同名提交项目 + 导入
    // ③组装分类表 + 关联 orf 注释，跳提交页自动展开（adoptSubmit）。
    (byKey.orfa && byKey.orfa.status === 'done' ? `
    <div class="pipe-templates" style="margin-bottom:8px">
      <span class="hint" style="margin:0">${t('pp.submitHint', '注释产物已就绪：')}</span>
      <button class="btn small" title="${t('pp.submitTip', '建/复用同名提交项目，导入组装分类表与注释产物；跳转后可直接生成 featuretable / 提交包')}"
              data-sample="${esc(d.sample)}"
              onclick="adoptSubmit(this.dataset.sample)">📦 ${t('pp.submitGo', '一键生成提交项目（NCBI）')}</button>
    </div>` : '');

  if (pre.length) {
    html += `<div class="group-title">${esc(gmap[pre[0].stage] || '')}</div>` +
      pre.map(s => card(s)).join(arrow()) + arrow();
  }
  if (byKey.kvsuite) {
    html += `<div class="group-title">${esc(gmap.kvsuite || '')}</div>` +
      card(byKey.kvsuite) + arrow();
  }
  if (byKey.assembly) {
    html += `<div class="group-title">${esc(gmap.assembly || '')}</div>` + card(byKey.assembly);
  }
  if (fanout.length) {
    html += `<div class="group-title">${esc(gmap.orf || '')} ·
        ${t('pp.branchHint', '并行分支（互不依赖，点节点展开参数与日志）')}</div>
      <div class="pipe-branch-row${fanout.length > 1 ? ' fanned' : ''}">
        ${fanout.map(s => `
          <div class="pipe-node ${s.status}" data-stage="${s.stage}"
               onclick="toggleBranch('${s.stage}')">
            <b>${esc(s.name)}</b>
            <span class="sbadge ${s.status}">${stMap[s.status] || s.status}</span>
            ${(s.status === 'ready' || s.status === 'done' ||
               s.status === 'skipped') ? `
              <button class="btn small" title="${t('pp.runThis')}"
                onclick="event.stopPropagation();if(curBranch!=='${s.stage}')toggleBranch('${s.stage}');
                         runStages('${esc(d.sample)}', ['${s.stage}'])">▶</button>` : ''}
          </div>`).join('')}
      </div>
      <div class="pipe-branch-details">
        ${fanout.map(s => card(s, true)).join('')}
      </div>` + arrow();
  }
  if (byKey.report) {
    html += `<div class="group-title">${esc(gmap.report || '')}</div>` + card(byKey.report);
  }
  extras.forEach(s => { html += arrow() + card(s); });
  $('pipe').innerHTML = html;
  hydrateKvLibSelects(saved.kv_lib || '');
}

/* ── ②b 鉴定库下拉：库列表来自 /api/kv_index_list（自动扫描平台内
   virusref_db/ + 外部登记），末尾带「✍ 自选目录…」可临时指定任意库目录
   （后端 resolve_kv_lib 统一解析库名/绝对路径）。pref = 重渲染前的选中值，
   保证翻卡片/重画管道时不丢用户已选的库。 */
const KVLIB_CUSTOM = '__custom__';
let _kvLibWired = new Set();

async function hydrateKvLibSelects(pref) {
  const sels = [...document.querySelectorAll('.kv-lib-select')];
  if (!sels.length) return;
  let libs = [];
  try {
    libs = (await (await fetch('/api/kv_index_list')).json()).libs || [];
  } catch (e) { /* 拉不到列表就只留默认项 + 自选 */ }
  const libVal = l => l.scope === 'platform' ? l.name : l.path;
  const usable = l => !!(l.salmon || l.minibwa || l.reference);
  const isDefault = l => l.scope === 'platform' && l.name === 'kv_index';
  const nonDefault = libs.filter(l => !isDefault(l));
  const def = libs.find(isDefault);
  sels.forEach(sel => {
    const cur = (pref !== undefined && pref !== null) ? pref : sel.value;
    sel.options.length = 0;
    sel.add(new Option(t('pp.kvLibDefault', '默认鉴定库 kv_index'), ''));
    nonDefault.forEach(l => {
      const tag = l.scope === 'external' ? t('pp.kvLibExternal', '（外部）')
                : l.scope === 'internal' ? t('pp.kvLibInternal', '（库根内）')
                : '';
      sel.add(new Option(
        l.name + tag +
        (l.salmon ? '' : ' · salmon✘') + (l.minibwa ? ' · minibwa✔' : ''),
        libVal(l)));
    });
    sel.add(new Option(t('pp.kvLibCustom', '✍ 自选目录…'), KVLIB_CUSTOM));
    let want = cur;
    if (want && want !== KVLIB_CUSTOM &&
        ![...sel.options].some(o => o.value === want)) {
      // 重渲染前选的是自选目录：作为动态选项保留
      sel.add(new Option('📂 ' + want + t('pp.kvLibCustomTmp', '（自选）'), want));
    }
    if (!want && !(def && usable(def))) {
      const first = nonDefault.find(usable);   // 默认库不可用 → 首个可用库
      if (first) want = libVal(first);
    }
    if (want && [...sel.options].some(o => o.value === want)) sel.value = want;
    if (!_kvLibWired.has(sel)) {
      _kvLibWired.add(sel);
      sel.addEventListener('change', () => {
        if (sel.value !== KVLIB_CUSTOM) return;
        const prev = sel.dataset.customPrev || '';
        pickDirWith(dir => {
          dir = String(dir || '').trim();
          if (!dir) { sel.value = prev; return; }
          let opt = [...sel.options].find(o => o.value === dir);
          if (!opt) {
            opt = new Option('📂 ' + dir + t('pp.kvLibCustomTmp', '（自选）'), dir);
            sel.add(opt);
          }
          sel.value = dir;
          sel.dataset.customPrev = dir;
        });
      });
    }
    sel.dataset.customPrev = sel.value;
  });
}

// 各阶段专属参数（渲染进对应卡片；label 为 i18n 键 pp.<key>）
const STAGE_PARAMS = {
  subsample: [
    { id: 'subsample_mode', label: 'pp.subMode', type: 'select', def: 'prefix',
      opts: [['prefix', 'prefix（前 N 对，快速冒烟）'],
             ['random', 'random（随机 N 对，等概率·seed 可复现）'],
             ['proportion', 'proportion（按比例随机，条数精确）']] },
    { id: 'subsample', label: 'pp.subPairs', type: 'number', def: 100000 },
    { id: 'subsample_prop', label: 'pp.subProp', type: 'number', def: 0.5 },
  ],
  fastp: [
    { id: 'fastp_dedup', label: 'pp.fastp_dedup', type: 'chk', def: false },
  ],
  fq2fa: [
    { id: 'do_fq2fa', label: 'pp.do_fq2fa', type: 'chk', def: true },
  ],
  host: [
    { id: 'db_host', label: 'pp.db_host', type: 'dir', def: '' },
  ],
  /* 病毒参考库目录挂在 ②b 上。以前它挂在已退役的 'virus' 阶段键下，
     而 pipeline_overview 永远不会产出 'virus' 卡片 → 输入框从不渲染，
     collectParams 取到 null，用户无法在管道页指定病毒库（只能用默认值）。
     kv_lib = ②b 鉴定库（salmon/minibwa 比对索引库，下拉可选/自选目录）；
     db_virus = kunpeng 分类库，供 ③b/④ 使用，与鉴定库正交。 */
  kvsuite: [
    { id: 'kv_lib', label: 'pp.kv_lib', type: 'kvlib', def: '' },
    { id: 'db_virus', label: 'pp.db_virus', type: 'dir', def: '' },
  ],
  assembly: [
    { id: 'assembly_mode', label: 'pp.assembly_mode', type: 'select', def: 'rnaviral',
      opts: [['rnaviral', 'rnaviral（RNA 病毒，默认）'], ['metaviral', 'metaviral'],
             ['rna', 'rna'], ['meta', 'meta'], ['isolate', 'isolate']] },
    { id: 'assembly_input', label: 'pp.assembly_input', type: 'select', def: 'virus',
      opts: [['virus', '② virus'], ['kept', '① host-free'], ['raw', 'raw']] },
    { id: 'memory', label: 'pp.memory', type: 'number', def: 64 },
    { id: 'min_contig_len', label: 'pp.min_contig_len', type: 'number', def: 200 },
  ],
  orf: [
    { id: 'min_orf_aa', label: 'pp.min_orf_aa', type: 'number', def: 100 },
  ],
  verify: [
    { id: 'verify_host', label: 'pp.verify_host', type: 'select', def: 'all',
      opts: [['all', 'all'], ['plants', 'plants'], ['fungi', 'fungi'],
             ['bacteria', 'bacteria'], ['invertebrates', 'invertebrates'],
             ['vertebrates', 'vertebrates'], ['algae', 'algae'],
             ['protists', 'protists'], ['archaea', 'archaea']] },
    { id: 'verify_methods', label: 'pp.verify_methods', type: 'multichk',
      def: ['blastx', 'bfvd'],
      opts: [['blastx', 'blastx'], ['bfvd', 'BFVD']] },
    { id: 'verify_combine', label: 'pp.verify_combine', type: 'select',
      def: 'union',
      opts: [['union', 'union'], ['intersection', 'intersection']] },
  ],
  phylo: [
    { id: 'do_trim', label: 'pp.do_trim', type: 'chk', def: true },
    { id: 'top_n_refs', label: 'pp.top_n_refs', type: 'number', def: 10 },
    { id: 'tree_tool', label: 'pp.tree_tool', type: 'select', def: 'fasttree',
      opts: [['fasttree', 'FastTree'], ['nj', 'NJ（快速）'], ['raxml-ng', 'RAxML-NG']] },
    { id: 'tree_sampling', label: 'pp.tree_sampling', type: 'select', def: 'blast',
      opts: [['blast', 'blast'], ['macro', 'macro'], ['genus', 'genus'],
             ['lineage', 'lineage']] },
    { id: 'ncbi_refs', label: 'pp.ncbi_refs', type: 'text', def: '' },
  ],
  primer: [
    { id: 'primer_mode', label: 'pp.primer_mode', type: 'select', def: 'conserved',
      opts: [['conserved', 'conserved'], ['plain', 'plain']] },
    { id: 'specificity', label: 'pp.specificity', type: 'chk', def: false },
  ],
  gbdraw: [
    { id: 'gbdraw_max', label: 'pp.gbdraw_max', type: 'number', def: 12 },
    { id: 'plot_engine', label: 'pp.plot_engine', type: 'select', def: 'auto',
      opts: [['auto', 'auto'], ['gbdraw', 'gbdraw'], ['dfv', 'DFV']] },
    { id: 'gbdraw_fasta', label: 'pp.gbdraw_fasta', type: 'file', def: '' },
    { id: 'gbdraw_ann', label: 'pp.gbdraw_ann', type: 'file', def: '' },
  ],
  report: [],
};

function renderStageParams(stage, saved) {
  const defs = STAGE_PARAMS[stage] || [];
  if (!defs.length) return '';
  const items = defs.map(p => {
    // 初值优先级：页面已填 > 设置页默认 > 代码缺省
    const dv = defVal(p.id, p.def);
    const cur = saved[p.id] !== undefined && saved[p.id] !== null ? saved[p.id] : dv;
    if (p.type === 'select') {
      const opts = p.opts.map(([v, txt]) =>
        `<option value="${v}" ${String(cur) === String(v) ? 'selected' : ''}>${txt}</option>`).join('');
      return `<div><label>${t(p.label)}</label><select id="${p.id}">${opts}</select></div>`;
    }
    if (p.type === 'kvlib') {
      // 鉴定库下拉：先渲染空壳，hydrateKvLibSelects() 异步填充库列表
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <select id="${p.id}" class="kv-lib-select"></select></div>`;
    }
    if (p.type === 'dir') {
      const ph = { db_host: t('pp.dbHostPh'), db_virus: t('pp.dbVirusPh') }[p.id]
                 || t('pp.defPh');
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <div class="filerow">
          <input id="${p.id}" type="text" value="${cur}" placeholder="${ph}">
          <button class="btn small" onclick="browseDir('${p.id}')">📁</button>
        </div></div>`;
    }
    if (p.type === 'file') {
      const ph = { gbdraw_fasta: t('pp.gFastaPh'),
                   gbdraw_ann: t('pp.gAnnPh') }[p.id] || '';
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <div class="filerow">
          <input id="${p.id}" type="text" value="${cur}" placeholder="${ph}">
          <button class="btn small" onclick="browse('${p.id}')">📁</button>
        </div></div>`;
    }
    if (p.type === 'chk') {
      const on = saved[p.id] !== undefined ? saved[p.id] : !!dv;
      return `<label class="chk" style="align-self:end;margin-bottom:8px">
        <input type="checkbox" id="${p.id}" ${on ? 'checked' : ''}> ${t(p.label)}</label>`;
    }
    if (p.type === 'multichk') {
      // 多选（如验证证据方法）：name=p.id，值取勾选项；全不勾时 collectParams 回传 undefined，
      // 服务端回退默认值。
      const vals = Array.isArray(cur) ? cur.map(String)
                 : String(cur || '').split(/[,;]/).filter(Boolean);
      const boxes = p.opts.map(([v, txt]) =>
        `<label class="chk" style="margin-right:14px">
          <input type="checkbox" name="${p.id}" value="${v}"
            ${vals.includes(String(v)) ? 'checked' : ''}> ${txt}</label>`).join('');
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <div>${boxes}</div></div>`;
    }
    if (p.type === 'text') {
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <input id="${p.id}" type="text" value="${cur}" placeholder="${t('pp.ncbi_ph')}"></div>`;
    }
    return `<div><label>${t(p.label)}</label>
      <input id="${p.id}" type="number" value="${cur}"></div>`;
  }).join('');
  return `<div class="stage-params">${items}</div>`;
}

function collectParams() {
  return {
    threads: +($('threads')?.value) || undefined,
    confidence: +($('confidence')?.value) || 0,
    chunk_dir: ($('chunk_dir')?.value || '').trim() || null,
    subsample: +($('subsample')?.value) || undefined,
    subsample_mode: $('subsample_mode')?.value || undefined,
    subsample_prop: +($('subsample_prop')?.value) || undefined,
    db_host: ($('db_host')?.value || '').trim() || null,
    db_virus: ($('db_virus')?.value || '').trim() || null,
    // ②b 鉴定库：''/自选占位未确认 = 默认 kv_index（resolve_kv_lib 解析）
    kv_lib: (() => {
      const v = ($('kv_lib')?.value || '').trim();
      return (!v || v === KVLIB_CUSTOM) ? null : v;
    })(),
    assembly_mode: $('assembly_mode')?.value,
    assembly_input: $('assembly_input')?.value,
    memory: +($('memory')?.value) || undefined,
    min_contig_len: +($('min_contig_len')?.value) || undefined,
    min_orf_aa: +($('min_orf_aa')?.value) || undefined,
    do_trim: $('do_trim') === null ? undefined : !!($('do_trim')?.checked),
    top_n_refs: +($('top_n_refs')?.value) || undefined,
    tree_tool: $('tree_tool')?.value,
    tree_sampling: $('tree_sampling')?.value,
    ncbi_refs: ($('ncbi_refs')?.value || '').trim() || null,
    primer_mode: $('primer_mode')?.value,
    specificity: !!($('specificity')?.checked),
    fastp_dedup: !!($('fastp_dedup')?.checked),
    do_fq2fa: $('do_fq2fa') === null ? undefined : !!($('do_fq2fa')?.checked),
    gbdraw_max: +($('gbdraw_max')?.value) || undefined,
    plot_engine: $('plot_engine')?.value,
    gbdraw_fasta: ($('gbdraw_fasta')?.value || '').trim() || null,
    gbdraw_ann: ($('gbdraw_ann')?.value || '').trim() || null,
    force: !!($('force')?.checked),
    verify_host: $('verify_host')?.value,
    verify_methods: (() => {
      const els = document.querySelectorAll('input[name="verify_methods"]:checked');
      return els.length ? [...els].map(e => e.value) : undefined;
    })(),
    verify_combine: $('verify_combine')?.value,
  };
}

async function runStages(sample, stages) {
  const body = Object.fromEntries(
    Object.entries(collectParams()).filter(([_k, v]) => v !== undefined));
  if (stages) body.stages = stages;
  if (body.force && !confirm(t('pp.forceConfirm',
      '强制重跑会忽略断点标记，所选步骤将重新执行。继续？'))) return;
  try {
    const r = await fetch('/api/pipeline/' + encodeURIComponent(sample) + '/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)});
    if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    const tb = $('tasks');
    if (tb) tb.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) { if (netFail(e)) setConnBanner(true); alert(t('c.connFail', '无法连接平台服务') + ': ' + e); }
}

function runUpTo(sample, stage) {
  const effective = curStages.filter(s => s.status !== 'unavailable');
  const firstUndone = effective.findIndex(s => s.status !== 'done');
  const effOrder = effective.map(s => s.stage);
  const start = firstUndone < 0 ? effOrder.indexOf(stage) : firstUndone;
  const end = effOrder.indexOf(stage);
  const stages = effOrder.slice(Math.min(start, end), end + 1);
  if (!stages.length) { alert(t('pp.alreadyDone', '该步骤已完成（可用「重跑此步」单独重跑）')); return; }
  runStages(sample, stages);
}

async function createSample(btn) {
  return withBtn(btn, async () => {
    const r1 = $('r1').value.trim();
    const errBox = $('cErr');
    errBox.textContent = '';
    if (!r1) { errBox.textContent = t('pp.needR1', '请选择 R1 FASTQ 文件'); return; }
    const body = {
      sample: $('sample').value.trim(),
      r1, r2: $('r2').value.trim() || null,
      project: $('sampleProject')?.value.trim() || null,
      subsample: ($('subsample_on')?.checked)
        ? (+($('subsample')?.value) || defVal('subsample', 100000)) : 0,
    };
    const res = await apiCreateSample(body);
    if (!res.ok) {
      if (res.conn) setConnBanner(true);
      errBox.textContent = res.error; return;
    }
    $('r1').value = ''; $('r2').value = ''; $('sample').value = '';
    // 后端可能把样品名规范化过（中文/空格 → '_'）：显式告知实际登记名，
    // 否则用户以为建的是「样品A」，列表里却是「A」，无从对应。
    if (res.data && res.data.note) {
      errBox.textContent = res.data.note;
    }
    selectSample(res.data.sample);
    loadSamples();
  }, '⏳ 创建中…');
}

