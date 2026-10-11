/* 植物病毒分析平台 · AI 助手对话坞（总览页）
 * 依赖：app.js 的 $/esc、i18n.js 的 t()（均先于本脚本加载）。
 * 流式：POST /api/ai/chat → SSE（data: {...}\n\n），fetch ReadableStream 逐块解析。
 * 服务端无状态：每次请求携带完整对话历史（客户端持有），上下文注入在服务端做。
 */
'use strict';

/* 防重复包含：本脚本随 _nav.html 注入每个页面，同一页只会加载一次；
 * 万一某页额外手滑再引一次，这里直接短路，避免重复绑定与双份事件流。 */
if (window.AI) {
  /* already loaded */
} else {

const AI_MSGS_KEY = 'vp_ai_msgs';
const AI_SEEN_KEY = 'vp_ai_lastseen';

const AI = {
  open: false, busy: false, ready: false, ctxLoaded: false,
  sample: '',            // 当前选中样品（'' = 不注入样品 digest）
  messages: [],          // {role:'user'|'assistant', content}
  _selFilled: false,
  lastEventId: 0,        // 本会话流水位（EventSource 重连自动带 Last-Event-ID）
  lastSeenId: 0,         // 跨刷新未读水位（localStorage）：看过的不再计未读
  unread: 0,             // 对话坞关闭期间的重要事件数
  _es: null,
  _ctrl: null,           // 当前生成的 AbortController
  _title0: null,         // 页面原始标题（未读计数前缀恢复用）
};

try { AI.lastSeenId = parseInt(localStorage.getItem(AI_SEEN_KEY) || '0', 10)
      || 0; } catch (e) {}

/* ---------- 开合与初始化 ---------- */

AI.toggle = function (force) {
  const want = (force === undefined) ? !AI.open : !!force;
  AI.open = want;
  $('aiDock').hidden = !want;
  $('aiFab').hidden = want;
  if (want) {
    AI.clearUnread();
    AI.ensureCtx();
    $('aiInput').focus();
  }
};

/* ---------- 主动播报事件流（三层驱动） ----------
 * EventSource 常驻订阅 /api/ai/events；重连自动带 Last-Event-ID 续传。
 * 坞关着时重要事件（AI 播报/失败/待决闸门）计未读、FAB 挂红点。 */
AI.startEvents = function () {
  if (AI._es || typeof EventSource === 'undefined') return;
  try { AI._es = new EventSource('/api/ai/events'); } catch (e) { return; }
  AI._es.onmessage = (m) => {
    let ev = null;
    try { ev = JSON.parse(m.data); } catch (e) { return; }
    if (!ev || ev.id === undefined) return;
    if (ev.id <= AI.lastEventId) return;      // 重连补发去重
    AI.lastEventId = ev.id;
    AI.renderEvent(ev);
  };
};
/* 离开页面时显式关掉 SSE：服务端 SSE 生成器没有客户端断连检测（写心跳到
 * 半关连接永远不会失败），导航后若不主动 close，旧连接的套接字两头都不释放
 * ——每次导航泄漏 ~2 条 ESTABLISHED，攒到 ~6 条（HTTP/1.1 同主机上限）后
 * 新 fetch 全部排队冻死（实测 2026-10-11：连续 6 页后 /api/samples 5s 超时、
 * kvsuite 选样弹窗打不开）。pagehide 在导航/刷新时都会触发，比 unload 可靠。 */
window.addEventListener('pagehide', function () {
  if (AI._es) { AI._es.close(); AI._es = null; }
});

AI.clearUnread = function () {
  AI.unread = 0;
  AI.lastSeenId = Math.max(AI.lastSeenId, AI.lastEventId);
  try { localStorage.setItem(AI_SEEN_KEY, String(AI.lastSeenId)); } catch (e) {}
  const b = $('aiFabBadge');
  if (b) { b.hidden = true; b.textContent = ''; }
  AI._titleSync();
};

AI._bumpUnread = function (ev) {
  if (AI.open) return;
  if (ev && ev.id <= AI.lastSeenId) return;   // 上一会话已看过（刷新重放）
  AI.unread += 1;
  const b = $('aiFabBadge');
  if (b) { b.hidden = false; b.textContent = AI.unread > 9 ? '9+' : String(AI.unread); }
  AI._titleSync();
};

// 浏览器标签页标题带未读数：用户在别的标签页也能看到有 AI 事件
AI._titleSync = function () {
  try {
    if (AI._title0 === null) AI._title0 = document.title.replace(/^\(\d+\)\s*/, '');
    document.title = AI.unread > 0 ? `(${AI.unread}) ${AI._title0}` : AI._title0;
  } catch (e) { /* 标题不可写则静默 */ }
};

AI.renderEvent = function (ev) {
  const s = ev.sample ? ` · ${ev.sample}` : '';
  if (ev.type === 'stage_start') {
    const n = document.createElement('div');
    n.className = 'ai-event ai-ev-start';
    n.innerHTML = `📖 <b>${esc(ev.name || ev.stage)}</b>${esc(s)}`
      + `<div class="ai-ev-note">${esc(ev.note || '')}</div>`;
    $('aiMsgs').appendChild(n);
    _scroll();
  } else if (ev.type === 'stage_done') {
    const n = document.createElement('div');
    n.className = 'ai-event ai-ev-done';
    n.textContent = `✔ ${ev.name || ev.stage} 完成${s}`;
    $('aiMsgs').appendChild(n);
    _scroll();
  } else if (ev.type === 'ai_read') {
    AI._bumpUnread(ev);
    const n = document.createElement('div');
    n.className = 'ai-msg ai-bot ai-read';
    n.innerHTML = `<div class="ai-read-head">✨ ${esc(t('ai.ev.read', 'AI 播报'))}`
      + ` · ${esc(ev.name || ev.stage || '')}</div>`
      + (ev.failure ? `<div class="ai-read-fail">✘ ${esc(t('ai.ev.failed', '阶段失败诊断'))}</div>` : '')
      + `<div class="ai-read-body">${aiFmt(ev.text || '')}</div>`;
    $('aiMsgs').appendChild(n);
    _scroll();
  } else if (ev.type === 'stage_failed') {
    AI._bumpUnread(ev);
    const n = document.createElement('div');
    n.className = 'ai-event ai-ev-fail';
    n.innerHTML = `✘ <b>${esc(ev.name || ev.stage)}</b> ${esc(t('ai.ev.failed', '阶段失败'))}${esc(s)}`
      + (ev.error ? `<div class="ai-ev-note">${esc(ev.error)}</div>` : '');
    $('aiMsgs').appendChild(n);
    _scroll();
  } else if (ev.type === 'run_done') {
    const n = document.createElement('div');
    n.className = 'ai-event ai-ev-done ai-ev-run';
    n.textContent = `🏁 ${t('ai.ev.runDone', '本轮分析全部完成')}${s}`;
    $('aiMsgs').appendChild(n);
    _scroll();
  } else if (ev.type === 'ai_gate') {
    AI._bumpUnread(ev);
    // 闸门会阻塞流水线：待决时自动弹出对话坞（唯一自动开坞的事件，
    // 必须有人拍板）。历史回放（backlog）不弹——刷新页面时旧闸门
    // 已决或已超时，弹坞纯属打扰
    if (ev.status === 'pending' && !AI.open && !ev.backlog) AI.toggle(true);
    const stTxt = {pending: t('ai.ev.gatePending', '等待你的决定'),
                   continue: t('ai.ev.gateContinue', '已确认继续'),
                   stop: t('ai.ev.gateStop', '已终止任务'),
                   timeout: t('ai.ev.gateTimeout', '超时自动继续')}[ev.status]
                || ev.status;
    // 同一样品同一阶段已有卡片：原地翻转状态（按钮移除、徽章更新），不再追加
    const old = [...document.querySelectorAll('.ai-gate')]
      .find(n => n.dataset.gate === `${ev.sample}|${ev.stage}`
                 && n.querySelector('.ai-gate-st.pending'));
    if (old && ev.status !== 'pending') {
      const b = old.querySelector('.ai-gate-st');
      b.textContent = stTxt;
      b.classList.replace('pending', 'done');
      const btns = old.querySelector('.ai-gate-btns');
      if (btns) btns.remove();
      _scroll();
      return;
    }
    const pending = ev.status === 'pending';
    const n = document.createElement('div');
    n.className = 'ai-msg ai-bot ai-gate';
    n.dataset.gate = `${ev.sample}|${ev.stage}`;
    n.innerHTML = `🚦 <b>${esc(t('ai.ev.gate', 'AI 决策闸门'))} · ${esc(ev.name || ev.stage)}</b>${esc(s)}`
      + `<span class="ai-gate-st ${pending ? 'pending' : 'done'}">${esc(stTxt)}</span>`
      + `<div class="ai-read-body">${aiFmt(ev.text || '')}</div>`
      + (pending ? `<div class="ai-gate-btns">
           <button class="btn small primary" onclick="AI.gate('${esc(ev.sample)}','${esc(ev.stage)}','continue')">${esc(t('ai.ev.gateBtnGo', '▶ 继续下一阶段'))}</button>
           <button class="btn small danger" onclick="AI.gate('${esc(ev.sample)}','${esc(ev.stage)}','stop')">${esc(t('ai.ev.gateBtnStop', '■ 终止任务'))}</button>
         </div>` : '');
    $('aiMsgs').appendChild(n);
    _scroll();
  } else if (ev.type === 'ai_action') {
    AI._bumpUnread(ev);
    // 执行单不打扰（不自动弹坞）——它不阻塞流水线，红点足够
    const stTxt = {pending: t('ai.act.pending', '待你批准'),
                   executing: t('ai.act.executing', '正在启动'),
                   launched: t('ai.act.launched', '✓ 已启动运行'),
                   rejected: t('ai.act.rejected', '已驳回'),
                   failed: t('ai.act.failed', '启动失败')}[ev.status]
                || ev.status;
    // 同一执行单已有卡片：原地翻转
    const oldA = [...document.querySelectorAll('.ai-action')]
      .find(n => n.dataset.aid === ev.aid
                 && n.querySelector('.ai-gate-st.pending'));
    if (oldA && ev.status !== 'pending') {
      const b = oldA.querySelector('.ai-gate-st');
      b.textContent = stTxt;
      b.classList.replace('pending', 'done');
      const btns = oldA.querySelector('.ai-gate-btns');
      if (btns) btns.remove();
      if (ev.error) {
        const e = document.createElement('div');
        e.className = 'ai-err';
        e.textContent = ev.error;
        oldA.appendChild(e);
      }
      _scroll();
      return;
    }
    const pendingA = ev.status === 'pending';
    const n = document.createElement('div');
    n.className = 'ai-msg ai-bot ai-action';
    n.dataset.aid = ev.aid;
    // 参数 diff：与上一轮基准不同的键高亮并加 *
    const base = ev.base_params || {};
    const paramsHtml = Object.entries(ev.params || {}).map(([k, v]) => {
      const changed = String(base[k]) !== String(v);
      return `<span class="ai-param ${changed ? 'changed' : ''}">`
        + `${esc(k)}=${esc(v)}${changed ? ' *' : ''}</span>`;
    }).join(' ') || `<span class="ai-hint">${esc(t('ai.act.sameParams', '（沿用上一轮参数）'))}</span>`;
    n.innerHTML = `🔧 <b>${esc(t('ai.act.title', 'AI 执行单'))} · ${esc(ev.title || ev.action)}</b>${esc(ev.sample ? ' · ' + ev.sample : '')}`
      + `<span class="ai-gate-st ${pendingA ? 'pending' : 'done'}">${esc(stTxt)}</span>`
      + (ev.stages && ev.stages.length
         ? `<div class="ai-act-stages">📑 ${esc((ev.stages || []).join(' → '))}</div>` : '')
      + (ev.action === 'download_runs'
         ? `<div class="ai-act-stages">⬇ ${esc((ev.runs || []).slice(0, 4).join(', '))}${((ev.runs || []).length > 4) ? ' …' : ''}（共 ${esc(String((ev.params || {}).runs_n || (ev.runs || []).length))} 条 run → 平台下载队列）</div>`
         : `<div class="ai-act-params">${paramsHtml}</div>`)
      + (ev.reason ? `<div class="ai-act-reason">${aiFmt(ev.reason)}</div>` : '')
      + (pendingA ? `<div class="ai-gate-btns">
           <button class="btn small primary" onclick="AI.action('${esc(ev.aid)}','approve')">${esc(t('ai.act.approve', '✓ 批准执行'))}</button>
           <button class="btn small danger" onclick="AI.action('${esc(ev.aid)}','reject')">${esc(t('ai.act.reject', '✗ 驳回'))}</button>
         </div>` : '');
    $('aiMsgs').appendChild(n);
    _scroll();
  }
};

AI.action = async function (aid, decision) {
  // 先本地置灰（乐观更新），结果状态由事件流原地翻转
  document.querySelectorAll('.ai-action').forEach(n => {
    if (n.dataset.aid === aid)
      n.querySelectorAll('.ai-gate-btns button').forEach(b => b.disabled = true);
  });
  try {
    const r = await fetch('/api/ai/actions/resolve', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({id: aid, decision})});
    if (!r.ok) {
      const msg = (await r.json()).error || `HTTP ${r.status}`;
      AI.sysNote(`⚠ ${msg}`);
      document.querySelectorAll('.ai-action').forEach(n => {
        if (n.dataset.aid === aid)
          n.querySelectorAll('.ai-gate-btns button')
            .forEach(b => b.disabled = false);
      });
    }
  } catch (e) { AI.sysNote(`⚠ ${e}`); }
};

AI.gate = async function (sample, stage, action) {
  // 同一样品同一阶段的旧闸门卡片立即置灰（新状态卡随后由事件流更新）
  document.querySelectorAll('.ai-gate').forEach(n => {
    if (n.dataset.gate === `${sample}|${stage}`)
      n.querySelectorAll('.ai-gate-btns button').forEach(b => b.disabled = true);
  });
  try {
    const r = await fetch('/api/ai/gate', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({sample, stage, action})});
    if (!r.ok) {
      const msg = (await r.json()).error || `HTTP ${r.status}`;
      AI.sysNote(`⚠ ${msg}`);
      document.querySelectorAll('.ai-gate').forEach(n => {
        if (n.dataset.gate === `${sample}|${stage}`)
          n.querySelectorAll('.ai-gate-btns button').forEach(b => b.disabled = false);
      });
    }
  } catch (e) { AI.sysNote(`⚠ ${e}`); }
};

AI.ensureCtx = async function () {
  try {
    const q = AI.sample ? `?sample=${encodeURIComponent(AI.sample)}` : '';
    const d = await (await fetch('/api/ai/context' + q)).json();
    AI.ready = !!d.ready;
    $('aiReadyDot').hidden = false;
    $('aiReadyDot').classList.toggle('ok', AI.ready);
    $('aiReadyDot').classList.toggle('no', !AI.ready);
    if (!AI._selFilled) {
      const sel = $('aiSampleSel');
      const opt0 = document.createElement('option');
      opt0.value = ''; opt0.textContent = t('ai.dock.any', '（不指定样品）');
      sel.appendChild(opt0);
      (d.samples || []).forEach(s => {
        const o = document.createElement('option');
        o.value = s.name;
        o.textContent = `${s.name} (${s.done}/${s.total})`;
        sel.appendChild(o);
      });
      // 默认选完成度最高的样品（对话最有的聊）
      const best = (d.samples || []).slice().sort((a, b) => b.done - a.done)[0];
      if (best && best.done > 0) { AI.sample = best.name; sel.value = best.name; }
      AI._selFilled = true;
    }
    if (!AI.ctxLoaded) {
      AI.ctxLoaded = true;
      if (!AI.ready) AI.notice(t('ai.dock.notice',
        'AI 助手尚未配置：请到 设置 → AI 助手 选择提供方并填写 API Key'
        + '（本地 Ollama 可完全离线使用）。'));
      // 有历史则恢复，没有才打招呼
      if (!AI._restoreHistory()) AI.greet();
    }
  } catch (e) {
    AI.notice(t('ai.dock.loadFail', '上下文加载失败') + ': ' + e);
  }
};

AI.onSampleChange = function (v) {
  AI.sample = v;
  if (v) AI.sysNote(`📎 ${t('ai.dock.sample', '样品')}: ${v}`);
};

/* ---------- 对话历史本地持久化 ----------
 * 刷新/换页不丢对话（localStorage，最多 30 条）。事件卡不持久化——
 * 它们由服务端环形缓冲在重连时重放。 */
AI._saveHistory = function () {
  try {
    localStorage.setItem(AI_MSGS_KEY,
      JSON.stringify(AI.messages.slice(-30)));
  } catch (e) { /* 存储不可用则静默 */ }
};

AI._restoreHistory = function () {
  try {
    const arr = JSON.parse(localStorage.getItem(AI_MSGS_KEY) || '[]');
    if (!Array.isArray(arr)) return false;
    for (const m of arr.slice(-30)) {
      if (!m || (m.role !== 'user' && m.role !== 'assistant')) continue;
      AI.messages.push({role: m.role, content: String(m.content || '')});
      const b = _bubble(m.role);
      if (m.role === 'user') b.textContent = String(m.content || '');
      else b.innerHTML = aiFmt(String(m.content || ''));
    }
    return AI.messages.length > 0;
  } catch (e) { return false; }
};

/* ---------- 消息渲染 ---------- */

// 极简 Markdown 渲染：先转义再放行少量标记（安全：不引入任何原始 HTML）
function aiFmt(s) {
  let h = esc(s);
  h = h.replace(/^### (.*)$/gm, '<b>$1</b>')
       .replace(/^## (.*)$/gm, '<b>$1</b>')
       .replace(/\*\*([^*\n]+)\*\*/g, '<b>$1</b>')
       .replace(/`([^`\n]+)`/g, '<code>$1</code>');
  // 连续的 "- " 行折叠成列表
  const lines = h.split('\n');
  let out = [], inList = false;
  for (const ln of lines) {
    const m = ln.match(/^\s*[-*]\s+(.*)$/);
    if (m) { if (!inList) { out.push('<ul>'); inList = true; } out.push(`<li>${m[1]}</li>`); }
    else {
      if (inList) { out.push('</ul>'); inList = false; }
      out.push(ln ? `<div>${ln}</div>` : '<div class="ai-br"></div>');
    }
  }
  if (inList) out.push('</ul>');
  return out.join('');
}

function _scroll() { const m = $('aiMsgs'); m.scrollTop = m.scrollHeight; }

function _bubble(role) {
  const wrap = document.createElement('div');
  wrap.className = `ai-msg ${role === 'user' ? 'ai-user' : 'ai-bot'}`;
  $('aiMsgs').appendChild(wrap);
  _scroll();
  return wrap;
}

AI.greet = function () {
  const name = AI.sample ? `"${AI.sample}"` : '';
  const b = _bubble('assistant');
  b.innerHTML = aiFmt(t('ai.dock.greeting',
    '你好！我是平台内置 AI 助手。我可以解读分析结果、解释参数、建议下一步。')
    + (AI.sample ? `\n\n📎 ${t('ai.dock.sample', '样品')}: **${AI.sample}**` : ''));
};

AI.notice = function (text) {
  const n = document.createElement('div');
  n.className = 'ai-notice';
  n.innerHTML = esc(text)
    + ` <a href="/settings#aiCard">${t('ai.dock.goSettings', '去配置')} ›</a>`;
  $('aiMsgs').appendChild(n);
  _scroll();
};

AI.sysNote = function (text) {
  const n = document.createElement('div');
  n.className = 'ai-sysnote';
  n.textContent = text;
  $('aiMsgs').appendChild(n);
  _scroll();
};

AI.clear = function () {
  AI.messages = [];
  try { localStorage.removeItem(AI_MSGS_KEY); } catch (e) {}
  $('aiMsgs').innerHTML = '';
  AI.greet();
};

/* ---------- 快捷指令 ---------- */

AI.chip = function (kind) {
  const q = {
    interpret: () => AI.sample
      ? t('ai.dock.qInterpret', '请解读这个样品的分析结果：检出了什么、各阶段质量如何、有什么要注意的。')
      : t('ai.dock.qNoSample', '平台目前有哪些样品，各自分析进行到哪一步？'),
    tasks: () => t('ai.dock.qTasks', '现在平台里有哪些分析任务？运行状态如何？'),
    next: () => AI.sample
      ? t('ai.dock.qNext', '基于这个样品目前的结果，下一步建议做什么？请给出可执行的具体步骤。')
      : t('ai.dock.qNextAny', '我刚接触这个平台，建议从哪里开始？'),
    params: () => AI.sample
      ? t('ai.dock.qParams', '如果要对这个样品重新分析（换更优参数），当前参数有哪些可以调整？请结合结果给出具体建议值。')
      : t('ai.dock.qParamsAny', '新建样品分析时，哪些参数最值得根据数据特点调整？'),
  }[kind];
  AI.send(q ? q() : '');
};

/* ---------- 发送与流式接收（可中断） ---------- */

// 发送 / 中断二合一：空闲时是发送键，生成中变成停止键
AI.sendOrStop = function () {
  if (AI.busy && AI._ctrl) { AI._ctrl.abort(); return; }
  AI.send();
};

AI.send = async function (preset) {
  if (AI.busy) return;
  const inp = $('aiInput');
  const text = (preset !== undefined ? preset : inp.value).trim();
  if (!text) return;
  if (preset === undefined) inp.value = '';
  AI.messages.push({role: 'user', content: text});
  AI._saveHistory();
  const ub = _bubble('user');
  ub.textContent = text;

  const bb = _bubble('assistant');
  bb.innerHTML = `<span class="ai-hint">${t('ai.dock.thinking', '思考中…')}</span>`;
  let acc = '';
  let toolsBox = null;
  let aborted = false;
  AI.busy = true;
  AI._ctrl = (typeof AbortController !== 'undefined') ? new AbortController() : null;
  const btn = $('aiSendBtn');
  btn.textContent = '■';
  btn.title = t('ai.dock.stop', '停止生成');
  btn.disabled = false;
  try {
    await AI._stream({messages: AI.messages, sample: AI.sample}, ev => {
      if (ev.type === 'delta') {
        acc += ev.content || '';
        bb.innerHTML = aiFmt(acc);
        _scroll();
      } else if (ev.type === 'tool') {
        if (!toolsBox) { toolsBox = document.createElement('div'); toolsBox.className = 'ai-tools'; bb.prepend(toolsBox); }
        const chip = document.createElement('span');
        chip.className = 'ai-toolchip';
        chip.textContent = `🔧 ${ev.name}(${Object.keys(ev.args || {}).length ? JSON.stringify(ev.args) : ''})`;
        toolsBox.appendChild(chip);
      } else if (ev.type === 'error') {
        const e = document.createElement('div');
        e.className = 'ai-err';
        e.innerHTML = esc(ev.message || 'error')
          + (String(ev.message || '').includes('设置')
             ? ` <a href="/settings#aiCard">${t('ai.dock.goSettings', '去配置')} ›</a>` : '');
        bb.appendChild(e);
      }
    }, AI._ctrl ? AI._ctrl.signal : undefined);
  } catch (e) {
    if (e && (e.name === 'AbortError' || String(e).includes('abort'))) {
      aborted = true;
      const s = document.createElement('div');
      s.className = 'ai-hint';
      s.textContent = `⏹ ${t('ai.dock.aborted', '已停止生成（保留以上内容）')}`;
      bb.appendChild(s);
    } else {
      const er = document.createElement('div');
      er.className = 'ai-err';
      er.textContent = t('ai.dock.busyFail', '请求失败') + ': ' + e;
      bb.appendChild(er);
    }
  } finally {
    AI.busy = false;
    AI._ctrl = null;
    btn.textContent = '➤';
    btn.title = t('ai.dock.send', '发送');
    if (!bb.querySelector('.ai-err') && !aborted && !acc) bb.innerHTML =
      `<span class="ai-hint">${t('ai.dock.empty2', '（无返回内容）')}</span>`;
    if (acc) { AI.messages.push({role: 'assistant', content: acc}); AI._saveHistory(); }
    _scroll();
  }
};

AI._stream = async function (body, onEvent, signal) {
  const r = await fetch('/api/ai/chat', {
    method: 'POST', headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body), signal});
  if (!r.ok) {
    let msg = `HTTP ${r.status}`;
    try { msg = (await r.json()).error || msg; } catch (e) {}
    onEvent({type: 'error', message: msg});
    return;
  }
  const reader = r.body.getReader();
  const dec = new TextDecoder();
  let buf = '';
  for (;;) {
    const {done, value} = await reader.read();
    if (done) break;
    buf += dec.decode(value, {stream: true});
    let i;
    while ((i = buf.indexOf('\n\n')) >= 0) {
      const frame = buf.slice(0, i).trim();
      buf = buf.slice(i + 2);
      if (!frame.startsWith('data:')) continue;
      try { onEvent(JSON.parse(frame.slice(5).trim())); } catch (e) {}
    }
  }
};

/* ---------- 启动（延迟到 load） ----------
 * 本脚本随 _nav.html 在页面顶部执行，而 $/esc/t 来自页面底部的
 * app.js / i18n.js —— 必须等 load 再绑定输入框与启动事件流，
 * 否则全局函数尚不存在（TDZ ReferenceError）。 */
(() => {
  const boot = () => {
    const inp = $('aiInput');
    if (!inp) return;                    // 本页无导航坞（理论不会发生）
    inp.addEventListener('keydown', e => {
      if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); AI.send(); }
    });
    inp.addEventListener('input', () => {
      inp.style.height = 'auto';
      inp.style.height = Math.min(inp.scrollHeight, 120) + 'px';
    });
    AI.startEvents();   // 常驻订阅主动播报（讲解/解读/闸门三层共用一条流）
  };
  if (document.readyState === 'complete') boot();
  else window.addEventListener('load', boot);
})();

window.AI = AI;

}  // end 防重复包含
