// ---------------- 任务轮询 + 结果预览 + 阶段卡片内嵌运行日志 + 上次结果恢复 ----------------
// 2026-09-13 自 app.js 拆出。普通 script 共享全局作用域：本文件在
// app.js 之后加载，可直接使用 $/esc/t/toast/jfetch 等核心工具；页面
// 内联脚本与其它 app-*.js 里的同名调用按全局解析（加载顺序见模板）。
// 注意：本块内引用的其它 app-*.js 函数均为调用期解析（事件/回调触发），
// 不存在加载期交叉依赖。
// ---------------- 任务轮询 + 结果预览 ----------------
let pollTimer = null;
const taskLogOpen = new Set();   // 展开日志的任务 id
const notifiedTasks = new Set(); // 已弹过结束通知的任务
let prevRunning = new Set();     // 上次轮询时运行中的任务

// 自适应轮询：有任务在跑时 2.5s（保持进度实时），全部空闲时降到 8s。
// 任务列表每 2.5s 全量回传，长时间挂着页面会产生大量无谓请求，
// 空闲降频可在不影响体验的前提下显著减轻服务端与浏览器负担。
const POLL_BUSY_MS = 2500;
const POLL_IDLE_MS = 8000;
let lastTaskList = [];
let connDown = false;

function startPolling() {
  if (pollTimer) clearTimeout(pollTimer);
  schedulePoll(0);
}

function schedulePoll(delay) {
  if (pollTimer) clearTimeout(pollTimer);
  pollTimer = setTimeout(async () => {
    // 递归调度不同于 setInterval：refreshTasks 抛异常会中断整条轮询链，
    // 这里必须兜住（原实现用 setInterval 不会有这个问题）
    try {
      await refreshTasks();
    } catch (e) {
      console.error('[poll] refreshTasks failed:', e);
    }
    const busy = lastTaskList.some(t => t.status === 'running') || connDown;
    schedulePoll(busy ? POLL_BUSY_MS : POLL_IDLE_MS);
  }, delay);
}

function taskStatusText(s) {
  return {running: t('c.st.running'), done: t('c.st.done'),
          failed: t('c.st.failed'), cancelled: t('c.st.cancelled')}[s] || s;
}

/* 任务结果预览面板（后端 result 字段为规范化预览结构） */
function taskResultHtml(res) {
  if (!res || !res.kind) return '';
  const chips = [];
  let links = [];
  if (res.kind === 'sample') {
    (res.stages || []).forEach(st => {
      if (st.status === 'done' && st.summary) {
        chips.push(`<span class="tr-chip" title="${esc(st.name)}">${esc(st.name)}：${esc(st.summary)}</span>`);
      }
    });
    if (res.report) {
      links.push(`<a class="btn small primary" href="${esc(res.report)}" target="_blank">${t('c.openReport')}</a>`);
      links.push(`<a class="btn small" href="/results" target="_blank">${t('nav.results')}</a>`);
    }
  } else if (res.kind === 'tool') {
    for (const [k, v] of Object.entries(res.stats || {})) {
      if (v === null || v === '') continue;
      const s = String(v);
      const isPath = /[\\/]/.test(s);
      const label = isPath ? s.split(/[\\/]/).pop() : s;
      chips.push(`<span class="tr-chip" title="${esc(k)}=${esc(s)}">${esc(k)}=${esc(label)}</span>`);
    }
    (res.files || []).forEach(f => {
      links.push(`<a class="btn small" href="/tool_runs/${esc(f.path)}" download title="${esc(f.path)}">📄 ${esc(f.path.split('/').pop())} <span class="hint">(${esc(f.size)})</span></a>`);
    });
  }
  if (!chips.length && !links.length) return '';
  return `<div class="task-result">
    <div class="tr-title">${t('c.resultTitle')}</div>
    ${chips.length ? `<div class="tr-stats">${chips.slice(0, 14).join('')}</div>` : ''}
    ${links.length ? `<div class="tr-links">${links.join('')}</div>` : ''}
  </div>`;
}

async function refreshTasks() {
  let tasks;
  try {
    const r = await fetch('/api/tasks');
    tasks = await r.json();
    lastTaskList = tasks || [];
    connFails = 0;
    connDown = false;
    setConnBanner(false);
    if (window.__rhOnTasks) { try { window.__rhOnTasks(tasks); } catch (e2) {} }
  } catch (e) {
    // 连接异常时保持短间隔重试，避免失联横幅要等满 3 次空闲周期才出现；
    // 仅真网络失败（TypeError）计入门槛，服务端 500/解析失败不谎报断连
    if (netFail(e)) {
      connDown = true;
      if (++connFails >= 3) setConnBanner(true);
    }
    return;
  }
  try {
  const box = $('tasks');
  // 页面没有「任务」汇总区（如工具页已移除底部任务面板）时跳过面板渲染；
  // 连接横幅、结束通知、工具卡内日志注入（injectStageLogs）仍照常执行。
  if (!box) {
    /* 无面板：仅保留通知与卡内日志注入 */
  } else if (!tasks.length) {
    box.innerHTML = `<p class="hint">${t('c.noTasks')}</p>`;
  } else {
  box.innerHTML = tasks.map(task => {
    const pct = Math.round((task.pct || 0) * 100);
    const stat = taskStatusText(task.status);
    const logHtml = taskLogOpen.has(task.id)
      ? `<pre data-autoscroll="${task.status === 'running' ? 1 : 0}">${esc((task.log || []).join('\n'))}</pre>` : '';
    const err = task.error ? `<div class="err">${esc(task.error)}</div>` : '';
    let timing = '';
    if (task.status === 'running') {
      const el = fmtDur(Date.now() / 1000 - (task.started || Date.now() / 1000));
      timing = `<span class="hint" style="margin-left:10px">${t('c.elapsed')} ${el}` +
        (task.eta != null ? ` · ${t('c.eta')} ~${fmtDur(task.eta)}` : '') + '</span>';
    } else if (task.finished && task.started) {
      timing = `<span class="hint" style="margin-left:10px">${t('c.used')} ${fmtDur(task.finished - task.started)}</span>`;
    }
    const preview = task.status === 'done' ? taskResultHtml(task.result) : '';
    return `<div class="task ${task.status}">
      <div style="display:flex;justify-content:space-between;align-items:center">
        <span class="tname">${esc(task.name)}</span>
        <span><span class="tstat ${task.status}">${stat}</span>${timing}
        ${task.status === 'running' ? `<button class="btn small danger" onclick="cancelTask('${task.id}')">${t('c.cancel')}</button>` : ''}
        </span>
      </div>
      <div class="tmsg">${esc(trTaskMsg(task.msg || ''))}</div>
      <div class="bar"><div style="width:${pct}%"></div></div>
      ${err}${preview}
      <button class="btn small" onclick="toggleLog('${task.id}')">${taskLogOpen.has(task.id) ? t('c.collapseLog') : t('c.expandLog')}</button>
      ${logHtml}
    </div>`;
  }).join('');
  }

  // 任务结束 → toast/系统通知 + 刷新管道视图
  const nowRunning = new Set(tasks.filter(x => x.status === 'running').map(x => x.id));
  for (const task of tasks) {
    if (task.status !== 'running' && prevRunning.has(task.id)
        && !notifiedTasks.has(task.id)) {
      notifiedTasks.add(task.id);
      if (task.status === 'done') {
        const hasReport = task.result && task.result.kind === 'sample' && task.result.report;
        toast(`${t('c.resultNotify')} · ${task.name}`,
              task.result && task.result.kind === 'sample' ? '' : (task.msg || ''),
              {actions: hasReport ? [
                {label: t('c.openReport'), primary: true,
                 onClick: () => window.open(task.result.report, '_blank')},
                {label: t('c.close'), onClick: () => {}},
              ] : [{label: t('c.close'), onClick: () => {}}]});
        sysNotify(t('c.resultNotify'), task.name);
      } else if (task.status === 'failed') {
        toast(`${t('c.failedNotify')} · ${task.name}`, task.error || '',
              {kind: 'failed', ttl: 12000, actions: [
                {label: t('c.close'), onClick: () => {}}]});
        sysNotify(t('c.failedNotify'), task.name);
      }
      loadSamples();
      if (curSample) selectSample(curSample);
      if (typeof loadToolRuns === 'function') loadToolRuns();
      break;   // 一次轮询只弹一条，避免轰炸
    }
  }
  prevRunning = nowRunning;
  } catch (e) { console.warn('任务面板渲染失败:', e); }
  injectStageLogs(tasks);
  autoscrollLogs();
}

// ---------------- 阶段卡片内嵌运行日志 ----------------
// 任务名按服务端语言生成；两种语言的标签都参与匹配
const STAGE_LABELS = {
  subsample:['预处理(子采样)', 'Prep (subsample)'],
  fastp:    ['⓪ Fastp 质控', '⓪ Fastp QC'],
  fq2fa:    ['⓪b 序列转换(FASTQ→FASTA)', '⓪b FASTQ→FASTA'],
  host:     ['① 宿主去除', '① Host removal'],
  kvsuite:  ['②b 已知病毒识别与定量', '②b Known virus ID & quantification'],
  assembly: ['③ 组装·分类·提取', '③ Assembly'],
  verify:   ['③b 候选序列验证', '③b Candidate verify'],
  consensus:['③c 共识序列与变异', '③c Consensus & variants'],
  hostana:  ['④ 宿主预测(ICTV)', '④ Host prediction (ICTV)'],
  orf:      ['⑥ ORF 预测', '⑥ ORF prediction'],
  orfa:     ['⑥b ORF 功能注释', '⑥b ORF annotation'],
  phylo:    ['⑦ 进化树与 SDT', '⑦ Phylogeny'],
  primer:   ['⑧ 引物设计', '⑧ Primer design'],
  gbdraw:   ['⑨ 基因组图', '⑨ Genome plots'],
  report:   ['⑩ 可视化报告', '⑩ Visual report'],
};

const TOOL_LABELS = {
  convert:  ['格式转换', 'Format convert'],   // 缺这条时运行中的转换卡日志区一直是空的
  fastp:    ['质控预处理', 'QC preprocess'],
  hostremoval: ['宿主去除与序列提取', 'Host removal'],
  hostpredict: ['宿主预测', 'Host prediction'],
  orf:      ['ORF 预测', 'ORF predict'],
  orfa:     ['功能注释', 'ORF annotate'],
  genoplot: ['基因组图谱', 'Genome plots'],
  primer:   ['引物设计', 'Primer design'],
  identify: ['病毒鉴定', 'Virus identify'],
  assemble: ['病毒组装', 'Virus assembly'],
  contigs:  ['contig分类', 'Contig classify'],
  verify:   ['候选序列验证', 'Candidate verify'],
  consensus: ['共识序列与变异', 'Consensus & variants'],
  // 一键流程两张卡此前没登记：injectStageLogs 找不到对应 toolrun 容器，
  // 卡片输出区一直是空的（2026-09-11 用户实测「一键分析（全流程）没有输出」）。
  kvsuite:  ['已知病毒识别与定量', 'Known virus suite'],
  kvchain:  ['病毒定量与共识·一键', 'KV chain'],
  virchain: ['病毒识别分类·一键', 'Virus chain'],
  ncbi:     ['NCBI下载', 'NCBI download'],
  gbdown:   ['GenBank下载', 'GenBank download'],
  gbimport: ['GenBank导入', 'GenBank import'],
  mirna:    ['miRNA 靶标预测', 'miRNA target'],
  // 以下条目 2026-10-26 补齐：/tasks 里这些任务名一直没有英文对照，
  // 英文界面下会原样漏出中文。
  prepchain: ['一键前处理', 'One-click prep'],
  pdprep:   ['进化动力学数据接入', 'Phylodynamics data prep'],
  phylodyn: ['⚡ 一键分析（定年 → 地理迁移 → 汇总溯源）',
             '⚡ One-click analysis (dating → migration → synthesis)'],
  phylogeo: ['系统地理分析', 'Phylogeography'],
  quicktree:['快速建树', 'Quick tree'],
  identity: ['NT+AA 同一性表', 'NT+AA identity table'],
  sdt:      ['SDT 精确分析', 'SDT exact analysis'],
  structcmp:['结构比较', 'Structure comparison'],
  align:    ['序列比对', 'Sequence alignment'],
};

// 数据库构建页：建库任务名 → 卡片内嵌日志容器（buildrun-<key>）。
// 两种语言的标签都参与匹配；refvirus/rvdb 用各自独有词区分。
const BUILD_TASK_LABELS = {
  taxonomy: ['Taxonomy', 'taxonomy'],
  host:     ['宿主库构建', 'Host DB build'],
  refvirus: ['RefSeq Viral'],
  rvdb:     ['RVDB 库构建', 'C-RVDB'],
  k2:       ['Kraken2 库转换', 'Kraken2 convert'],
};

function injectStageLogs(tasks) {
  snapshotLogScroll();
  document.querySelectorAll('.stage-log').forEach(el => { el.innerHTML = ''; });
  document.querySelectorAll('.toolrun').forEach(el => { el.innerHTML = ''; });
  document.querySelectorAll('.buildrun').forEach(el => { el.innerHTML = ''; });
  // 每个容器只吃**第一个**匹配任务：/api/tasks 是新→旧序（order.insert(0)），
  // 若不去重，同卡的旧任务（数组靠后）会后写覆盖新任务——实测连跑两次
  // 格式转换，第二次完成后输出区仍显示第一次的运行日志（2026-09-30 用户报告）。
  const claimed = new Set();
  for (const task of tasks) {
    if (!['running', 'done', 'failed'].includes(task.status)) continue;
    // 样品管道任务：按后端上报的当前阶段（task.stage）精确归属——
    // 日志只进"正在跑的这一张"阶段卡，不再整条链每张卡都重复一份。
    // stage 尚未上报（刚启动）时退化为名称匹配，且仅当名称只命中
    // 一个阶段时注入，避免多阶段链一次性铺满所有卡。
    // 另按样品名过滤：其它样品的管道任务不出现在当前样品的卡上。
    const nmStage = task.name || '';
    const sampleOk = !curSample || nmStage.includes(curSample);
    if (sampleOk) {
      let hitKeys = [];
      if (task.stage && STAGE_LABELS[task.stage]) {
        hitKeys = [task.stage];
      } else if (!task.stage) {
        hitKeys = Object.entries(STAGE_LABELS)
          .filter(([, labels]) => labels.some(label => nmStage.includes(label)))
          .map(([key]) => key);
        if (hitKeys.length > 1) hitKeys = [];   // 多阶段链未上报 → 等首个进度
      }
      for (const key of hitKeys) {
        const boxId = 'log-' + key;
        if (claimed.has(boxId)) continue;
        const box = $(boxId);
        if (box) { box.innerHTML = stageLogHtml(task); claimed.add(boxId); }
      }
    }
    // 专项分析工具任务（按工具标签匹配 → 卡内嵌 toolrun 容器）
    for (const [key, labels] of Object.entries(TOOL_LABELS)) {
      const nm = task.name || '';
      if (labels.some(label => nm.includes(label))
          && (nm.includes('工具·') || nm.includes('Tool·')
              || nm.includes('NCBI') || nm.includes('GenBank')
              || nm.includes('同属') || nm.includes('Synteny'))) {
        const boxId = 'toolrun-' + key;
        if (claimed.has(boxId)) continue;
        const box = $(boxId);
        if (box) { box.innerHTML = toolRunHtml(task, boxId); claimed.add(boxId); }
      }
    }
    // 建库任务 → 数据库构建页对应卡片内嵌日志（.buildrun 容器，按任务名匹配）
    for (const [key, labels] of Object.entries(BUILD_TASK_LABELS)) {
      const nm = task.name || '';
      if (labels.some(label => nm.includes(label))) {
        const boxId = 'buildrun-' + key;
        if (claimed.has(boxId)) continue;
        const box = $(boxId);
        if (box) { box.innerHTML = stageLogHtml(task); claimed.add(boxId); }
      }
    }
  }
  autoscrollLogs();
}

// 输出区标签状态（boxId → 'result' | 'log' | 'files'）。
// 输出区已改为纵向平铺（日志/结果/文件同时展示），不再渲染标签页按钮；
// 这里保留状态与 setOutTab 仅为向后兼容（防御性，不再被触发）。
const OUT_TAB_STATE = {};
const OUT_TASKS = {};

/* 工具产物文件名 → 语义化下载名 + 按钮配色。
   用于把原始英文文件名（run.log / contigs.filtered.fasta 等）翻译成
   「运行日志 / 组装结果」这类可读名称，并按主结果/辅助产物配色。 */
const DL_LABEL_SPECS = [
  // [文件名包含, 语义化名, 按钮色]（先匹配先得，结果类放前）
  // —— 组装 / 分类主结果 ——
  ['contigs.filtered.fasta', '组装结果 (contigs)', 'primary'],
  ['transcripts.fasta', '组装结果 (transcripts)', 'primary'],
  ['viral_contigs.fasta', '病毒序列 (viral contigs)', 'primary'],
  ['virus_classification.tsv', '病毒分类结果', 'primary'],
  ['viral_ids.tsv', '分类序列 ID', 'primary'],
  ['viral_sequences.', '病毒序列', 'primary'],
  ['contigs.fasta', '组装结果 (contigs)', 'primary'],
  ['scaffolds.fasta', '组装结果 (scaffolds)', 'primary'],
  ['pyrodigal.faa', 'ORF 蛋白', 'primary'],
  ['pyrodigal.ffn', 'ORF 核苷酸', 'primary'],
  ['host_prediction.tsv', '宿主预测结果', 'primary'],
  // —— 报告 / 汇总 ——
  ['report.html', '分析报告', 'green'],
  ['kreport', '分类报告 (kreport)', 'green'],
  ['keep_R1.fastq.gz', '保留 reads R1', 'green'],
  ['keep_R2.fastq.gz', '保留 reads R2', 'green'],
  ['kept_R1.fastq.gz', '保留 reads R1', 'green'],
  ['kept_R2.fastq.gz', '保留 reads R2', 'green'],
  ['taxburst', 'Krona 旭日图', 'green'],
  // —— 图 / 结构文件 ——
  ['.svg', '基因组图 (SVG)', 'green'],
  ['.gfa', '组装图 (GFA)', 'amber'],
  ['before_rr.fasta', '纠错前序列', 'amber'],
  ['fastp.', '质控报告', 'amber'],
  // —— 日志 / 配置 ——
  ['run.log', '运行日志', 'gray'],
  ['spades.log', 'SPAdes 日志', 'gray'],
  ['warnings.log', '警告日志', 'gray'],
  ['spades.sh', 'SPAdes 脚本', 'gray'],
  ['spades.yaml', 'SPAdes 配置', 'gray'],
  ['run_spades.sh', 'SPAdes 脚本', 'gray'],
  ['run_spades.yaml', 'SPAdes 配置', 'gray'],
  ['params.txt', '参数表', 'gray'],
  ['input_dataset.yaml', '输入数据描述', 'gray'],
  ['dataset.info', '数据描述', 'gray'],
  ['.json', '结果数据 (JSON)', 'gray'],
  ['summary.', '阶段汇总', 'gray']
];
function dlLabel(base, path) {
  const low = String(base || path || '').toLowerCase();
  for (const [key, label, color] of DL_LABEL_SPECS) {
    if (low.includes(key)) return { label, color };
  }
  return { label: base, color: '' };   // 未知 → 原始文件名
}

/* 输出目录行：结果产出位置随结果一起明示（完整路径 + 复制 + 打开）。
   路径经 data-* 属性传递再由 dataset 读回：Windows 路径里的 \t \r 这类
   序列若直接拼进 onclick 的 JS 字符串字面量会被当成转义符。 */
function outDirLine(outDir) {
  if (!outDir) return '';
  const od = String(outDir);
  return `<div style="display:flex;gap:6px;align-items:center;flex-wrap:wrap;margin:8px 0 2px;font-size:12.5px">
    <span style="font-weight:700">${t('c.outdir', '输出目录')}</span>
    <code class="mono" style="user-select:all;word-break:break-all" title="${esc(od)}">${esc(od)}</code>
    <button class="btn small" style="padding:1px 8px" data-d="${esc(od)}" onclick="copyText(this.dataset.d)" title="${t('c.copyPath', '复制路径')}">⧉</button>
    <button class="btn small" style="padding:1px 8px" data-d="${esc(od)}" onclick="openOutDir(this.dataset.d)" title="${t('c.openDir', '打开目录')}">📂</button>
  </div>`;
}

/* 在资源管理器中打开平台内的输出目录（/api/open_dir 做平台内校验） */
async function openOutDir(p) {
  try {
    const r = await fetch('/api/open_dir', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ path: p }) });
    if (!r.ok) {
      const d = await r.json().catch(() => ({}));
      alert((d.error || t('ui.common.openFail')) + ' (HTTP ' + r.status + ')');
    }
  } catch (e) { alert(t('ui.common.openFail') + ': ' + e); }
}

function outFilesHtml(task) {
  const files = (task.result && task.result.files) || [];
  if (!files.length) {
    return `<p class="hint" style="margin:6px 0">${
      task.status === 'running' ? t('c.filesRunning') : t('c.noFiles')}</p>`;
  }
  // 链接必须含运行名（/tool_runs/<run>/<file>），task.result.run 由后端注入；
  // download 属性带运行名前缀：不同运行的同名产物（viral_ids.tsv 等）
  // 下载到本地不再互相覆盖。
  const run = (task.result && task.result.run) || '';
  const btns = files.map(f => {
    const normPath = String(f.path).replace(/\\/g, '/');
    const segs = normPath.split('/').map(encodeURIComponent).join('/');
    const base = String(f.path).split(/[\\/]/).pop();
    const dl = run ? `${run}_${base}` : base;
    const href = run ? `/tool_runs/${encodeURIComponent(run)}/${segs}` : '#';
    const { label, color } = dlLabel(base, f.path);
    const cls = color ? ` dlbtn ${color}` : ' dlbtn';
    const low = base.toLowerCase();
    // 图件（SVG/PNG）额外给一个「整理」入口：同源 URL 直接喂给内置 svguest。
    const edit = (run && /\.(svg|png)$/i.test(base))
      ? ' ' + svgEditLink(href, '整理', 'dlbtn') : '';
    // 基因组图 (SVG)：直接内联预览，不再只给下载按钮。
    if (low.endsWith('.svg') && run) {
      const img = `<img src="/tool_runs/${encodeURIComponent(run)}/${segs}" alt="${esc(base)}" style="display:block;width:100%;height:220px;object-fit:contain;border:1px solid var(--line-100);border-radius:8px;margin:8px 0;background:#fff" loading="lazy">`;
      return `<div class="svg-prev" style="border-radius:8px;padding:8px;background:#fff">${img}<a class="${cls.trim()}" href="${href}" download="${esc(dl)}" title="${esc(f.path)}【${esc(f.size)}】">⬇ 下载 ${esc(label)}</a>${edit}</div>`;
    }
    return `<a class="${cls.trim()}" href="${href}" download="${esc(dl)}" title="${esc(f.path)}【${esc(f.size)}】">⬇ ${esc(label)}</a>${edit}`;
  }).join('');
  return `<div class="dlbtn-row">${btns}</div>`;
}

/* 工具运行的关键数字 chips（仅标量；路径取文件名；run 名已在状态徽标里，不重复） */
const TOOL_STAT_LABELS = {
  n_classified: ['已分类序列', 'Classified seqs'],
  n_extracted:  ['提取病毒序列', 'Extracted viral seqs'],
  n_contigs:    ['过滤后 contigs', 'Contigs kept'],
  n_viral:      ['病毒 contigs', 'Viral contigs'],
  total_bp:     ['总长度(bp)', 'Total bp'],
};

function toolStatChips(res) {
  if (!res || res.kind !== 'tool') return '';
  const zh = (typeof VP_LANG === 'undefined' || VP_LANG !== 'en');
  const chips = [];
  for (const [k, v] of Object.entries(res.stats || {})) {
    if (v === null || v === '' || k === 'run' || k === 'out_dir') continue;
    const lab = (TOOL_STAT_LABELS[k] || [k, k])[zh ? 0 : 1];
    const s = String(v);
    const shown = /[\\/]/.test(s) ? s.split(/[\\/]/).pop() : s;
    chips.push(`<span class="tr-chip" title="${esc(k)}=${esc(s)}">${esc(lab)}=${esc(shown)}</span>`);
  }
  return chips.length ? `<div class="tr-stats">${chips.slice(0, 14).join('')}</div>` : '';
}

function toolRunHtml(task, boxId) {
  boxId = boxId || '';
  if (boxId) { OUT_TASKS[boxId] = task; if (!OUT_TAB_STATE[boxId]) OUT_TAB_STATE[boxId] = 'result'; }
  const pct = Math.round((task.pct || 0) * 100);
  const stat = taskStatusText(task.status);
  const err = task.error ? `<div class="err" style="margin:6px 0 0">${esc(task.error)}</div>` : '';
  const chips = task.status === 'done' ? toolStatChips(task.result) : '';
  // 输出目录：任务创建即带（运行中也能看到将写到哪），完成后结果里再确认一次
  const outdir = task.out_dir || (task.result && task.result.out_dir) || '';
  const lines = esc((task.log || []).join('\n'));
  // 纵向平铺：状态/进度 → 输出目录 → 关键数字 → 运行日志 → Downloads（下载）。
  // 产物文件统一在 Downloads 一处列出（分类报告面板内不再重复）。
  return `<div style="margin:10px 0 2px">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:4px">
      <span class="sbadge ${task.status}">${stat} · ${esc(task.name)}</span>
      ${task.status === 'running' ? `<button class="btn small danger" onclick="cancelTask('${task.id}')">${t('c.cancel')}</button>` : ''}
    </div>
    <div class="bar"><i style="width:${pct}%"></i></div>
    ${err}${outdir ? outDirLine(outdir) : ''}${chips}
    <div style="font-size:12.5px;font-weight:700;margin:10px 0 4px">${t('c.runlog')}</div>
    <pre class="logbox" style="max-height:180px;overflow-y:auto" data-autoscroll="${task.status === 'running' ? 1 : 0}">${lines || t('c.noLog')}</pre>
    <div style="font-size:12.5px;font-weight:700;margin:10px 0 4px">${t('c.downloads')}</div>
    ${outFilesHtml(task)}
  </div>`;
}

function setOutTab(boxId, tab) {
  OUT_TAB_STATE[boxId] = tab;
  const task = OUT_TASKS[boxId];
  const box = $(boxId);
  if (task && box) box.innerHTML = toolRunHtml(task, boxId);
}
document.addEventListener('click', ev => {
  const btn = ev.target.closest?.('.out-tab');
  if (!btn) return;
  const holder = btn.closest('.toolrun');
  if (holder) setOutTab(holder.id, btn.dataset.tab);
});

/* ============================================================
 * 全局「上次运行结果」恢复（页面加载 / 模块切换时自动执行）
 * ------------------------------------------------------------
 * 工具卡的结果原本只在任务完成那一刻渲染进 DOM：跳到别的页面再回来，
 * 结果区就空了，只能重跑。这里在页面加载/模块切换后，从磁盘上的
 * tool_runs/（/api/tool/runs，按 mtime 倒序）找出每张**可见**卡片最近一次
 * 运行并自动恢复：
 *   - 通用卡：渲染「上次运行 · <run>」+ 产物下载块（outFilesHtml 复用）；
 *   - 富结果卡（SDT / 同一性）：调用其专用加载函数，恢复交互热图与矩阵表。
 * 有任务在跑、或卡片已有结果时不覆盖（交给轮询/用户操作）。
 * ============================================================ */
const CARD_RUN_PREFIX = {
  convert: 'convert', fastp: 'fastp', identify: 'identify', assemble: 'assemble',
  contigs: 'contigs', verify: 'verify', virchain: 'virchain',
  consensus: 'consensus', kvsuite: 'kvsuite', kvchain: 'kvchain',
  align: 'align', treebuild: 'quicktree', sdt: 'sdt',
  orf: 'orf', orfa: 'orfa', genoplot: 'genoplot', primer: 'primer',
  hostremoval: 'hostremoval'
};
/* 富结果卡的恢复函数（键 = 卡片 id 去掉 t- 前缀） */
const TOOL_RESTORE = {
  sdt: run => { sdRun = run; sdtExactLoad(); },
  identity: run => { if (typeof identityLoad === 'function') identityLoad(run); },
  mirna: run => { if (typeof loadMirnaResult === 'function') loadMirnaResult(run); }
};
/* 卡片「已有结果」探测：避免覆盖正在跑 / 刚跑完的结果 */
function cardResultFilled(key) {
  const probe = { sdt: 'sdExact', identity: 'idtResult' }[key];
  if (probe) {
    const el = $(probe);
    return !!(el && el.innerHTML.trim() && el.offsetParent !== null);
  }
  const box = $('toolrun-' + key);
  return !!(box && box.innerHTML.trim());
}
function runRestoreHtml(run) {
  const fake = { status: 'done', out_dir: run.out_dir || '',
                 result: { run: run.name, files: run.files || [] } };
  return `<div style="margin:10px 0 2px">
      <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:4px">
        <span class="sbadge done">${t('c.restoredRun', '上次运行')} · ${esc(run.name)}</span>
        <span class="hint">${t('c.restoredHint', '结果已从磁盘恢复，无需重跑')}</span>
      </div>
      ${outDirLine(run.out_dir)}
      <div style="font-size:12.5px;font-weight:700;margin:10px 0 4px">${t('c.downloads')}</div>
      ${outFilesHtml(fake)}
    </div>`;
}
let _restoreBusy = false;
async function restoreLastResults() {
  if (_restoreBusy) return;
  const cards = [...document.querySelectorAll('section.card[id^="t-"]')]
    .filter(c => c.offsetParent !== null);           // /tools 一次只显示一个二级模块
  if (!cards.length) return;
  _restoreBusy = true;
  try {
    const tasks = await (await fetch('/api/tasks?logs=0')).json();
    if ((tasks || []).some(x => x.status === 'running')) return;   // 有任务在跑 → 交给轮询
    const runs = await (await fetch('/api/tool/runs')).json();
    const newest = {};
    for (const r of (runs || [])) {
      const m = String(r.name).match(/^(.+?)_\d{8}_\d{6}/);
      const p = m ? m[1] : String(r.name);
      if (!newest[p]) newest[p] = r;                 // 接口已按 mtime 倒序
    }
    for (const card of cards) {
      const key = card.id.slice(2);
      const prefix = CARD_RUN_PREFIX[key];
      if (!prefix || cardResultFilled(key)) continue;
      const run = newest[prefix];
      if (!run) continue;
      if (TOOL_RESTORE[key]) {
        try { TOOL_RESTORE[key](run.name); } catch (e) { /* 恢复失败不打扰用户 */ }
      } else {
        const box = $('toolrun-' + key);
        if (box) box.innerHTML = runRestoreHtml(run);
      }
    }
  } catch (e) { /* 静默：恢复失败不影响正常使用 */ }
  finally { _restoreBusy = false; }
}
/* 页面加载完成 & /tools 模块切换（hash 变化，不重载页面）后各触发一次 */
window.addEventListener('load', () => setTimeout(restoreLastResults, 900));
window.addEventListener('hashchange', () => setTimeout(restoreLastResults, 600));

function stageLogHtml(task) {
  const stat = taskStatusText(task.status);
  const lines = esc((task.log || []).join('\n'));
  const closed = logClosedSet().has(task.id);
  return `
    <div style="display:flex;justify-content:space-between;align-items:center;margin:8px 0 4px">
      <span class="sbadge ${task.status}">${stat} · ${t('c.liveLog')}</span>
      <span>
        <button class="btn small" onclick="logToggle('${task.id}')" title="折叠/展开日志">${closed ? '▸' : '▾'} ${t('c.logBtn', '日志')}</button>
        <button class="btn small" onclick="copyStageLog(this)">${t('c.copyLog')}</button>
        ${task.status === 'running' ? `<button class="btn small danger" onclick="cancelTask('${task.id}')">${t('c.cancel')}</button>` : ''}
      </span>
    </div>
    <pre class="logbox" style="margin:0;${closed ? 'display:none' : ''}" data-tid="${esc(task.id)}" data-autoscroll="${task.status === 'running' ? 1 : 0}">${lines || t('c.noLog')}</pre>`;
}

function copyStageLog(btn) {
  const pre = btn.closest('.stage-log')?.querySelector('pre.logbox');
  if (!pre) return;
  const text = pre.textContent;
  (navigator.clipboard?.writeText(text) ?? Promise.reject())
    .then(() => { btn.textContent = t('c.copied'); setTimeout(() => btn.textContent = t('c.copyLog'), 1500); })
    .catch(() => {
      const ta = document.createElement('textarea');
      ta.value = text; document.body.appendChild(ta); ta.select();
      document.execCommand('copy'); ta.remove();
      btn.textContent = t('c.copied'); setTimeout(() => btn.textContent = t('c.copyLog'), 1500);
    });
}

/* 日志滚动策略：重绘前记录每个日志框 scrollTop 与是否贴底；重绘后——
   运行中且原本贴底 → 跟随新日志；其余（含已结束任务、用户上翻阅读中）
   → 精确恢复原位置，不再被 2.5s 轮询重绘拉回顶部。 */
const _logScrollMemo = new WeakMap();
function snapshotLogScroll() {
  document.querySelectorAll('.toolrun, .stage-log, .buildrun').forEach(box => {
    const st = [...box.querySelectorAll('pre.logbox')].map(pre => ({
      top: pre.scrollTop,
      atBottom: pre.scrollHeight - pre.scrollTop - pre.clientHeight < 16
    }));
    _logScrollMemo.set(box, st);
  });
}
function restoreLogScroll() {
  document.querySelectorAll('.toolrun, .stage-log, .buildrun').forEach(box => {
    const memo = _logScrollMemo.get(box) || [];
    [...box.querySelectorAll('pre.logbox')].forEach((pre, i) => {
      const st = memo[i];
      if (pre.dataset.autoscroll === '1' && (!st || st.atBottom)) {
        pre.scrollTop = pre.scrollHeight;
      } else if (st) {
        pre.scrollTop = st.top;
      }
    });
  });
}
function autoscrollLogs() { restoreLogScroll(); }

function toggleLog(tid) {
  taskLogOpen.has(tid) ? taskLogOpen.delete(tid) : taskLogOpen.add(tid);
  refreshTasks();
}
function logClosedSet() {
  try { return new Set(JSON.parse(sessionStorage.getItem('vp_logclosed') || '[]')); }
  catch (e) { return new Set(); }
}
function logToggle(tid) {
  const s = logClosedSet();
  s.has(tid) ? s.delete(tid) : s.add(tid);
  try { sessionStorage.setItem('vp_logclosed', JSON.stringify([...s])); } catch (e) {}
  refreshTasks();
}

async function cancelTask(tid) {
  await fetch('/api/task/' + tid + '/cancel', {method: 'POST'});
  refreshTasks();
}

