// ---------------- 模块历史运行（折叠 · 上下游衔接 · 删除） ----------------
// 2026-09-13 自 app.js 拆出。普通 script 共享全局作用域：本文件在
// app.js 之后加载，可直接使用 $/esc/t/toast/jfetch 等核心工具；页面
// 内联脚本与其它 app-*.js 里的同名调用按全局解析（加载顺序见模板）。
// 注意：本块内引用的其它 app-*.js 函数均为调用期解析（事件/回调触发），
// 不存在加载期交叉依赖。
// ---------------- 模块历史运行（折叠 · 上下游衔接 · 删除） ----------------
/* 每个工具卡/独立模块页的持久历史：数据来自 /api/tool/runs（tool_runs/ 目录，
   重启不丢）。条目可展开文件下载、一键复用为下游模块输入（跨页暂存回填）、
   可删除。默认折叠，避免过多。 */
const RH_PAGES = {
  f_r1: '/tools#t-fastp', f_r2: '/tools#t-fastp',
  hr_r1: '/hostremoval', hr_r2: '/hostremoval',
  i_input: '/tools#t-identify', i_input2: '/tools#t-identify',
  a_r1: '/tools#t-assemble', a_r2: '/tools#t-assemble',
  c_fa: '/tools#t-contigs',
  of_fa: '/orf', oa_fa: '/annotation', oa_run: '/annotation',
  gp_fa: '/genome', gp_ann: '/genome', pr_fa: '/primer',
  al_fa: '/tools#t-align', qt_fa: '/tools#t-treebuild',
  tv_file: '/tools#t-treebuild', sd_fa: '/tools#t-sdt'
};

function rhPick(files, sub, excl) {
  return ((files || []).filter(f => (f.path || '').includes(sub)
                            && !(excl && (f.path || '').includes(excl)))[0] || {}).path;
}
function rhFqPair(files) {
  const fqs = (files || []).map(f => f.path).filter(p => /\.fastq\.gz$/.test(p));
  if (!fqs.length) return null;
  const num = p => { const b = p.replace(/\.fastq\.gz$/, '').split(/[\/]/).pop();
                     return b.endsWith('_1') || b.endsWith('.1') || b.endsWith('_R1') ? 1
                          : b.endsWith('_2') || b.endsWith('.2') || b.endsWith('_R2') ? 2 : 0; };
  const r1 = fqs.find(p => num(p) === 1), r2 = fqs.find(p => num(p) === 2);
  if (r1 && r2) return [r1, r2];
  return [fqs[0], null];                       // 单端
}
const RH_DEFS = [
  { box: 'rh-convert', prefixes: ['convert_'],
    chains: files => {
      const out = [];
      const fa = rhPick(files, '.fa.gz') || rhPick(files, '.fasta.gz');
      const fq = rhFqPair(files);
      if (fq) out.push(
        { label: '→ 质控', page: RH_PAGES.f_r1,
          fills: { f_r1: fq[0], ...(fq[1] ? { f_r2: fq[1] } : {}) } },
        { label: '→ 宿主去除', page: RH_PAGES.hr_r1,
          fills: { hr_r1: fq[0], ...(fq[1] ? { hr_r2: fq[1] } : {}) } },
        { label: '→ 病毒识别', page: RH_PAGES.i_input,
          fills: { i_input: fq[0], ...(fq[1] ? { i_input2: fq[1] } : {}) } });
      if (fa) out.push(
        { label: '→ 识别(FASTA)', page: RH_PAGES.i_input,
          fills: { i_input: fa, i_type: 'fasta' } },
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } });
      return out;
    } },
  { box: 'rh-fastp', prefixes: ['fastp_'],
    chains: files => {
      const p = f => rhPick(files, f);
      const r1 = p('fastp_R1'), r2 = p('fastp_R2');
      if (!r1) return [];
      return [
        { label: '→ 宿主去除', page: RH_PAGES.hr_r1,
          fills: { hr_r1: r1, hr_r2: r2 || '' } },
        { label: '→ 病毒识别', page: RH_PAGES.i_input,
          fills: { i_input: r1, i_input2: r2 || '' } },
        { label: '→ 组装', page: RH_PAGES.a_r1,
          fills: { a_r1: r1, a_r2: r2 || '' } }];
    } },
  { box: 'rh-hostremoval', prefixes: ['hostremoval_'],
    chains: files => {
      const r1 = rhPick(files, 'kept_R1'), r2 = rhPick(files, 'kept_R2');
      if (!r1) return [];
      return [
        { label: '→ 病毒识别', page: RH_PAGES.i_input,
          fills: { i_input: r1, i_input2: r2 || '' } },
        { label: '→ 组装', page: RH_PAGES.a_r1,
          fills: { a_r1: r1, a_r2: r2 || '' } },
        { label: '→ 质控', page: RH_PAGES.f_r1,
          fills: { f_r1: r1, f_r2: r2 || '' } }];
    } },
  { box: 'rh-identify', prefixes: ['identify_'],
    chains: files => {
      const fa = rhPick(files, 'viral_sequences');
      if (!fa) return [];
      return [
        { label: '→ 组装', page: RH_PAGES.a_r1, fills: { a_r1: fa } },
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } }];
    } },
  { box: 'rh-assemble', prefixes: ['assemble_'],
    chains: files => {
      const fa = rhPick(files, 'contigs.filtered.fasta');
      if (!fa) return [];
      return [
        { label: '→ contigs 分类', page: RH_PAGES.c_fa, fills: { c_fa: fa } },
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } },
        { label: '→ 图谱', page: RH_PAGES.gp_fa, fills: { gp_fa: fa } },
        { label: '→ 引物', page: RH_PAGES.pr_fa, fills: { pr_fa: fa } }];
    } },
  { box: 'rh-contigs', prefixes: ['contigs_'],
    chains: files => {
      const fa = rhPick(files, 'viral_contigs.fasta');
      const out = [];
      if (fa) out.push(
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 功能注释', page: RH_PAGES.oa_fa, fills: { oa_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } },
        { label: '→ 图谱', page: RH_PAGES.gp_fa, fills: { gp_fa: fa } },
        { label: '→ 引物', page: RH_PAGES.pr_fa, fills: { pr_fa: fa } });
      /* 「→ 宿主预测」链随 /hostpredict 孤儿页删除（2026-09-27）一并移除：
         宿主预测现在由 contigs 分类后自动运行（tools.html autoHostIfMissing），
         不再需要手动载入分类表。 */
      return out;
    } },
  { box: 'rh-orf', prefixes: ['orf_'],
    chains: (files, run) => {
      const out = [];
      const faa = rhPick(files, 'pyrodigal.faa');
      const nt = rhPick(files, 'pyrodigal.ffn');
      if (faa) out.push(
        { label: '→ 比对(蛋白)', page: RH_PAGES.al_fa, fills: { al_fa: faa } },
        { label: '→ SDT(AA)', page: RH_PAGES.sd_fa, fills: { sd_fa: faa } });
      if (nt) out.push(
        { label: '→ SDT(NT)', page: RH_PAGES.sd_fa, fills: { sd_fa: nt } });
      out.push({ label: '→ 功能注释(该运行)', page: '/annotation',
                 fills: { oa_run: run.name } });
      return out;
    } },
  { box: 'rh-orfa', prefixes: ['orfa_'],
    chains: (files, run) => {
      const gff = rhPick(files, 'orf_annotation.gff3');
      const fa = rhPick(files, 'viral_contigs.fasta');
      const out = [];
      if (gff && fa) out.push(
        { label: '→ 基因组图谱(注释)', page: '/genome',
          fills: { gp_ann: gff, gp_fa: fa } });
      const realRun = (run.name || '').replace(/^orfa_/, 'orf_');
      out.push({ label: '→ 功能注释(orf 运行)', page: '/annotation',
                 fills: { oa_run: realRun } });
      return out;
    } },
  { box: 'rh-verify', prefixes: ['verify_'],
    chains: files => {
      // 验证通过/候选的序列 → 下游注释分析（历史运行此前没登记，这个盒子一直是空的）
      const fa = rhPick(files, 'virus_input.fasta')
        || rhPick(files, 'viroid_input.fasta');
      if (!fa) return [];
      return [
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 功能注释', page: RH_PAGES.oa_fa, fills: { oa_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } }];
    } },
  { box: 'rh-genoplot', prefixes: ['genoplot_'], chains: () => [] },
  { box: 'rh-primer', prefixes: ['primer_'], chains: () => [] },
  { box: 'rh-align', prefixes: ['align_'],
    chains: files => {
      const aln = rhPick(files, 'aln.trim.fasta') || rhPick(files, 'aln.fasta');
      if (!aln) return [];
      return [
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: aln } },
        { label: '→ SDT(已比对)', page: RH_PAGES.sd_fa,
          fills: { sd_fa: aln, sd_aligned: true } }];
    } },
  { box: 'rh-treebuild', prefixes: ['structcmp_', 'quicktree_'],
    chains: files => {
      const out = [];
      const aln = rhPick(files, 'aln.fasta');
      const nwk = rhPick(files, 'nj.nwk') || rhPick(files, 'tree.nwk');
      if (aln) out.push(
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: aln } },
        { label: '→ SDT(已比对)', page: RH_PAGES.sd_fa,
          fills: { sd_fa: aln, sd_aligned: true } });
      if (nwk) out.push(
        { label: '→ 树查看', page: RH_PAGES.tv_file, fills: { tv_file: nwk } });
      return out;
    } },
  { box: 'rh-sdt', prefixes: ['sdt_', 'identity_'],
    chains: files => {
      const out = [];
      const csv = rhPick(files, 'sdt_matrix.csv');
      if (csv) out.push({ label: '⬇ 矩阵 CSV', page: '', fills: {}, dl: csv });
      return out;
    } }
];

let RH_DATA = null;
const RH_OPEN = new Set();

async function rhRefresh() {
  try {
    const r = await fetch('/api/tool/runs');
    RH_DATA = await r.json();
  } catch (e) { RH_DATA = []; }
  for (const def of RH_DEFS) {
    const box = $(def.box);
    if (!box) continue;
    const runs = RH_DATA.filter(x => def.prefixes.some(p => x.name.startsWith(p)));
    const n = runs.length;
    const head = box.querySelector('.rh-head');
    if (head) head.querySelector('.rh-count').textContent = n;
    const body = box.querySelector('.rh-body');
    if (!body) continue;
    if (!n) { body.innerHTML = '<p class="hint">（暂无历史运行）</p>'; continue; }
    body.innerHTML = runs.map(run => {
      const chains = (def.chains || (() => []))(run.files || [], run) || [];
      const chainBtns = chains.map((c, i) => c.dl
        ? `<a class="btn small" href="/tool_runs/${encodeURIComponent(run.name)}/${c.dl}" download>${esc(c.label)}</a> `
        : `<button class="btn small" onclick="rhGo(this, '${esc(run.name)}', '${esc(def.box)}', ${i})">${esc(c.label)}</button> `).join('');
      const open = RH_OPEN.has(run.name);
      // 关键产物：只展示主结果/报告/序列（primary/green），其余收进“全部文件”
      const keyed = (run.files || []).map(f => {
        const base = String(f.path).split('/').pop();
        const { label, color } = dlLabel(base, f.path);
        return { f, base, label, color };
      });
      const keys = keyed.filter(x => x.color === 'primary' || x.color === 'green');
      const rest = keyed.filter(x => x.color !== 'primary' && x.color !== 'green');
      const keyHtml = keys.length ? keys.map(x => {
        const cls = x.color ? ` dlbtn ${x.color}` : ' dlbtn';
        const low = x.base.toLowerCase();
        // 基因组图 (SVG)：直接内联预览
        if (low.endsWith('.svg')) {
          const normPath = String(x.f.path).replace(/\\/g, '/');
          const segs = normPath.split('/').map(encodeURIComponent).join('/');
          return `<div style="background:#fff;border-radius:8px;margin:6px 0;padding:8px"><img src="/tool_runs/${encodeURIComponent(run.name)}/${segs}" alt="${esc(x.base)}" style="display:block;width:100%;height:220px;object-fit:contain;border:1px solid var(--line-100);border-radius:8px;background:#fff" loading="lazy"><a class="${cls.trim()}" href="/tool_runs/${encodeURIComponent(run.name)}/${encodeURIComponent(normPath)}" title="${esc(x.f.path)}【${esc(x.f.size)}】">⬇ 下载 ${esc(x.label)}</a></div>`;
        }
        return `<a class="${cls.trim()}" href="/tool_runs/${encodeURIComponent(run.name)}/${encodeURIComponent(String(x.f.path).replace(/\\/g, '/'))}" title="${esc(x.f.path)}【${esc(x.f.size)}】">⬇ ${esc(x.label)}</a>`;
      }).join('') : '<p class="hint">（暂无关键产物）</p>';
      const restHtml = rest.length ? `<details class="rpt-sec" style="margin:8px 0 0"><summary>全部文件（${rest.length}）</summary><div class="rpt-sec-body" style="display:flex;flex-wrap:wrap;gap:4px">${rest.map(x => {
        const cls = x.color ? ` dlbtn ${x.color}` : ' dlbtn';
        return `<a class="${cls.trim()}" href="/tool_runs/${encodeURIComponent(run.name)}/${esc(x.f.path)}" title="${esc(x.f.path)}【${esc(x.f.size)}】">⬇ ${esc(x.label)}</a>`;
      }).join('')}</div></details>` : '';
      const fileHtml = open ? `<div class="dlbtn-row">${keyHtml}</div>${restHtml}` : '';
      const rptBtn = (run.name.startsWith('identify_') || run.name.startsWith('contigs_'))
        ? `<button class="btn small" onclick="rhOpenReport('${esc(run.name)}')" title="查看该运行的分类报告（桑基图/分类表/明细）">📊 报告</button> `
        : '';
      return `<div style="border-bottom:1px dashed var(--line-100);padding:5px 0">
        <div style="display:flex;gap:6px;align-items:center;flex-wrap:wrap">
          <b class="mono" style="font-size:12px;cursor:pointer" onclick="rhFiles('${esc(run.name)}')" title="点击展开/收起文件列表">${esc(run.name)}</b>
          <span class="hint">${(run.files || []).length} 文件</span>
          ${chainBtns}
          ${rptBtn}
          <button class="btn small" onclick="rhFiles('${esc(run.name)}')">📁 文件</button>
          <button class="btn small" data-dir="${esc(run.out_dir || '')}" onclick="rhOpenDir(this.dataset.dir)" title="在资源管理器中打开该运行目录（窗口置顶）">📂 打开</button>
          <button class="btn small danger" onclick="rhDel(this, '${esc(run.name)}')" title="删除该运行目录">🗑</button>
        </div>
        ${run.out_dir ? `<div class="hint mono" style="font-size:11px;margin:2px 0 0;word-break:break-all" title="结果输出位置">📍 ${esc(run.out_dir)}</div>` : ''}
        ${open ? (run.name.startsWith('contigs_') ? `<div id="rhcontig-${esc(run.name)}" style="margin-top:8px"><p class="hint">加载病毒序列分类明细…</p></div>` : '') + `<div style="margin-top:4px" class="dlbtn-row">${fileHtml}</div>` : ''}
      </div>`;
    }).join('');
  }
}

async function rhOpenDir(dir) {
  if (!dir) return;
  try {
    const r = await fetch('/api/open_dir', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ path: dir }) });
    if (!r.ok) {
      const d = await r.json().catch(() => ({}));
      toast(t('c.openFail', '打开失败'), d.error || (r.status + ''), {kind: 'failed', ttl: 12000});
    }
  } catch (e) {
    toast(t('c.openFail', '打开失败'), String(e), {kind: 'failed', ttl: 12000});
  }
}

function rhToggle(boxId) {
  const body = $(boxId).querySelector('.rh-body');
  body.style.display = body.style.display === 'none' ? '' : 'none';
}
function rhFiles(name) {
  if (RH_OPEN.has(name)) RH_OPEN.delete(name); else RH_OPEN.add(name);
  rhRefresh();
  if (RH_OPEN.has(name) && name.startsWith('contigs_')) rhLoadContigDetail(name);
}

/* 历史运行点开 → 展示该 run 的病毒序列分类（contig 明细）关键信息。
   读取 virus_classification.tsv，展示每条病毒 contig 的分类/长度/得分。 */
async function rhLoadContigDetail(run) {
  const box = document.getElementById('rhcontig-' + run);
  if (!box) {   // rhRefresh 异步渲染可能未完成，稍后重试
    setTimeout(() => rhLoadContigDetail(run), 300);
    return;
  }
  try {
    const r = await fetch(`/tool_runs/${encodeURIComponent(run)}/virus_classification.tsv`);
    if (!r.ok) { box.innerHTML = '<p class="hint">（无病毒分类明细）</p>'; return; }
    const text = await r.text();
    const lines = text.split(/\r?\n/).filter(Boolean);
    if (lines.length < 2) { box.innerHTML = '<p class="hint">（无病毒分类明细）</p>'; return; }
    const head = lines[0].split('\t');
    const idx = k => head.indexOf(k);
    const ic = idx('contig'), it = idx('taxon'), il = idx('length'), isc = idx('score'), ih = idx('host');
    const rows = lines.slice(1).map(l => l.split('\t'));
    const show = rows.slice(0, 12);
    const n = rows.length;
    const trs = show.map(row => {
      const c = ic >= 0 ? row[ic] : '';
      const t = it >= 0 ? row[it] : '';
      const l = il >= 0 ? row[il] : '';
      const s = isc >= 0 ? row[isc] : '';
      const h = ih >= 0 ? row[ih] : '';
      return `<tr><td class="mono" style="font-size:11px">${esc(c)}</td><td>${esc(t)}</td><td class="num">${esc(l)}</td><td class="num">${esc(s)}</td><td>${esc(h)}</td></tr>`;
    }).join('');
    box.innerHTML = `<div style="font-size:12.5px;font-weight:700;margin-bottom:6px;color:#1a5276">病毒序列分类（contig 明细 · ${n} 条）</div>
      <div style="max-height:260px;overflow:auto;border:1px solid #eee;border-radius:6px">
      <table class="tb"><thead><tr><th>Contig</th><th>分类</th><th class="num">长度</th><th class="num">得分</th><th>宿主</th></tr></thead>
      <tbody>${trs}</tbody></table></div>
      ${n > show.length ? `<p class="hint" style="margin-top:4px">共 ${n} 条，仅显示前 ${show.length} 条（点“📊 报告”看完整）</p>` : ''}`;
  } catch (e) { box.innerHTML = '<p class="hint">明细加载失败</p>'; }
}
async function rhDel(btn, name) {
  if (!confirm(t('ui.rh.delConfirm').replace('{name}', name))) return;
  try {
    const r = await fetch(`/api/tool/runs/${encodeURIComponent(name)}/delete`,
                          { method: 'POST' });
    if (!r.ok) { alert((await r.json()).error || t('ui.common.delFail')); return; }
    RH_OPEN.delete(name);
    rhRefresh();
  } catch (e) { alert(t('ui.common.connShort') + ': ' + e); }
}
function rhGo(btn, runName, boxId, idx) {
  const def = RH_DEFS.find(d => d.box === boxId);
  const run = (RH_DATA || []).find(x => x.name === runName);
  if (!def || !run) return;
  const chain = (def.chains(run.files || [], run) || [])[idx];
  if (!chain) return;
  // 暂存全部回填值 → 跳目标页 → rhApply 落位（同页目标立即填充）
  let fills = Object.assign({}, chain.fills);
  if (chain.page) {
    const samePage = !chain.page.startsWith('/') ||
      location.pathname === chain.page.split('#')[0];
    const missing = Object.keys(fills).filter(k => !$(k));
    if (samePage && !missing.length) {
      rhApplyFills(fills);
      toast('已回填下游输入', Object.keys(fills).join(', '), { ttl: 3000 });
      if (chain.page.includes('#')) location.hash = chain.page.split('#')[1];
      return;
    }
    try { sessionStorage.setItem('vp_runfill', JSON.stringify(fills)); } catch (e) {}
    location.href = chain.page;
  } else {
    rhApplyFills(fills);
  }
}
function rhApplyFills(fills) {
  for (const [id, val] of Object.entries(fills)) {
    const el = $(id);
    if (!el) continue;
    if (el.type === 'checkbox') el.checked = !!val;
    else if (el.tagName === 'SELECT') {
      if ([...el.options].some(o => o.value === val)) el.value = val;
    } else el.value = val;
  }
}
function rhApplyStashed() {
  let fills = null;
  try { fills = JSON.parse(sessionStorage.getItem('vp_runfill') || 'null'); } catch (e) {}
  if (!fills) return;
  const applied = [];
  for (const id of Object.keys(fills)) {
    if ($(id)) { rhApplyFills({ [id]: fills[id] }); applied.push(id); delete fills[id]; }
  }
  if (applied.length) {
    if (Object.keys(fills).length)
      try { sessionStorage.setItem('vp_runfill', JSON.stringify(fills)); } catch (e) {}
    else sessionStorage.removeItem('vp_runfill');
    toast('已回填历史运行的衔接输入', applied.join(', '), { ttl: 4000 });
  }
}
window.__rhOnTasks = tasks => {
  // 任何任务刚结束（running → 终态）时刷新历史区
  const done = new Set(tasks.filter(t => ['done', 'failed', 'cancelled'].includes(t.status))
                            .map(t => t.id));
  const newly = [...done].filter(id => !RH_SEEN.has(id));
  if (RH_SEEN.size && newly.length) rhRefresh();
  RH_SEEN.clear(); tasks.forEach(t => RH_SEEN.add(t.id));
};
const RH_SEEN = new Set();
function rhInitContainers() {
  for (const def of RH_DEFS) {
    const box = $(def.box);
    if (!box || box.dataset.init) continue;
    box.dataset.init = '1';
    box.innerHTML = `<div class="rh-head" style="display:flex;gap:8px;align-items:center;margin-top:8px">
        <button class="btn small" onclick="rhToggle('${def.box}')">🕘 历史运行（<span class="rh-count">…</span>）</button>
        <span class="hint">重启不丢 · 可复用为下游输入 · 可删除（默认折叠）</span>
      </div>
      <div class="rh-body" style="display:none;margin-top:4px"></div>`;
  }
}
document.addEventListener('DOMContentLoaded', () => {
  rhInitContainers();
  rhRefresh();
  setTimeout(rhApplyStashed, 350);   // 等各页自身 init（下拉等）先就绪
});

