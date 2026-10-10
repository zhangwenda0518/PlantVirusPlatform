// ---------------- 文件浏览（记忆上次目录） + 路径输入历史 ----------------
// 2026-09-13 自 app.js 拆出。普通 script 共享全局作用域：本文件在
// app.js 之后加载，可直接使用 $/esc/t/toast/jfetch 等核心工具；页面
// 内联脚本与其它 app-*.js 里的同名调用按全局解析（加载顺序见模板）。
// 注意：本块内引用的其它 app-*.js 函数均为调用期解析（事件/回调触发），
// 不存在加载期交叉依赖。
// ---------------- 文件浏览（记忆上次目录） ----------------
let browseTarget = null;
let browseMode = 'file';   // 'file' 选数据文件 | 'dir' 选目录
let browseCwd = '.';
// 目录选择回调（无输入框可回填的场景：下拉「✍ 自选目录」等）。
// pickDir() 确认时若非空则把所选目录交给回调；browse()/browseDir() 会清掉它，
// 避免上一次未消费的回调串场。
let browseDirCb = null;

function browse(inputId) {
  browseMode = 'file';
  browseTarget = inputId;
  browseDirCb = null;
  $('dlgMask').style.display = 'flex';
  const dp = $('dlgPath');
  if (dp) wirePathInput(dp);              // 路径栏：粘贴清洗 + 跳转历史
  const ttl = $('dlgTitle');
  if (ttl) ttl.textContent = t('c.chooseFile');
  const last = localStorage.getItem('vp_browse_cwd');
  loadBrowse(last || '.');
}

function browseDir(inputId) {
  browseMode = 'dir';
  browseTarget = inputId;
  browseDirCb = null;
  $('dlgMask').style.display = 'flex';
  const dp = $('dlgPath');
  if (dp) wirePathInput(dp);              // 路径栏：粘贴清洗 + 跳转历史
  const ttl = $('dlgTitle');
  if (ttl) ttl.textContent = t('c.chooseDir', '选择目录');
  const last = localStorage.getItem('vp_browse_cwd');
  loadBrowse(last || '.');
}

/* 打开「选目录」对话框，确认后把所选目录交给 cb（不回填任何输入框）。
   供下拉框「自选目录…」项使用：选中后把目录作为动态选项插回下拉。 */
function pickDirWith(cb) {
  if (typeof cb !== 'function') return;
  browseMode = 'dir';
  browseTarget = null;
  browseDirCb = cb;
  $('dlgMask').style.display = 'flex';
  const dp = $('dlgPath');
  if (dp) wirePathInput(dp);
  const ttl = $('dlgTitle');
  if (ttl) ttl.textContent = t('c.chooseDir', '选择目录');
  const last = localStorage.getItem('vp_browse_cwd');
  loadBrowse(last || '.');
}

function pickDir() {
  if (browseCwd === '此电脑') { alert(t('ui.br.needDir')); return; }
  if (browseDirCb) {
    const cb = browseDirCb;
    browseDirCb = null;
    closeDlg();
    cb(browseCwd);
    return;
  }
  if (browseTarget && $(browseTarget)) {
    $(browseTarget).value = browseCwd;
    pushPathHist(browseTarget, browseCwd);
  }
  closeDlg();
}

function _joinBrowse(cwd, name) {
  /* 拼接子路径：
     - '此电脑' 层级的子项就是盘符绝对路径（如 D:\），直接返回；
     - 盘根目录（D:\ 带尾反斜杠）直接拼名字；
     - 其余用 / 连接（服务端 normpath 兼容正斜杠）。 */
  if (cwd === '此电脑') return name;
  if (/[\\/]$/.test(cwd)) return cwd + name;
  return cwd + '/' + name;
}

async function loadBrowse(path, _isRetry) {
  try {
    // 选目录模式带 all=1：服务端不限文件扩展名，方便查看目录里有什么
    const r = await fetch('/api/browse?path=' + encodeURIComponent(path) +
                          (browseMode === 'dir' ? '&all=1' : ''));
    if (!r.ok) {
      // 记忆的上次目录可能已被删除/改名（U 盘拔出、运行目录清理等）：
      // 退回「此电脑」重试一次，只有根视图也打不开才报错。
      if (!_isRetry && path && path !== '.') {
        try { localStorage.removeItem('vp_browse_cwd'); } catch (e) {}
        loadBrowse('.', true);
        return;
      }
      alert(t('ui.br.openDirFail'));
      return;
    }
    const d = await r.json();
    browseCwd = d.cwd;
    try { localStorage.setItem('vp_browse_cwd', d.cwd === '此电脑' ? '' : d.cwd); } catch (e) {}
    $('dlgPath').value = d.cwd;
    const pick = $('dlgPickDir');
    if (pick) pick.style.display = browseMode === 'dir' ? '' : 'none';
    const items = [];
    if (d.parent !== '') {
      items.push(`<div class="fitem dir" onclick="loadBrowse(${JSON.stringify(d.parent).replace(/"/g, '&quot;')})">⬆ ..</div>`);
    } else if (d.cwd !== '此电脑') {
      items.push(`<div class="fitem dir" onclick="loadBrowse('此电脑')">🖥 此电脑</div>`);
    }
    for (const dir of d.dirs) {
      const target = _joinBrowse(d.cwd, dir);
      items.push(`<div class="fitem dir" onclick="loadBrowse(${JSON.stringify(target).replace(/"/g, '&quot;')})">📂 ${esc(dir)}${browseMode === 'dir' ? ` <span class="hint">${t('c.dirHint', '双击进入 · 或到上级选此目录')}</span>` : ''}</div>`);
    }
    if (browseMode === 'file') {
      for (const f of d.files) {
        const full = _joinBrowse(d.cwd, f.name);
        items.push(`<div class="fitem" onclick="pickFile(${JSON.stringify(full).replace(/"/g, '&quot;')})">
          🧬 ${esc(f.name)} <span class="hint">${esc(f.size)}</span></div>`);
      }
    } else if (d.files.length) {
      // 目录模式：文件仅作内容预览（不可选，选中的是目录本身）
      for (const f of d.files) {
        items.push(`<div class="fitem" style="cursor:default">
          📄 ${esc(f.name)} <span class="hint">${esc(f.size)}</span></div>`);
      }
    }
    if (browseMode === 'dir' && !d.dirs.length && !d.files.length) {
      items.push(`<p class="hint">${t('c.emptyDir', '（空目录）')}</p>`);
    }
    $('dlgList').innerHTML = items.join('') || `<p class="hint">${t('c.noMatchFile')}</p>`;
  } catch (e) {
    if (netFail(e)) setConnBanner(true);
  }
}

function pickFile(fullPath) {
  if (browseTarget && $(browseTarget)) {
    $(browseTarget).value = fullPath;
    pushPathHist(browseTarget, fullPath);
  }
  closeDlg();
}

/* 路径文本清洗：资源管理器「复制文件地址」带首尾引号、浏览器地址栏是
   file:/// 链接、表格单元格粘贴可能带换行——统一裁成可直接使用的路径。
   只在文本“长得像路径”时才接管（引号 / 盘符 / file:// / 双反斜杠 / 根斜杠
   开头），避免干扰普通文本输入。 */
function cleanPathText(s) {
  const raw = String(s == null ? '' : s);
  const x = raw.trim();
  if (!x) return raw;
  if (!/^(["']|file:\/\/|[A-Za-z]:[\\/]|\\\\|\/)/.test(x)) return raw;
  let out = x.replace(/^file:\/\/\/?/i, '')
             .replace(/^["']+/, '').replace(/["']+$/, '').trim();
  out = out.split(/\r?\n/)[0].trim();
  return out || raw;
}

/* 路径栏粘贴/输入 + 回车（或点「跳转」）：
   - 指向文件：文件模式 → 直接选中；目录模式 → 跳到其所在目录；
   - 指向目录：跳进该目录（目录模式再用「选择此目录 ✔」确认）。 */
async function jumpBrowse() {
  const box = $('dlgPath');
  if (!box) return;
  const raw = cleanPathText(box.value);
  if (!raw) return;
  box.value = raw;
  try {
    const r = await fetch('/api/path_info?path=' + encodeURIComponent(raw));
    const d = await r.json();
    if (!d.exists) {
      toast(t('c.pathMissing', '路径不存在'), raw, {kind: 'failed', ttl: 5000});
      return;
    }
    pushPathHist('dlgPath', raw);
    if (!d.is_dir && browseMode === 'file') { pickFile(d.path); return; }
    loadBrowse(d.is_dir ? d.path : (d.parent || d.path));
  } catch (e) {
    if (netFail(e)) setConnBanner(true);
  }
}

// ---------------- 路径输入历史（点输入框 → 弹最近用过的路径 → 点选回填） ----------------
/* 按输入框 id 存 localStorage（vp_path_hist），每 id 最多 10 条、去重、最新在前。
   入史时机：change（手输失焦 / 粘贴清洗后）、pickFile/pickDir、拖拽上传成功、
   对话框路径栏跳转。示例按钮等程序化回填不触发 change，故意不入史。 */
const PATH_HIST_KEY = 'vp_path_hist';
const PATH_HIST_MAX = 10;
let _histBox = null, _histInput = null;

function getPathHist(id) {
  try {
    const all = JSON.parse(localStorage.getItem(PATH_HIST_KEY) || '{}') || {};
    return Array.isArray(all[id]) ? all[id] : [];
  } catch (e) { return []; }
}
function pushPathHist(id, p) {
  const v = String(p || '').trim();
  if (!id || !v || !$(id)) return;          // 占位目标（如测试用 __probe__）不入史
  try {
    const all = JSON.parse(localStorage.getItem(PATH_HIST_KEY) || '{}') || {};
    if ((all[id] || [])[0] === v) return;   // 已是最新，免写
    const arr = (all[id] || []).filter(x => x !== v);
    arr.unshift(v);
    all[id] = arr.slice(0, PATH_HIST_MAX);
    localStorage.setItem(PATH_HIST_KEY, JSON.stringify(all));
  } catch (e) {}
}
function dropPathHist(id, p) {
  try {
    const all = JSON.parse(localStorage.getItem(PATH_HIST_KEY) || '{}') || {};
    all[id] = (all[id] || []).filter(x => x !== p);
    localStorage.setItem(PATH_HIST_KEY, JSON.stringify(all));
  } catch (e) {}
}
function clearPathHist(id) {
  try {
    const all = JSON.parse(localStorage.getItem(PATH_HIST_KEY) || '{}') || {};
    delete all[id];
    localStorage.setItem(PATH_HIST_KEY, JSON.stringify(all));
  } catch (e) {}
}

function closePathHist() {
  if (_histBox) { _histBox.remove(); _histBox = null; _histInput = null; }
}

function showPathHist(input) {
  closePathHist();
  const id = input.id;
  if (!id) return;
  const items = getPathHist(id);
  if (!items.length) return;
  const box = document.createElement('div');
  box.className = 'pathhist';
  box.innerHTML =
    `<div class="ph-head"><span>${t('c.recentPaths', '最近使用')}</span>` +
    `<button type="button" class="ph-clear">${t('c.histClear', '清空')}</button></div>` +
    items.map(p =>
      `<div class="ph-item" data-p="${esc(p)}">` +
      `<span class="ph-path" title="${esc(p)}">${esc(p)}</span>` +
      `<button type="button" class="ph-del" title="${t('c.histDel', '删除这条')}">✕</button></div>`
    ).join('');
  box.addEventListener('pointerdown', e => {
    e.preventDefault();                     // 防止输入框因失焦抢先关闭下拉
    const item = e.target.closest('.ph-item');
    if (!item) return;
    if (e.target.closest('.ph-del')) {
      dropPathHist(id, item.dataset.p);
      showPathHist(input);                  // 删除后原地重画
      return;
    }
    input.value = item.dataset.p;
    input.dispatchEvent(new Event('input', { bubbles: true }));
    input.dispatchEvent(new Event('change', { bubbles: true }));
    closePathHist();
  });
  box.querySelector('.ph-clear').addEventListener('pointerdown', e => {
    e.preventDefault();
    clearPathHist(id);
    closePathHist();
    toast(t('c.histClear', '清空'), t('c.histCleared', '已清空该输入框的路径历史'), { ttl: 2500 });
  });
  const host = input.parentElement;
  if (!host) return;
  host.style.position = 'relative';
  host.appendChild(box);
  _histBox = box; _histInput = input;
}
// 点输入框/下拉以外任意处收起（capture：抢在页面其它 handler 之前）
document.addEventListener('pointerdown', e => {
  if (!_histBox) return;
  if (e.target === _histInput || _histBox.contains(e.target)) return;
  closePathHist();
}, true);

function closeDlg() { $('dlgMask').style.display = 'none'; }

