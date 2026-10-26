// ---------------- 磁盘水位与存储面板 ----------------
// 2026-09-13 自 app.js 拆出。普通 script 共享全局作用域：本文件在
// app.js 之后加载，可直接使用 $/esc/t/toast/jfetch 等核心工具；页面
// 内联脚本与其它 app-*.js 里的同名调用按全局解析（加载顺序见模板）。
// 注意：本块内引用的其它 app-*.js 函数均为调用期解析（事件/回调触发），
// 不存在加载期交叉依赖。
// ---------------- 磁盘水位与存储面板 ----------------
// 磁盘剩余每次实时取（毫秒级）；目录占用要遍历 6 万+ 文件，后端异步统计，
// 前端轮询等待结果，避免首页被 47GB 的扫描拖住。
async function checkStorageBanner() {
  try {
    const r = await fetch('/api/storage');
    if (!r.ok) return;
    const d = await r.json();
    const low = (d.disks || []).filter(x => x.warn);
    let b = $('storage-banner');
    if (!low.length) { if (b) b.remove(); return; }
    // 本标签页会话内点过「知道了」就不再打扰；关掉浏览器/开新会话会再提醒
    if (sessionStorage.getItem('stBannerDismissed') === '1') {
      if (b) b.remove();
      return;
    }
    if (!b) {
      b = document.createElement('div');
      b.id = 'storage-banner';
      b.style.cssText = 'position:sticky;top:0;z-index:98;background:#f9a825;'
        + 'color:#3e2723;padding:10px 16px;font-size:14px;line-height:1.6;';
      const cb = $('conn-banner');
      if (cb && cb.parentNode) cb.parentNode.insertBefore(b, cb.nextSibling);
      else document.body.prepend(b);
    }
    b.innerHTML = '<b>⚠ ' + t('st.lowTitle', '磁盘空间不足') + '</b><br>'
      + low.map(x => t('ui.st.lbl.' + x.label, x.label) + t('st.disk', '盘')
        + t('ui.st.parenL') + esc(x.path) + t('ui.st.parenR')
        + t('st.onlyLeft', ' 仅剩 ') + esc(x.free_h) + t('ui.st.comma')
        + t('st.suggest', '建议保留 ') + d.threshold_gb + 'GB '
        + t('st.above', '以上')).join('<br>')
      + '<br>' + t('st.lowHint',
          '可到「设置 → 存储与磁盘」查看占用明细；或把数据库目录改到空间更大的盘。')
      + '<div style="margin-top:6px">'
      + '<button class="btn small" onclick="dismissStorageBanner()">'
      + t('st.gotIt', '知道了，本次不再提示') + '</button>'
      + '<button class="btn small" style="margin-left:6px" '
      + 'onclick="sessionStorage.removeItem(\'stBannerDismissed\');'
      + 'location.href=\'/settings\'">'
      + t('st.gotoStorage', '去清理 / 查看占用') + '</button></div>';
  } catch (e) { /* 磁盘检查失败不应影响正常使用 */ }
}

function dismissStorageBanner() {
  try { sessionStorage.setItem('stBannerDismissed', '1'); } catch (e) {}
  const b = $('storage-banner');
  if (b) b.remove();
}

function _storageRow(x, indent) {
  const pad = indent ? 'padding-left:' + (indent * 16) + 'px;' : '';
  const kindTxt = ({ database: t('st.kDb', '数据库'),
                     output: t('st.kOut', '产物'),
                     input: t('st.kIn', '输入'),
                     build: t('st.kBuild', '打包产物') })[x.kind] || x.kind;
  return '<tr>'
    + '<td style="' + pad + '">' + (indent ? '└ ' : '') + esc(x.name) + '</td>'
    + '<td>' + esc(kindTxt) + '</td>'
    + '<td style="text-align:right">' + esc(x.size_h) + '</td>'
    + '<td style="text-align:right">' + (x.files || 0).toLocaleString() + '</td>'
    + '</tr>';
}

async function loadStoragePanel(times) {
  const box = $('storagePanel');
  if (!box) return;
  let d;
  try {
    const r = await fetch('/api/storage?detail=1');
    if (!r.ok) {
      box.innerHTML = '<p class="hint">' + t('st.loadFail', '加载失败') + '</p>';
      return;
    }
    d = await r.json();
  } catch (e) {
    box.innerHTML = '<p class="hint">' + t('st.loadFail', '加载失败') + '</p>';
    return;
  }

  const disks = (d.disks || []).map(x => {
    const pct = x.used_pct || 0;
    const col = x.warn ? '#c62828' : (pct > 85 ? '#f9a825' : '#2e7d32');
    return '<div style="margin-bottom:8px">'
      + '<div style="display:flex;justify-content:space-between;font-size:13px">'
      + '<span>' + esc(t('ui.st.lbl.' + x.label, x.label)) + t('st.disk', '盘') + ' · '
      + '<span style="color:#666">' + esc(x.path) + '</span></span>'
      + '<span>' + esc(x.free_h) + ' ' + t('st.free', '可用') + ' / '
      + esc(x.total_h) + '</span></div>'
      + '<div style="height:8px;background:#eee;border-radius:4px;overflow:hidden">'
      + '<div style="height:100%;width:' + pct + '%;background:' + col + '"></div>'
      + '</div></div>';
  }).join('');

  const items = d.scan.items || [];
  let rows = '';
  for (const it of items) {
    rows += _storageRow(it, 0);
    for (const c of (it.children || [])) rows += _storageRow(c, 1);
  }
  const table = items.length
    ? '<table class="tbl" style="width:100%;font-size:13px;margin-top:6px">'
      + '<thead><tr><th>' + t('st.colName', '目录') + '</th><th>'
      + t('st.colKind', '类别') + '</th><th style="text-align:right">'
      + t('st.colSize', '占用') + '</th><th style="text-align:right">'
      + t('st.colFiles', '文件数') + '</th></tr></thead>'
      + '<tbody>' + rows + '</tbody></table>'
    : '';

  const status = d.scan.ready
    ? t('st.at', '统计于 ') + esc(d.scan.scanned_at || '')
    : (d.scan.scanning
        ? t('st.scanning', '正在统计（数据库较大，约需几十秒）…')
        : t('st.waiting', '等待统计…'));

  box.innerHTML = disks
    + '<div style="display:flex;align-items:center;gap:8px;margin:10px 0 4px">'
    + '<b style="font-size:13px">' + t('st.breakdown', '目录占用') + '</b>'
    + '<span class="hint">' + status + '</span>'
    + '<button class="btn small" onclick="rescanStorage()" style="margin-left:auto">'
    + t('st.rescan', '重新统计') + '</button></div>'
    + (d.scan.error
        ? '<p class="hint" style="color:#c62828">' + esc(d.scan.error) + '</p>' : '')
    + table
    + '<p class="hint" style="margin-top:8px">' + t('st.tip',
        '提示：dist（打包产物）、tool_runs（工具运行）、logs（日志）可安全清理；'
        + 'databases 下的库请不要手删——重建宿主库代价很高。'
        + '空间紧张时建议把「数据库目录」改到更大的盘。') + '</p>';

  if (!d.scan.ready && times > 0) {
    setTimeout(function () { loadStoragePanel(times - 1); }, 3000);
  }
}

async function rescanStorage() {
  try { await fetch('/api/storage/rescan', { method: 'POST' }); } catch (e) {}
  loadStoragePanel(40);
}

checkStorageBanner();
loadStoragePanel(0);

