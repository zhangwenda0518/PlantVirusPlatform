/* 主导航交互（2026-10-21 布局重构：顶部横条 → 左侧固定导航栏）
 *
 * 职责只有三件：
 *   1. 分组展开/收起（点击 .sn-caret 或组链接右侧的箭头区）；
 *   2. 窄屏抽屉开合（☰ / 遮罩 / Esc / 点链接后自动关）；
 *   3. 顶栏标题填当前分组名（侧栏收起时主区仍知道自己在哪）。
 *
 * 不负责：链接跳转、当前页高亮（服务端渲染时已定）、上下文胶囊（app.js）。
 * 与 app.js / i18n.js 的契约：本脚本先于它们执行，只做 DOM 绑定，
 * 不依赖 $ / t / esc 等全局函数；所有查询都是原生 API。
 */
(function () {
  'use strict';

  var STORE_KEY = 'vp_nav_open';       // 记住展开的分组（跨页面保持）

  function byId(id) { return document.getElementById(id); }

  /* ---------------- 1. 分组展开 ---------------- */

  function loadOpen() {
    try {
      var raw = localStorage.getItem(STORE_KEY);
      var arr = raw ? JSON.parse(raw) : null;
      return Array.isArray(arr) ? arr : [];
    } catch (e) { return []; }
  }

  function saveOpen() {
    try {
      var ids = [];
      document.querySelectorAll('.sn-group.open').forEach(function (g) {
        ids.push(g.dataset.group);
      });
      localStorage.setItem(STORE_KEY, JSON.stringify(ids));
    } catch (e) { /* 隐私模式下 localStorage 不可写，忽略 */ }
  }

  function setOpen(group, open) {
    group.classList.toggle('open', !!open);
    var link = group.querySelector('a.sn-link');
    if (link) link.setAttribute('aria-expanded', open ? 'true' : 'false');
  }

  function initGroups() {
    var groups = document.querySelectorAll('.sn-group');
    if (!groups.length) return;

    // 恢复上次展开状态：服务端已按当前组打上 open 的保持不动，
    // 其余按 localStorage 补开（同一页面里可能有多个组被用户手动展开过）
    var remembered = loadOpen();
    groups.forEach(function (g) {
      if (g.classList.contains('open')) return;         // 当前组由服务端标记
      if (remembered.indexOf(g.dataset.group) >= 0) setOpen(g, true);
    });

    groups.forEach(function (g) {
      // 整行标题即开关：点一次展开、再点收起（热区 = 整行，不再依赖小箭头）。
      // href 指向组概览页，这里拦下来只做开合；想进组概览点展开后的
      // 「🏠 组概览」那一项（保持"点标题=开合"的单一语义，不猜用户意图）。
      var link = g.querySelector('a.sn-link');
      if (link) {
        link.addEventListener('click', function (e) {
          e.preventDefault();
          setOpen(g, !g.classList.contains('open'));
          saveOpen();
        });
      }
      // 箭头只是视觉指示（span），点击落到整行标题上即可；这里再兜一层，
      // 万一日后样式把箭头做成独立可点元素也不会漏。
      var caret = g.querySelector('.sn-caret');
      if (caret) {
        caret.addEventListener('click', function (e) {
          e.preventDefault();
          e.stopPropagation();
          setOpen(g, !g.classList.contains('open'));
          saveOpen();
        });
      }
    });
  }

  /* ---------------- 2. 窄屏抽屉 ---------------- */

  function initDrawer() {
    var nav = byId('sideNav');
    var toggle = byId('snToggle');
    var scrim = byId('snScrim');
    if (!nav || !toggle) return;

    function setOpen(open) {
      nav.classList.toggle('open', open);
      toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
      if (scrim) scrim.hidden = !open;
      document.body.classList.toggle('sn-drawer-open', open);
    }

    toggle.addEventListener('click', function () {
      setOpen(!nav.classList.contains('open'));
    });
    var closeBtn = byId('snClose');
    if (closeBtn) closeBtn.addEventListener('click', function () { setOpen(false); });
    if (scrim) scrim.addEventListener('click', function () { setOpen(false); });

    document.addEventListener('keydown', function (e) {
      if (e.key === 'Escape' && nav.classList.contains('open')) setOpen(false);
    });

    // 抽屉里点了链接就关（窄屏下侧栏盖住了主区，不关会挡住目标页）
    nav.addEventListener('click', function (e) {
      var a = e.target.closest && e.target.closest('a');
      if (a && window.matchMedia('(max-width: 1024px)').matches) setOpen(false);
    });

    // 拉宽到桌面断点时清掉抽屉态，否则 body 上残留的滚动锁定会生效
    var mq = window.matchMedia('(max-width: 1024px)');
    var onChange = function (ev) { if (!ev.matches) setOpen(false); };
    if (mq.addEventListener) mq.addEventListener('change', onChange);
    else if (mq.addListener) mq.addListener(onChange);
  }

  /* ---------------- 3. 顶栏标题 ---------------- */

  function initTitle() {
    var el = byId('tbTitle');
    if (!el) return;
    // 用服务端标记的 data-current（当前组）而不是「第一个 .open」：
    // .open 还包含用户手动展开过、靠 localStorage 恢复的其它组，
    // 拿第一个会让标题显示成别的组（实测 /results 显示「样本处理」）。
    var cur = document.querySelector('.sn-group[data-current] > a.sn-link');
    var text = cur ? (cur.textContent || '').trim() : '';
    if (!text) {
      // 非分组页（总览 / 任务中心 / 结果中心 / 设置 / 手册）：
      // 取侧栏底部区或首页链接里 active 的那条
      var act = document.querySelector('.sn-foot a.active')
             || document.querySelector('.sn-links > a.active');
      text = act ? (act.textContent || '').trim() : '';
    }
    el.textContent = text;
  }

  function boot() {
    initGroups();
    initDrawer();
    initTitle();
    // i18n 落地后再算一次标题（见 i18n.js 的 vp:i18n 广播）
    document.addEventListener('vp:i18n', initTitle);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else {
    boot();
  }
})();
