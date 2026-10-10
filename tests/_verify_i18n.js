/* i18n 缺失键实证：加载真实 i18n.js，用最小 DOM 桩跑 applyI18n，
 * 断言「字典缺失的键会把元素内容写成 key 字面量」；
 * 并扫源码文本查**重复键定义**（JS 对象字面量后者静默覆盖前者，解析后看不见）。
 *
 * 用法: node tests/_verify_i18n.js [--verbose]
 * 退出码 0 = 无缺失键 / 无重复键；1 = 存在（打印清单）。
 *
 * 设计（2026-10-01 重写）
 * ----------------------
 * 旧版把待检键**硬编码**在脚本里（tk.* / vx.ch* / sc.*）。2026-09-15 体检
 * 报告里记的 vx.chYearT 等 5 个缺陷后来被修掉了（tools.html 现在一个
 * data-i18n="vx.*" 锚点都没有），但 fixture 没跟着改 —— 于是这个脚本
 * 变成**永久红灯**，跑一次红一次，反而训练人忽略它。
 *
 * 现在改为**从真实模板自动提取**锚点：扫 webapp/templates/*.html 里所有
 * data-i18n / data-i18n-ph / data-i18n-title 的值，逐个喂给 applyI18n。
 * 模板改了什么，这里就查什么，不会再过期。
 */
'use strict';
const fs = require('fs');
const path = require('path');

const ROOT = path.dirname(__dirname);
const TPL_DIR = path.join(ROOT, 'webapp', 'templates');
const SRC = fs.readFileSync(path.join(ROOT, 'webapp', 'static', 'i18n.js'), 'utf8');
const VERBOSE = process.argv.includes('--verbose');

/* —— 从模板提取锚点 ——
 * 三类属性对应 applyI18n 的三段：
 *   data-i18n        -> el.innerHTML  = t(key)
 *   data-i18n-ph     -> el.placeholder = t(key)
 *   data-i18n-title  -> el.title      = t(key)
 * 含 Jinja 语法（{{ }} / {% %} / ${}）的键是动态拼的，静态查不了，跳过并计数。
 */
const ATTRS = [
  ['data-i18n', 'i18n', 'innerHTML'],
  ['data-i18n-ph', 'i18nPh', 'placeholder'],
  ['data-i18n-title', 'i18nTitle', 'title'],
];

function extractAnchors() {
  const found = new Map();          // key -> {attr, prop}
  const skipped = [];
  for (const fn of fs.readdirSync(TPL_DIR)) {
    if (!fn.endsWith('.html')) continue;
    const html = fs.readFileSync(path.join(TPL_DIR, fn), 'utf8');
    for (const [attr, prop, target] of ATTRS) {
      const re = new RegExp(attr + '="([^"]*)"', 'g');
      let m;
      while ((m = re.exec(html)) !== null) {
        const key = m[1];
        if (!key) continue;
        if (/[{}$%]/.test(key)) { skipped.push(`${fn}: ${key}`); continue; }
        if (!found.has(key)) found.set(key, { prop, target, file: fn });
      }
    }
  }
  return { anchors: found, skipped };
}

const { anchors, skipped } = extractAnchors();
if (!anchors.size) {
  console.error('✘ 没从模板里提取到任何 data-i18n 锚点 —— 路径或属性名可能变了');
  process.exit(1);
}

/* —— 最小 DOM 桩：只实现 applyI18n 用到的 querySelectorAll + dataset/属性 —— */
function mkEl(prop, key) {
  const el = { dataset: {}, innerHTML: '', placeholder: '', title: '' };
  // dataset 键名转换：data-i18n-ph -> dataset.i18nPh
  el.dataset[prop] = key;
  return el;
}

const els = [];
for (const [key, meta] of anchors) els.push(mkEl(meta.prop, key));

const byProp = { innerHTML: 'i18n', placeholder: 'i18nPh', title: 'i18nTitle' };
const documentStub = {
  querySelector: () => null,
  querySelectorAll: (sel) => {
    // 按选择器把对应属性的元素挑出来 —— 一个元素只有一个锚点，故一一对应
    const prop = sel === '[data-i18n]' ? 'innerHTML'
      : sel === '[data-i18n-ph]' ? 'placeholder'
        : sel === '[data-i18n-title]' ? 'title' : null;
    if (!prop) return [];
    return els.filter(e => e.dataset[byProp[prop]] !== undefined
      && anchors.get(e.dataset[byProp[prop]]).target === prop);
  },
  addEventListener: () => {},          // i18n.js 末尾挂了 DOMContentLoaded 钩子
};
const localStorageStub = { getItem: () => null, setItem: () => {} };

const factory = new Function('document', 'localStorage', 'window',
  SRC + '\nreturn { I18N_DICT, t, applyI18n, setVpLang: (l) => { VP_LANG = l; } };');
const m = factory(documentStub, localStorageStub, globalThis);

const langs = ['zh', 'en'];
const missingByLang = {};
let bad = 0;
for (const lang of langs) {
  // 预置哨兵值：如果 applyI18n 没写，值会保持哨兵；若写了缺失键，值会变成 key 本身。
  // ⚠ 判据必须是「val === key」，不能用「val 没变」——
  //   t(key, null) 在缺键时返回的是 **key 字面量**（见 i18n.js:4030-4036
  //   的 `return fallback != null ? fallback : key`），所以 applyI18n 的
  //   `if (v != null)` 恒真、一定会赋值。踩过：用「没被赋值」当判据 ⇒ 恒过。
  const SENTINEL = '\u0000ORIGINAL\u0000';
  els.forEach(e => {
    e.innerHTML = SENTINEL;
    e.placeholder = SENTINEL;
    e.title = SENTINEL;
  });
  m.setVpLang(lang);
  m.applyI18n(documentStub);

  const broken = [];
  for (const e of els) {
    const key = e.dataset.i18n || e.dataset.i18nPh || e.dataset.i18nTitle;
    const meta = anchors.get(key);
    const val = e[meta.target];
    if (val === key) broken.push({ key, ...meta });
  }
  missingByLang[lang] = broken;
  bad += broken.length;
  console.log(`[${lang}] 模板锚点 ${els.length} 个，字典缺失 ${broken.length} 个`);
  if (VERBOSE) {
    broken.forEach(b => console.log(`   ✘ ${b.key}  (${b.file} · ${b.prop})`));
  }
}

// 中英必须成对：只缺一种语言同样是坏文案
const zhSet = new Set(missingByLang.zh.map(b => b.key));
const enSet = new Set(missingByLang.en.map(b => b.key));
const asymmetric = [...new Set([...zhSet, ...enSet])]
  .filter(k => zhSet.has(k) !== enSet.has(k));
if (asymmetric.length) {
  console.log(`\n⚠ 中英不成对（只有一种语言缺）：${asymmetric.length} 个`);
  asymmetric.forEach(k => {
    console.log(`   ${zhSet.has(k) ? 'zh 缺' : 'en 缺'} ${k}`);
  });
}

if (skipped.length) {
  console.log(`\n（跳过含模板语法的动态键 ${skipped.length} 个，静态查不了）`);
  if (VERBOSE) skipped.forEach(s => console.log(`   ~ ${s}`));
}

/* —— 重复键检测（2026-10-08 新增）——
 * 为什么不能只看 I18N_DICT：JS 对象字面量**允许重复键**，后者静默覆盖前者，
 * 解析出来的对象里根本看不到重复 ⇒ 必须扫**源码文本**。
 * 危害：同一语义两份实现，改前一处**静默无效** —— 本仓最稳定的缺陷模式。
 * 口径陷阱：i18n.js **一行可放多个键**（'nav.home': …, 'dl.title': …），
 * 必须**全行匹配**；只按行首匹配会漏（实测 zh 少算 244 条：2152 vs 2396）。
 */
function findDuplicateKeys(src) {
  const lines = src.split(/\r?\n/);
  const zi = lines.findIndex((l) => l.trim() === 'zh: {');
  const ei = lines.findIndex((l) => l.trim() === 'en: {');
  const xi = lines.findIndex((l, i) => i > ei && l.trim() === '};');
  if (zi < 0 || ei < 0 || xi < 0) {
    // 范围自证：定位不到分段就显式失败，别静默"通过"
    console.error('✘ 无法定位 I18N_DICT 的 zh: / en: 段 —— i18n.js 结构变了，'
      + '重复键检查已失效（请同步更新 _verify_i18n.js）');
    process.exit(1);
  }
  const KEY = /(['"])([A-Za-z][A-Za-z0-9_.\-]*)\1\s*:/g;
  const dups = [];
  const scan = (a, b, lang) => {          // 区间 [a, b)，与上方 eval 切片一致
    const seen = new Map();               // key -> 首次出现的行下标
    for (let i = a; i < b; i++) {
      KEY.lastIndex = 0;
      let m;
      while ((m = KEY.exec(lines[i])) !== null) {
        const k = m[2];
        if (seen.has(k)) {
          dups.push({ lang, key: k, first: seen.get(k) + 1, again: i + 1 });
        } else {
          seen.set(k, i);
        }
      }
    }
  };
  scan(zi + 1, ei - 1, 'zh');             // 跳过 'zh: {' 与其收尾的 '},'
  scan(ei + 1, xi - 1, 'en');
  return dups;
}

const dups = findDuplicateKeys(SRC);
if (dups.length) {
  console.log(`\n⚠ 重复键定义（同一语义多份实现）：${dups.length} 处`);
  dups.forEach((d) => console.log(
    `   ✘ [${d.lang}] ${d.key}  L${d.first} 与 L${d.again} 重复`
    + '（后者覆盖前者，前一处是死代码）'));
} else {
  console.log('\n重复键定义 0 处 ✔');
}

const total = bad + asymmetric.length + dups.length;
if (total) {
  console.log(`\n结论：存在 ${bad} 处「显示为原始 key」+ ${asymmetric.length} 处中英不成对`
    + ` + ${dups.length} 处重复键`);
  if (!VERBOSE) console.log('（加 --verbose 看完整清单）');
} else {
  console.log(`\n结论：模板锚点 ${els.length} 个全部有定义，中英成对，无重复键`);
}
process.exit(total ? 1 : 0);
