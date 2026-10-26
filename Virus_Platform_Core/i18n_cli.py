# -*- coding: utf-8 -*-
"""CLI 与运行日志的中英双语层（2026-10-26）。

背景：网页端的中英切换走 `webapp/static/i18n.js`（前端字典），但命令行帮助、
`--cli` 子命令的状态行以及写入 `run/logs/*.log` 的运行日志正文，全部由 Python
侧直接 `print` / `logger.log` 产生，前端字典管不到。本模块在**输出边界**统一
翻译，避免给 900+ 处消息各写两份文案：

* 语言解析优先级：``VP_CLI_LANG`` / ``VP_LANG`` 环境变量 > 命令行 ``--lang``
  > ``platform.json`` 的 ``language``（即「设置 → 语言」，网页与 CLI 共用同一
  开关）> ``zh``。
* :func:`tr` 三段式：**精确表** → **模板正则**（f-string 渲染后的消息用占位符
  还原成正则）→ **短语表**（兜底，处理中英混排与未登记的新消息）。英文模式下
  反复最多 3 轮，让「嵌套消息」（如 ``⚠ 阶段 ① 完成: 文件不存在``）逐层英文化。
* :func:`localize_parser` 就地改写 argparse 的 description/help/epilog，因此
  `add_argument` 处不需要写两份文案。
* 目录（中文→英文）由 ``i18n_cli_catalog.py`` 提供，可由脚本重新生成；缺失时
  本模块退化为「原样输出」，不会影响功能。

调用点：``utils.TaskLogger.log``（所有 logger 的唯一出口）、``main.py`` 的
``_say``、以及 :func:`localize_parser`。
"""
from __future__ import annotations

import os
import re
import sys

__all__ = ['lang', 'set_lang', 'enabled', 'tr', 'tr_argv_hint',
           'localize_parser', 'residual_cjk', 'install_print_hook']

_LANG = None                 # None=未解析；'zh' | 'en'
_EXACT = {}                  # 无占位符的整串映射
_RULES = []                  # [(head, compiled_regex, repl)]
_BY_CHAR = {}                # 首字符 -> [rule, ...]，用于快速筛选
_PHRASES = []                # [(zh, en)]，长串优先
_CJK = re.compile(r'[\u4e00-\u9fff]')
_PH = re.compile(r'\{[^{}]*\}')
_MARKERS = ' \t\r\n·:：,，;；→-—✔✘△⚠⏳✅❌'


# ---------------------------------------------------------------- 语言解析
def _env_lang():
    for k in ('VP_CLI_LANG', 'VP_LANG'):
        v = (os.environ.get(k) or '').strip().lower()
        if v in ('zh', 'en'):
            return v
    return None


def _argv_lang(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    for i, a in enumerate(argv):
        if a == '--lang' and i + 1 < len(argv):
            v = argv[i + 1].strip().lower()
            if v in ('zh', 'en'):
                return v
        if a.startswith('--lang='):
            v = a.split('=', 1)[1].strip().lower()
            if v in ('zh', 'en'):
                return v
    return None


def _cfg_lang():
    try:
        from .config import get_config
        return get_config().lang
    except Exception:
        return None


def _resolve():
    return (_env_lang() or _argv_lang() or _cfg_lang() or 'zh')


def set_lang(value=None):
    """显式指定语言；``None`` 表示按环境变量/配置重新解析。"""
    global _LANG
    if value is None:
        _LANG = None
    else:
        v = str(value).strip().lower()
        _LANG = v if v in ('zh', 'en') else None
    return lang()


def lang():
    global _LANG
    if _LANG is None:
        _LANG = _resolve()
        _ensure_catalog()
    return _LANG


def enabled():
    return lang() == 'en'


# ---------------------------------------------------------------- 目录加载
def _rule_head(tpl):
    """取模板里第一个够长的中文字面块，作为快速命中索引。"""
    for part in _PH.split(tpl):
        p = part.strip(_MARKERS)
        if len(p) >= 2 and _CJK.search(p):
            return p
    m = _CJK.search(tpl)
    return m.group(0) if m else ''


def _index_char(head, tpl=''):
    """索引键：head 里第一个中文字（head 可能以 ASCII 开头，如 'ICTV 库目录'）。"""
    m = _CJK.search(head) or _CJK.search(tpl)
    return m.group(0) if m else ''


_PH_KEY = re.compile(r"'[^']*'|\"[^\"]*\"")   # 占位符内部引号字面量（可被翻译）


def _ph_keys(tpl):
    """占位符的可比骨架序列（引号内字面量归一），用于中英模板对应。"""
    return [_PH_KEY.sub("''", x) for x in _PH.findall(tpl)]


def _compile_rule(tpl, en):
    """把 ``阶段 {k} 完成`` / ``stage {k} done`` 之类还原成正则 + 替换串。

    英文语序可以与中文不同（占位符顺序变化）：捕获组按 **中文模板** 的位置
    命名 g0..gn，替换串再按英文模板里每个占位符对应的中文表达式取用；
    占位符里带引号字面量时（``{'预演' if dry else '完成'}``）比较骨架，
    允许英文侧翻译这些字面量。
    """
    tpl, en = tpl.strip(), en.strip()
    zh_ph, en_ph = _PH.findall(tpl), _PH.findall(en)
    if not zh_ph or sorted(_ph_keys(tpl)) != sorted(_ph_keys(en)):
        return None
    zh_key = _ph_keys(tpl)
    en_key = _ph_keys(en)
    # 中文占位符（按骨架）-> 组名；同一表达式重复出现时按出现次序区分
    zh_map, seen = {}, {}
    for i, ph in enumerate(zh_key):
        k = (ph, seen.get(ph, 0))
        seen[ph] = k[1] + 1
        zh_map[k] = 'g%d' % i
    pieces, names = [], []
    for i, part in enumerate(_PH.split(tpl)):
        pieces.append(re.escape(part))
        if i < len(zh_ph):
            name = 'g%d' % i
            names.append(name)
            pieces.append('(?P<%s>.+?)' % name)
    # 最后一个捕获组放宽到行尾，避免惰性匹配把尾巴丢掉
    tail_lazy = '(?P<%s>.+?)' % names[-1]
    pattern = ''.join(pieces)
    if pattern.endswith(tail_lazy):
        pattern = pattern[:-len(tail_lazy)] + '(?P<%s>.+)' % names[-1]
    repl, seen2 = '', {}
    en_parts = _PH.split(en)
    for i, part in enumerate(en_parts):
        repl += part.replace('\\', '\\\\')
        if i < len(en_ph):
            ph = en_key[i]
            idx = seen2.get(ph, 0)
            seen2[ph] = idx + 1
            repl += '\\g<%s>' % zh_map.get((ph, idx), 'g%d' % i)
    try:
        return re.compile('^' + pattern + '$', re.S), repl
    except re.error:
        return None


def _ensure_catalog():
    global _EXACT, _RULES, _BY_CHAR, _PHRASES
    if _RULES or _EXACT:
        return
    try:
        from . import i18n_cli_catalog as cat
    except Exception:
        return
    _EXACT = dict(getattr(cat, 'EXACT', {}) or {})
    _PHRASES = list(getattr(cat, 'PHRASES', []) or [])
    rules = []
    for tpl, en in (getattr(cat, 'RULES', {}) or {}).items():
        if not en or tpl == en:
            continue
        if _PH.search(tpl):
            got = _compile_rule(tpl, en)
            if got:
                rules.append((_rule_head(tpl), got[0], got[1]))
        else:
            _EXACT.setdefault(tpl, en)
    _RULES = rules
    _BY_CHAR = {}
    for head, rx, repl in rules:
        ch = _index_char(head, rx.pattern)
        if ch:
            _BY_CHAR.setdefault(ch, []).append((head, rx, repl))


# ---------------------------------------------------------------- 翻译主体
def _apply_rules(text):
    if not _CJK.search(text):
        return text
    cands = []
    for ch in set(_CJK.findall(text)):
        cands.extend(_BY_CHAR.get(ch, ()))
    if not cands:
        return text
    out = text
    for head, rx, repl in cands:
        if head not in out:
            continue
        out = rx.sub(repl, out)
    return out


def _apply_phrases(text):
    for zh, en in _PHRASES:
        if zh in text:
            text = text.replace(zh, en)
    return text


# 语言相关的时间/数量单位（2026-10-26）：这些值由 fmt_eta / 各阶段用 f-string
# 直接拼进句子里（"预估 ~1秒"、"保留 12 条"），模板翻译管不到，最后统一收口。
_UNITS = [
    (re.compile(r'(\d+)\s*(?:小时|时)\s*(\d+)\s*分(?:钟)?'), r'\1h \2m'),
    (re.compile(r'(\d+)\s*分(?:钟)?\s*(\d+)\s*秒'), r'\1m \2s'),
    (re.compile(r'(\d+(?:\.\d+)?)\s*秒'), r'\1s'),
    (re.compile(r'(\d+(?:\.\d+)?)\s*毫秒'), r'\1ms'),
    (re.compile(r'(\d+)\s*分钟'), r'\1min'),
    (re.compile(r'(\d+)\s*个?小时'), r'\1h'),
    (re.compile(r'(\d+)\s*天'), r'\1d'),
    (re.compile(r'(\d+)\s*对'), r'\1 pairs'),
    (re.compile(r'(\d+)\s*科'), r'\1 families'),
    (re.compile(r'(\d+)\s*属'), r'\1 genera'),
    (re.compile(r'(\d+)\s*种'), r'\1 species'),
    (re.compile(r'(\d+)\s*行'), r'\1 rows'),
    (re.compile(r'(\d+)\s*列'), r'\1 cols'),
    (re.compile(r'(\d+)\s*片'), r'\1 parts'),
    (re.compile(r'(\d+)\s*次'), r'\1x'),
    (re.compile(r'(\d+)\s*倍'), r'\1x'),
    (re.compile(r'(\d+(?:\.\d+)?)\s*万'), r'\1e4'),
    (re.compile(r'(\d+)\s*条'), r'\1'),
    (re.compile(r'(\d+)\s*张'), r'\1'),
    (re.compile(r'(\d+)\s*个'), r'\1'),
]
_PUNCT_EN = [('，', ', '), ('；', '; '), ('：', ': '), ('（', ' ('), ('）', ')'),
             ('、', ', '), ('。', '. '), ('！', '!'), ('？', '?'), ('～', '~')]
# 英文输出里不该出现全角标点（2026-10-26）：逐词翻译后可能出现
# "identity not recorded（no host_db.json manifest）" 这种残留。
_FULLWIDTH = re.compile(r'[，；：（）【】〔〕、「」。！？～]')


def _units_pass(text):
    """把数字后面残留的中文单位换成英文（只处理「数字+单位」形态）。"""
    out = text
    for rx, repl in _UNITS:
        out = rx.sub(repl, out)
    return out


def _punct_pass(text):
    """英文模式下把残留的全角标点换成半角，并收口标点两侧的空格。

    只收敛紧挨 ASCII 标点的多余空格，表格靠连续空格对齐的地方保持原样。
    """
    out = text
    for zh, en in _PUNCT_EN:
        out = out.replace(zh, en)
    out = re.sub(r'[ ]{2,}(?=[,;:.!?)])', '', out)
    out = re.sub(r'(?<=[,;:.!?])[ ]{2,}', ' ', out)
    out = re.sub(r'[ ]{2,}(?=\()', ' ', out)
    return out


def _tr_one(body):
    hit = _EXACT.get(body)
    if hit is not None:
        return _punct_pass(hit) if _FULLWIDTH.search(hit) else hit
    out = body
    for _ in range(3):
        prev = out
        out = _apply_rules(out)
        if out == prev:
            break
    if _CJK.search(out) and _PHRASES:
        out = _apply_phrases(out)
    out = _units_pass(out)
    out = _punct_pass(out)
    return out


def tr(text):
    """把一条消息翻成当前语言的英文（zh 模式原样返回）。

    保留首尾空白（CLI 靠缩进对齐），逐行处理多行文本；未登记的动态串原样
    返回，绝不吞内容。
    """
    if not isinstance(text, str) or not text:
        return text
    if not enabled():
        return text
    if not _CJK.search(text):
        # 已经是英文，但仍可能夹着拼接留下的全角标点：只做标点收口。
        return _punct_pass(text) if _FULLWIDTH.search(text) else text
    if not _EXACT and not _RULES:
        _ensure_catalog()
    body = text.strip()
    if not body:
        return text
    lead = text[:len(text) - len(text.lstrip())]
    tail = text[len(text.rstrip()):]
    # 先整块匹配（帮助文本/多行 banner 在目录里是整串登记的），
    # 整块没命中的长文本再逐行处理，让每行都能命中行级模板。
    out = _tr_one(body)
    if '\n' in out and _CJK.search(out):
        out = '\n'.join(_tr_one(x) if _CJK.search(x) else x
                        for x in out.split('\n'))
    return lead + out + tail


def residual_cjk(text):
    """测试用：返回消息里残留的中文片段列表。"""
    if not isinstance(text, str):
        return []
    return re.findall(r'[\u4e00-\u9fff][\u4e00-\u9fff\uff00-\uffef]*', text)


def tr_argv_hint(text):
    """argparse 的 usage 提示（英文 argparse 自带英文，这里只处理中文文案）。"""
    return tr(text)


_PRINT_HOOKED = False


def install_print_hook():
    """把 builtins.print 包一层：核心模块里的直接 print 也跟随语言开关。

    CLI 与网页进程启动时各调一次即可。英文模式且含中文才翻译，纯 ASCII 走
    快速路径，开销可以忽略；未登记的动态串原样输出。
    """
    global _PRINT_HOOKED
    import builtins
    if _PRINT_HOOKED or getattr(builtins.print, '_vp_i18n', False):
        _PRINT_HOOKED = True
        return
    _orig = builtins.print

    def _print(*args, **kw):
        if args and enabled():
            args = tuple(tr(a) if isinstance(a, str) else a for a in args)
        _orig(*args, **kw)

    _print._vp_i18n = True          # type: ignore[attr-defined]
    builtins.print = _print
    _PRINT_HOOKED = True


# ---------------------------------------------------------------- argparse
def localize_parser(p):
    """就地英文化 argparse 解析器；zh 模式直接返回。

    只改文案字段（description/help/epilog/_choices_actions），不动 options、
    默认值与 dest，因此行为完全一致。
    """
    if p is None or not enabled():
        return p
    try:
        import argparse
    except Exception:            # pragma: no cover
        return p
    if getattr(p, 'description', None):
        p.description = tr(p.description)
    if getattr(p, 'epilog', None):
        p.epilog = tr(p.epilog)
    for a in getattr(p, '_actions', ()):  # noqa: SLF001 - argparse 无公开 API
        if getattr(a, 'help', None) and a.help is not argparse.SUPPRESS:
            a.help = tr(a.help)
        if isinstance(a, argparse._SubParsersAction):
            for ch in a._choices_actions:
                if ch.help:
                    ch.help = tr(ch.help)
            for _name, sub in a.choices.items():
                if getattr(sub, 'description', None):
                    sub.description = tr(sub.description)
                if getattr(sub, 'epilog', None):
                    sub.epilog = tr(sub.epilog)
                for sa in getattr(sub, '_actions', ()):  # noqa: SLF001
                    if getattr(sa, 'help', None) and sa.help is not argparse.SUPPRESS:
                        sa.help = tr(sa.help)
    return p
