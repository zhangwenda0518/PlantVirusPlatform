# -*- coding: utf-8 -*-
"""★ 对拍：我们的治理层 vs 上游 `MMPV-RNA/virome_phylo_pipeline/utils/metadata_governance.py`。

为什么这么验：移植最怕"看着像、口径悄悄不一样"。这里**直接 import 上游模块当参照**，
对同一批输入逐条比 `iso / precision / decimal / issue / normalized / usable`，
不一致就报出来 —— 差异要么是我改错，要么是**有意差异**（必须写明理由）。

用法: python tests/_check_govern_vs_upstream.py
退出码: 0 全一致（或差异都在白名单里）/ 1 有未解释的差异
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
# 外部参照仓（延伸基因组/MMPV-RNA 工作区）：默认仍是维护者本机路径，
# 可用环境变量 VP_EXT_ROOT 覆盖（指向 MMPV-RNA 工作区根）；
# 目录不存在时跳过（返回 0，不判失败）——上游只是对拍参照，不是运行依赖。
EXT_ROOT = os.environ.get('VP_EXT_ROOT', r'D:\桌面\延伸基因组\MMPV-RNA')
UP = os.path.join(EXT_ROOT, 'virome_phylo_pipeline')
if not os.path.isdir(UP):
    print(f'[SKIP] 上游参照仓不存在: {UP}'
          '（如已迁移，请设 VP_EXT_ROOT 指向 MMPV-RNA 工作区根）')
    sys.exit(0)
sys.path.insert(0, UP)

from Virus_Platform_Core import phylodyn_govern as og      # noqa: E402

PASS, FAIL, DIFF = [], [], []


def chk(cond, msg):
    (PASS if cond else FAIL).append(msg)
    print(('  ok   ' if cond else '  FAIL ') + msg, flush=True)


DATES = [
    # 常规
    '2013-05-01', '2013-05', '2013', '1999-12-31', '2000-02-29', '1900-01-01',
    # GenBank / 欧洲写法 / Mon-YYYY
    '01-May-2013', '5-Dec-2008', 'May-2013', '15/05/2013', '15.05.2013',
    # 小数年（≥4 位）与「年.月」（2 位）—— 消歧判据
    '2007.569473', '2018.49589', '2014.05', '2022.11',
    '2019.5', '2019.75',          # ← 与上游**有意不同**，见 INTENTIONAL
    # 占位符（形似合法）
    'Yyyy-Mm-Dd', 'yyyy', 'YYYY-MM', 'XX.XX N XXX.XX E', 'Country:Region',
    'unknown', 'NA', '', '未',
    # 非法 / 越界
    '2020-13-01', '2020-02-30', '1900-00-01', 'abcd', '13/13/2020', '3000-01-01',
]
LOCS = [
    'China:Ningxia', 'China, Ningxia, Yinchuan_AI', 'China, Unknown, Yinchuan_AI',
    'China, Unknown, Unknown_AI', 'Unknown', 'Unknown:Unknown', 'unknown',
    'Russia: Stavropol', 'Russia', 'USA', 'Germany: Bavaria, Munich',
    'XX.XX N XXX.XX E', 'Country:Region', '', 'China:Beijing', 'Macao, China',
]

# 有意差异白名单：输入 → 理由（只在"结果语义等价、只是标签/措辞不同"时登记）
INTENTIONAL = {
    # 两者都**拒绝**（iso=None），只是 issue 标签不同：我们区分"越界"与"认不出"，
    # 上游一律记 unparsable。我们的更细，且平台既有 `validate_date_str` 也用
    # "年份超范围"这种措辞 → 保留我们的标签。
    '3000-01-01': '超出 1900–2100：我们标 out_of_range（更细），上游标 unparsable；行为同为拒绝',
    # 上游只认「小数位 ≥4」为小数年（他们的数据源是 BioAider 的「年.月」写法）；
    # **本平台的数据源是 explorer 的 date 列**，那里 `2019.5` 就是小数年
    # （平台既有 `decimal_to_iso('2019.5') → 2019-07-02`，且有测试钉住）。
    # 照抄上游会把它读成"5 月"、偏半年 → 我们改判据为「带前导零的两位才是年.月」。
    '2019.5': '小数年位数判据与上游不同：我们的更宽（数据源不同），语义上我们的对',
    '2019.75': '同上（上游会把 2019.75 读成"年.75"→ 越界拒绝）',
}


def main():
    try:
        from utils import metadata_governance as up
    except Exception as e:                                          # noqa: BLE001
        print(f'⚠️ 无法 import 上游模块（{type(e).__name__}: {e}）→ 跳过对拍')
        return 0

    print('[1] 日期逐条对拍（iso / precision / decimal / issue）')
    n_same = 0
    for v in DATES:
        u = up.parse_date(v)
        o = og.parse_date(v)
        # 上游 issue 用 None 表示无问题
        ui, oi = (u.issue or ''), (o.issue or '')
        ui_eff = 'out_of_range' if ui in ('out_of_range',) else ui
        same = (u.iso == o.iso and u.precision == o.precision
                and (u.decimal is None) == (o.decimal is None)
                and (u.decimal is None or abs(u.decimal - o.decimal) < 1e-9)
                and (ui_eff or '') == (oi or ''))
        if not same and v in INTENTIONAL:
            same = True          # 有意差异：白名单登记过，**打印出来保持可见**
            print(f'  ok*  {v!r:18s} 有意差异（上游 {u.iso}/{u.decimal} vs '
                  f'我们 {o.iso}/{o.decimal}）：{INTENTIONAL[v]}')
        if same:
            n_same += 1
            print(f'  ok   {v!r:18s} iso={o.iso!r:14s} prec={o.precision:6s} '
                  f'dec={o.decimal}')
        else:
            DIFF.append((v, ui, oi, u.iso, o.iso, u.decimal, o.decimal))
            print(f'  ✗    {v!r:18s} 上游(iso={u.iso!r} prec={u.precision} '
                  f'dec={u.decimal} issue={ui!r}) vs 我们(iso={o.iso!r} '
                  f'prec={o.precision} dec={o.decimal} issue={oi!r})')
    chk(n_same == len(DATES), f'日期 {n_same}/{len(DATES)} 条与上游一致')

    print('\n[2] 地理逐条对拍（normalized / dropped_levels / issue）')
    n_same2 = 0
    for v in LOCS:
        u = up.parse_location(v)
        o = og.parse_location(v)
        ui, oi = (u.issue or ''), (o.issue or '')
        same = ((u.normalized or '') == (o.normalized or '')
                and sorted(u.dropped_levels or []) == sorted(o.dropped_levels or [])
                and ui == oi)
        if not same and v in INTENTIONAL:
            same = True
            print(f'  ok*  {v!r:34s} 白名单：{INTENTIONAL[v]}')
        if same:
            n_same2 += 1
            print(f'  ok   {v!r:34s} → {o.normalized!r:22s} '
                  f'issue={oi!r} usable={o.usable}')
        else:
            DIFF.append((v, ui, oi, u.normalized, o.normalized,
                         u.dropped_levels, o.dropped_levels))
            print(f'  ✗    {v!r:34s} 上游(norm={u.normalized!r} '
                  f'drop={u.dropped_levels} issue={ui!r}) vs '
                  f'我们(norm={o.normalized!r} drop={o.dropped_levels} issue={oi!r})')
    chk(n_same2 == len(LOCS), f'地理 {n_same2}/{len(LOCS)} 条与上游一致')

    print('\n[3] 不变量（无论与上游是否逐字一致，这几条必须成立）')
    chk(og.parse_date('2020-13-01').issue == 'out_of_range',
        '越界月份判 out_of_range（**不夹取**，否则 13 月会被当好数据放行）')
    chk(og.parse_date('2020-02-30').issue == 'out_of_range', '越界日期判 out_of_range')
    chk(og.parse_date('2007.569473').fmt == 'decimal_year'
        and og.parse_date('2014.05').precision == 'month',
        '小数年（≥4 位）与「年.月」（2 位）消歧正确')
    chk(og.parse_location('China, Unknown, Unknown_AI').usable is False,
        '只剩国家名 → usable=False（不赋国家中心点）')
    chk(og.parse_location('China:Ningxia').usable is True
        and og.parse_location('China:Ningxia').normalized == 'China, Ningxia',
        '冒号也拆（修「只拆逗号」的缺口）')
    chk(og.decimal_year('2013', mode='mid') > og.decimal_year('2013', mode='start'),
        'mode=mid 在「只给年」时晚于 mode=start（口径显式，不藏在实现里）')
    chk(abs(og.decimal_year('2013-05-01') - (2013 + (4 + 0 / 31) / 12)) < 1e-12,
        '完整日期口径 = yr+(mo-1+(dy-1)/dim)/12（与上游一致）')
    chk(og.is_placeholder('XX.XX N XXX.XX E') and not og.is_placeholder('China:Ningxia'),
        '占位符识别不误伤真实地理（上游曾用泛化正则误伤 72 条）')

    print('\n' + '=' * 52)
    if DIFF:
        print(f'与上游存在 {len(DIFF)} 处差异：')
        for d in DIFF:
            print('  -', d)
    if FAIL:
        print(f'FAILED: {len(FAIL)} 项')
        for m in FAIL:
            print('  - ' + m)
        return 1
    print(f'GOVERN CROSSCHECK PASSED（{len(PASS)} 项；与上游日期/地理逐条一致）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
