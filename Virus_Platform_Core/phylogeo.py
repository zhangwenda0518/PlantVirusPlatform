# -*- coding: utf-8 -*-
"""系统地理分析（本地实现）——地理/时间元数据驱动的迁移重构 + 分子钟定年。

对标 39.106.101.94 远端 Phylogeography（Dash 应用）的本地化版本：
  序列 FASTA（MAFFT 比对）→ 建树（复用平台 NJ / FastTree）→
  Fitch 最大简约把「地理区划」性状重构到内部节点 →
  输出：迁移矩阵（区间转移计数）、带区划标注的树、逐叶状态表。

本模块另外提供三个**可选**统计层（默认关闭，不拖慢常规运行）：
  · `region_permutation_test` 区域随机化检验（简约法）—— 区划是否真有信号
  · `rssp_bootstrap` 根状态不确定性 —— 给迁移数配 bootstrap 置信区间
  · `lsd2_dating` LSD2 分子钟定年 —— 速率 + tMRCA + 定年树（真正的时间标度）
`root_to_tip` 是免依赖的时间信号初筛（根到尾回归），与 LSD2 不同层次：
前者只回答「有没有时间信号」，后者给绝对年代。

元数据来源：CSV（列 seq_id + 区划列）或 FASTA 头 `>acc|区域|年份`。
⚠️ 头回退解析有个陷阱：头若写成 `acc|<整条描述>`（平台示例集就是这样），
「区划」会变成整条描述 → 每条样本一个独立区划、迁移矩阵与置换检验全部失去意义。
analyze() 会在这种情况下**显式告警**（见返回值的 warnings）。
"""

import json
import math
import os
import re
import shutil
import subprocess
from collections import Counter, OrderedDict


# ---------------------------------------------------------------- 地理坐标（迁移弧线地图）

# 内置质心表：区域/国家（英/中别名）→ (纬度, 经度)。
# 用途：迁移弧线地图的**兜底坐标**（用户坐标表优先，见 merge_coords）。
# 口径：国家取首都附近的示意点、大洲取人口/采样重心一类的**近似质心**，
# 只用于把区划落到世界地图上，不要当作精确地理中心。
# 未收录的名字请在「区域坐标表」里显式给坐标（也兼容 VirPhyKit
# Spatial.txt 的逐序列坐标，逐条取该区划的质心）。
BUILTIN_REGION_COORDS = {
    # 大洲 / 泛区域
    'Asia': (34.0, 100.0), '亚洲': (34.0, 100.0),
    'Europe': (50.0, 10.0), '欧洲': (50.0, 10.0),
    'Africa': (0.0, 25.0), '非洲': (0.0, 25.0),
    'North America': (45.0, -100.0), '北美洲': (45.0, -100.0),
    'South America': (-15.0, -60.0), '南美洲': (-15.0, -60.0),
    'Americas': (10.0, -80.0), '美洲': (10.0, -80.0),
    'Oceania': (-25.0, 140.0), '大洋洲': (-25.0, 140.0),
    'East Asia': (35.0, 115.0), '东亚': (35.0, 115.0),
    'Southeast Asia': (10.0, 105.0), '东南亚': (10.0, 105.0),
    'South Asia': (20.0, 78.0), '南亚': (20.0, 78.0),
    'Central Asia': (45.0, 65.0), '中亚': (45.0, 65.0),
    'West Asia': (30.0, 45.0), '西亚': (30.0, 45.0),
    'Middle East': (30.0, 45.0), '中东': (30.0, 45.0),
    'Eurasia': (45.0, 60.0), '欧亚大陆': (45.0, 60.0),
    'Latin America': (0.0, -70.0), '拉美': (0.0, -70.0),
    'Caribbean': (18.0, -75.0), '加勒比': (18.0, -75.0),
    'Other': (30.0, 80.0), '其他': (30.0, 80.0),
    # 东亚
    'China': (35.0, 105.0), '中国': (35.0, 105.0),
    'Japan': (36.0, 138.0), '日本': (36.0, 138.0),
    'Korea': (37.0, 127.5), '韩国': (37.0, 127.5), 'South Korea': (37.0, 127.5),
    'North Korea': (40.0, 127.0), '朝鲜': (40.0, 127.0),
    'Taiwan': (23.7, 121.0), '台湾': (23.7, 121.0),
    'Beijing': (39.9, 116.4), '北京': (39.9, 116.4),
    'Ningxia': (37.3, 106.0), '宁夏': (37.3, 106.0),
    'Neimenggu': (44.0, 113.0), '内蒙古': (44.0, 113.0),
    # 南亚 / 东南亚
    'India': (20.0, 78.0), '印度': (20.0, 78.0),
    'Pakistan': (30.0, 70.0), '巴基斯坦': (30.0, 70.0),
    'Bangladesh': (24.0, 90.0), '孟加拉国': (24.0, 90.0),
    'Sri Lanka': (7.0, 81.0), '斯里兰卡': (7.0, 81.0),
    'Nepal': (28.0, 84.0), '尼泊尔': (28.0, 84.0),
    'Indonesia': (-2.0, 118.0), '印度尼西亚': (-2.0, 118.0), '印尼': (-2.0, 118.0),
    'Thailand': (15.0, 101.0), '泰国': (15.0, 101.0),
    'Vietnam': (16.0, 108.0), '越南': (16.0, 108.0),
    'Philippines': (13.0, 122.0), '菲律宾': (13.0, 122.0),
    'Malaysia': (4.0, 102.0), '马来西亚': (4.0, 102.0),
    'Singapore': (1.35, 103.8), '新加坡': (1.35, 103.8),
    # 北美
    'USA': (39.0, -98.0), 'United States': (39.0, -98.0), '美国': (39.0, -98.0),
    'Canada': (56.0, -106.0), '加拿大': (56.0, -106.0),
    'Mexico': (23.0, -102.0), '墨西哥': (23.0, -102.0),
    'Cuba': (22.0, -80.0), '古巴': (22.0, -80.0),
    # 拉美
    'Brazil': (-10.0, -52.0), '巴西': (-10.0, -52.0),
    'Argentina': (-35.0, -65.0), '阿根廷': (-35.0, -65.0),
    'Chile': (-33.0, -71.0), '智利': (-33.0, -71.0),
    'Peru': (-10.0, -76.0), '秘鲁': (-10.0, -76.0),
    'Colombia': (4.0, -72.0), '哥伦比亚': (4.0, -72.0),
    'Venezuela': (7.0, -66.0), '委内瑞拉': (7.0, -66.0),
    'Ecuador': (-1.8, -78.0), '厄瓜多尔': (-1.8, -78.0),
    'Bolivia': (-17.0, -65.0), '玻利维亚': (-17.0, -65.0),
    'Uruguay': (-33.0, -56.0), '乌拉圭': (-33.0, -56.0),
    'Paraguay': (-23.0, -58.0), '巴拉圭': (-23.0, -58.0),
    # 欧洲
    'United Kingdom': (54.0, -2.0), 'UK': (54.0, -2.0), '英国': (54.0, -2.0),
    'France': (46.0, 2.0), '法国': (46.0, 2.0),
    'Germany': (51.0, 10.0), '德国': (51.0, 10.0),
    'Italy': (42.5, 12.5), '意大利': (42.5, 12.5),
    'Spain': (40.0, -4.0), '西班牙': (40.0, -4.0),
    'Portugal': (39.5, -8.0), '葡萄牙': (39.5, -8.0),
    'Netherlands': (52.2, 5.3), '荷兰': (52.2, 5.3),
    'Belgium': (50.5, 4.5), '比利时': (50.5, 4.5),
    'Poland': (52.0, 19.0), '波兰': (52.0, 19.0),
    'Czech Republic': (50.0, 15.0), '捷克': (50.0, 15.0),
    'Austria': (47.5, 14.5), '奥地利': (47.5, 14.5),
    'Switzerland': (47.0, 8.0), '瑞士': (47.0, 8.0),
    'Sweden': (62.0, 15.0), '瑞典': (62.0, 15.0),
    'Norway': (62.0, 9.0), '挪威': (62.0, 9.0),
    'Finland': (64.0, 26.0), '芬兰': (64.0, 26.0),
    'Denmark': (56.0, 10.0), '丹麦': (56.0, 10.0),
    'Hungary': (47.0, 19.0), '匈牙利': (47.0, 19.0),
    'Romania': (46.0, 25.0), '罗马尼亚': (46.0, 25.0),
    'Serbia': (44.0, 21.0), '塞尔维亚': (44.0, 21.0),
    'Greece': (39.0, 22.0), '希腊': (39.0, 22.0),
    'Bulgaria': (42.7, 25.5), '保加利亚': (42.7, 25.5),
    'Croatia': (45.1, 15.2), '克罗地亚': (45.1, 15.2),
    'Ireland': (53.0, -8.0), '爱尔兰': (53.0, -8.0),
    'Russia': (60.0, 90.0), '俄罗斯': (60.0, 90.0),
    'Ukraine': (49.0, 32.0), '乌克兰': (49.0, 32.0),
    # 西亚 / 北非 / 非洲
    'Turkey': (39.0, 35.0), '土耳其': (39.0, 35.0),
    'Iran': (32.0, 53.0), '伊朗': (32.0, 53.0),
    'Iraq': (33.0, 44.0), '伊拉克': (33.0, 44.0),
    'Israel': (31.5, 35.0), '以色列': (31.5, 35.0),
    'Saudi Arabia': (24.0, 45.0), '沙特阿拉伯': (24.0, 45.0),
    'Egypt': (26.0, 30.0), '埃及': (26.0, 30.0),
    'Morocco': (32.0, -6.0), '摩洛哥': (32.0, -6.0),
    'Kenya': (0.0, 38.0), '肯尼亚': (0.0, 38.0),
    'Nigeria': (10.0, 8.0), '尼日利亚': (10.0, 8.0),
    'Ethiopia': (9.0, 40.0), '埃塞俄比亚': (9.0, 40.0),
    'South Africa': (-29.0, 24.0), '南非': (-29.0, 24.0),
    # 大洋洲
    'Australia': (-25.0, 134.0), '澳大利亚': (-25.0, 134.0), '澳洲': (-25.0, 134.0),
    'New Zealand': (-42.0, 173.0), '新西兰': (-42.0, 173.0),
    # 补充国家（2026-09-16：Explorer 公共库实测未命中项，按行数排序补齐；
    # 坐标同样是国家/岛屿地理中心附近的**近似**值，只用于底图定位与弧线）
    'Tunisia': (34.0, 9.0), '突尼斯': (34.0, 9.0),
    'Ghana': (8.0, -1.0), '加纳': (8.0, -1.0),
    'Tanzania': (-6.0, 35.0), '坦桑尼亚': (-6.0, 35.0),
    'Uganda': (1.0, 32.0), '乌干达': (1.0, 32.0),
    'Sudan': (15.0, 30.0), '苏丹': (15.0, 30.0),
    'Senegal': (14.0, -14.0), '塞内加尔': (14.0, -14.0),
    'Cameroon': (6.0, 12.0), '喀麦隆': (6.0, 12.0),
    "Cote d'Ivoire": (7.5, -5.5), '科特迪瓦': (7.5, -5.5),
    'Ivory Coast': (7.5, -5.5),
    'DR Congo': (-3.0, 23.0), 'Democratic Republic of the Congo': (-3.0, 23.0),
    '刚果民主共和国': (-3.0, 23.0),
    'Republic of the Congo': (-1.0, 15.0),
    'Madagascar': (-19.0, 47.0), '马达加斯加': (-19.0, 47.0),
    'Reunion': (-21.1, 55.5), '留尼汪': (-21.1, 55.5),
    'Zimbabwe': (-19.0, 30.0), '津巴布韦': (-19.0, 30.0),
    'Zambia': (-14.0, 27.0), '赞比亚': (-14.0, 27.0),
    'Mozambique': (-18.0, 35.0), '莫桑比克': (-18.0, 35.0),
    'Namibia': (-22.0, 17.0), '纳米比亚': (-22.0, 17.0),
    'Botswana': (-22.0, 24.0), '博茨瓦纳': (-22.0, 24.0),
    'Mali': (17.0, -4.0), '马里': (17.0, -4.0),
    'Niger': (17.0, 8.0), '尼日尔': (17.0, 8.0),
    'Chad': (15.0, 19.0), '乍得': (15.0, 19.0),
    'Guinea': (11.0, -11.0), '几内亚': (11.0, -11.0),
    'Sierra Leone': (8.5, -11.8), '塞拉利昂': (8.5, -11.8),
    'Liberia': (6.5, -9.5), '利比里亚': (6.5, -9.5),
    'Togo': (8.0, 1.2), '多哥': (8.0, 1.2),
    'Benin': (9.3, 2.3), '贝宁': (9.3, 2.3),
    'Burkina Faso': (12.3, -1.5), '布基纳法索': (12.3, -1.5),
    'Gabon': (-1.0, 11.8), '加蓬': (-1.0, 11.8),
    'Malawi': (-13.5, 34.0), '马拉维': (-13.5, 34.0),
    'Rwanda': (-2.0, 29.9), '卢旺达': (-2.0, 29.9),
    'Burundi': (-3.4, 29.9), '布隆迪': (-3.4, 29.9),
    'Somalia': (5.0, 46.0), '索马里': (5.0, 46.0),
    'Libya': (27.0, 17.0), '利比亚': (27.0, 17.0),
    'Algeria': (28.0, 3.0), '阿尔及利亚': (28.0, 3.0),
    'Angola': (-12.0, 18.0), '安哥拉': (-12.0, 18.0),
    'Mauritius': (-20.3, 57.6), '毛里求斯': (-20.3, 57.6),
    'Seychelles': (-4.7, 55.5),
    'Cape Verde': (16.0, -24.0),
    'Eritrea': (15.3, 39.0),
    'Djibouti': (11.8, 42.6),
    'Jordan': (31.0, 36.0), '约旦': (31.0, 36.0),
    'Lebanon': (33.9, 35.9), '黎巴嫩': (33.9, 35.9),
    'Syria': (35.0, 38.0), '叙利亚': (35.0, 38.0),
    'Oman': (21.0, 57.0), '阿曼': (21.0, 57.0),
    'Yemen': (15.5, 48.0), '也门': (15.5, 48.0),
    'Kuwait': (29.4, 47.7), '科威特': (29.4, 47.7),
    'Qatar': (25.3, 51.2), '卡塔尔': (25.3, 51.2),
    'United Arab Emirates': (24.0, 54.0), '阿联酋': (24.0, 54.0),
    'Bahrain': (26.0, 50.5), '巴林': (26.0, 50.5),
    'Afghanistan': (33.9, 66.0), '阿富汗': (33.9, 66.0),
    'Georgia': (42.0, 43.5), '格鲁吉亚': (42.0, 43.5),
    'Armenia': (40.1, 45.0), '亚美尼亚': (40.1, 45.0),
    'Azerbaijan': (40.3, 47.9), '阿塞拜疆': (40.3, 47.9),
    'Cyprus': (35.1, 33.2), '塞浦路斯': (35.1, 33.2),
    'Slovenia': (46.1, 14.8), '斯洛文尼亚': (46.1, 14.8),
    'Slovakia': (48.7, 19.7), '斯洛伐克': (48.7, 19.7),
    'Lithuania': (55.2, 23.9), '立陶宛': (55.2, 23.9),
    'Latvia': (56.9, 24.6), '拉脱维亚': (56.9, 24.6),
    'Estonia': (58.6, 25.0), '爱沙尼亚': (58.6, 25.0),
    'Belarus': (53.7, 28.0), '白俄罗斯': (53.7, 28.0),
    'Moldova': (47.4, 28.4), '摩尔多瓦': (47.4, 28.4),
    'Albania': (41.2, 20.2), '阿尔巴尼亚': (41.2, 20.2),
    'Bosnia and Herzegovina': (44.2, 17.8),
    'North Macedonia': (41.6, 21.7),
    'Montenegro': (42.7, 19.4),
    'Iceland': (65.0, -18.0), '冰岛': (65.0, -18.0),
    'Luxembourg': (49.8, 6.1),
    'Malta': (35.9, 14.4),
    'Myanmar': (21.9, 96.0), '缅甸': (21.9, 96.0),
    'Cambodia': (12.6, 105.0), '柬埔寨': (12.6, 105.0),
    'Laos': (19.9, 102.5), '老挝': (19.9, 102.5),
    'Mongolia': (46.9, 103.8), '蒙古': (46.9, 103.8),
    'Kazakhstan': (48.0, 67.0), '哈萨克斯坦': (48.0, 67.0),
    'Uzbekistan': (41.4, 64.6), '乌兹别克斯坦': (41.4, 64.6),
    'Turkmenistan': (39.0, 59.5),
    'Kyrgyzstan': (41.2, 74.8),
    'Tajikistan': (38.9, 71.3),
    'Bhutan': (27.5, 90.4),
    'Maldives': (3.2, 73.2),
    'Brunei': (4.5, 114.7),
    'Guatemala': (15.8, -90.2), '危地马拉': (15.8, -90.2),
    'Honduras': (14.8, -86.2),
    'El Salvador': (13.8, -88.9),
    'Nicaragua': (12.9, -85.2),
    'Costa Rica': (9.7, -83.8), '哥斯达黎加': (9.7, -83.8),
    'Panama': (8.5, -80.8), '巴拿马': (8.5, -80.8),
    'Dominican Republic': (19.0, -70.2),
    'Haiti': (19.0, -72.3),
    'Jamaica': (18.1, -77.3), '牙买加': (18.1, -77.3),
    'Trinidad and Tobago': (10.5, -61.3),
    'Belize': (17.2, -88.5),
    'Guyana': (5.0, -58.9),
    'Suriname': (4.0, -56.0),
    'French Guiana': (4.0, -53.0), '法属圭亚那': (4.0, -53.0),
    'Bermuda': (32.3, -64.8),
    'Puerto Rico': (18.2, -66.5),
    'Papua New Guinea': (-6.0, 145.0),
    'Fiji': (-18.0, 178.0), '斐济': (-18.0, 178.0),
    'Samoa': (-13.8, -172.1),
    'Tonga': (-20.0, -175.0),
    'New Caledonia': (-21.3, 165.5), '新喀里多尼亚': (-21.3, 165.5),
    'Guam': (13.4, 144.8),
    'French Polynesia': (-17.6, -149.4),
    'Solomon Islands': (-9.6, 160.2),
    'Vanuatu': (-15.4, 166.9),
    # ---- 补录：Explorer 全库扫描后仍缺质心的主权国家/属地（2026-09-16）----
    'Central African Republic': (6.6, 20.9),
    'Comoros': (-11.6, 43.3),
    'Guadeloupe': (16.3, -61.6),
    'Mayotte': (-12.8, 45.2),
    'Martinique': (14.6, -61.0),
    'Grenada': (12.1, -61.7),
    'Barbados': (13.2, -59.5),
    'South Sudan': (7.9, 30.2),
    'Gambia': (13.4, -15.4),
    'Sao Tome and Principe': (0.2, 6.6),
    'Lesotho': (-29.6, 28.2),
    'Palau': (7.5, 134.6),
    'Equatorial Guinea': (1.7, 10.3),
    'Saint Vincent and the Grenadines': (13.3, -61.2),
    'Saint Kitts and Nevis': (17.4, -62.8),
    'Antigua and Barbuda': (17.1, -61.8),
    'Niue': (-19.1, -169.9),
    'Cook Islands': (-21.2, -159.8),
    'Wallis and Futuna': (-13.3, -176.2),
    'Norfolk Island': (-29.0, 168.0),
    'Antarctica': (-75.0, 0.0),
    'Eswatini': (-26.5, 31.5),
    'Timor-Leste': (-8.9, 125.7),
    'Palestine': (31.9, 35.2),
    'Federated States of Micronesia': (7.4, 150.5),
    'Kosovo': (42.6, 20.9),
}

# 区划名别名：库里 / 服务器表里常见的旧称、简称、异写 → 质心表的键。
# 只做「同一地点换个写法」的归一，不含糊地跨地区合并；历史实体（Yugoslavia
# 一类）按继承国质心落点，并在结果里以原名显示，不会假装是当代采样点。
REGION_ALIASES = {
    'Czechia': 'Czech Republic',
    'Macedonia': 'North Macedonia',
    'The former Yugoslav Republic of Macedonia': 'North Macedonia',
    'Swaziland': 'Eswatini',
    'East Timor': 'Timor-Leste',
    'Micronesia': 'Federated States of Micronesia',
    'West Bank': 'Palestine',
    'State of Palestine': 'Palestine',
    'Gaza Strip': 'Palestine',
    'Yugoslavia': 'Serbia',
    'Serbia and Montenegro': 'Serbia',
    'Czechoslovakia': 'Czech Republic',
    'USSR': 'Russia',
    'Russian Federation': 'Russia',
    'Republic of Korea': 'South Korea',
    'Korea': 'South Korea',
    'Viet Nam': 'Vietnam',
    'United States of America': 'United States',
}

# ⚠️ 2026-10-03 删：原常量 `NON_PLACE_REGIONS = {'Unknown', 'missing', …}` 全仓
#    **无任何代码引用**（M3 采样分布图早已改用别的口径），只留下"排除伪区划"
#    的注释性常量 → 死代码，删除。判"没有区划信息"的统一实现见
#    `phylodyn_kit.is_no_region`（与 `phylodyn_govern._is_unknown_level` 对齐）。

_LAT_ALIASES = ('latitude', 'lat', '纬度', 'y')
_LON_ALIASES = ('longitude', 'lon', 'lng', 'long', '经度', 'x')
_NAME_ALIASES = ('region', 'name', 'locality', 'location', 'state', 'area',
                 'province', 'country', '地区', '区划', '区域', '地名', '省',
                 # 逐样本坐标表（每序列一行）常用第一列名 —— 不认它们的话，
                 # 表头会被当数据行开头处理，整表错位
                 'sample', 'seq_id', 'seqid', 'accession')


def _read_coord_rows(path):
    """坐标表解析内核（`parse_coord_table` / `parse_coord_rows` 共用）。

    兼容三种常见形态（自动识别；分隔符取 tab/逗号/分号中最多者，都没有则按空白切）：
      1. 表头式（推荐）：首行含 Region/Latitude/Longitude（或 name/lat/lon 等别名）；
      2. VirPhyKit `Spatial.txt` 式：同 1（Region/Latitude/Longitude + 其它列），
         **逐序列**多行同区划；
      3. 无表头式：前三列 = 区划 / 纬度 / 经度。
    坏行**不静默**：缺列 / 空名 / 非数值 / 坐标越界（|lat|>90 或 |lon|>180）
    逐行记入 skipped（人可读原因 + 行号），调用方必须展示出来。
    返回 (rows[(name, lat, lon, 行号)], skipped)。
    """
    import csv
    try:
        with open(path, encoding='utf-8-sig', errors='replace', newline='') as f:
            text = f.read()
    except OSError as e:
        return [], [f'坐标表读取失败：{e}']
    lines = [ln for ln in text.splitlines()
             if ln.strip() and not ln.lstrip().startswith('#')]
    if not lines:
        return [], ['坐标表是空的（没有可用数据行）']
    counts = {d: text.count(d) for d in ('\t', ',', ';')}
    best = max(counts, key=counts.get)
    if counts[best] > 0:
        def _split(ln):
            return next(csv.reader([ln], delimiter=best))
    else:
        def _split(ln):
            return ln.split()
    rows = [_split(ln) for ln in lines]
    header = [c.strip().lower() for c in rows[0]]
    i_name = i_lat = i_lon = None
    for i, h in enumerate(header):
        if i_name is None and h in _NAME_ALIASES:
            i_name = i
        if i_lat is None and h in _LAT_ALIASES:
            i_lat = i
        if i_lon is None and h in _LON_ALIASES:
            i_lon = i
    if i_name is not None and i_lat is not None and i_lon is not None:
        data, off = rows[1:], 1
    else:
        i_name, i_lat, i_lon = 0, 1, 2
        data, off = rows, 0
    need = max(i_name, i_lat, i_lon) + 1
    out, skipped = [], []
    for k, row in enumerate(data, off + 1):
        if len(row) < need:
            skipped.append(f'第{k}行：列数不足（需要区划/纬度/经度三列）')
            continue
        # 头一个单元格常残留 BOM（utf-8 无 sig 打开的机器上）——不算区划名一部分
        name = row[i_name].strip().lstrip('\ufeff').strip()
        if not name:
            skipped.append(f'第{k}行：区划名为空')
            continue
        try:
            lat = float(row[i_lat])
            lon = float(row[i_lon])
        except ValueError:
            skipped.append(
                f'第{k}行：纬度/经度不是数值（{row[i_lat]!r} / {row[i_lon]!r}）')
            continue
        if (not math.isfinite(lat) or not math.isfinite(lon)
                or not (-90.0 <= lat <= 90.0) or not (-180.0 <= lon <= 180.0)):
            skipped.append(f'第{k}行：坐标越界或非有限值（{lat}, {lon}）')
            continue
        out.append((name, lat, lon, k))
    return out, skipped


def _circular_mean_lon(lons):
    """经度序列 → **环形均值**（单位向量合成：mean(cosθ)/mean(sinθ) → atan2）。

    ⚠️ 2026-10-03 修：经度是**环形量**，跨 ±180° 日界线的锚点
    （178 / −179.5 / 179 / −178）用普通均值会翻到地球另一侧（≈0）→ 弧线横穿
    整张地图；环形均值给 ≈179.9，仍贴近日界线。纬度不是环形量（−90..90 无环），
    故仍用普通均值/中位数。返回 [−180, 180] 内的度数。
    退化情形（正负完全抵消，如 0 与 180 对点）向量长度≈0、方向无意义 →
    退回普通均值并保持有限值（不返回 NaN）。
    """
    vals = [float(v) for v in lons]
    if not vals:
        return None
    s = sum(math.sin(math.radians(v)) for v in vals)
    c = sum(math.cos(math.radians(v)) for v in vals)
    if abs(s) < 1e-12 and abs(c) < 1e-12:
        return sum(vals) / len(vals)
    return math.degrees(math.atan2(s, c))


def coord_table_centroid(rows):
    """逐行坐标条目 → {名称: [lat, lon]}（同名多行取均值 = 质心）。

    两种逐行条目都收：`_read_coord_rows` 的 (name, lat, lon, 行号) 元组，
    与 `parse_coord_rows` 的 {name, lat, lon, line} 字典。

    ⚠️ 纬度普通均值；**经度改环形均值**（`_circular_mean_lon`）—— 原实现逐轴
    均值，跨 ±180° 的同区划锚点会翻到地球另一侧（审计实测
    178/−179.5/179/−178 → lon ≈ 0）。A2 的 `merge_coords → parse_coord_table`
    弧线落点链路都吃这里，修一处即全链一致。
    """
    acc = {}
    for r in rows:
        if isinstance(r, dict):
            name, lat, lon = r.get('name'), r.get('lat'), r.get('lon')
        else:
            name, lat, lon = r[0], r[1], r[2]
        a = acc.setdefault(name, ([], []))
        a[0].append(float(lat))
        a[1].append(float(lon))
    return {nm: [round(sum(la) / len(la), 5),
                 round(_circular_mean_lon(lo), 5)]
            for nm, (la, lo) in acc.items()}


def parse_coord_table(path):
    """区域坐标表 → ({区划: [lat, lon]}, skipped)。

    兼容形态见 `_read_coord_rows`；**同名多行按质心（均值）合并**
    （VirPhyKit Spatial.txt 每序列一行 —— 逐样本点层请用 `parse_coord_rows`）。
    """
    rows, skipped = _read_coord_rows(path)
    return coord_table_centroid(rows), skipped


def parse_coord_rows(path):
    """坐标表**逐行**条目（不做区划聚合）→ ([{name, lat, lon, line}], skipped)。

    与 `parse_coord_table` 的区别：那个把同名行聚合成「区划质心」供弧线地图用；
    这个保留每一行 —— 逐样本坐标表（每序列一行 lat/lon）要的就是它，
    供「样本点层」与「区划位置 = 样本坐标中位数」使用。
    """
    rows, skipped = _read_coord_rows(path)
    return [{'name': nm, 'lat': lat, 'lon': lon, 'line': k}
            for nm, lat, lon, k in rows], skipped


def match_coord(region, coords):
    """在坐标表里找区划坐标：先精确，再去空白后不分大小写。

    仍找不到时再查 `REGION_ALIASES`（旧称/简称/异写）——**只对内置质心表生效**，
    用户自己给的坐标表不做别名替换（用户写了什么就是什么，不许替他改）。
    找不到返回 None。
    """
    if not region or not coords:
        return None
    if region in coords:
        return coords[region]
    low = {str(k).strip().lower(): v for k, v in coords.items()}
    hit = low.get(str(region).strip().lower())
    if hit is not None:
        return hit
    if coords is BUILTIN_REGION_COORDS:
        alias = REGION_ALIASES.get(str(region).strip())
        if alias:
            return BUILTIN_REGION_COORDS.get(alias)
    return None


def region_name_candidates(region):
    """区划名的候选写法（**精度从高到低**），供坐标匹配逐级兜底。

    为什么需要：本平台元数据的地点列惯用「国家: 子区」写法
    （`Russia: Stavropol`、`China: Ningxia`），而内置质心表只到国家一级 ——
    不做前缀兜底时这些区划会**全部**落进 `coord_missing`，迁移弧线地图近乎空白
    （实测 10/10 条真实地点都是"全名不命中、国家命中"）。
    候选顺序保证**先精确后粗**：① 原名 ② 「国家: 子区」的国家部分 ③ 去括号后的名字。
    兜底命中在 `merge_coords` 的 src 里标 `user_prefix` / `builtin_prefix` ——
    **不许当成精确命中**：落点是国家（或粗一级名字）的近似质心，比数据说的粗。
    """
    r = str(region or '').strip()
    out = [r] if r else []
    for sep in (':', '：'):
        if sep in r:
            head = r.split(sep, 1)[0].strip()
            if head and head not in out:
                out.append(head)
    base = re.sub(r'[（(][^)）]*[)）]', '', r).strip()
    if base and base not in out:
        out.append(base)
    return out


def merge_coords(regions, user_coords):
    """把「用户坐标表 + 内置质心表」合到分析结果的区划集合上。

    返回 (coords, src, missing)：
      coords  {区划: [lat, lon]}（用户表优先，其次内置质心表）
      src     {区划: 'user' | 'builtin' | 'user_prefix' | 'builtin_prefix'}
              —— 图上要有出处，别把近似质心当用户数据；`*_prefix` 表示是
                 退到「国家: 子区」的**国家部分**才命中的，落点比数据说的粗一级
                 （候选顺序见 `region_name_candidates`）。
      missing 没有任何坐标的区划列表 —— 前端必须显式提示，不许静默少画点。
    """
    coords, src, missing = {}, {}, []
    for r in regions:
        cands = region_name_candidates(r)
        hit = None
        # 用户表优先：**先把自己给的坐标表按候选逐级试完**，再退回内置表 ——
        # 否则「用户表只有国家、内置表恰好有同名子区」会拿内置的粗坐标顶掉
        # 用户明确给的坐标，等于悄悄改了用户的意思。
        for label, table in (('user', user_coords),
                             ('builtin', BUILTIN_REGION_COORDS)):
            for i, cand in enumerate(cands):
                c = match_coord(cand, table)
                if c is not None:
                    hit = ([float(c[0]), float(c[1])],
                           label if i == 0 else label + '_prefix')
                    break
            if hit is not None:
                break
        if hit is None:
            missing.append(r)
            continue
        coords[r], src[r] = hit
    return coords, src, missing


# ---------------------------------------------------- 逐样本经纬度（点层 / 距离）

_EARTH_R_KM = 6371.0088


def haversine_km(a, b):
    """两经纬度点 [lat, lon] 的大圆距离（km；地球平均半径 6371.0088 km）。

    为什么不用「经纬度差的平方和」：度数欧氏距离在高纬会把经度方向严重拉伸
    （60°N 处 1° 经度只有约 55.8 km，赤道是 111.3 km）——**跨纬度比较会失真**。
    迁移事件的公里级分类（Direct/Indirect/Distant）与任何距离展示都走这个函数，
    口径只有一个。
    """
    lat1, lon1 = math.radians(float(a[0])), math.radians(float(a[1]))
    lat2, lon2 = math.radians(float(b[0])), math.radians(float(b[1]))
    h = (math.sin((lat2 - lat1) / 2.0) ** 2
         + math.cos(lat1) * math.cos(lat2)
         * math.sin((lon2 - lon1) / 2.0) ** 2)
    return 2.0 * _EARTH_R_KM * math.asin(min(1.0, math.sqrt(h)))


def _pick_col(keys, aliases):
    """在列名集合里按别名找列（大小写/首尾空白容错）。找不到 → None。"""
    for k in keys:
        if str(k).strip().lower() in aliases:
            return k
    return None


def _coord_pair(lat_raw, lon_raw):
    """两单元格 → [lat, lon]（5 位小数）；非数值/越界/非有限 → None。"""
    try:
        lat = float(str(lat_raw).strip())
        lon = float(str(lon_raw).strip())
    except (TypeError, ValueError):
        return None
    if (not math.isfinite(lat) or not math.isfinite(lon)
            or not (-90.0 <= lat <= 90.0) or not (-180.0 <= lon <= 180.0)):
        return None
    return [round(lat, 5), round(lon, 5)]


def match_sample_coords(leaf_names, safe2header, meta, coord_rows):
    """逐样本经纬度：元数据 lat/lon 列 → 坐标表**逐行**条目。

    匹配口径与 `match_states_and_years` **逐步一致**（净化名→原头→前缀/包含、
    元数据按插入序取第一个命中键）——两处若各写一套，会出现「区划命中的是 A 行、
    坐标命中的是 B 行」这种自相矛盾的结果。

    返回 dict：
      points  {叶名: [lat, lon]}
      src     {叶名: 'meta' | 'coords'}（图上/表里如实标出处）
      n_meta / n_table / n_bad
              （n_bad＝元数据里经纬度列取到了值但不是有效数值/越界的条数，
                坏值不静默——由调用方报出）
    """
    points, src = {}, {}
    n_meta = n_table = n_bad = 0
    items = list(meta.items()) if meta else []
    for n in leaf_names:
        h = safe2header.get(n)
        if h is None:
            h = next((hh for k, hh in safe2header.items()
                      if n.startswith(k) or k in n), n)
        pt, hit = None, ''
        for k, v in items:
            if k and (h.startswith(k) or k in h or n.startswith(k)):
                lk = _pick_col(v.keys(), _LAT_ALIASES)
                ok_ = _pick_col(v.keys(), _LON_ALIASES)
                if lk and ok_ and str(v.get(lk) or '').strip() \
                        and str(v.get(ok_) or '').strip():
                    pt = _coord_pair(v.get(lk), v.get(ok_))
                    if pt is None:
                        n_bad += 1
                break                       # 与区划匹配同口径：取第一个命中行
        if pt is not None:
            hit = 'meta'
        elif coord_rows:
            for r in coord_rows:
                nm = r['name']
                if (nm == n or nm == h or h.startswith(nm) or nm in h
                        or n.startswith(nm)):
                    pt = [r['lat'], r['lon']]
                    hit = 'coords'
                    break
        if pt is not None:
            points[n], src[n] = pt, hit
            if hit == 'meta':
                n_meta += 1
            else:
                n_table += 1
    return {'points': points, 'src': src, 'n_meta': n_meta,
            'n_table': n_table, 'n_bad': n_bad}


def resolve_region_positions(regions, states, sample_points, user_coords):
    """区划落点优先级：**样本坐标（纬度中位数 / 经度环形均值）> 用户坐标表 > 内置质心表**。

    样本坐标口径：该区划下所有有坐标样本的**纬度中位数** + **经度环形均值**
    （环形量不能逐轴平均，跨 ±180° 日界线会翻到地球另一侧 —— 2026-10-03 修，
    见 `_circular_mean_lon`）。比「区划质心」更贴数据（样本散布时，质心会把
    弧线端点拖到没有样本的位置）。
    回退部分直接复用 `merge_coords`（用户表优先 → 内置表 → missing），
    保持与既有口径完全一致。

    返回 (coords, src{'samples'|'user'|'builtin'|'user_prefix'|'builtin_prefix'},
          missing, n_samples_used)。
    """
    by_region = {}
    for leaf, st in (states or {}).items():
        pt = (sample_points or {}).get(leaf)
        if pt and st and st != 'Unknown':
            by_region.setdefault(st, []).append(pt)
    coords, src, missing = {}, {}, []
    n_used = 0
    rest = []
    for r in regions:
        pts = by_region.get(r)
        if not pts:
            rest.append(r)
            continue
        lats = sorted(p[0] for p in pts)
        lons = [p[1] for p in pts]
        m = len(pts) // 2
        if len(pts) % 2:
            lat = lats[m]
        else:
            lat = (lats[m - 1] + lats[m]) / 2.0
        lon = _circular_mean_lon(lons)     # 经度：环形均值（非逐轴中位数）
        coords[r], src[r] = [round(lat, 5), round(lon, 5)], 'samples'
        n_used += len(pts)
    c2, s2, m2 = merge_coords(rest, user_coords)
    coords.update(c2)
    src.update(s2)
    missing.extend(m2)
    return coords, src, missing, n_used


# ---------------------------------------------------------------- 元数据

def load_metadata(meta_path):
    """读元数据 CSV/TSV：第一列（或 seq_id 列）为序列标识，其余为性状列。"""
    import csv
    meta = {}
    with open(meta_path, encoding='utf-8-sig', errors='replace', newline='') as f:
        sample = f.read(4096)
        f.seek(0)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
        rdr = csv.reader(f, delimiter=delim)
        rows = list(rdr)
    if not rows:
        return meta
    header = [h.strip().lstrip('#') for h in rows[0]]
    key_i = 0
    for i, h in enumerate(header):
        if h.lower() in ('seq_id', 'seqid', 'accession', 'sample', 'name'):
            key_i = i
            break
    cols = [(i, h) for i, h in enumerate(header) if i != key_i and h]
    for row in rows[1:]:
        if not row or not row[key_i].strip():
            continue
        key = row[key_i].strip()
        meta[key] = {h: row[i].strip() for i, h in cols if i < len(row)}
    return meta


def _header_traits(header):
    """FASTA 头回退解析 → `{region:…, year:…}`。

    认两种命名：
      ① 本平台约定 `>acc|区域|年份`（竖线分隔）
      ② **VirPhyKit 约定** `区域_登录号_十进制年`，如 `CCD_AB622861_2008.58197`、
         `AS_KC430335_2010.6685` —— 真实世界的病毒系统地理数据集常用这个
         （VirPhyKit 自己的 Example 全是这种），区域在**前**、年份是**十进制**。
         不认这一种的话，这类数据集进来会 0 命中（实测 Geosubsampler 的
         RSV_209CP.fasta 209 条全判未命中）。

    ⚠️ ② 是会**误判**的兜底：像 `SRR12805583_OR489165.1_2021` 这种「run_登录号_年」
    会被当成 region=SRR12805583。所以它只在 ① 不成立时启用，且要求
    恰好三段的末段是 4 位年（decimal 可带）。数据集命名不规律时建议直接给
    元数据 CSV，别靠头解析。
    """
    parts = header.split('|')
    if len(parts) >= 2:
        out = {'region': parts[1].strip()}
        if len(parts) >= 3:
            out['year'] = parts[2].strip()
        return out
    tok = (header.split() or [header])[0]
    m = re.match(r'^([A-Za-z][A-Za-z0-9]*)_([A-Za-z0-9.]+)_(\d{4}(?:\.\d+)?)$', tok)
    if m:
        return {'region': m.group(1), 'year': m.group(3)}
    return {}


def match_states(names, meta, trait):
    """把元数据映射到序列名（前缀/包含匹配）。返回 {name: state}（未匹配→''）。"""
    out = {}
    for n in names:
        st = ''
        if meta:
            if n in meta:
                st = meta[n].get(trait, '')
            else:
                for k, v in meta.items():
                    if k and (n.startswith(k) or k in n):
                        st = v.get(trait, '')
                        break
        out[n] = st.strip() or 'Unknown'
    return out


def safe_first_token(h):
    """FASTA 头 → 建树时用的净化名（与 phylo._newick_safe 同口径）。"""
    tok = h.split()[0] if h.split() else h
    return re.sub(r'[(),:;\[\]\'"]+', '_', tok) or 'seq'


def match_states_and_years(leaf_names, safe2header, meta, trait='region',
                           date_trait='year'):
    """把元数据/FASTA 头匹配到叶名。**analyze 与 ID 对齐报告共用的唯一实现**
    —— 两处若各写一套，报告说的「命中情况」就会与分析实际用的不一致。

    匹配规则（与既有行为逐字一致，别改）：
      1. 先按净化名↔原头映射找回原头；找不到就按「前缀/包含」找第一个
      2. 有元数据表：按 `dates` 的**插入序**取第一个满足
         `头.startswith(键) or 键 in 头 or 叶名.startswith(键)` 的行
      3. 元数据没给该列 → 回退解析 FASTA 头（`>acc|区划|年份` 约定）
      4. 还取不到 → 区划记 'Unknown'（unresolved+1）；年份留空

    2026-09-21 修：trait 留空时**自动识别**区划列——取元数据值字典里第一个
    命中 ('region','location','country','geo_location') 的键；都没有再退回
    'region'（FASTA 头约定）。此前 trait='' 会传到 v.get('') → 所有节点
    区划判空 → 全 Unknown（系统地理卡示例实测踩中）。
    返回 (states{叶名:区划}, leaf_year{叶名:int年}, unresolved, rows[明细])。
    """
    if not trait and meta:
        first_val = next(iter(meta.values()), {})
        for cand in ('region', 'location', 'country', 'geo_location'):
            if cand in first_val:
                trait = cand
                break
        else:
            trait = 'region'
    states, leaf_year, rows = {}, {}, []
    unresolved = 0
    items = list(meta.items()) if meta else []
    for n in leaf_names:
        h = safe2header.get(n)
        if h is None:
            h = next((hh for k, hh in safe2header.items()
                      if n.startswith(k) or k in n), n)
        st, yr_raw, hit_key = '', '', ''
        for k, v in items:
            if k and (h.startswith(k) or k in h or n.startswith(k)):
                st = v.get(trait, '')
                yr_raw = v.get(date_trait, '')
                hit_key = k
                break
        src = 'meta' if hit_key else ''
        traits = _header_traits(h) or {}
        if not st:
            st = traits.get(trait, '')
            if st:
                src, hit_key = 'header', h
        if not yr_raw:
            yr_raw = traits.get(date_trait, '')
        yr = None
        if yr_raw:
            _ym = re.match(r'(\d{4})', str(yr_raw))
            if _ym:
                yr = int(_ym.group(1))
                leaf_year[n] = yr
        if not st:
            st = 'Unknown'
            unresolved += 1
            src = src or 'none'
        states[n] = st
        rows.append({'seq': n, 'header': h, 'matched_key': hit_key,
                     'region': st, 'year': yr, 'source': src or 'meta',
                     'status': 'unmatched' if st == 'Unknown' else 'matched'})
    return states, leaf_year, unresolved, rows


def id_alignment_report(aln_path, meta_path=None, trait='region',
                        date_trait='year'):
    """SeqIDRename + SeqGrouper 的核心：报告每条序列对到了什么、分成几组。

    解决的实际痛点：区划列没匹配上时会**静默回退**去读 FASTA 头的第二段
    （`acc|X|年份`）—— 若头是 `acc|<整条描述>`，每条序列都变成"独立区划"，
    迁移矩阵与置换检验随之失去意义，而界面上看不出来。这里把命中/未命中、
    分组计数、以及**不确定**的情形全部摆出来。

    ⚠️ 键统一用 `safe_first_token(头)`（净化名）—— 与 analyze 从树里拿到的
    叶名同口径，否则「报告说命中」而分析实际用的是另一套键。

    返回 dict：n_seqs / n_matched / n_unmatched / groups{区划:条数} /
    n_meta_keys / dup_ids / rows[明细] / warnings[]。
    """
    from Virus_Platform_Core.utils import load_alignment
    aln = load_alignment(aln_path)
    headers = _fasta_headers(aln_path)
    safe2header, dup = {}, {}
    for h in headers:
        k = safe_first_token(h)
        if k in safe2header:
            dup.setdefault(k, [safe2header[k]]).append(h)
        safe2header.setdefault(k, h)
    meta = load_metadata(meta_path) if (meta_path and os.path.isfile(meta_path)) else {}
    states, leaf_year, unresolved, rows = match_states_and_years(
        list(safe2header), safe2header, meta, trait=trait, date_trait=date_trait)
    for r in rows:                      # 明细里补上原始头，便于人工核对
        r['raw_header'] = safe2header.get(r['seq'], '')
    groups = Counter(states.values())
    warnings = []
    known = [v for v in states.values() if v and v != 'Unknown']
    if len(known) >= 5 and len(set(known)) == len(known):
        warnings.append(f'{len(known)} 条序列有 {len(set(known))} 个各不相同的'
                        f'{trait} —— 多为该列没匹配上、回退成了整条 FASTA 描述；'
                        '请检查元数据键或把 FASTA 头写成 `名称|区划|年份`。')
    if unresolved:
        warnings.append(f'{unresolved} 条序列没有 {trait}（分析时按缺失数据处理）。')
    if dup:
        warnings.append('净化后同名的序列（建树时会撞名）：'
                        + '；'.join(f'{k}（{len(v)} 条）'
                                   for k, v in list(dup.items())[:5]))
    return {
        'n_seqs': len(aln),
        'n_matched': sum(1 for r in rows if r['status'] == 'matched'),
        'n_unmatched': unresolved,
        'groups': dict(groups.most_common()),
        'n_meta_keys': len(meta),
        'dup_ids': {k: v for k, v in list(dup.items())[:20]},
        'rows': rows,
        'warnings': warnings,
        'trait': trait,
        'date_trait': date_trait,
    }


def write_normalized_fasta(aln_path, out_path, meta_path=None, trait='region',
                           date_trait='year'):
    """导出**ID 规范化**的 FASTA：头统一成 `净化名|区划|年份`。

    这样后续任何一步（建树 / 系统地理 / 定年）都能直接从头的 `|` 段解析，
    不必再依赖元数据表的键匹配。返回 (写出条数, 报告 dict)。
    """
    from Virus_Platform_Core.utils import load_alignment
    rep = id_alignment_report(aln_path, meta_path, trait, date_trait)
    aln = load_alignment(aln_path)
    by_key = {r['seq']: r for r in rep['rows']}
    n = 0
    d = os.path.dirname(os.path.abspath(out_path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8', newline='\n') as f:
        for raw_name, seq in aln.items():
            key = safe_first_token(raw_name)
            r = by_key.get(key) or {}
            reg = r.get('region') or 'Unknown'
            yr = r.get('year')
            # 基名取第一个 `|` 之前那一段：原头若本来就是 `acc|区划|年份`，
            # 直接用整串会拼成 `acc|区划|年份|区划|年份`（重复）。
            # 这样写也保证本函数**幂等**（重复跑不会越拼越长）。
            base = key.split('|')[0] or 'seq'
            f.write(f'>{base}|{reg}{"|" + str(yr) if yr else ""}\n')
            f.write(seq + '\n')
            n += 1
    return n, rep


def apply_rename_table(aln_path, table_path, out_path, delim=None):
    """按「旧名 → 新名」对照表批量改名（VirPhyKit **SeqIDRenamer** 的核心能力）。

    VirPhyKit 的 `Example/SeqIDRenamer/SeqIDRenamer.txt` 是**无表头的两列 TSV**
    （`PL086098<TAB>Seq01_trait01_date01`）。这里对分隔符不敏感：显式 `delim`
    优先；否则按「TAB > 逗号 > 空白」自动判别。

    匹配规则：与平台其他地方一致 —— **先精确、再前缀/包含**（FASTA 头常带
    `|` 描述段，而对照表里只写登录号）。表里没命中的序列**原样保留**并计入
    `unmatched`（不静默丢弃）。

    返回 (n_out, n_renamed, unmatched[列表], n_extra_table_rows)。
    """
    from Virus_Platform_Core.utils import load_alignment
    aln = load_alignment(aln_path)

    pairs = []
    with open(table_path, encoding='utf-8-sig', errors='replace') as f:
        for ln in f.read().splitlines():
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            if delim:
                parts = ln.split(delim)
            elif '\t' in ln:
                parts = ln.split('\t')
            elif ',' in ln:
                parts = ln.split(',')
            else:
                parts = ln.split()
            parts = [p.strip() for p in parts]
            if len(parts) >= 2 and parts[0] and parts[1]:
                pairs.append((parts[0], parts[1]))

    exact = {a: b for a, b in pairs}
    renamed, unmatched = 0, []
    d = os.path.dirname(os.path.abspath(out_path))
    if d:
        os.makedirs(d, exist_ok=True)
    used = set()
    with open(out_path, 'w', encoding='utf-8', newline='\n') as f:
        for raw, seq in aln.items():
            base = safe_first_token(raw)
            new = exact.get(base) or exact.get(raw)
            if new:
                used.add(base)
            else:
                # 前缀/包含回退（表里可能只写了登录号，而头里是 `acc|描述`）
                hit = next(((a, b) for a, b in pairs
                            if base.startswith(a) or a in base), None)
                if hit:
                    new = hit[1]
                    used.add(hit[0])
            if new:
                renamed += 1
            else:
                new = base
                unmatched.append(base)
            f.write(f'>{new}\n' + seq + '\n')
    return (len(aln), renamed, unmatched,
            max(0, len(pairs) - len(used)))


# ---------------------------------------------------------------- Newick


def parse_newick(text):
    """极简 Newick 解析 → 根节点 dict {name, length, children[]}。

    支持括号树、内部节点命名/长度、叶名与长度；不处理 NHX 注释。

    ⚠️ 解析不完就**抛错**，不许「截断成子树」返回：括号写错（少/多一个 `)`）
    或顶层出现多个逗号分隔子树时，旧实现会把剩余部分静默丢掉 —— 于是下游
    （Fitch / MOTP / 时间树）拿到的是一棵**少了若干样品的树**，结果照常算出来却
    只覆盖子集，谁都不会发现。这里改成 ValueError，附带位置与片段，让调用方
    在入口就炸掉。
    """
    s = text.strip().rstrip(';')
    pos = [0]
    _head = s[:60].replace('\n', ' ')

    def _name_len():
        buf = ''
        while pos[0] < len(s) and s[pos[0]] not in ',():':
            buf += s[pos[0]]
            pos[0] += 1
        name = buf.strip()
        length = 0.0
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            buf = ''
            while pos[0] < len(s) and s[pos[0]] not in ',()':
                buf += s[pos[0]]
                pos[0] += 1
            try:
                length = float(buf)
            except ValueError:
                length = 0.0
        return name, length

    def _node():
        node = {'name': '', 'length': 0.0, 'children': []}
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1
            while True:
                node['children'].append(_node())
                if pos[0] < len(s) and s[pos[0]] == ',':
                    pos[0] += 1
                    continue
                if pos[0] < len(s) and s[pos[0]] == ')':
                    pos[0] += 1
                    break
                # ⚠️ 读到结尾（或意外字符）都没等到 ',' / ')' —— 括号不配对。
                # 没有这道护栏时 while 会无限递归 append(_node())，以
                # MemoryError 收场（2026-09-16 实测截断树 '((a:1,b:1)' 触发；
                # 文档一直写"解析不完就抛错"，但此前只覆盖"顶层剩余字符"，
                # 漏了"括号没闭合"这条路径）。
                raise ValueError(
                    'Newick 括号不配对：在位置 %d 读到 %s，预期 "," 或 ")"。'
                    '树开头：%s'
                    % (pos[0],
                       '结尾' if pos[0] >= len(s) else repr(s[pos[0]]),
                       _head))
            node['name'], node['length'] = _name_len()
        else:
            node['name'], node['length'] = _name_len()
        return node

    root = _node()
    rest = s[pos[0]:].strip()
    if rest:
        raise ValueError(
            'Newick 解析未走完：在位置 %d 处还剩 %d 个字符没解析（片段：%r）——'
            '多半是括号不配对或顶层有多个逗号分隔的子树；继续下去会把剩余分支'
            '静默丢掉，故直接报错。树开头：%s'
            % (pos[0], len(rest), rest[:40], _head))
    return root


def write_annotated(node, states):
    """导出带区划标注的 Newick：叶名后缀 `|状态`，内部节点写内部名。"""
    counter = [0]

    def _rec(n):
        if not n['children']:
            st = states.get(n['name'], 'Unknown')
            name = n['name'] if n['name'] else 'leaf'
            return f"{name}|{st}:{n['length']}"
        kids = ','.join(_rec(c) for c in n['children'])
        counter[0] += 1
        name = n['name'] or f'Node{counter[0]}'
        return f"({kids}){name}:{n['length']}"

    return _rec(node) + ';'


# ---------------------------------------------------------------- Fitch 简约

def _iter_post(n):
    """后序遍历（先子后父）。"""
    for c in n['children']:
        yield from _iter_post(c)
    yield n


def _iter_pre(n):
    """先序遍历（先父后子）。"""
    yield n
    for c in n['children']:
        yield from _iter_pre(c)


def _fitch(root, states):
    """Fitch 最大简约：离散性状（区划）重构。**纯函数，不写回节点**。

    经典两遍（Felsenstein 无序性状）：
      up-pass（后序）—— 每个节点求**候选状态集**：
          叶 = {观测态}；内部 = 各子集两两求交，交为空则取并集
      down-pass（先序）—— 每个节点**定一个具体状态**：
          根 = 候选集的规范首元；其余 = 父状态若落在自身候选集内则沿用，否则取规范首元

    返回 (state_by_node {id(node): state}, transitions [(pname, cname, from, to)])。

    ⚠️ 修的是这里（2026-09-15 AUDIT 记的缺陷）：原实现在 up-pass 里
    **也给内部节点写了 `_state`**，于是 down-pass 的 `if '_state' not in c`
    恒假、**自顶向下那一遍从未执行**，内部节点状态退化成「候选集字典序最小者」。
    后果实测：系统性**高估迁移数**，并把内部状态推向字母序靠前的区划
    （Africa 占 59–64%，而叶里真实只 ~25%）。现在 up-pass 只算集合、不写值。

    ⚠️ 另一处同源缺陷一并修：原来把「元数据缺失」当成字面状态 `'Unknown'`，
    于是每个缺元数据的叶都凭空造出一次迁移（父态 → Unknown）。Fitch 对
    **缺失数据**的正确处理是「任意状态皆可」——这里给这类叶赋全体观测态集合，
    因此它们不再强行制造迁移（缺多少条仍由 analyze() 的 `unresolved` 如实报出）。

    注：Fitch 的**最少变化数**唯一，但**具体悬挂方式不唯一**（多个赋值同样最优）
    → 迁移矩阵依赖 tie-break 选择。这正是要配 bootstrap 置信区间的原因
    （见 phylogeo.rssp_bootstrap）。
    """
    universe = {v for v in states.values() if v and v != 'Unknown'}
    missing = universe or {'Unknown'}          # 缺数据的叶：任意状态皆可
    leaf_set = {}

    def _leaf_states(name):
        st = states.get(name, '')
        if not st or st == 'Unknown':
            return missing
        return {st}

    # ---- up-pass：候选状态集（只存集合，不写节点）----
    sets = {}
    for n in _iter_post(root):
        if not n['children']:
            s = _leaf_states(n['name'])
            leaf_set[id(n)] = s
            sets[id(n)] = set(s)
            continue
        # 多元（含二叉、多叉）Fitch —— 对每个候选状态 x 统计"候选集包含 x 的
        # 子节点个数" c(x)，祖先候选集 = argmax_x c(x)。
        #
        # 推导：节点取状态 x 时，与子节点 c 的连接代价 = 0（x ∈ S_c）或 1（x ∉ S_c），
        # 故该节点总代价 = m − c(x)（m = 子节点数）；取最小即 max c(x)。
        # 并列最大值的状态**全部保留**在候选集里（悬挂方式不唯一，交给 down-pass
        # 的 tie-break 决定，这是有意的）。
        #
        # ⚠️ 2026-10-03 修：原实现是"两两求交，交空取并集"—— 二叉时与上式等价
        # （交集非空 → 每个 x 的 c=2、argmax=交集；交集空 → 每个 x 至多属一个子集、
        # argmax=并集），但**多叉节点**会漏掉"某状态被多数子集共享、交集却已空"
        # 的情况 → 祖先候选集偏大 → 系统性高估最少迁移步数。审计反例
        # `((B,A),(A,B,B),(B,A))`：旧 4 步、真值 3 步（修后 3 步，见回归用例）。
        cand = {}
        for c in n['children']:
            for x in sets[id(c)]:
                cand[x] = cand.get(x, 0) + 1
        if not cand:
            sets[id(n)] = {'Unknown'}
            continue
        best = max(cand.values())
        sets[id(n)] = {x for x, v in cand.items() if v == best} or {'Unknown'}

    def _canon(s):
        """候选集的规范取值：优先非 Unknown（确定性 tie-break，便于复现）。"""
        known = sorted(x for x in s if x != 'Unknown')
        return known[0] if known else 'Unknown'

    # ---- down-pass：定值 + 数转移 ----
    #
    # ⚠️ 状态只放在**局部字典**里，绝不写回节点。原因：置换检验会对同一棵树
    # 反复调用本函数（几百次），若写回 `n['_state']`，等置换跑完节点上留的是
    # **最后一次置换**的状态，随后 write_annotated() 就会把标注树写错。
    # 写回只由 fitch_mugration() 负责，且只在「用观测标签算完」之后做一次。
    state_by_node = {id(root): _canon(sets[id(root)] or {'Unknown'})}
    transitions = []
    child_ids = []                             # 与 transitions 平行：子节点 id
    for n in _iter_pre(root):
        pstate = state_by_node[id(n)]
        for c in n['children']:
            cs = sets.get(id(c)) or {'Unknown'}
            # 父状态在子候选集内 → 沿用（不产生转移）；否则取候选集规范首元
            cstate = pstate if pstate in cs else _canon(cs)
            state_by_node[id(c)] = cstate
            if cstate != pstate:
                transitions.append((n.get('name') or '(root)',
                                    c.get('name') or '(node)',
                                    pstate, cstate))
                child_ids.append(id(c))
    return state_by_node, transitions, child_ids


def fitch_mugration(root, states):
    """Fitch 重构并**把状态写回节点**（`_state`，供标注树使用）。

    返回 (node_state {id(node): state}, transitions [(parent, child, from, to)],
          n_changes)。
    """
    state_by_node, transitions, _child_ids = _fitch(root, states)
    for n in _iter_pre(root):
        n['_state'] = state_by_node[id(n)]
    return state_by_node, transitions, len(transitions)


def count_transitions(root, states):
    """只数迁移次数，**不改动节点** —— 供置换检验/bootstrap 反复调用。

    必须与 fitch_mugration 走同一套 `_fitch`，否则置换零分布与实测量不同口径。
    """
    return len(_fitch(root, states)[1])


def region_permutation_test(root, states, n_perm=1000, seed=0):
    """区域随机化检验（简约法 / BaTS 式）：区划与系统树是否显著关联。

    把 tip 的区划标签在叶之间随机重排 n_perm 次（**保持各区划条数不变**，
    即只打乱「谁属于哪个区划」），每次重跑同一套 Fitch 迁移计数 → 零分布。

    统计量＝Fitch 最少迁移数。地理结构与树相关 ⇒ 迁移数**少于**随机 ⇒
    关心的是**左尾**。两个尾都给，别只看一个：

      p_low  = P(零分布 ≤ 实测)   小 ⇒ 实测显著**少**于随机 ⇒ **地理结构显著**
      p_high = P(零分布 ≥ 实测)   小 ⇒ 实测显著**多**于随机 ⇒ 区划过度分散

    ⚠️ **与 VirPhyKit 的 RRT 不是同一个东西，别混称**：VirPhyKit 的
    `RRT/function_rrt.py` 用的是**贝叶斯后验概率**口径 —— 取「最可能区划的后验
    概率最大值」作统计量，零分布来自把区划随机化后**重跑贝叶斯分析**得到的 N 组
    MCC 树（`max_prob > max(随机组的 max_prob)` 才判显著）。那需要 BEAST 类 MCMC。
    这里是**无需模型的简约法替代**（BaTS 的思路），纯 Python、秒级，能回答同一个
    科学问题（区划是否有信号），但**不能声称是 VirPhyKit 的实现**。
    两者都是"与时间信号无关"的另一件事，别互相替代。

    返回 dict：observed / n_perm / null_mean / null_sd / null_min / null_max /
    p_low / p_high / seed / n_regions / interpretation。
    """
    import random as _random

    labels = sorted({v for v in states.values() if v and v != 'Unknown'})
    # 只打乱「有区划」的叶；缺数据的叶保持缺数据（Fitch 里=任意状态，不影响计数）
    keys = [k for k, v in states.items() if v and v != 'Unknown']
    observed = count_transitions(root, states)
    if not keys or len(labels) < 2:
        return {'observed': observed, 'n_perm': 0, 'p_low': None,
                'p_high': None, 'null_mean': None, 'null_sd': None,
                'null_min': None, 'null_max': None, 'seed': seed,
                'n_regions': len(labels), 'null_counts': [],
                'interpretation': '可用区划 < 2 或无可打乱样本，未做置换'}

    pool = [states[k] for k in keys]           # 保持各区划条数不变的打乱池
    rng = _random.Random(seed)
    counts = []
    for _ in range(int(n_perm)):
        rng.shuffle(pool)                      # 对多重集洗牌 → 均匀于全部排列
        perm = dict(states)
        for k, v in zip(keys, pool):
            perm[k] = v
        counts.append(count_transitions(root, perm))

    n_le = sum(1 for x in counts if x <= observed)
    n_ge = sum(1 for x in counts if x >= observed)
    n = len(counts)
    mean = sum(counts) / n
    var = (sum((x - mean) ** 2 for x in counts) / n) if n > 1 else 0.0
    p_low = round((n_le + 1) / (n + 1), 5)
    p_high = round((n_ge + 1) / (n + 1), 5)
    if p_low < 0.05 and p_low <= p_high:
        interp = '地理结构显著（迁移数显著少于随机）'
    elif p_high < 0.05:
        interp = '区划显著过度分散（迁移数显著多于随机）'
    else:
        interp = '未检出显著地理结构'
    return {
        'observed': observed,
        'n_perm': n,
        'null_mean': round(mean, 3),
        'null_sd': round(var ** 0.5, 3),
        'null_min': min(counts),
        'null_max': max(counts),
        'p_low': p_low,
        'p_high': p_high,
        'seed': seed,
        'n_regions': len(labels),
        'interpretation': interp,
        # 零分布全量（写 rrt_permutation.tsv 用；summary.json 前会被 pop 掉）
        'null_counts': counts,
    }


# ---------------------------------------------------------------- 引入 vs 本地传播

def _contiguous_group(child, node_states, states, region):
    """跨区边（祖先区 → region）的后代里，与 region **连续同区**的叶数。

    从跨区边的子节点出发，只沿「重构状态 == region」的节点往下走；一遇到
    别的区划就停（那是另一次跨区边的事）。组大小只数**观测到 region 的叶**：
    缺元数据（'Unknown'）的叶不计入，也不额外造组 —— 与 `_fitch` 的
    「缺失=任意状态、不凭空造迁移」同一口径（见 tests/_check_phylogeo_fitch.py）。
    """
    tips = 0
    stack = [child]
    while stack:
        n = stack.pop()
        if not n['children']:
            if states.get(n['name']) == region:
                tips += 1
            continue
        for c in n['children']:
            if node_states.get(id(c)) == region:
                stack.append(c)
    return tips


def import_export_balance(root, states, meta=None, node_states=None):
    """引入 vs 本地传播（import / export）—— 1/group_size 口径的纯 Python 实现。

    算法出处（先读后重实现，**不抄代码**）：SARS-CoV-2_phylogeo 的
    `docs/Analyses/R_functions/Relative_load_importation_local_transmission.R`
    第 15–30 行（对每个后代组取 `group_size = Result$Lineage_sizes[group_no]`，
    逐叶累加 `importations += 1/group_size`、`local_transmissions += 1 − 1/group_size`），
    用法见同仓 `docs/Analyses/Omicron/01_Mapping_Omicron.Rmd` 207–216 行
    （LineageHomology：每个焦点区支系都始于一次引入，其后的分叉都是本地传播）。
    该仓 GPL-3.0 → 只学方法、不复制代码。

    本实现的语义（对 R 脚本逐叶摊派的**求和等价形式**）：
      · 每条**跨区边**（祖先区 A → 后代区 B，A≠B）＝ 1 次引入：记到 B 的
        importations、同时记 1 次 exportations 到 A（方向性：祖先区→后代区）；
      · 该边后代中与 B 连续同区的叶组大小 k → 组内贡献 k−1 次本地传播；
        逐叶看即每叶 1/k 次引入、1−1/k 次本地传播（与 R 脚本一致）。
      ⇒ per_region[r] = importations / local_transmissions / exportations /
        net(=importations−exportations)；全局 Σ net = 0（一进一出）。
      ⇒ 组大小 k=0 的跨区边（后代没有观测到 B 的叶）不计权重，只计入
        n_empty_groups —— 与 R 脚本「组里没有叶就不会进循环」一致。

    缺失状态（元数据缺 → 'Unknown'）沿用修复版 `_fitch` 语义：既不凭空造迁移，
    也不把缺数据的叶算进任何组（`_contiguous_group` 只认观测态）。

    参数：
      meta       可选 {'focal': 区划名}（或直接传区划名字符串）→ 额外给出
                 「焦点区划」的引入/本地传播拆分（R 脚本默认就是这个单焦点口径）
      node_states 可选，`_fitch` 已算好的 {id(节点): 状态}（analyze 里复用，
                 避免第二遍重构）；不给则内部调用 `_fitch`（纯函数，不改树）

    返回 JSON 友好 dict：per_region / top_corridors / flow_partition /
    n_crossings / n_empty_groups / focal / method ...
    """
    if node_states is None:
        node_states, _tr, _cids = _fitch(root, states)
    if isinstance(meta, str):
        meta = {'focal': meta}
    focal_region = (meta or {}).get('focal') if isinstance(meta, dict) else None

    per = {}

    def _row(r):
        return per.setdefault(r, {'n_tips': 0, 'n_import_groups': 0,
                                  'importations': 0.0, 'local_transmissions': 0.0,
                                  'exportations': 0.0, 'net': 0.0})

    n_missing = 0
    for v in states.values():
        if v and v != 'Unknown':
            _row(v)['n_tips'] += 1
        else:
            n_missing += 1

    corridors = {}
    groups = []
    n_cross = 0
    n_empty = 0
    for n in _iter_pre(root):
        pstate = node_states.get(id(n))
        for c in n['children']:
            cstate = node_states.get(id(c))
            if cstate == pstate:
                continue
            n_cross += 1
            k = _contiguous_group(c, node_states, states, cstate)
            if k < 1:
                # 该跨区边下没有任何观测到 cstate 的叶 → 不进权重（R 脚本同口径）
                n_empty += 1
                continue
            _row(cstate)['importations'] += 1.0
            _row(cstate)['local_transmissions'] += k - 1.0
            _row(cstate)['n_import_groups'] += 1
            _row(pstate)['exportations'] += 1.0
            co = corridors.setdefault((pstate, cstate),
                                      {'weight': 0.0, 'n_edges': 0, 'n_tips': 0})
            co['weight'] += 1.0
            co['n_edges'] += 1
            co['n_tips'] += k
            groups.append({'from': pstate, 'to': cstate, 'group_size': k,
                           'import_weight': round(1.0 / k, 6),
                           'local_weight': round(1.0 - 1.0 / k, 6),
                           'child': c.get('name') or '(node)'})

    # 四舍五入 + net；排序：引入多的区划排前面（表格/前端直接按序渲染）
    per_region = OrderedDict()
    for r in sorted(per, key=lambda x: (-per[x]['importations'],
                                        -per[x]['exportations'], x)):
        v = per[r]
        v = {'n_tips': v['n_tips'], 'n_import_groups': v['n_import_groups'],
             'importations': round(v['importations'], 6),
             'local_transmissions': round(v['local_transmissions'], 6),
             'exportations': round(v['exportations'], 6),
             'net': round(v['importations'] - v['exportations'], 6)}
        per_region[r] = v

    tot_imp = sum(v['importations'] for v in per_region.values())
    tot_local = sum(v['local_transmissions'] for v in per_region.values())
    denom = tot_imp + tot_local
    top_corridors = [{'from': a, 'to': b, 'weight': round(v['weight'], 6),
                      'n_edges': v['n_edges'], 'n_tips': v['n_tips']}
                     for (a, b), v in sorted(
                         corridors.items(),
                         key=lambda kv: (-kv[1]['weight'], -kv[1]['n_tips'],
                                         kv[0][0], kv[0][1]))]
    focal = None
    if focal_region:
        fv = per_region.get(focal_region)
        if fv:
            fd = fv['importations'] + fv['local_transmissions']
            focal = {'region': focal_region,
                     'importations': fv['importations'],
                     'local_transmissions': fv['local_transmissions'],
                     'import_frac': (round(fv['importations'] / fd, 4)
                                     if fd else None),
                     'local_frac': (round(fv['local_transmissions'] / fd, 4)
                                    if fd else None)}
    return {
        'per_region': per_region,
        'top_corridors': top_corridors,
        'flow_partition': {
            # 加权口径：组内分支 = 本地传播（Σ k−1），跨区边 = 引入（Σ 1）
            'within_region': round(tot_local, 6),
            'between_region': round(tot_imp, 6),
            'within_frac': round(tot_local / denom, 4) if denom else None,
            'between_frac': round(tot_imp / denom, 4) if denom else None,
            'basis': '按 1/group_size 加权：组内分支记本地传播，跨区边记引入；'
                     '分母＝两者之和（R 脚本的 sum_cases 口径）',
        },
        'n_crossings': n_cross,
        'n_empty_groups': n_empty,
        'n_weighted_crossings': n_cross - n_empty,
        'n_missing_tips': n_missing,
        'n_tips': sum(v['n_tips'] for v in per_region.values()),
        'n_regions': len(per_region),
        'groups': groups,
        'focal': focal,
        'method': ('LineageHomology / 1/group_size 口径：跨区边 = 1 次引入'
                   '（记到后代区）+ 1 次输出（记自祖先区）；该边下连续同区叶组 k 个'
                   '→ 本地传播 k−1 次（逐叶即 1/k 引入、1−1/k 本地）；'
                   '缺元数据的叶不计入组、也不造迁移。'
                   '来源：SARS-CoV-2_phylogeo（GPL-3.0，只学方法）'
                   'Relative_load_importation_local_transmission.R 15–30 行。'
                   '口径注意：节点状态来自 Fitch 简约重构（非 MCC 注释 MAP 态），'
                   '跨区计数与 MOT/注释后验口径不可直接对比'
                   '（同一棵 RSV MCC 树实测 54 vs 130 条边；'
                   '见 docs/方法对比_VirPhyKit_Example_20260916.md）。'),
    }


# ---------------------------------------------------------------- 迁移事件分类

# 迁移事件分类：主标签 = 距离分带（三分类，对齐 phymapr 的三档距离），
# 顺序即表格/排序优先级；Unresolved 居末（没坐标，不是距离档）。
_CLASS_ORDER = ('Direct', 'Indirect', 'Distant', 'Unresolved')


def _subtree_tips(root):
    """节点 id → (后代叶名列表, 节点对象)；另给 子 id → 父 id。

    分类要按「边」把两侧的实际样本摆出来（最近样本对距离、两侧采样密度），
    所以需要**后代叶集合**与节点本体（拿枝长）。内部节点名常常是空的
    （NJ/FastTree 输出），不能靠名字定位——一律用 id()。
    """
    tips, node_of, parent_of = {}, {}, {}

    def walk(n):
        names = []
        for c in n['children']:
            parent_of[id(c)] = id(n)
            names.extend(walk(c))
        tips[id(n)] = names or [n['name'] or '']
        node_of[id(n)] = n
        return tips[id(n)]

    walk(root)
    return tips, node_of, parent_of


def _median(vals):
    """中位数（空 → None）。"""
    if not vals:
        return None
    s = sorted(vals)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2.0


def classify_transitions(root, states, node_states, coords, sample_points,
                         leaf_year, aln_len):
    """把每条跨区迁移边按**地理距离**分类：Direct / Indirect / Distant。

    分类口径（自研；**不沿用 phymapr 的固定公里阈值**——固定阈值换数据集就失真）：
      · 两端区划锚点（样本坐标中位数 > 用户表 > 内置质心表）的 haversine
        距离分**分位带**：≤q33 → Direct、≥q67 → Distant、中间 → Indirect；
        阈值随本次运行自适应；跨区边 <6 条或分位退化（q67≤q33）时退回
        500/3000 km，并在 thresholds.src 如实标注；
      · 任一端没有坐标 → Unresolved（如实报，不硬塞档）。
    另给两个正交字段（不进主分类，避免把距离带吞掉）：
      intro      该边是否为「引入边」——与 import/export 同一判定
                 （`_contiguous_group` 在子侧数到 ≥1 条观测叶）。注意：带观测叶
                 的跨区边几乎都满足它，所以**不当主标签**，否则三分类形同虚设；
      group_size k：该边下连续同区观测叶组大小（k=1 孤立引入；k>1 引入后又
                 本地传播 k−1 次；k=0 后代没有观测叶 → 记 note=no-observed-tips）。
    四通道证据（全部进表，不只做主分类）：
      geo_km   两端锚点距离（主分类轴）
      step_km  跨边**实测样本对**的最近距离（锚点是中位数落点，这个更贴边）
      edge_snp 边替换数 ≈ 枝长 × 比对长度（低分化远程 → 疑似人为/运输带毒）
      dt_years 两侧样本年份中位数差（同期远程 → 同上）
      采样密度 两侧局部叶数 n_p/n_c（任一侧 <3 记 sparse-sampling，分类不稳要提示）
    返回 dict：rows[] / counts{} / thresholds{q33,q67,src} / flags{} /
    n_edges / n_intro / n_no_tips / n_unknown_state / n_state_changes / snp_median。
    """
    tips_of, _node_of, parent_of = _subtree_tips(root)
    rows = []
    n_unknown = 0
    n_changes = 0
    for n in _iter_pre(root):
        pstate = node_states.get(id(n))
        for c in n['children']:
            cstate = node_states.get(id(c))
            if not pstate or not cstate or pstate == cstate:
                continue
            n_changes += 1
            if 'Unknown' in (pstate, cstate):
                n_unknown += 1
                continue
            cid = id(c)
            names_c = [x for x in tips_of.get(cid, []) if x]
            pid = parent_of.get(cid)
            set_c = set(names_c)
            names_p = [x for x in tips_of.get(pid, []) if x and x not in set_c]
            geo = None
            if pstate in coords and cstate in coords:
                geo = round(haversine_km(coords[pstate], coords[cstate]), 1)
            pc = [sample_points[x] for x in names_c if x in sample_points]
            pp = [sample_points[x] for x in names_p if x in sample_points]
            step = None
            if pc and pp:
                step = round(min(haversine_km(a, b) for a in pp for b in pc), 1)
            yc = [leaf_year[x] for x in names_c if x in leaf_year]
            yp = [leaf_year[x] for x in names_p if x in leaf_year]
            dt = None
            if yc and yp:
                dt = round(_median(yc) - _median(yp), 3)
            rows.append({
                'parent': n['name'] or '(node)',
                'child': c['name'] or '(node)',
                'from': pstate, 'to': cstate,
                'geo_km': geo, 'step_km': step,
                'edge_snp': int(round((c['length'] or 0.0) * aln_len)),
                'dt_years': dt,
                'n_p': len(names_p), 'n_c': len(names_c),
                # 引入判定与 import_export 同一口径（同一函数、同一 node_states）
                'group_size': _contiguous_group(c, node_states, states,
                                                cstate),
            })
    for r in rows:
        r['intro'] = r['group_size'] >= 1

    geos = [r['geo_km'] for r in rows if r['geo_km'] is not None]
    if len(geos) >= 6:
        q33 = _quantile(geos, 1.0 / 3.0)
        q67 = _quantile(geos, 2.0 / 3.0)
        tsrc = 'run-quantiles'
        if q67 <= q33:                       # 距离全相同/退化 → 固定阈值兜底
            q33, q67, tsrc = 500.0, 3000.0, 'default(degenerate)'
    else:
        q33, q67, tsrc = 500.0, 3000.0, 'default(n<6)'
    snp_med = _median([r['edge_snp'] for r in rows])
    for r in rows:
        if r['geo_km'] is None:
            label = 'Unresolved'
        elif r['geo_km'] <= q33:
            label = 'Direct'
        elif r['geo_km'] >= q67:
            label = 'Distant'
        else:
            label = 'Indirect'
        notes = []
        if not r['group_size']:
            # 该边后代没有观测到该区的叶：R 脚本同口径——不记权重，也不能当实证迁移
            notes.append('no-observed-tips')
        if r['geo_km'] is not None and r['geo_km'] >= q67:
            if snp_med is not None and r['edge_snp'] <= snp_med:
                notes.append('low-divergence')
            if r['dt_years'] is not None and abs(r['dt_years']) <= 1:
                notes.append('same-period')
        if min(r['n_p'], r['n_c']) < 3:
            notes.append('sparse-sampling')
        r['label'] = label
        r['note'] = ';'.join(notes)
    order = {k: i for i, k in enumerate(_CLASS_ORDER)}
    rows.sort(key=lambda r: (order.get(r['label'], 9),
                             -(r['geo_km'] if r['geo_km'] is not None else -1)))
    return {
        'rows': rows,
        'counts': dict(Counter(r['label'] for r in rows)),
        'thresholds': {'q33': q33, 'q67': q67, 'src': tsrc},
        'flags': dict(Counter(t for r in rows
                              for t in r['note'].split(';') if t)),
        'n_edges': len(rows),
        'n_intro': sum(1 for r in rows if r['intro']),
        'n_no_tips': sum(1 for r in rows if not r['group_size']),
        'n_unknown_state': n_unknown,
        'n_state_changes': n_changes,
        'snp_median': snp_med,
        'n_with_step': sum(1 for r in rows if r['step_km'] is not None),
        'method': ('分类口径（自研）：主标签 = 两端区划锚点的 haversine 距离'
                   '**运行内分位带**（≤q33 Direct / ≥q67 Distant / 中间 Indirect）。'
                   '阈值随数据自适应，不沿用 phymapr 的固定公里阈值；'
                   '跨区边 <6 条或分位退化时退回 500/3000 km（thresholds.src 标注）。'
                   '低分化远程（edge_snp ≤ 中位数）与同期远程（|dt| ≤1 年）'
                   '只作证据标记，不改分类。另给 intro/group_size（是否为引入边'
                   '及其连续同区叶组大小；与 import/export 同一判定，因几乎每条'
                   '带叶跨区边都满足，故不作主标签）。'),
    }


def corridor_classes(rows):
    """按走廊（from→to）聚合分类计数 → [{from,to,label,counts,n}]。

    地图上一条弧线往往对应多条边：label = 计数最多者，同票按 `_CLASS_ORDER`
    （近距离 → 远距离）定夺；intro/n_intro 单独带上（引出与距离是两条轴）。
    """
    acc = {}
    for r in rows:
        e = acc.setdefault((r['from'], r['to']),
                           {'counts': Counter(), 'n_intro': 0,
                            'max_group': 0})
        e['counts'][r['label']] += 1
        if r.get('intro'):
            e['n_intro'] += 1
        e['max_group'] = max(e['max_group'], r.get('group_size') or 0)
    order = {k: i for i, k in enumerate(_CLASS_ORDER)}
    out = []
    for (a, b), e in acc.items():
        c = e['counts']
        label = min(c.items(), key=lambda kv: (-kv[1], order.get(kv[0], 9)))[0]
        out.append({'from': a, 'to': b, 'label': label,
                    'counts': dict(c), 'n': sum(c.values()),
                    'n_intro': e['n_intro'], 'max_group': e['max_group']})
    out.sort(key=lambda x: (-x['n'], x['from'], x['to']))
    return out


# ---------------------------------------------------------------- RSPP / 不确定性

def _clade_key(node):
    """节点名下全部叶名（排序元组）——bootstrap 复本之间靠它对齐同一个 clade。"""
    out = []

    def _rec(n):
        if not n['children']:
            out.append(n['name'])
            return
        for c in n['children']:
            _rec(c)

    _rec(node)
    return tuple(sorted(out))


def _quantile(vals, p):
    """线性插值分位数（numpy 口径）；空返回 None。"""
    if not vals:
        return None
    vs = sorted(vals)
    k = (len(vs) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(vs) - 1)
    return round(vs[lo] + (vs[hi] - vs[lo]) * (k - lo), 3)


def rssp_bootstrap(aln, states, work_dir, method='nj', n_boot=100, seed=0,
                   progress=None, observed=None):
    """根状态不确定性（RSPP）：bootstrap 重建 + Fitch，量化「内部节点区划」
    与「迁移数」的不确定性 —— 补 AUDIT 记的「迁移计数无置信区间」。

    做法：对**比对列**有放回重抽样 → 重建树 → 在**同一套 tip 区划**下重跑 Fitch：
      · 每个内部节点按**后代叶集合（clade）**跨复本对齐 → 各区划出现频率 ⇒ 后验式概率
      · 每个复本得一个迁移数 ⇒ 迁移数的经验置信区间（2.5% / 97.5%）

    为什么按 clade 而不是节点序号对齐：bootstrap 之间拓扑不同，节点序号毫无对应；
    只有「后代叶集合」是跨拓扑可比的身份。

    ⚠️ 这不等于贝叶斯后验概率：零分布来自**重抽样 + 简约法**，没有模型、没有先验。
    定位是「给点估计配一个不确定性」，别当成 BEAST 的后验。

    返回 dict：n_boot / n_ok / failed / transitions{mean,sd,p2_5,p97_5,min,max} /
    corridor_support[]（各方向对 from→to 的复本支持率，按支持复本数降序）/
    clades[]（按出现复本数降序）。产物：work_dir/rssp/bs_*.nwk（复本树，留证）。
    """
    import random as _random
    from collections import Counter
    from Virus_Platform_Core.phylo import _run_nj, _run_fasttree

    names = list(aln.keys())
    if len(names) < 3:
        raise RuntimeError('RSPP 需要 ≥3 条序列')
    widths = {len(s) for s in aln.values()}
    if len(widths) != 1:
        raise RuntimeError(f'RSPP 需要等长（已比对）序列，当前长度有 {sorted(widths)}')
    L = widths.pop() if widths else 0
    if L < 1:
        raise RuntimeError('RSPP：比对长度为 0')

    bs_dir = os.path.join(work_dir, 'rssp')
    os.makedirs(bs_dir, exist_ok=True)
    rng = _random.Random(seed)
    counts = []
    clade_hits = {}
    # 走廊（from→to）支持率：每个复本里出现过的迁移**方向对**各记一次。
    # 出处：polio-wpv1-phylodynamics 的 transition_support —— 单报「迁移数」
    # 会掩盖方向不稳定：总数一致但 A→B / B→A 来回换的复本，平均下来看不出来。
    # 支持率 = 出现过该走廊的复本数 / 有效复本数，是"方向有多稳"的直接度量。
    corridor_reps = {}
    failed = 0
    for b in range(int(n_boot)):
        cols = [rng.randrange(L) for _ in range(L)]
        # ⚠️ 重抽样的比对**保留不删**：它正是复本树 bs_<i>.nwk 的输入，
        # 删掉就没法复现「这条树是怎么来的」。同时这也避免了在循环里做
        # 几十上百次 os.remove —— 在某些受管环境里删除操作会被安全策略
        # **阻塞**（不是报错，是挂住），足以让整个分析卡死（2026-09-16 实测）。
        fa = os.path.join(bs_dir, f'bs_{b}.fasta')
        with open(fa, 'w', encoding='utf-8', newline='\n') as f:
            for n, s in aln.items():
                f.write(f'>{n}\n' + ''.join(s[i] for i in cols) + '\n')
        nwk = os.path.join(bs_dir, f'bs_{b}.nwk')
        try:
            used = (_run_fasttree(fa, nwk) if method == 'fasttree'
                    else _run_nj(fa, nwk))
            rb = parse_newick(open(used, encoding='utf-8',
                                   errors='replace').read())
            st_map, _tr, nch = fitch_mugration(rb, states)
            counts.append(nch)
            for _p, _c, frm, to in _tr:            # _tr 已只含真迁移边
                corridor_reps.setdefault((frm, to), set()).add(b)
            for node in _iter_pre(rb):
                if not node['children']:
                    continue
                key = _clade_key(node)
                if len(key) >= 2:                  # 至少有 2 片叶才谈得上共祖
                    clade_hits.setdefault(key, Counter())[st_map[id(node)]] += 1
        except Exception as e:                     # 单棵复本失败不该毁掉整体
            failed += 1
            if progress:
                progress('rssp', 0.1,
                         f'复本 {b + 1} 失败（{type(e).__name__}: {e}）— 跳过')
        if progress:
            progress('rssp', 0.1 + 0.8 * (b + 1) / max(1, int(n_boot)),
                     f'bootstrap 复本 {b + 1}/{n_boot}')

    clades = []
    for key, cnt in clade_hits.items():
        tot = sum(cnt.values())
        top_state, top_n = cnt.most_common(1)[0]
        clades.append({
            'clade': list(key),
            'n_tips': len(key),
            'n_rep': tot,
            'top_state': top_state,
            'top_prob': round(top_n / tot, 3),
            'states': {k: round(v / tot, 3) for k, v in cnt.most_common()},
        })
    # 出现复本越多、clade 越大者越靠前（这些是最稳的共祖推断）
    clades.sort(key=lambda x: (-x['n_rep'], -x['n_tips']))
    mean = sum(counts) / len(counts) if counts else None
    sd = ((sum((x - mean) ** 2 for x in counts) / len(counts)) ** 0.5
          if counts and len(counts) > 1 else 0.0)
    lo, hi = _quantile(counts, 0.025), _quantile(counts, 0.975)
    # 点估计是否落在 bootstrap 区间内 —— 落外说明「全量比对的那棵树」与
    # 「重抽样得到的树族」给出不同的迁移数（常见于拓扑接近平局时），
    # 此时迁移数不稳，别只报点估计。这一条是刻意报出来的，不是异常。
    in_range = (None if (observed is None or lo is None)
                else bool(lo <= observed <= hi))
    note = None
    if in_range is False:
        note = ('全量比对的迁移数（%s）落在 bootstrap 区间 [%s, %s] 之外 —— '
                '该数据集的迁移数对拓扑选择敏感，解读时以区间为准'
                % (observed, lo, hi))
    elif not counts:
        note = '所有 bootstrap 复本都失败，未得到区间'
    # 走廊支持率：出现该方向对的复本数 / 有效复本数（只列 ≥1 个复本的）。
    # 与 `transitions.mean`（平均迁移**条数**）互补——条数均值掩盖方向摇摆。
    n_ok = len(counts)
    corridor_support = [
        {'from': a, 'to': b, 'n_rep': len(reps),
         'support': round(len(reps) / n_ok, 3) if n_ok else 0.0}
        for (a, b), reps in sorted(corridor_reps.items(),
                                   key=lambda kv: -len(kv[1]))]
    return {
        'n_boot': int(n_boot),
        'n_ok': len(counts),
        'failed': failed,
        'seed': seed,
        'method': method,
        'observed': observed,
        'observed_in_range': in_range,
        'transitions': {
            'mean': round(mean, 3) if mean is not None else None,
            'sd': round(sd, 3),
            'p2_5': lo,
            'p97_5': hi,
            'min': min(counts) if counts else None,
            'max': max(counts) if counts else None,
        },
        'corridor_support': corridor_support,
        'note': note,
        'clades': clades,
    }


# 后验权重分带：方向支持率 → 四档（弧线线宽 / 表格 / 报告共用一套阈值）。
# 阈值取 0.95 / 0.70 / 0.50：0.70 是平台 MOT 后验阈值（docs/方法对比_*），
# 0.50 = 「复本里过半」（也是弧线虚线阈值 M5）—— 三处共用同一根线，避免同一份
# 数据在图上「稳」、在表里「不稳」。
WEIGHT_BAND_DEFS = (
    ('high', '高置信', 0.95, 1.0),
    ('mid', '中高', 0.70, 0.95),
    ('weak', '偏弱', 0.50, 0.70),
    ('low', '低（复本不过半）', 0.0, 0.50),
)


def band_of(support):
    """支持率 → 分带 key（边界含下端：0.95 属 high、0.70 属 mid、0.50 属 weak）。
    `None`（未评估）返回 None，**不当作 0** —— 没测过 ≠ 测出 0。"""
    if support is None:
        return None
    for key, _label, lo, _hi in WEIGHT_BAND_DEFS:   # 阈值降序，取首个命中
        if support >= lo:
            return key
    return None      # 非有限值（NaN 等）→ 调用方按最低档兜底


def posterior_weight_bands(corridor_support, counts,
                           driver='rssp.corridor_support'):
    """按**后验式权重**给走廊分带（驱动源：RSPP bootstrap 方向支持率）。

    为什么要有这一层：`corridor_support` 只回答「这条方向对在复本里多稳」，
    `matrix` 只回答「树上有几条边」—— 两者分开看会得出相反印象（支持率 100%
    但只有 1 条边、或 20 条边但方向来回换）。分带把两者合到**一个权重**上：
    `weight = support`（稳定性）、`exp_events = n_events × support`（按稳定性
    折减后的期望事件数），再按 weight 切档 —— 这就是 MAPLE 弧线「粗细 = 后验
    权重」的本地版（我们无 BSSVS 的 BF，驱动源换成 RSPP 支持率）。

    ⚠️ **不是贝叶斯后验**：零分布来自列重抽样 + 简约法，没有模型/先验
    （同 rssp_bootstrap 的告诫）。`driver_label` 与前端文案必须沿用同一措辞。

    corridor_support: rssp_bootstrap 的 corridor_support（[{from,to,support,…}]）
    counts: {(from, to): n_events}（Fitch 重构的跨区迁移计数，即 analyze 的 matrix）

    返回 dict：bands[]（四档全给，空档也在，前端表格布局才不跳）、
    n_corr / n_events / exp_events 合计、n_zero_rep（树上出现但复本里一次都没
    出现的走廊数 —— 它们是**最不稳**的一档，单独计数免被"未评估"混同）、
    n_uneval（列了方向对却没给支持率 → **未评估**，不进分带也不当 0）、
    corridors[]（全部走廊按期望事件数降序）、driver / driver_label / note。
    """
    sup = {}
    for c in corridor_support or []:
        s = c.get('support')
        sup[(c.get('from'), c.get('to'))] = \
            (None if s is None else float(s))
    corridors = []
    n_zero_rep = 0
    n_uneval = 0
    for (a, b), n in (counts or {}).items():
        if (a, b) in sup:
            s = sup[(a, b)]
            if s is None:      # 列了方向对、没给支持率＝未评估（≠ 测出 0）
                n_uneval += 1
                continue
            zero = False
        else:
            # 树上出现、复本里 0 次 → 支持率 0（这是**测量结果**，不是"未测"）
            s, zero = 0.0, True
            n_zero_rep += 1
        corridors.append({
            'from': a, 'to': b, 'support': s, 'n_events': int(n),
            'exp_events': round(int(n) * s, 3),
            'band': band_of(s) or 'low',
            'zero_rep': zero,
        })
    corridors.sort(key=lambda c: (-c['exp_events'], -c['n_events'],
                                  str(c['from']), str(c['to'])))
    bands = []
    for key, label, lo, hi in WEIGHT_BAND_DEFS:
        rows = [c for c in corridors if c['band'] == key]
        bands.append({
            'band': key, 'label': label, 'lo': lo, 'hi': hi,
            'n_corr': len(rows),
            'n_events': sum(c['n_events'] for c in rows),
            'exp_events': round(sum(c['exp_events'] for c in rows), 3),
            'corridors': rows,
        })
    note = None
    if n_zero_rep or n_uneval:
        parts = []
        if n_zero_rep:
            parts.append('%d 条走廊在树上出现、但在任何 bootstrap 复本里都没出现'
                         '（支持率 0，已归入最低档）——这类方向对最不稳，引用前'
                         '先看 leaf_states/transitions 的原始边。' % n_zero_rep)
        if n_uneval:
            parts.append('%d 条走廊列在方向对表里但没有支持率（未评估），'
                         '未纳入分带。' % n_uneval)
        note = ' '.join(parts)
    return {
        'driver': driver,
        'driver_label': 'RSPP bootstrap 方向支持率（非贝叶斯后验）',
        'edges': [d[3] for d in WEIGHT_BAND_DEFS],
        'bands': bands,
        'corridors': corridors,
        'n_corr': len(corridors),
        'n_events': sum(c['n_events'] for c in corridors),
        'exp_events': round(sum(c['exp_events'] for c in corridors), 3),
        'n_zero_rep': n_zero_rep,
        'n_uneval': n_uneval,
        'note': note,
        # 期望事件数这一列的口径（前端脚注/报告直接取用）：两个 driver 下
        # 这个列算的不是同一件事，措辞必须跟着 driver 走，不能一套话打天下。
        'exp_note': ('期望事件数 = 事件数 × 支持率（按方向稳定性折减）；'
                     '支持率 = 该方向对出现的 bootstrap 复本比例。'),
    }


# BEAST 驱动版分带：驱动量来自**真实 BEAST 产物**（三种口径，见 beast_io 模块头）。
# 阈值仍与 WEIGHT_BAND_DEFS 同一套（0.95/0.70/0.50）—— 同一根线在四处含义一致；
# 只有最低档的说明词按口径分三种：「样本树不过半」「后验不过半」「枝端后验不过半」
# 不是同一件事，标签必须跟着口径走（前端 i18n 也按 driver 分三套键）。
BEAST_BAND_DEFS = (
    ('high', '高置信', 0.95, 1.0),
    ('mid', '中高', 0.70, 0.95),
    ('weak', '偏弱', 0.50, 0.70),
    ('low', '低（后验不过半）', 0.0, 0.50),
)

# 三种口径各自的措辞（driver_label / 最低档标签 / 期望事件数口径 / 告诫）。
# ⚠️ 三者措辞**不许互串**：tree_sample 是后验频率、jump 是后验跳转计数、
# mcc 是保守下界（不是该走廊的后验概率）。前端 bandLabel 按 driver 取键。
BEAST_DRIVERS = {
    'beast.tree_sample': {
        'label': 'BEAST 后验树样本走廊频率（逐树是否含该走廊）',
        'low': '低（样本树不过半）',
        'unit': '棵样本树',
        'support_note': '支持率 = 含该走廊的后验样本树占比',
        'exp_note': ('期望事件数 = 逐树事件数的均值（每棵树先数事件、再对树'
                     '取均值）；支持率 = 含该走廊的样本树比例。'),
        'caveat': ('走廊按每棵树各自的节点状态**最可能值**判定，不是对每条枝的'
                   '后验积分；同一枝上 a→b→a 的来回按各次分别计。'),
    },
    'beast.jump': {
        'label': 'BEAST 后验跳转计数 P(跳转 > 0)',
        'low': '低（后验不过半）',
        'unit': '个后验采样点',
        'support_note': '支持率 = P(该走廊跳转数 > 0)',
        'exp_note': ('期望事件数 = 后验跳转计数均值（BEAST 在采样历史上积分的'
                     '期望跳转次数，含多跳路径）；支持率 = 跳转数 > 0 的'
                     '后验样本比例。'),
        'caveat': ('跳转计数是**采样历史下的期望跳转次数**，与 Fitch 观测事件数'
                   '是两条口径；本平台不算贝叶斯因子、不做模型比较。'),
    },
    'beast.mcc': {
        'label': 'BEAST MCC 树枝端后验下界（非该走廊的后验概率）',
        'low': '低（枝端后验不过半）',
        'unit': '棵 MCC 树（单树，无后验样本）',
        'support_note': '支持率 = 枝端节点后验的**下界**',
        'exp_note': ('期望事件数 = Σ P(父=a)·P(子=b)（枝端后验乘积，端独立'
                     '近似）；支持率 = 承载该走廊的各枝端节点后验的**下界**'
                     '（其中最不确定的那个），真值 ≥ 该数。'),
        'caveat': ('单棵树只给**保守下界**：真值 ≥ 该数，不是该走廊的后验概率；'
                   '没有逐样本 → 没有可信区间、也不做逐树动画'
                   '（要这些请导入 .trees 或 .log）。'),
    },
}


def _normalize_beast_corridors(corridors):
    """三种口径的走廊表 → 统一形（support 优先取 support、其次 p_used）。"""
    out = []
    for c in corridors or []:
        s = c.get('support')
        if s is None:
            s = c.get('p_used')
        if s is None:
            continue
        out.append({
            'from': c.get('from'), 'to': c.get('to'),
            'support': round(float(s), 4),
            'mean_events': round(float(c.get('mean_events') or 0.0), 4),
            'ci95': c.get('ci95'), 'max_events': c.get('max_events'),
            'n_events': c.get('n_events'),
        })
    return out


def beast_weight_bands(driver, corridors, counts=None, n_units=0):
    """按**BEAST 驱动量**给走廊分带（口径由 driver 决定，见 BEAST_DRIVERS）。

    走廊全集 = BEAST 表 ∪ Fitch 观测，缺的一侧按 0 计并计入
    n_zero_rep / n_extra，两边都不丢；驱动量口径三选一：
      `beast.tree_sample` 支持率 = 含该走廊的样本树占比（真后验频率）
      `beast.jump`        支持率 = P(跳转数 > 0)（后验跳转计数）
      `beast.mcc`         支持率 = 枝端节点后验的**下界**（保守下界）

    ⚠️ 与 RSPP 版的关键差别：期望事件数**直取**各口径的 mean_events，
    不再乘支持率 —— 这三个量本身就是按事件数定义的（RSPP 那里 n_events
    是 Fitch 边数、必须乘支持率折减，口径不同，别混用话术）。

    driver: 三者之一；counts: {(from, to): n_events}（Fitch 计数，即 matrix）
    n_units: 该口径的样本单位数（样本树数 / 采样点数 / 1），仅供展示
    """
    info = BEAST_DRIVERS.get(driver)
    if not info:
        raise RuntimeError(f'未知的 BEAST 驱动口径：{driver}')
    defs = tuple((k, (info['low'] if k == 'low' else lab), lo, hi)
                 for k, lab, lo, hi in BEAST_BAND_DEFS)
    post = {}
    for c in corridors or []:
        post[(c.get('from'), c.get('to'))] = c
    obs = dict(counts or {})
    corridors_out = []
    n_zero_rep = 0
    n_extra = 0
    for key in set(post) | set(obs):
        a, b = key
        c = post.get(key)
        n_obs = int(obs.get(key, 0) or 0)
        if c is None:
            # BEAST 产物里一次没出现 → 支持率 0（**测量结果**，不是"未测"）
            p, mean_ev, zero, ci, mx = 0.0, 0.0, True, None, None
            n_zero_rep += 1
        else:
            p = float(c.get('support') or 0.0)
            mean_ev = float(c.get('mean_events') or 0.0)
            zero = (p <= 0.0)
            ci, mx = c.get('ci95'), c.get('max_events')
        if n_obs <= 0:
            n_extra += 1
        corridors_out.append({
            'from': a, 'to': b,
            # 键名沿用 'support'（各 driver 的驱动量同名同位），口径由
            # driver_label / support_note 说明。
            'support': round(p, 4),
            'n_events': n_obs,
            'exp_events': round(mean_ev, 3),
            'band': band_of(p) or 'low',
            'zero_rep': zero,
            'ci95': ci, 'max_events': mx,
        })
    corridors_out.sort(key=lambda c: (-c['exp_events'], -c['support'],
                                      str(c['from']), str(c['to'])))
    bands = []
    for key, label, lo, hi in defs:
        rows = [c for c in corridors_out if c['band'] == key]
        bands.append({
            'band': key, 'label': label, 'lo': lo, 'hi': hi,
            'n_corr': len(rows),
            'n_events': sum(c['n_events'] for c in rows),
            'exp_events': round(sum(c['exp_events'] for c in rows), 3),
            'corridors': rows,
        })
    parts = []
    if n_zero_rep:
        parts.append('%d 条走廊在 Fitch 重构里出现、但在 BEAST 产物里一次都'
                     '没出现（支持率 0，已归入最低档）。' % n_zero_rep)
    if n_extra:
        parts.append('%d 条方向对在 BEAST 产物里出现、但 Fitch 最简约重构没给出'
                     '（事件数列为 0，期望事件数来自 BEAST 口径）——两者不矛盾：'
                     '简约法不给"多一跳"的历史，模型法按速率先验允许它出现。'
                     % n_extra)
    return {
        'driver': driver,
        'driver_label': info['label'],
        'unit': info['unit'],
        'n_units': int(n_units or 0),
        'edges': [d[3] for d in defs],
        'bands': bands,
        'corridors': corridors_out,
        'n_corr': len(corridors_out),
        'n_events': sum(c['n_events'] for c in corridors_out),
        'exp_events': round(sum(c['exp_events'] for c in corridors_out), 3),
        'n_zero_rep': n_zero_rep,
        'n_extra': n_extra,
        'n_uneval': 0,
        'support_note': info['support_note'],
        'exp_note': info['exp_note'],
        'caveat': info['caveat'],
        'note': (' '.join(parts) if parts else None),
    }


# ---------------------------------------------------------------- BEAST 导入

# summary 里逐列诊断的截断上限（全量始终落 beast_log.json）—— 与 TC_MAX /
# SP_MAX 同一纪律：截断就如实报总数，不静默。
_BEAST_COL_MAX = 40


def _beast_state_check(beast_states, meta_states):
    """BEAST 状态名 vs 元数据区划名：两侧清单 + 交集 + 仅单侧（不猜同义）。"""
    bs = sorted({s for s in (beast_states or []) if s})
    ms = sorted({s for s in (meta_states or []) if s})
    inter, only_b, only_m = [], [], []
    for s in bs:
        (inter if s in ms else only_b).append(s)
    for s in ms:
        if s not in bs:
            only_m.append(s)
    return {'beast_states': bs, 'meta_states': ms, 'matched': inter,
            'only_in_beast': only_b, 'only_in_meta': only_m,
            'n_matched': len(inter), 'aligned': bool(bs) and bs == ms}


def _json_safe(obj):
    """NaN/±Inf → null（这两个产物是**原样发文件**给浏览器的：
    `r.json()`（JSON.parse）遇到裸 NaN 直接报错 → 动画静默不出结果。
    与 tool_jobs._dump_json 同一理由，这里给独立导入工具也备一份）。"""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def import_beast(beast_log=None, beast_trees=None, burnin=0.1, min_prob=0.0,
                 max_trees=2000, out_dir=None, meta_states=(), progress=None):
    """导入 BEAST 产物（.log / .trees / MCC 树）→ (beast_summary, driver, corridors)。

    - `.log`   → 逐列后验 + ESS 诊断（log_summary）+ 跳转走廊表（jump_corridors）
    - `.trees` → ≥2 棵树＝后验样本（mode='sample'：逐树事件、走廊频率、
                 时间直方图、动画数据）；1 棵树＝MCC（mode='mcc'：枝端后验
                 下界，保守口径）
    全量落 `<out_dir>/beast_log.json`、`<out_dir>/beast_trees.json`；
    beast_summary 是给 summary.json 的精简版（逐列诊断截断 `_BEAST_COL_MAX`、
    逐树事件不放进 summary —— 前端按需拉产物文件）。

    驱动优先序（同一份运行里只用**一个** driver 给弧线分带，否则「粗」在图上
    不再可比）：`tree_sample` > `jump` > `mcc`。真实 BEAST 产物优先于 RSPP
    —— 它来自用户为这次分析跑的真实模型（见 analyze 的分带块）。

    返回 driver=None 表示三类口径都拿不到**或状态名与元数据对不上**
    （后者见 state_check；弧线分带仍退回 RSPP，BEAST 面板原样展示）。
    """
    from Virus_Platform_Core import beast_io
    def _p(f, m):
        if progress:
            progress(f, m)

    out = {'log': None, 'jump': None, 'jump_error': None, 'trees': None,
           'driver': None, 'driver_label': None, 'n_units': 0,
           'state_check': None, 'artifacts': [], 'warnings': [],
           'log_file': (os.path.basename(beast_log) if beast_log else None),
           'trees_file': (os.path.basename(beast_trees) if beast_trees
                          else None)}
    corr = []
    driver = None

    if beast_log:
        _p(0.93, '导入 BEAST 日志（逐列诊断 + 跳转计数）')
        ls = beast_io.log_summary(beast_log, burnin=burnin)
        out['warnings'] += list(ls.get('warnings') or [])
        try:
            jc = beast_io.jump_corridors(beast_log, burnin=burnin)
            out['warnings'] += list(jc.get('warnings') or [])
        except RuntimeError as e:
            # 只有扁平 `*indicators*` 列时的明确拒绝理由（不猜顺序）
            jc = None
            out['jump_error'] = str(e)
            out['warnings'].append('BEAST 日志里没有按走廊命名的跳转计数列：'
                                   + str(e))
        cols = ls.get('columns') or []
        out['log'] = {
            'file': ls.get('file'), 'n_samples': ls.get('n_samples'),
            'n_skipped': ls.get('n_skipped'), 'n_total': ls.get('n_total'),
            'burnin': ls.get('burnin'), 'ci': ls.get('ci'),
            'header_kind': ls.get('header_kind'), 'reader': ls.get('reader'),
            'n_columns': len(cols),
            'columns': cols[:_BEAST_COL_MAX],
            'columns_truncated': len(cols) > _BEAST_COL_MAX,
            'group_counts': {k: len(v) for k, v in
                             (ls.get('groups') or {}).items() if v},
            'tree_height': ls.get('tree_height'),
            'ess_min': ls.get('ess_min'), 'n_low_ess': ls.get('n_low_ess'),
            'n_nonnumeric': ls.get('n_nonnumeric'),
            'warnings': list(ls.get('warnings') or []),
        }
        if jc is not None:
            out['jump'] = {
                'file': jc.get('file'), 'n_samples': jc.get('n_samples'),
                'burnin': jc.get('burnin'), 'ci': jc.get('ci'),
                'total_mean': jc.get('total_mean'),
                'total_source': jc.get('total_source'),
                'corridors': jc.get('corridors'),
                'dwell': jc.get('dwell'),
                'warnings': list(jc.get('warnings') or []),
            }
        if out_dir:
            p = os.path.join(out_dir, 'beast_log.json')
            with open(p, 'w', encoding='utf-8', newline='\n') as f:
                json.dump(_json_safe({'file': os.path.basename(beast_log),
                                      'burnin': float(burnin),
                                      'log': ls, 'jump': jc}), f,
                          ensure_ascii=False, indent=1)
            out['artifacts'].append('beast_log.json')

    if beast_trees:
        # 无 log 的树导入 = 未做收敛诊断（审查 P1-10）：后验树自相关/有效
        # 样本量无法从树数本身判断，走廊频率只能当探索性结果。显式落进
        # warnings 与 trees 摘要，防止被读成"稳定的后验支持"。
        if not beast_log:
            out['warnings'].append(
                '仅导入 BEAST 树（未提供 .log）：未做 ESS/收敛诊断，'
                '走廊/频率结果仅为探索性，正式报告需配套日志诊断')
        _p(0.95, '导入 BEAST 树（后验样本 / MCC 单树）')
        tb = beast_io.tree_corridors(
            beast_trees, burnin=burnin, max_trees=max_trees,
            min_prob=min_prob,
            progress=(lambda f, m: _p(0.95 + 0.03 * min(f, 1.0), m)))
        out['warnings'] += list(tb.get('warnings') or [])
        out['trees'] = {k: tb.get(k) for k in (
            'file', 'mode', 'n_trees', 'n_trees_used', 'n_branches',
            'n_branches_per_tree', 'n_unresolved', 'n_same_state',
            'n_lowprob', 'labels', 'corridors', 'tree_height', 'time_hist',
            'anim_truncated', 'max_anim_trees', 'min_prob', 'warnings')}
        # P1-10：收敛是否被核验（有 log 才有 ESS 诊断），前端/报告可见
        out['trees']['convergence_checked'] = bool(beast_log)
        # 逐树事件（动画数据）不放进 summary：前端按需拉产物文件
        out['trees']['n_anim_trees'] = len(tb.get('sample_events') or [])
        out['trees']['has_events'] = bool(tb.get('sample_events'))
        if out_dir:
            p = os.path.join(out_dir, 'beast_trees.json')
            with open(p, 'w', encoding='utf-8', newline='\n') as f:
                json.dump(_json_safe({'file': os.path.basename(beast_trees),
                                      'burnin': float(burnin),
                                      'min_prob': float(min_prob), **tb}), f,
                          ensure_ascii=False, indent=1)
            out['artifacts'].append('beast_trees.json')

    # ---- 状态名核对 + 驱动选择（三口径优先序，见 docstring）----
    bstates = []
    if out['trees']:
        bstates = list(out['trees'].get('labels') or [])
    if not bstates and out['jump']:
        bstates = sorted({c['from'] for c in out['jump']['corridors']}
                         | {c['to'] for c in out['jump']['corridors']})
    out['state_check'] = _beast_state_check(bstates, meta_states)
    if ((out['trees'] or out['jump'])
            and out['state_check']['meta_states']
            and not out['state_check']['n_matched']):
        out['warnings'].append(
            'BEAST 产物里的状态名（%s）与元数据区划（%s）**没有交集** —— '
            '两侧的走廊表按各自命名展示；弧线分带不采用 BEAST 驱动'
            '（覆盖范围对不上）。若是同一套区划的不同写法，请统一命名后重导。'
            % ('/'.join(out['state_check']['beast_states'][:8]),
               '/'.join(out['state_check']['meta_states'][:8])))
    elif (out['state_check']['n_matched']
          and not out['state_check']['aligned']):
        out['warnings'].append(
            'BEAST 状态名与元数据区划只是部分一致：BEAST 多出 %s、元数据多出 %s '
            '—— 只交集部分参与分带对齐，两侧多出的走廊按各自命名展示。'
            % ('/'.join(out['state_check']['only_in_beast'][:6]) or '无',
               '/'.join(out['state_check']['only_in_meta'][:6]) or '无'))

    # 没有元数据可对齐时（独立导入工具）照样给 driver —— 它只用来标口径；
    # 有元数据但一条都对不上时 driver=None（覆盖范围对不上，弧线分带会误导）。
    if out['state_check']['n_matched'] or not out['state_check']['meta_states']:
        # 优先序：树样本（真后验频率，含拓扑不确定性）> 跳转计数 > MCC 单树
        if out['trees'] and out['trees'].get('corridors') \
                and out['trees'].get('mode') == 'sample':
            driver = 'beast.tree_sample'
            corr = out['trees']['corridors']
            out['n_units'] = int(out['trees'].get('n_trees_used') or 0)
        elif out['jump'] and out['jump'].get('corridors'):
            driver = 'beast.jump'
            corr = out['jump']['corridors']
            out['n_units'] = int(out['jump'].get('n_samples') or 0)
        elif out['trees'] and out['trees'].get('corridors'):
            driver = 'beast.mcc'
            corr = out['trees']['corridors']
            out['n_units'] = 1
    out['driver'] = driver
    out['driver_label'] = (BEAST_DRIVERS[driver]['label'] if driver else None)
    out['n_corr'] = len(corr or [])
    out['support_kind'] = driver
    # 口径措辞**原样带出**：前端面板/PDF 报告直接渲染这三个字段（不另写一套
    # 文案）。缺了它们界面会印出"（；）"并把最关键的那句告诫（MCC 单树只给
    # 保守下界）吞掉 —— 措辞纪律要求这句话必须跟在数字旁边。
    _di = BEAST_DRIVERS.get(driver) or {}
    out['support_note'] = _di.get('support_note')
    out['exp_note'] = _di.get('exp_note')
    out['caveat'] = _di.get('caveat')
    out['unit'] = _di.get('unit')
    return out, driver, _normalize_beast_corridors(corr)


# ⚠️ 2026-10-03 删：`ltt_and_skyline`（LSD2 定年树的 LTT + 经典 skyline 点估计）
#    全仓**无任何代码引用** —— A1 于 2026-09-17 收敛到 TreeTime 后，LSD2 通路
#    未接线，LTT 的唯一实现是 `phylodyn_trees.ltt_from_tree`（t-pdltt 卡在用）。
#    该函数只被 docs/*.md 的历史对比文档提到过；保留它会让"哪套 LTT 才是活的"
#    继续含糊。若将来恢复 LSD2 通路（`lsd2_dating` 仍在），按 git 历史取回即可。


# ---------------------------------------------------------------- MOTP 迁移随时间

def parse_lsd2_dated_tree(path):
    """从 LSD2 的 `<tree>.result.date.nexus` 取**带节点年代**的树。

    NEXUS 的 tree 行形如 `tree 1 = (...)[&date="1998.33"];` —— **每个节点（含内部）
    都带 `[&date="..."]`**，所以节点年代是白拿的，不必用枝长÷速率去换算。

    实测佐证「`.result.nwk` 是**替换率尺度**而非时间尺度」：合成例里 C4 采样 2010
    距根 0.04629 = (2010 − tMRCA 1998.33) × rate 0.00396777，逐位吻合。

    返回 (root, {id(node): float 年代})；解析不出返回 (None, {})。
    """
    try:
        txt = open(path, encoding='utf-8', errors='replace').read()
    except OSError:
        return None, {}
    m = re.search(r'tree\s+\S+\s*=\s*(.+?);', txt, re.S | re.I)
    if not m:
        return None, {}
    try:
        root = parse_newick(m.group(1))
    except ValueError:
        # 候选文件里出现畸形/非 LSD2 的树（如把 BEAST mcc.nexus 指到这儿：
        # 那串以 `[&R] (((…)))` 开头）→ 按本函数契约「解析不出返回 (None, {})」
        # 继续试下一个候选，让调用方给「没有节点年代」的明确告警，而不是硬崩。
        return None, {}
    dates = {}

    def _rec(n):
        nm = n.get('name') or ''
        dm = re.search(r'\[&date="([^"]+)"\]', nm)
        if dm:
            try:
                dates[id(n)] = float(dm.group(1))
            except ValueError:
                pass
            n['name'] = nm[:dm.start()].strip()      # 剥掉注释，名字还原
        for c in n['children']:
            _rec(c)

    _rec(root)
    return (root, dates) if dates else (None, {})


def tree_root_distances(root):
    """每个节点到根的距离（枝长累加）→ {id(node): float}。

    TreeDater 的 `treedater_dated.nwk` **枝长＝年但没有节点年代注释**（年代写在
    `node_heights.csv` 里，且索引是 R 自己的 ape 编号，与 newick 节点对不上）。
    与其去猜 R 的编号，不如用枝长自己算：`root_h = max(dist)` 就是最深叶的采样
    时刻，`root_h − dist` 就是「距今多少年」，逐个节点白拿且与 R 的
    `age_before_present` 同口径。
    """
    out = {}

    def _rec(n, d):
        out[id(n)] = d
        for c in n['children']:
            try:
                ln = float(c.get('length') or 0.0)
            except (TypeError, ValueError):
                ln = 0.0
            _rec(c, d + ln)

    _rec(root, 0.0)
    return out


def write_rtt_dated_nexus(tree_path, out_nexus, slope, intercept):
    """RTT 回归（slope/intercept）→ LSD2 兼容的带年代 NEXUS。

    年代 = intercept + slope × 距根遗传距离（与 root_to_tip 同口径），
    每个节点在名字位带 `[&date="年"]` —— 与 parse_lsd2_dated_tree 的读取
    格式一致，作为「时间信号与定年 → 系统地理」的 dated_tree 格式桥。
    2026-09-21 新增（此前 RTT 卡只产 newick，年代不随树走，
    系统地理卡的 dated_tree 入参读不到年代 → 接力断点）。
    """
    txt = open(tree_path, encoding='utf-8', errors='replace').read().strip()
    m = re.search(r'tree\s+\S+\s*=\s*(.+?);', txt, re.S | re.I)
    body = (m.group(1) if m else txt).strip().rstrip(';')
    if body.startswith('[&R]'):
        body = body[4:].strip()
    root = parse_newick(body)

    dist = tree_root_distances(root)

    def _emit(n):
        y = intercept + slope * dist.get(id(n), 0.0)
        tag = '[&date="%.4f"]' % y
        head = (n.get('name') or '') + tag
        bl = ':%.6g' % (n.get('length') or 0.0)
        if n['children']:
            return '(' + ','.join(_emit(c) for c in n['children']) + ')' + head + bl
        return head + bl

    out = _emit(root)
    with open(out_nexus, 'w', encoding='utf-8') as f:
        f.write('#NEXUS\nBEGIN TREES;\n\ttree 1 = ' + out + ';\nEND;\n')
    return out_nexus


def time_tree_segments(root, x_of, axis_label='采样年（日历年）', source=''):
    """（定年）树 → **时间轴直角树**的 plotly 单 trace 折线坐标（T3 时间树）。

    x_of: {id(node): 时间坐标}（LSD2 用 `parse_lsd2_dated_tree` 的年代；TreeDater
    用 `tree_root_distances` 换出的「距今年」）。与 `ltt_from_tree` 同一套输入
    约定 —— 传进来的必须是**同一棵树对象**上算出的字典（id() 才有效）。

    形态：叶按 DFS 顺序占 y=0..n-1，内部节点 y = 子节点 y 均值；每条边画
    「父处竖线 + 到子处横线」，全部塞进一个 'lines' trace（段间用 None 断开，
    这是 plotly 的标准做法，比逐个节点画省 trace 数）。

    ⚠️ 内部节点缺时间坐标时**不静默**：回退用父节点坐标（视觉上＝零长枝），根节点
    缺坐标才退到子节点均值，计数进 `n_backfilled` 并写进 note —— 定年工具本该给全
    节点时间，缺了就是异常，允许画但必须在图上说明，不能假装是定年结果。

    返回 dict：x/y（plotly trace 坐标）/ tips[{name,x,y}] / min_x / max_x /
    n_tips / n_internal / n_backfilled（缺时间坐标的节点总数）/ n_back_tips
    （其中是叶的条数，叶缺时间比内部节点严重，单独计数并在 note 里点名）/
    axis_label / source / note。
    """
    tips, internal = [], []
    for n in _iter_pre(root):
        (tips if not n['children'] else internal).append(n)
    if len(tips) < 3:
        return {'x': [], 'y': [], 'tips': [], 'n_tips': len(tips),
                'n_internal': 0, 'n_backfilled': 0, 'note': '叶数不足（<3），不画时间树'}
    if len(x_of) < max(3, len(tips) // 2):
        return {'x': [], 'y': [], 'tips': [], 'n_tips': len(tips),
                'n_internal': 0, 'n_backfilled': 0,
                'note': ('带时间坐标的节点太少'
                         f'（{len(x_of)}/{len(tips) + len(internal)}）——'
                         '先确认定年结果里是否真有节点年代')}

    y_of = {id(t): float(i) for i, t in enumerate(tips)}

    def _y(n):
        if id(n) not in y_of:
            vals = [_y(c) for c in n['children']]
            y_of[id(n)] = sum(vals) / len(vals) if vals else 0.0
        return y_of[id(n)]

    _y(root)
    n_back = [0]

    def _x(n, parent_x):
        v = x_of.get(id(n))
        if v is not None:
            return float(v)
        if parent_x is not None:
            return parent_x                          # 缺坐标 → 塌成零长枝
        kids = [_x(c, None) for c in n['children']]
        return sum(kids) / len(kids) if kids else None

    xs = {}
    n_back_tip = [0]

    def _walk(n, px):
        xs[id(n)] = _x(n, px)
        if x_of.get(id(n)) is None:                  # 计数只在这里做一次
            n_back[0] += 1
            if not n['children']:
                n_back_tip[0] += 1                   # 叶也缺时间 → 必须单独说清
        for c in n['children']:
            _walk(c, xs[id(n)])

    _walk(root, None)

    X, Y = [], []

    def _emit(node):
        y0 = y_of[id(node)]
        x0 = xs[id(node)]
        for c in node['children']:
            y1 = y_of[id(c)]
            x1 = xs[id(c)]
            if x0 is None or x1 is None:
                continue
            if y1 != y0:
                X.extend([x0, x0, None])
                Y.extend([y0, y1, None])
            X.extend([x0, x1, None])
            Y.extend([y1, y1, None])
            _emit(c)

    _emit(root)
    tip_rows = [{'name': (t['name'] or ''),
                 'x': round(xs[id(t)], 6) if xs.get(id(t)) is not None else None,
                 'y': y_of[id(t)]} for t in tips]
    vals = [v for v in xs.values() if v is not None]
    n_internal_back = n_back[0] - n_back_tip[0]
    if not n_back[0]:
        _scale = '坐标全部取自定年结果的节点时间。'
    else:
        _parts = []
        if n_internal_back:
            _parts.append(f'其中 {n_internal_back} 个内部节点无时间坐标，'
                          '已按父节点位置塌成零长枝。')
        if n_back_tip[0]:
            _parts.append(f'⚠ {n_back_tip[0]} 个叶（tip）没有时间坐标（定年结果缺'
                          '这些样品的年代）——它们的时间位置按父节点塌陷，别当真实'
                          '采样时刻解读。')
        _scale = ''.join(_parts)
    note = (f'时间树（{axis_label}）：{len(tips)} 叶，' + _scale + '未做任何额外插值。')
    return {
        'x': [round(v, 6) if v is not None else None for v in X],
        'y': Y,
        'tips': tip_rows,
        'min_x': min(vals) if vals else None,
        'max_x': max(vals) if vals else None,
        'n_tips': len(tips), 'n_internal': len(internal),
        'n_backfilled': n_back[0],
        'n_back_tips': n_back_tip[0],
        'axis_label': axis_label, 'source': source, 'note': note,
    }


def leaf_years_from_states(root, years):
    """给每个节点算「后代叶采样年」的中位/范围（无定年树时的近似年代）。

    years: {叶名: int 年份}。返回 (median_by_id, span_by_id)。
    """
    med, span = {}, {}

    def _rec(n):
        if not n['children']:
            y = years.get(n['name'])
            vals = [y] if y else []
        else:
            vals = []
            for c in n['children']:
                _rec(c)
                vals += _leaf_vals[id(c)]
        _leaf_vals[id(n)] = vals
        if vals:
            vs = sorted(vals)
            med[id(n)] = vs[len(vs) // 2]
            span[id(n)] = (vs[0], vs[-1])
        return vals

    _leaf_vals = {}
    _rec(root)
    return med, span


def migration_over_time(root, states, node_dates=None, years=None, bin_width=1.0,
                        bin_start=None):
    """MOTP：把 Fitch 的每条转移按**子节点年代**分箱 → 随时间变化的转移计数。

    时间轴来源（按优先级，取到什么就在返回值里注明用的是哪种）：
      1. `node_dates`（LSD2 定年树 `[&date=...]` 的节点年代）—— 严格口径
      2. `years`（叶采样年）→ 用**后代叶采样年中位**近似节点年代 —— 近似口径
         （⚠️ 这不是节点年代，只是「这些叶什么时候采的」的中点；要严格年代先做 LSD2）

    返回 dict：axis（'dated_tree' | 'tip_years' | None）/ bin_width /
    bins[{lo,hi,n}]（按时间升序）/ rows[{lo,hi,from,to,count}] / n_undated /
    series[{from,to,points[[bin_center,count]...]}]（供前端画堆叠折线）/
    flow_series{centers,total,regions[{region,in,out,net}]}（区划视角净流量）。
    """
    _st, transitions, child_ids = _fitch(root, states)
    med, span = ({}, {})
    if node_dates:
        axis = 'dated_tree'
        t_of = lambda cid: node_dates.get(cid)          # noqa: E731
    elif years:
        axis = 'tip_years'
        med, span = leaf_years_from_states(root, years)
        t_of = lambda cid: med.get(cid)                 # noqa: E731
    else:
        return {'axis': None, 'bin_width': bin_width, 'bins': [], 'rows': [],
                'n_undated': len(transitions), 'series': [], 'flow_series': None,
                'note': '既没有定年树也没有采样年，无法按时间分箱'}

    timed = []
    n_undated = 0
    for (pname, cname, frm, to), cid in zip(transitions, child_ids):
        t = t_of(cid)
        if t is None:
            n_undated += 1
            continue
        timed.append((float(t), pname, cname, frm, to))
    if not timed:
        return {'axis': axis, 'bin_width': bin_width, 'bins': [], 'rows': [],
                'n_undated': n_undated, 'series': [], 'flow_series': None,
                'note': '所有转移都没有可用时间，未分箱'}

    tmin = min(t for t, *_ in timed)
    tmax = max(t for t, *_ in timed)
    bw = float(bin_width) if bin_width and float(bin_width) > 0 else 1.0
    lo0 = float(bin_start) if bin_start else math.floor(tmin / bw) * bw
    if tmax < lo0:                                       # 时间轴倒退时的兜底
        lo0 = math.floor(tmax / bw) * bw

    def _bin_of(t):
        return int(math.floor((t - lo0) / bw))

    nbin = max(1, _bin_of(tmax) + 1)
    rows = Counter()
    bins = [0] * nbin
    for t, _p, _c, frm, to in timed:
        b = max(0, min(_bin_of(t), nbin - 1))
        rows[(b, frm, to)] += 1
        bins[b] += 1

    rows_out = [{'lo': round(lo0 + b * bw, 3), 'hi': round(lo0 + (b + 1) * bw, 3),
                 'from': frm, 'to': to, 'count': n}
                for (b, frm, to), n in sorted(rows.items(),
                                              key=lambda x: (x[0][0], -x[1]))]
    bins_out = [{'lo': round(lo0 + i * bw, 3), 'hi': round(lo0 + (i + 1) * bw, 3),
                 'n': bins[i]} for i in range(nbin)]
    pairs = sorted({(r['from'], r['to']) for r in rows_out})
    points_by_pair = {pr: [] for pr in pairs}
    for r in rows_out:
        points_by_pair[(r['from'], r['to'])].append(
            [round((r['lo'] + r['hi']) / 2, 3), r['count']])
    series = [{'from': a, 'to': b, 'points': pts}
              for (a, b), pts in points_by_pair.items()]
    # 净流量时间曲线：每个区划逐时间窗的「进 / 出 / 净（进−出）」。
    # 与 series 互补 —— series 是路线视角（谁迁到谁），这是区划视角（谁在净输出、
    # 谁在净接收，随时间怎么变），都从同一份 rows 派生，口径必然一致。
    regions_flow = sorted({r['from'] for r in rows_out} | {r['to'] for r in rows_out})
    _in = {r: [0] * nbin for r in regions_flow}
    _out = {r: [0] * nbin for r in regions_flow}
    for (b, frm, to), n in rows.items():
        _in[to][b] += n
        _out[frm][b] += n
    flow_series = {
        'centers': [round((bins_out[i]['lo'] + bins_out[i]['hi']) / 2, 3)
                    for i in range(nbin)],
        'total': list(bins),
        'regions': [{'region': r, 'in': _in[r], 'out': _out[r],
                     'net': [_in[r][i] - _out[r][i] for i in range(nbin)]}
                    for r in regions_flow],
    }
    note = None
    if axis == 'tip_years':
        note = ('时间轴用的是**后代叶采样年的中位**（近似），不是节点年代 —— '
                '要严格年代请先跑「时间信号与定年」里的 LSD2，再把 `.date.nexus` 填进来。')
    return {'axis': axis, 'bin_width': bw, 'bins': bins_out, 'rows': rows_out,
            'n_undated': n_undated, 'series': series, 'note': note,
            'flow_series': flow_series,
            'n_transitions': len(transitions), 'n_timed': len(timed)}




def resolve_dated_tree(path):
    """用户给的定年树文件 → 找到**真带节点年代**那个。返回 (path, root, dates)。

    为什么需要：LSD2 一次跑出三个同名兄弟文件 ——
      `<x>.result`（文本结果，含 rate/tMRCA）
      `<x>.result.nwk`（树，但是**替换率尺度**，没有年代注释）
      `<x>.result.date.nexus`（**时间树**，每个节点带 `[&date="…"]`）
    用户（和测试）很容易指到 `.result.nwk`，但年代只在 `.nexus` 里。
    这里按「原文件 → 补 `.date.nexus` → 同目录下同前缀的 `*date.nexus`」依次试，
    避免 MOTP 静默退回近似时间轴。
    """
    cands = [path]
    if path:
        base = path
        for suf in ('.result.nwk', '.result', '.nwk'):
            if base.endswith(suf):
                base = base[:-len(suf)]
                break
        d = os.path.dirname(path)
        cands += [path + '.date.nexus', base + '.result.date.nexus',
                  base + '.date.nexus']
        try:
            for f in os.listdir(d or '.'):
                if f.endswith('date.nexus') and f.startswith(
                        os.path.basename(base)):
                    cands.append(os.path.join(d, f))
        except OSError:
            pass
    seen = set()
    for c in cands:
        if not c or c in seen or not os.path.isfile(c):
            continue
        seen.add(c)
        root, dates = parse_lsd2_dated_tree(c)
        if root is not None and dates:
            return c, root, dates
    return None, None, {}


# ---------------------------------------------------------------- LSD2 分子钟定年

def find_lsd2():
    """定位 lsd2.exe：platform.json 的 tools.lsd2 优先，其次 3rd/tools/lsd2/。"""
    from Virus_Platform_Core.config import get_config, PLATFORM_ROOT
    try:
        p = get_config().tool('lsd2')
        if p and os.path.isfile(p):
            return p
    except Exception:
        pass
    for cand in ('3rd/tools/lsd2/lsd2.exe', '3rd/bin/lsd2.exe'):
        p = os.path.join(PLATFORM_ROOT, *cand.split('/'))
        if os.path.isfile(p):
            return p
    return None


def write_lsd2_dates(path, dates):
    """写 LSD2 的日期约束文件（`-d`）。返回写入的约束条数。

    格式（LSD2 v2.4.1 `-h` 实测）：第 1 行＝约束个数；其后每行 `名称 取值`，
    空格分隔，取值 = `固定值` / `l(下界)` / `u(上界)` / `b(下界,上界)`，
    名称可为单个 tip 或 `mrca(a,b,c)` —— 与 YR-MPE
    `plugins/components/methods/lsd2_logistics.py` 的实现逐字一致。

    ⚠️ **Windows 下必须 CRLF**：LSD2 win32 版期望 `\\r\\n`（YR-MPE 的插件专门按
    `platform.system()` 分支处理过这件事）。写成 `\\n` 会读成乱码。
    """
    items = [(k, str(v).strip()) for k, v in dates.items() if str(v).strip()]
    lines = [str(len(items))] + [f'{k} {v}' for k, v in items]
    with open(path, 'w', encoding='utf-8', newline='') as f:
        f.write('\r\n'.join(lines) + '\r\n')
    return len(items)


def lsd2_dating(tree_path, dates, seq_len, out_dir, exe=None, extra_args=None,
                prefix=None, progress=None, timeout=1800, rooting='auto'):
    """LSD2 最小二乘定年 → 速率 / tMRCA / 定年树（真正的时间标度）。

    tree_path  含全部叶名的树（本平台 NJ/FastTree 产物即可，叶名已净化）。
    dates      {树中叶名: 日期}；未匹配的叶按「无日期」处理。
    seq_len    比对长度（LSD2 的 `-s`，必需）——把枝长换算成每位点替换数。
    extra_args 追加参数：`-f 100`（bootstrap 给置信区间）、`-e 3`（Z>3 判离群剔除）、
               `-a <根日期>`、`-b <方差参数>`。
    rooting    **定根方法**（`-r`）：'auto'（默认）| 'a' | 'l' | 'k' | 'as' | None。

    ⚠️ **`-r` 不是钟模型，是定根方法**（LSD2 `-h` 原文）：
       `-r l` 局部重估根、`-r a` 在**全部枝上搜索**根位置（忽略原根）、
       `-r as` 约束模式全枝搜索、`-r k` 同枝重估。
       —— 而**严格钟 / 松弛钟**在 LSD2 里是 `-b varianceParameter`
       （越小越接近严格钟，越大越松弛）。曾把 `-r a` 当成"松弛钟"是**错的**。
    'auto' 的行为：先按原样跑；**满足任一条就自动加 `-r a`（全枝搜根）重跑**：
       ① LSD2 自报 `input trees are not rooted`（根节点三分叉的裸无根树
          —— VirPhyKit `h3n2_na_500.nwk` 实测走这条）；
       ② 结果里速率**顶到 `-t` 下界**（根节点恰是**两分叉**但根位乱放的树
          —— 本平台 NJ 产物实测走这条：LSD2 **不报错**，只把 rate 压成
          1e-10 再吐一个假 tMRCA=-3.1e7，见下方结果文件自报行）。
       本平台的树来自 NJ / FastTree（**无根**），所以这条回退几乎总会触发；
       之所以不默认直接加 `-r a`，是为了不悄悄丢掉用户**有意给定**的根。
       实际用了哪种会记在返回值的 `rooting_used` / `rooted_retry` 里。

    产物落 out_dir：`<prefix>.result`（主结果）、`<prefix>.result.nwk`（**定年树**）、
    `<prefix>.result.date.nexus`（**时间树**，节点带 `[&date=]`，MOTP/LTT 要用它）。
    返回 dict：rate / tmrca / objective / result_file / dated_tree / nexus /
    seq_len / n_constraints / cmd / rooting_used / rooted_retry / log_tail。

    ⚠️ 为什么用 `-d` 文件而不是把日期塞进 tip 名：实测（2026-09-16）不给 `-d` 时
    LSD2 输出 **tMRCA=0**（退回相对定年 T[root]=0/T[tips]=1），它**不解析** tip 名里的
    日期。这也顺带说明「平台建树会净化 tip 名」不构成障碍。
    ⚠️ LSD2 把产物写在**输入树所在目录** → 先把树拷进 out_dir 再跑，产物自然落在那儿。
    """
    exe = exe or find_lsd2()
    if not exe:
        raise RuntimeError('未找到 lsd2.exe（请登记 platform.json 的 tools.lsd2，'
                           '或放到 3rd/tools/lsd2/）')
    if not os.path.isfile(tree_path):
        raise RuntimeError(f'输入树不存在：{tree_path}')
    try:
        seq_len = int(float(seq_len))
    except (TypeError, ValueError):
        raise RuntimeError(f'seq_len 必须是数字（LSD2 -s）: {seq_len!r}')
    if seq_len <= 0:
        raise RuntimeError('seq_len 必须为正（LSD2 -s 用于把枝长换算成每位点替换数）')

    os.makedirs(out_dir, exist_ok=True)
    prefix = prefix or os.path.splitext(os.path.basename(tree_path))[0]
    local_tree = os.path.join(out_dir, f'{prefix}.nwk')
    # 输入树本来就在 out_dir（常见：调用方先在该目录建树）时不能 copyfile，
    # 同路径会抛 SameFileError。
    if os.path.abspath(tree_path) != os.path.abspath(local_tree):
        shutil.copyfile(tree_path, local_tree)
    date_file = os.path.join(out_dir, f'{prefix}.dates.txt')
    n_con = write_lsd2_dates(date_file, dates)
    if n_con == 0:
        raise RuntimeError('没有任何带日期的叶 —— LSD2 定年至少需要一条时间约束')

    base = [exe, '-i', local_tree, '-d', date_file, '-s', str(seq_len)]
    extra = [str(x) for x in (extra_args or [])]

    def _run(extra_args_):
        cmd_ = base + extra_args_
        if progress:
            progress('lsd2', 0.2, f'LSD2 定年（{n_con} 条时间约束）')
        try:
            # 按字节捕获再 decode_output（utf-8→GBK 兜底）：LSD2 把含中文
            # 路径的提示行按系统 ANSI 打印，text=True+utf-8 会在捕获时就
            # 把那些行打成 U+FFFD（不可逆），log_tail 进报告后无法挽回
            p = subprocess.run(cmd_, cwd=out_dir, capture_output=True,
                               timeout=timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f'LSD2 超时（>{timeout}s）')
        from Virus_Platform_Core.utils import decode_output
        out = decode_output((p.stdout or b'') + (p.stderr or b'')).strip()
        return cmd_, p.returncode, out

    used_rooting = (f'-r {rooting}' if rooting not in ('auto', None)
                    else ('-r (用户自定)' if any(a == '-r' for a in extra)
                          else None))

    def _result_path():
        rf = local_tree + '.result'
        if not os.path.isfile(rf):        # 兜底：按后缀扫目录
            cands = [os.path.join(out_dir, f) for f in os.listdir(out_dir)
                     if f.startswith(prefix) and f.endswith('.result')]
            if cands:
                rf = cands[0]
        return rf

    def _read_result(path):
        rate_ = tmrca_ = objective_ = None
        if os.path.isfile(path):
            for line in open(path, encoding='utf-8', errors='replace'):
                if 'tMRCA' in line:        # 结果行：rate …, tMRCA …, objective …
                    m = re.search(r'\brate\s+([-\d.eE+]+)', line)
                    if m:
                        rate_ = float(m.group(1))
                    m = re.search(r'tMRCA\s+([-\d.eE+]+)', line)
                    if m:
                        tmrca_ = float(m.group(1))
                    m = re.search(r'objective function\s+([-\d.eE+]+)', line)
                    if m:
                        objective_ = float(m.group(1))
                    break
        return rate_, tmrca_, objective_

    cmd, rc, log = _run(extra)
    rooted_retry = False
    retry_why = None
    if rc != 0 and rooting == 'auto' and 'not rooted' in log and '-r' not in extra:
        # 无根树：LSD2 要求 -g（外群）或 -r（搜索根位置）。这里自动补 -r a 重跑。
        retry_why = '输入树无根'
    elif rc == 0 and rooting == 'auto' and '-r' not in extra:
        # ⚠️ 静默垃圾兜底（2026-09-16 实测 LSD2 v2.4.1）：无根/根位乱放的树
        # 它**不报错**，而是把速率压到 `-t` 下界再返回一个可解析的假结果
        # （实测 rate=1e-10、tMRCA=-3.14e7，最早采样年才 1968；-r a 后
        # rate=1.4e-3、tMRCA=1960.8 才正常）。它自己在结果文件里写
        # "The estimated rate reaches the given lower bound." —— 拿这条自报
        # 当信号。不查会怎样：假值照样流进 LTT / 时间树（x=-3e7），无人察觉。
        try:
            with open(_result_path(), encoding='utf-8', errors='replace') as _f:
                if 'reaches the given lower bound' in _f.read():
                    retry_why = '速率顶到下界（多半是根位不对）'
        except OSError:
            pass
    if retry_why:
        extra2 = extra + ['-r', 'a']
        cmd, rc, log = _run(extra2)
        used_rooting, rooted_retry = f'-r a（自动补：{retry_why}）', True
    if rc != 0:
        raise RuntimeError(f'LSD2 退出码 {rc}：'
                           + ' / '.join(log.splitlines()[-3:]))

    result_file = _result_path()
    rate, tmrca, objective = _read_result(result_file)
    dated_tree = result_file + '.nwk'
    nexus = result_file + '.date.nexus'
    # LSD2 按系统 ANSI 写自产结果文件（.result 里回显输入路径，中文路径
    # 即 GBK）：跑完就地转 UTF-8，与平台其余产物编码统一
    from Virus_Platform_Core.utils import to_utf8_file
    for _p in (result_file, dated_tree, nexus):
        to_utf8_file(_p)
    if progress:
        progress('lsd2', 1.0, f'定年完成：rate={rate} tMRCA={tmrca}')
    return {
        'engine': 'lsd2',
        'rate': rate,
        'tmrca': tmrca,
        'objective': objective,
        'result_file': result_file if os.path.isfile(result_file) else None,
        'dated_tree': dated_tree if os.path.isfile(dated_tree) else None,
        'nexus': nexus if os.path.isfile(nexus) else None,
        'seq_len': seq_len,
        'n_constraints': n_con,
        'cmd': ' '.join(cmd),
        'rooting_used': used_rooting,
        'rooted_retry': rooted_retry,
        'log_tail': '\n'.join(log.splitlines()[-6:]),
    }


# ---------------------------------------------------------------- 年份解析

_DAYS_BEFORE = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)


def resolve_years(leaf_names, dates, date_trait='year'):
    """把叶名解析成采样年代：{叶名: int 年 或 float 十进制年}（解析不出则不入表）。

    规则（与 root_to_tip 原先的内联实现**同口径**，抽出是为了让 RTT 与
    LSD2 定年共用一套，避免两处各写一套而悄悄漂移）：
      对每个叶名，优先**精确匹配**键名；没有精确键时，再按 dates 的**插入序**
      取第一个前缀/包含匹配的键（精确优先是 2026-09-16 修掉的静默错配：
      `t10` 会被更短的 `t1` 前缀抢先命中，拿到别人的年份且毫无提示）。
      取该键的值（dict 取 date_trait 列，否则本身就是日期串），
      再用开头 4 位数字解析年份；取不到就放弃该叶
      （注意：命中但不含年份也会放弃，与旧实现一致）。

    ⚠️ 2026-09-16 提升（十进制年）：日期带月/日时，旧实现把 `2020-07-01` 直接
    截成 `2020` —— 月/日级采样（同一年的不同时点）全部压成同一 x 值，RTT
    的回归与 LSD2 的 tip 时间约束等于丢掉全部年内信息（这是 zika_Vietnam
    管道 `date2years` 十进制年口径要解决的问题）。现在：
      · 只有 4 位年份 → 返回 **int**（与旧行为逐位一致，老数据零变化）；
      · `YYYY-MM[-DD]` / `YYYY/MM[/DD]` → 返回 **float 十进制年**
        = 年 + (年内第几天 − 1)/365.25（**1 月 1 日 = 年首整点**，与 BEAST
        日期换算口径一致；只有月份时取当月 15 日）。
    ⚠️ 只认 `-` 与 `/` 两种分隔符：`2008.58197` 这类**本身就是十进制年**的值
    （本项目 tip 名里实测出现过）用 `.` 分隔，若按"点号日期"解读会把 .58 当
    月号 → 被钳到 12 月、年份漂移近一年。`.` 形式一律维持旧口径（取前 4 位）。
    LSD2 的 `-d` 约束文件接受十进制年（浮点时间轴），实测验证见
    tests/_check_dyn_fixes.py 的 H 段。
    """
    out = {}
    items = list(dates.items())
    for nm in leaf_names:
        if not nm:
            continue
        hit = next(((k, v) for k, v in items if nm == k), None)
        if hit is None:
            hit = next(((k, v) for k, v in items
                        if nm.startswith(k) or k in nm), None)
        if hit is None:
            continue
        raw = hit[1].get(date_trait) if isinstance(hit[1], dict) else hit[1]
        if raw:
            m = re.match(
                r'(\d{4})(?:[-/](\d{1,2})(?:[-/](\d{1,2}))?)?',
                str(raw).strip())
            if m:
                y = int(m.group(1))
                if not m.group(2):
                    out[nm] = y
                else:
                    mo = min(12, max(1, int(m.group(2))))
                    d = min(31, max(1, int(m.group(3) or 15)))
                    doy = _DAYS_BEFORE[mo - 1] + d
                    out[nm] = round(y + (doy - 1) / 365.25, 5)
    return out


# ---------------------------------------------------------------- 主入口

def analyze(aln_path, tree_path, meta_path=None, trait='region',
            method='nj', out_dir=None, progress=None,
            rrt_perm=0, seed=0, rssp_bs=0,
            dated_tree=None, motp_bin=0.0, date_trait='year',
            import_export=True, coords_path=None,
            beast_log=None, beast_trees=None, beast_burnin=0.1,
            beast_min_prob=0.0, beast_max_trees=2000):
    """主入口：比对 + 树 + 元数据 → 迁移重构产物。

    meta 缺省时从 FASTA 头解析（`>acc|区域|年份`）。
    返回 summary dict；产物（标注树 / 迁移矩阵 / 叶表）写入 out_dir。

    可选（默认关，避免拖慢常规运行；import_export 例外——默认开、纯计数）：
      rrt_perm   区域随机化检验的置换次数（0=不做；常用 1000）
      rssp_bs    根状态不确定性的 bootstrap 复本数（0=不做；常用 100）
      seed       两者共用的随机种子（结果可复现）
      dated_tree LSD2 的 `.date.nexus`（有则 MOTP 用**真正的节点年代**）
      motp_bin   迁移随时间的分箱宽度（年）；>0 才做 MOTP
      date_trait 日期列名（无定年树时 MOTP 用它做近似时间轴）
      import_export 引入 vs 本地传播拆分（默认开；纯计数、毫秒级，不拖慢运行，
                 结果同时是返回值 import_export 与产物 import_export.tsv）
                 迁移事件分类（transition_classes：Import/Direct/Indirect/Distant/
                 Unresolved，四通道证据）跟随同一开关，产物 transition_classes.tsv
      coords_path 坐标表（可选）：逐样本点层（元数据 lat/lon 列优先，其次本表
                 逐行条目）+ 区划落点（样本坐标中位数 > 本表质心 > 内置质心表）；
                 命中/坏行/缺坐标全部进返回值，不静默
      beast_log  BEAST `.log`（服务器端跑完的产物；有则导入逐列诊断 + 跳转计数）
      beast_trees BEAST `.trees` / MCC 树（≥2 棵＝后验样本；1 棵＝MCC 单树）
      beast_burnin 导 BEAST 产物时丢的开头比例（日志行 / 树；默认 0.1）
      beast_min_prob MCC/样本树里只统计"父子两端后验都 ≥ 该阈值"的枝（0=不过滤）
      beast_max_trees `.trees` 最多用多少棵（超出如实上报，不静默）
    """
    from Virus_Platform_Core.phylo import _run_nj, _run_fasttree
    from Virus_Platform_Core.utils import load_alignment

    def prog(f, m):
        if progress:
            progress('phylogeo', f, m)

    aln = load_alignment(aln_path)
    names = list(aln.keys())
    if len(names) < 3:
        raise RuntimeError('系统地理分析需要 ≥3 条序列')
    meta = load_metadata(meta_path) if (meta_path and os.path.isfile(meta_path)) else {}

    prog(0.15, f'建树（{method}）')
    if method == 'fasttree':
        tree_path_used = _run_fasttree(aln_path, tree_path)
    else:
        tree_path_used = _run_nj(aln_path, tree_path)
    text = open(tree_path_used, encoding='utf-8', errors='replace').read()
    root = parse_newick(text)

    # 叶名对齐：Newick 里的名字做过安全替换（phylo._newick_safe），
    # 按净化名 ↔ 原 FASTA 头第一词回映射（树拓扑不保证保持输入顺序）
    leaves = []

    def _leaves(n):
        if not n['children']:
            leaves.append(n)
        for c in n['children']:
            _leaves(c)

    _leaves(root)
    leaf_names = [x['name'] for x in leaves]

    # 元数据匹配：优先 CSV 精确/前缀匹配，否则 FASTA 头（读原比对头）
    prog(0.45, '匹配区划元数据')
    headers = _fasta_headers(aln_path)

    def _safe_first(h):
        tok = h.split()[0] if h.split() else h
        return re.sub(r'[(),:;\[\]\'"]+', '_', tok) or 'seq'

    safe2header = {}
    for h in headers:
        safe2header.setdefault(safe_first_token(h), h)
    # 匹配逻辑与「ID 对齐报告」共用同一实现（match_states_and_years），
    # 否则报告说的命中情况会和这里实际用的不一致。
    states, leaf_year, unresolved, _rows = match_states_and_years(
        leaf_names, safe2header, meta, trait=trait, date_trait=date_trait)

    # ---- 数据质量告警（只报事实，不改行为）----
    # 踩过的坑：元数据里没有区划列时，_header_traits 会去读 FASTA 头的第二个字段
    # （`acc|X|年份` 约定取 X）。若头是 `acc|<整条描述>`（例如平台示例集
    # `NC_001440.1|Cucumber_mosaic_virus_RNA3_RefSeq`），X 就是整条描述 →
    # **每条序列都成"独立区划"**，迁移矩阵与区域随机化检验随之全部失去意义
    # （RRT 会老实地报 p≈1.0，但看不出是数据问题）。这里显式点出来。
    warnings = []
    known = [v for v in states.values() if v and v != 'Unknown']
    n_unk = len(states) - len(known)
    if len(known) >= 5 and len(set(known)) == len(known):
        warnings.append(
            f'{len(known)} 条样本有 {len(set(known))} 个各不相同的区划 —— 多为区划列没匹配上、'
            '回退成了整条 FASTA 描述。请检查元数据 CSV，或把 FASTA 头写成 `名称|区划|年份`；'
            '否则迁移矩阵与区域随机化检验都无意义。')
    if n_unk:
        warnings.append(f'{n_unk} 条样本没有区划 —— 按缺失数据处理（Fitch 中不产生迁移），'
                        '不影响拓扑但会削弱区划信号的检出力。')
    # ---- 方法学警告（条件触发；结论出处：drifting_into_nowhere 仿真研究 ——
    # 论文 "Can Bayesian phylogeography reconstruct migrations and expansions in
    # human history?"（Neureiter et al.）的全部仿真/评估代码，27 仓解读报告的
    # E1 段有摘要（报告在姊妹平台：植物病毒分析平台/docs/git-repo全量解读_20260916.md，
    # 本仓 docs/ 下没有副本）。结论：方向性扩散（跨区扩张）越强、且没有历史/古样本时，
    # 根位置与迁移方向的重建误差越大，HPD 覆盖率也会偏离名义水平。这里只做
    # 「数据条件检查」，不替用户下结论、也不编造具体数值。）----
    if len(set(known)) < 2:
        warnings.append(
            f'可用区划只有 {len(set(known))} 个 —— 区间迁移矩阵与区域随机化检验'
            '都失去意义（至少要 2 个区划才谈得上「区间迁移」）。')
    if not leaf_year:
        warnings.append(
            '样本里没有任何采样年份 —— 无历史/古样本时，根位置与迁移方向的重建'
            '误差更大；若数据是方向性扩散（跨区扩张），根位与方向尤其不可轻信。')

    prog(0.7, 'Fitch 迁移重构')
    node_states, transitions, n_changes = fitch_mugration(root, states)

    # 区划计数 + 迁移矩阵
    region_counts = Counter(states.values())
    matrix = Counter()
    for _p, _c, frm, to in transitions:
        matrix[(frm, to)] += 1
    regions = sorted(region_counts)

    # ---- 核心产物**先落盘** ----
    # 顺序很重要：RRT/RSPP 是可选统计，万一它失败/超时/被环境卡住，
    # 迁移矩阵与标注树也必须已经在磁盘上（2026-09-16 实测踩过：把统计放前面，
    # 统计一卡就连一个核心产物都没写出来）。
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, 'tree_annotated.nwk'), 'w',
                  encoding='utf-8', newline='\n') as f:
            f.write(write_annotated(root, states) + '\n')
        with open(os.path.join(out_dir, 'migration_matrix.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('from\tto\tcount\n')
            for (frm, to), n in sorted(matrix.items(), key=lambda x: -x[1]):
                f.write(f'{frm}\t{to}\t{n}\n')
        with open(os.path.join(out_dir, 'leaf_states.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('sample\t' + trait + '\n')
            for n in leaf_names:
                f.write(f'{n}\t{states[n]}\n')
        with open(os.path.join(out_dir, 'transitions.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('parent\tchild\tfrom\tto\n')
            for p, c, frm, to in transitions:
                f.write(f'{p}\t{c}\t{frm}\t{to}\n')

    # ---- 逐样本经纬度 + 区划落点（样本点层 / 三分类距离的底料）----
    # 样本点来源：元数据 lat/lon 列 → 坐标表逐行条目（逐样本坐标表）；
    # 区划落点 = 该区划样本坐标中位数 > 坐标表质心 > 内置质心表（src 如实标出处）。
    prog(0.75, '逐样本经纬度 / 区划落点')
    coord_rows, coord_skips = [], []
    user_coords = {}
    if coords_path and os.path.isfile(coords_path):
        coord_rows, coord_skips = parse_coord_rows(coords_path)
        user_coords = coord_table_centroid(coord_rows)
    sc = match_sample_coords(leaf_names, safe2header, meta, coord_rows)
    sample_points = sc['points']
    coords, coord_src, coord_missing, n_pts_used = resolve_region_positions(
        regions, states, sample_points, user_coords)
    if sc['n_bad']:
        warnings.append(
            f"元数据里有 {sc['n_bad']} 条样本的经纬度不是有效数值/越界 —— 已跳过，"
            '不参与样本点层与区划落点计算。')
    if out_dir and coords:
        with open(os.path.join(out_dir, 'coords_resolved.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('region\tlatitude\tlongitude\tsource\n')
            for r_, c_ in sorted(coords.items()):
                f.write(f"{r_}\t{c_[0]}\t{c_[1]}\t{coord_src.get(r_, '')}\n")
    if out_dir and sample_points:
        with open(os.path.join(out_dir, 'sample_coords.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('sample\tlatitude\tlongitude\tregion\tsource\n')
            for n_ in sorted(sample_points):
                pt = sample_points[n_]
                f.write(f"{n_}\t{pt[0]}\t{pt[1]}\t{states.get(n_, '')}"
                        f"\t{sc['src'].get(n_, '')}\n")

    # ---- 引入 vs 本地传播（import/export；默认开，纯计数）----
    # 复用上面 fitch_mugration 的 node_states（同一棵树、同一批 id），不重跑重构。
    ie = None
    if import_export:
        prog(0.76, '引入 / 本地传播拆分')
        ie = import_export_balance(root, states, node_states=node_states)
        if out_dir:
            with open(os.path.join(out_dir, 'import_export.tsv'), 'w',
                      encoding='utf-8', newline='') as f:
                f.write('region\timportations\tlocal_transmissions\t'
                        'exportations\tnet\tn_tips\tn_import_groups\n')
                for r_, v_ in ie['per_region'].items():
                    f.write(f"{r_}\t{v_['importations']}\t"
                            f"{v_['local_transmissions']}\t"
                            f"{v_['exportations']}\t{v_['net']}\t"
                            f"{v_['n_tips']}\t{v_['n_import_groups']}\n")

    # ---- 迁移事件分类（Direct/Indirect/Distant/Import；四通道证据）----
    # 与 import_export 同量级成本（同一 node_states、同一引入判定），默认开；
    # 距离带依赖上面的 coords/sample_points，没有坐标时除 Import 外如实记
    # Unresolved（不硬塞档），并发一条告警。
    prog(0.77, '迁移事件分类')
    cls = classify_transitions(root, states, node_states, coords,
                               sample_points, leaf_year,
                               max((len(v) for v in aln.values()), default=0))
    cls_corr = corridor_classes(cls['rows'])
    if cls['n_edges'] and not any(r['geo_km'] is not None for r in cls['rows']):
        warnings.append(
            f"{cls['n_edges']} 条跨区边全部没有两端坐标 —— 迁移事件分类除 Import 外"
            '都是 Unresolved；要出距离分带请提供坐标表（coords）或元数据 lat/lon 列。')
    if out_dir:
        with open(os.path.join(out_dir, 'transition_classes.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('parent\tchild\tfrom\tto\tlabel\tgeo_km\tstep_km\t'
                    'edge_snp\tdt_years\tn_p\tn_c\tgroup_size\tintro\tnote\n')

            def _s(v):
                return '' if v is None else v

            for r_ in cls['rows']:
                f.write(f"{r_['parent']}\t{r_['child']}\t{r_['from']}\t"
                        f"{r_['to']}\t{r_['label']}\t{_s(r_['geo_km'])}\t"
                        f"{_s(r_['step_km'])}\t{r_['edge_snp']}\t"
                        f"{_s(r_['dt_years'])}\t{r_['n_p']}\t{r_['n_c']}\t"
                        f"{r_['group_size']}\t"
                        f"{1 if r_['intro'] else 0}\t{r_['note']}\n")

    # ---- 定年树基准（LSD2 `.date.nexus`）：时间树与 MOTP 共用 ----
    # 只解析一次：droot/ndates 是**同一棵树对象**上的年代字典（id() 对齐），
    # 谁要用都得在这棵树上重算，不能拿未定年树的节点 id 去查。
    dpath, droot, ndates = (None, None, {})
    if dated_tree:
        # 用户可能指到 `.result.nwk`（替换率尺度、无年代），
        # resolve_dated_tree 会去同目录找 `.date.nexus`
        dpath, droot, ndates = resolve_dated_tree(dated_tree)
        if droot is None:
            warnings.append(
                f'定年树里没有节点年代（{os.path.basename(dated_tree)}）——'
                'LSD2 的年代写在 `.date.nexus` 里，`.result.nwk` 是替换率尺度；'
                '时间树与 MOTP 已退回「后代叶采样年」近似时间轴。')

    # ---- 时间树（T3）：x ＝ 定年结果的节点年代 ----
    # 没有定年树就**不给**这张图：把「后代叶采样年中位数」当节点年代画出来会得到
    # 一棵与定年无关的假时间树（近似轴只给 MOTP 用，并有 axis 字段明示）。
    timetree = None
    if droot is not None:
        prog(0.78, '时间树（定年结果）')
        timetree = time_tree_segments(
            droot, ndates, axis_label='采样年（日历年）',
            source=os.path.basename(dpath or dated_tree or ''))
        if out_dir and timetree.get('x'):
            with open(os.path.join(out_dir, 'timetree_tips.tsv'), 'w',
                      encoding='utf-8', newline='') as f:
                f.write('sample\tyear\n')
                for t in timetree['tips']:
                    f.write(f"{t['name']}\t{t['x']}\n")
            with open(os.path.join(out_dir, 'timetree_nodes.tsv'), 'w',
                      encoding='utf-8', newline='') as f:
                f.write('node\tdate\n')
                for n in _iter_pre(droot):
                    v = ndates.get(id(n))
                    if v is not None:
                        f.write(f"{n['name'] or 'Node'}\t{round(v, 4)}\n")
        if timetree.get('note'):
            prog(0.79, timetree['note'])

    # ---- 可选项③：迁移随时间（MOTP）----
    motp = None
    if float(motp_bin or 0) > 0:
        prog(0.80, '迁移随时间（MOTP）分箱')
        if droot is not None:
            # 在**定年树本身**上重跑 Fitch —— 不能拿 root 的节点 id 去查定年树的
            # 年代：两棵树是不同对象、id 无关（哪怕拓扑同源）。叶名一致，所以在
            # 定年树上重算一遍转移，年代天然对齐。
            dnames = [x['name'] for x in _iter_pre(droot) if not x['children']]
            dstates = {nm: states.get(nm, 'Unknown') for nm in dnames}
            motp = migration_over_time(droot, dstates, node_dates=ndates,
                                       bin_width=float(motp_bin))
            motp['source'] = os.path.basename(dpath or dated_tree)
        else:
            motp = migration_over_time(root, states, years=leaf_year or None,
                                       bin_width=float(motp_bin))
        if out_dir and motp.get('rows'):
            with open(os.path.join(out_dir, 'motp.tsv'), 'w',
                      encoding='utf-8', newline='') as f:
                f.write('bin_start\tbin_end\tfrom\tto\tcount\n')
                for r in motp['rows']:
                    f.write(f"{r['lo']}\t{r['hi']}\t{r['from']}\t"
                            f"{r['to']}\t{r['count']}\n")

    # ---- 可选项①：区域随机化检验（RRT）----
    # 顺序无妨：_fitch 是纯函数（状态只写局部字典），置换几百次不会污染节点，
    # 上面写出的标注树仍是观测标签的结果。
    perm = None
    null_counts = []
    if int(rrt_perm or 0) > 0:
        prog(0.82, f'区域随机化检验（置换 {int(rrt_perm)} 次）')
        perm = region_permutation_test(root, states, n_perm=int(rrt_perm),
                                       seed=int(seed))
        null_counts = perm.pop('null_counts', [])
        # 跑了检验但没判出显著结构 → 条件提示（不显著时不加这句反而会误导）。
        # 出处同上：方向性扩散 + 无历史样本时，根位/方向推断的风险更高。
        if perm.get('p_low') is not None and not (
                perm['p_low'] < 0.05 and perm['p_low'] <= perm['p_high']):
            warnings.append(
                f"区域随机化检验未判出显著地理结构（p_low={perm['p_low']}）——"
                '区划与树的关系不能排除随机，迁移方向与根位置（进而 import/export '
                '的方向）只能当**探索性**结果看。')

    # ---- 可选项②：RSPP bootstrap 不确定性 ----
    rssp = None
    if int(rssp_bs or 0) > 0:
        prog(0.88, f'RSPP bootstrap（{int(rssp_bs)} 个复本）')
        rssp = rssp_bootstrap(aln, states, out_dir or os.path.dirname(
            os.path.abspath(tree_path)), method=method, n_boot=int(rssp_bs),
            seed=int(seed), progress=progress, observed=n_changes)

    # ---- BEAST 产物导入（.log / .trees / MCC；服务器端跑、本地读）----
    # 本地不跑 BEAST（见 beast_io 模块头）：只导入服务器端产物。三口径
    # （tree_sample / jump / mcc）与 RSPP 一样**只出一个 driver**
    # 给弧线分带 —— 同一张图上一条弧线只能有一个线宽口径。
    beast = None
    beast_driver = None
    beast_corr = []
    if beast_log or beast_trees:
        prog(0.92, '导入 BEAST 产物（.log / .trees）')
        try:
            beast, beast_driver, beast_corr = import_beast(
                beast_log=beast_log, beast_trees=beast_trees,
                burnin=float(beast_burnin or 0.1),
                min_prob=float(beast_min_prob or 0.0),
                max_trees=int(beast_max_trees or 0), out_dir=out_dir,
                meta_states=sorted({s for s in states.values() if s}),
                progress=prog)
        except Exception as e:                      # noqa: BLE001 — 透出原文
            # 解析失败不该毁掉主分析：树已建好、矩阵/统计都在，BEAST 面板
            # 给红灯原文，前端显式说明（错误原文透出，不吞）。
            beast = {'log': None, 'jump': None, 'jump_error': None,
                     'trees': None, 'driver': None, 'driver_label': None,
                     'n_units': 0, 'state_check': None, 'artifacts': [],
                     'warnings': [], 'n_corr': 0,
                     'log_file': (os.path.basename(beast_log) if beast_log
                                  else None),
                     'trees_file': (os.path.basename(beast_trees) if beast_trees
                                    else None),
                     'error': f'{type(e).__name__}: {e}'}
            warnings.append('BEAST 产物导入失败（其余分析结果不受影响）：'
                            + beast['error'])
        else:
            warnings += beast['warnings']

    # ---- 后验权重分带（驱动源优先序：BEAST 产物 → RSPP）----
    # 真实 BEAST 产物（用户为这次分析跑的外部模型）优先于 RSPP 支持率
    # （复本频率，非后验口径）。
    # 都没有 → None：**不降级成计数口径冒充"权重"**，前端据此显示「未评估」
    # （与弧线虚线的 M5 处理同一纪律）。多个都跑了也只出一份分带：同一张图上
    # 一条弧线只能有一个线宽口径，driver 混着用会让「粗」在图上不再可比。
    # BEAST 状态名与元数据区划对不上时 driver=None（见 import_beast 的核对）
    # —— 覆盖范围对不上的分带比没有分带更误导。
    weight_bands = None
    if beast_driver:
        weight_bands = beast_weight_bands(beast_driver, beast_corr,
                                         dict(matrix),
                                         n_units=(beast or {}).get('n_units'))
    elif rssp and rssp.get('n_ok'):
        weight_bands = posterior_weight_bands(rssp['corridor_support'],
                                              dict(matrix))
    if weight_bands and out_dir:
        with open(os.path.join(out_dir, 'weight_bands.json'), 'w',
                  encoding='utf-8', newline='\n') as f:
            json.dump(weight_bands, f, ensure_ascii=False, indent=1)

    if out_dir and null_counts:
        with open(os.path.join(out_dir, 'rrt_permutation.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('rep\tn_transitions\n')
            for i, v in enumerate(null_counts):
                f.write(f'{i + 1}\t{v}\n')

    return {
        'n_seqs': len(names),
        'trait': trait,
        'method': method,
        'regions': dict(region_counts),
        'n_transitions': n_changes,
        'transitions_top': [{'from': a, 'to': b, 'count': n}
                            for (a, b), n in matrix.most_common(12)],
        'unresolved': unresolved,
        'matrix': [{'from': a, 'to': b, 'count': n}
                   for (a, b), n in sorted(matrix.items(), key=lambda x: -x[1])],
        'tree': tree_path_used,
        # 带采样年份的叶数（前端据此提示「树无时间信息」的方法学风险）
        'n_years': len(leaf_year),
        # 数据质量告警（前端直接展示；空列表＝无告警）
        'warnings': warnings,
        # 引入 vs 本地传播（默认开；未开启时为 None，前端据此隐藏面板）
        'import_export': ie,
        # 迁移事件分类（Direct/Indirect/Distant/Import；默认开，纯计数）
        # 同类走廊聚合（地图弧线按主导类上色用）
        'transition_classes': cls,
        'corridor_classes': cls_corr,
        # 可选结果（未开启时为 None，前端据此隐藏对应面板）
        'motp': motp,
        'rrt': perm,
        'rssp': rssp,
        # 后验权重分带（驱动源：BEAST 产物 > RSPP；都没有为 None）
        'weight_bands': weight_bands,
        # BEAST 产物导入（.log / .trees / MCC；未导入时为 None。全量落
        # beast_log.json / beast_trees.json，summary 是精简版）
        'beast': beast,
        # 时间树（x＝定年结果的节点年代；无定年树时为 None）
        'timetree': timetree,
        # 逐样本经纬度 + 区划落点（样本点层；弧线地图/三分类距离的底料）
        'coords': coords, 'coord_src': coord_src,
        'coord_missing': coord_missing, 'coord_skips': coord_skips,
        'coord_table_regions': len(user_coords),
        'coord_table_rows': len(coord_rows),
        'sample_points': [{'name': n_, 'lat': sample_points[n_][0],
                           'lon': sample_points[n_][1],
                           'region': states.get(n_, ''),
                           'src': sc['src'].get(n_, '')}
                          for n_ in sorted(sample_points)],
        'n_sample_coords': len(sample_points),
        'n_sample_coords_used': n_pts_used,
    }




def root_to_tip(tree_path, dates, date_trait='year'):
    """根到尾回归（TreeTime RTT 的免依赖实现，含最优重根搜索）。

    NJ/FastTree 输出的 Newick 根位置是任意的——错误根位会摧毁时间信号
    （RTT 的经典陷阱）。这里把树视为无根图，逐枝尝试中点重根，
    取 R² 最高的根位（TreeTime 的 reroot 策略同思路）。

    dates: {样本名或原头首词: 日期串/含日期的 dict}（取 4 位年份）。
    返回 dict：slope（替换/位点/年）、intercept、r2、n、points、root_branch。
    """
    text = open(tree_path, encoding='utf-8', errors='replace').read()
    root = parse_newick(text)

    # 有根树 → 无根邻接图 {id: [(id2, length)]}，叶名表
    adj = {}
    names = {}

    def _reg(n):
        i = id(n)
        names[i] = n['name'] or ''
        adj.setdefault(i, [])
        return i

    def _walk(n):
        i = _reg(n)
        for c in n['children']:
            ci = _walk(c)
            ln = c['length'] or 0.0
            adj[i].append((ci, ln))
            adj[ci].append((i, ln))
        return i

    _walk(root)      # 只为建 adj/names（返回值 _walk 里即 _reg(root) 的 id，
                     #  从未被使用 → F841，2026-10-03 删去赋值）

    # 删掉辅助键，收集叶（度=1）与其名字
    leaves = [i for i, nb in adj.items() if len(nb) == 1]
    name_of = {i: names.get(i, '') for i in adj}
    # 年份**预解析一次**：原实现把查找写在下面的重根循环里，于是每换一个候选根位
    # 都要对每片叶再遍历一遍 dates（O(枝数 × 叶数 × |dates|)）。规则不变，
    # 只是挪到循环外 —— 顺便让 RTT 与 LSD2 定年共用同一套口径（resolve_years）。
    year_by_leaf = resolve_years(list(name_of.values()), dates, date_trait)

    def _tip_dists(a, block):
        """从 a 出发（不走 block 方向）的各叶距离。"""
        out = {}
        stack = [(a, 0.0, block)]
        while stack:
            cur, d, frm = stack.pop()
            nb = [(x, l) for x, l in adj[cur] if x != frm]
            kids = [x for x, _ in nb]
            if not kids:
                out[cur] = d
                continue
            for x, l in nb:
                stack.append((x, d + l, cur))
        return out

    best = None
    seen_branch = set()
    for a in list(adj):
        for b, ln in adj[a]:
            key = (min(a, b), max(a, b))
            if key in seen_branch:
                continue
            seen_branch.add(key)
            half = ln / 2.0
            # 以该枝中点为根：a 侧叶距 = d(t,a)+half；b 侧 = d(t,b)+half
            da = _tip_dists(a, b)
            db = _tip_dists(b, a)
            pts = []
            for leaf in leaves:
                nm = name_of.get(leaf, '')
                ds = year_by_leaf.get(nm)
                if ds is None:
                    continue
                d0 = da.get(leaf)
                if d0 is None:
                    d0 = db.get(leaf)
                    if d0 is None:
                        continue
                    pts.append((nm, d0 + half, float(ds)))
                else:
                    pts.append((nm, d0 + half, float(ds)))
            if len(pts) < 5:
                continue
            n = len(pts)
            mx = sum(p[2] for p in pts) / n
            my = sum(p[1] for p in pts) / n
            sxx = sum((p[2] - mx) ** 2 for p in pts)
            sxy = sum((p[2] - mx) * (p[1] - my) for p in pts)
            syy = sum((p[1] - my) ** 2 for p in pts)
            if sxx == 0 or syy == 0:
                continue
            slope = sxy / sxx
            r2 = (sxy ** 2) / (sxx * syy)
            if best is None or r2 > best['r2']:
                best = {'slope': slope, 'intercept': my - slope * mx,
                        'r2': round(r2, 4), 'n': n, 'root_branch': key,
                        'points': [{'name': nm, 'dist': round(d, 6),
                                    'year': int(yr)} for nm, d, yr in pts]}
    if best is None:
        raise RuntimeError('有年份的样本不足（≥5）或采样年份无变化：'
                           '无法做根到尾回归')
    best['r2'] = round(best['r2'], 4)
    return best


def _fasta_headers(path):
    heads = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                heads.append(line[1:].strip())
    return heads
