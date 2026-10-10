# -*- coding: utf-8 -*-
"""默认鉴定库覆盖审计：主要植物病毒/类病毒清单 vs kv_index 参考库。

清单覆盖：粮油/蔬菜/果树/葡萄/经济作物的主要病毒 + EPPO A1/A2 关注种 +
主要类病毒。命中判定：Species_NCBI 列包含关键词（大小写不敏感）。
"""
import csv
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
# 平台根：默认自动推导（tests/bench/ 上溯 3 层 = 平台根），可用 VP_ROOT 覆盖
PLATFORM = (os.environ.get('VP_ROOT')
            or os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))))
KV_INFO = os.path.join(PLATFORM, 'databases', 'virusref_db', 'kv_index',
                       'reference.ref_info.tsv')

# (组, 显示名, 库内 Species 关键词列表 —— 任一命中即算覆盖)
CHECKLIST = [
    ('粮食/薯类', '番茄斑萎病毒 TSWV', ['tomato spotted wilt']),
    ('粮食/薯类', '马铃薯Y病毒 PVY', ['potato virus y']),
    ('粮食/薯类', '马铃薯X病毒 PVX', ['potato virus x']),
    ('粮食/薯类', '马铃薯S病毒 PVS', ['potato virus s']),
    ('粮食/薯类', '马铃薯A病毒 PVA', ['potato virus a']),
    ('粮食/薯类', '马铃薯卷叶病毒 PLRV', ['potato leafroll']),
    ('粮食/薯类', '烟草花叶病毒 TMV', ['tobacco mosaic virus']),
    ('粮食/薯类', '番茄花叶病毒 ToMV', ['tomato mosaic virus']),
    ('粮食/薯类', '黄瓜花叶病毒 CMV', ['cucumber mosaic virus']),
    ('粮食/薯类', '番茄褪绿病毒 ToCV(传毒介体)', ['tomato chlorosis']),
    ('粮食/薯类', '番茄 souhaitez 斑萎病毒属其它种', ['orthotospovirus']),
    ('蔬菜', '番茄黄化曲叶病毒 TYLCV', ['yellow leaf curl']),
    ('蔬菜', '番茄褐色皱果病毒 ToBRFV', ['tomato brown rugose fruit', 'tobamovirus']),
    ('蔬菜', '辣椒轻斑驳病毒 PMMoV', ['mild mottle']),
    ('蔬菜', '番茄顶曲病毒 ToTV', ['tomato torrado']),
    ('蔬菜', '番茄黄化矮缩病毒 TYDV', ['yellow dwarf']),
    ('蔬菜', '蚕豆萎蔫病毒 BBWV', ['broad bean wilt']),
    ('蔬菜', '莴苣传染性黄化 LIYV', ['lettuce infectious yellows']),
    ('蔬菜', '甜菜曲顶病毒 BCTV', ['beet curly top']),
    ('葫芦科', '小西葫芦黄化花叶 ZYMV', ['zucchini yellow mosaic']),
    ('葫芦科', '西瓜花叶 WMV', ['watermelon mosaic']),
    ('葫芦科', '南瓜花叶 SqMV', ['squash mosaic']),
    ('葫芦科', '甜瓜坏死斑点 MNSV', ['melon necrotic spot']),
    ('葫芦科', '番木瓜环斑 PRSV', ['papaya ringspot']),
    ('果树/柑橘', '柑橘衰退病毒 CTV', ['citrus tristeza']),
    ('果树/果树', '李痘病毒 PPV', ['plum pox']),
    ('果树/果树', '李矮化 PDV', ['prune dwarf']),
    ('果树/果树', '李属坏死环斑 PNRSV', ['necrotic ringspot']),
    ('果树/果树', '苹果褪绿叶斑 ACLSV', ['chlorotic leaf spot']),
    ('果树/果树', '苹果茎沟 ASGV', ['apple stem grooving']),
    ('果树/果树', '苹果茎痘 ASPV', ['apple stem pitting']),
    ('果树/果树', '苹果花叶 ApMV', ['apple mosaic']),
    ('果树/葡萄', '葡萄扇叶 GFLV', ['fanleaf']),
    ('果树/葡萄', '阿拉伯花叶 ArMV', ['arabis mosaic']),
    ('果树/葡萄', '葡萄卷叶相关 GLRaV-1', ['grapevine leafroll-associated virus 1']),
    ('果树/葡萄', '葡萄卷叶相关 GLRaV-2', ['grapevine leafroll-associated virus 2']),
    ('果树/葡萄', '葡萄卷叶相关 GLRaV-3', ['grapevine leafroll-associated virus 3']),
    ('果树/葡萄', '葡萄卷叶相关 GLRaV-4', ['grapevine leafroll-associated virus 4']),
    ('果树/葡萄', '葡萄病毒A GVA', ['grapevine virus a']),
    ('果树/葡萄', '葡萄病毒B GVB', ['grapevine virus b']),
    ('果树/葡萄', '葡萄病毒F GVF', ['grapevine virus f']),
    ('果树/葡萄', '葡萄钉状皮 GRSPaV', ['rupestris stem pitting']),
    ('果树/葡萄', '葡萄斑点 GFkV', ['grapevine fleck']),
    ('果树/核果', '樱桃绿环斑 CRLV', ['cherry rasp leaf']),
    ('果树/核果', '樱桃卷叶 CLRV', ['cherry leaf roll']),
    ('果树/核果', '樱桃小果 LChV1', ['little cherry virus 1']),
    ('果树/核果', '樱桃小果 LChV2', ['little cherry virus 2']),
    ('果树/核果', '核果坏死环斑 PBNSPaV', ['plum bark necrosis']),
    ('热带', '香蕉线条 BSV', ['banana streak']),
    ('热带', '木薯曲叶 ACMV', ['african cassava mosaic']),
    ('热带', '木薯褐条 CBSV', ['cassava brown streak']),
    ('热带', '番薯羽状斑 SPFMV', ['sweet potato feathery mottle']),
    ('热带', '甘蔗花叶 SCMV', ['sugarcane mosaic']),
    ('热带', '柑橘裂皮 CEVd(类病毒)', ['exocortis']),
    ('热带', '香蕉束顶 BBTV', ['bunchy top']),
    ('谷类', '大麦黄矮 BYDV', ['barley yellow dwarf']),
    ('谷类', '小麦条纹花叶 WSMV', ['wheat streak mosaic']),
    ('谷类', '水稻东格鲁 RTBV', ['rice tungro']),
    ('谷类', '水稻黑条矮缩 RBSDV', ['black streaked dwarf']),
    ('谷类', '大麦黄色花叶 BaYMV', ['barley yellow mosaic']),
    ('豆类', '大豆花叶 SMV', ['soybean mosaic']),
    ('豆类', '豌豆耳突花叶 PEMV', ['pea enation']),
    ('豆类', '豇豆轻斑驳 CMPMV', ['cowpea mild mottle']),
    ('类病毒', '啤酒花矮化 HSVd', ['hop stunt']),
    ('类病毒', '葡萄黄点 GYSVd-1', ['yellow speckle viroid 1', 'yellow speckle']),
    ('类病毒', '葡萄黄点 GYSVd-2', ['yellow speckle viroid 2']),
    ('类病毒', '桃潜隐花叶 PLMVd', ['peach latent mosaic']),
    ('类病毒', '苹果锈果 ASSVd', ['apple scar skin']),
    ('类病毒', '柑橘曲叶 CVd/CDVd', ['citrus dwarfing', 'citrus viroid']),
    ('其他重要属', '线虫传多面体属 Nepovirus', ['nepovirus']),
    ('其他重要属', '马铃薯Y病毒属 Potyvirus', ['potyvirus']),
    ('其他重要属', '番茄斑萎病毒属 Orthotospovirus', ['orthotospovirus', 'tospovirus']),
    ('其他重要属', '双生病毒属 Begomovirus', ['begomovirus']),
    ('其他重要属', '陷阱病毒属/毛形病毒属 Trichovirus', ['trichovirus']),
    ('其他重要属', '凹陷病毒属 Foveavirus', ['foveavirus']),
    ('其他重要属', '葡萄病毒属 Vitivirus', ['vitivirus']),
    ('其他重要属', '蚕豆病毒属 Fabavirus', ['fabavirus']),
    ('其他重要属', '绿萝病毒属 Potexvirus', ['potexvirus']),
    ('其他重要属', '香石竹斑驳 Carmovirus', ['carmovirus']),
    ('其他重要属', '耳突花叶属 Umbravirus', ['umbravirus']),
    ('其他重要属', '纤细病毒属 Tenuivirus', ['tenuivirus']),
    ('其他重要属', '巨脉病毒属 Varicosavirus', ['varicosavirus']),
]


def main():
    rows = []
    with open(KV_INFO, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rows.append(((r.get('Species_NCBI') or '') + '|' +
                         (r.get('Virus name(s)') or '')).lower())
    print(f'库参考总数: {len(rows)}\n')
    print(f"{'组':<10}{'目标':<38}{'状态':<6}命中示例")
    miss_total = 0
    miss_list = []
    for grp, name, terms in CHECKLIST:
        hits = []
        for t in terms:
            hits = [x for x in rows if t in x]
            if hits:
                break
        if hits:
            sample = hits[0].split('|')[0][:36]
            print(f"{grp:<10}{name:<38}{'✓':<6}{sample}")
        else:
            miss_total += 1
            miss_list.append((grp, name))
            print(f"{grp:<10}{name:<38}{'✗ 缺失'}")
    print(f'\n审计结论: {len(CHECKLIST) - miss_total}/{len(CHECKLIST)} 覆盖, '
          f'缺失 {miss_total} 项')
    if miss_list:
        print('\n缺失清单:')
        for grp, name in miss_list:
            print(f'  - [{grp}] {name}')


if __name__ == '__main__':
    main()
