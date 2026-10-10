# -*- coding: utf-8 -*-
"""宿主库检查：库身份、元数据一致性，以及"库与物种真的对得上"。

为什么需要这个检查（2026-09-11 实际事故）：
  1) 宿主分类库原本是单一固定槽位 host-db/host/classify，建第二个宿主只能
     覆盖/混进同一目录，且库本体不记录物种；
  2) kunpeng 的 add-library 往 seqid2taxid.map **追加**，而 build_db 的替换式
     清理漏了该文件与 taxo.k2d → 陈旧 seqid→taxid 长期累积；
  3) 实测该库的 map 里同时有 4081（番茄 Solanum lycopersicum）与 112863
     （枸杞 Lycium barbarum），taxo.k2d 也带上 Solanum —— 而 ①宿主去除
     只按路径取库、不核对物种，用错宿主是**静默**的。
现在改为按物种建目录 + host_db.json 清单 + active_host_db 指针；本脚本守住
"声明与内容一致"，并用一段真实序列做端到端探针确认库确实认这个物种。

用法:
    python tests/_check_hostdb.py             # 身份检查 + 分类探针
    python tests/_check_hostdb.py --no-probe  # 只做身份检查（不调 kunpeng）
    python tests/_check_hostdb.py --rebuild   # 追加"同目录换 taxid 重建"语义测试
"""
import io
import json
import os
import shutil
import sys

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Virus_Platform_Core import config as C          # noqa: E402

fails = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        fails.append(msg)


def section(t):
    print(f'\n[{t}]')


def probe_genome(genome, out_fa, bp=20000):
    """取参考基因组第一条序列的前 bp 个碱基做探针。"""
    head, chunks = None, []
    total = 0
    with io.open(genome, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if head is not None:
                    break
                head = line[1:].strip()
                continue
            if head is None:
                continue
            chunks.append(line.strip())
            total += len(line.strip())
            if total >= bp:
                break
    seq = ''.join(chunks)[:bp]
    os.makedirs(os.path.dirname(out_fa), exist_ok=True)
    with io.open(out_fa, 'w', encoding='utf-8') as f:
        f.write(f'>probe_{head}\n')
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + '\n')
    return head, len(seq)


def rebuild_semantics():
    """同目录换 taxid 重建：map 里只应剩最后一次的 taxid（替换语义）。"""
    from Virus_Platform_Core import kunpeng
    work = os.path.join('host-db', '_check_hostdb_probe')
    db = os.path.join(work, 'db')
    shutil.rmtree(work, ignore_errors=True)

    def mk(path, names):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with io.open(path, 'w', encoding='utf-8') as f:
            for n in names:
                f.write(f'>{n}\n')
                unit = 'ACGTACGGTTACCGATCGATCGATCGGCTA' * 1700
                for i in range(0, len(unit), 70):
                    f.write(unit[i:i + 70] + '\n')

    fa_a = os.path.join(work, 'a.fa')
    fa_b = os.path.join(work, 'b.fa')
    mk(fa_a, ['seqA1', 'seqA2'])
    mk(fa_b, ['seqB1', 'seqB2'])
    try:
        kunpeng.build_host_db(fa_a, 3702, db_dir=db, clean_mid=True)
        check(sorted(C.host_db_taxids(db)) == [3702],
              '首次建库 map 只含本次 taxid(3702)')
        kunpeng.build_host_db(fa_b, 112863, db_dir=db, rebuild=True,
                              clean_mid=True)
        t = sorted(C.host_db_taxids(db))
        check(t == [112863],
              f'同目录换 taxid 重建后 map 只剩 112863（实际 {t}）——'
              f'替换式重建未把陈旧 seqid→taxid 留下')
        man = C.host_db_manifest(db)
        check(man.get('taxid') == 112863,
              '清单 taxid 随重建更新')
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    args = set(sys.argv[1:])
    do_probe = '--no-probe' not in args
    cfg = C.get_config()
    active = C.get_active_host_db()

    section('配置与目录')
    check(bool(active) and os.path.isdir(active),
          f'当前宿主库存在: {active or "(未配置/未找到)"}')
    check(os.path.normcase(cfg.databases['host']) == os.path.normcase(active or ''),
          'databases["host"] 与 active_host_db 解析一致')
    if not active:
        return 1
    info = C.host_db_info(active)
    from Virus_Platform_Core.kunpeng import db_ready
    check(db_ready(active), f'库结构完整（hash + opts/taxo/hash_config）: '
                            f'{info["name"]}')

    section('库身份')
    man = C.host_db_manifest(active)
    check(bool(man), '库内有 host_db.json 清单')
    taxids = sorted(C.host_db_taxids(active))
    if man:
        print(f'       声明: taxid={man.get("taxid")} '
              f'species={man.get("species")!r} '
              f'source={os.path.basename(str(man.get("source_genome")))}')
    print(f'       seqid2taxid.map 内 taxid: {taxids}')
    check(len(taxids) == 1,
          f'map 内只有 1 个 taxid（实际 {len(taxids)} 个: {taxids}）——'
          f'多于 1 个说明库元数据被历史建库污染，应换新目录重建')
    if man and taxids:
        check(man.get('taxid') in taxids,
              '清单声明的 taxid 与 map 内容一致')

    section('旧布局')
    legacy = C.legacy_host_db_dir()
    if os.path.isdir(legacy):
        li = C.host_db_info(legacy)
        print(f'       legacy 槽位存在: {legacy}')
        print(f'       taxid={li.get("taxid")} species={li.get("species")!r} '
              f'map={li.get("taxids_in_map")}')
        if li.get('conflicted'):
            print('       （已知污染：该目录混入过 4081 番茄建库；'
                  '新库不受影响，未设为当前即不参与分析）')
    else:
        print('       （无 legacy 槽位）')

    if do_probe:
        section('分类探针（库是否真的认这个物种）')
        genome = str(man.get('source_genome') or '') or (
            C.current_host_genome() or '')
        if not genome or not os.path.isfile(genome):
            check(False, f'拿不到源基因组做探针: {genome or "(空)"}')
        else:
            from Virus_Platform_Core import kunpeng
            pf = os.path.join('run', '_check_hostdb', 'probe.fa')
            head, n = probe_genome(genome, pf)
            out = os.path.join('run', '_check_hostdb', 'out')
            try:
                res = kunpeng.classify(active, [pf], out, threads=4)
                leaves = []
                with io.open(res['kreport'], encoding='utf-8',
                             errors='replace') as f:
                    for line in f:
                        p = line.rstrip('\n').split('\t')
                        if len(p) >= 6:
                            leaves.append((p[3], p[4], p[5].strip()))
                assigned = leaves[-1][1] if leaves else None
                check(assigned == str(man.get('taxid')),
                      f'探针（>{head} 前 {n}bp）判为 {assigned} '
                      f'= 声明的 taxid {man.get("taxid")}')
                # 原断言写死「不得出现 Solanum」——那是 2026-09-11 枸杞库被
                # 番茄污染的事故口径；番茄库（taxid=4081）的本名出现是正确
                # 行为。普适化：任何命名节点的 taxid 都必须等于声明 taxid
                # （宿主库单物种，异种节点 = 污染）。
                _want = str(man.get('taxid'))
                # 污染判据 = 出现「另一个物种级节点」：kreport 的
                # R/K/P/C/O/F/G 等层级是声明物种的祖先链（单物种库必然
                # 带，合法）；只有 S/S1… 物种级出现声明之外 taxid 才是
                # 他种污染（2026-09 枸杞库混入番茄即此形态）
                bad = [x for x in leaves
                       if (x[0] or '').upper().startswith('S')
                       and x[1] not in ('0', '-', _want)]
                check(not bad,
                      f'kreport 出现声明物种之外的节点（污染）: {bad[:3]}')
            finally:
                shutil.rmtree(os.path.join('run', '_check_hostdb'),
                              ignore_errors=True)

    if '--rebuild' in args:
        section('替换式重建语义（--rebuild）')
        rebuild_semantics()

    print('\n' + '=' * 52)
    if fails:
        print(f'FAILED: {len(fails)} 项')
        for f in fails:
            print('  - ' + f)
        return 1
    print('宿主库检查通过 ✔')
    return 0


if __name__ == '__main__':
    sys.exit(main())
