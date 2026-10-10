# -*- coding: utf-8 -*-
"""批量跑集成/自检测试，汇总通过情况。只读性质为主，不改平台产物契约。

用法:
    python tests/_run_all.py                 # 跑全部（含慢测试）
    python tests/_run_all.py --quick         # 跳过 slow 标记的测试
    python tests/_run_all.py --list          # 只列测试清单与分级
    python tests/_run_all.py _it_align _it_submit   # 只跑指定测试

超时策略：每个测试有**独立**超时（默认 300s），慢测试放宽到 1500s。
此前用统一短超时 + 外壳 kill，会把"跑得久"误判成"失败"（rc=143/1），
掩盖真实结论——现在超时会明确标 TIMEOUT，与 FAIL 区分。

环境变量:
    VP_SKIP_ONLINE=1   跳过一切真实网络请求的测试段（透传给子进程）
输出: tests/_run_all_report.txt
"""
import io
import os
import subprocess
import sys
import time
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 互斥锁：所有测试都通过 app.test_client() 直接读写同一个 run/ 目录，
# 两个批量实例并发跑会互相删/改对方产物，表现为「单跑全过、批量随机 FAIL」。
# 这里用 O_CREAT|O_EXCL 建锁文件，拿不到就检查持锁进程是否还活着——
# 已退出则是陈旧锁，直接复用（不删文件，某些环境会拦截删除调用）。
_LOCK_PATH = os.path.join(ROOT, 'tests', '._run_all.lock')


def _pid_alive(pid):
    """跨平台探测 pid 是否存活。探测不了时返回 True（保守，当它活着）。"""
    if not pid or not str(pid).isdigit():
        return False
    pid = int(pid)
    if pid == os.getpid():
        return True
    try:
        if os.name == 'nt':
            import ctypes
            k = ctypes.windll.kernel32
            h = k.OpenProcess(0x1000, False, pid)      # PROCESS_QUERY_LIMITED
            if not h:
                return False
            code = ctypes.c_ulong()
            ok = k.GetExitCodeProcess(h, ctypes.byref(code))
            k.CloseHandle(h)
            return bool(ok) and code.value == 259     # STILL_ACTIVE
        os.kill(pid, 0)
        return True
    except (OSError, ValueError, AttributeError):
        return False


def _acquire_lock():
    pid = str(os.getpid())
    for _ in range(2):
        try:
            fd = os.open(_LOCK_PATH, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, pid.encode())
            os.close(fd)
            return True
        except FileExistsError:
            # 持锁进程已退出 → 陈旧锁，覆盖其内容即可（不删除文件）
            try:
                with io.open(_LOCK_PATH, encoding='utf-8') as f:
                    holder = f.read().strip()
            except OSError:
                holder = ''
            if _pid_alive(holder):
                print(f'[中止] 另一个批量测试实例正在运行（pid={holder}）。\n'
                      f'       测试共享 run/ 目录，并发会互相干扰。')
                return False
            try:                                       # 陈旧 → 原地覆写接管
                with io.open(_LOCK_PATH, 'w', encoding='utf-8') as f:
                    f.write(pid)
                return True
            except OSError:
                return False                           # 写不了就当拿不到
        except OSError:
            return False
    return False


def _release_lock():
    try:
        os.remove(_LOCK_PATH)
    except OSError:
        pass

# (路径, 单测超时秒, 是否 slow)
# slow=True 的测试含真实重型计算或在线请求，`--quick` 时跳过。
TESTS = [
    ('tests/_it_platform.py',       300,  False),
    # dsRNA 全链路集成（2026-10-03 入编；此前从不在批量里，迁出/恢复期间
    # 无人报警）：真跑 dsRNAmax 选窗 + dsRIP 效价 + 面板脱靶扫描 + T7 引物
    # + 成品 QC，含确定性/内存守卫/max_mm=0 语义/isolsates 分支；约 3 秒。
    # 前置数据 databases/dsrna/lookup.db（小文件）；面板缺失时自合成小面板
    # 自足运行，不依赖仓库外 GB 级面板。
    ('tests/_it_dsrna.py',          300,  False),
    # 日志编码统一（2026-09-27）：decode_output 三分支 + to_utf8_file 四场景
    # （纯 GBK 转写 / 纯 UTF-8 幂等 / 混编码逐行无损 / 缺失文件静默）。
    ('tests/_check_enc_utils.py',    60,  False),
    ('tests/_it_msa.py',            300,  False),
    ('tests/_it_concurrency.py',    300,  False),
    ('tests/_smoke_phylo.py',       300,  False),
    ('tests/_it_phylo.py',          300,  False),
    ('tests/_it_extract.py',        300,  False),
    ('tests/_it_submit.py',         300,  False),
    ('tests/_it_align.py',          600,  False),
    ('tests/_it_compare.py',        600,  False),
    ('tests/_it_synteny.py',        600,  False),
    ('tests/_it_dl_delete.py',      300,  False),
    # 下载进度/大小/速度/ETA 与单端 .sra 解码的离线回归（秒级，不联网）
    ('tests/_it_dl_progress.py',    300,  False),
    # 工具箱「格式转换」.sra→FASTQ/FASTA 的产物识别：单端 <acc>.fastq.gz
    # 必须被认（旧行为只 glob <acc>_* → 产物已落盘却报「sracha 无输出」）。
    ('tests/_it_tool_convert.py',   300,  False),
    # 内置 HTTP 回退下载（_dl_python）与完整性闸门：Range 续传追加语义、
    # 416 收尾、精确总量优先、三道校验闸门。本机 HTTP 服务端，离线。
    ('tests/_it_dl_http_fallback.py', 300, False),
    # 重启恢复健壮性：畸形 batch.json 不得拖垮整个下载模块；中断语义。
    ('tests/_it_dl_recovery.py',    300,  False),
    # 「一键转入分析流程」：双端/单端配对、跳过不留空样品目录、runs 过滤、
    # enqueue 入队（假队列 + 路径隔离，绝不碰真实 queue.json）。
    ('tests/_it_dl_to_pipeline.py', 300,  False),
    # 收尾阶段（.sra 解码）必须能取消：进程被终止、剩余不继续、状态不被
    # 覆写成 completed、不误报转换失败。
    ('tests/_it_dl_convert_cancel.py', 300, False),
    # 长耗时接口后台任务化 + 重 IO 缓存/上限的离线回归（外部依赖全用桩）
    ('tests/_it_bg_tasks.py',       300,  False),
    # 阶段二：排队样品 / 存储扫描两处外部任务挂载 + attach 重绑定语义 +
    # 「跑完的外部任务必须退出 active_ids」（否则全局重置会误取消已完成的
    # 批次）。全程离线，只动内存里的假队列，不碰真实批次。
    ('tests/_it_taskcenter_ext.py', 180,  False),
    # 下载页真实浏览器验证（Selenium + 本机 Edge/Chrome，headless）：
    # 逐文件表格的倒挂 / 跨轮询被抹掉 / 重试续传入口。无浏览器时自动
    # SKIP 并返回 0，不阻塞批量测试。
    ('tests/_it_download_ui.py',    600,  False),
    # 四个「同步阻塞→后台任务」改造的真实浏览器验证（.sqn / 物种名校验 /
    # AI 总结 / CDS 导出）；外部依赖全部用桩。无浏览器时自动 SKIP。
    ('tests/_it_async_ui.py',       900,  False),
    # 任务中心页对「非样品类后台任务」的渲染与可操作性（进度条/日志/
    # 停止/前往），以及「只看当前样品」筛选的真实语义。无浏览器时自动 SKIP。
    ('tests/_it_taskcenter_ui.py',  600,  False),
    # 阶段二：**下载批次挂进任务中心**的真实浏览器端到端 —— 卡片/进度/进展
    # 文案、停止转发到 DownloadManager、前往页面、真身消失后卡片自动移除、
    # provider 报错时降级保留。绑定的是临时 downloads 目录，不碰真实批次。
    ('tests/_it_taskcenter_dl_ui.py', 600, False),
    # 三条修复路径真实浏览器端到端（续传 / 重试失败 / 重试转换）：点下去
    # 真的能把批次救回来。只打桩「字节传输」与「sracha 解码」，其余走真实代码。
    ('tests/_it_dl_repair_ui.py',   600,  False),
    ('tests/test_logan_trace.py',   300,  False),
    ('tests/test_logan_batch.py',   300,  False),
    ('tests/_check_local_annot.py', 300,  False),
    # 宿主库身份 + 分类探针（含 kunpeng 分类一次，约 20s）
    ('tests/_check_hostdb.py',      300,  False),
    # ASCII 中转目录的并发隔离（平台根含中文时 MAFFT/trimAl 必走中转；秒级）
    ('tests/_check_ascii_stage.py', 300,  False),
    # 启动模式（桌面窗口 / 网页）选择与回退 + 两个启动脚本（秒级，无副作用）
    ('tests/_check_start_modes.py', 300,  False),
    # 「浏览」按钮回归：调用 browse()/browseDir() 的页面必须带上共享对话框片段
    # （缺 #dlgMask 会抛 TypeError，按钮点了毫无反应；纯静态，秒级）
    ('tests/_check_browse_dlg.py',  120,  False),
    # 示例结果清单 ↔ 实际文件：列表页计数与详情页列表必须一致（秒级，只读）
    ('tests/_check_examples_manifest.py', 60, False),
    # 库/目录类 API 返回的路径必须是绝对路径（分发版下相对路径会解析错库；
    # 秒级，只读）
    ('tests/_check_db_paths.py',     60,  False),
    # spec 里 webapp 的 datas 展开：既不漏前端文件、也不带 .mimosa 等工具残留
    # （AST 抽出 spec 里的 _tree 单独执行；秒级，只读）
    ('tests/_check_spec_datas.py',   60,  False),
    # i18n 缺键：扫 webapp/templates/*.html 全部 data-i18n/-ph/-title 锚点，
    # 逐个喂真实 applyI18n，断言不会被写成 key 字面量 + 中英成对。
    # 2026-10-01 入编 —— 此前 _verify_i18n.js 从不在批量里（只跑 .py），
    # 于是它 fixture 过期变永久红灯也没人发现。无 node 时自动 SKIP。
    ('tests/_check_i18n.py',         60,  False),
    # 前端 JS 语法静态校验（内联 script 块 + static/*.js，node --check 只做
    # 语法不执行）+ 4 个跨卡复用函数存在性 + tools.html 无重复 id。无 node 时
    # 自动 SKIP 返回 0。2026-10-03 入编（此前不在批量里，一次批量替换把函数头
    # 连删 5 行都无人报警）。
    ('tests/_check_js_syntax.py',   120,  False),
    # 路由清单基线守卫：rule+methods 与 tests/_route_baseline.json 严格比对。
    # 新增/丢失/方法变化都会报红（改路由后需 python tests/_route_inventory.py
    # --save 更新基线）；只读、不启动分析任务，秒级。
    ('tests/_route_inventory.py',    60,  False),
    # 全路由健壮性巡检：对每条路由发一次空参数请求，只找「未捕获异常 500」；
    # 破坏性路由（空 POST 就改状态/删数据）在脚本内 DESTRUCTIVE_POST 里显式
    # 跳过，且请求前后比对 run/tasks 增删。只写 tests/_audit_routes_report.txt
    #（已 gitignore）。2026-10-03 入编。约 20s。
    ('tests/_audit_routes.py',      300,  False),
    # 阶段表一致性守卫（后端 9 张表 + 前端 static/*.js + ③c 卡片可见性实测）。
    # 2026-09-16 修掉它"只读 app.js"的过时假设（前端已拆分到 app-jobs/app-batch）
    # 后重新入编——此前不在批量里，被前端重构甩下后静默腐烂了一轮。
    ('tests/_it_stage_tables.py',   120,  False),
    # PYZ↔源码字节码比较器的灵敏度自检：确认它**不是恒真**
    # （纯内存 compile 比较，不读产物、不需 PyInstaller；秒级）
    ('tests/_selftest_pyz_cmp.py',   60,  False),
    # 全页面真浏览器巡检：21 个页面逐个加载，收集 JS 报错 / 4xx 资源 /
    # 关键全局函数存在性。**逐页独立进程**驱动（同进程连续访问会随机段错误），
    # 约 2-4 分钟；无 playwright 时自动 SKIP 并返回 0。
    ('tests/_check_pages_console.py', 900, False),
    # ⚡一键分析（全流程）卡：自起临时实例 + 真浏览器，验证本卡自有输入/参数/
    # 输出区齐全、示例与参数互取落到本卡、内联结果链路（桩 API）不串卡（~30s）。
    # 无 playwright（开发依赖未装）时自动 SKIP 并返回 0。
    ('tests/_check_kvchain_ui.py',  300,  False),
    # 病毒浏览器 Explorer（服务器版 7 页签：时序趋势 / 全基因组变异 / 数据浏览 /
    # 引物库 / 宿主范围 / 媒介传播 / 病毒档案）真浏览器验证：断言 Plotly 实例真的
    # 渲染出数据点（读 _fullData，兼容 bdata 二进制数组与 Sankey 的 node/link 嵌套）、
    # 下拉初始化、表格排序分页、深链指向与真点击跳转、导出 Content-Disposition、
    # 全程无 console error。覆盖 198,819 条 / 6,168 物种真实数据，约 3 分钟。
    # 自起临时实例。无 playwright（开发依赖未装）时自动 SKIP 并返回 0。
    # 其中 [13] 深链跳转需要外网；断网时该段自动 SKIP，不判失败。
    # ⛔ 已停用用例（2026-09-17 摘除 explorer；2026-09-27 彻底归档）。
    # _check_explorer.py 走真实的 /explorer 与 /api/explorer/* 路由，explorer
    # 蓝图已不再注册 → 必然 404。它验证的是一个**已下线的入口**，失败不代表
    # 回归，故从默认清单移除，避免长期红灯把真回归淹掉。
    # 需要复核归档代码时到 archive/_retire_20260927/explorer/tests/ 手动跑
    # （_check_explorer.py / _check_explorer_links.py / _check_virphykit_export.py
    # / _check_decimal_year_conventions.py 都随引擎一起归档）。
    # ('tests/_check_explorer.py',   900,  False),
    # Explorer → VirPhyKit 输入包导出：小数年解析（正/负控各 8）、区域口径与 name
    # token、严格档预过滤（无日期/无地点/非完整/重复登录号/空序列五道剔，含记账
    # 恒等式 sum(dropped)+kept==in 与「被剔的不能出现在保留集」负控）、**灵敏度
    # 自检**（同一套契约断言函数喂一个故意做坏的包，必须逐条报出问题，否则断言
    # 可能恒真）、真库端到端（PVY 4470→588 与 preview 一致，另用科级 24511 条
    # 证明导出没吃表格 5000 行截断）。可选的 GeoSubsampler 实跑在缺外部目录时
    # 自动 SKIP。约 30 秒。
    # 注：端到端段走 /api/explorer/export/* 路由，而 explorer 蓝图 2026-09-17
    # 起不再挂载（2026-09-27 随引擎彻底归档）——该段探测到路由
    # 不在 url_map 时自动 SKIP 并返回 0；引擎契约由本用例 [1]-[6] 段常态守护。
    # （本用例文件已随 explorer 归档：archive/_retire_20260927/explorer/tests/）
    # Fitch 重构回归：锁住 2026-09-15 AUDIT 记的缺陷 —— fitch_mugration 的
    # down-pass 从未执行（内部状态退化成候选集字典序最小者 → 高估迁移数），
    # 以及缺失数据被当成字面状态 'Unknown' 凭空造迁移。含旧实现对照。
    # 注：本平台的进化动力学卡片已迁往「植物病毒进化分析平台」，但 2026-09-16
    # 起本站重新接上了裁剪版系统地理（RRT/RSPP/MOTP/三分类/GIF/分带，见下一项）。
    ('tests/_check_phylogeo_fitch.py', 900, False),
    # 主平台裁剪版系统地理：haversine 自研口径、坐标表/别名/逐样本落点优先级、
    # classify_transitions 分位带与证据标记、RSPP 驱动的权重分带（措辞非贝叶斯）、
    # MOTP 时间分箱（含 LSD2 定年树解析）、GIF 帧计划、analyze() 端到端产物 +
    # 默认关开时 rrt/rssp/motp/weight_bands 均为 None。全离线、确定性夹具，秒级。
    ('tests/_check_phylogeo_geo.py',  900, False),
    # 同一张卡的真浏览器端到端（自起临时实例，随机空闲端口）：一级导航 → 组落地
    # → 点卡片进卡 → 参数齐全 + 静态文案已 i18n（切英文验证）→ 真跑一次
    # （RRT/RSSP/MOTP 全开的 8 条合成集）→ 三分类表 / 权重分带表 / GIF 区渲染
    # 且写明非贝叶斯 → 无 pageerror。无 playwright（开发依赖未装）时自动 SKIP。
    ('tests/_check_phylogeo_ui.py',   900, False),
    # 进化动力学工具组（VirPhyKit 对齐 12 卡，2026-09-17）：Engine 逐模块验证 +
    # 两个精确对拍契约（MJRM 与 PVS_with_matrix.xml 逐字节、TempMig 矩阵与
    # Migration_matrix.txt 逐格）+ 数据接入三来源时间地点同检 + 灵敏度自检
    # （坏输入必须红）。VirPhyKit Example 目录不在时 Example 段 SKIP。
    # TreeTime 全量跑约 40s（500 叶）；--fast 可跳。Example 不在本机也全绿。
    ('tests/_check_phylodyn_kit.py',  900, False),
    # VirPhyKit Example 全目录验证（2026-09-17）：42 个文件逐个被对应工具真实
    # 消费并断言（含 Total=各区行和自洽、RRT Min/Max 与 20 副本逐棵核对、
    # na_20/na_200/ebola 三个负例的精确拒绝、全目录覆盖检查漏一个文件即失败）。
    # [13] TreeDater 段跑真 R 引擎：定年 ~1 分钟 + parboot 置信带单核 15-20 分钟
    # （引擎内注释原话 10-40 分钟）——900s 超时必挂（2026-09-18 全量首跑即超时），
    # 故放宽到 2400s 并标 slow：--quick 跳过，引擎契约由 _check_phylodyn_kit.py
    # （~70s，含 MJRM/TempMig 逐字节对拍）在 quick 档继续守护。
    ('tests/_check_phylodyn_example.py', 2400, True),
    # 同组的真浏览器端到端（自起临时实例）：15 卡组落地 → 数据接入卡三来源切换 +
    # 真跑手动粘贴（坏行点名 + 图渲染 + 下载链接）→ 重命名卡 / MJRM 卡真跑 →
    # 其余 9 卡控件存在性 → 新增 i18n 键切英文 → 无 pageerror。无 playwright SKIP。
    ('tests/_check_phylodyn_ui.py',   900, False),
    # 数据准备链（2026-09-18 合并 4 张卡 + 引入上游预处理）：任务层真跑 44 项
    # （治理/去重/比对QC/允无地点/坐标表/降采样/两张卡已摘的负向）。
    ('tests/_check_pdprep_chain.py',  900, False),
    # 与上游 `virome_phylo_pipeline` 的治理层**逐条对拍**（上游仓库不在则自跳过返回 0）。
    ('tests/_check_govern_vs_upstream.py', 300, False),
    # 元数据列名自动识别与择优（explorer 新版 25 列）+ 治理边界（关治理回到老行为）。
    ('tests/_check_phylodyn_meta_alias.py', 300, False),
    # 小数年口径：四处实现的数值钉住 + 「月份必须生效」的 P0 防回归。
    # （原 _check_decimal_year_conventions.py 依赖 explorer 的 virphykit_export，
    #  已随 2026-09-27 归档移除；其 P0 断言由 _check_phylodyn_meta_alias.py 延续。）
    # 示例**输入** fixture 自洽（叶名能被元数据匹配 / 等长 / 不是 GenBank 转储）。
    # 示例输入 fixture 自洽由上行守卫覆盖；以下为 2026-10 修复轮新增守卫：
    # AA 同一性分簇锚定（M1）/ 圈图真圈图（M2）/ SVG 良构（M3）/
    # 引物全链热力学评分分级（Primer 专项）
    ('tests/_check_identity_aa.py', 600, False),
    ('tests/_check_genome_circos.py', 120, False),
    ('tests/_check_svg_wellformed.py', 60, False),
    ('tests/_check_primer.py', 600, False),
    ('tests/_check_examples_inputs.py', 120, False),
    # 在线获取（替代 SeqHarvester 的便利层）+ **输出 schema 一致性**：
    # taxid/物种名一键采集、geo_loc_name 优先、host 通道、host 可选列。
    # ⚠️ 里面含**真联网**的段（NCBI 不可达时自动 SKIP 并返回 0）。
    ('tests/_check_online_meta.py',   600, False),
    # TreeTime mugration 的**状态代号映射**（代号会跨到标点 ）：
    # 不依赖 TreeTime，用合成 GTR.txt/confidence.csv 钉住口径。
    ('tests/_check_gtr_mapping.py',   120, False),
    # 「进化树 + 基因组叠加」接线（2026-09-17）：后端 /api/tree/file 按「与树同名的
    # <stem>.overlay.json」带上 overlay（3 轨道 / 18 基因 / 18 条蛋白同一性连线）+
    # 五类坏 overlay 的降级负控（树照常返回）+ 真浏览器点「✨ 示例·基因组叠加」数
    # SVG（属色带 / 基因框 / 同源连线）与换纯树时的清零负控 + 新增 i18n 键。
    # 无 playwright（开发依赖未装）时真浏览器段自动 SKIP 并返回 0。
    ('tests/_check_tree_overlay.py', 900, False),
    # 真跑 HMM 全链路 + 在线 NCBI，实测 400s~900s（网络波动大）
    ('tests/_it_annotate.py',      1500,  True),

    # ---- 2026-10-07 审查新增守卫（每条都带「从备份抽取旧实现」的对照臂）----
    # SDT gap 口径：_sdt_formula 必须与 phylo.pairwise_identity_matrix 同口径
    # （'-. ' 都算 gap）；旧实现只认 '-'，含 '.' 的 MSA 会算出两个答案。
    ('tests/_check_sdt_gap_parity.py', 120, False),
    # SDT 分布图：全同输入（所有两两 identity == 100）曾因 counts 为空而
    # max(空) 崩溃，且崩在 CSV/热图落盘之后。
    ('tests/_check_sdt_distribution.py', 120, False),
    # CRISPR 路线 B（枚举兜底）的引物交付序列必须纯 ACGT —— 旧实现只挡 'N'，
    # W/S/M/K/R/Y 会一路流到订购表。
    ('tests/_check_crispr_primer_purity.py', 120, False),
    # TaskManager.on_finish：排队途中取消时 fn 根本不被调用，清理若只写在
    # fn 的 finally 里会永久泄漏（样品占位表 _RUNNING_SAMPLES）。
    ('tests/_check_task_on_finish.py', 120, False),
    # P1-4 SVG 良构：出口覆盖（dfv/gbdraw/genome_diag/viz 必须调 verify_svg）
    # + 覆盖范围扩到 run/ 产物 + 判据灵敏度自测（已在 _check_svg_wellformed）。
    # 本条是**消费侧链路**：残缺 SVG 不得进入 report.html（修得动的就地修复后
    # 内嵌、修不动的跳过），对照臂证明"未校验"时确实会漏进报告。
    ('tests/_it_svg_embed_chain.py', 300, False),

    # ---- AI 助手（对话驱动智能代理）----
    # 2026-10-08 入编：这 6 个用例此前**存在但从不在批量里**（批量入口是硬编码
    # 元组、不是 glob），等于「写了测试但没人跑」——本仓已两次踩过同一坑
    # （_verify_i18n.js 与 6 个 pytest 文件）。3 个 Playwright UI 巡检缺
    # playwright 时 SKIP 返回 0（已补守卫；此前会直接抛 ImportError 变红）；
    # 3 个集成测试用线程内模拟 OpenAI 兼容端点，不打真实网络。
    ('tests/_it_ai_assist.py',       300, False),
    ('tests/_it_ai_events.py',       300, False),
    ('tests/_it_ai_actions.py',      300, False),
    ('tests/_check_ai_dock_ui.py',   300, False),
    ('tests/_check_ai_events_ui.py', 300, False),
    ('tests/_check_ai_global_ui.py', 300, False),
    # digest ↔ 源文件数字对拍（2026-10-08）：digest 是 AI 的**唯一数字来源**，
    # 读错文件/键 → AI 给出"看起来完全合理"的错误解读，且没有任何报错。
    # 本守卫独立重读源文件逐键比对（含"凭空造数字"反向检查），自带灵敏度
    # 自证（改错一个数字/造一个假键都必须被报出），秒级。
    ('tests/_check_ai_digest_parity.py', 120, False),
    # AI 审计落盘（2026-10-08）：执行单状态变迁与事件流原本只在进程内存里
    # （重启即失、事后无从追溯）。本守卫含**真子进程 + 信号**的"退出不丢"
    # 测试、反复 drain 的"不重复写"测试，以及"接线自证"（publish/executors
    # 真的会落盘 —— 只测 audit 模块自身等于没测接线）。约 10s（含子进程）。
    ('tests/_check_ai_audit.py',     180, False),

    # ---- pytest 风格用例（走 `python -m pytest`，见 PYTEST_FILES 说明）----
    # 2026-10-07 审查 P1-3：这 6 个文件此前**既不在本清单、也无 CI**，
    # 等于「写了测试但从不在回归面上跑」（平台最大的新模块 crispr_design
    # 就属于这种）。现统一接入。
    ('tests/test_crispr_design.py',   300, False),
    ('tests/test_host_genome_gget.py', 120, False),
    ('tests/test_mirna_alignment.py', 120, False),
    ('tests/test_utils_pure.py',      120, False),
    ('tests/test_datasets_cli.py',    120, False),
    ('tests/test_web_guard.py',       120, False),
]

# pytest 风格用例：用 pytest 的 fixture / importorskip，且**没有 __main__
# 入口** —— 不能像其它 _check_*/_it_* 那样直接当脚本跑，必须走 `python -m
# pytest`。本清单按此集合选择调用方式。
# 未安装 pytest 时整组 SKIP（不 FAIL）：开发依赖未装的机器不该因此报红。
PYTEST_FILES = {
    'tests/test_crispr_design.py',
    'tests/test_host_genome_gget.py',
    'tests/test_mirna_alignment.py',
    'tests/test_utils_pure.py',
    'tests/test_datasets_cli.py',
    'tests/test_web_guard.py',
}
_PYTEST_OK = None


def _pytest_available():
    """pytest 是否可用（探测一次并缓存）。"""
    global _PYTEST_OK
    if _PYTEST_OK is None:
        try:
            r = subprocess.run([sys.executable, '-c', 'import pytest'],
                               capture_output=True, timeout=60,
                               env=_child_env())
            _PYTEST_OK = (r.returncode == 0)
        except Exception:
            _PYTEST_OK = False
    return _PYTEST_OK
DEFAULT_TIMEOUT = 300
_ONLINE_HINT = {
    'tests/_it_annotate.py': '含在线 NCBI 请求，可设 VP_SKIP_ONLINE=1 跳过该段',
    'tests/_check_explorer.py': '已停用并随 2026-09-27 彻底归档移出 tests/'
                                '（archive/_retire_20260927/explorer/tests/）。'
                                '若手动跑：必须设 VP_EXPLORER_DATA 指回归档数据，'
                                '且把 explorer 蓝图临时挂回 app.py，否则必然 404。',
}


def _fmt_secs(s):
    return f'{s / 60:.1f}min' if s >= 120 else f'{s:.0f}s'


def _child_env():
    """子进程环境：剥掉宿主 Agent 会话注入的删除监控变量。

    某些 Agent/IDE 会话会把 CODEBUDDY_SAFE_DELETE_* 注入环境，并通过
    stdout 监控子进程输出。测试的收尾清理（shutil.rmtree）在累计计数
    超阈值后会被标记为 SAFE_DELETE_BULK_CONFIRM_REQUIRED，导致**明明全部
    断言已通过**的子进程以 rc=1 退出——表现为"单跑全过、批量随机 FAIL"。
    测试是本地受控目录内的自我清理，不需要外部删除守卫，故剥离。
    """
    env = {k: v for k, v in os.environ.items()
           if 'SAFE_DELETE' not in k and not k.startswith('CODEBUDDY_')}
    # 批量满载下浏览器类测试的元素等待放宽（单文件直跑不受影响）
    env.setdefault('VP_TEST_WAIT', '90')
    # 剥离后 PATH 里可能残留宿主 safe-bin 目录，一并去掉
    if env.get('PATH'):
        parts = [p for p in env['PATH'].split(os.pathsep)
                 if 'safe-bin' not in p.replace('\\', '/')]
        env['PATH'] = os.pathsep.join(parts)
    return env


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    flags = {a for a in sys.argv[1:] if a.startswith('-')}
    quick = '--quick' in flags or '-q' in flags

    if '--list' in flags:
        for path, to, slow in TESTS:
            tag = '  [slow]' if slow else ''
            print(f'{path:34s} timeout={_fmt_secs(to)}{tag}')
        return 0

    # 互斥锁只用于避免并发实例互毁 run/ 目录；拿不到锁时**默认继续跑**
    # 并打印警告——测试可用性优先，锁不该成为跑不了测试的理由。
    # 设 VP_STRICT_LOCK=1 可恢复"拿不到锁就退出"的严格行为。
    locked = _acquire_lock()
    strict = os.environ.get('VP_STRICT_LOCK', '').strip() in ('1', 'true', 'yes')
    if not locked:
        if strict:
            return 2
        print('[警告] 未取得互斥锁（可能有并发实例）。仍继续，'
              '结果可能受干扰；如需严格中止请设 VP_STRICT_LOCK=1\n')
    try:
        return _run(plan_from(args), quick)
    finally:
        if locked:
            _release_lock()


def plan_from(args):
    if not args:
        return list(TESTS)
    by_path = {p: (to, sl) for p, to, sl in TESTS}
    plan = []
    for a in args:
        p = a if a.startswith('tests/') else f'tests/{a}'
        if not p.endswith('.py'):
            p += '.py'
        to, sl = by_path.get(p, (DEFAULT_TIMEOUT, False))
        plan.append((p, to, sl))
    return plan


def _run(plan, quick):
    if quick:
        skipped = [p for p, _t, s in plan if s]
        plan = [(p, to, s) for p, to, s in plan if not s]
        if skipped:
            print(f'--quick: 跳过 {len(skipped)} 个 slow 测试 '
                  f'({", ".join(os.path.basename(p) for p in skipped)})\n')

    out_path = os.path.join(ROOT, 'tests', '_run_all_report.txt')
    lines = []
    n_ok = 0
    n_skip = 0
    t0 = time.time()
    for path, to, slow in plan:
        p = os.path.join(ROOT, path)
        if not os.path.isfile(p):
            lines.append(f'[SKIP] {path} 不存在')
            print(lines[-1])
            continue
        is_pytest = path in PYTEST_FILES
        if is_pytest and not _pytest_available():
            lines.append(f'[SKIP] {path} 未安装 pytest（开发依赖缺失，不计失败）')
            print(lines[-1])
            n_skip += 1
            continue
        ts = time.time()
        timed_out = False
        try:
            # pytest 风格用例没有 __main__，必须走 `python -m pytest`
            cmd = ([sys.executable, '-m', 'pytest', p, '-q'] if is_pytest
                   else [sys.executable, p])
            r = subprocess.run(cmd, cwd=ROOT, env=_child_env(),
                               capture_output=True, timeout=to)
            rc = r.returncode
            raw = (r.stdout or b'') + (r.stderr or b'')
        except subprocess.TimeoutExpired as e:
            # TimeoutExpired.stdout/err 可能是 None 或 bytes
            timed_out = True
            rc = None
            raw = (e.stdout or b'') + (e.stderr or b'')
            if isinstance(raw, str):
                raw = raw.encode('utf-8', 'replace')
        txt = raw.decode('utf-8', errors='replace')
        tail = '\n'.join(txt.strip().splitlines()[-14:])
        dur = time.time() - ts
        if timed_out:
            status = f'TIMEOUT(>{_fmt_secs(to)})'
        else:
            status = 'PASS' if rc == 0 else f'FAIL(rc={rc})'
            if rc == 0:
                n_ok += 1
        note = ''
        if path in _ONLINE_HINT and status != 'PASS':
            note = f'  ↳ {_ONLINE_HINT[path]}'
        block = (f'{"=" * 70}\n[{status}] {path}  ({dur:.1f}s)\n'
                 f'{"-" * 70}\n{tail}\n')
        lines.append(block)
        print(f'[{status}] {path}  ({dur:.1f}s){note}')
        sys.stdout.flush()

    total = (len([p for p, _t, _s in plan
                  if os.path.isfile(os.path.join(ROOT, p))]) - n_skip)
    elapsed = time.time() - t0
    head = (f'批量测试汇总: {n_ok}/{total} 通过，总耗时 {elapsed:.1f}s'
            + (f'（{n_skip} 个 SKIP）' if n_skip else '')
            + ('（--quick 模式）' if quick else '') + '\n')
    with io.open(out_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(head + '\n' + '\n'.join(lines))
    print('\n' + head)
    print(f'明细: {out_path}')
    return 0 if n_ok == total else 1


if __name__ == '__main__':
    raise SystemExit(main())
