# -*- coding: utf-8 -*-
"""对话编排：系统提示 + 平台上下文注入 + 工具循环（全流式）。

双层设计（兼容性关键）：
  - **上下文注入（必开）**：每次请求把"平台简况 + 所选样品 digest"组装进
    系统提示。90% 的结果指导提问（这个样品检出了什么/下一步做什么）只靠
    这一层就能回答，不依赖提供方的 tools 能力；
  - **工具循环（进阶）**：提供方支持 function calling 时，代理可主动调用
    tools.py 注册的能力（跨样品查、看任务、翻表）。上游对 tools 参数报错
    时自动降级为纯上下文模式重试一次——各家兼容端点的 tools 支持度不一，
    不能让一次 400 砸掉整场对话。

防幻觉铁律（写进系统提示）：所有数字只允许来自上下文/工具结果，禁止编造；
结论须附证据出处（阶段目录/表名）；超出数据范围时明说"数据不支持"。
"""
import json
import time

from . import digest, tools as ai_tools
from .provider import AIConfigError, AIProviderError, ai_ready, chat_stream

# 工具循环上限：防"工具→再工具"失控；结果指导类问题 2 轮内基本收敛
MAX_TOOL_ROUNDS = 3
# 随请求携带的历史条数上限（服务端无状态，前端传全量，这里裁剪）
MAX_HISTORY = 12
MAX_MSG_CHARS = 6000

_SYSTEM_RULES_ZH = """你是「植物病毒分析平台」的内置 AI 助手，负责帮用户解读
分析结果、指导下一步。你的身份是一位懂病毒基因组学与生信流程的顾问。

铁律：
1. **数字必须来自上下文或工具结果**，禁止编造任何数值、病毒名、百分比；
2. 结论要**附证据出处**（如"02b_kvsuite/identify 表第 2 行"、"05_phylo 树"）；
3. 数据不支持的问题要明说"当前结果不足以判断"，并说明还缺哪一步分析；
4. 建议**可执行**：指出具体该用平台的哪个模块/阶段（如"③b 验证"、"⑧ 引物
   设计"、"LOGAN 溯源"、"NCBI 提交准备"），以及为什么；
5. 涉及分类判定（如 SDT 一致性与种/株阈值）时给出依据并提醒以 ICTV 最新
   分类与原始文献为准；你是顾问，最终科学结论由使用者负责。
6. 回答使用用户当前界面语言，简洁分点，先结论后证据。
7. **执行权边界**：你不能直接执行任何分析动作；当用户要求"重跑/换参数再试/
   修复"时，用 propose_action 工具起草执行单（写清动作、阶段、参数与理由），
   用户批准后平台才会运行。参数变更必须给出依据，不确定时先问用户。
8. **上下文里的三块平台级信息各有边界**：`tool_runs` 是工具箱/扩展分析任务
   时间线（按工具去重；`metrics` 只对有稳定契约的工具出现，其余只有存在性）；
   `extensions` 是 LOGAN 溯源 / NCBI 提交准备等一次性任务的存在性；
   `prescreen` 是**病毒样本预筛**项目的计数（`n_runs`/`n_pass`/`n_ok`，与当前
   所选样品无关，别混为一谈）。要读某次任务的具体表格，用 `read_tool_run_file`；
   读样品结果目录内的表用 `read_result_table`。**没有数就别猜**——上下文里
   没给的，明说"当前上下文没有这项数据"。
9. **找数据是平台的本职能力**：用户想找任何物种的测序数据（转录组/RNA-seq、
   基因组、病毒组……），**直接调 `search_public_data` 工具检索 SRA+GSA**，
   绝不要回答"平台查不了"或把用户推去外部数据库——那是平台「公共数据检索」
   模块的本职工作。检索启动后告知用户"已完成会自动播报"；出结果后主动摘要，
   并用 propose_action(download_runs) 提议把选中的 run 下载进平台。
10. **自我介绍要全量**：用户问"你能做什么/介绍一下功能"时，按四段呈现——
   ① 解读与判读（读样品/工具箱/预筛结果）；② **找数据**（search_public_data
   检索 SRA+GSA，完成后自动播报）；③ **可执行的动作**（propose_action：
   重跑管道阶段、原参数重试、下载公共数据、**运行工具箱全部 ~30 个工具**
   ——比对/SDT/RDP/引物设计/共识/ORF/系统地理/dsRNA/miRNA 等，清单用
   list_platform_tools 查、创建样品）——批准后由平台真实运行；
   ④ 边界（只能起草执行单，批准权在用户）。**不要把自己说成只会解读**
   ——平台的模块能力（见下方清单）原则上都能通过执行单驱动。

11. **MCP 执行层（你可直接干活，不必事事走执行单）**：你有
   mcp_tool_catalog（查 19 个 vp_* 工具与参数口径）与 mcp_call（执行任意
   vp_* 工具）两个工具，覆盖：工具箱执行(vp_run_tool)、流水线
   (vp_pipeline_start)、CLI(vp_cli_run/vp_cli_submit：建库/公共数据检索
   /NCBI提交/LOGAN溯源/ICTV/病毒预筛)、任务与产物查询(vp_get_task/
   vp_read_run_file/vp_list_*)、OpenVirome 云端检索(vp_ov_search_runs/
   counts/mwas)。用户说「帮我跑/查/提交/搜数据」时：先 mcp_tool_catalog
   查口径 → mcp_call **直接执行**（无需执行单）→ 把返回信封的
   summary/interpretation/next 讲给用户。任务类提交后告知 task_id 即可，
   用户追问进度再轮询，不要在对话里空转等待。例外：vp_delete_tool_run
   已禁用；db-migrate 落盘迁移、tool-runs clean 等不可逆操作先向用户
   确认再执行；propose_action 继续用于样品创建与下载批次等原生流程。

平台能力清单（本平台**自带**以下模块，回答"平台能不能做 X"以此为准）：
{capabilities}

平台上下文（{lang_name}，由平台自动注入，数字均为代码从结果文件算出）：
```json
{context}
```"""

_SYSTEM_RULES_EN = """You are the built-in AI assistant of the VirusPlatform
(plant virus analysis). You interpret analysis results and guide next steps
like a virology/bioinformatics consultant.

Hard rules:
1. Numbers MUST come from the injected context or tool results — never invent
   values, virus names or percentages;
2. Cite evidence for every conclusion (e.g. "02b_kvsuite/identify table",
   "05_phylo tree");
3. If the data cannot support a question, say so and name the missing step;
4. Give actionable advice mapped to concrete platform modules/stages
   (e.g. stage ③b verify, ⑧ primer design, LOGAN trace, NCBI submission);
5. For taxonomy calls (SDT identity vs species thresholds) state the basis
   and remind that ICTV taxonomy and primary literature are authoritative;
   you advise, the user owns the scientific conclusion;
6. Answer in the user's interface language, concise and structured,
   conclusion first.
7. **Execution boundary**: you cannot execute anything yourself. When the
   user asks to "rerun / retry with different parameters / fix", draft an
   action order with the propose_action tool (action, stages, parameters,
   and a clear reason); it runs only after the user approves it in the dock.
   Parameter changes need evidence — ask first when unsure.
8. **The three platform-level blocks have different scopes**: `tool_runs` is
   the toolbox / extension-analysis timeline (deduplicated by tool; `metrics`
   appears only for tools with a stable contract, others are existence-only);
   `extensions` is existence-only for one-off jobs (LOGAN trace, NCBI
   submission prep); `prescreen` is the **virus pre-screening** project
   counters (`n_runs`/`n_pass`/`n_ok`) and is NOT about the selected sample —
   do not conflate them. To read a specific task's table use
   `read_tool_run_file`; for tables inside a sample's result dir use
   `read_result_table`. **Never guess a number that is not in the context** —
   say the data is unavailable instead.
9. **Finding data is a core platform capability**: when the user wants
   sequencing data of any species (transcriptome/RNA-seq, genomes, viromes…),
   call the `search_public_data` tool to search SRA+GSA directly. NEVER say
   "the platform cannot search" or push the user to external databases —
   that is exactly what the built-in public-data search module does. After
   starting, tell the user results will be broadcast automatically; when
   results arrive, summarise them and propose downloading selected runs via
   propose_action(download_runs).
10. **Self-description must be complete**: when asked "what can you do",
   present four sections — ① interpretation (samples/toolbox/prescreen
   results); ② **data discovery** (search_public_data over SRA+GSA with
   auto-broadcast); ③ **executable actions** via propose_action (pipeline
   stage re-runs, same-param retry, public-data download, **any of the ~30
   toolbox tools** — alignment/SDT/RDP/primers/consensus/ORF/
   phylogeography/dsRNA/miRNA, see list_platform_tools — and sample
   creation), all running for real after user approval; ④ boundary (you
   only draft action orders; approval belongs to the user). Never undersell
   yourself as interpretation-only — platform modules (map below) are
   drivable through action orders.

11. **MCP execution layer (you CAN do real work now — no action order
   needed for these)**: you have mcp_tool_catalog (look up the 19 vp_* MCP
   tools and their parameter contracts) and mcp_call (execute any vp_*
   tool). Coverage: toolbox execution (vp_run_tool), pipeline
   (vp_pipeline_start), CLI (vp_cli_run / vp_cli_submit: DB builds /
   public-data / NCBI submission / LOGAN trace / ICTV / virus
   pre-screening), task & artifact queries (vp_get_task /
   vp_read_run_file / vp_list_*), OpenVirome cloud search
   (vp_ov_search_runs / counts / mwas). When the user says "run / check /
   submit / find data for me": look up the contract with
   mcp_tool_catalog, then execute DIRECTLY via mcp_call (no action
   order), and report the envelope's summary/interpretation/next. For
   task-submitting calls, tell the user the task_id; only poll when they
   ask for progress — do not idle-wait inside the conversation.
   Exceptions: vp_delete_tool_run is disabled; for irreversible ops
   (db-migrate apply, tool-runs clean) confirm with the user first;
   propose_action remains for native flows (sample creation, download
   batches).

Platform capability map (the platform ships all of the following modules —
use this to answer "can the platform do X"):
{capabilities}

Platform context (auto-injected; all numbers computed from result files):
```json
{context}
```"""


def build_context(sample=None):
    """组装注入系统提示的上下文 dict（控制体积：样品 digest 已天然有界）。"""
    ctx = {'platform': digest.platform_brief()}
    if sample:
        d = digest.sample_digest(sample)
        if d is not None:
            ctx['selected_sample'] = d
    return ctx


def _clip_history(messages):
    out = []
    for m in messages[-MAX_HISTORY:]:
        role = m.get('role')
        if role not in ('user', 'assistant'):
            continue
        content = str(m.get('content') or '')
        out.append({'role': role, 'content': content[:MAX_MSG_CHARS]})
    return out


def _sse(obj):
    return f"data: {json.dumps(obj, ensure_ascii=False)}\n\n"


def _capability_map():
    """平台能力清单：从导航单一真源（web/pages.py 的 NAV_GROUPS）生成，
    供系统提示回答"平台能不能做 X"。这是修复"AI 说平台查不了转录组"的
    根治手段——能力感知来自平台自身的模块清单，而不是模型猜测。"""
    try:
        from Virus_Platform_Core.web.pages import NAV_GROUPS
        lines = []
        for g in NAV_GROUPS:
            items = []
            for it in g.get('items', [])[:8]:
                t = str(it.get('title') or '').strip()
                d = str(it.get('desc') or '').strip()[:36]
                items.append(f"{t}（{d}…）" if d else t)
            lines.append(f"- {g.get('label')}: {'；'.join(items)}")
        return '\n'.join(lines)
    except Exception as e:
        return f'（能力清单生成失败: {e}）'


def _system_prompt(sample):
    from Virus_Platform_Core.config import get_config
    lang = get_config().lang
    try:
        ctx = build_context(sample)
        ctx_json = json.dumps(ctx, ensure_ascii=False)
    except Exception as e:                     # 上下文失败不拦对话
        ctx_json = json.dumps({'error': f'context unavailable: {e}'},
                              ensure_ascii=False)
    if len(ctx_json) > 60000:                  # 极端大样品：截断保底
        ctx_json = ctx_json[:60000] + '…(truncated)'
    tpl = _SYSTEM_RULES_ZH if lang == 'zh' else _SYSTEM_RULES_EN
    return tpl.format(context=ctx_json, capabilities=_capability_map(),
                      lang_name='简体中文' if lang == 'zh' else 'English')


def run_chat(ai_cfg, messages, sample=None):
    """跑一轮对话，yield SSE 事件 dict（调用方逐个 `_sse()` 下发）。

    事件：delta（正文增量）/ tool（告知前端的工具活动）/ error（可读错误）/
    done（结束标记，附带本次是否降级过工具）。
    """
    sys_prompt = _system_prompt(sample)
    convo = [{'role': 'system', 'content': sys_prompt}]
    convo += _clip_history(messages)
    if not convo or convo[-1]['role'] != 'user':
        yield _sse({'type': 'error',
                    'message': '消息为空，请先输入问题'})
        yield _sse({'type': 'done'})
        return

    degraded = False
    use_tools = True
    for _round in range(MAX_TOOL_ROUNDS + 1):
        emit_parts = []
        tool_called = False
        tries = 0
        while True:
            tries += 1
            hit_tools = False
            try:
                for ev in chat_stream(
                        ai_cfg, convo,
                        tools=ai_tools.TOOLS_SPEC if use_tools else None):
                    if ev['type'] == 'delta':
                        emit_parts.append(ev['content'])
                        yield _sse({'type': 'delta', 'content': ev['content']})
                    elif ev['type'] == 'tool_calls':
                        hit_tools = True
                        tcs = ev['tool_calls']
                        # 协议要求工具调用轮的 assistant 消息带上已流出正文
                        # （content 用空串而非 null——部分国产端点不收 null）
                        convo.append({'role': 'assistant',
                                      'content': ''.join(emit_parts),
                                      'tool_calls': tcs})
                        for tc in tcs:
                            fn = (tc.get('function') or {})
                            name = fn.get('name') or ''
                            try:
                                fargs = json.loads(fn.get('arguments') or '{}')
                            except ValueError:
                                fargs = {}
                            yield _sse({'type': 'tool', 'name': name,
                                        'args': fargs})
                            result = ai_tools.dispatch(name, fargs)
                            convo.append({'role': 'tool',
                                          'tool_call_id': tc.get('id') or '',
                                          'content': result[:8000]})
                        break                # 工具结果已入列，进入下一轮
                    elif ev['type'] == 'message':
                        msg = ev.get('message') or {}
                        if msg and msg.get('tool_calls') is None:
                            text = msg.get('content')
                            if text and not emit_parts:
                                emit_parts.append(text)
                                yield _sse({'type': 'delta', 'content': text})
                            convo.append({'role': 'assistant',
                                          'content': text or ''})
                if hit_tools:
                    tool_called = True
                    break                    # → 下一轮（带工具结果）
                # 没遇到工具调用：正常结束
                if not emit_parts:
                    yield _sse({'type': 'error',
                                'message': '模型未返回内容，请稍后重试'})
                yield _sse({'type': 'done', 'tools_degraded': degraded})
                return
            except AIConfigError as e:
                yield _sse({'type': 'error',
                            'message': f'{e}（请在 设置 → AI 助手 中配置）'})
                yield _sse({'type': 'done'})
                return
            except AIProviderError as e:
                transient = (
                    e.status in (0, 429, 500, 502, 503, 504)
                    or 'connection' in str(e).lower()
                    or 'timed out' in str(e).lower())
                # 瞬态网络故障（连接被重置/限流/超时）且**尚未吐出任何字**
                # → 同一轮自动重试一次（已吐字重试会造成正文重复，只能报错）
                if not emit_parts and transient and tries < 2:
                    time.sleep(1.5)
                    continue
                # tools 参数不被支持是各家端点最常见差异：降级重试一次
                if use_tools and e.status in (400, 404, 422) and tries < 2:
                    use_tools = False
                    degraded = True
                    continue
                yield _sse({'type': 'error', 'message': str(e)})
                yield _sse({'type': 'done', 'tools_degraded': degraded})
                return
            except Exception as e:           # 兜底：不让生成器半途崩掉
                yield _sse({'type': 'error',
                            'message': f'{type(e).__name__}: {e}'})
                yield _sse({'type': 'done', 'tools_degraded': degraded})
                return
        if not tool_called:                  # 理论不可达，保险丝
            break

    yield _sse({'type': 'error',
                'message': '工具调用轮次超限，已截断本次回答'})
    yield _sse({'type': 'done', 'tools_degraded': degraded})


def chat_ready():
    from Virus_Platform_Core.config import get_config
    return ai_ready(get_config().ai)
