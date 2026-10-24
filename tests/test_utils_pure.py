# -*- coding: utf-8 -*-
"""utils 纯函数单测：路径校验 / 输出解码 / 步骤名白名单 / run_cmd 超时。

这些函数零外部依赖、是全平台安全与稳定性的地基（check_path 挡路径穿越，
run_cmd 管进程生命周期），却长期没有测试兜底——2026-09-13 全面体检后补上。
"""
import os
import sys

import pytest

from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.utils import (TaskLogger, check_path, check_step,
                                       decode_output, run_cmd)


# ── decode_output：GBK/UTF-8 双编码无损 ──────────────────────
def test_decode_output_utf8():
    raw = '中文路径'.encode('utf-8')
    assert decode_output(raw) == '中文路径'


def test_decode_output_gbk_fallback():
    raw = '中文路径'.encode('gbk')          # 非 UTF-8 字节串
    assert decode_output(raw) == '中文路径'


def test_decode_output_passthrough_str():
    assert decode_output('already str') == 'already str'


# ── check_path：穿越拒绝 / 平台内限定 / 存在性 ────────────────
def test_check_path_rejects_dotdot():
    with pytest.raises(ValueError):
        check_path(os.path.join('run', '..', 'etc'))


def test_check_path_rejects_dotdot_absolute():
    with pytest.raises(ValueError):
        check_path(os.path.join(str(PLATFORM_ROOT), 'run', '..', 'x'))


def test_check_path_resolves_relative_to_platform():
    p = check_path(os.path.join('run', 'x.txt'))
    assert os.path.isabs(p)
    assert p.startswith(str(PLATFORM_ROOT))


def test_check_path_in_platform_rejects_outside():
    outside = os.path.join(os.path.abspath(os.sep), 'windows', 'temp', 'x')
    with pytest.raises(ValueError):
        check_path(outside, in_platform=True)


def test_check_path_must_exist_missing():
    with pytest.raises(FileNotFoundError):
        check_path(os.path.join('run', 'no_such_file_xyz_12345'),
                   must_exist=True)


# ── check_step：步骤名白名单 ─────────────────────────────────
def test_check_step_accepts_plain():
    assert check_step('fastp') == 'fastp'
    assert check_step('02b-kvsuite') == '02b-kvsuite'


@pytest.mark.parametrize('bad', ['../evil', 'a b', 'a;b', 'a|b', '', 'x/y'])
def test_check_step_rejects_bad(bad):
    with pytest.raises(ValueError):
        check_step(bad)


# ── run_cmd：正常路径 + 默认超时看门狗 ────────────────────────
def test_run_cmd_success():
    rc = run_cmd([sys.executable, '-c', 'print("hello")'], silent=True)
    assert rc == 0


def test_run_cmd_failure_raises():
    with pytest.raises(RuntimeError):
        run_cmd([sys.executable, '-c', 'import sys; sys.exit(3)'], silent=True)


def test_run_cmd_timeout_kills_and_reports():
    """挂死命令被看门狗击杀，且报错明确说「超时」而非「退出码<0」。"""
    with pytest.raises(RuntimeError) as ei:
        run_cmd([sys.executable, '-c', 'import time; time.sleep(60)'],
                timeout=1, silent=True)
    assert '超时' in str(ei.value)


def test_run_cmd_timeout_zero_disables_watchdog():
    """timeout=0 显式关闭看门狗（快命令正常完成即视为关闭成功）。"""
    rc = run_cmd([sys.executable, '-c', 'print("ok")'], timeout=0, silent=True)
    assert rc == 0


# ── TaskLogger：logging.Logger 风格便捷级鸭子兼容 ─────────────
def test_tasklogger_logger_style_levels(tmp_path):
    """kvsuite 子模块按 logging.Logger 习惯调 logger.warning(...)，而流水线
    传入 TaskLogger——缺便捷级会 AttributeError 崩掉整个阶段（2026-10-13
    实测：salmon BAM 无 QUAL 告警分支每次必踩，②b 全线崩）。"""
    # TaskLogger 经 safe_open，只允许平台内路径 → 用 run/ 下临时目录
    tdir = os.path.join(PLATFORM_ROOT, 'run', '_pytest_tasklogger')
    os.makedirs(tdir, exist_ok=True)
    tlog = os.path.join(tdir, 't.log')
    lg = TaskLogger(tlog)
    lg.info('i')
    lg.warning('w')
    lg.error('e')
    lg.debug('d')                       # 静默，不抛即合规
    lg.close()
    with open(tlog, encoding='utf-8') as f:
        text = f.read()
    os.remove(tlog)
    try:
        os.rmdir(tdir)
    except OSError:
        pass
    assert '[INFO] i' in text and '[WARN] w' in text and '[ERROR] e' in text
    assert '[DEBUG]' not in text        # debug 不落盘
