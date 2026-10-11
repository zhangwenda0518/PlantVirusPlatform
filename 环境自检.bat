@echo off
chcp 65001 >nul
title 环境自检（模块 / 工具 / 数据库 / 磁盘）
cd /d "%~dp0"
rem 优先使用项目内置的绿色版 Python（3rd\python），不存在时回退系统 Python
set "PYTHONNOUSERSITE=1"
set "PYTHON=python"
if exist "%~dp03rd\python\python.exe" set "PYTHON=%~dp03rd\python\python.exe"
rem 依赖兜底（2026-10-11 起依赖包不随仓库分发，仅解释器随包）：缺包自动装一次
"%PYTHON%" -c "import flask, Bio, pandas, numpy, matplotlib, scipy, yaml" >nul 2>&1
if errorlevel 1 (
    echo [环境自检] 检测到缺少运行依赖，开始自动安装（requirements.lock，约几分钟，清华镜像）...
    "%PYTHON%" -m pip install -r requirements.lock -i https://pypi.tuna.tsinghua.edu.cn/simple
)
"%PYTHON%" main.py selfcheck
echo.
"%PYTHON%" Virus_Platform_Core\mem_check.py
pause
