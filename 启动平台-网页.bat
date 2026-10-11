@echo off
chcp 65001 >nul
title 植物病毒分析平台（网页模式）
cd /d "%~dp0"
echo.
echo   ╔══════════════════════════════════════════╗
echo   ║   植物病毒分析平台 · 网页模式（浏览器）     ║
echo   ╚══════════════════════════════════════════╝
echo.
echo   服务启动后自动在默认浏览器打开；地址也会打印在下方，
echo   可复制到其它浏览器/设备调试。关闭本窗口即退出平台。
echo.
rem 优先使用项目内置的绿色版 Python（3rd\python），不存在时回退系统 Python
set "PYTHONNOUSERSITE=1"
set "PYTHON=python"
if exist "%~dp03rd\python\python.exe" set "PYTHON=%~dp03rd\python\python.exe"
rem 依赖兜底（2026-10-11 起依赖包不随仓库分发，仅解释器随包）：缺包自动装一次
"%PYTHON%" -c "import flask, Bio, pandas, numpy, matplotlib, scipy, yaml" >nul 2>&1
if errorlevel 1 (
    echo   首次运行：正在安装运行依赖（requirements.lock 锁定版本，约几分钟，清华镜像）...
    "%PYTHON%" -m pip install -r requirements.lock -i https://pypi.tuna.tsinghua.edu.cn/simple
)
"%PYTHON%" app.py --web
pause
