@echo off
chcp 65001 >nul
title 拖拽代理（桌面 → 平台）
cd /d "%~dp0"
echo.
echo   ╔══════════════════════════════════════════╗
echo   ║   拖拽代理 · 桌面悬浮窗                    ║
echo   ╚══════════════════════════════════════════╝
echo.
echo   从资源管理器把文件拖到悬浮窗，真实路径直送平台，
echo   自动填入「最后点过的输入框」（零拷贝，不复制文件）。
echo   需先启动平台；依赖 tkinterdnd2（缺失时自动安装）。
echo.
rem 优先使用项目内置的绿色版 Python（3rd\python），不存在时回退系统 Python
set "PYTHONNOUSERSITE=1"
set "PYTHON=python"
if exist "%~dp03rd\python\python.exe" set "PYTHON=%~dp03rd\python\python.exe"
"%PYTHON%" -c "import tkinterdnd2" 2>nul || (
  echo   首次运行：安装 tkinterdnd2 …
  "%PYTHON%" -m pip install tkinterdnd2
)
"%PYTHON%" scripts\drag_proxy.py %*
echo.
echo   拖拽代理已退出。
pause
