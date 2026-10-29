@echo off
rem ============================================================
rem  MEGAHIT v1.2.9 for Windows -- launcher
rem
rem  Python lookup order:
rem    1) python\python.exe under this folder  (embedded portable)
rem    2) ..\..\python\python.exe              (platform 3rd\python)
rem    3) python.exe on PATH                   (system Python)
rem
rem  Keep this file ASCII-only: cmd.exe mis-parses UTF-8 batch files.
rem ============================================================
setlocal enableextensions
set "MEGAHIT_HOME=%~dp0"
set "PYEXE="

if exist "%MEGAHIT_HOME%python\python.exe" set "PYEXE=%MEGAHIT_HOME%python\python.exe"
if not defined PYEXE if exist "%MEGAHIT_HOME%..\..\python\python.exe" set "PYEXE=%MEGAHIT_HOME%..\..\python\python.exe"
if not defined PYEXE for /f "delims=" %%I in ('where python.exe 2^>nul') do if not defined PYEXE set "PYEXE=%%I"

if not defined PYEXE (
    echo [megahit] No Python interpreter found.
    echo   Option 1: install Python 3.6+ and add it to PATH
    echo   Option 2: unzip a portable Python into "%MEGAHIT_HOME%python\"
    echo   Option 3: use megahit-1.2.9-win64-with-python.zip instead
    exit /b 9009
)

"%PYEXE%" "%MEGAHIT_HOME%bin\megahit.py" %*
exit /b %ERRORLEVEL%
