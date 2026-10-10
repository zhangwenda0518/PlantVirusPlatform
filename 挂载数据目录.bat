@echo off
chcp 65001 >nul
title 挂载数据目录（源码平台 → 分发版平台）
rem 将源码平台的测序数据目录挂载进分发版平台（NTFS 目录联接，不复制文件）。
rem 用法：双击本文件；重复运行安全（已挂载则跳过）。
rem 约定：分发包与源码平台同级，即 ..\PVAP\VirusPlatform；
rem       若放别处，请编辑下方 PKG 一行指向实际的 VirusPlatform 目录。
setlocal
set "SRC=%~dp0run\fastq"
set "PKG=%~dp0..\PVAP\VirusPlatform"

if not exist "%PKG%\VirusPlatform.exe" (
  echo 未找到分发版平台: %PKG%
  echo 请确认 PVAP\VirusPlatform 与源码平台在同一个父目录下，
  echo 或编辑本文件把 PKG 改成实际路径后重试。
  echo.
  pause
  exit /b 1
)

set "DST=%PKG%\run\fastq"
if exist "%DST%" (
  echo 已存在: %DST%
  echo （若这里面看不到你的测序数据，说明之前挂到了别处，可删除该目录后重跑）
) else (
  mklink /J "%DST%" "%SRC%" >nul
  if errorlevel 1 (
    echo 挂载失败（需要 NTFS 文件系统；确认 %PKG%\run 目录存在）
    echo.
    pause
    exit /b 1
  )
  echo 已挂载: %DST%  -^>  %SRC%
)

echo.
echo 之后在分发版平台「分析管道」新建样品时，点 ■ 进入 fastq 目录
echo 即可直接选用源码平台里的测序数据（不占用额外磁盘空间）。
pause
