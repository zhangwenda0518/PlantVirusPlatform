# MEGAHIT v1.2.9 for Windows（平台版，命令行，不内置 Python）

本目录是 [MEGAHIT](https://github.com/voutcn/megahit) v1.2.9 在 Windows 上从源码编译出的版本，
已放入 `3rd\tools\`，**不内置 Python**，直接复用平台的 `3rd\python`（Python 3.12.10）。

- 上游版本：v1.2.9（commit `c77d45c`）
- 编译器：MinGW-w64 g++ 16.2.0 (UCRT, posix threads, seh)
- 四个可执行文件全部**静态链接**，除 Windows 自带的 UCRT 外无任何外部 DLL 依赖
- 本目录**只有命令行入口**，供平台脚本调用
- 需要图形界面时用独立便携包：`D:\桌面\文章写作2\megahit-1.2.9-win64-with-python.zip`
  （内置 Python 3.12.10 + 图形界面，可整体拷到没装 Python 的机器上用）

## 入口

| 入口 | 用途 |
| --- | --- |
| `megahit.bat` | 命令行入口，适合批量、脚本化调用 |

## 目录结构

```
megahit-1.2.9-win64\
├─ megahit.bat                ← 命令行入口
├─ bin\
│  ├─ megahit.py              ← 官方 Python 驱动（已做 Windows 适配）
│  ├─ megahit_core.exe        ← 运行时自动选择（POPCNT+BMI2）
│  ├─ megahit_core_popcnt.exe
│  ├─ megahit_core_no_hw_accel.exe
│  └─ megahit_toolkit.exe     ← 辅助工具（统计、转换、取子图等）
├─ share\megahit\test_data\   ← 官方自检数据
├─ README_中文.md
└─ 编译与移植说明.md
```

## Python 是怎么找到的

`megahit.bat` 按下面顺序找解释器，用第一个存在的：

1. `<本目录>\python\python.exe`（如果你把 Python 放进了本目录）
2. `<本目录>\..\..\python\python.exe`（即平台 `3rd\python`）
3. PATH 里的 `python.exe`

都找不到时提示一行信息并退出（`exit /b 9009`），不会静默失败。

## 快速验证

```bat
cd /d D:\桌面\植物病毒分析平台\3rd\tools\megahit-1.2.9-win64
megahit.bat --version
megahit.bat --test
```

`--test` 通过时输出 `Test passed.`，并在当前目录生成 `megahit_test_out`（可删）。

## 常用命令

```bat
:: 双端
megahit.bat -1 reads_1.fq.gz -2 reads_2.fq.gz -o result -t 8

:: 交错双端
megahit.bat --12 interleaved.fq.gz -o result -t 8

:: 单端 / 长读长
megahit.bat -r reads.fq.gz -o result -t 8

:: 多组文库（不同 insert size）
megahit.bat -1 lib1_1.fq.gz -2 lib1_2.fq.gz -1 lib2_1.fq.gz -2 lib2_2.fq.gz -o result -t 8
```

常用参数：

| 参数 | 说明 |
| --- | --- |
| `-t/--num-cpu-threads` | 线程数，默认用满所有逻辑核 |
| `-m/--memory` | 内存上限：`0~1` 之间表示占总内存的比例（如 `-m 0.8`），`>=1` 表示字节数（如 `-m 32e9`） |
| `--min-contig-len` | 输出最短 contig，默认 200 |
| `--k-list` | 自定义 k 值，如 `--k-list 21,29,39,59,79,99,119` |
| `--min-count` | 最小 k-mer 计数，默认 2；单物种/高覆盖度数据可用 1 |
| `--no-mercy` | 不做低覆盖度剔除（更灵敏、更慢） |
| `--presets` | 只有 `meta-sensitive`、`meta-large` 两个；默认 meta 模式不加该参数 |
| `-r` / `--12` | 单端或长读长 / 交错双端 |

完整参数：`megahit.bat -h`

## 在 Python 里调用

平台的其他脚本可以直接调用，不用经过 `.bat`：

```python
import subprocess, os

TOOL = r"D:\桌面\植物病毒分析平台\3rd\tools\megahit-1.2.9-win64"
PY = r"D:\桌面\植物病毒分析平台\3rd\python\python.exe"

cmd = [PY, os.path.join(TOOL, "bin", "megahit.py"),
       "-1", r"D:\data\R1.fq.gz", "-2", r"D:\data\R2.fq.gz",
       "-o", r"D:\out\megahit", "-t", "8", "-m", "0.8"]
subprocess.run(cmd, cwd=TOOL, check=True)
```

## 输入格式

- FASTA / FASTQ：`.fa/.fna/.fq/.fasta/.fastq`
- **gzip（`.gz`）由核心程序直接读**，无需先解压
- **bzip2（`.bz2`）由 Python 驱动先解压**到输出目录 `tmp\` 下（Windows 无 `mkfifo`，故走临时文件）
- 路径支持中文和空格（已实测），但仍建议尽量用英文路径

## 注意事项

1. **不要直接跑 `bin\megahit_core.exe`**，它是底层子命令程序，请用 `megahit.bat` 或调用 `bin\megahit.py`。
2. 内存：复杂宏基因组占用高，建议显式 `-m`（如 `-m 0.8`），否则默认按可用内存自动计算。
3. 中间文件在输出目录 `tmp\`，成功后可直接删除；中断后重跑同一命令会**断点续跑**。
4. 杀毒软件首次运行可能扫描 `megahit_core.exe`，有几秒延迟，属正常。
5. 许可证：MEGAHIT 为 GPLv3，本目录同样遵循 GPLv3，源码见 https://github.com/voutcn/megahit 。

## 引用

Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 2015, 31(10): 1674-1676.

Li D, Luo R, Liu CM, et al. MEGAHIT v1.0: A fast and scalable metagenome assembler driven by advanced methodologies and community practices. *Methods*, 2016, 102: 3-11.
