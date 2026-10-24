# 第三方组件清单（THIRD PARTY NOTICES）

本文件是**索引**而非法律意见：各组件的权威许可证以其目录内随附的
LICENSE/COPYING 与上游声明为准；Python 包以安装元数据（dist-info
METADATA）的 License 字段为准。重新生成：
`python dev_tools/third_party_notices.py`。

## 一、随包分发的第三方工具（33 项，3rd/bin + 3rd/tools）

| 组件（随包路径） | 目录内随附的许可/说明文件 |
|---|---|
| 3rd/bin/aria2c.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/bin/crabz.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/bin/kun_peng-v0.7.12-x86_64-pc-windows-msvc.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/bin/minibwa.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/bin/muscle-win64.v5.3.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/bin/seqkit.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/bin/sracha.exe | 3rd/bin/（单文件工具，许可见上游 Release 说明） |
| 3rd/tools/Blast/ | 3rd/tools/Blast/doc/README.txt |
| 3rd/tools/FastTree/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/SPAdes-Windows/ | 3rd/tools/SPAdes-Windows/README-WINDOWS.txt; 3rd/tools/SPAdes-Windows/python/LICENSE.txt; 3rd/tools/SPAdes-Windows/share/spades/LICENSE; 3rd/tools/SPAdes-Windows/share/spades/README.md |
| 3rd/tools/bcftools/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/datasets/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/diamond/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/dsrnamax/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/fastp/ | 3rd/tools/fastp/LICENSE; 3rd/tools/fastp/README.md |
| 3rd/tools/iQtree/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/jre-snpeff/ | 3rd/tools/jre-snpeff/legal/java.base/LICENSE; 3rd/tools/jre-snpeff/legal/java.datatransfer/LICENSE; 3rd/tools/jre-snpeff/legal/java.desktop/LICENSE; 3rd/tools/jre-snpeff/legal/java.logging/LICENSE |
| 3rd/tools/lsd2/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/mafft-win/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/minimap2/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/mmseqs/ | 3rd/tools/mmseqs/LICENSE.md; 3rd/tools/mmseqs/README.md |
| 3rd/tools/pandepth/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/raxml-ng/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/rdp5/ | 3rd/tools/rdp5/README_3seqTable.txt |
| 3rd/tools/salmon2/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/samtools/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/snpeff/ | 3rd/tools/snpeff/snpEff/LICENSE.md |
| 3rd/tools/snpgenie/ | 3rd/tools/snpgenie/LICENSE.txt; 3rd/tools/snpgenie/README.md |
| 3rd/tools/srna-target/ | 3rd/tools/srna-target/LICENSES/NOTICE.txt; 3rd/tools/srna-target/miranda/README.txt; 3rd/tools/srna-target/pita/README.txt; 3rd/tools/srna-target/pita/perl/README.txt |
| 3rd/tools/strawberry-perl/ | 3rd/tools/strawberry-perl/README.txt; 3rd/tools/strawberry-perl/licenses/License.rtf; 3rd/tools/strawberry-perl/licenses/fribidi/AUTHORS; 3rd/tools/strawberry-perl/licenses/fribidi/COPYING |
| 3rd/tools/table2asn/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/trimAl/ | （目录内未随附许可文件，见上游项目） |
| 3rd/tools/viral_consensus/ | （目录内未随附许可文件，见上游项目） |

例外：`3rd/tools/rdp5/3seqTable`（105MB，超 GitHub 单文件限制）经
Release（tag: external-tools）分发；RDP5 会话运行态文件（LastSave.rdp5 / tempfile）
不随版本库与插件包分发。

## 二、内置 Python 运行时的第三方包（129 项，3rd/python）

完整版本快照见 `requirements.lock`；下表为包名 / 版本 / 元数据声明的许可证。

| 包 | 版本 | 许可证（元数据） |
|---|---|---|
| PySocks | 1.7.1 | BSD |
| adjustText | 1.4.0 | MIT |
| altgraph | 0.17.5 | MIT |
| annotated-types | 0.8.0 | MIT |
| anyio | 4.15.1 | MIT |
| archspec | 0.2.6 | Apache-2.0 OR MIT |
| attrs | 26.1.0 | MIT |
| bcbio-gff | 0.7.1 | Biopython License |
| biopython | 1.88 | LicenseRef-Biopython-License-Agreement |
| blinker | 1.9.0 | （元数据未声明） |
| bottle | 0.13.4 | MIT |
| certifi | 2026.7.22 | MPL-2.0 |
| cffi | 2.1.1 | MIT-0 |
| charset-normalizer | 3.5.1 | MIT |
| choreographer | 1.3.0 | # MIT License |
| click | 8.5.0 | BSD-3-Clause |
| clr_loader | 0.3.1 | （元数据未声明） |
| colorama | 0.4.6 | （元数据未声明） |
| configs | 3.0.3 | MIT |
| contourpy | 1.3.3 | BSD 3-Clause License |
| cryptography | 50.0.2 | Apache-2.0 OR BSD-3-Clause |
| cycler | 0.12.1 | Copyright (c) 2015, matplotlib project |
| Cython | 3.3.0 | Apache-2.0 |
| defusedxml | 0.7.1 | PSFL |
| distinctipy | 1.3.4 | MIT License |
| dna_features_viewer | 3.1.5 | MIT |
| et_xmlfile | 2.0.0 | MIT |
| Flask | 3.1.3 | BSD-3-Clause |
| fonttools | 4.65.0 | MIT |
| gbdraw | 0.14.0b0 | MIT |
| git-filter-repo | 2.47.0 | MIT |
| graphviz | 0.21 | MIT |
| greenlet | 3.5.6 | MIT AND PSF-2.0 |
| h11 | 0.16.0 | MIT |
| httpcore | 1.0.9 | BSD-3-Clause |
| httpcore2 | 2.12.0 | BSD-3-Clause |
| httpx | 0.28.1 | BSD-3-Clause |
| httpx2 | 2.12.0 | BSD-3-Clause |
| httpx-sse | 0.4.3 | MIT |
| idna | 3.19 | BSD-3-Clause |
| iniconfig | 2.3.0 | MIT |
| itsdangerous | 2.2.0 | （元数据未声明） |
| Jinja2 | 3.1.6 | （元数据未声明） |
| jiter | 0.16.0 | MIT |
| jsonschema | 4.26.0 | MIT |
| jsonschema-specifications | 2025.9.1 | MIT |
| kaleido | 1.4.0 | The MIT License (MIT) |
| kiwisolver | 1.5.1 | ========================= |
| logistro | 2.0.1 | MIT License |
| MarkupSafe | 3.0.3 | BSD-3-Clause |
| matplotlib | 3.11.1 | License agreement for matplotlib versions 1.3.0 and later |
| mcp | 1.30.0 | MIT |
| mcp-types | 2.3.0 | MIT |
| narwhals | 2.26.0 | MIT |
| numpy | 2.5.3 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |
| openai | 3.13.0 | Apache-2.0 |
| openpyxl | 3.1.5 | MIT |
| opentelemetry-api | 1.45.1 | Apache-2.0 |
| orfipy | 0.0.4 | （元数据未声明） |
| orjson | 3.12.0 | MPL-2.0 AND (Apache-2.0 OR MIT) |
| outcome | 1.3.0.post0 | MIT OR Apache-2.0 |
| packaging | 26.3 | Apache-2.0 OR BSD-2-Clause |
| pandas | 3.0.5 | BSD 3-Clause License |
| pefile | 2024.8.26 | MIT |
| phylo-treetime | 0.12.1 | MIT |
| pillow | 12.3.0 | MIT-CMU |
| pip | 26.2.1 | MIT |
| platformdirs | 4.11.8 | MIT |
| playwright | 1.62.0 | Apache-2.0 |
| plotly | 7.0.0 | MIT |
| pluggy | 1.6.0 | MIT |
| primer3-py | 2.2.0 | GPLv2 |
| progress | 1.6.1 | ISC |
| proxy_tools | 0.1.0 | MIT |
| psutil | 7.2.2 | BSD-3-Clause |
| pyahocorasick | 2.3.1 | BSD-3-Clause and Public-Domain |
| pyCirclize | 1.10.1 | MIT |
| pycparser | 3.0 | BSD-3-Clause |
| pydantic | 2.13.5 | MIT |
| pydantic_core | 2.46.5 | MIT |
| pydantic-settings | 2.16.0 | MIT |
| pyee | 13.0.1 | MIT |
| pyfastx | 2.3.1 | MIT |
| Pygments | 2.21.0 | BSD-2-Clause |
| pyhmmer | 0.12.3 | MIT AND BSD-3-Clause AND BSD-2-Clause |
| pyinstaller | 6.22.2 | GPLv2-or-later with a special exception which allows to use PyInstaller to build and distribute n… |
| pyinstaller-hooks-contrib | 2026.7 | （元数据未声明） |
| PyJWT | 2.15.1 | MIT |
| pyparsing | 3.3.2 | MIT |
| pyrodigal | 3.7.1 | GPL-3.0-or-later |
| pyrodigal_rv | 0.1.0 | GPL-3.0-only |
| pytest | 8.4.2 | MIT |
| python-dateutil | 2.9.0.post0 | Dual License |
| python-dotenv | 1.2.4 | BSD-3-Clause |
| python-multipart | 0.0.32 | Apache-2.0 |
| pythonnet | 3.1.0 | MIT |
| pywebview | 6.2.1 | BSD 3-Clause License |
| pywin32 | 312 | PSF |
| pywin32-ctypes | 0.2.3 | BSD-3-Clause |
| PyYAML | 6.0.3 | MIT |
| referencing | 0.37.0 | MIT |
| requests | 2.34.2 | Apache-2.0 |
| rpds-py | 2026.9.1 | MIT |
| scipy | 1.18.1 | Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers. |
| seaborn | 0.13.2 | （元数据未声明） |
| selenium | 4.49.0 | Apache-2.0 |
| setuptools | 84.0.0 | MIT |
| simplejson | 4.1.2 | MIT OR AFL-2.1 |
| six | 1.17.0 | MIT |
| sniffio | 1.3.1 | MIT OR Apache-2.0 |
| sortedcontainers | 2.4.0 | Apache 2.0 |
| sse-starlette | 3.5.0 | BSD-3-Clause |
| starlette | 1.8.0 | BSD-3-Clause |
| svgwrite | 1.4.3 | MIT License |
| taxburst | 0.3.2 | （元数据未声明） |
| tqdm | 4.70.1 | MPL-2.0 AND MIT |
| trio | 0.34.0 | MIT OR Apache-2.0 |
| trio-websocket | 0.12.2 | （元数据未声明） |
| truststore | 0.10.4 | MIT |
| typing_extensions | 4.16.0 | PSF-2.0 |
| typing-inspection | 0.4.4 | MIT |
| tzdata | 2026.3 | Apache-2.0 |
| urllib3 | 2.7.0 | MIT |
| uvicorn | 0.54.0 | BSD-3-Clause |
| ViennaRNA | 2.7.2 | Disclaimer and Copyright |
| websocket-client | 1.9.2 | Apache-2.0 |
| Werkzeug | 3.1.8 | BSD-3-Clause |
| wsproto | 1.3.2 | MIT |
| xmltodict | 1.0.4 | MIT |

## 三、分发口径备注

- 平台自身代码随仓库分发；第三方工具本体**克隆即用**（3rd/ 全量纳管，
  例外见上），或经「设置 → 插件管理」从插件源按需安装（32 个工具包，
  逐包 sha256 校验）。
- GPL 系组件（RDP5、SPAdes、snpEff、snipit、svguest 等）以**独立进程/  独立文件**形式随包，不与其源码链接；其许可证文本随各自目录分发。
- 本清单由构建脚本扫描生成，**不构成法律意见**；商用再分发前请按
  各组件许可证自行核对（尤其 GPL 系的传染条款与商标条款）。
