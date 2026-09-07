# CLAUDE.md

本文件为 Claude Code (claude.ai/code) 在本仓库中工作时提供指引。

## 项目概述

BLINK 是 HXMT HE（慧眼硬X射线调制望远镜，高能探测器）饱和分析与 TGF（地球伽马射线闪）探测工具包。用于处理卫星遥测数据，检测并重建探测器饱和（FIFO 缓冲区溢出）期间丢失的数据，以及搜索瞬变伽马射线事件。

## 构建与运行

Rust 工作区，使用 Cargo，resolver v3，edition 2024：

```bash
cargo build --release             # 构建所有 crate
cargo build -p blink          # 仅构建 CLI
cargo run -p blink -- <参数>  # 运行 CLI
```

Web 前端（Next.js + React 19 + TypeScript）位于 `web/`，使用 pnpm。

Python 脚本位于 `scripts/` 和 `analysis/`，用于绘图和数据分析（matplotlib、astropy、python-pptx）。使用 `uv` 管理 Python 依赖，运行脚本前需要激活 venv：

```bash
uv venv                          # 创建虚拟环境（首次）
source .venv/bin/activate         # 激活
uv pip install numpy matplotlib  # 安装依赖
python3 scripts/plot_compare.py ...  # 运行脚本
```

## CLI 用法

```
blink sat <COMMAND>

  # burst-centric (TRIGGER 必填; 1B 小时从 TRIGGER 自动派生)
  report <TRIGGER> --before <s> --after <s> -o <DIR>
                          # 完整诊断数据 pack: events_obs/rec/1k.csv + resets.csv + manifest.json
  detect <TRIGGER> --before <s> --after <s> [--box a|b|c]
                          # 列出 FIFO reset
  reconstruct <TRIGGER> --before <s> --after <s> [--box a|b|c] [--bin <s>]
                          # 跨 box 重建光变曲线 (1B 观测 + gap-fill)
  extract <TRIGGER> --before <s> --after <s> [--box a|b|c] [--source 1b|1k]
                          # 输出窗口内逐事件 (默认 1B 原始)
  compare <TRIGGER> --before <s> --after <s> [--box a|b|c]
                          # 1B vs 1K 互校 (1s 粗 bin, 0.1s 细 bin + 1ms cross-correlation)

  # 离线扫描 (无 TRIGGER)
  scan --epoch <YYYY-MM-DDTHH> [--box a|b|c]
                          # 扫描整小时 1B 找 FIFO reset

  # 低层诊断
  dump <SUB> --epoch <YYYY-MM-DDTHH> ...
    times | packets | events | hist | diag    # window-based, 同样接受 <TRIGGER> + --before/--after
    ptime | check-offset                       # 包索引范围 [pkt_min, pkt_max]

blink search <FROM> <TO>    # 在日期范围内扫描 TGF 候选信号
blink wwlln                 # 候选富集（挂 associated + coincidence_probability + train{neighbors_10min,is_train}，不筛选/不丢弃）
blink acd-audit <CSV> -o <CSV> [--scint csi|nai]  # 离线复算候选窗 ACD 符合计数（输入需 start/stop 列；csi=搜索同款 keep，nai=电子截止层正控制；需 1K 档案）
blink catalog [tgfs.json] -o <CSV>  # 目录生成：池级清洁（is_train 摘除）→ 论文判选 fa≤1e-5 ∪ (fa≤1 ∩ assoc)；逐步计数打印，无静默过滤
```

- TRIGGER 接 MET 数字或 UTC 字符串 (`2026-06-01T19:12:49.900`); CLI 内部转 MET 并把 1B 小时下取整。
- 跨小时窗口会 warning 但只加载 trigger 所在小时。
- 报告 pack 用法: 见 `scripts/plot_burst_report.py --pack <DIR> -o <PNG>`

## 架构

### 工作区结构 (`crates/`)

**Core（核心）** — 共享逻辑，不涉及具体仪器：
- `blink_core` — 核心类型、trait、FITS 文件读写（通过 `fitsio`）
- `blink_algorithms` — 统计算法（泊松分析、光变曲线、快照步进）
- `blink_solar` — 太阳几何计算
- `blink_lightning` — 闪电数据库查询（SQLite，通过 `rusqlite`）

**Instruments（仪器）** — 探测器专用实现：
- `blink_hxmt_he` — 主要焦点：HXMT HE 饱和检测/重建、1B/1K 数据读写
- `blink_fermi_gbm` — Fermi GBM 支持
- `blink_svom_grm` — SVOM GRM 支持

**Workflows（工作流）** — CLI 工具与处理管线：
- `blink` — 主入口，饱和分析命令
- `blink_search` — TGF 搜索算法
- `blink_filter` — TGF 过滤与闪电关联
- `blink_load` — 数据加载工具
- `blink_workflow` — 通用工作流工具

### 数据处理流程

1. **Solve** — 从 Level 1B 原始遥测包中提取事件 MET 时间
2. **Detect** — 识别饱和区间（FIFO 复位）
3. **Reconstruct** — 通过跨 Box 插值填充缺口
4. **Compare** — 与 Level 1K 管线输出进行互相关验证

### 关键概念

- **1B vs 1K**：Level 1B 是原始遥测数据；Level 1K 是标准管线产品。本工具从 1B 重建事件，恢复 1K 处理中因饱和丢失的数据。
- **Box A/B/C**：HXMT HE 的三个独立探测器箱体，各有独立的 FIFO 缓冲区。
- **FIFO 复位**：探测器缓冲区溢出时的硬件事件，标志饱和边界。
- **静默丢数（Silent drops）**：FIFO 满时 FPGA 静默丢弃事件，无硬件标记。已分析确认影响可忽略，不做检测。
- **MET**：Mission Elapsed Time（任务经过时间），HXMT 事件的主时间基准。

## Git 规则

- CLAUDE.md 不要提交到 git，也不要加入 .gitignore
- commit message 中不要出现 Claude、Anthropic、Co-Authored-By 等痕迹
- 不要在代码注释中提及 AI 辅助

## 主要依赖

- `fitsio` — FITS 科学数据格式读写
- `clap`（derive 模式）— CLI 参数解析
- `chrono` — 日期时间处理
- `statrs` — 统计函数
- `uom` — 类型安全的物理量单位
