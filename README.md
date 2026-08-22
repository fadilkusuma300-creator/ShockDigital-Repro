# 🏢 SME 数字化转型渠道分析与企业复苏

> 三波纵贯研究中，冲击条件性数字化渠道对中小企业（SME）复苏的影响效应分析。

<input type="radio" name="lang" id="zh" checked>
<input type="radio" name="lang" id="en">

<div class="tabs">
  <label for="zh" class="tab">🇨🇳 简体中文</label><label for="en" class="tab">🇬🇧 English</label>
</div>

<div class="zh">

# 📖 项目简介

本仓库是一套完整的因果分析流水线代码，用于研究**冲击条件性数字化渠道**如何影响中小企业（SME）的复苏。流水线完成以下工作：

- 🔗 构建纵向三波样本（疫情前基线 + 两轮追踪）
- ⚖️ 估计双重稳健（Doubly Robust）处理效应
- 📈 刻画初始销售冲击强度上的连续效应异质性
- 🛡️ 执行稳健性分析、安慰剂检验与阴性对照
- 🧪 进行半合成数据验证
- 📊 生成论文级图表

---

# 🗂️ 数据来源

World Bank Enterprise Surveys（世界银行企业调查）微观数据为第三方数据，**不随本仓库分发**。

官方来源：

- 🌐 企业调查数据门户：https://www.enterprisesurveys.org/en/data
- 🦠 COVID-19 企业调查与追踪调查信息：https://www.enterprisesurveys.org/en/covid-19

注册用户可免费下载企业层面数据与调查文档。请将相关数据放入：

```text
data/raw/baseline/    # 疫情前基线调查文件
data/raw/followup/    # COVID-19 追踪调查文件
```

支持格式：`.dta`、`.sav`、`.csv`、`.txt`、`.rds`。原始变量名通过 `config/variable_map.yml` 统一映射；若调查版本使用不同的原始列名，请在该配置文件中添加对应的候选列表。

---

# 🧱 项目结构

```text
📁 R/                         分析函数（编号分模块）
📁 scripts/                   分阶段入口脚本
📄 config/analysis.yml        分析设置与超参数
📄 config/variable_map.yml    世界银行变量映射
📄 config/figures.yml         论文图表设置
🐍 python/causalpfn_bridge.py CausalPFN 验证桥接
🐍 python/fig1_framework.py   研究框架示意图
▶️  run_all.R                 端到端流水线
```

生成文件写入 `data/derived/` 与 `results/`，这两个目录已被版本控制排除。

---

# ⚙️ 安装依赖

### 1️⃣ R 依赖

```bash
Rscript install.R
```

### 2️⃣ Python 依赖（研究框架图）

```bash
python -m pip install -r requirements-figures.txt
```

> 🪶 轻量：仅包含 `matplotlib`、`PyYAML`。

### 3️⃣ Python 依赖（CausalPFN 验证，需单独环境）

```bash
python -m pip install -r requirements-causalpfn.txt
```

> 🐻 较重：包含 `torch`、`causalpfn` 等，建议使用独立 Python 环境，避免只为一张图而安装几 GB 的深度学习栈。

---

# 🚀 运行分析

从项目根目录执行：

```bash
Rscript run_all.R
```

流水线包含 6 个阶段：

1. 🗃️ 数据清洗与相邻三波样本构建
2. 💰 主要销售复苏分析
3. 📊 次要复苏结果分析
4. 🛡️ 稳健性、安慰剂与阴性对照分析
5. 🧪 半合成验证
6. 📈 论文图表生成

也可单独运行各阶段：

```bash
Rscript scripts/01_prepare_data.R
Rscript scripts/02_primary_sales.R
Rscript scripts/03_secondary_outcomes.R
Rscript scripts/04_robustness.R
Rscript scripts/05_validation.R
Rscript scripts/06_figures_tables.R
```

---

# 📦 主要输出

关键输出写入 `results/`：

- 📋 `results/tables/`：样本流转、效应估计、稳健性结果及图表数值输入
- 🖼️ `results/figures/Fig1.pdf` ~ `Fig4.pdf`：论文配图
- 🧪 `results/validation/`：半合成验证结果
- 🔍 `results/diagnostics/`：平衡性、样本流转与数据质量摘要

图表外观由 `config/figures.yml` 统一控制，修改呈现样式无需改动估计代码。

</div>

<div class="en">

# 📖 About

This repository provides a complete causal-analysis pipeline to study how **shock-conditional digital channels** affect SME recovery. The pipeline:

- 🔗 Constructs an adjacent three-wave longitudinal sample (pre-pandemic baseline + two follow-ups)
- ⚖️ Estimates doubly robust treatment effects
- 📈 Traces continuous effect heterogeneity over initial sales-shock intensity
- 🛡️ Runs robustness, placebo, and negative-control analyses
- 🧪 Performs semi-synthetic validation
- 📊 Generates publication-ready figures and tables

---

# 🗂️ Data

The World Bank Enterprise Surveys microdata are third-party data and are **not distributed** with this repository.

Official sources:

- 🌐 Enterprise Surveys data portal: https://www.enterprisesurveys.org/en/data
- 🦠 COVID-19 Enterprise Surveys and follow-up survey information: https://www.enterprisesurveys.org/en/covid-19

Registered users can download firm-level data and documentation for free. Place files in:

```text
data/raw/baseline/    # pre-pandemic baseline files
data/raw/followup/    # COVID-19 follow-up files
```

Supported formats: `.dta`, `.sav`, `.csv`, `.txt`, `.rds`. Source-variable names are harmonized through `config/variable_map.yml`; if a survey release uses a different raw column name, add it to the relevant candidate list there.

---

# 🧱 Project structure

```text
📁 R/                         analysis functions
📁 scripts/                   stage-by-stage entry points
📄 config/analysis.yml        analysis settings and hyperparameters
📄 config/variable_map.yml    World Bank variable mappings
📄 config/figures.yml         publication figure settings
🐍 python/causalpfn_bridge.py CausalPFN validation bridge
🐍 python/fig1_framework.py   study framework figure
▶️  run_all.R                  end-to-end pipeline
```

Generated files are written to `data/derived/` and `results/`; both are excluded from version control.

---

# ⚙️ Installation

### 1️⃣ R dependencies

```bash
Rscript install.R
```

### 2️⃣ Python dependencies (study framework figure)

```bash
python -m pip install -r requirements-figures.txt
```

> 🪶 Lightweight: `matplotlib`, `PyYAML` only.

### 3️⃣ Python dependencies (CausalPFN validation — separate environment)

```bash
python -m pip install -r requirements-causalpfn.txt
```

> 🐻 Heavier: `torch`, `causalpfn`, etc. Use a dedicated Python environment so you don't install a multi-GB deep-learning stack just to render one figure.

---

# 🚀 Run the analysis

From the project root:

```bash
Rscript run_all.R
```

The pipeline runs 6 stages:

1. 🗃️ Data harmonization and adjacent three-wave sample construction
2. 💰 Primary sales-recovery analysis
3. 📊 Secondary recovery outcomes
4. 🛡️ Robustness, placebo, and negative-control analyses
5. 🧪 Semi-synthetic validation
6. 📈 Publication figures and tables

Individual stages can also be run directly:

```bash
Rscript scripts/01_prepare_data.R
Rscript scripts/02_primary_sales.R
Rscript scripts/03_secondary_outcomes.R
Rscript scripts/04_robustness.R
Rscript scripts/05_validation.R
Rscript scripts/06_figures_tables.R
```

---

# 📦 Main outputs

Key outputs are written under `results/`:

- 📋 `results/tables/`: sample flow, effect estimates, robustness results, and numerical inputs to figures
- 🖼️ `results/figures/Fig1.pdf` ~ `Fig4.pdf`: publication figures
- 🧪 `results/validation/`: semi-synthetic validation results
- 🔍 `results/diagnostics/`: balance, sample-flow, and data-quality summaries

Figure appearance is controlled centrally through `config/figures.yml`, so presentation changes require no changes to the estimation code.

</div>

<style>
  input#zh, input#en { display: none; }
  .tabs { border-bottom: 2px solid #d0d7de; margin: 16px 0; }
  .tab { display: inline-block; padding: 10px 24px; margin-bottom: -2px; cursor: pointer; font-weight: 600; color: #57606a; border-bottom: 2px solid transparent; user-select: none; }
  .tab:hover { color: #24292f; }
  .zh, .en { display: none; }
  #zh:checked ~ .tabs .tab[for="zh"],
  #en:checked ~ .tabs .tab[for="en"] { color: #1f4e79; border-bottom-color: #1f4e79; }
  #zh:checked ~ .zh { display: block; }
  #en:checked ~ .en { display: block; }
</style>
