<div align="center">

# OmicVerse.Skills

> *"我又不是不会分析，就是不知道用哪个函数，参数怎么填，跑完怎么检查，数据有问题怎么办……"*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://python.org)
[![PyPI](https://img.shields.io/pypi/v/omicverse-skills.svg)](https://pypi.org/project/omicverse-skills/)
[![Claude Code](https://img.shields.io/badge/Claude%20Code-Skill-blueviolet)](https://claude.ai/code)
[![Skills](https://img.shields.io/badge/Skills-43-green)](#skill-目录)

<br>

单细胞分析跑到一半，不知道 batch correction 该用 Harmony 还是 scVI？<br>
SCENIC 跑完 regulon 是空的，不知道是数据库没下载还是参数写错了？<br>
bulk DEG 做完想接 GSEA，却不知道 geneset 格式从哪里来？<br>
spatial 数据想做 deconvolution，三个方法不知道怎么选、怎么装？<br>

**把每一个分析流程固化成可调用的 Skill，让 AI Agent 真正接管你的组学分析！**

<br>

安装后配合 `ov.Agent` 或 Claude Code 即可使用<br>
43 个 Skill 覆盖单细胞、bulk RNA-seq、空间转录组、数据处理全流程<br>
每个 Skill 包含完整工作流、API 细节、防御性校验和真实报错处理

[Skill 目录](#skill-目录) · [安装](#安装) · [使用](#使用) · [贡献](#贡献)

</div>

---

## Skill 目录

### 单细胞分析

| Skill | 功能描述 |
|-------|---------|
| `single-cell-preprocessing` | QC、归一化、HVG、PCA、UMAP 完整预处理流程 |
| `single-cell-clustering-backends` | Leiden、Louvain、scICE、GMM 聚类后端选择 |
| `single-cell-annotation` | SCSA、MetaTiME、CellVote、GPTAnno、KNN 标注 |
| `single-cell-batch-integration` | Harmony、scVI、BBKNN、Combat 批次校正 |
| `single-cell-scenic` | RegDiffusion GRN、cisTarget regulon、AUCell 评分 |
| `single-cell-trajectory-inference` | PAGA、Palantir、VIA 轨迹推断 |
| `single-cell-rna-velocity` | scVelo、dynamo、latentvelo、graphvelo RNA 速率 |
| `single-cell-sctour-trajectory` | scTour 连续时间轨迹建模 |
| `single-cell-cellphonedb-communication` | CellPhoneDB v5 配体-受体通讯分析 |
| `single-cell-differential-expression` | Wilcoxon、t-test、memento-DE 单细胞差异表达 |
| `single-cell-differential-abundance` | scCODA、Milopy、Milo 细胞组成差异分析 |
| `single-cell-cytotrace2` | CytoTRACE2 发育潜能预测 |
| `single-cell-cnmf-program-discovery` | cNMF 基因表达程序发现 |
| `single-cell-lda-topic-clustering` | LDA 主题模型聚类 |
| `single-cell-kb-alignment` | kallisto/bustools 单细胞比对定量 |
| `single-cellfate-analysis` | CellFateGenie 伪时间基因发现 |
| `single-downstream-analysis` | AUCell、scDrug、NOCD 下游分析 |
| `single-multiomics` | MOFA、GLUE、SIMBA 多组学整合 |
| `single-popv-annotation` | PopV 10 算法共识细胞类型标注 |
| `single-to-spatial-mapping` | Single2Spatial scRNA 到空间转录组映射 |
| `cross-modal-celltype-transfer` | 跨模态细胞类型迁移 |
| `reference-label-transfer` | 参考图谱标签迁移 |

### Bulk RNA-seq

| Skill | 功能描述 |
|-------|---------|
| `bulk-deg-analysis` | 基因 ID 映射、DESeq2 归一化、差异表达、火山图 |
| `bulk-deseq2-analysis` | PyDESeq2 差异表达、fold-change 筛选、GSEA |
| `bulk-combat-correction` | pyComBat 批次效应校正 |
| `bulk-wgcna-analysis` | WGCNA 共表达网络、模块检测、hub 基因 |
| `bulk-stringdb-ppi` | STRING PPI 网络构建与可视化 |
| `bulk-to-single-deconvolution` | Bulk2Single beta-VAE 解卷积 |
| `bulk-trajblend-interpolation` | BulkTrajBlend 轨迹插值 |
| `gsea-enrichment` | GSEA 富集分析、基因集格式处理 |

### 空间转录组

| Skill | 功能描述 |
|-------|---------|
| `spatial-tutorials` | Visium/HD、Stereo-seq 预处理、Tangram、Starfysh 解卷积 |

### 数据处理与可视化

| Skill | 功能描述 |
|-------|---------|
| `data-io-loading` | h5ad、10x、Visium 数据读写 |
| `data-transform` | pandas/numpy 数据清洗与变换 |
| `data-stats-analysis` | scipy、statsmodels 统计检验 |
| `data-viz-plots` | matplotlib/seaborn 发表级图表 |
| `data-export-excel` | Excel 报告导出 |
| `data-export-pdf` | PDF 报告生成 |
| `datasets-loading` | OmicVerse 内置数据集加载 |
| `plotting-visualization` | 火山图、Venn、嵌入图、热图等 |
| `biocontext-knowledge` | UniProt、STRING、GO、PubMed 基因注释查询 |
| `bulk-fastq-quantification` | Bulk RNA-seq FASTQ→counts：SRA 下载、fastp、STAR+featureCounts 或 kb-python `technology='BULK'`，并衔接 `ov.bulk.pyDEG` |
| `single-cell-foundation-model` | scGPT、Geneformer、scFoundation、UCE、CellPLM —— 由 `ov.llm.SCLLMManager` 统一驱动 |
| `tcga-preprocessing` | TCGA bulk RNA-seq 预处理与生存分析 |

---

## 安装

```bash
pip install omicverse-skills
```

配合 OmicVerse Agent 使用时会自动加载；也可手动克隆到 `.claude/skills/`：

```bash
git clone https://github.com/Starlitnightly/omicverse-skills .claude/skills/omicverse-skills
```

---

## 使用

### Python API

```python
from omicverse_skills import skill_root, list_skills, load_skill_text

# 列出所有 skill
for skill in list_skills():
    print(skill["slug"], "—", skill["description"])

# 加载指定 skill 内容
text = load_skill_text("single-cell-preprocessing")
print(text)
```

### 配合 ov.Agent

```python
import omicverse as ov

agent = ov.Agent()
agent.run("对这个 AnnData 做单细胞预处理，然后跑 Leiden 聚类")
```

Agent 会自动通过 `skill_lookup` 匹配并加载对应 Skill 的执行指引。

---

## Skill 结构

每个 Skill 目录包含：

```
single-cell-preprocessing/
├── SKILL.md          # 执行指引（工作流、API、校验、分支选择）
└── references/       # 深度参考（函数签名、notebook 映射、兼容性说明）
```

`SKILL.md` 供 Agent 匹配与加载；`references/` 供 Agent 运行时按需读取。

---

## 贡献

欢迎提交新 Skill 或改进现有 Skill：

1. 在 `src/omicverse_skills/skills/` 下新建目录
2. 按照现有格式编写 `SKILL.md`（需含 `name`、`slug`、`description` frontmatter）
3. 提交 PR

---

## Star History

<a href="https://www.star-history.com/?repos=Starlitnightly%2Fomicverse-skills&type=date&legend=top-left">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/image?repos=Starlitnightly/omicverse-skills&type=date&theme=dark&legend=top-left" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/image?repos=Starlitnightly/omicverse-skills&type=date&legend=top-left" />
   <img alt="Star History Chart" src="https://api.star-history.com/image?repos=Starlitnightly/omicverse-skills&type=date&legend=top-left" />
 </picture>
</a>

---

<div align="center">

MIT License © [Starlitnightly](https://github.com/Starlitnightly)

</div>
