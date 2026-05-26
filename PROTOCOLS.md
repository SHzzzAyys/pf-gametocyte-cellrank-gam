# 单细胞测序分析 — 模型构建 Protocol 资源汇总

搜集自 GitHub，按分析阶段分类，聚焦**轨迹推断、RNA velocity、命运概率、基因调控网络**等模型构建方向。

---

## 一、综合最佳实践 / 教程

| 仓库 | 简介 | 链接 |
|------|------|------|
| **theislab/single-cell-best-practices** ★1198 | Theis Lab 官方单细胞分析最佳实践书籍（配套代码），覆盖预处理→聚类→轨迹→多组学全流程 | https://github.com/theislab/single-cell-best-practices |
| **scRNA-tools/scRNA-tools** ★339 | 单细胞 RNA-seq 工具数据库，按功能分类索引 500+ 工具 | https://github.com/scRNA-tools/scRNA-tools |
| **kpatel427/YouTubeTutorials** ★378 | 配套 YouTube 教程的 R/Python 单细胞分析代码合集（Seurat / Scanpy / DESeq2 等）| https://github.com/kpatel427/YouTubeTutorials |
| **agitter/single-cell-pseudotime** ★439 | 伪时间算法大汇总：列出并评比 50+ 种伪时间/轨迹推断方法，附论文与代码链接 | https://github.com/agitter/single-cell-pseudotime |

---

## 二、上游处理 / 定量

| 仓库 | 简介 | 链接 |
|------|------|------|
| **nf-core/scrnaseq** ★328 | Nextflow 单细胞 RNA-seq 上游流程（10x / Smart-seq），支持 STARsolo / Alevin / kb-python | https://github.com/nf-core/scrnaseq |
| **COMBINE-lab/salmon** ★885 | 转录本定量工具（selective alignment），支持 alevin 单细胞模式 | https://github.com/COMBINE-lab/salmon |
| **COMBINE-lab/alevin-fry** ★207 | alevin 的高速继任者，单细胞 UMI 处理，内存占用极低 | https://github.com/COMBINE-lab/alevin-fry |
| **pachterlab/kallistobustools** ★121 | kallisto \| bustools 单细胞预处理流程（kb-python），速度极快 | https://github.com/pachterlab/kallistobustools |
| **sdparekh/zUMIs** ★294 | 支持 UMI 的单细胞 RNA-seq 流程，处理 SMART-seq/Drop-seq/10x | https://github.com/sdparekh/zUMIs |

---

## 三、核心分析框架

| 仓库 | 简介 | 链接 |
|------|------|------|
| **scverse/scanpy** ★2468 | Python 单细胞分析标准库（预处理/聚类/可视化），可扩展到亿级细胞 | https://github.com/scverse/scanpy |
| **scverse/anndata** ★745 | 单细胞数据格式标准（h5ad），Scanpy/scVelo/CellRank 的数据基础 | https://github.com/scverse/anndata |
| **scverse/muon** ★268 | 多模态组学框架（RNA + ATAC + Protein），扩展自 AnnData | https://github.com/scverse/muon |

---

## 四、RNA Velocity 模型

| 仓库 | 简介 | 链接 |
|------|------|------|
| **velocyto-team/velocyto.py** ★176 | RNA velocity 原始实现（La Manno 2018），从 BAM 文件生成 spliced/unspliced loom | https://github.com/velocyto-team/velocyto.py |
| **theislab/scvelo** | scVelo：随机/动力学 RNA velocity（Bergen 2020），`recover_dynamics` + `latent_time`，本项目核心依赖 | https://github.com/theislab/scvelo |
| **basilkhuder/Seurat-to-RNA-Velocity** ★142 | Seurat 对象接入 RNA velocity 的完整指南（R + Python 混合流程）| https://github.com/basilkhuder/Seurat-to-RNA-Velocity |
| **GuangyuWangLab2021/cellDancer** ★72 | 深度学习 RNA velocity（逐细胞动力学参数估计，适合异质性强的数据）| https://github.com/GuangyuWangLab2021/cellDancer |
| **aristoteleo/dynamo-release** ★498 | Dynamo：向量场重建 + 微分几何分析，支持代谢标记 scRNA-seq，命运预测最完整 | https://github.com/aristoteleo/dynamo-release |

---

## 五、细胞命运 / 轨迹推断

| 仓库 | 简介 | 链接 |
|------|------|------|
| **theislab/cellrank** | **CellRank 1.5.x**（本项目版本）：VelocityKernel → GPCCA → 吸收概率，API 见 `cellrank.tl.*` | https://github.com/theislab/cellrank |
| **scverse/cellrank_notebooks** ★8 | CellRank 官方教程 notebooks（多种 kernel、多种生物场景）| https://github.com/scverse/cellrank_notebooks |
| **broadinstitute/wot** ★162 | Waddington OT：最优传输推断跨时间点的细胞命运转变（时间序列数据专用）| https://github.com/broadinstitute/wot |
| **WPZgithub/CEFCON** ★31 | 从 scRNA-seq 解码细胞命运决策的驱动调控因子（GRN + 命运网络）| https://github.com/WPZgithub/CEFCON |

---

## 六、基因调控网络（GRN）

| 仓库 | 简介 | 链接 |
|------|------|------|
| **aertslab/pySCENIC** ★606 | SCENIC：单细胞转录因子调控网络推断（TF motif + co-expression），本论文也用了该方法 | https://github.com/aertslab/pySCENIC |
| **smorabit/hdWGCNA** ★475 | 高维 WGCNA：单细胞/空间转录组共表达网络模块分析 | https://github.com/smorabit/hdWGCNA |
| **aertslab/SCope** ★79 | SCENIC 结果的交互可视化工具（配套 pySCENIC 使用）| https://github.com/aertslab/SCope |

---

## 七、细胞类型注释

| 仓库 | 简介 | 链接 |
|------|------|------|
| **ZJUFanLab/scCATCH** ★241 | 基于 marker gene 字典的自动细胞类型注释（支持 100+ 组织）| https://github.com/ZJUFanLab/scCATCH |
| **cafferychen777/mLLMCelltype** ★641 | 多 LLM 共识投票的细胞类型注释（GPT-4/Claude/Gemini）| https://github.com/cafferychen777/mLLMCelltype |
| **CelVoxes/ceLLama** ★153 | 本地 LLM（Ollama）细胞类型注释，隐私安全，支持自定义报告 | https://github.com/CelVoxes/ceLLama |

---

## 八、可视化

| 仓库 | 简介 | 链接 |
|------|------|------|
| **dtm2451/dittoSeq** ★215 | 色盲友好的单细胞/bulk RNA-seq 可视化（R，Seurat/SCE 兼容）| https://github.com/dtm2451/dittoSeq |
| **zhanghao-njmu/SCP** ★644 | 端到端单细胞 pipeline + 高质量可视化，R 包，输出发表级图 | https://github.com/zhanghao-njmu/SCP |
| **pwwang/scplotter** ★298 | 基于 plotthis 的单细胞/空间可视化 R 包，API 简洁 | https://github.com/pwwang/scplotter |

---

## 九、本仓库流程定位

```
上游定量            中游分析                    下游建模
──────────          ─────────────────────────   ──────────────────────────
nf-core/scrnaseq    scverse/scanpy              theislab/scvelo
velocyto.py    ───► scverse/anndata        ───► theislab/cellrank  ───► GAM (pygam)
alevin-fry          ZJUFanLab/scCATCH           aertslab/pySCENIC       本仓库实现
                                                aristoteleo/dynamo
```

本仓库（`pf-gametocyte-cellrank-gam`）对应 **CellRank 1.5.x + GAM** 阶段，输入为已预处理的 h5ad，输出为命运概率 UMAP 和基因趋势图。

---

## 参考阅读

| 论文 | 对应工具 |
|------|---------|
| La Manno et al., *Nature* 2018 | velocyto（RNA velocity 提出）|
| Bergen et al., *Nature Biotechnology* 2020 | scVelo（动力学模型）|
| Lange et al., *Nature Methods* 2022 | CellRank 1.x |
| Brody et al., *Nature* 2023 | CellRank 2.x |
| Aibar et al., *Nature Methods* 2017 | SCENIC / pySCENIC |
| Mohammed et al., *Nature Communications* 2024 | 本仓库复现目标 |
