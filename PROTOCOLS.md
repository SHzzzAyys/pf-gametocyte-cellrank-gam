# 单细胞测序分析 — 模型构建 Protocol 资源汇总

> 搜集自 GitHub，逐仓库核实信息（stars / 语言 / topics / 最近更新），按分析阶段分类，附详细功能注释。
> 最后更新：2026-05-26

---

## 分析流程总览

```
原始数据（FASTQ）
    │
    ▼ ── A. 上游定量 ──────────────────────────────────────────────
    │   nf-core/scrnaseq（流程编排）
    │   salmon / alevin-fry / kallisto|bustools（定量工具）
    │   velocyto.py（spliced/unspliced loom）
    │
    ▼ ── B. 核心分析框架 ──────────────────────────────────────────
    │   scanpy + anndata（Python 主力）
    │   muon（多模态）
    │
    ▼ ── C. RNA Velocity ─────────────────────────────────────────
    │   scVelo（动力学模型 → latent_time）
    │   dynamo（向量场重建 → 势能景观）
    │   cellDancer（深度学习 velocity）
    │
    ▼ ── D. 细胞命运 / 轨迹推断 ───────────────────────────────────
    │   CellRank（本项目核心）
    │   WOT（最优传输，时间序列）
    │   CEFCON（驱动调控因子）
    │
    ▼ ── E. 基因调控网络（GRN）────────────────────────────────────
    │   pySCENIC（TF regulon 推断）
    │   hdWGCNA（共表达模块）
    │
    ▼ ── F. 细胞类型注释 ─────────────────────────────────────────
    │   scCATCH（marker 数据库）
    │   mLLMCelltype（多 LLM 共识）
    │   ceLLama（本地 LLM）
    │
    ▼ ── G. 可视化 ────────────────────────────────────────────────
        SCP / scplotter / dittoSeq / SCope
```

---

## A. 综合教程 / 最佳实践

### 1. theislab/single-cell-best-practices
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 1,198 |
| 语言 | Jupyter Notebook |
| 维护状态 | 活跃（2026-05-26 更新）|
| 官网 | https://www.sc-best-practices.org |
| 链接 | https://github.com/theislab/single-cell-best-practices |

**功能注释：**
- Theis Lab（慕尼黑亥姆霍兹中心）出品的单细胞分析最佳实践**在线书籍 + 配套代码**
- 覆盖完整流程：原始数据处理 → QC → 标准化 → 降维 → 聚类 → 注释 → 轨迹 → 多组学
- 每一章都有可运行的 Jupyter Notebook，理论与代码并行
- **适用场景**：入门学习、实验室标准化流程制定、方法选型参考
- **推荐程度**：⭐⭐⭐⭐⭐（新手和专家都适用的核心参考）

---

### 2. agitter/single-cell-pseudotime
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 439 |
| 语言 | Markdown（文档型）|
| 维护状态 | 活跃（2026-04-28 更新）|
| 链接 | https://github.com/agitter/single-cell-pseudotime |

**功能注释：**
- 纯文档型仓库：系统整理 **50+ 种伪时间 / 轨迹推断算法**，每种附论文、代码链接、适用场景
- 按方法类型分类：扩散映射、MST、贝叶斯、深度学习等
- 包含方法比较表：准确性、可扩展性、是否需要已知起点等维度
- **适用场景**：选轨迹方法前的系统调研，快速找到目标场景下的最优工具
- **推荐程度**：⭐⭐⭐⭐（方法选型必读清单）

---

### 3. scRNA-tools/scRNA-tools
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 339 |
| 语言 | R |
| 维护状态 | 活跃（2026-05-13 更新）|
| 链接 | https://github.com/scRNA-tools/scRNA-tools |

**功能注释：**
- 单细胞工具**数据库 + 分类索引**，收录 500+ 个工具
- 按功能标签分类（聚类、轨迹、批次校正、空间、多组学等），可筛选语言（R/Python）
- 每个工具附发表论文、GitHub 链接、最后更新时间
- **适用场景**：查找特定任务工具、横向比较同类工具、跟踪领域新工具
- **推荐程度**：⭐⭐⭐⭐（工具检索手册，适合反复查阅）

---

### 4. kpatel427/YouTubeTutorials
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 378 |
| 语言 | R |
| 维护状态 | 活跃（2026-05-26 更新）|
| 链接 | https://github.com/kpatel427/YouTubeTutorials |

**功能注释：**
- 配套 YouTube 频道的单细胞分析代码库，教程驱动
- 主要覆盖：Seurat（聚类/注释）、Monocle3（轨迹）、DESeq2（差异表达）、scType 等
- 代码简洁，直接可运行，适合快速上手
- **适用场景**：视频 + 代码同步学习，R 用户入门首选
- **推荐程度**：⭐⭐⭐（教学向，适合初学者）

---

## B. 上游定量 / 预处理

### 5. nf-core/scrnaseq
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 328 |
| 语言 | Nextflow |
| 维护状态 | 活跃（2026-05-21 更新）|
| 链接 | https://github.com/nf-core/scrnaseq |

**功能注释：**
- nf-core 标准化单细胞 RNA-seq 上游流程，支持 **10x / Drop-seq / Smart-seq**
- 内置多种比对器：STARsolo / Alevin / Alevin-fry / kb-python / CellRanger
- 自动空液滴过滤（EmptyDrops），生成标准化的 h5ad 输出
- **适用场景**：生产环境批量处理，HPC/云计算部署，多样本并行
- **输入**：原始 FASTQ + 参考基因组  **输出**：counts matrix（h5ad/mtx）
- **推荐程度**：⭐⭐⭐⭐⭐（上游处理首选，可重复性强）

---

### 6. COMBINE-lab/salmon ★885
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 885 |
| 语言 | C++ |
| 维护状态 | 活跃（2026-05-14 更新）|
| 链接 | https://github.com/COMBINE-lab/salmon |

**功能注释：**
- 转录本定量工具，使用 **selective alignment**（速度 vs 准确率平衡优）
- 单细胞模式：`alevin`，处理 10x / Drop-seq barcode，输出 UMI counts
- 相比 STAR 等比对器速度提升 10-100×，内存占用低
- **适用场景**：bulk RNA-seq 定量 + 单细胞 alevin 模式（已被 alevin-fry 继承）
- **推荐程度**：⭐⭐⭐⭐（基础工具，alevin-fry 的上游依赖）

---

### 7. COMBINE-lab/alevin-fry ★207
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 207 |
| 语言 | Rust |
| 维护状态 | 活跃（2026-05-22 更新）|
| 链接 | https://github.com/COMBINE-lab/alevin-fry |

**功能注释：**
- alevin 的**高速 Rust 重写版**，单细胞 UMI 定量
- 内存占用比 alevin 降低 10×，速度提升 5× 以上
- 支持 10x Chromium / BD Rhapsody / SPLiT-seq / 多模态 feature barcoding
- **适用场景**：大规模单细胞数据（10万+细胞），内存受限环境
- **推荐程度**：⭐⭐⭐⭐（高效上游处理，适合规模化分析）

---

### 8. pachterlab/kallistobustools ★121
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 121 |
| 语言 | Python / Shell |
| 维护状态 | 较活跃 |
| 链接 | https://github.com/pachterlab/kallistobustools |

**功能注释：**
- `kallisto | bustools` 工作流（kb-python 封装）
- 最快的单细胞预处理方案之一，原理：pseudoalignment
- 支持 10x v2/v3、Smart-seq3、SPLiT-seq 等协议，直接输出 h5ad
- **适用场景**：快速原型验证，笔记本电脑上跑大数据集
- **推荐程度**：⭐⭐⭐（速度极快但准确率略低于 STARsolo）

---

### 9. velocyto-team/velocyto.py ★176
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 176 |
| 语言 | Python |
| 维护状态 | 仍可用（2026-05-10 有更新）|
| 链接 | https://github.com/velocyto-team/velocyto.py |

**功能注释：**
- RNA velocity 的**原始数据生成工具**：从 BAM 文件中计算 spliced / unspliced / ambiguous reads
- 两种模式：`velocyto run`（通用）、`velocyto run10x`（10x 专用）
- 输出 `.loom` 文件，供 scVelo / CellRank 使用
- **本项目关联**：P. falciparum 论文的 loom 文件即由此生成
- **适用场景**：任何需要 RNA velocity 的分析的第一步（BAM → loom）
- **注意**：速度较慢，大数据集建议改用 STARsolo 的 velocity 输出
- **推荐程度**：⭐⭐⭐⭐（velocity 分析必备前置步骤）

---

## C. 核心分析框架

### 10. scverse/scanpy ★2468
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 2,468 |
| 语言 | Python |
| 维护状态 | 高度活跃（2026-05-26 更新）|
| 官网 | https://scanpy.scverse.org |
| 链接 | https://github.com/scverse/scanpy |

**功能注释：**
- Python 单细胞分析**事实标准**，scverse 生态核心
- 覆盖完整分析链：`sc.pp.*`（预处理）→ `sc.tl.*`（分析）→ `sc.pl.*`（可视化）
- 可扩展至亿级细胞（基于 AnnData + sparse 矩阵）
- **关键 API**：`neighbors`、`umap`、`leiden`/`louvain`、`dpt`、`rank_genes_groups`
- **本项目角色**：AnnData 构建、neighbors 图、DPT 伪时间计算
- **推荐程度**：⭐⭐⭐⭐⭐（Python 单细胞分析必装）

---

### 11. scverse/anndata ★745
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 745 |
| 语言 | Python |
| 维护状态 | 高度活跃（2026-05-25 更新）|
| 链接 | https://github.com/scverse/anndata |

**功能注释：**
- 单细胞数据**容器标准**（h5ad 格式的 Python 实现）
- 统一存储：`X`（表达矩阵）、`obs`（细胞元数据）、`var`（基因元数据）、`obsm`/`obsp`/`uns`
- 所有 scverse 工具（Scanpy/scVelo/CellRank/muon）共用此格式
- **本项目角色**：所有数据操作的基础容器
- **推荐程度**：⭐⭐⭐⭐⭐（生态基础，无需选择直接用）

---

### 12. scverse/muon ★268
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 268 |
| 语言 | Python |
| 维护状态 | 活跃（2026-05-18 更新）|
| 链接 | https://github.com/scverse/muon |

**功能注释：**
- 多模态单细胞数据框架（MuData 格式扩展 AnnData）
- 支持 RNA + ATAC + Protein（CITE-seq）联合分析
- 内置加权最近邻（WNN）降维、跨模态注意力机制
- **适用场景**：多组学单细胞数据（10x Multiome、CITE-seq）
- **推荐程度**：⭐⭐⭐⭐（多模态首选，单模态用 scanpy 即可）

---

## D. RNA Velocity 模型

### 13. theislab/scvelo ★501
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 501 |
| 语言 | Python |
| 维护状态 | 活跃（2026-05-20 更新）|
| 链接 | https://github.com/theislab/scvelo |

**功能注释：**
- RNA velocity 的**动力学模型实现**（Bergen et al., *Nature Biotechnology* 2020）
- 三种模式：`steady_state`（快）→ `stochastic`（中）→ `dynamical`（最准）
- `recover_dynamics`：从剪接动力学方程恢复基因特异性参数，输出 `latent_time`
- **本项目核心依赖**：`latent_time` 作为 GAM 趋势图的 X 轴
- **关键 API**：`scv.pp.moments` → `scv.tl.velocity` → `scv.tl.recover_dynamics` → `scv.tl.latent_time`
- **版本锁定原因**：0.2.4 与 CellRank 1.5.1 API 兼容，0.3+ 有 breaking changes
- **推荐程度**：⭐⭐⭐⭐⭐（velocity 分析必备）

---

### 14. basilkhuder/Seurat-to-RNA-Velocity ★142
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 142 |
| 语言 | R / Python |
| 维护状态 | 偶发更新（2026-04-18）|
| 链接 | https://github.com/basilkhuder/Seurat-to-RNA-Velocity |

**功能注释：**
- 完整指南：将 **Seurat 对象转换并接入 scVelo / velocyto** 的全流程
- 解决 R/Python 互操作痛点：barcode 格式统一、loom 合并、AnnData 转换
- 覆盖工具链：Seurat → loom → scVelo → CellRank
- **适用场景**：从 Seurat 工作流过渡到 velocity 分析的 R 用户
- **推荐程度**：⭐⭐⭐（R 转 Python velocity 的桥接教程）

---

### 15. GuangyuWangLab2021/cellDancer ★72
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 72 |
| 语言 | Jupyter Notebook / Python |
| 维护状态 | 活跃（2026-05-25 更新）|
| 链接 | https://github.com/GuangyuWangLab2021/cellDancer |

**功能注释：**
- **深度学习 RNA velocity**：用神经网络逐细胞估计动力学参数（splicing/degradation rates）
- 克服 scVelo 假设局限：无需全局稳定假设，处理高异质性数据更准
- 基于 PyTorch，支持 GPU 加速
- **适用场景**：细胞类型混杂、发育轨迹复杂的数据；scVelo velocity 方向错误时的替代方案
- **推荐程度**：⭐⭐⭐（新方法，适合 scVelo 结果不理想时尝试）

---

### 16. aristoteleo/dynamo-release ★498
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 498 |
| 语言 | Python |
| 维护状态 | 高度活跃（2026-05-25 更新）|
| 链接 | https://github.com/aristoteleo/dynamo-release |

**功能注释：**
- RNA velocity 的**理论最完整实现**，超越 scVelo 的向量场框架
- 核心功能：
  - 向量场重建（神经网络插值，光滑连续向量场）
  - **微分几何分析**：散度、旋度、曲率、加速度
  - **势能景观**（Waddington landscape）可视化
  - 细胞命运谱系预测、转分化路径分析
  - 支持代谢标记 scRNA-seq（sci-fate / scNT-seq / SLAM-seq）
- **适用场景**：需要量化细胞命运转变概率、势能景观的深度分析
- **推荐程度**：⭐⭐⭐⭐（velocity 进阶分析首选，学习成本较高）

---

## E. 细胞命运 / 轨迹推断

### 17. theislab/cellrank ★450
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 450 |
| 语言 | Python |
| 维护状态 | 高度活跃（2026-05-26 更新）|
| 链接 | https://github.com/theislab/cellrank |

**功能注释：**
- 细胞命运分析框架，核心：**马尔可夫链 + GPCCA**
- 版本差异（本项目用 1.5.1）：
  - `1.5.x`：`cellrank.tl.kernels.VelocityKernel` → `GPCCA` → `compute_absorption_probabilities`
  - `2.x`：API 重构，`cellrank.kernels.*`，方法名变化，不兼容
- 主要输出：终态概率 → 吸收概率（命运概率）→ lineage driver 基因
- **关键 kernel**：VelocityKernel（RNA velocity 驱动）、ConnectivityKernel（图拓扑驱动）
- **本项目核心**：整套命运分析流程的主体
- **推荐程度**：⭐⭐⭐⭐⭐（细胞命运分析必装，注意版本）

---

### 18. scverse/cellrank_notebooks ★8
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 8 |
| 语言 | Jupyter Notebook |
| 维护状态 | 2026-03-09 更新 |
| 链接 | https://github.com/scverse/cellrank_notebooks |

**功能注释：**
- CellRank **官方教程 notebooks 集合**（对应 CellRank 2.x 文档）
- 覆盖多种场景：VelocityKernel、ConnectivityKernel、RealTimeKernel（时间序列）、CytoTRACEKernel
- 每个 notebook 对应一个具体生物问题（造血、胰腺内分泌发育等）
- **适用场景**：照着跑通 CellRank 2.x 标准流程（注意：本项目用 1.5.x，API 不同）
- **推荐程度**：⭐⭐⭐⭐（官方教程，新项目建议参考 2.x）

---

### 19. broadinstitute/wot ★162
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 162 |
| 语言 | Jupyter Notebook / Python |
| 维护状态 | 偶发更新（2026-04-21）|
| 链接 | https://github.com/broadinstitute/wot |

**功能注释：**
- **Waddington OT**：用最优传输（Optimal Transport）分析跨**时间点**的细胞命运
- 不依赖 RNA velocity，通过时间序列快照推断细胞转变概率
- 输出：转变矩阵（任意两个时间点间细胞的命运流向）
- **适用场景**：有时间标签的单细胞数据（如发育时间序列、药物处理时间点）
- **与 CellRank 的区别**：WOT 用时间点信息，CellRank 用 velocity/拓扑
- **推荐程度**：⭐⭐⭐⭐（时间序列数据的命运分析首选）

---

### 20. WPZgithub/CEFCON ★31
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 31 |
| 语言 | Jupyter Notebook / Python |
| 维护状态 | 活跃（2026-05-11 更新）|
| 链接 | https://github.com/WPZgithub/CEFCON |

**功能注释：**
- 从 scRNA-seq 中解码**细胞命运决策的驱动调控因子**
- 方法：图神经网络（GNN）+ 网络控制理论 + 注意力机制
- 将 GRN 推断与命运决策结合，找出控制分化方向的 master regulator TF
- **输出**：核心调控因子排序、命运决策调控网络
- **适用场景**：想知道"哪些 TF 驱动了这条分化路径"
- **推荐程度**：⭐⭐⭐（新颖方法，适合深度机制分析）

---

## F. 基因调控网络（GRN）

### 21. aertslab/pySCENIC ★606
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 606 |
| 语言 | Python |
| 维护状态 | 活跃（2026-05-26 更新）|
| 链接 | https://github.com/aertslab/pySCENIC |

**功能注释：**
- SCENIC 的 Python 实现：单细胞**转录因子调控网络推断**
- 三步流程：
  1. `GRNBoost2`：基因共表达模块（TF → target gene）
  2. `cisTarget`：motif 富集过滤（保留有 DNA binding 证据的连接）
  3. `AUCell`：每个细胞的 regulon 活性打分
- **本论文相关**：Mohammed et al. 2024 用 SCENIC 分析 P. falciparum TF 活性（图中 Regulon 列）
- **适用场景**：推断转录因子调控网络、细胞类型特异性 regulon 识别
- **推荐程度**：⭐⭐⭐⭐⭐（GRN 分析金标准）

---

### 22. smorabit/hdWGCNA ★475
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 475 |
| 语言 | R |
| 维护状态 | 活跃（2026-05-25 更新）|
| 链接 | https://github.com/smorabit/hdWGCNA |

**功能注释：**
- **单细胞 / 空间转录组版 WGCNA**：高维加权基因共表达网络分析
- 在伪细胞（metacell）水平构建共表达模块，克服单细胞稀疏噪声
- 支持 scRNA-seq 和 Visium 空间转录组
- 输出：基因模块（hub gene 识别）、模块特异性细胞类型、模块与表型的关联
- **适用场景**：无先验 TF motif 数据时的 GRN 分析；功能模块发现
- **与 pySCENIC 的区别**：SCENIC 基于 motif 数据库，hdWGCNA 纯数据驱动
- **推荐程度**：⭐⭐⭐⭐（R 生态 GRN 分析首选）

---

### 23. aertslab/SCope ★79
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 79 |
| 语言 | Python / React.js |
| 维护状态 | 偶发更新（2026-04-22）|
| 链接 | https://github.com/aertslab/SCope |

**功能注释：**
- pySCENIC 配套的**交互可视化平台**（Web 应用）
- 支持 loom 格式大规模单细胞数据的浏览：UMAP、基因表达、regulon 活性同时展示
- 部署在云端（AWS），也可本地运行
- **适用场景**：SCENIC 结果的探索性可视化，生成论文展示图
- **推荐程度**：⭐⭐⭐（与 pySCENIC 配套使用，独立价值有限）

---

## G. 细胞类型注释

### 24. ZJUFanLab/scCATCH ★241
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 241 |
| 语言 | R |
| 维护状态 | 活跃（2026-05-04 更新）|
| 链接 | https://github.com/ZJUFanLab/scCATCH |

**功能注释：**
- 基于 **marker gene 数据库**的自动细胞类型注释
- 内置数据库覆盖 100+ 种组织、400+ 种细胞类型（人 / 小鼠）
- 方法：对 Seurat/Scanpy 聚类结果，用统计打分匹配 marker
- **输出**：每个 cluster 的细胞类型预测 + 置信度评分
- **适用场景**：无参考数据集时的标准组织注释（人类常见组织首选）
- **局限**：寄生虫（P. falciparum）等非标准物种没有数据库支持
- **推荐程度**：⭐⭐⭐⭐（人/鼠组织注释快速可靠）

---

### 25. cafferychen777/mLLMCelltype ★641
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 641 |
| 语言 | Python |
| 维护状态 | 高度活跃（2026-05-25 更新）|
| 链接 | https://github.com/cafferychen777/mLLMCelltype |

**功能注释：**
- **多 LLM 共识投票**细胞类型注释：GPT-4o / Claude 3.5 / Gemini Pro 等多模型并行推断
- 根据 marker gene 列表让多个 LLM 独立预测，用 Delphi 共识算法整合
- 支持 Seurat（R）和 Scanpy（Python）
- **优势**：减少单一 LLM 幻觉，置信度可量化
- **适用场景**：非标准物种/组织注释（LLM 有更广泛的知识库）、常规注释加速
- **推荐程度**：⭐⭐⭐⭐（LLM 注释工具中最成熟的方案）

---

### 26. CelVoxes/ceLLama ★153
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 153 |
| 语言 | R |
| 维护状态 | 2026-02-24 更新 |
| 链接 | https://github.com/CelVoxes/ceLLama |

**功能注释：**
- 用**本地 LLM（Ollama）**进行细胞类型注释，数据不出本地
- 支持 llama3 / mistral / gemma 等开源模型
- 生成包含推理过程的自定义报告
- **适用场景**：数据隐私敏感场景（临床数据、未发表研究）、无 API 费用预算
- **与 mLLMCelltype 的区别**：本地运行 vs 云端 API，准确率略低但完全私有
- **推荐程度**：⭐⭐⭐（隐私优先场景下的 LLM 注释方案）

---

## H. 可视化

### 27. zhanghao-njmu/SCP ★644
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 644 |
| 语言 | R |
| 维护状态 | 高度活跃（2026-05-26 更新）|
| 链接 | https://github.com/zhanghao-njmu/SCP |

**功能注释：**
- **端到端单细胞 pipeline**：从 Seurat/AnnData 到发表级图一站完成
- 覆盖全流程：QC → 降维 → 聚类 → 注释 → 差异分析 → 轨迹 → 空间
- 内置 Shiny 交互界面，可视化功能极丰富（35+ 种图形类型）
- 支持 scRNA-seq / scATAC-seq / 空间转录组
- **适用场景**：需要快速出高质量图的研究，不想手写可视化代码
- **推荐程度**：⭐⭐⭐⭐⭐（R 用户综合 pipeline 首选）

---

### 28. pwwang/scplotter ★298
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 298 |
| 语言 | HTML / R |
| 维护状态 | 活跃（2026-05-23 更新）|
| 链接 | https://github.com/pwwang/scplotter |

**功能注释：**
- 基于 `plotthis` 的单细胞**可视化 R 包**，API 极简洁
- 支持 Seurat、Giotto（空间）对象
- 擅长：细胞比例图、marker 热图、轨迹可视化、空间基因表达图
- **适用场景**：快速生成论文质量可视化，补充 Seurat 默认图的不足
- **推荐程度**：⭐⭐⭐⭐（R 用户可视化补充工具）

---

### 29. dtm2451/dittoSeq ★215
| 项目 | 值 |
|------|-----|
| ⭐ Stars | 215 |
| 语言 | R |
| 维护状态 | 活跃（2026-05-22 更新）|
| 链接 | https://github.com/dtm2451/dittoSeq |

**功能注释：**
- **色盲友好**的单细胞 / bulk RNA-seq 可视化（基于 Bioconductor）
- 支持 Seurat 对象和 SingleCellExperiment 对象（双平台兼容）
- 内置 20+ 种图形：UMAP、violin、dot plot、热图、条形图等
- 配色方案经色盲测试，适合论文发表
- **适用场景**：论文图表制作，需要色盲可访问性的场合
- **推荐程度**：⭐⭐⭐⭐（Bioconductor 生态可视化首选）

---

## 汇总对比表

| # | 仓库 | ⭐ | 语言 | 阶段 | 活跃度 | 推荐度 |
|---|------|-----|------|------|--------|--------|
| 1 | theislab/single-cell-best-practices | 1198 | Jupyter | 教程 | ✅ | ⭐⭐⭐⭐⭐ |
| 2 | agitter/single-cell-pseudotime | 439 | 文档 | 方法综述 | ✅ | ⭐⭐⭐⭐ |
| 3 | scRNA-tools/scRNA-tools | 339 | R | 工具索引 | ✅ | ⭐⭐⭐⭐ |
| 4 | kpatel427/YouTubeTutorials | 378 | R | 教程 | ✅ | ⭐⭐⭐ |
| 5 | nf-core/scrnaseq | 328 | Nextflow | 上游流程 | ✅ | ⭐⭐⭐⭐⭐ |
| 6 | COMBINE-lab/salmon | 885 | C++ | 定量 | ✅ | ⭐⭐⭐⭐ |
| 7 | COMBINE-lab/alevin-fry | 207 | Rust | 定量 | ✅ | ⭐⭐⭐⭐ |
| 8 | pachterlab/kallistobustools | 121 | Python | 定量 | 较活跃 | ⭐⭐⭐ |
| 9 | velocyto-team/velocyto.py | 176 | Python | Velocity 前处理 | ✅ | ⭐⭐⭐⭐ |
| 10 | scverse/scanpy | 2468 | Python | 核心框架 | ✅ | ⭐⭐⭐⭐⭐ |
| 11 | scverse/anndata | 745 | Python | 数据格式 | ✅ | ⭐⭐⭐⭐⭐ |
| 12 | scverse/muon | 268 | Python | 多模态框架 | ✅ | ⭐⭐⭐⭐ |
| 13 | theislab/scvelo | 501 | Python | RNA Velocity | ✅ | ⭐⭐⭐⭐⭐ |
| 14 | basilkhuder/Seurat-to-RNA-Velocity | 142 | R/Python | Velocity 桥接 | 偶发 | ⭐⭐⭐ |
| 15 | GuangyuWangLab2021/cellDancer | 72 | Python | RNA Velocity（DL）| ✅ | ⭐⭐⭐ |
| 16 | aristoteleo/dynamo-release | 498 | Python | Velocity 进阶 | ✅ | ⭐⭐⭐⭐ |
| 17 | theislab/cellrank | 450 | Python | 细胞命运 | ✅ | ⭐⭐⭐⭐⭐ |
| 18 | scverse/cellrank_notebooks | 8 | Jupyter | CellRank 教程 | 较活跃 | ⭐⭐⭐⭐ |
| 19 | broadinstitute/wot | 162 | Python | 时间序列命运 | 偶发 | ⭐⭐⭐⭐ |
| 20 | WPZgithub/CEFCON | 31 | Python | 命运驱动因子 | ✅ | ⭐⭐⭐ |
| 21 | aertslab/pySCENIC | 606 | Python | GRN 推断 | ✅ | ⭐⭐⭐⭐⭐ |
| 22 | smorabit/hdWGCNA | 475 | R | GRN 共表达 | ✅ | ⭐⭐⭐⭐ |
| 23 | aertslab/SCope | 79 | Python | GRN 可视化 | 偶发 | ⭐⭐⭐ |
| 24 | ZJUFanLab/scCATCH | 241 | R | 细胞注释 | ✅ | ⭐⭐⭐⭐ |
| 25 | cafferychen777/mLLMCelltype | 641 | Python | 细胞注释（LLM）| ✅ | ⭐⭐⭐⭐ |
| 26 | CelVoxes/ceLLama | 153 | R | 细胞注释（本地 LLM）| 偶发 | ⭐⭐⭐ |
| 27 | zhanghao-njmu/SCP | 644 | R | 综合可视化 | ✅ | ⭐⭐⭐⭐⭐ |
| 28 | pwwang/scplotter | 298 | R | 可视化 | ✅ | ⭐⭐⭐⭐ |
| 29 | dtm2451/dittoSeq | 215 | R | 可视化 | ✅ | ⭐⭐⭐⭐ |
