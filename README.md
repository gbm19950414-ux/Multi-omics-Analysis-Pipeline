# Multi-omics Analysis Pipeline

RNA 与脂质组的数据处理、质量检查和派生分析脚本。RNA-seq and lipidomics processing scripts, with inspectable synthetic examples.

开发来源是 multiomics_mech，是科研三仓系统的组学处理部分。本导出包括通用 RNA 预处理／合并／差异分析和脂质格式、QC、PQN／批次处理、矩阵与图形脚本；排除真实数据、课题特异的机制轴／指标、EphB1 下游推断、论文修订绘图及本机环境锁。它是开发脚本集，不是已经验证全部科学结果的产品。

## 最小例子（只处理虚构数据）

需要 Python 3、R 与 readr、dplyr、tidyr、stringr、purrr。不需要 FASTQ、真实样本、网络 API 或生物学参考文件。

```sh
python3 examples/run_demo.py --rscript Rscript
```

例子在 demo-output 下创建虚构 count 表和脂质长表，调用实际 07_merge_counts.R 与 04_within_batch_pqn.R，再检查基因／样本合并及已知强度比例。它不证明批次效应已经在真实数据中消除。

## 其他脚本

RNA shell 脚本需要 Bash 4+、FastQC、fastp、STAR 或 HISAT2、samtools、subread／featureCounts 和 MultiQC；R 下游按对应 library 声明安装 DESeq2／ashr／ggplot2／pheatmap 等。脂质可选模块还依赖 YAML、limma 及各脚本列明的包。不要把安装例子当作已跑通 FASTQ 全流程的证据。

默认路径由 scripts/rna/00_config.env 从导出根目录解析，也允许 PROJECT_DIR 环境变量覆盖；脂质脚本以调用者当前工作目录为工作空间。00_config.env 和 YAML 是需要用户按数据结构审读的配置模板，分组／批次推断和 QC 阈值不是普适科学标准。保持现有方法语义，未在本轮重新验收差异分析、交互项、机制解释或全部图形。99_run_all_batches.sh 含旧产物清理逻辑，运行真实数据前须逐步检查输入与备份；公开例子不调用它。

## 科研系统的关系

[Research Automation Workflow](https://github.com/gbm19950414-ux/Research-Automation-Workflow)负责协调与实验执行；本仓负责组学处理；具体 EphB1 开发源保留本地；[EphB1 analysis code](https://github.com/gbm19950414-ux/EphB1-multiomics-analysis-code)是固定论文修订快照。快照不随这里的通用代码更新而重新认定为科学验收。

PUBLIC_EXPORT.json 保存选中来源和公开路径适配的哈希。开发仍在原来源维护，公开代码从同一来源更新，不整推其历史。原公开仓许可保留；真实科研素材未随代码发布。
