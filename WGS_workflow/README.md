# 基因组重测序流程使用指南

## 目录

- [基因组重测序流程使用指南](#基因组重测序流程使用指南)
  - [目录](#目录)
  - [简介](#简介)
  - [环境依赖](#环境依赖)
  - [使用方法](#使用方法)
    - [参数说明](#参数说明)
    - [示例命令](#示例命令)
  - [流程步骤](#流程步骤)
    - [步骤 0：设置工作目录和准备索引文件](#步骤-0设置工作目录和准备索引文件)
    - [步骤 1：使用 Fastp 进行质量控制和剪切](#步骤-1使用-fastp-进行质量控制和剪切)
    - [步骤 2：使用 BWA-MEM 进行比对和排序](#步骤-2使用-bwa-mem-进行比对和排序)
    - [步骤 3：计算数据指标](#步骤-3计算数据指标)
    - [步骤 4：去除重复序列](#步骤-4去除重复序列)
    - [步骤 5：生成 MultiQC 报告](#步骤-5生成-multiqc-报告)
    - [步骤 6：使用 FreeBayes 进行变异检测](#步骤-6使用-freebayes-进行变异检测)
    - [步骤 7：使用 bcftools 过滤变异](#步骤-7使用-bcftools-过滤变异)
  - [输出文件说明](#输出文件说明)
  - [注意事项](#注意事项)

---

## 简介

一个用于基因组重测序分析的自动化流程脚本。该脚本从原始测序数据（Fastq 格式）开始，经过质量控制、比对、变异检测和过滤等步骤，最终生成高质量的变异结果文件（VCF 格式）。


## 环境依赖

在运行该脚本之前，请确保已安装以下软件并配置好环境：

- **Bash**：Unix Shell 脚本解释器
- **Fastp**：0.20.1 或以上版本
- **BWA**：0.7.17 或以上版本
- **Samtools**：1.9 或以上版本
- **Qualimap**：2.2.1 或以上版本
- **Picard**：2.20.8 或以上版本
- **MultiQC**：1.7 或以上版本
- **FreeBayes**：1.3.1 或以上版本
- **Bcftools**：1.9 或以上版本
- **Conda/Mamba**：用于环境管理

**注意**：建议使用 Conda 或 Mamba 创建独立的环境来管理上述软件，以避免版本冲突。该脚本使用 conda run 解决Picard与其他软件之间的环境冲突。请根据自己的系统环境进行修改。

---

## 使用方法

### 参数说明

```bash
bash script.sh [选项] -f reference -s sample -1 fq1 -2 fq2
```

- `-f, --reference`：**必需**，参考基因组文件（FASTA 格式）
- `-s, --sample`：**必需**，样本名称
- `-1, --fq1`：**必需**，前向 reads 文件（FASTQ 格式）
- `-2, --fq2`：**必需**，反向 reads 文件（FASTQ 格式）
- `-t, --threads`：可选，使用的线程数，默认值为 5
- `-h, --help`：显示帮助信息并退出

### 示例命令

```bash
bash script.sh -f reference.fa -s sample1 -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -t 8
```

---

## 流程步骤

### 步骤 0：设置工作目录和准备索引文件

**目的**：设置工作目录，检查并准备参考基因组的索引文件。

**操作**：

1. **设置工作目录**：脚本将在 `~/WGS_workflow/` 下创建以样本名命名的目录，用于存放中间文件和结果文件。
   
2. **检查 BWA 索引文件**：
   - 脚本会检查参考基因组目录下是否存在以下 BWA 索引文件：
     - `*.bwt`、`*.ann`、`*.amb`、`*.pac`、`*.sa`
   - 如果不存在，将自动运行 `bwa index` 创建索引。

3. **检查 Samtools 索引文件**：
   - 脚本会检查是否存在参考基因组的 `.fai` 索引文件。
   - 如果不存在，将自动运行 `samtools faidx` 创建索引。

### 步骤 1：使用 Fastp 进行质量控制和剪切

**目的**：对原始测序数据进行质量过滤和剪切，去除低质量数据和接头序列，提高后续分析的准确性。

**命令**：

```bash
fastp --thread [线程数] \
      --qualified_quality_phred 15 \
      --unqualified_percent_limit 40 \
      --n_base_limit 10 \
      --length_required 50 \
      --detect_adapter_for_pe \
      -i [前向 reads] -I [反向 reads] \
      -o [前向 clean reads] -O [反向 clean reads] \
      -h [HTML 报告] -j [JSON 报告]
```

**参数说明**：

- `--qualified_quality_phred 15`：设定合格的 Phred 质量值阈值为 15。
- `--unqualified_percent_limit 40`：允许每个 read 中低质量碱基的最大比例为 40%。
- `--n_base_limit 10`：允许每个 read 中 N 碱基的最大数量为 10。
- `--length_required 50`：过滤掉长度小于 50 bp 的 reads。
- `--detect_adapter_for_pe`：自动检测并剪切双端测序的接头序列。

**输出**：

- 清洗后的前向 reads：`[样本名]_1.clean.fq.gz`
- 清洗后的反向 reads：`[样本名]_2.clean.fq.gz`
- 质量控制报告（HTML 和 JSON 格式）

### 步骤 2：使用 BWA-MEM 进行比对和排序

**目的**：将清洗后的 reads 比对到参考基因组上，并对比对结果进行排序和索引。

**命令**：

```bash
bwa mem -t [线程数] -M \
        -R "@RG\tID:G\tLB:[样本名]\tPL:ILLUMINA\tSM:[样本名]" \
        [参考基因组] [前向 clean reads] [反向 clean reads] | \
    samtools sort -@ [线程数] -o [样本名].sorted.bam

samtools index [样本名].sorted.bam
samtools flagstat [样本名].sorted.bam > [样本名]_flagstat.txt
```

**参数说明**：

- `-M`：标记比对中的次优比对，以兼容 Picard。
- `-R`：指定 Read Group 信息，包含：
  - `ID`：Read Group 标识，设为 `G`。
  - `LB`：文库（Library）名称，设为样本名。
  - `PL`：测序平台，设为 `ILLUMINA`。
  - `SM`：样本名。

**输出**：

- 排序后的 BAM 文件：`[样本名].sorted.bam`
- BAM 索引文件：`[样本名].sorted.bam.bai`
- 比对统计信息：`[样本名]_flagstat.txt`

### 步骤 3：计算数据指标

**目的**：评估 BAM 文件的质量，包括测序深度、覆盖度等指标。

**命令**：

```bash
qualimap bamqc -bam [样本名].sorted.bam -outdir [样本名]_bamqc
```

**输出**：

- 质量评估报告目录：`[样本名]_bamqc`

### 步骤 4：去除重复序列

**目的**：使用 Picard 标记并去除 PCR 或光学重复的 reads，避免影响后续变异调用的准确性。

**命令**：

```bash
mamba run -n picard picard MarkDuplicates \
    --INPUT [样本名].sorted.bam \
    --OUTPUT [样本名]_dup.bam \
    --METRICS_FILE [样本名]_dup_metrics.txt
```

**输出**：

- 去除重复后的 BAM 文件：`[样本名]_dup.bam`
- 重复序列指标文件：`[样本名]_dup_metrics.txt`

### 步骤 5：生成 MultiQC 报告

**目的**：整合各步骤的质量控制报告，生成综合性的报告，便于查看和分析。

**命令**：

```bash
multiqc [工作目录] \
    --module fastp --module picard --module samtools \
    -o [样本名]_multiqc_report
```

**输出**：

- MultiQC 报告目录：`[样本名]_multiqc_report`

### 步骤 6：使用 FreeBayes 进行变异检测

**目的**：对去除重复后的 BAM 文件进行变异检测，生成原始的 VCF 文件。

**命令**：

```bash
freebayes -f [参考基因组] [样本名]_dup.bam > [样本名]_raw.vcf
```

**输出**：

- 原始变异结果文件：`[样本名]_raw.vcf`

### 步骤 7：使用 bcftools 过滤变异

**目的**：对原始 VCF 文件进行质量过滤，保留高可信度的变异。

**命令**：

```bash
bcftools filter \
    --include 'QUAL>30 && INFO/DP>10' \
    --SnpGap 5 --IndelGap 10 \
    -O v -o [样本名]_filter.vcf [样本名]_raw.vcf
```

**参数说明**：

- `QUAL>30`：仅保留质量值大于 30 的变异。
- `INFO/DP>10`：仅保留测序深度大于 10 的变异。
- `--SnpGap 5`：过滤与其他 SNP 距离小于 5 bp 的 SNP。
- `--IndelGap 10`：过滤与其他 Indel 距离小于 10 bp 的 Indel。

**输出**：

- 过滤后的变异结果文件：`[样本名]_filter.vcf`

---

## 输出文件说明

| 文件类型                   | 文件名格式                     | 描述                         |
|----------------------------|--------------------------------|------------------------------|
| 清洗后的fq序列文件            | `[样本名]_1.clean.fq.gz`       | 清洗后的前向 reads 文件       |
|                            | `[样本名]_2.clean.fq.gz`       | 清洗后的反向 reads 文件       |
| **fastp 报告**             | `[样本名]_fastp.html`          | fastp 生成的 HTML 报告        |
|                            | `[样本名]_fastp.json`          | fastp 生成的 JSON 报告        |
| **排序后的 BAM 文件及索引** | `[样本名].sorted.bam`          | 排序后的 BAM 文件             |
|                            | `[样本名].sorted.bam.bai`      | BAM 文件的索引               |
| **比对统计信息**           | `[样本名]_flagstat.txt`        | Samtools 生成的比对统计信息   |
| **Qualimap 报告**          | `[样本名]_bamqc` 目录          | Qualimap 生成的 BAM 质量评估  |
| **去除重复后的 BAM 文件**  | `[样本名]_dup.bam`             | 去除重复序列后的 BAM 文件     |
| **重复序列指标文件**       | `[样本名]_dup_metrics.txt`     | Picard 生成的重复序列指标报告 |
| **MultiQC 报告**           | `[样本名]_multiqc_report` 目录 | MultiQC 整合的质量控制报告    |
| **原始变异结果文件**       | `[样本名]_raw.vcf`             | FreeBayes 生成的原始 VCF 文件 |
| **过滤后的变异结果文件**   | `[样本名]_filter.vcf`          | 过滤后的高质量 VCF 文件       |

---

## 注意事项

1. **参考基因组**：
   - **指定方式**：通过命令行参数 `-f` 或 `--reference` 指定参考基因组文件（FASTA 格式）。
   - **索引文件**：脚本会自动检查参考基因组所在目录是否存在所需的 BWA 和 Samtools 索引文件，如果不存在则自动创建。这可能会耗费一些时间，尤其是对于较大的基因组。

2. **软件环境**：建议使用 Conda 或 Mamba 创建独立环境，安装所需软件包，避免版本冲突。

3. **参数调整**：根据实际数据情况，适当调整脚本中的参数，如质量阈值、线程数等。

4. **数据安全**：处理完毕后，脚本会删除临时目录，请确保重要文件已备份。

5. **脚本错误**：在执行过程中，如果遇到错误，请检查输入参数是否正确，或者查看日志信息获取更多详情。

6. **输入文件格式**：确保输入的 FASTQ 文件为标准格式，且与参考基因组序列匹配。

7. **比对参数**：在 BWA 比对步骤中，脚本默认使用了特定的参数（如 `-M`），可根据需要进行调整。

---

**注**：在使用脚本前，请仔细阅读并理解各步骤，根据自己的数据和需求进行调整。