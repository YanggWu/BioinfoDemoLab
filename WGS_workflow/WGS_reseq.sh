#!/bin/bash
############################################################################################
# 上海百兆瑞生物科技有限公司
# Description: Genome resequencing pipeline
# Author: wuyang
# Date: 2024.11.07
# Version: 1.0
# Last Update: 2024.11.08
#############################################################################################
# 检查是否在 Bash 中运行, 如果不是则切换到 Bash。
if [ -z "$BASH_VERSION" ]; then
  echo "Switching to bash..."
  exec /bin/bash "$0" "$@"
  exit
fi
# 帮助信息
function usage() {
  echo "Usage: bash $0 [options] -f reference -s sample -1 fq1 -2 fq2"
  echo
  echo "Options:"
  echo "  -f, --reference  Reference genome (fasta file)"
  echo "  -s, --sample     Sample name"
  echo "  -1, --fq1        Forward reads (fastq file)"
  echo "  -2, --fq2        Reverse reads (fastq file)"
  echo "  -t, --threads    Number of threads (default: 5)"
  echo "  -h, --help       Show this help message and exit"
  exit 1
}

# Default parameters
nt=5 # Number of threads

# Parse command-line arguments
# 解析命令行参数
# 解析命令行参数
while [[ "$#" -gt 0 ]]; do
  case $1 in
  -s | --sample)
    sample="$2"
    shift
    ;;
  -1 | --fq1)
    fq1="$2"
    shift
    ;;
  -2 | --fq2)
    fq2="$2"
    shift
    ;;
  -f | --reference)
    reference="$2"
    shift
    ;;
  -t | --threads)
    nt="$2"
    shift
    ;;
  -h | --help) usage ;;
  *)
    echo "未知参数: $1"
    usage
    ;;
  esac
  shift
done

# Check required parameters
if [[ -z "$sample" || -z "$fq1" || -z "$fq2" || -z "$reference" ]]; then
  echo -e "\033[31mError\033[0m: Sample name, input fastq files, and reference genome are required." >&2
  usage
fi

i=${sample}
# ******************************************
# 0. Setup	设置工作目录,准备索引文件和比对信息。
# ******************************************
echo -e "\n\n\033[34m Step 0 \033[0m: Setup work directory and prepare index files\n"
#### ------工作目录------ ####
workdir=~/WGS_workflow/${i}
[ ! -d $workdir ] && mkdir -p $workdir
tmpdir=$workdir/temp
[ ! -d $tmpdir ] && mkdir -p $tmpdir
cd $workdir
echo -e "\n\033[32mWorking directory: $workdir\033[0m"



# 获取参考基因组所在目录和文件名
refdir=$(dirname "$reference")
refname=$(basename "$reference")
refbase=${refname%.*}

REF="$reference"
BWA_INDEX="$REF"

# 检查并创建 BWA 索引文件
echo -e "\n\033[34mChecking BWA index files...\033[0m"

if [ -e "${BWA_INDEX}.bwt" ] && [ -e "${BWA_INDEX}.ann" ] && [ -e "${BWA_INDEX}.amb" ] && [ -e "${BWA_INDEX}.pac" ] && [ -e "${BWA_INDEX}.sa" ]; then
    echo -e "\n\033[32mBWA index files found.\033[0m"
else
    echo -e "\n\033[33mBWA index files not found. Creating BWA index...\033[0m"
    bwa index "$BWA_INDEX"
fi

# 检查并创建 Samtools faidx 索引文件
if [ -e "${REF}.fai" ]; then
    echo -e "\033[32mSamtools faidx index found.\033[0m"
else
    echo -e "\033[33mSamtools faidx index not found. Creating faidx index...\033[0m"
    samtools faidx "$REF"
fi

# bam 文件比对信息，避免后续变异调用时出现问题
group="G"
platform="ILLUMINA"
mq=30

#### ------工作目录------ ####
workdir=~/WGS_workflow/${i}
[ ! -d $workdir ] && mkdir -p $workdir
tmpdir=$workdir/temp
[ ! -d $tmpdir ] && mkdir -p $tmpdir
cd $workdir

# ******************************************
# 1. 利用 Fastp 对原始测序数据过滤
# ******************************************
echo -e "\n\n\033[34m Step 1 \033[0m: Fastp filter and trim\n"

# Output files
fq1_clean=${sample}_1.clean.fq.gz
fq2_clean=${sample}_2.clean.fq.gz
html=${sample}_fastp.html
json=${sample}_fastp.json

fastp --thread ${nt} \
  --qualified_quality_phred 15 \
  --unqualified_percent_limit 40 \
  --n_base_limit 10 \
  --length_required 50 \
  --detect_adapter_for_pe \
  -i ${fq1} -I ${fq2} \
  -o ${fq1_clean} -O ${fq2_clean} \
  -h ${html} -j ${json}

# ******************************************
# 2. 利用 BWA-MEM 进行比对并排序
# ******************************************
echo -e "\n\n\033[34m Step 2 \033[0m: BWA-MEM mapping and sorting\n"

# Output file
sortedbam=$i.sorted.bam

bwa mem \
  -t ${nt} -M \
  -R "@RG\tID:${group}\tLB:${i}\tPL:${platform}\tSM:${i}" \
  ${BWA_INDEX} ${fq1_clean} ${fq1_clean} |
  samtools sort -@ ${nt} >$sortedbam && samtools index $sortedbam
# Output mapping statistics
samtools flagstat $sortedbam >${i}_flagstat.txt

# ******************************************
# 3. Calculate data metrics
# ******************************************
echo -e "\n\n\033[34m Step 3 \033[0m: Calculate data metrics\n"

# Output directory
bamqcdir=${i}_bamqc # qualimap 输出目录
# 测序深度与覆盖度统计
qualimap bamqc \
  -bam $sortedbam \
  -outdir $bamqcdir

# ******************************************
# 4. 去除 Duplicate Reads
# ******************************************
echo -e "\n\n\033[34m Step 4 \033[0m: Remove Duplicate Reads\n"

mamba run -n picard picard MarkDuplicates \
  --INPUT $sortedbam \
  --OUTPUT ${i}_dup.bam \
  --METRICS_FILE ${i}_dup_metrics.txt

# ******************************************
# 5. MultiQC 数据整合
# ******************************************
echo -e "\n\n\033[34m Step 5 \033[0m: MultiQC report\n"

multiqc $workdir \
  --module fastp --module picard --module samtools \
  -o $workdir/${i}_multiqc_report

# ******************************************
# 6. freebayes 变异检测
# ******************************************
echo -e "\n\n\033[34m Step 6 \033[0m: Variant calling with Freebayes\n"

freebayes -f $REF ${i}_dup.bam >${i}_raw.vcf

# ******************************************
# 7. bcftools filter 过滤变异
# ******************************************
echo -e "\n\n\033[34m Step 7 \033[0m: Filter variants with bcftools\n"

bcftools filter \
  --include 'QUAL>30 && INFO/DP>10' \
  --SnpGap 5 --IndelGap 10 \
  -O v -o ${i}_filter.vcf ${i}_raw.vcf

rm -rf $tmpdir
echo -e "\n\n\033[35m Pipeline completed successfully. \033[0m\n"
