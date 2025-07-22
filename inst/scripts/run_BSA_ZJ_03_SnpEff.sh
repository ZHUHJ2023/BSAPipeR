#!/bin/bash


usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -sample_mut    Sample name for mutation calling (required)"
    echo "  -workdir      Working directory (default: current directory)"
    echo "  -GFF3          Reference genome in GFF3 format (full path)"
    echo "  -FA            Reference genome in FASTA format (full path)"
    echo "  -snpEff_path   Path to snpEff directory (optional)"
    echo "  -mdepth        depth of coverage for mutation calling (default: 0)"
    echo "  -h, --help     Show this help message"
    exit 1
}

# 参数解析
while [[ $# -gt 0 ]]; do
    case "$1" in
        -sample_mut) sample_mut="$2"; shift 2 ;;
        -workdir) workdir="$2"; shift 2 ;;
        -GFF3) GFF3="$2"; shift 2 ;;
        -FA) FA="$2"; shift 2 ;;
        -snpEff_path) snpEff_path="$2"; shift 2 ;;
        -mdepth) mdepth="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# 参数验证
[[ -z "$sample_mut" ]] && { echo "Error: Sample name for mutation calling is required" >&2; usage; exit 1; }
[[ -z "$GFF3" ]] && { echo "Error: GFF3 file is required" >&2; usage; exit 1; }
[[ -z "$FA" ]] && { echo "Error: FASTA file is required" >&2; usage; exit 1; }
[[ -z "$snpEff_path" ]] && { echo "snpEff path not specified" >&2; usage; exit 1; }
[[ -z "$workdir" ]] && workdir=$(pwd)  # 如果未指定工作目录，则使用当前目录

: ${mdepth:=0}  # 设置默认深度值

# 创建主日志目录
if [ ! -d "${workdir}/logs/snpEff/${sample_mut}" ]; then
    mkdir -p "${workdir}/logs/snpEff/${sample_mut}"
fi

# 日志函数 - 使用子shell隔离每个步骤的日志
run_with_logging() {
    local step_name=$1
    local sample=$2
    shift 2
    local cmd=("$@")

    local log_dir="${workdir}/logs/snpEff/${sample}"
    if [ ! -d "$log_dir" ]; then
        mkdir -p "$log_dir"
    fi
    local LOG_FILE="${log_dir}/${step_name}.log"
    local ERR_FILE="${log_dir}/${step_name}.err"

    # 在子shell中执行命令并重定向输出
    (
        echo "===== STARTING STEP: $step_name ====="
        echo "Command: ${cmd[*]}"
        echo "Timestamp: $(date)"

        # 执行实际命令
        "${cmd[@]}"

        local exit_code=$?
        echo "Exit code: $exit_code"
        echo "Timestamp: $(date)"
        echo "===== FINISHED STEP: $step_name ====="
        exit $exit_code
    ) > "$LOG_FILE" 2> "$ERR_FILE"

    return $?
}

#######################################
# 获取snpEff路径
# 输出:
#   设置全局变量SNP_EFF_DIR
# 返回:
#   0-成功 1-失败
#######################################

SNP_EFF_DIR="$snpEff_path"
# 验证路径是否存在
if [[ ! -d "$SNP_EFF_DIR" ]]; then
    echo "Error: snpEff directory not found at $SNP_EFF_DIR" >&2
    return 1
fi

# 验证必要的文件是否存在
if [[ ! -f "$SNP_EFF_DIR/snpEff.jar" ]]; then
    echo "Error: snpEff.jar not found in $SNP_EFF_DIR" >&2
    return 1
fi

#######################################
# 构建snpEff数据库
#######################################
build_snpEff_database() {
    local GFF3="$1"
    local FA="$2"
    local genome_version
    local gffname=$(basename "$GFF3")
    genome_version=$(echo "$gffname" | awk -F. '{print $1}')
    local data_dir="$SNP_EFF_DIR/data/$genome_version"

    # 检查数据库是否已存在
    if [[ -f "$data_dir/snpEffectPredictor.bin" ]]; then
        echo "snpEff database for $genome_version already exists at $data_dir"
        return 0
    fi

    echo "Building snpEff database for $genome_version..."

    # 准备数据库目录
    if [ ! -d "$data_dir" ]; then
        echo "Creating data directory: $data_dir"
        mkdir -p "$data_dir" || return 1
    fi


    # 准备GFF文件(包含FASTA)
    echo "Preparing GFF+FASTA file..."
    local combined_gff="$SNP_EFF_DIR/data/${genome_version}.gff3"
    cp "$GFF3" "$combined_gff" || return 1

    echo "##FASTA" >> "$combined_gff"
    cat "$FA" >> "$combined_gff" || return 1

    # 更新snpEff配置文件
    echo "Updating snpEff config..."
    local config_file="$SNP_EFF_DIR/snpEff.config"

    # 检查是否已存在该基因组配置
    if ! grep -q "^$genome_version.genome" "$config_file"; then
        echo "$genome_version.genome : $genome_version" >> "$config_file"
        echo "data_dir=${data_dir}" >> "$config_file"
    fi

    # 移动文件到基因组目录
    mv "$combined_gff" "$data_dir/genes.gff" || return 1

    # 构建数据库
    echo "Building database (this may take a while)..."
    (cd "$SNP_EFF_DIR" && \
     $java -XX:ParallelGCThreads=1 -Xmx4g -jar snpEff.jar \
        build -gff3 -v "$genome_version") || return 1

    echo "snpEff database for $genome_version built successfully at $data_dir"
}

#######################################
# 使用snpEff进行变异注释
#######################################
run_snpEff_annotation() {
    local GFF3="$1"
    local FA="$2"
    local VCF="$3"
    local outputdir="$4"
    local SNP_EFF_DIR="$5"

    # 获取基因组版本名
    local gffname=$(basename "$GFF3")
    local genome_version=$(echo "$gffname" | awk -F. '{print $1}')
    local data_dir="$SNP_EFF_DIR/data/$genome_version"

    # 检查数据库是否存在，不存在则构建
    if [[ ! -f "$data_dir/snpEffectPredictor.bin" ]]; then
        echo "snpEff database not found, building now..."
        build_snpEff_database "$GFF3" "$FA" || return 1
    fi

    # 准备输出目录
    if [ ! -d "$outputdir" ]; then
        mkdir -p "$outputdir" || return 1
    fi


    # 准备输出文件名
    local vcfname=$(basename "$VCF")
    local outputname="$(echo $vcfname | sed 's/vcf/snpEff.vcf/g')"
    local statsname="$(echo $vcfname | sed 's/vcf.gz/snpEff.html/g')"

    # 运行注释
    echo "Annotating $VCF using snpEff..."
    (cd "$SNP_EFF_DIR" && \
     $java -XX:ParallelGCThreads=4 -Xmx8g -jar snpEff.jar \
        ann -c snpEff.config \
        "$genome_version" -i vcf -o vcf "$VCF" \
        > "$outputdir/$outputname" \
        -stats "$outputdir/$statsname") || return 1

    echo "Annotation completed successfully!"
    echo "Output files:"
    echo "  - Annotated VCF: $outputdir/$outputname"
    echo "  - Stats report: $outputdir/$statsname"
}

#######################################
# 主snpEff流程
#######################################
run_snpeff_pipeline() {
    local sample_mut="$1"
    local gff="$2"
    local ref="$3"
    local workdir="$4"
    local SNP_EFF_DIR="$5"
    local mdepth="${6:-0}"

    local depth_prefix=""
    [[ "$mdepth" -gt 0 ]] && depth_prefix="${mdepth}x."

    # 创建输出目录
    local output_dir="$workdir/06_SnpEff/$sample_mut"
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
    fi

    # 处理SNP注释
    local snp_vcf="$workdir/05_snpIndel/${depth_prefix}${sample_mut}.snp.VQSR.vcf.gz"
    if [[ -f "$snp_vcf" ]]; then
        run_with_logging "${depth_prefix}snp_annotation" "$sample_mut" \
            run_snpEff_annotation "$gff" "$ref" "$snp_vcf" "$output_dir" "$SNP_EFF_DIR" || return 1
    else
        echo "Warning: SNP VCF file not found: $snp_vcf" >&2
    fi

    # 处理Indel注释
    local indel_vcf="$workdir/05_snpIndel/${depth_prefix}${sample_mut}.indel.VQSR.vcf.gz"
    if [[ -f "$indel_vcf" ]]; then
        run_with_logging "${depth_prefix}indel_annotation" "$sample_mut" \
            run_snpEff_annotation "$gff" "$ref" "$indel_vcf" "$output_dir" "$SNP_EFF_DIR" || return 1
    else
        echo "Warning: Indel VCF file not found: $indel_vcf" >&2
    fi
}

# 主执行部分
run_with_logging "main" "$sample_mut" \
    run_snpeff_pipeline "$sample_mut" "$GFF3" "$FA" "$workdir" "$snpEff_path" "$mdepth"
