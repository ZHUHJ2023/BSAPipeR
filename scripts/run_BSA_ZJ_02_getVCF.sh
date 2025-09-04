#!/bin/bash

usage() {
    echo "Usage: $0 -ref ref.fa -refbwa ref.fa(bwa index) -workdir workpath -sample_mut -sample_wt -datadir fqpath"
    echo "Options:"
    echo "  -ref ref.fa          reference fasta file"
    echo "  -refbwa ref.fa(bwa index)  reference fasta file(bwa index)"
    echo "  -workdir workpath    working directory"
    echo "  -sample_mut sample_mut_name  sample name of mutant"
    echo "  -sample_wt sample_wt_name    sample name of wild type"
    echo "  -depthFT False/True  whether to specify sequencing depth"
    echo "  -mdepth mdepth      target sequencing depth (>0)"
    echo "  -h, --help           show this help message and exit"
    exit 1
}


# Parameter parsing
while [ "$#" -gt 0 ]; do
    case "$1" in
        -ref) ref="$2" ; shift 2;;
        -refbwa) refbwa="$2" ; shift 2;;
        -workdir) workdir="$2" ; shift 2;;
        -sample_mut) sample_mut="$2" ; shift 2;;
        -sample_wt) sample_wt="$2" ; shift 2;;
        -depthFT) depthFT="$2" ; shift 2;;
        -mdepth) mdepth="$2" ; shift 2;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; exit 1;;
    esac
done

# Validate required parameters
if [ -z $ref ] || [ -z $refbwa ] || [ -z $workdir ] || [ -z $sample_mut ] || [ -z $sample_wt ] ; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi
# 日志函数 - 使用子shell隔离每个步骤
run_step_with_logging() {
    local step_name=$1
    local sample=$2
    shift 2
    local cmd=("$@")

    local log_dir="${workdir}/logs/${sample}"
    mkdir -p "${log_dir}"
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

# Validate depth parameters
if [ $depthFT == "True" ]; then
    if [ -z $mdepth ]; then
        echo "The sequencing depth is not specified.(mdepth)"
        usage
        exit 1
    fi
    if [ "$mdepth" == "0" ]; then
        echo "The sequencing depth is 0."
        usage
        exit 1
    fi
fi

# Create main log directory
if [ ! -d "${workdir}/logs" ]; then
    mkdir -p "${workdir}/logs"
fi

run_get_x_depth_bam() {
    local sample_mut=""
    local sample_wt=""
    local mdepth=""
    local workdir=""
    local threads=8

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --sample_mut)
                sample_mut="$2"
                shift 2
                ;;
            --sample_wt)
                sample_wt="$2"
                shift 2
                ;;
            --mdepth)
                mdepth="$2"
                shift 2
                ;;
            --workdir)
                workdir="$2"
                shift 2
                ;;
            --threads)
                threads="$2"
                shift 2
                ;;
            *)
                echo "Unknown parameter: $1"
                return 1
                ;;
        esac
    done

    if [[ -z $sample_mut || -z $sample_wt || -z $mdepth || -z $workdir ]]; then
        echo "Error: Missing required arguments."
        return 1
    fi

    # 创建深度计算目录
    if [ ! -d $workdir/03_gatk_depth ]; then
        mkdir -p $workdir/03_gatk_depth
    fi

    # 处理突变样本
    if [ ! -s $workdir/03_gatk_depth/${sample_mut}.mapped.bam.chr.stat.gz ]; then
        # 创建映射的BAM文件
        if [ ! -s $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam ]; then
            run_step_with_logging "create_mapped_bam_mut" "${sample_mut}" \
                $samtools view -@ 8 -b -F 4 $workdir/03_gatk/${sample_mut}/${sample_mut}.BQSR.bam > $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam
        fi

        # 计算深度
        run_step_with_logging "depth_calculation_mut" "${sample_mut}" \
            $pandepth -i $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam -t $threads -o $workdir/03_gatk_depth/${sample_mut}.mapped.bam
    fi

    # 计算降采样率
    local alldepth_mut=$(less $workdir/03_gatk_depth/${sample_mut}.mapped.bam.chr.stat.gz| awk 'NR>1 && ($1 ~ /^chr/ || $1 ~ /^[0-9]/) {sum += $NF; count++} END {if (count > 0) print sum/count}' )
    local p_mut=$( echo "scale=2; $mdepth/$alldepth_mut" |bc)
    p_mut=$(echo "0$p_mut")
    echo "${sample_mut} ${mdepth}x Downsample rate: $p_mut"

    # 执行降采样
    if [ ! -s $workdir/03_gatk/${sample_mut}/${mdepth}x.${sample_mut}.BQSR.bam ]; then
        run_step_with_logging "picard_downsample_mut_${mdepth}" "${sample_mut}" \
            bash -c "echo '${sample_mut} ${mdepth}x Downsample rate: $p_mut' && \
            $picard DownsampleSam \
                   --INPUT $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam \
                   --OUTPUT $workdir/03_gatk/${sample_mut}/${mdepth}x.${sample_mut}.BQSR.bam \
                   --PROBABILITY $p_mut \
                   --RANDOM_SEED $mdepth \
                   --CREATE_INDEX TRUE "
    fi

    # 处理野生型样本
    if [ ! -s $workdir/03_gatk_depth/${sample_wt}.mapped.bam.chr.stat.gz ]; then
        # 创建映射的BAM文件
        if [ ! -s $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam ]; then
            run_step_with_logging "create_mapped_bam_wt" "${sample_wt}" \
                bash -c "$samtools view -@ 8 -b -F 4 $workdir/03_gatk/${sample_wt}/${sample_wt}.BQSR.bam > $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam"
        fi

        # 计算深度
        run_step_with_logging "depth_calculation_wt" "${sample_wt}" \
            $pandepth -i $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam -t $threads -o $workdir/03_gatk_depth/${sample_wt}.mapped.bam
    fi

    # 计算降采样率
    local alldepth_wt=$(less $workdir/03_gatk_depth/${sample_wt}.mapped.bam.chr.stat.gz| awk 'NR>1 && ($1 ~ /^chr/ || $1 ~ /^[0-9]/) {sum += $NF; count++} END {if (count > 0) print sum/count}' )
    local p_wt=$( echo "scale=2; $mdepth/$alldepth_wt" |bc)
    p_wt=$(echo "0$p_wt")
    echo "${sample_wt} ${mdepth}x Downsample rate: $p_wt"

    # 执行降采样
    if [ ! -s $workdir/03_gatk/${sample_wt}/${mdepth}x.${sample_wt}.BQSR.bam ]; then
        run_step_with_logging "picard_downsample_wt_${mdepth}" "${sample_wt}" \
            bash -c "echo '${sample_wt} ${mdepth}x Downsample rate: $p_wt' && \
            $picard DownsampleSam \
                    --INPUT $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam \
                    --OUTPUT $workdir/03_gatk/${sample_wt}/${mdepth}x.${sample_wt}.BQSR.bam \
                    --PROBABILITY $p_wt \
                    --RANDOM_SEED $mdepth \
                    --CREATE_INDEX TRUE "
    fi
}

run_gatk_pipeline() {
    local sample_wt=$1
    local sample_mut=$2
    local ref=$3
    local workdir=$4
    local mdepth=${5:-0}

    local depth_prefix=""
    if [ "$mdepth" -gt 0 ]; then
        depth_prefix="${mdepth}x."
    fi

    # 处理每个样本的HaplotypeCaller
    process_sample() {
        local sample=$1
        local depth_prefix=$2

        # 运行HaplotypeCaller
        if [ ! -s "$workdir/03_gatk/${sample}/${depth_prefix}${sample}.MAPQ20.BQSR.g.vcf.gz" ]; then
            run_step_with_logging "${depth_prefix}haplotype_caller_${sample}" "${sample}" \
                bash -c "$samtools index '$workdir/03_gatk/${sample}/${depth_prefix}${sample}.BQSR.bam' && \
                $gatk HaplotypeCaller \
                    -R $ref \
                    -I $workdir/03_gatk/${sample}/${depth_prefix}${sample}.BQSR.bam \
                    -ERC GVCF \
                    -O $workdir/03_gatk/${sample}/${depth_prefix}${sample}.MAPQ20.BQSR.g.vcf.gz \
                    --native-pair-hmm-threads 20"
        fi
    }

    # 处理样本
    process_sample "$sample_wt" "$depth_prefix"
    process_sample "$sample_mut" "$depth_prefix"

    # 合并GVCFs
    if [ ! -d "$workdir/04_finalVCF/${sample_mut}" ]; then
        mkdir -p "$workdir/04_finalVCF/${sample_mut}"
    fi

    if [ ! -s "$workdir/04_finalVCF/${sample_mut}/${depth_prefix}gatk.${sample_mut}.${sample_wt}.genotypeGVCFs.vcf.gz" ]; then
        # 合并GVCF
        if [ ! -s "$workdir/03_gatk/${sample_mut}/${depth_prefix}gatk.${sample_mut}.${sample_wt}.combine.g.vcf.gz" ]; then
            run_step_with_logging "${depth_prefix}combine_gvcfs" "${sample_mut}" \
                $gatk CombineGVCFs \
                    -R "$ref" \
                    -V "$workdir/03_gatk/${sample_mut}/${depth_prefix}${sample_mut}.MAPQ20.BQSR.g.vcf.gz" \
                    -V "$workdir/03_gatk/${sample_wt}/${depth_prefix}${sample_wt}.MAPQ20.BQSR.g.vcf.gz" \
                    -O "$workdir/03_gatk/${sample_mut}/${depth_prefix}gatk.${sample_mut}.${sample_wt}.combine.g.vcf.gz"
        fi

        # 基因分型
        run_step_with_logging "${depth_prefix}genotype_gvcfs" "${sample_mut}" \
            $gatk GenotypeGVCFs \
                -R "$ref" \
                -V "$workdir/03_gatk/${sample_mut}/${depth_prefix}gatk.${sample_mut}.${sample_wt}.combine.g.vcf.gz" \
                -O "$workdir/04_finalVCF/${sample_mut}/${depth_prefix}gatk.${sample_mut}.${sample_wt}.genotypeGVCFs.vcf.gz"
    fi

    # VQSR处理
    vqsr_process() {
        local mode=$1
        local depth_prefix=$2

        local input_vcf="$workdir/04_finalVCF/${sample_mut}/${depth_prefix}gatk.${sample_mut}.${sample_wt}.genotypeGVCFs.vcf.gz"
        local output_prefix="$workdir/04_finalVCF/${sample_mut}/${depth_prefix}${sample_mut}"

        if [ "$mode" == "SNP" ] && [ ! -s "${output_prefix}.snp.VQSR.vcf.gz" ]; then
            run_step_with_logging "${depth_prefix}vqsr_snp" "${sample_mut}" \
                bash -c "cd '$workdir/04_finalVCF/${sample_mut}' && \
                $gatk VariantRecalibrator \
                    -R '$ref' \
                    --variant '$input_vcf' \
                    -resource:samtools_bcftools,training=true,truth=true \
                    '$workdir/02_vcf/merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf' \
                    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS \
                    -an SOR -an DP -tranche 100.0 -tranche 99.9 \
                    -tranche 99.0 -tranche 95.0 -tranche 90.0 \
                    -mode SNP -O '${output_prefix}.snp.recal' \
                    --tranches-file '${output_prefix}.snp.tranches' \
                    --rscript-file '${output_prefix}.snp.plots.R'
                sed -i 's/space=\"rgb\"/space=\"Lab\"/g' '${output_prefix}.snp.plots.R'
                Rscript '${output_prefix}.snp.plots.R'
                $gatk ApplyVQSR \
                    -R '$ref' \
                    --variant '$input_vcf' \
                    --tranches-file '${output_prefix}.snp.tranches' \
                    --recal-file '${output_prefix}.snp.recal' \
                    -mode SNP -O '${output_prefix}.snp.VQSR.vcf.gz'"
        fi

        if [ "$mode" == "INDEL" ] && [ ! -s "${output_prefix}.VQSR.vcf.gz" ]; then
            run_step_with_logging "${depth_prefix}vqsr_indel" "${sample_mut}" \
                bash -c "cd '$workdir/04_finalVCF/${sample_mut}' && \
                $gatk VariantRecalibrator \
                    -R '$ref' \
                    --variant '${output_prefix}.snp.VQSR.vcf.gz' \
                    --max-gaussians 6 \
                    -resource:samtools_bcftools,training=true,truth=true \
                    '$workdir/02_vcf/merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf' \
                    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
                    --tranches-file '${output_prefix}.snp.indel.tranches' \
                    --rscript-file '${output_prefix}.snp.indel.plots.R' \
                    -mode INDEL -O '${output_prefix}.snp.indel.recal'
                sed -i \"s/space=\"rgb\"/space=\"Lab\"/g\" '${output_prefix}.snp.indel.plots.R'
                Rscript '${output_prefix}.snp.indel.plots.R'
                $gatk ApplyVQSR \
                    -R $ref \
                    --variant '${output_prefix}.snp.VQSR.vcf.gz' \
                    --tranches-file '${output_prefix}.snp.indel.tranches' \
                    --recal-file '${output_prefix}.snp.indel.recal' \
                    -mode INDEL -O '${output_prefix}.VQSR.vcf.gz' "
        fi
    }

    vqsr_process "SNP" "$depth_prefix"
    vqsr_process "INDEL" "$depth_prefix"

    # 过滤变异
    filter_variants() {
        local variant_type=$1
        local depth_prefix=$2

        local input_file=""
        local output_file=""

        if [ "$variant_type" == "snp" ]; then
            input_file="$workdir/04_finalVCF/${sample_mut}/${depth_prefix}${sample_mut}.snp.VQSR.vcf.gz"
            output_file="$workdir/05_snpIndel/${depth_prefix}${sample_mut}.snp.VQSR.vcf.gz"
        else
            input_file="$workdir/04_finalVCF/${sample_mut}/${depth_prefix}${sample_mut}.VQSR.vcf.gz"
            output_file="$workdir/05_snpIndel/${depth_prefix}${sample_mut}.indel.VQSR.vcf.gz"
        fi

        if [ ! -d "$workdir/05_snpIndel" ]; then
            mkdir -p "$workdir/05_snpIndel"
        fi

        if [ ! -s "$output_file" ]; then
            if [ "$variant_type" == "snp" ]; then
                cmd="$bcftools filter -e 'QD < 2.0 || MQ < 40.0 || QUAL < 30' \"$input_file\" | \
                     $bcftools view -m2 -M2 -v snps - | gzip > \"$output_file\""
            else
                cmd="$bcftools filter -e 'QD < 2.0 || MQ < 40.0 || QUAL < 30' \"$input_file\" | \
                     $bcftools view -m2 -M2 -v indels - | gzip > \"$output_file\""
            fi

            run_step_with_logging "${depth_prefix}filter_${variant_type}" "${sample_mut}" \
                bash -c "$cmd"
        fi
        # if [ ! -s "$output_file" ]; then
        #     run_step_with_logging "${depth_prefix}filter_${variant_type}" "${sample_mut}" \
        #         if [ "$variant_type" == "snp" ]; then \
        #             $bcftools filter -e 'QD < 2.0 || MQ < 40.0 || QUAL < 30' "$input_file" | \
        #             $bcftools view -m2 -M2 -v snps - | gzip > "$output_file"; \
        #         else \
        #             $bcftools filter -e 'QD < 2.0 || MQ < 40.0 || QUAL < 30' "$input_file" | \
        #             $bcftools view -m2 -M2 -v indels - | gzip > "$output_file"; \
        #         fi
        # fi
    }

    filter_variants "snp" "$depth_prefix"
    filter_variants "indel" "$depth_prefix"

    echo "GATK pipeline completed for ${sample_wt} and ${sample_mut} with depth ${mdepth}x."
}

# 主执行部分
if [ $depthFT == "True" ]; then
    # 创建深度计算目录
    if [ ! -d $workdir/03_gatk_depth ]; then
        mkdir -p $workdir/03_gatk_depth
    fi

    # 计算突变样本深度
    if [ ! -s $workdir/03_gatk_depth/${sample_mut}.mapped.bam.chr.stat.gz ]; then
        # 创建映射的BAM文件
        if [ ! -s $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam ]; then
            run_step_with_logging "create_mapped_bam_mut" "${sample_mut}" \
                bash -c "$samtools view -@ 8 -b -F 4 $workdir/03_gatk/${sample_mut}/${sample_mut}.BQSR.bam > $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam"
        fi

        # 计算深度
        run_step_with_logging "initial_depth_mut" "${sample_mut}" \
            $pandepth -i $workdir/03_gatk/${sample_mut}/${sample_mut}.mapped.bam -t 8 -o $workdir/03_gatk_depth/${sample_mut}.mapped.bam
    fi

    # 计算野生型样本深度
    if [ ! -s $workdir/03_gatk_depth/${sample_wt}.mapped.bam.chr.stat.gz ]; then
        # 创建映射的BAM文件
        if [ ! -s $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam ]; then
            run_step_with_logging "create_mapped_bam_wt" "${sample_wt}" \
                bash -c "$samtools view -@ 8 -b -F 4 $workdir/03_gatk/${sample_wt}/${sample_wt}.BQSR.bam > $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam"
        fi

        # 计算深度
        run_step_with_logging "initial_depth_wt" "${sample_wt}" \
            $pandepth -i $workdir/03_gatk/${sample_wt}/${sample_wt}.mapped.bam -t 8 -o $workdir/03_gatk_depth/${sample_wt}.mapped.bam
    fi

    # 验证深度要求
    alldepth_mut=$(less $workdir/03_gatk_depth/${sample_mut}.mapped.bam.chr.stat.gz| awk 'NR>1 && ($1 ~ /^chr/ || $1 ~ /^[0-9]/) {sum += $NF; count++} END {if (count > 0) print sum/count}' )
    alldepth_wt=$(less $workdir/03_gatk_depth/${sample_wt}.mapped.bam.chr.stat.gz| awk 'NR>1 && ($1 ~ /^chr/ || $1 ~ /^[0-9]/) {sum += $NF; count++} END {if (count > 0) print sum/count}' )

    if (( $(echo "$mdepth > $alldepth_mut" | bc -l) )) || (( $(echo "$mdepth > $alldepth_wt" | bc -l) )); then
        echo "Error: mdepth is larger than total sequencing depth."
        exit 1
    fi

    # 执行降采样和GATK流程
    run_step_with_logging "run_get_x_depth_bam" "${sample_mut}_${sample_wt}_${mdepth}" \
        run_get_x_depth_bam --sample_mut $sample_mut --sample_wt $sample_wt --mdepth ${mdepth} --workdir $workdir

    run_step_with_logging "run_gatk_pipeline" "${sample_mut}_${sample_wt}_${mdepth}" \
        run_gatk_pipeline $sample_wt $sample_mut $ref $workdir $mdepth
elif [ $depthFT == "False" ]; then
    # 直接执行GATK流程
    run_step_with_logging "run_gatk_pipeline" "${sample_mut}_${sample_wt}" \
        run_gatk_pipeline $sample_wt $sample_mut $ref $workdir
fi
