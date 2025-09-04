#!/bin/bash

usage() {
    echo "Usage: $0 -ref ref.fa -refbwa ref.fa(bwa index) -workdir workpath -sample_mut -sample_wt -datadir fqpath"
    echo "Options:"
    echo "  -ref ref.fa          reference fasta file"
    echo "  -refbwa ref.fa(bwa index)  reference fasta file(bwa index)"
    echo "  -workdir workpath    working directory"
    echo "  -sample_mut sample_mut_name  sample name of mutant"
    echo "  -sample_wt sample_wt_name    sample name of wild type"
    echo "  -mutfq1               mutant fastq file 1"
    echo "  -mutfq2               mutant fastq file 2"
    echo "  -wtfq1                wild type fastq file 1"
    echo "  -wtfq2                wild type fastq file 2"
    echo "  -h, --help           show this help message and exit"
    exit 1
}


# Parse parameters
while [ "$#" -gt 0 ]; do
    case "$1" in
        -ref) ref="$2" ; shift 2;;
        -refbwa) refbwa="$2" ; shift 2;;
        -workdir) workdir="$2" ; shift 2;;
        -sample_mut) sample_mut="$2" ; shift 2;;
        -sample_wt) sample_wt="$2" ; shift 2;;
        -mutfq1) mutfq1="$2" ; shift 2;;
        -mutfq2) mutfq2="$2" ; shift 2;;
        -wtfq1) wtfq1="$2" ; shift 2;;
        -wtfq2) wtfq2="$2" ; shift 2;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; exit 1;;
    esac
done

# 日志函数 - 使用子shell隔离每个步骤
run_step_with_logging() {
    local step_name=$1
    local sample=$2
    shift 2
    local cmd=("$@")

    local log_dir="${workdir}/logs/${sample}"
    if [ ! -d "${log_dir}" ]; then
        mkdir -p "${log_dir}"
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
# Create main log directory
if [ ! -d "${workdir}/logs" ]; then
    mkdir -p "${workdir}/logs"
fi

# 01. bwa index for ref
if [ ! -s $refbwa.bwt ]; then
    run_step_with_logging "bwa_index" "reference" \
        $bwa index $refbwa
fi
samplelist=($sample_mut $sample_wt)
for sample in ${samplelist[@]}
do
    # 02. fastp.QC
    if [ ! -s $workdir/03_gatk/${sample}/${sample}.BQSR.bam ]; then
        if [ ! -d ${workdir} ]; then mkdir -p ${workdir} ;fi
        cd ${workdir}
        fpout=${workdir}/00_fastp
        mkdir -p $fpout/data/${sample}
        mkdir $fpout/outfile
        if [ $sample == $sample_mut ];then
            fq1=$mutfq1
            fq2=$mutfq2
        else
            fq1=$wtfq1
            fq2=$wtfq2
        fi
        if [ ! -s $fpout/outfile/${sample}.html ];then
            run_step_with_logging "fastp_QC" "${sample}" \
                $fastp -q 20 -u 50 -l 40 -n 7 -e 20 -3 \
                --thread 8 \
                -i $fq1 \
                -I $fq2 \
                -o $fpout/data/${sample}/filter.${sample}.1.fq.gz \
                -O $fpout/data/${sample}/filter.${sample}.2.fq.gz \
                -h $fpout/outfile/${sample}.html \
                -j $fpout/outfile/${sample}.json
        fi

        # 03. BWA alignment
        if [ ! -s $workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.dup.bam ]; then
            if [ ! -d ${workdir}/01_bwa/${sample} ]; then mkdir -p ${workdir}/01_bwa/${sample} ;fi
            if [ ! -s $workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.bam ];then
                # fq1
                run_step_with_logging "bwa_aln_1" "${sample}" \
                    bash -c "$bwa aln -t 8 $refbwa $workdir/00_fastp/data/${sample}/filter.${sample}.1.fq.gz > $workdir/01_bwa/${sample}/${sample}.1.sai "

                # fq2
                run_step_with_logging "bwa_aln_2" "${sample}" \
                    bash -c "$bwa aln -t 8 $refbwa $workdir/00_fastp/data/${sample}/filter.${sample}.2.fq.gz > $workdir/01_bwa/${sample}/${sample}.2.sai "

                # BWA pe
                run_step_with_logging "bwa_sampe" "${sample}" \
                    bash -c "$bwa sampe -a 1000 -r \"@RG\tID:${sample}\tPL:DNBseq\tLB:${sample}\tSM:${sample}\" $refbwa $workdir/01_bwa/${sample}/${sample}.1.sai $workdir/01_bwa/${sample}/${sample}.2.sai $workdir/00_fastp/data/${sample}/filter.${sample}.1.fq.gz $workdir/00_fastp/data/${sample}/filter.${sample}.2.fq.gz > $workdir/01_bwa/${sample}/${sample}.sam "
                # sam to bam
                run_step_with_logging "sam_to_bam" "${sample}" \
                    bash -c "\
                    $samtools view -@ 8 -bS $workdir/01_bwa/${sample}/${sample}.sam > $workdir/01_bwa/${sample}/${sample}.bam && \
                    $samtools sort -@ 8 $workdir/01_bwa/${sample}/${sample}.bam -o $workdir/01_bwa/${sample}/${sample}.sorted.bam && \
                    $samtools view -@ 8 -q 20 -b $workdir/01_bwa/${sample}/${sample}.sorted.bam > $workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.bam && \
                    rm -f $workdir/01_bwa/${sample}/${sample}.sam $workdir/01_bwa/${sample}/${sample}.bam $workdir/01_bwa/${sample}/${sample}.sorted.bam $workdir/01_bwa/${sample}/${sample}.1.sai $workdir/01_bwa/${sample}/${sample}.2.sai \
                    "
            fi

            # samtools index & picard markDuplicate
            run_step_with_logging "samtools_index" "${sample}" \
                $samtools index -@ 8 $workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.bam

            run_step_with_logging "picard_mark_duplicates" "${sample}" \
                $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT \
                I=$workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.bam \
                O=$workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.dup.bam \
                M=$workdir/01_bwa/${sample}/${sample}.MAPQ20.dup.matrics.txt && \
                rm -f $workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.bam
        fi

        # bcftools mpileup + call SNP
        if [ ! -d ${workdir}/02_vcf ]; then mkdir -p ${workdir}/02_vcf ;fi
        if [ ! -s $workdir/02_vcf/${sample}.raw.bedtools.vcf.gz ];then
            run_step_with_logging "bcftools_mpileup_call" "${sample}" \
                bash -c "$bcftools mpileup --threads 8 --max-depth 100000 -E -C50 -Q20 -q20 -Ou -f $ref \
                $workdir/01_bwa/${sample}/${sample}.MAPQ20.sorted.dup.bam | \
                $bcftools call --threads 8 -mv -Oz - > $workdir/02_vcf/${sample}.raw.bedtools.vcf.gz"
        fi

        # combine vcf files
        cd $workdir/02_vcf
        if [ ! -s sorted.${sample}.raw.bedtools.vcf.gz.tbi ];then
            run_step_with_logging "bcftools_sort_vcf" "${sample}" \
                $bcftools sort -Oz -o sorted.${sample}.raw.bedtools.vcf.gz ${sample}.raw.bedtools.vcf.gz && \
                $bcftools index sorted.${sample}.raw.bedtools.vcf.gz
        fi
    fi
done

# Merge vcf for each mutant line
cd $workdir/02_vcf
if [ ! -s merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf ];then
    run_step_with_logging "bcftools_index_mut" "${sample_mut}" \
        $bcftools index sorted.${sample_mut}.raw.bedtools.vcf.gz

    run_step_with_logging "bcftools_index_wt" "${sample_wt}" \
        $bcftools index sorted.${sample_wt}.raw.bedtools.vcf.gz

    run_step_with_logging "bcftools_merge_vcf" "${sample_mut}_${sample_wt}" \
        $bcftools merge sorted.${sample_mut}.raw.bedtools.vcf.gz sorted.${sample_wt}.raw.bedtools.vcf.gz \
        -Oz -o merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf.gz && \
        gunzip merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf.gz
fi

run_step_with_logging "gatk_create_dict" "reference" \
    $gatk CreateSequenceDictionary -R $ref

# 06. Use gatk bqsr module to generate xx.BQSR.bam
run_step_with_logging "gatk_index_vcf_wt" "${sample_wt}" \
    $gatk IndexFeatureFile -I $workdir/02_vcf/sorted.${sample_wt}.raw.bedtools.vcf.gz
# create data.table for WT
if [ ! -d ${workdir}/03_gatk/${sample_wt} ]; then mkdir -p ${workdir}/03_gatk/${sample_wt} ;fi
if [ ! -s $workdir/03_gatk/${sample_wt}/${sample_wt}.data.table ];then
    run_step_with_logging "gatk_bqsr_wt" "${sample_wt}" \
        $gatk BaseRecalibrator \
        -R $ref \
        -I $workdir/01_bwa/${sample_wt}/${sample_wt}.MAPQ20.sorted.dup.bam \
        --known-sites $workdir/02_vcf/sorted.${sample_wt}.raw.bedtools.vcf.gz \
        -O $workdir/03_gatk/${sample_wt}/${sample_wt}.data.table
fi

if [ ! -s $workdir/03_gatk/${sample_wt}/${sample_wt}.BQSR.bam ];then
    run_step_with_logging "gatk_apply_bqsr_wt" "${sample_wt}" \
        $gatk ApplyBQSR \
        -R $ref \
        -I $workdir/01_bwa/${sample_wt}/${sample_wt}.MAPQ20.sorted.dup.bam \
        -bqsr $workdir/03_gatk/${sample_wt}/${sample_wt}.data.table \
        -O $workdir/03_gatk/${sample_wt}/${sample_wt}.BQSR.bam
fi

# For mutant sample
run_step_with_logging "gatk_index_vcf_mut" "${sample_mut}" \
    $gatk IndexFeatureFile -I $workdir/02_vcf/merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf

# create data.table for Mut
if [ ! -d ${workdir}/03_gatk/${sample_mut} ]; then mkdir -p ${workdir}/03_gatk/${sample_mut} ;fi
if [ ! -s $workdir/03_gatk/${sample_mut}/${sample_mut}.data.table ];then
    run_step_with_logging "gatk_bqsr_mut" "${sample_mut}" \
        $gatk BaseRecalibrator \
        -R $ref \
        -I $workdir/01_bwa/${sample_mut}/${sample_mut}.MAPQ20.sorted.dup.bam \
        --known-sites $workdir/02_vcf/merged.${sample_mut}.${sample_wt}.raw.bedtools.vcf \
        -O $workdir/03_gatk/${sample_mut}/${sample_mut}.data.table
fi

if [ ! -s $workdir/03_gatk/${sample_mut}/${sample_mut}.BQSR.bam ]; then
    run_step_with_logging "gatk_apply_bqsr_mut" "${sample_mut}" \
        $gatk ApplyBQSR \
        -R $ref \
        -I $workdir/01_bwa/${sample_mut}/${sample_mut}.MAPQ20.sorted.dup.bam \
        -bqsr $workdir/03_gatk/${sample_mut}/${sample_mut}.data.table \
        -O $workdir/03_gatk/${sample_mut}/${sample_mut}.BQSR.bam
fi
