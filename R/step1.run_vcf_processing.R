#' run_vcf_processing
#'
#' @param ref
#' @param refbwa
#' @param outdir
#' @param sample_mut
#' @param sample_wt
#' @param mutfq1
#' @param mutfq2
#' @param wtfq1
#' @param wtfq2
#' @param depthFT
#' @param mdepth
#'
#' @returns
#' @export
#'
#' @examples set_software_paths(bwa = "/path/to/bwa",fastp = "/path/to/fastp",samtools = "/path/to/samtools",picard = "/path/to/picard.jar",bcftools = "/path/to/bcftools",gatk = "/path/to/gatk.jar",snpEff = "/path/to/snpEff") run_vcf_processing(ref = "/path/to/reference.fa",refbwa = "/path/to/bwa_index_prefix",workdir = "/path/to/workdir",sample_mut = "MUT",sample_wt = "WT",mutfq1 = "/path/to/mut_1.fq.gz",mutfq2 = "/path/to/mut_2.fq.gz",wtfq1 = "/path/to/wt_1.fq.gz",wtfq2 = "/path/to/wt_2.fq.gz")
run_vcf_processing <- function(ref, refbwa, outdir, sample_mut, sample_wt,GFF3,
                               mutfq1, mutfq2, wtfq1, wtfq2, depthFT = "False", mdepth = 0) {

# #软件绝对路径
# bwa=/data/software/bwa-0.7.17/bwa
# fastp=/data/software/fastp
# samtools=/data/software/samtools-1.17/bin/samtools
# picard="/data/software/jdk-18.0.2.1/bin/java -XX:ParallelGCThreads=8 -Xmx30g -jar /data/home/hjzhustu/home/software/picard-3.1.1/picard.jar"
# bcftools=/data/home/hjzhustu/anaconda3/envs/BSAtools/bin/bcftools
# gatk="/data/software/jdk-18.0.2.1/bin/java -XX:ParallelGCThreads=8 -Xmx30g -jar /data/software/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"  #使用了java -jar 所以gatk后面的--java-options '-XX:ParallelGCThreads=20 -Xmx100G'需要删掉
# snpEff_path="/data/home/hjzhustu/home/software/snpEff-5.0"
 # 获取配置的软件路径
  software_paths <- get_software_paths()

  # 设置环境变量传递给shell脚本
  env_vars <- c(
    paste0("bwa=", software_paths$bwa),
    paste0("fastp=", software_paths$fastp),
    paste0("samtools=", software_paths$samtools),
    # 对于包含空格的命令，使用引号包裹
    paste0("picard='", software_paths$picard, "'"),
    paste0("bcftools=", software_paths$bcftools),
    paste0("gatk='", software_paths$gatk, "'"),
    paste0("snpEff_path=", software_paths$snpEff),
    paste0("pandepth=", software_paths$pandepth),
    paste0("java=", software_paths$java)
  )


  workdir <- paste0(outdir, "/Processing")
  # 获取脚本路径
  shell_01script <- system.file("scripts/run_BSA_ZJ_01_getBAM.sh", package = "BSAliulab")
  if (!file.exists(shell_01script)) {
    stop("Shell script not found in package resources")
  }

  # 执行第一个脚本
  system(paste(paste(env_vars, collapse = " "), "bash", shell_01script,
               "-ref", shQuote(ref),
               "-refbwa", shQuote(refbwa),
               "-workdir", shQuote(workdir),
               "-sample_mut", shQuote(sample_mut),
               "-sample_wt", shQuote(sample_wt),
               "-mutfq1", shQuote(mutfq1),
               "-mutfq2", shQuote(mutfq2),
               "-wtfq1", shQuote(wtfq1),
               "-wtfq2", shQuote(wtfq2)))

  # 第二个脚本
  shell_02script <- system.file("scripts/run_BSA_ZJ_02_getVCF.sh", package = "BSAliulab")
  if (!file.exists(shell_02script)) {
    stop("Shell script not found in package resources")
  }

  if (depthFT == "True") {
    system(paste(paste(env_vars, collapse = " "), "bash", shell_02script,
                 "-ref", shQuote(ref),
                 "-refbwa", shQuote(refbwa),
                 "-workdir", shQuote(workdir),
                 "-sample_mut", shQuote(sample_mut),
                 "-sample_wt", shQuote(sample_wt),
                 "-depthFT", depthFT,
                 "-mdepth", mdepth))
  } else {
    system(paste(paste(env_vars, collapse = " "), "bash", shell_02script,
                 "-ref", shQuote(ref),
                 "-refbwa", shQuote(refbwa),
                 "-workdir", shQuote(workdir),
                 "-sample_mut", shQuote(sample_mut),
                 "-sample_wt", shQuote(sample_wt),
                 "-depthFT", depthFT))
  }

  # 第三个脚本
  shell_03script <- system.file("scripts/run_BSA_ZJ_03_SnpEff.sh", package = "BSAliulab")
  if (!file.exists(shell_03script)) {
    stop("Shell script not found in package resources")
  }

  if (depthFT == "True") {
    system(paste(paste(env_vars, collapse = " "), "bash", shell_03script,
                 "-sample_mut", shQuote(sample_mut),
                 "-workdir", shQuote(workdir),
                 "-GFF3", shQuote(GFF3),
                 "-FA", shQuote(ref),
                 "-snpEff_path", shQuote(software_paths$snpEff),
                 "-mdepth", mdepth))
  } else {
    system(paste(paste(env_vars, collapse = " "), "bash", shell_03script,
                 "-sample_mut", shQuote(sample_mut),
                 "-workdir", shQuote(workdir),
                 "-GFF3", shQuote(GFF3),
                 "-FA", shQuote(ref),
                 "-snpEff_path", shQuote(software_paths$snpEff)))
  }
}
