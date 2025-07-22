#' Bulked Segregant Analysis (BSA) Pipeline
#'
#' 该函数实现完整的 BSA 分析流程，包括原始数据管理、变异检测和结果生成。
#' 主要步骤：1) 创建原始数据链接文件 2) 运行变异检测流程 3) 生成结果。
#'
#' @param ref 参考基因组 FASTA 文件路径
#' @param refbwa BWA 索引参考基因组路径(与 `ref` 相同,如果没有创建索引会自动创建，如果有建好的索引可以直接写建好索引的路径：path/ref.fa)
#' @param outdir 输出结果的主目录路径
#' @param sample_mut 突变型样本名称(用于结果标识)
#' @param sample_wt 野生型样本名称(用于结果标识)
#' @param GFF3 基因组注释文件 (GFF3格式) 路径
#' @param mutfq1 突变样本的 FASTQ 文件路径(read 1)
#' @param mutfq2 突变样本的 FASTQ 文件路径(read 2)
#' @param wtfq1 野生型样本的 FASTQ 文件路径(read 1)
#' @param wtfq2 野生型样本的 FASTQ 文件路径(read 2)
#' @param depthFT 是否启用深度过滤(默认 "False")，可选 "True" 启用
#' @param mdepth 最小深度阈值(默认 0，需指定不为0的值)，当 depthFT="True" 时生效
#' @param GT_mut 目标突变基因型(默认 "aa")，如 "aa" (纯合突变), "Aa" (杂合)
#' @param mindp 最小深度阈值(默认 4)，用于过滤变异位点
#'
#' @return 无直接返回值。结果文件将写入 `outdir` 目录，包括：
#' * rawdata/: 原始 FASTQ 链接文件
#' * Processing/: 中间比对文件以及VCF结果
#' * Results/: 最终结果
#'
#' @export
#'
#' @examples
#' \dontrun{
#' BSAPipeR(
#'   ref = "ref_genome.fasta",
#'   refbwa = "ref_genome",
#'   outdir = "Project_sample_mut/",
#'   sample_mut = "mutant",
#'   sample_wt = "wildtype",
#'   GFF3 = "annotation.gff3",
#'   mutfq1 = "mut_R1.fq.gz",
#'   mutfq2 = "mut_R2.fq.gz",
#'   wtfq1 = "wt_R1.fq.gz",
#'   wtfq2 = "wt_R2.fq.gz",
#'   depthFT = "True",
#'   mdepth = 10,
#'   GT_mut = "aa"
#' )
#' }

BSAPipeR <- function(ref, refbwa, outdir, sample_mut, sample_wt,GFF3,
                      mutfq1, mutfq2, wtfq1, wtfq2, depthFT = "False",
                      mdepth = 0,GT_mut="aa",mindp=4) {

  # 查看当前配置的路径
  get_software_paths()

  pddepth(depthFT, mdepth)

  # 为原始数据创建超链接到rawdata文件下
  rawdatadir <- file.path(outdir, "rawdata")
  if (!dir.exists(rawdatadir)) {
    dir.create(rawdatadir, recursive = TRUE)
  }
  cmd <- paste("ln -s", mutfq1, rawdatadir,sep = " ")
  system(cmd)
  cmd <- paste("ln -s", mutfq2, rawdatadir,sep = " ")
  system(cmd)
  cmd <- paste("ln -s", wtfq1, rawdatadir,sep = " ")
  system(cmd)
  cmd <- paste("ln -s", wtfq2, rawdatadir,sep = " ")
  system(cmd)

  mutfq1 <- file.path(rawdatadir, basename(mutfq1))
  mutfq2 <- file.path(rawdatadir, basename(mutfq2))
  wtfq1 <- file.path(rawdatadir, basename(wtfq1))
  wtfq2 <- file.path(rawdatadir, basename(wtfq2))

  # 运行分析流程
  run_vcf_processing(
    ref = ref,
    refbwa = refbwa,
    outdir = outdir,
    sample_mut = sample_mut,
    sample_wt = sample_wt,
    mutfq1 = mutfq1,
    mutfq2 = mutfq2,
    wtfq1 = wtfq1,
    wtfq2 = wtfq2,
    GFF3 = GFF3,
    depthFT = depthFT,
    mdepth = mdepth
  )
  # 创建结果
  create_results(
    outdir = outdir,
    sample_mut = sample_mut,
    sample_wt = sample_wt,
    depthFT = depthFT,
    mdepth = mdepth,
    GT_mut = GT_mut,
    mindp = mindp
  )

}


#' Title
#'
#' @param depthFT
#' @param mdepth
#'
#' @returns
#' @export
#'
#' @examples
pddepth <- function(depthFT, mdepth) {
  if (depthFT == "True") {
    if (mdepth <= 0) {
      stop("当 depthFT = 'True' 时，mdepth 必须大于 0")
    }
    return(mdepth)
  } else {
    return(0)
  }
}

