# 从vcfgztxt文件获取得到xlsx
#' Title get_table
#' @description 从vcfgztxt文件获取得到xlsx
#' @param file get_txt得到的txt文件路径
#' @param chrlist 需要提取的染色体列表
#' @param mindp 最低测序深度，默认4
#' @param minSNPindexmut 最低SNPindex-M，默认0.3
#' @param maxSNPindexmut 最高SNPindex-M，默认0.7
#'
#' @returns Null
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom data.table fread
#' @importFrom utils write.table
#'
#' @examples get_table("B3304.snp.VQSR.snpEff.vcf.gz.txt", paste("chr", c(1:10), sep = ""), mindp = 4)
get_table <- function(file, chrlist, mindp = 4, minSNPindexmut = 0.3, maxSNPindexmut = 0.7,outdir,M) {
  df <- data.table::fread(file)
  df <- dplyr::filter(df, df$CHROM %in% chrlist)
  # 1.过滤掉WT杂合的突变
  fse <- dplyr::filter(df, !GT_WT %in% c("1/0", "1|0", "0/1", "0|1"))
  # 2.过滤掉WT与M相同的基因型
  fse <- dplyr::filter(fse, gsub("/", "|", GT_MUT) != gsub("/", "|", GT_WT))
  # #2.过滤掉WT中的aa
  # fse <- filter(fse,!GT_WT %in% c("1/1" ,"1|1"))
  # 3.过滤掉M中的aa----因为这次的突变体全部aa致死
  fse <- dplyr::filter(fse, !GT_MUT %in% c("0/0", "0|0"))
  fse <- dplyr::filter(fse, !GT_MUT %in% c("1/1", "1|1"))
  # 4.过滤测序深度低于4的SNP
  DP1 <- data.table::fread(text = fse$AD_WT, header = F, sep = ",", data.table = F)
  DP2 <- data.table::fread(text = fse$AD_MUT, header = F, sep = ",", data.table = F)
  index <- DP1[, 1] + DP1[, 2] >= mindp & DP2[, 1] + DP2[, 2] >= mindp
  index <- which(index)
  fse1 <- fse[index, ]
  # data.table::fwrite(fse,gsub("vcf.gz.txt","filter.mutation_effect.xls",file),sep = "\t",quote = F,row.names = F,col.names = T)
  # print filter SNP in mutation_effect
  fse <- dplyr::filter(fse1, grepl("splice_acceptor_variant|splice_donor_variant|splice_region_variant|stop_lost|start_lost|stop_gained|missense_variant|coding_sequence_variant|inframe_insertion|disruptive_inframe_insertion|inframe_deletion|disruptive_inframe_deletion|exon_variant|exon_loss_variant|exon_loss_variant|duplication|inversion|frameshift_variant|feature_ablation|duplication|gene_fusion|bidirectional_gene_fusion|rearranged_at_DNA_level|miRNA|initiator_codon_variant|start_retained", fse1$mutation_effect))
  # data.table::fwrite(fse,gsub("vcf.gz.txt","final.mutation_effect.xls",file),sep = "\t",quote = F,row.names = F,col.names = T)
  # diuqi <- data.table::fread("../data/Maa_WTAA.snpindel.txt")
  # diuqi$ID <- paste(diuqi$CHROM,diuqi$POS,sep = "_")
  f <- dplyr::filter(fse, SNP_index_MUT > minSNPindexmut & SNP_index_MUT < maxSNPindexmut)
  f$ID <- paste(f$CHROM, f$POS, sep = "_")
  # ff <- filter(f,!f$ID %in% diuqi$ID)
  # library(openxlsx)
  ## Create a new workbook
  wb <- openxlsx::createWorkbook()
  ## Add 3 worksheets
  openxlsx::addWorksheet(wb, sheetName = "00_raw", tabColour = "white") # 自定义sheet名称
  openxlsx::addWorksheet(wb, sheetName = paste0("01_MAa_WTAA_DP>", mindp), tabColour = "grey") # 自定义sheet名称
  openxlsx::addWorksheet(wb, sheetName = "02_mutation_effect", tabColour = "grey") ## 给sheet标签添加颜色
  openxlsx::addWorksheet(wb, sheetName = paste0("03_", minSNPindexmut, "<SNPindexM<", maxSNPindexmut), tabColour = "#A593E0")

  ## Write data
  openxlsx::writeData(wb, sheet = 1, df)
  openxlsx::writeData(wb, sheet = 2, fse1)
  openxlsx::writeData(wb, sheet = 3, fse)
  openxlsx::writeData(wb, sheet = 4, f)

  ## 保存工作簿
  openxlsx::saveWorkbook(wb, paste0(outdir,"/Aa.",M,".CandidateGene.xlsx"), overwrite = TRUE)
}
