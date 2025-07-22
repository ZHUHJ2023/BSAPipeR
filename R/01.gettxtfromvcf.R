#' Title
#'
#' @param sample_wt 野生型的样本名称
#' @param sample_mut 突变体的样本名称
#' @param vcffile vcf文件的路径
#'
#' @importFrom dplyr filter select %>%
#' @importFrom vcfR read.vcfR extract.gt
#'
#' @returns null
#' @export
#'
#' @examples get_txt("B73QL", "B3304", "B3304.snp.VQSR.snpEff.vcf.gz")
get_txt <- function(sample_wt, sample_mut, vcffile) {
  # library(dplyr)
  # vcffile <- paste0("../largedata/",sample_mut,".snp.VQSR.snpEff.vcf.gz")
  vcf <- vcfR::read.vcfR(vcffile)

  df <- as.data.frame(vcf@fix)
  df$ID <- paste(df$CHROM, df$POS, sep = "_")
  AD <- vcfR::extract.gt(vcf, "AD")
  GT <- vcfR::extract.gt(vcf, "GT")

  # > colnames(df)
  # [1] "CHROM"  "POS"    "ID"     "REF"    "ALT"    "QUAL"   "FILTER" "INFO"

  df$mutation_effect <- lapply(strsplit(df$INFO, "|", fixed = TRUE), function(x) {
    x[2]
  }) %>% unlist()

  df$CDS_change <- lapply(strsplit(df$INFO, "|", fixed = TRUE), function(x) {
    x[10]
  }) %>% unlist()

  df$protein_change <- lapply(strsplit(df$INFO, "|", fixed = TRUE), function(x) {
    x[11]
  }) %>% unlist()

  df$Gene <- lapply(strsplit(df$INFO, "|", fixed = TRUE), function(x) {
    x[4]
  }) %>% unlist()

  # > colnames(df)
  #  [1] "CHROM"           "POS"             "ID"              "REF"
  #  [5] "ALT"             "QUAL"            "FILTER"          "INFO"
  #  [9] "mutation_effect" "CDS_change"      "protein_change"  "Gene"

  newdf <- df %>% dplyr::select(CHROM, POS, ID, REF, ALT, mutation_effect, CDS_change, protein_change, Gene)


  if (length(colnames(GT)) == 2) {
    newdf <- cbind(newdf, GT)
    colnames(newdf)[c(10, 11)] <- paste0("GT_", colnames(newdf)[c(10, 11)])
    newdf <- cbind(newdf, AD)
    colnames(newdf)[c(12, 13)] <- paste0("AD_", colnames(newdf)[c(12, 13)])
    # 仅保留有一个逗号的行，去掉有2个逗号的行；匹配"[1-9],[1-9]"
    index <- which(grepl(",", newdf$ALT))
    if (length(index) > 0) {
      AD <- AD[-index, ]
      newdf <- newdf[-index, ]
      GT <- GT[-index, ]
    }

    x <- AD[, 1]
    DP1 <- data.table::fread(text = x, header = F, sep = ",", data.table = F)
    colnames(DP1) <- c("ref_Depth", "alt_Depth")
    x <- AD[, 2]
    DP2 <- data.table::fread(text = x, header = F, sep = ",", data.table = F)
    colnames(DP2) <- c("ref_Depth", "alt_Depth")

    SNP_index <- DP1$alt_Depth / (DP1$ref_Depth + DP1$alt_Depth)
    newdf <- cbind(newdf, SNP_index)
    colnames(newdf)[c(14)] <- paste0("SNP_index_", colnames(AD)[1])
    SNP_index <- DP2$alt_Depth / (DP2$ref_Depth + DP2$alt_Depth)
    newdf <- cbind(newdf, SNP_index)
    colnames(newdf)[c(15)] <- paste0("SNP_index_", colnames(AD)[2])
    colnames(newdf) <- gsub(sample_wt, "WT", colnames(newdf))
    colnames(newdf) <- gsub(sample_mut, "MUT", colnames(newdf))
  }

  write.table(newdf, file = paste0(vcffile, ".txt"), sep = "\t", row.names = F, quote = F)
}
