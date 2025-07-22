#' Title plot
#'
#' @description plot
#' @param file get_table得到的xlsx文件路径
#'
#' @importFrom dplyr filter
#' @importFrom openxlsx
#' @importFrom grDevices pdf png dev.off
#' @importFrom graphics par
#'
#' @returns null
#' @export
#'
#' @examples get_graphs("B3304.snp.VQSR.snpEff.xlsx")
get_graphs <- function(file,outdir,sample_mut) {
  # library(ggplot2)
  # library(CMplot)
  # library(dplyr)
  # source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")
  df <- openxlsx::read.xlsx(file, sheet = "01_MAa_WTAA_DP>4")
  df <- as.data.frame(df)
  # df <- filter(df,df$SNP_index_MUT < 0.8 & df$SNP_index_MUT > 0.2)
  chrom <- df[, 1]
  pos <- df[, 2]
  window <- 2e6
  # colnames(df)[1] <- 'CHROM'
  # colnames(df)[2] <- 'POS'
  res <- sapply(1:nrow(df), function(i) {
    if (i %% 1000 == 0) {
      print(i)
    }
    ii <- which(chrom == chrom[i] & abs(pos - pos[i]) < window)
    # N0.7 = sum(df$SNP_index_MUT[ii] >= 0.7)
    # N0.8 = sum(df$SNP_index_MUT[ii] >= 0.8)
    # N0.85 = sum(df$SNP_index_MUT[ii] >= 0.85)
    # N0.8 = sum(df$SNP_index_MUT[ii] >= 0.8)
    # N0.9 = sum(df$SNP_index_MUT[ii] >= 0.9)
    # N0.95 = sum(df$SNP_index_MUT[ii] >= 0.95)
    N0.8 <- sum(abs(df$SNP_index_MUT[ii] - 0.5) <= 0.2)
    N0.9 <- sum(abs(df$SNP_index_MUT[ii] - 0.5) <= 0.1)
    N0.95 <- sum(abs(df$SNP_index_MUT[ii] - 0.5) <= 0.05)

    c(length(ii), mean(df$SNP_index_MUT[ii]), mean(df$SNP_index_WT[ii]), N0.8, N0.9, N0.95)
    # c(length(ii),mean(df$SNP_index_MUT[ii]),mean(df$SNP_index_WT[ii]),N0.7,N0.8,N0.85)
  })
  res <- t(res)
  colnames(res) <- c("N", "Mean_SNP-index", "Mean_SNP-index_WT", "N0.8", "N0.9", "N0.95")
  # colnames(res) = c("N","Mean_SNP-index","Mean_SNP-index_WT","N0.7","N0.8","N0.85")

  d <- cbind(df, res)
  d$INFO <- NULL
  d <- dplyr::filter(d, d$CHROM %in% paste("chr", c(1:10), sep = ""))
  # dx = d[,c("ID","CHROM","POS","SNP_index_MUT","Mean_SNP-index","N0.7","N0.8","N0.85")]
  dx <- d[, c("ID", "CHROM", "POS", "SNP_index_MUT", "Mean_SNP-index", "N0.8", "N0.9", "N0.95")]

  # pdf(file=gsub(".xlsx",".CMplot.pdf",file),width = 10,height = 3)
  outfile=paste0(outdir,"/Aa.",gsub(".VQSR.snpEff.xlsx",".CMplot.png",basename(file)))
  png(filename = outfile, width = 1000, height = 1000)
  par(ps = 10, mar = c(3, 4, 1, 0.5), mgp = c(1, 0.5, 0), mfcol = c(5, 1), cex = 2)
  CMplot(dx,
    LOG10 = F, plot.type = "m", chr.labels = paste("Chr", c(1:10), sep = ""), r = 0.4, cir.axis = TRUE,
    outward = FALSE, cir.axis.col = "black", cir.chr.h = 1.3, chr.den.col = "black", file = "jpg",
    file.name = "", dpi = 300, file.output = F, verbose = TRUE, width = 10, height = 3
  )
  dev.off()
  outfile=paste0(outdir,"/Aa.",gsub(".VQSR.snpEff.xlsx",".CMplot.pdf",basename(file)))
  pdf(file = outfile, width = 10, height = 10)
  par(ps = 10, mar = c(3, 4, 1, 0.5), mgp = c(1, 0.5, 0), mfcol = c(5, 1), cex = 2)
  CMplot(dx,
    LOG10 = F, plot.type = "m", chr.labels = paste("Chr", c(1:10), sep = ""), r = 0.4, cir.axis = TRUE,
    outward = FALSE, cir.axis.col = "black", cir.chr.h = 1.3, chr.den.col = "black", file = "jpg",
    file.name = "", dpi = 300, file.output = F, verbose = TRUE, width = 10, height = 3
  )
  dev.off()
}
