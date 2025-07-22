#' Title
#'
#' @param outdir
#' @param sample_mut
#' @param sample_wt
#' @param depthFT
#' @param mdepth
#' @param GT_mut
#'
#' @returns
#' @export
#'
#' @examples
create_results <- function(outdir, sample_mut, sample_wt, depthFT = "False", mdepth = 0,GT_mut = "aa", mindp = 4) {
  workdir <- paste0(outdir, "/Processing")
  # snp_vcf="$workdir/05_snpIndel/${depth_prefix}${sample_mut}.snp.VQSR.vcf.gz"
  # indel_vcf="$workdir/05_snpIndel/${depth_prefix}${sample_mut}.indel.VQSR.vcf.gz"
  # local depth_prefix=""
  # [[ "$mdepth" -gt 0 ]] && depth_prefix="${mdepth}x."
  if (depthFT == "True") {
    depth_prefix <- paste0(mdepth, "x.")
  } else {
    depth_prefix <- ""
  }
  snp_vcf <- paste0(workdir, "/06_SnpEff/",sample_mut,"/" ,depth_prefix, sample_mut, ".snp.VQSR.snpEff.vcf.gz")
  indel_vcf <- paste0(workdir, "/06_SnpEff/",sample_mut,"/" ,depth_prefix, sample_mut, ".indel.VQSR.snpEff.vcf.gz")
  snp_html <- paste0(workdir, "/06_SnpEff/",sample_mut,"/" ,depth_prefix, sample_mut, ".snp.VQSR.snpEff.html")
  indel_html <- paste0(workdir, "/06_SnpEff/",sample_mut,"/" ,depth_prefix, sample_mut, ".indel.VQSR.snpEff.html")
  # Create results directory if it doesn't exist
  results_dir <- paste0(outdir, "/Results/01_raw_vcf")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  # Cp results files
  system(paste0("cp ", snp_vcf, " ", results_dir, "/"))
  system(paste0("cp ", indel_vcf, " ", results_dir, "/"))
  system(paste0("cp ", snp_html, " ", results_dir, "/"))
  system(paste0("cp ", indel_html, " ", results_dir, "/"))
  # Create a table from the html files
  snp_vcf <- paste0(results_dir, "/", basename(snp_vcf))
  indel_vcf <- paste0(results_dir, "/", basename(indel_vcf))
  snp_html <- paste0(results_dir, "/", basename(snp_html))
  indel_html <- paste0(results_dir, "/", basename(indel_html))
  get_table_from_html(snp_html)
  get_table_from_html(indel_html)
  # CK results (我们的BSA流程的结果)
  results_ck_dir <- paste0(outdir, "/Results/02_BSA_results")
  if (!dir.exists(results_ck_dir)) {
    dir.create(results_ck_dir, recursive = TRUE)
  }
  get_results_ck(vcffile = snp_vcf,outdir = results_ck_dir,sample_wt = sample_wt,
                 sample_mut = sample_mut,GT_mut = GT_mut, mindp = mindp)
  # ML model
  results_ml_dir <- paste0(outdir, "/Results/03_ML_results")
  if (!dir.exists(results_ml_dir)) {
    dir.create(results_ml_dir, recursive = TRUE)
  }
  get_results_ml(vcffile = snp_vcf, outdir = results_ml_dir, sample_wt = sample_wt,
                 sample_mut = sample_mut, GT_mut = "aa",mindp=mindp)
  # 输出结果说明
  print("Results have been generated in the following directories:")
  print(paste0("1. Raw VCF files: ", results_dir))
  print(paste0("2. BSA results: ", results_ck_dir))
  print(paste0("3. Model prediction result: ", results_ml_dir))
}


#' Extract Tables from HTML File and Save to Excel
#'
#' This function parses an HTML file containing tables, extracts specific tables of interest,
#' and saves them to an Excel workbook with multiple worksheets.
#'
#' @param htmlfile Path to the input HTML file containing tables to extract.
#'
#' @importFrom XML htmlParse getNodeSet readHTMLTable
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom dplyr filter %>%
#'
#' @return Invisibly returns NULL. The main output is an Excel file saved to disk.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' get_table_from_html("AeOM/AeOM.snp.VQSR.20x.snpEff.html")
#' }
get_table_from_html <- function(htmlfile) {
  # Parse the HTML file
  doc <- XML::htmlParse(htmlfile)

  # Get all tables from the HTML document
  total_table <- XML::getNodeSet(doc, "//table")

  # Read all tables into a list
  tablelist <- lapply(total_table, XML::readHTMLTable)

  # Extract specific tables of interest
  df3 <- tablelist[[3]]  # Variants rate details
  df4 <- tablelist[[4]]  # Number variants by type
  df5 <- tablelist[[5]]  # Number of effects by impact
  df6 <- tablelist[[6]]  # Number of effects by functional class
  df7 <- tablelist[[7]]  # Number of effects by type and region

  # Process table 7 to separate type and region data
  df <- data.frame(
    "Type" = df7$V1[-1:-3],
    "Count" = df7$V3[-1:-3],
    "Percent" = df7$V4[-1:-3]
  )

  # Safely find the separator row
  n <- which(df$Type %in% "Type (alphabetical order)")
  # print(n)  # Debugging line to check separator index
  # Handle case where separator isn't found
  if (length(n) == 0) {
    df7_1 <- df  # use all data as type data
    df7_2 <- data.frame()  # empty region data
  } else {
    n <- as.numeric(n)
    df7_1 <- df[1:(n-1), ]  # type data
    df7_2 <- df[(n+1):nrow(df), ]  # region data
  }

  # Create Excel workbook and add worksheets
  wb <- openxlsx::createWorkbook()

  # Add sheets with descriptive names
  openxlsx::addWorksheet(wb, sheetName = "01_Variants_rate_details")
  openxlsx::writeData(wb, sheet = "01_Variants_rate_details", df3)

  openxlsx::addWorksheet(wb, sheetName = "02_Number_variants_by_type")
  openxlsx::writeData(wb, sheet = "02_Number_variants_by_type", df4)

  openxlsx::addWorksheet(wb, sheetName = "03_Number_of_effects_by_impact")
  openxlsx::writeData(wb, sheet = "03_Number_of_effects_by_impact", df5)

  openxlsx::addWorksheet(wb, sheetName = "04_N_by_functional_class")
  openxlsx::writeData(wb, sheet = "04_N_by_functional_class", df6)

  openxlsx::addWorksheet(wb, sheetName = "05_Number_of_effects_by_type")
  openxlsx::writeData(wb, sheet = "05_Number_of_effects_by_type", df7_1)

  openxlsx::addWorksheet(wb, sheetName = "06_Number_of_effects_by_region")
  openxlsx::writeData(wb, sheet = "06_Number_of_effects_by_region", df7_2)

  # Save the Excel file with same name as HTML but .xlsx extension
  output_file <- gsub(".html", ".xlsx", htmlfile)
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

  invisible(NULL)
}

#' Process VCF File for SNP Analysis
#'
#' This function analyzes VCF files to identify candidate genes based on SNP indices,
#' generates summary tables, and creates visualization plots. It handles both homozygous
#' (aa) and heterozygous (Aa) mutation cases differently.
#'
#' @param vcffile Path to the input VCF file
#' @param outdir Output directory for results
#' @param sample_wt Name of the wild-type sample in VCF
#' @param sample_mut Name of the mutant sample in VCF
#' @param GT_mut Genotype of mutant, either "aa" (homozygous) or "Aa" (heterozygous)
#'
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom dplyr filter mutate %>%
#' @importFrom data.table fread
#' @importFrom furrr future_map_dfr
#' @importFrom future plan multisession availableCores
#' @importFrom ggplot2 ggplot geom_smooth scale_x_continuous ylab xlab theme_bw theme
#' @importFrom ggplot2 element_blank element_text element_line unit margin
#' @importFrom purrr map2_dfr
#' @importFrom utils write.table
#' @importFrom stats na.omit
#'
#' @return Invisibly returns NULL. Outputs are saved to files in the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' get_results_ck("input.vcf", "output_dir", "WT_sample", "MUT_sample", "aa")
#' }
get_results_ck <- function(vcffile, outdir, sample_wt, sample_mut, GT_mut, mindp = 4) {

  # Internal function to process chromosome data
  process_chr <- function(chr_data, chr, window_size, step_size) {
    chr_min <- min(chr_data$POS, na.rm = TRUE)
    chr_max <- max(chr_data$POS, na.rm = TRUE)

    if (is.finite(chr_min) && is.finite(chr_max) && chr_max > chr_min) {
      starts <- seq(chr_min, chr_max - window_size, by = step_size)
      ends <- starts + window_size

      purrr::map2_dfr(starts, ends, function(start, end) {
        window_data <- chr_data %>% dplyr::filter(POS >= start & POS <= end)
        if (nrow(window_data) > 0) {
          boot_result <- c(mean(window_data$SNP_index, na.rm = TRUE),
                           mean(window_data$SNP_index_WT, na.rm = TRUE))
          data.frame(CHROM = chr, Start = start, End = end,
                     Mean_SNP_index = boot_result[1],
                     Mean_SNP_index_WT = boot_result[2])
        } else {
          NULL
        }
      })
    } else {
      return(NULL)
    }
  }

  # Internal function to process VCF and generate plots (for aa genotype)
  plot_vcf <- function(vcffile, M, WT, outdir) {
    # Read and process VCF
    suppressWarnings({
      vcf <- vcfR::read.vcfR(vcffile)
    })
    df <- as.data.frame(vcf@fix)
    gt <- as.data.frame(vcf@gt)
    df$ID <- paste(df$CHROM, df$POS, sep = "_")

    GT <- vcfR::extract.gt(vcf, "GT")
    AD <- vcfR::extract.gt(vcf, "AD")

    # Filter for mutant genotype (not 0/0 or 0|0)
    index <- !(GT[, M] %in% c("0/0", "0|0"))
    df <- df[index, ]
    GT <- GT[index, ]
    AD <- AD[index, ]

    # Calculate SNP indices
    x <- AD[, M]
    DP1 <- data.table::fread(text = x, header = FALSE, sep = ",", data.table = FALSE)
    colnames(DP1) <- c("mut_ref_Depth", "mut_alt_Depth")
    df$SNP_index <- DP1[, 2] / (DP1[, 1] + DP1[, 2])

    x <- AD[, WT]
    DP2 <- data.table::fread(text = x, header = FALSE, sep = ",", data.table = FALSE)
    colnames(DP2) <- c("wt_ref_Depth", "wt_alt_Depth")
    df$SNP_index_WT <- DP2[, 2] / (DP2[, 1] + DP2[, 2])

    # Filter for minimum depth of 4
    x1 <- AD[, M]
    DP1 <- data.table::fread(text = x1, header = FALSE, sep = ",", data.table = FALSE)
    x2 <- AD[, WT]
    DP2 <- data.table::fread(text = x2, header = FALSE, sep = ",", data.table = FALSE)
    index <- (DP1[, 1] + DP1[, 2] >= mindp) & (DP2[, 1] + DP2[, 2] >= mindp)
    df <- df[index, ]
    GT <- GT[index, ]
    AD <- AD[index, ]

    # Calculate ED (Euclidean Distance)
    ad_flt <- AD[, c(WT, M)]
    df$ED <- apply(ad_flt, 1, function(x) {
      count <- as.numeric(unlist(strsplit(x, ",", fixed = TRUE, useBytes = TRUE)))
      depth1 <- count[1] + count[2]
      depth2 <- count[3] + count[4]

      ED <- sqrt((count[3] / depth2 - count[1] / depth1)^2 +
                   (count[4] / depth2 - count[2] / depth1)^2)
      return(ED^4)
    })

    # Filter for candidate genes
    res <- df %>%
      dplyr::filter(SNP_index >= 0.9 & SNP_index_WT <= 0.1) %>%
      dplyr::filter(grepl(paste0("splice_acceptor_variant|splice_donor_variant|",
                                 "splice_region_variant|stop_lost|start_lost|",
                                 "stop_gained|missense_variant|coding_sequence_variant|",
                                 "inframe_insertion|disruptive_inframe_insertion|",
                                 "inframe_deletion|disruptive_inframe_deletion|",
                                 "exon_variant|exon_loss_variant|duplication|",
                                 "inversion|frameshift_variant|feature_ablation|",
                                 "gene_fusion|bidirectional_gene_fusion|",
                                 "rearranged_at_DNA_level|miRNA|initiator_codon_variant|",
                                 "start_retained"), INFO))
    outfile <- paste0(outdir, "/", gsub(".VQSR.snpEff.vcf.gz", ".CandidateGene.txt", basename(vcffile)))
    utils::write.table(res, outfile,
                       row.names = FALSE, quote = FALSE, sep = "\t")

    # Filter for different genotypes between WT and mutant
    index <- !(gsub('/', '|', GT[, WT]) == gsub('/', '|', GT[, M]))
    df <- df[index, ]
    GT <- GT[index, ]
    AD <- AD[index, ]

    # Process chromosome data
    # 筛选主要染色体
    df$CHROM <- gsub("chr", "", df$CHROM)
    # 去掉chr，pos为na的行
    df <- df %>%
      dplyr::filter(!is.na(POS) & !is.na(CHROM))
    # 定义主要染色体模式
    maize_pattern <- "^\\d{1,2}$"  # 匹配纯数字染色体 (1-10)
    wheat_pattern <- "^[1-7][ABD]$|^Un$"  # 匹配小麦染色体 (1A,1B,1D,...,7A,7B,7D) 和 Un

    # 检测染色体模式
    is_maize <- all(grepl(maize_pattern, unique(df$CHROM)[1:min(10, length(unique(df$CHROM)))]))
    is_wheat <- all(grepl(wheat_pattern, unique(df$CHROM)[1:min(21, length(unique(df$CHROM)))]))

    if (is_maize) {
      # 玉米: 选择 1-10 号染色体
      main_chroms <- as.character(1:10)
      df <- df %>% dplyr::filter(CHROM %in% main_chroms)
      chrom_levels <- main_chroms
    } else if (is_wheat) {
      # 小麦: 选择 1A-7D 和 Un 染色体
      wheat_chroms <- c(paste0(rep(1:7, each = 3), c('A','B','D')))
      df <- df %>% dplyr::filter(CHROM %in% wheat_chroms)
      chrom_levels <- wheat_chroms
    } else {
      # 其他物种: 使用所有染色体
      chrom_levels <- sort(unique(df$CHROM))
      warning("无法自动识别物种，使用所有染色体")
    }

    # Calculate window means in parallel
    df$POS <- as.numeric(as.character(df$POS))
    window_size <- 2e6
    step_size <- 10e3
    # 查询可用的工作线程数
    t <- future::availableCores()
    future::plan(future::multisession, workers = t)
    results <- furrr::future_map_dfr(unique(df$CHROM), function(chr) {
      chr_data <- df %>% dplyr::filter(CHROM == chr)
      process_chr(chr_data, chr, window_size, step_size)
    })
    future::plan(future::sequential)

    colnames(results) <- c("CHROM", "Start", "End", "Mean_SNP-index", "Mean_SNP-index_WT")

    # Generate plot
    d <- df
    d$POS <- as.numeric(as.character(d$POS))
    results$Start <- as.numeric(as.character(results$Start))
    results$End <- as.numeric(as.character(results$End))
    results$`Mean_SNP-index` <- as.numeric(results$`Mean_SNP-index`)
    d$CHROM <- factor(d$CHROM, levels = chrom_levels)
    results$CHROM <- factor(results$CHROM, levels = chrom_levels)

    p <- ggplot2::ggplot() +
      ggplot2::facet_grid(~CHROM, scales = "free_x") +
      ggplot2::geom_smooth(data = results,
                           aes(x = (Start + End)/2/1e6,
                               y = round(`Mean_SNP-index`, 4)),
                           size = 0.8, color = "#74C1DD",
                           method = "loess", span = 0.3) +
      ggplot2::scale_x_continuous(
        breaks = seq(0, max(d$POS/1e6), 100),
        labels = ceiling(seq(0, max(d$POS/1e6), 100)/100)) +
      ggplot2::ylab(paste0(M, " SNP index")) +
      ggplot2::xlab("Chromosome position (x100Mb)") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.spacing = unit(0, "cm", data = NULL),
        strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        strip.text.x = element_text(size = 8, angle = 0)
      )
    outfile=paste0(outdir, "/", gsub(".VQSR.snpEff.vcf.gz",".SNP_index.pdf",basename(vcffile)))
    ggplot2::ggsave(outfile,p, width = 8, height = 3.5)
  }

  # Main function logic
  if (GT_mut == "aa") {
    plot_vcf(vcffile = vcffile,
             M = sample_mut,
             WT = sample_wt,
             outdir = outdir)
  } else if (GT_mut == "Aa") {
    get_txt(sample_wt, sample_mut, vcffile)
    get_table(file = paste0(vcffile, ".txt"),
              chrlist = paste("chr", c(1:10), sep = ""),
              mindp = mindp,
              minSNPindexmut = 0.3,
              maxSNPindexmut = 0.7,
              M = sample_mut,
              outdir = outdir)
    get_graphs(file = paste0(outdir, "/Aa.", sample_mut, ".CandidateGene.xlsx"),
               sample_mut = sample_mut,
               outdir = outdir)
  }

  invisible(NULL)
}


#' Process VCF File Using Machine Learning Model
#'
#' This function analyzes VCF files using a pre-trained machine learning model to identify
#' candidate mutations. It handles data processing, feature engineering, model prediction,
#' and result visualization.
#'
#' @param vcffile Path to the input VCF file (required)
#' @param outdir Output directory for results (required)
#' @param sample_wt Name of the wild-type sample in VCF (required)
#' @param sample_mut Name of the mutant sample in VCF (required)
#' @param GT_mut Genotype of mutant ("aa" only, currently supported)
#' @param model_path Optional path to custom model file (default: uses built-in model)
#'
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom dplyr select filter bind_cols %>% mutate if_else
#' @importFrom data.table fread
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook read.xlsx
#' @importFrom ggplot2 ggplot aes geom_point facet_grid scale_color_gradientn labs
#' @importFrom ggplot2 theme_bw theme element_blank element_line element_text unit margin
#' @importFrom workflowsets extract_workflow
#' @importFrom stats predict na.omit
#' @importFrom utils str
#'
#' @return Invisibly returns NULL. Outputs are saved to files in the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' get_results_ml("input.vcf", "output_dir", "WT_sample", "MUT_sample", "aa")
#' }
get_results_ml <- function(vcffile, outdir, sample_wt, sample_mut, GT_mut,mindp=4, model_path = NULL) {
    # Parameter validation
    if (missing(vcffile) || !file.exists(vcffile)) {
      stop("VCF file is required and must exist")
    }

    if (missing(outdir)) {
      stop("Output directory is required")
    }

    if (missing(sample_wt) || missing(sample_mut)) {
      stop("Sample names for wild-type and mutant are required")
    }

    if (GT_mut != "aa") {
      stop("ML model currently only supports 'aa' mutation type")
    }

    # Create output directory if it doesn't exist
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    # Check for required packages
    required_pkgs <- c("vcfR", "data.table", "dplyr", "tidymodels", "openxlsx",
                       "caret", "class", "ggplot2", "tidymodels")
    missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))

    if (length(missing_pkgs) > 0) {
      stop(paste("The following required packages are missing:",
                 paste(missing_pkgs, collapse = ", ")))
    }

    # Read and process VCF file
    df <- get_vcf_features(vcffile = vcffile,
                           sample_wt = sample_wt,
                           sample_mut = sample_mut) %>%
      mutate(Sample = sample_mut) %>%
      mutate(
        DP_WT_norm = DP_WT / mean(DP_WT, na.rm = TRUE),
        DP_MUT_norm = DP_MUT / mean(DP_MUT, na.rm = TRUE) ) %>%
      filter(DP_WT >= mindp & DP_MUT >= mindp) #修改为使用测序深度大于等于4的位点
    df$ID <- paste0(df$CHROM, "_", df$POS, "_", df$Sample)
    df <- df %>%
      mutate(
        DP_WT_norm = DP_WT / mean(DP_WT, na.rm = TRUE),
        DP_MUT_norm = DP_MUT / mean(DP_MUT, na.rm = TRUE)
      )
    df$detaAF <- df$AF_MUT - df$AF_WT
    df$FILTER <- ifelse(df$FILTER == "PASS", 1, 0)
    df$GT_MUT <- gsub("/","|",df$GT_MUT)
    df$GT_WT <- gsub("/","|",df$GT_WT)
    gxdata <- df
    if ( "mutation_effect_level" %in% colnames(gxdata) ) {
      gxdata$mutation_effect_level <- factor(gxdata$mutation_effect_level,levels=c("HIGH","MODERATE","MODIFIER","LOW"))
    }
    gxdata$GT_MUT <- factor(gxdata$GT_MUT, levels = c("1|1", "0|1", "0|0"))
    gxdata$GT_WT <- factor(gxdata$GT_WT, levels = c("0|0", "0|1", "1|1"))
    gxdata$Pr_change_is <- ifelse(gxdata$Pr_change_is %in% "1", "TRUE", "FALSE")
    gxdata$Pr_change_is <- factor(gxdata$Pr_change_is,levels=c("TRUE","FALSE"))
    gxdata$EMS_is <- ifelse(gxdata$EMS_is %in% "1", "TRUE", "FALSE")
    gxdata$EMS_is <- factor(gxdata$EMS_is,levels=c("TRUE","FALSE"))
    # Load model
    if (is.null(model_path)) {
      model_file <- system.file("extdata", "final.lightgbm_fit.rds",
                                package = "BSAPipeR")
      if (!file.exists(model_file)) {
        stop("Built-in model file not found in package")
      }
    } else if (!file.exists(model_path)) {
      stop("Specified model file does not exist")
    } else {
      model_file <- model_path
    }
    lightgbm_fit <- readRDS(model_file) %>% workflowsets::extract_workflow()
    pred_lightgbm <- gxdata %>%
      bind_cols(predict(lightgbm_fit, gxdata, type = "prob")) %>%
      bind_cols(predict(lightgbm_fit, gxdata, type = "class"))
    f_lightgbm <- dplyr::filter(pred_lightgbm,.pred_class==TRUE)
    #06.保存结果

    # library(openxlsx)
    ## Create a new workbook
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "00_raw_dataframe", tabColour = "white")
    openxlsx::addWorksheet(wb, sheetName = "01_Pred_LightGBM", tabColour = "grey")
    openxlsx::addWorksheet(wb, sheetName = "02_Pred_LightGBM_TRUE", tabColour = "#A593E0")

    openxlsx::writeData(wb, sheet = 1, gxdata)
    openxlsx::writeData(wb, sheet = 2, pred_lightgbm)
    openxlsx::writeData(wb, sheet = 3, f_lightgbm)
    outfile <- paste0(outdir, "/", gsub(".VQSR.snpEff.vcf.gz",".Predict_lightgbm.xlsx",basename(vcffile)))
    openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)

    # Generate visualization
    Pred_Final <- openxlsx::read.xlsx(outfile, sheet = "01_Pred_LightGBM")
    Pred_Final$CHROM <-  gsub("chr","",Pred_Final$CHROM)
    Pred_Final$POS <- as.numeric(Pred_Final$POS)
    Pred_Final$CHROM <-  factor(Pred_Final$CHROM,levels = seq(1:12))
    Pred_Final$level <- ifelse(Pred_Final$mutation_effect_level == "HIGH",4,ifelse(Pred_Final$mutation_effect_level == "MODERATE",3,ifelse(Pred_Final$mutation_effect_level == "MODIFIER",2,1)))

    p4 <- ggplot2::ggplot(Pred_Final, ggplot2::aes(x = .data$POS/1e6, y = .data$.pred_TRUE,
                                                   color = .data$level)) +
      ggplot2::geom_point(size = 0.5) +
      ggplot2::facet_grid(~.data$CHROM, scales = "free_x", space = "free_x", switch = "x") +
      ggplot2::scale_color_gradientn(colours = c("#f9d423", "#f83600", "#020f75")) +
      ggplot2::labs(title = "", x = "", y = "LightGBM") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 8, color = "black"),
        strip.background = ggplot2::element_blank(),
        strip.placement = "outside",
        axis.line.x = ggplot2::element_line(linetype = 1, color = "black", linewidth = 0.5),
        axis.line.y = ggplot2::element_line(linetype = 1, color = "black", linewidth = 0.5),
        panel.spacing = ggplot2::unit(0, "cm"),
        axis.title.y = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 8),
        axis.title = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 8),
        legend.position = "right",
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = ggplot2::unit(0.2, "cm"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
        plot.margin = ggplot2::margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
      )

    pdf_file <- gsub(".xlsx", ".pdf", outfile)
    grDevices::pdf(file = pdf_file, width = 8, height = 3)
    print(p4)
    grDevices::dev.off()

    invisible(NULL)
}

#' Extract and Annotate Variant Features from VCF File
#'
#' Processes a VCF file to extract genomic variants, compute allele frequencies,
#' add functional annotations, and calculate quality metrics for mutation analysis.
#' Designed for paired wild-type (WT) and mutant (MUT) samples.
#'
#' @param vcffile Path to the input VCF file (compressed or uncompressed).
#' @param sample_wt Name of the wild-type sample as it appears in the VCF header.
#' @param sample_mut Name of the mutant sample as it appears in the VCF header.
#'
#' @return A data frame containing annotated variant features with the following columns:
#' \itemize{
#'   \item \strong{Basic Variant Info:} CHROM, POS, ID (chromosome_position), REF, ALT, FILTER
#'   \item \strong{Functional Annotation:} mutation_effect, mutation_effect_level, CDS_change, protein_change, Gene
#'   \item \strong{Genotype Calls:} GT_WT, GT_MUT (genotype strings)
#'   \item \strong{Depth Metrics:} AD_WT, AD_MUT (allelic depths), DP_WT, DP_MUT (total depths)
#'   \item \strong{Quality Metrics:} PL_WT, PL_MUT (phred-scaled likelihoods)
#'   \item \strong{Frequency Metrics:} AF_WT, AF_MUT (allele frequencies), detaAF (AF_MUT - AF_WT)
#'   \item \strong{Variant Characteristics:} RTA (reference>alternate), EMS_is (1 if EMS-type mutation), Pr_change_is (1 if protein-changing)
#'   \item \strong{Statistical Metrics:} ED4 (Euclidean distance^4 of allele frequencies)
#'   \item \strong{Regional Density:} SNP_count_2Mb (variants in ±1Mb window)
#' }
#'
#' @details
#' Key processing steps include:
#' \enumerate{
#'   \item Parse VCF using \code{vcfR} and extract fixed fields + genotype matrices (GT/AD/PL)
#'   \item Split INFO field to extract functional annotations (VEP/Snpeff-style)
#'   \item Calculate allele frequencies and depths from AD fields
#'   \item Filter multi-allelic variants and missing/invalid genotypes
#'   \item Compute variant-specific features:
#'     \itemize{
#'       \item EMS mutation flag (G>A or C>T transitions)
#'       \item Protein change impact (excluding synonymous)
#'       \item Allelic divergence metric (ED4)
#'     }
#'   \item Annotate regional variant density using parallel computation
#' }
#'
#' @note
#' \itemize{
#'   \item Requires exactly 2 samples in VCF (WT and MUT)
#'   \item INFO field must contain pipe-delimited annotations (standard in VEP/Snpeff)
#'   \item Uses 8-core parallelization for regional density calculation (modify with \code{registerDoParallel})
#'   \item Multi-allelic sites are automatically excluded
#' }
#'
#' @section Warning:
#' Parallel computation uses \code{foreach} and \code{doParallel}. Ensure sufficient memory for large VCFs.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' vcf_path <- system.file("extdata", "sample.vcf.gz", package = "vcfR")
#' features <- get_vcf_features(
#'   vcffile = vcf_path,
#'   sample_wt = "P1_WT",
#'   sample_mut = "P1_MUT"
#' )
#' head(features)
#' }
#'
#' @import vcfR
#' @import dplyr
#' @import data.table
#' @import foreach
#' @import doParallel
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom data.table fread
#' @export

get_vcf_features <- function(vcffile, sample_wt, sample_mut) {
  vcf <- read.vcfR(vcffile)
  df <- as.data.frame(vcf@fix)
  df$ID <- paste(df$CHROM,df$POS,sep="_")
  AD <- extract.gt(vcf, "AD")
  GT <- extract.gt(vcf, "GT")
  PL <- extract.gt(vcf, "PL")
  # > colnames(df)
  # [1] "CHROM"  "POS"    "ID"     "REF"    "ALT"    "QUAL"   "FILTER" "INFO"

  df$mutation_effect <- lapply(strsplit(df$INFO, "|",fix=TRUE), function(x) {
    x[2]
  }) %>% unlist()
  df$mutation_effect_level <- lapply(strsplit(df$INFO, "|",fix=TRUE), function(x) {
    x[3]
  }) %>% unlist()

  df$CDS_change <- lapply(strsplit(df$INFO, "|",fix=TRUE), function(x) {
    x[10]
  }) %>% unlist()

  df$protein_change <- lapply(strsplit(df$INFO, "|",fix=TRUE), function(x) {
    x[11]
  }) %>% unlist()

  df$Gene <- lapply(strsplit(df$INFO, "|",fix=TRUE), function(x) {
    x[4]
  }) %>% unlist()

  # > colnames(df)
  #  [1] "CHROM"           "POS"             "ID"              "REF"
  #  [5] "ALT"             "QUAL"            "FILTER"          "INFO"
  #  [9] "mutation_effect" "CDS_change"      "protein_change"  "Gene"

  newdf <- df %>% select(CHROM,POS,ID,REF,ALT,FILTER,mutation_effect,mutation_effect_level,CDS_change,protein_change,Gene)

  if ( length(colnames(GT)) == 2 ) {
    n <- dim(newdf)[2]
    newdf <- cbind(newdf,GT)
    colnames(newdf)[c(n+1,n+2)] <- paste0("GT_",colnames(newdf)[c(n+1,n+2)])
    newdf <- cbind(newdf,AD)
    colnames(newdf)[c(n+3,n+4)] <- paste0("AD_",colnames(newdf)[c(n+3,n+4)])
    newdf <- cbind(newdf,PL)
    colnames(newdf)[c(n+5,n+6)] <- paste0("PL_",colnames(newdf)[c(n+5,n+6)])
    #仅保留有一个逗号的行，去掉有2个逗号的行；匹配"[1-9],[1-9]"
    index <- which(grepl(",", newdf$ALT))
    if (length(index) > 0) {
      AD <- AD[-index,]
      newdf <- newdf[-index,]
      GT <- GT[-index,]
    }

    x=AD[,1]
    DP1 = fread(text = x,header=F,sep=",",data.table=F)
    colnames(DP1) = c("ref_Depth","alt_Depth")
    x=AD[,2]
    DP2 = fread(text = x,header=F,sep=",",data.table=F)
    colnames(DP2) = c("ref_Depth","alt_Depth")

    SNP_index <- DP1$alt_Depth/(DP1$ref_Depth + DP1$alt_Depth )
    newdf <- cbind(newdf,SNP_index)
    colnames(newdf)[c(n+7)] <- paste0("AF_",colnames(AD)[1])
    SNP_index <- DP2$alt_Depth/(DP2$ref_Depth + DP2$alt_Depth )
    newdf <- cbind(newdf,SNP_index)
    colnames(newdf)[c(n+8)] <- paste0("AF_",colnames(AD)[2])
    newdf$DP_1 <- DP1$ref_Depth + DP1$alt_Depth
    colnames(newdf)[c(n+9)] <- paste0("DP_",colnames(AD)[1])
    newdf$DP_2 <- DP2$ref_Depth + DP2$alt_Depth
    colnames(newdf)[c(n+10)] <- paste0("DP_",colnames(AD)[2])
    colnames(newdf) <- gsub(sample_wt,"WT",colnames(newdf))
    colnames(newdf) <- gsub(sample_mut,"MUT",colnames(newdf))
    newdf$PL_MUT <- sapply(newdf$PL_MUT, function(x) {
      vals <- as.numeric(unlist(strsplit(x, ",")))
      sort(vals)[2]
    })
    newdf$PL_WT <- sapply(newdf$PL_WT, function(x) {
      vals <- as.numeric(unlist(strsplit(x, ",")))
      sort(vals)[2]
    })
    newdf$RTA <- paste0(newdf$REF, ">", newdf$ALT)
    newdf$EMS_is = ifelse(newdf$RTA %in% c("G>A", "C>T"), 1, 0)
    newdf$Pr_change_is = ifelse(newdf$protein_change %in% "", 0, 1)
    newdf$Pr_change_is[newdf$mutation_effect %in% "synonymous_variant"] <- 0
    ad_flt <- newdf[,c("AD_WT","AD_MUT")]
    ad_flt[ad_flt == "."] = "1,0"
    newdf$ED4 <- apply(ad_flt, 1, function(x){
      count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
      depth1 <- count[1] + count[2]
      depth2 <- count[3] + count[4]

      ED <- sqrt((count[3] / depth2 - count[1] / depth1)^2 +
                   (count[4] / depth2- count[2] /depth1)^2)
      return(ED^4)

    })
    newdf$detaAF <- newdf$AF_MUT - newdf$AF_WT
    #去掉有空值或NA的行
    df <- newdf %>% tidyr::drop_na()
    df <- df %>% filter(!GT_MUT %in% "" )
    df <- df %>% filter(!GT_WT %in% "" )
    # 计算2Mb内的SNP数量
    # df$SNP_count_2Mb <- sapply(1:nrow(df), function(i) {
    #   chrom <- df$CHROM[i]
    #   pos <- df$POS[i] %>% as.numeric()
    #   start_pos <- max(1, pos - 1000000)
    #   end_pos <- pos + 1000000
    #   count <- sum(df$CHROM == chrom & as.numeric(df$POS) >= start_pos & as.numeric(df$POS) <= end_pos)
    #   return(count)
    # })
    library(doParallel)
    library(foreach)

    # 注册并行后端
    registerDoParallel(cores = 8)

    # 并行计算
    df$SNP_count_2Mb <- foreach(i = 1:nrow(df), .combine = c) %dopar% {
      chrom <- df$CHROM[i]
      pos <- as.numeric(df$POS[i])
      start_pos <- max(1, pos - 1000000)
      end_pos <- pos + 1000000
      sum(df$CHROM == chrom & as.numeric(df$POS) >= start_pos & as.numeric(df$POS) <= end_pos)
    }
    # 停止集群
    stopImplicitCluster()
    return(df)
  }

}



# 使用旧机器学习模型
# get_results_ml <- function(vcffile, outdir, sample_wt, sample_mut, GT_mut,mindp=4, model_path = NULL) {
#   # Parameter validation
#   if (missing(vcffile) || !file.exists(vcffile)) {
#     stop("VCF file is required and must exist")
#   }
#
#   if (missing(outdir)) {
#     stop("Output directory is required")
#   }
#
#   if (missing(sample_wt) || missing(sample_mut)) {
#     stop("Sample names for wild-type and mutant are required")
#   }
#
#   if (GT_mut != "aa") {
#     stop("ML model currently only supports 'aa' mutation type")
#   }
#
#   # Create output directory if it doesn't exist
#   if (!dir.exists(outdir)) {
#     dir.create(outdir, recursive = TRUE)
#   }
#
#   # Check for required packages
#   required_pkgs <- c("vcfR", "data.table", "dplyr", "tidymodels", "openxlsx",
#                      "caret", "class", "ggplot2", "tidymodels")
#   missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
#
#   if (length(missing_pkgs) > 0) {
#     stop(paste("The following required packages are missing:",
#                paste(missing_pkgs, collapse = ", ")))
#   }
#
#   # Read and process VCF file
#   vcf <- vcfR::read.vcfR(vcffile)
#   df <- as.data.frame(vcf@fix)
#   df$ID <- paste(df$CHROM, df$POS, sep = "_")
#   AD <- vcfR::extract.gt(vcf, "AD")
#   GT <- vcfR::extract.gt(vcf, "GT")
#
#   # Extract INFO fields
#   extract_info_field <- function(info, field_num) {
#     sapply(strsplit(info, "|", fixed = TRUE), function(x) {
#       if (length(x) >= field_num) x[field_num] else NA_character_
#     })
#   }
#
#   df <- df %>%
#     dplyr::mutate(
#       mutation_effect = extract_info_field(.data$INFO, 2),
#       mutation_effect_level = extract_info_field(.data$INFO, 3),
#       Gene = extract_info_field(.data$INFO, 4),
#       CDS_change = extract_info_field(.data$INFO, 10),
#       protein_change = extract_info_field(.data$INFO, 11)
#     )
#
#   # Select relevant columns
#   newdf <- df %>%
#     dplyr::select(CHROM, POS, ID, REF, ALT, mutation_effect,
#                   mutation_effect_level, CDS_change, protein_change, Gene)
#
#   # Process genotype data
#   if (ncol(GT) == 2) {
#     # Add GT and AD data
#     newdf <- cbind(
#       newdf,
#       GT = GT,
#       AD = AD
#     )
#     colnames(newdf)[(ncol(newdf)-3):ncol(newdf)] <- c(
#       paste0("GT_", colnames(GT)),
#       paste0("AD_", colnames(AD))
#     )
#
#     # Filter out variants with multiple ALT alleles
#     multi_alt <- grepl(",", newdf$ALT)
#     if (any(multi_alt)) {
#       newdf <- newdf[!multi_alt, ]
#       AD <- AD[!multi_alt, ]
#       GT <- GT[!multi_alt, ]
#     }
#
#     # Calculate SNP indices
#     calc_snp_index <- function(ad_col) {
#       dp <- data.table::fread(text = ad_col, header = FALSE, sep = ",",
#                               data.table = FALSE, fill = TRUE)
#       dp <- dp[, 1:2]  # Ensure only two columns
#       colnames(dp) <- c("ref_Depth", "alt_Depth")
#       dp$alt_Depth / (dp$ref_Depth + dp$alt_Depth)
#     }
#
#     newdf <- newdf %>%
#       dplyr::mutate(
#         SNP_index_WT = calc_snp_index(.data[[paste0("AD_", sample_wt)]]),
#         SNP_index_MUT = calc_snp_index(.data[[paste0("AD_", sample_mut)]])
#       )
#
#     # Rename columns
#     colnames(newdf) <- gsub(sample_wt, "WT", colnames(newdf))
#     colnames(newdf) <- gsub(sample_mut, "MUT", colnames(newdf))
#   }
#
#   # Load model
#   if (is.null(model_path)) {
#     model_file <- system.file("extdata", "final.lightgbm_fit.rds",
#                               package = "BSAliulab")
#     if (!file.exists(model_file)) {
#       stop("Built-in model file not found in package")
#     }
#   } else if (!file.exists(model_path)) {
#     stop("Specified model file does not exist")
#   } else {
#     model_file <- model_path
#   }
#
#   lightgbm_fit <- readRDS(model_file) %>% workflowsets::extract_workflow()
#
#   # Prepare data for prediction
#   df <- newdf
#
#   ad_flt <- df[,c("AD_WT","AD_MUT")]
#   ad_flt[ad_flt == "."] = "1,0"
#   df$ED4 <- apply(ad_flt, 1, function(x){
#     count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
#     depth1 <- count[1] + count[2]
#     depth2 <- count[3] + count[4]
#
#     ED <- sqrt((count[3] / depth2 - count[1] / depth1)^2 +
#                  (count[4] / depth2- count[2] /depth1)^2)
#     return(ED^4)
#
#   })
#   df$DP_WT <- apply(ad_flt, 1, function(x){
#     count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
#     depth1 <- count[1] + count[2]
#     return(depth1)
#   })
#   df$DP_MUT <- apply(ad_flt, 1, function(x){
#     count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
#     depth2 <- count[3] + count[4]
#     return(depth2)
#   })
#   df <- filter(df, DP_WT >= 4 & DP_MUT >= 4)
#   df <- na.omit(df)
#   df$ProteinChange <- ifelse(df$protein_change %in% "", "FALSE", "TRUE")
#   df$GT_MUT <- gsub("/","|",df$GT_MUT)
#   df$GT_WT <- gsub("/","|",df$GT_WT)
#   gxdata <- df
#   gxdata$mutation_effect_level <- factor(gxdata$mutation_effect_level,levels=c("HIGH","MODERATE","MODIFIER","LOW"))
#   gxdata$gt_mut <- gxdata$GT_MUT
#   gxdata$gt_mut <- gsub("0|0","AA",gxdata$gt_mut,fixed=TRUE)
#   gxdata$gt_mut <- gsub("1|1","aa",gxdata$gt_mut,fixed=TRUE)
#   gxdata$gt_mut <- gsub("0|1","Aa",gxdata$gt_mut,fixed=TRUE)
#   gxdata$gt_mut <- gsub("1|0","Aa",gxdata$gt_mut,fixed=TRUE)
#   gxdata$gt_wt <- gxdata$GT_WT
#   gxdata$gt_wt <- gsub("0|0","AA",gxdata$gt_wt,fixed=TRUE)
#   gxdata$gt_wt <- gsub("1|1","aa",gxdata$gt_wt,fixed=TRUE)
#   gxdata$gt_wt <- gsub("0|1","Aa",gxdata$gt_wt,fixed=TRUE)
#   gxdata$gt_wt <- gsub("1|0","Aa",gxdata$gt_wt,fixed=TRUE)
#   gxdata$GT_MUT_WT <- paste(gxdata$gt_mut,gxdata$gt_wt,sep="|")
#   # "AA|aa" "aa|aa" "Aa|aa" "Aa|AA" "aa|AA" "AA|AA" "AA|Aa" "Aa|Aa" "aa|Aa"
#   gxdata$GT_MUT_WT <- factor(gxdata$GT_MUT_WT,levels=c("aa|AA","aa|Aa","Aa|AA","Aa|Aa","Aa|aa","AA|aa","AA|Aa","aa|aa","AA|AA"))
#   gxdata$ProteinChange <- factor(gxdata$ProteinChange,levels=c("TRUE","FALSE"))
#
#   df <- gxdata
#   # Make predictions
#   pred_lightgbm <- df %>%
#     dplyr::bind_cols(
#       stats::predict(lightgbm_fit, ., type = "prob"),
#       stats::predict(lightgbm_fit, ., type = "class")
#     )
#
#   f_lightgbm <- pred_lightgbm %>%
#     dplyr::filter(.data$.pred_class == TRUE)
#
#   # Save results to Excel
#   wb <- openxlsx::createWorkbook()
#   openxlsx::addWorksheet(wb, sheetName = "00_raw_dataframe", tabColour = "white")
#   openxlsx::addWorksheet(wb, sheetName = "01_Pred_LightGBM", tabColour = "grey")
#   openxlsx::addWorksheet(wb, sheetName = "02_Pred_LightGBM_TRUE", tabColour = "#A593E0")
#
#   openxlsx::writeData(wb, sheet = 1, df)
#   openxlsx::writeData(wb, sheet = 2, pred_lightgbm)
#   openxlsx::writeData(wb, sheet = 3, f_lightgbm)
#   outfile <- paste0(outdir, "/", gsub(".VQSR.snpEff.vcf.gz",".Predict_lightgbm.xlsx",basename(vcffile)))
#   openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
#
#   # Generate visualization
#   Pred_Final <- openxlsx::read.xlsx(outfile, sheet = "01_Pred_LightGBM")
#   Pred_Final$CHROM <-  gsub("chr","",Pred_Final$CHROM)
#   Pred_Final$POS <- as.numeric(Pred_Final$POS)
#   Pred_Final$CHROM <-  factor(Pred_Final$CHROM,levels = seq(1:12))
#   Pred_Final$level <- ifelse(Pred_Final$mutation_effect_level == "HIGH",4,ifelse(Pred_Final$mutation_effect_level == "MODERATE",3,ifelse(Pred_Final$mutation_effect_level == "MODIFIER",2,1)))
#
#   p4 <- ggplot2::ggplot(Pred_Final, ggplot2::aes(x = .data$POS/1e6, y = .data$.pred_TRUE,
#                                                  color = .data$level)) +
#     ggplot2::geom_point(size = 0.5) +
#     ggplot2::facet_grid(~.data$CHROM, scales = "free_x", space = "free_x", switch = "x") +
#     ggplot2::scale_color_gradientn(colours = c("#f9d423", "#f83600", "#020f75")) +
#     ggplot2::labs(title = "", x = "", y = "LightGBM") +
#     ggplot2::theme_bw() +
#     ggplot2::theme(
#       text = ggplot2::element_text(size = 8, color = "black"),
#       strip.background = ggplot2::element_blank(),
#       strip.placement = "outside",
#       axis.line.x = ggplot2::element_line(linetype = 1, color = "black", linewidth = 0.5),
#       axis.line.y = ggplot2::element_line(linetype = 1, color = "black", linewidth = 0.5),
#       panel.spacing = ggplot2::unit(0, "cm"),
#       axis.title.y = ggplot2::element_text(size = 8),
#       axis.text.x = ggplot2::element_blank(),
#       axis.text = ggplot2::element_text(size = 8),
#       axis.title = ggplot2::element_text(size = 8),
#       plot.title = ggplot2::element_text(size = 8),
#       legend.text = ggplot2::element_text(size = 8),
#       legend.title = ggplot2::element_text(size = 8),
#       legend.position = "none",
#       panel.border = ggplot2::element_blank(),
#       panel.grid.major = ggplot2::element_blank(),
#       panel.grid.minor = ggplot2::element_blank(),
#       axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
#       axis.ticks.length = ggplot2::unit(0.2, "cm"),
#       axis.ticks.x = ggplot2::element_blank(),
#       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
#       plot.margin = ggplot2::margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
#     )
#
#   pdf_file <- gsub(".xlsx", ".pdf", outfile)
#   grDevices::pdf(file = pdf_file, width = 4, height = 1.5)
#   print(p4)
#   grDevices::dev.off()
#
#   invisible(NULL)
# }
