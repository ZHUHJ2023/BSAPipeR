#' Set or get software paths for BSA analysis
#'
#' @param bwa Path to bwa executable
#' @param fastp Path to fastp executable
#' @param samtools Path to samtools executable
#' @param picard Path to picard jar file
#' @param bcftools Path to bcftools executable
#' @param gatk Path to gatk jar file
#' @param snpEff Path to snpEff directory
#' @param java Path to java executable
#' @param pandepth Path to pandepth executable
#'
#' @return A list of software paths (invisibly)
#' @export
set_software_paths <- function(bwa = NULL, fastp = NULL, samtools = NULL,
                               picard = NULL, bcftools = NULL, gatk = NULL,
                               snpEff = NULL, java = NULL,pandepth = NULL) {

  # 使用包环境而不是命名空间
  if (!exists(".BSApaths", envir = .GlobalEnv)) {
    assign(".BSApaths", new.env(), envir = .GlobalEnv)
  }

  # 获取当前路径
  current_paths <- if (exists("paths", envir = .BSApaths)) {
    get("paths", envir = .BSApaths)
  } else {
    list(
      bwa = "",
      fastp = "",
      samtools = "",
      picard = "",
      bcftools = "",
      gatk = "",
      snpEff = "",
      java = "",
      pandepth = ""
    )
  }

  # 更新路径
  if (!is.null(bwa)) current_paths$bwa <- bwa
  if (!is.null(fastp)) current_paths$fastp <- fastp
  if (!is.null(samtools)) current_paths$samtools <- samtools
  if (!is.null(picard)) current_paths$picard <- picard
  if (!is.null(bcftools)) current_paths$bcftools <- bcftools
  if (!is.null(gatk)) current_paths$gatk <- gatk
  if (!is.null(snpEff)) current_paths$snpEff <- snpEff
  if (!is.null(java)) current_paths$java <- java
  if (!is.null(pandepth)) current_paths$pandepth <- pandepth

  # 存储更新后的路径
  assign("paths", current_paths, envir = .BSApaths)

  invisible(current_paths)
}

#' Get currently configured software paths
#'
#' @return A list of software paths
#' @export
get_software_paths <- function() {
  if (exists(".BSApaths", envir = .GlobalEnv)) {
    if (exists("paths", envir = .BSApaths)) {
      return(get("paths", envir = .BSApaths))
    }
  }
  stop("Software paths not configured. Please use set_software_paths() first.")
}
