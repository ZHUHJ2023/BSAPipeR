#!/usr/bin/env Rscript
#
# BSAPipeR Analysis Pipeline - Parameterized Version
# Supports both FASTQ-to-VCF pipeline and VCF-to-results only modes
#
# Usage examples:
# 1. For FASTQ input (full pipeline):
#    echo "bwa = '/path/to/bwa'
#    fastp = '/path/to/fastp'
#    samtools = '/path/to/samtools'
#    picard = '/path/to/java -XX:ParallelGCThreads=8 -Xmx30g -jar /path/to/picard.jar'
#    bcftools = '/path/to/bcftools'
#    gatk = '/path/to/java -XX:ParallelGCThreads=8 -Xmx30g -jar /path/to/gatk.jar', # gatk-4.4.0.0
#    snpEff = '/path/to/snpEff', #SnpEff directory, SnpEff-4.3t
#    java = '/path/to/java', # jdk-18.0.2.1
#    pandepth = '/path/to/pandepth' " > software_paths.config
#    # Run BSAPipeR
#    Rscript bsapipe.r --input_type=fastq --ref=/path/to/reference.fa --refbwa=/path/to/bwa_index_prefix --mutfq1=mutant_R1.fq.gz --mutfq2=mutant_R2.fq.gz --wtfq1=wildtype_R1.fq.gz --wtfq2=wildtype_R2.fq.gz --gff3=/path/to/annotation.gff3 --outdir=results --sample_mut=mutant --sample_wt=wildtype --GT_mut=aa --mindp=4 --config=software_paths.config
#
# 2. For VCF input (results only):
#    Rscript bsapipe.r --input_type vcf --outdir results \
#           --sample_mut mutant --sample_wt wildtype

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(BSAPipeR)
})

# Configure software paths (modify these paths as needed)

## software_paths.config
# # BSAPipeR Software Paths Configuration File
# # Edit the paths below to match your system installation
#
# # Alignment and mapping tools
# bwa = "/path/to/bwa"
# samtools = "/path/to/samtools"
#
# # Quality control and preprocessing
# fastp = "/path/to/fastp"
#
# # Variant calling and processing
# bcftools = "/path/to/bcftools"
#
# # Java-based tools (full command line including Java options)
# java = "/path/to/java"
# picard = "/path/to/java -XX:ParallelGCThreads=8 -Xmx30g -jar /path/to/picard.jar"
# gatk = "/path/to/java -XX:ParallelGCThreads=8 -Xmx30g -jar /path/to/gatk.jar"
#
# # Annotation tools
# snpEff = "/path/to/snpEff"
#
# # Depth calculation
# pandepth = "/path/to/pandepth"

# Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    # Input type selection
    make_option(c("-t", "--input_type"), type = "character", default = "fastq",
                help = "Input type: 'fastq' or 'vcf' [default: %default]",
                metavar = "TYPE"),

    # Common parameters for both modes
    make_option(c("-o", "--outdir"), type = "character", default = "BSAPipeR_results",
                help = "Output directory [default: %default]",
                metavar = "DIR"),
    make_option(c("-m", "--sample_mut"), type = "character", default = "mutant",
                help = "Mutant sample name [default: %default]",
                metavar = "NAME"),
    make_option(c("-w", "--sample_wt"), type = "character", default = "wildtype",
                help = "Wildtype sample name [default: %default]",
                metavar = "NAME"),
    make_option(c("-d", "--depthFT"), type = "character", default = "False",
                help = "Specify sequencing depth (True/False) [default: %default]",
                metavar = "BOOL"),
    make_option(c("--mdepth"), type = "integer", default = 10,
                help = "Minimum depth threshold [default: %default]",
                metavar = "INT"),
    make_option(c("-gt", "--GT_mut"), type = "character", default = "aa",
                help = "Mutant genotype [default: %default]",
                metavar = "GT"),
    make_option(c("--mindp"), type = "integer", default = 4,
                help = "Minimum depth for variant calling [default: %default]",
                metavar = "INT"),
    make_option(c("--config"), type = "character", default = "software_paths.config",
                help = "Configuration file for software paths",
                metavar = "FILE"),

    # FASTQ-specific parameters
    make_option(c("-r", "--ref"), type = "character", default = NULL,
                help = "Reference genome FASTA file (required for fastq mode)",
                metavar = "FILE"),
    make_option(c("-b", "--refbwa"), type = "character", default = NULL,
                help = "BWA index prefix (required for fastq mode)",
                metavar = "PREFIX"),
    make_option(c("--gff3"), type = "character", default = NULL,
                help = "Annotation GFF3 file (required for fastq mode)",
                metavar = "FILE"),
    make_option(c("--mutfq1"), type = "character", default = NULL,
                help = "Mutant R1 FASTQ file (required for fastq mode)",
                metavar = "FILE"),
    make_option(c("--mutfq2"), type = "character", default = NULL,
                help = "Mutant R2 FASTQ file (required for fastq mode)",
                metavar = "FILE"),
    make_option(c("--wtfq1"), type = "character", default = NULL,
                help = "Wildtype R1 FASTQ file (required for fastq mode)",
                metavar = "FILE"),
    make_option(c("--wtfq2"), type = "character", default = NULL,
                help = "Wildtype R2 FASTQ file (required for fastq mode)",
                metavar = "FILE"),

    # VCF-specific parameters
    make_option(c("--snp_vcf"), type = "character", default = NULL,
                help = "SNP VCF file (alternative to auto-detection for vcf mode)",
                metavar = "FILE"),

    # Optional parameters
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Print verbose output")

  )

  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

  # Parse arguments with error handling
  args <- tryCatch({
    parse_args(parser)
  }, error = function(e) {
    # If there's an error, print help and exit
    print_help(parser)
    stop(e$message)
  })

  return(args)
}


# Main function
main <- function() {
  cat(paste(rep("=", 80), collapse=""), "\n")
  cat("BSAPipeR Analysis Pipeline\n")
  cat(paste(rep("=", 80), collapse=""), "\n\n")

  # Parse command line arguments
  args <- parse_arguments()

  # Validate arguments based on input type
  if (args$input_type == "fastq") {
    cat("Running in FASTQ mode (full pipeline)...\n")

    # Check required arguments for FASTQ mode
    required_fastq <- c("ref", "refbwa", "mutfq1", "mutfq2", "wtfq1", "wtfq2")
    missing <- sapply(required_fastq, function(x) is.null(args[[x]]))
    if (any(missing)) {
      missing_args <- paste(required_fastq[missing], collapse = ", ")
      stop(paste("Missing required arguments for fastq mode:", missing_args))
    }

    # Run full BSAPipeR pipeline
    cat("Starting full analysis pipeline...\n")
    configure_software_paths(args$config)
    BSAPipeR(
      ref = args$ref,
      refbwa = args$refbwa,
      outdir = args$outdir,
      sample_mut = args$sample_mut,
      sample_wt = args$sample_wt,
      GFF3 = args$gff3,
      mutfq1 = args$mutfq1,
      mutfq2 = args$mutfq2,
      wtfq1 = args$wtfq1,
      wtfq2 = args$wtfq2,
      depthFT = args$depthFT,
      mdepth = args$mdepth,
      GT_mut = args$GT_mut,
      mindp = args$mindp
    )

  } else if (args$input_type == "vcf") {
    cat("Running in VCF mode (results generation only)...\n")

    # Run create_results function
    create_results_from_vcf(
      outdir = args$outdir,
      sample_mut = args$sample_mut,
      sample_wt = args$sample_wt,
      depthFT = args$depthFT,
      mdepth = args$mdepth,
      GT_mut = args$GT_mut,
      mindp = args$mindp,
      snp_vcf = args$snp_vcf
    )

  } else {
    stop(paste("Invalid input type:", args$input_type, ". Must be 'fastq' or 'vcf'."))
  }
}

# Run the main function with error handling
tryCatch({
  main()
}, error = function(e) {
  writeLines(c(
    paste(rep("*", 80), collapse=""),
    paste("ERROR:", e$message),
    paste(rep("*", 80), collapse=""),
    "Use --help for usage information"
  ))
  quit(status = 1)
})
