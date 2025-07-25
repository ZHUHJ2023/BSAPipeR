% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step1.run_vcf_processing.R
\name{run_vcf_processing}
\alias{run_vcf_processing}
\title{run_vcf_processing}
\usage{
run_vcf_processing(
  ref,
  refbwa,
  outdir,
  sample_mut,
  sample_wt,
  GFF3,
  mutfq1,
  mutfq2,
  wtfq1,
  wtfq2,
  depthFT = "False",
  mdepth = 0
)
}
\arguments{
\item{mdepth}{}
}
\description{
run_vcf_processing
}
\examples{
set_software_paths(bwa = "/path/to/bwa",fastp = "/path/to/fastp",samtools = "/path/to/samtools",picard = "/path/to/picard.jar",bcftools = "/path/to/bcftools",gatk = "/path/to/gatk.jar",snpEff = "/path/to/snpEff") run_vcf_processing(ref = "/path/to/reference.fa",refbwa = "/path/to/bwa_index_prefix",workdir = "/path/to/workdir",sample_mut = "MUT",sample_wt = "WT",mutfq1 = "/path/to/mut_1.fq.gz",mutfq2 = "/path/to/mut_2.fq.gz",wtfq1 = "/path/to/wt_1.fq.gz",wtfq2 = "/path/to/wt_2.fq.gz")
}
