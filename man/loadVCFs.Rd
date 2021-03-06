% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadVCFs.R
\name{loadVCFs}
\alias{loadVCFs}
\title{loadVCFs}
\usage{
loadVCFs(
  vcf.files,
  sample.names = NULL,
  cnvs.gr,
  min.total.depth = 10,
  regions.to.exclude = NULL,
  vcf.source = NULL,
  ref.support.field = NULL,
  alt.support.field = NULL,
  list.support.field = NULL,
  homozygous.range = c(90, 100),
  heterozygous.range = c(28, 72),
  exclude.indels = TRUE,
  genome = "hg19",
  exclude.non.canonical.chrs = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{vcf.files}{vector of VCFs paths. Both .vcf and .vcf.gz extensions are allowed.}

\item{sample.names}{Sample names vector containing sample names for each \code{vcf.files}. If NULL, sample name will be obtained from the VCF sample column.  (Defaults to NULL)}

\item{cnvs.gr}{\code{GRanges} object containg CNV calls. Call \code{loadCNVcalls} to obtain it. Only those variants in regions affected by CNVs will be loaded to speed up the load.}

\item{min.total.depth}{Minimum total depth. Variants under this value will be excluded. (Defaults to 10)}

\item{regions.to.exclude}{A \code{GRanges} object defining the regions for which the variants should be excluded.
Useful for defining known difficult regions like pseudogenes where the allele frequency is not trustable. (Defaults to NULL)}

\item{vcf.source}{VCF source, i.e., the variant caller used to generate the VCF file. If set, the \code{loadSNPsFromVCF} function will not try to recognize the source. (Defaults to NULL)}

\item{ref.support.field}{Reference allele depth field. (Defaults to NULL)}

\item{alt.support.field}{Alternative allele depth field. (Defaults to NULL)}

\item{list.support.field}{Allele support field in a list format: reference allele, alternative allele. (Defaults to NULL)}

\item{homozygous.range}{Homozygous range. Variants not in the homozygous/heterozygous intervals will be excluded. (Defaults to \code{c(90, 100)})}

\item{heterozygous.range}{Heterozygous range. Variants not in the homozygous/heterozygous intervals will be excluded. (Defaults to \code{c(28, 72)})}

\item{exclude.indels}{Whether to exclude indels when loading the variants. TRUE is the recommended value given that indels frequency varies in a different way than SNVs. (Defaults to TRUE)}

\item{genome}{The name of the genome. (Defaults to "hg19")}

\item{exclude.non.canonical.chrs}{Whether to exclude non canonical chromosomes (Defaults to TRUE)}

\item{verbose}{Whether to show information messages. (Defaults to TRUE)}
}
\value{
A list where names are the sample names, and values are the \code{GRanges} objects for each sample.
}
\description{
Loads VCFs files
}
\details{
Loads VCF files and computes alt allele frequency for each variant. It uses
\code{\link{loadSNPsFromVCF}} function load the data and identify the
correct VCF format for allele frequency computation.

If sample.names is not provided, the sample names included in the VCF itself
will be used. Both single-sample and multi-sample VCFs are accepted, but when
multi-sample VCFs are used, sample.names parameter must be NULL.

If vcf is not compressed with bgzip, the function compresses it and generates
the .gz file. If .tbi file does not exist for a given VCF file, the function
also generates it. All files are generated in a temporary folder.
}
\note{
Important: Compressed VCF must be compressed with
[bgzip ("block gzip") from Samtools htslib](http://www.htslib.org/doc/bgzip.html)
and not using the standard Gzip utility.
}
\examples{
# Load CNVs data (required by loadVCFs to speed up the load process)
cnvs.file <- system.file("extdata", "DECoN.CNVcalls.csv", package = "CNVfilteR", mustWork = TRUE)
cnvs.gr <- loadCNVcalls(cnvs.file = cnvs.file, chr.column = "Chromosome", start.column = "Start", end.column = "End", cnv.column = "CNV.type", sample.column = "Sample")

# Load VCFs data
vcf.files <- c(system.file("extdata", "variants.sample1.vcf.gz", package = "CNVfilteR", mustWork = TRUE),
               system.file("extdata", "variants.sample2.vcf.gz", package = "CNVfilteR", mustWork = TRUE))
vcfs <- loadVCFs(vcf.files, cnvs.gr = cnvs.gr)


}
