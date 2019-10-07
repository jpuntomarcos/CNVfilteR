#' loadCNVcalls
#'
#' @description
#' Loads CNV calls from a csv/tsv file
#'
#' @details
#' Loads a csv/tsv file containing CNV calls, and transform it into a GRanges with \code{cnv} and \code{sample} metadata columns.
#'
#'
#' @param cnvs.file Path to csv/tsv file containing the CNV calls.
#' @param chr.column Which column stores the chr location of the CNV.
#' @param start.column Which column stores the start location of the CNV.
#' @param end.column Which column stores the end location of the CNV.
#' @param coord.column CNV location in the chr:start-end format. Example: "1:538001-540000". If NULL, \code{chr.column},
#' \code{start.column} and {end.column} columns will be used. (Defaults to NULL)
#' @param cnv.column Which column stores the type of CNV (deletion or duplication).
#' @param sample.column Which column stores the sample name.
#' @param gene.column Which columns store the gene or genes affected (optional). (Defaults to NULL)
#' @param deletion Text used in the \code{cnv.column} to represent deletion CNVs. (Defaults to "deletion")
#' @param duplication Text used in the \code{cnv.column} to represent duplication CNVs. (Defaults to "duplication")
#' @param sep Separator symbol to load the csv/tsv file. (Defaults to "\\t")
#' @param skip Number of rows that shoud be skipped when reading the csv/tsv file. (Defaults to 0)
#' @param genome The name of the genome. (Defaults to "hg19")
#' @param exclude.non.canonical.chrs Whether to exclude non canonical chromosomes (Defaults to TRUE)
#'
#' @return A \code{GRanges} with a range per each CNV and the metadata columns:
#'  - \code{cnv}: type of CNV, "duplication" or "deletion"
#'  - \code{sample}: sample name
#'
#' @examples
#' # Load CNVs data
#' cnvs.file <- system.file("extdata", "DECoN.CNVcalls.csv", package = "CNVfilteR", mustWork = TRUE)
#' cnvs.gr <- loadCNVcalls(cnvs.file = cnvs.file, chr.column = "Chromosome", start.column = "Start", end.column = "End", cnv.column = "CNV.type", sample.column = "Sample")
#'
#'
#' @import assertthat
#' @importFrom regioneR toGRanges filterChromosomes
#' @importFrom utils read.csv
#'
#' @export loadCNVcalls
#'
loadCNVcalls <- function(cnvs.file, chr.column, start.column, end.column, coord.column = NULL, cnv.column, sample.column,
                         gene.column = NULL, deletion = "deletion", duplication = "duplication",
                         sep = "\t", skip = 0, genome = "hg19", exclude.non.canonical.chrs = TRUE){

  # intial vars
  CNV_DF_COLUMNS <- c("chr", "start", "end", "cnv", "sample")

  # Check input
  assertthat::assert_that(assertthat::is.string(cnvs.file))
  assertthat::assert_that(assertthat::is.string(chr.column))
  assertthat::assert_that(assertthat::is.string(start.column))
  assertthat::assert_that(assertthat::is.string(end.column))
  assertthat::assert_that(assertthat::is.string(coord.column) || is.null(coord.column))
  assertthat::assert_that(assertthat::is.string(cnv.column))
  assertthat::assert_that(assertthat::is.string(sample.column))
  assertthat::assert_that(assertthat::is.string(deletion))
  assertthat::assert_that(assertthat::is.string(duplication))
  assertthat::assert_that(assertthat::is.string(sep))
  assertthat::assert_that(assertthat::is.number(skip))
  assertthat::assert_that(assertthat::is.string(genome))
  assertthat::assert_that(is.logical(exclude.non.canonical.chrs))

  # Read data
  cnvs.df <- utils::read.csv(cnvs.file, sep=sep, header=TRUE, stringsAsFactors = FALSE, skip = skip)

  # Rename and select needed columns
  colnames(cnvs.df)[which(names(cnvs.df) == cnv.column)] <- "cnv"
  colnames(cnvs.df)[which(names(cnvs.df) == sample.column)] <- "sample"

  if (is.null(coord.column)){
    colnames(cnvs.df)[which(names(cnvs.df) == chr.column)] <- "chr"
    colnames(cnvs.df)[which(names(cnvs.df) == start.column)] <- "start"
    colnames(cnvs.df)[which(names(cnvs.df) == end.column)] <- "end"
  } else {
    parts <- do.call(rbind, strsplit(cnvs.df$coordinates, ":|-"))
    cnvs.df$chr <- parts[,1]
    cnvs.df$start <- parts[,2]
    cnvs.df$end <- parts[,3]
  }


  # depending on gene.column provided or not...
  if (assertthat::is.string(gene.column)){
    colnames(cnvs.df)[which(names(cnvs.df) == gene.column)] <- "gene"
    CNV_DF_COLUMNS <- c(CNV_DF_COLUMNS, "gene")
  }

  # Select only desired columns
  cnvs.df <- cnvs.df[, CNV_DF_COLUMNS]

  # reassing deletion / duplication text values
  cnvs.df[cnvs.df$cnv == deletion, "cnv"] <- "deletion"
  cnvs.df[cnvs.df$cnv == duplication, "cnv"] <- "duplication"

  # Transform to GRanges
  cnvs.gr <- regioneR::toGRanges(cnvs.df, genome = genome)

  # Exclude non cannonical chromosomes to avoid conflicts when comparing later with GRanges with different non canonical chromosomes
  if (exclude.non.canonical.chrs){
    cnvs.gr <- regioneR::filterChromosomes(cnvs.gr, organism = genome)
  }


  return(cnvs.gr)
}

