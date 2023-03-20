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
#' @param sample.name Sample name for all CNVs defined in \code{cnvs.file}. If set, \code{sample.column} is ignored (Defaults to NULL)
#' @param gene.column Which columns store the gene or genes affected (optional). (Defaults to NULL)
#' @param deletion Text used in the \code{cnv.column} to represent deletion CNVs. Multiple values are also allowed, for example: c("CN0", "CN1"). (Defaults to "deletion")
#' @param duplication Text used in the \code{cnv.column} to represent duplication CNVs.  Multiple values are also allowed, for example: c("CN3", "CN4") (Defaults to "duplication")
#' @param ignore.unexpected.rows Whether to ignore the rows which CNV \code{cnv.column} value is different to \code{deletion} or \code{duplication} values (Defaults to FALSE). It is useful for processing output from callers like LUMPY or Manta (they call also events that are not CNVs)
#' @param sep Separator symbol to load the csv/tsv file. (Defaults to "\\t")
#' @param skip Number of rows that should be skipped when reading the csv/tsv file. (Defaults to 0)
#' @param genome The name of the genome. (Defaults to "hg19")
#' @param exclude.non.canonical.chrs Whether to exclude non canonical chromosomes (Defaults to TRUE)
#' @param check.names.cnvs.file Whether to check \code{cnvs.file} names or not (Defaults to FALSE).  If TRUE then column names in the \code{cnvs.file} are checked to ensure that they are syntactically valid variable names. If necessary they are adjusted (by make.names) so that they are, and also to ensure that there are no duplicates
#'
#' @return A \code{GRanges} with a range per each CNV and the metadata columns:
#' \itemize{
#'  \item \code{cnv}: type of CNV, "duplication" or "deletion"
#'  \item \code{sample}: sample name
#' }
#'
#' Returns NULL if \code{cnvs.file} has no CNVs
#'
#' @examples
#' # Load CNVs data
#' cnvs.file <- system.file("extdata", "DECoN.CNVcalls.csv", package = "CNVfilteR", mustWork = TRUE)
#' cnvs.gr <- loadCNVcalls(cnvs.file = cnvs.file, chr.column = "Chromosome", start.column = "Start", end.column = "End", cnv.column = "CNV.type", sample.column = "Sample")
#'
#'
#' @import assertthat
#' @importFrom Biostrings head
#' @importFrom regioneR toGRanges filterChromosomes
#' @importFrom utils read.csv
#'
#' @export loadCNVcalls
#'
loadCNVcalls <- function(cnvs.file, chr.column, start.column, end.column,
                         coord.column = NULL, cnv.column, sample.column,
                         sample.name = NULL, gene.column = NULL,
                         deletion = "deletion", duplication = "duplication",
                         ignore.unexpected.rows = FALSE,
                         sep = "\t", skip = 0, genome = "hg19",
                         exclude.non.canonical.chrs = TRUE,  check.names.cnvs.file=FALSE){

  # intial vars
  CNV_DF_COLUMNS <- c("chr", "start", "end", "cnv", "sample")

  # Check input
  assertthat::assert_that(assertthat::is.string(cnvs.file))
  if (is.null(coord.column)){
    assertthat::assert_that(assertthat::is.string(chr.column))
    assertthat::assert_that(assertthat::is.string(start.column))
    assertthat::assert_that(assertthat::is.string(end.column))
  } else {
    assertthat::assert_that(assertthat::is.string(coord.column))
  }
  if (is.null(sample.name)){
    assertthat::assert_that(assertthat::is.string(sample.column))
  }
  assertthat::assert_that(assertthat::is.string(cnv.column))
  assertthat::assert_that(is.character(deletion))
  assertthat::assert_that(is.character(duplication))
  assertthat::assert_that(assertthat::is.string(sep))
  assertthat::assert_that(assertthat::is.number(skip))
  assertthat::assert_that(assertthat::is.string(genome))
  assertthat::assert_that(is.logical(exclude.non.canonical.chrs))
  assertthat::assert_that(is.logical(check.names.cnvs.file))

  # Read data
  cnvs.df <- utils::read.csv(cnvs.file, sep=sep, header=TRUE,
                             check.names = check.names.cnvs.file,
                             stringsAsFactors = FALSE, skip = skip)

  if (nrow(cnvs.df) == 0){
    message(paste("Warning: ", cnvs.file," has no rows"))
    return(NULL)
  }

  # Rename and select needed columns
  colnames(cnvs.df)[which(names(cnvs.df) == cnv.column)] <- "cnv"
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

  # Set sample name column
  if (is.null(sample.name)){
    colnames(cnvs.df)[which(names(cnvs.df) == sample.column)] <- "sample"
  } else {
    cnvs.df$sample <- sample.name
  }

  # Force sample to be character: it will avoid problems when in future steps
  cnvs.df$sample <- as.character(cnvs.df$sample)

  # Assert CNV type values are correct
  cnvValuesCheck <- unique(cnvs.df$cnv) %in%  c(deletion, duplication)
  if(!identical(unique(cnvValuesCheck), TRUE)) {
    if (!ignore.unexpected.rows){
      stop(paste0("Values found in ", cnv.column, " column are ",
                 toString(Biostrings::head(unique(cnvs.df$cnv))), " but expected values are only '",
                 toString(c(deletion, duplication)),
                 "'. Please, use the duplication and deletion parameters to select the expected values for the ",
                 cnv.column ," column, or set ignore.unexpected.rows parameter to TRUE"))
    }
  }


  # depending on gene.column provided or not...
  if (assertthat::is.string(gene.column)){
    colnames(cnvs.df)[which(names(cnvs.df) == gene.column)] <- "gene"
    CNV_DF_COLUMNS <- c(CNV_DF_COLUMNS, "gene")
  }

  # Select only desired columns
  cnvs.df <- cnvs.df[, CNV_DF_COLUMNS]

  # reassign deletion / duplication text values
  cnvs.df[cnvs.df$cnv %in% deletion, "cnv"] <- "deletion"
  cnvs.df[cnvs.df$cnv %in% duplication, "cnv"] <- "duplication"

  # Discard those rows with a CNV type different to the expected ones
  cnvs.df <- cnvs.df[cnvs.df$cnv %in% c("deletion", "duplication"),]

  # Transform to GRanges
  cnvs.gr <- regioneR::toGRanges(cnvs.df, genome = genome)

  # Exclude non cannonical chromosomes to avoid conflicts when comparing later with GRanges with different non canonical chromosomes
  if (exclude.non.canonical.chrs){
    cnvs.gr <- regioneR::filterChromosomes(cnvs.gr, organism = genome)
  }


  return(cnvs.gr)
}

