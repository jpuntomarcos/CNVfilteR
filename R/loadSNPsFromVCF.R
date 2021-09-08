#' loadSNPsFromVCF
#'
#' @description
#' Loads SNPs (SNVs/indels) from a VCF file
#'
#' @details
#' Given a VCF file path, the function recognizes the variant caller source to decide which fields should be used to calculate
#' ref/alt support and allelic frequency (see \code{return}). Current supported variant callers are VarScan2, Strelka/Strelka2, freebayes,
#' HaplotypeCaller and UnifiedGenotyper.
#'
#' Optionally, the fields where the data is stored can be manually set by using the parameters \code{ref.support.field},
#' \code{alt.support.field} and \code{list.support.field}
#'
#' Requirement: a TabixFile (.tbi) should exists in the same directory of the VCF file.
#'
#' @param vcf.file VCF file path
#' @param vcf.source VCF source, i.e., the variant caller used to generate the VCF file. If set, the function will not try to recognize the source. (Defaults to NULL)
#' @param ref.support.field Reference allele depth field. (Defaults to NULL)
#' @param alt.support.field Alternative allele depth field. (Defaults to NULL)
#' @param list.support.field Allele support field in a list format: reference allele, alternative allele. (Defaults to NULL)
#' @param regions.to.filter The regions to which limit the VCF import. It can be used to speed up the import process. (Defaults to NULL)
#' @param genome The name of the genome (Defaults to "hg19")
#' @param exclude.non.canonical.chrs Whether to exclude non canonical chromosomes (Defaults to TRUE)
#' @param verbose Whether to show information messages. (Defaults to TRUE)
#'
#' @return A list where names are sample names, and values are \code{GRanges} objects containing the variants for each sample, including the following metadata columns:
#' \itemize{
#'  \item \code{ref.support}: Reference allele depth field
#'  \item \code{alt.support}: Alternative allele depth field
#'  \item \code{alt.freq}: allelic frequency
#'  \item \code{total.depth}: total depth
#' }
#'
#' @examples
#' vcf.file <- system.file("extdata", "variants.sample1.vcf.gz", package = "CNVfilteR", mustWork = TRUE)
#' vcf <- loadSNPsFromVCF(vcf.file)
#'
#' @import assertthat
#' @import VariantAnnotation
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevels
#' @importFrom regioneR toGRanges filterChromosomes
#' @importFrom Rsamtools TabixFile headerTabix
#' @importFrom methods is
#' @importFrom IRanges reduce
#' @importFrom GenomicRanges mcols seqnames
#' @importFrom SummarizedExperiment rowRanges
#' @export loadSNPsFromVCF
#'
loadSNPsFromVCF <- function(vcf.file, vcf.source = NULL, ref.support.field = NULL, alt.support.field = NULL,
                            list.support.field = NULL, regions.to.filter = NULL, genome = "hg19",
                            exclude.non.canonical.chrs = TRUE, verbose = TRUE) {

  # Check input
  assertthat::assert_that(assertthat::is.string(vcf.file))
  assertthat::assert_that(assertthat::is.string(vcf.source) || is.null(vcf.source))
  assertthat::assert_that(assertthat::is.string(ref.support.field) || is.null(ref.support.field))
  assertthat::assert_that(assertthat::is.string(alt.support.field) || is.null(alt.support.field))
  assertthat::assert_that(assertthat::is.string(list.support.field) || is.null(list.support.field))
  assertthat::assert_that(methods::is(regions.to.filter, "GRanges") || is.null(regions.to.filter))
  assertthat::assert_that(assertthat::is.string(genome))
  assertthat::assert_that(is.logical(exclude.non.canonical.chrs))
  assertthat::assert_that(is.logical(verbose))

  if (verbose) message("Scanning file ", vcf.file, "...")

  if(!is.null(regions.to.filter))
    regions.to.filter <- tryCatch(regioneR::toGRanges(regions.to.filter), error = function(e){stop("regions.to.filter must be in a valid format accepted by toGRanges.\n ", e)})


  ##  GET VCF SOURCE  ##

  # Get VCF source
  vcf.source <- auxGetVcfSource(vcf.source, vcf.file)

  # Check if vcf source was finally recognized
  supported.tools <- c("VarScan2", "strelka", "freeBayes", "HaplotypeCaller", "UnifiedGenotyper")
  msg <- ""
  if (!vcf.source %in% supported.tools){
    if ( (is.null(ref.support.field) | is.null(alt.support.field)) & (is.null(list.support.field)) )  {
      stop(paste("VCF source was not recognized, and ref.support.field/alt.support.field/list.support.field were not provided. Expected options are:", utils::str(supported.tools)))
    } else {
      if (!is.null(list.support.field)) {
        msg <- paste("VCF source", vcf.source ,"is not supported, but list.support.field was provided. ")
      } else {
        msg <- paste("VCF source", vcf.source ,"is not supported, but ref.support.field/alt.support.field were provided. ")
      }
    }
  } else {
    msg <- paste(vcf.source, "was found as source in the VCF metadata,")
  }


  ##  GET VCF FIELDS  ##

  # Detect format: recognise fields where allelic depth and frequency are stored
  support.field.is.list <- FALSE
  ref.and.forward.fields <- FALSE

  if (vcf.source == "VarScan2") {
    if (is.null(ref.support.field)) ref.support.field <- "RD"
    if (is.null(alt.support.field)) alt.support.field <- "AD"
    msg <- paste(msg, ref.support.field, "will be used as ref allele depth field,", alt.support.field, "will be used as alt allele depth field.")
  } else if (vcf.source %in% c("strelka", "freeBayes", "HaplotypeCaller", "UnifiedGenotyper"))  {
    if (is.null(list.support.field))
      list.support.field <- "AD"
    msg <- paste(msg, list.support.field, "will be used as allele support field in a list format: ref allele, alt allele.")
  }

  # Set support.field.is.list
  if (!is.null(list.support.field)){
    support.field.is.list <- TRUE
  }

  # set genoField
  if (is.null(ref.support.field)){
    genoField <- c("GT", list.support.field)
  } else {
    genoField <- c("GT", ref.support.field, alt.support.field)
  }

  if (verbose) message(msg)


  ##  PROCESS VCF VARIANTS  ##

  # Get variants according to detected desired regions.to.filter
  if (!is.null(regions.to.filter)) {

    # reduce (join) all regions to avoid duplicated positions when using which command
    regions.to.filter <- IRanges::reduce(regions.to.filter)

    # Get seqnames from vcf
    availableSeqs <- Rsamtools::headerTabix(vcf.file)$seqnames

    # Convert regions.to.filter to Ensembl / UCSC chr style if necessary
    if (all(grepl("chr", seqlevels(regions.to.filter))) & all(!grepl("chr", availableSeqs))){
      GenomeInfoDb::seqlevelsStyle(regions.to.filter) <- "Ensembl"
    } else if (all(!grepl("chr", seqlevels(regions.to.filter))) & all(grepl("chr", availableSeqs))){
      GenomeInfoDb::seqlevelsStyle(regions.to.filter) <- "UCSC"
    }

    # filter regions.to.filter by available seqnames to prevent an error when calling readVcf
    regions.to.filter <- regions.to.filter[as.character(GenomicRanges::seqnames(regions.to.filter)) %in% availableSeqs]

    # create ScanVcfParam
    scan.vcf.param <- VariantAnnotation::ScanVcfParam(info=NA, geno = genoField, which = regions.to.filter)

  } else {
    scan.vcf.param <- VariantAnnotation::ScanVcfParam(info=NA, geno = genoField)
  }


  # load variants
  vars <- VariantAnnotation::readVcf(file=Rsamtools::TabixFile(vcf.file), genome = genome, param = scan.vcf.param)
  if (length(vars) > 0) {
    GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(SummarizedExperiment::rowRanges(vars))) <- "UCSC"
  }

  # process each sample
  samples <- colnames(vars)
  res <- list()
  for(s in samples) {

    v <- vars[,s]

    # Filter: only variants with GT denoting that the variant was found
    ids.to.select <- sapply(VariantAnnotation::geno(v)[["GT"]],
                            function(gt) gt %in% c("0/1", "0|1", "1/1", "1|1"))
    if (length(ids.to.select) > 0) {
      v <- v[ids.to.select]
    }


    if (support.field.is.list) {

      # Get alt depth
      depths.list <- VariantAnnotation::geno(v)[[list.support.field]]

      # process it expecting 4 or 2 items
      if (ref.and.forward.fields) {
        # Discard those variants with unexpected number of ref/alt values
        ids.to.select <- sapply(depths.list, function(i) length(i) == 4)
        v <- v[ids.to.select]
        depths.list <- depths.list[ids.to.select]

        # Get depths
        depths.ref1 <- unlist(lapply(depths.list, "[", 1))
        depths.ref2 <- unlist(lapply(depths.list, "[", 2))
        depths.alt3 <- unlist(lapply(depths.list, "[", 3))
        depths.alt4 <- unlist(lapply(depths.list, "[", 4))
        depths.ref <- depths.ref1 + depths.ref2
        depths.alt <- depths.alt3 + depths.alt4
      } else {

        # Discard those variants with unexpected number of ref/alt values
        ids.to.select <- sapply(depths.list, function(i) length(i) == 2)
        if (length(ids.to.select) > 0) {
          v <- v[ids.to.select]
          depths.list <- depths.list[ids.to.select]
        }

        # Get depths
        depths.ref <- unlist(lapply(depths.list, "[", 1))
        depths.alt <- unlist(lapply(depths.list, "[", 2))
      }

    } else {
      depths.alt <- unlist(VariantAnnotation::geno(v)[[alt.support.field]])
      depths.ref <- unlist(VariantAnnotation::geno(v)[[ref.support.field]])
    }

    # calculate allelic freq
    freqs.alt <- round(depths.alt / (depths.alt + depths.ref) * 100, 4)
    freqs.alt[is.nan(freqs.alt)] <- 0

    # calculate total depth
    total.depth <- depths.alt + depths.ref

    #Build the GRanges
    ref <- VariantAnnotation::ref(v)
    alt <- unlist(VariantAnnotation::alt(v))
    v <- SummarizedExperiment::rowRanges(v)
    mdata <- data.frame(ref = ref, alt = alt, ref.support = as.vector(depths.ref), alt.support = as.vector(depths.alt),
                        alt.freq = as.vector(freqs.alt), total.depth = as.vector(total.depth))
    GenomicRanges::mcols(v) <- mdata

    # Exclude non cannonical chromosomes to avoid conflicts when comparing later with GRanges with different non canonical chromosomes
    if (exclude.non.canonical.chrs){
      v <- regioneR::filterChromosomes(v, organism = genome)
    }

    # save variants
    res[[s]] <- v
  }

  return(res)
}


#' auxGetVcfSource
#'
#' @description
#' Obtains VCF source from a given VCF file path. Auxiliar function used by \code{loadSNPsFromVCF}.
#'
#' @param vcf.source VCF source. Leave NULL to allow the function to recognize it. Otherwise, the function will not try to recognize the source. (Defaults to NULL)
#' @param vcf.file VCF file path
#'
#' @return VCF source
#'
#'
#' @import assertthat
#' @import VariantAnnotation
#' @importFrom utils str
#'
auxGetVcfSource <- function(vcf.source = NULL, vcf.file){

  # Call scanVcfHeader() but controlling if knwon error occurs
  tryCatch(
    {
      vcf.header <- VariantAnnotation::scanVcfHeader(vcf.file)
    },
    error = function(cond){
      if (cond$message == "subscript out of bounds"){
        stop(paste0("There was an error when calling scanVcfHeader() with vcf.file=",
                   vcf.file, ": \"", cond$message,
                   "\". It's probably due to a VCF file with no variants."))
      } else {
        stop(cond$message)
      }
    }
  )

  vcf.header.meta <- VariantAnnotation::meta(vcf.header)

  # Get VCF source if necessary
  if (is.null(vcf.source)) {
    vcf.source <- vcf.header.meta$source[,1]
    if (!is.null(vcf.source) && is.na(vcf.source))
      vcf.source <- NULL
  }

  # if not found, look at alternative header fields
  if (is.null(vcf.source)){

    # Try to find GATK caller in raw header
    vcf.lines <- readLines(vcf.file)
    for (l in vcf.lines){
      if (grepl("GATKCommandLine.HaplotypeCaller", l)) {
        vcf.source <- "HaplotypeCaller"
        break
      } else if (grepl("GATKCommandLine.UnifiedGenotyper", l)) {
        vcf.source <- "UnifiedGenotyper"
        break
      } else if (grepl("#CHROM", l))
        break
    }
  }

  # other vcf checks
  if (is.null(vcf.source)) {
    stop("No VCF source was found at VCF file and not VCF source was defined as input param.")
  } else if (grepl("freeBayes", vcf.source)){
    vcf.source <- "freeBayes"
  }

  return(vcf.source)
}
