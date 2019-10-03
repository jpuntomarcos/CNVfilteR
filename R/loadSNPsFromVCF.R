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
#' @param verbose Whether to show information messages. (Defaults to TRUE)
#'
#' @return A list where names are sample names, and values are \code{GRanges} objects containing the variants for each sample, including the following metadata columns:
#'   - \code{ref.support}: Reference allele depth field
#'   - \code{alt.support}: Alternative allele depth field
#'   - \code{alt.freq}: allelic frequency
#'   - \code{total.depth}: total depth
#'
#' @examples
#' vcf.path <- system.file("extdata", "variants.sample1.vcf.gz", package = "CNVfilteR", mustWork = TRUE)
#' vcf <- loadSNPsFromVCF(vcf.path)
#'
#' @import assertthat
#' @import VariantAnnotation
#' @importFrom GenomeInfoDb seqlevelsStyle seqinfo
#' @importFrom regioneR toGRanges filterChromosomes
#' @importFrom Rsamtools TabixFile headerTabix
#' @importFrom IRanges subsetByOverlaps reduce
#' @importFrom GenomicRanges mcols seqnames
#' @importFrom SummarizedExperiment rowRanges
#' @export loadSNPsFromVCF
#'
loadSNPsFromVCF <- function(vcf.file, vcf.source = NULL, ref.support.field = NULL, alt.support.field = NULL,
                            list.support.field = NULL, regions.to.filter = NULL, genome = "hg19",
                            exclude.non.canonical.chrs = TRUE, verbose = TRUE) {

  if (verbose) message("Scanning file ", vcf.file, "...")

  if(!is.null(regions.to.filter))
    regions.to.filter <- tryCatch(regioneR::toGRanges(regions.to.filter), error = function(e){stop("regions.to.filter must be in a valid format accepted by toGRanges.\n ", e)})


  ##  GET VCF SOURCE  ##

  # Get VCF source
  vcf.source <- auxGetVcfSource(vcf.source, vcf.file)


  ##  GET VCF FIELDS  ##

  # Detect format: recognise fields where allelic depth and frequency are stored
  support.field.is.list <- FALSE
  ref.and.forward.fields <- FALSE
  msg <- paste(vcf.source, "was found as source in the VCF metadata,")

  if (vcf.source == "VarScan2") {
    if (is.null(ref.support.field)) ref.support.field <- "RD"
    if (is.null(alt.support.field)) alt.support.field <- "AD"
    msg <- paste(msg, ref.support.field, "will be used as ref allele depth field,", alt.support.field, "will be used as alt allele depth field.")
  } else if (vcf.source %in% c("strelka", "freeBayes", "HaplotypeCaller", "UnifiedGenotyper")) {
    if (is.null(list.support.field)) ref.support.field <- alt.support.field <- list.support.field <- "AD"
    support.field.is.list <- TRUE
    msg <- paste(msg, list.support.field, "will be used as allele support field in a list format: ref allele, alt allele.")
  }

  if (verbose) message(msg)


  ##  PROCESS VCF VARIANTS  ##

  # Get variants according to detected desired regions.to.filter
  if (!is.null(regions.to.filter)) {

    # reduce (join) all regions to avoid duplicated positions when using which command
    regions.to.filter <- reduce(regions.to.filter)

    # Get seqnames from vcf
    availableSeqs <- headerTabix(vcf.file)$seqnames

    # Convert regions.to.filter to Ensembl / UCSC chr style if necessary
    if (grepl("chr", seqlevels(regions.to.filter)) & !grepl("chr", availableSeqs)){
      GenomeInfoDb::seqlevelsStyle(regions.to.filter) <- "Ensembl"
    } else if (!grepl("chr", seqlevels(regions.to.filter)) & grepl("chr", availableSeqs)){
      GenomeInfoDb::seqlevelsStyle(regions.to.filter) <- "UCSC"
    }

    # filter regions.to.filter by available seqnames to prevent an error when calling readVcf
    regions.to.filter <- regions.to.filter[as.character(seqnames(regions.to.filter)) %in% availableSeqs]

    # create ScanVcfParam
    scan.vcf.param <- ScanVcfParam(info=NA, geno = c(alt.support.field, ref.support.field), which = regions.to.filter)

  } else {
    scan.vcf.param <- ScanVcfParam(info=NA, geno = c(alt.support.field, ref.support.field))
  }


  # load variants
  vars <- readVcf(file=TabixFile(vcf.file), genome = genome, param = scan.vcf.param)
  if (length(vars) > 0) {
    GenomeInfoDb::seqlevelsStyle(vars) <- "UCSC"
  }

  # Exclude non cannonical chromosomes to avoid conflicts when comparing later with GRanges with different non canonical chromosomes
  if (exclude.non.canonical.chrs){
    vars <- regioneR::filterChromosomes(vars, organism = genome(vars)[1])
  }

  # process each sample
  samples <- colnames(vars)
  res <- list()
  for(s in samples) {

    v <- vars[,s]
    if (support.field.is.list) {

      # Get alt depth
      depths.list <- geno(v)[[alt.support.field]]

      # process it expecting 4 or 2 items
      if (ref.and.forward.fields) {
        depths.ref1 <- unlist(lapply(depths.list, "[", 1))
        depths.ref2 <- unlist(lapply(depths.list, "[", 2))
        depths.alt3 <- unlist(lapply(depths.list, "[", 3))
        depths.alt4 <- unlist(lapply(depths.list, "[", 4))
        depths.ref <- depths.ref1 + depths.ref2
        depths.alt <- depths.alt3 + depths.alt4
      } else {
        depths.ref <- unlist(lapply(depths.list, "[", 1))
        depths.alt <- unlist(lapply(depths.list, "[", 2))
      }

    } else {
      depths.alt <- geno(v)[[alt.support.field]]
      depths.ref <- geno(v)[[ref.support.field]]
    }


    # calculate allelic freq
    freqs.alt <- round(depths.alt / (depths.alt + depths.ref) * 100, 4)
    freqs.alt[is.nan(freqs.alt)] <- NA

    # calculate total depth
    total.depth <- depths.alt + depths.ref


    #Build the GRanges
    ref <- ref(v)
    alt <- unlist(alt(v))
    if(length(ref) != length(alt))
      stop("Multiple alt fields are not supported.")
    v <- SummarizedExperiment::rowRanges(v)
    mdata <- data.frame(ref = ref, alt = alt, ref.support = as.vector(depths.ref), alt.support = as.vector(depths.alt),
                        alt.freq = as.vector(freqs.alt), total.depth = as.vector(total.depth))
    GenomicRanges::mcols(v) <- mdata
    res[[s]] <- v
  }

  return(res)
}


#' auxGetVcfSource
#'
#' @description
#' Obtains VCF source from a given VCF file path and checks if the source is supported. Auxiliar function used by \code{loadSNPsFromVCF}.
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
  vcf.header <- scanVcfHeader(vcf.file)
  vcf.header.meta <- meta(vcf.header)

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

  # Check if vcf source was finally recognized
  supported.tools <- c("VarScan2", "strelka", "freeBayes", "HaplotypeCaller", "UnifiedGenotyper")
  if (!vcf.source %in% supported.tools)
    stop(paste("VCF source was not recognized. Expected options are:", utils::str(supported.tools)))

  return(vcf.source)
}
