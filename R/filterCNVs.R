#' filterCNVs
#'
#' @description
#' Identifies those copy number calls that can be filtered out
#'
#' @details
#' Checks all the variants (SNV and optionally INDELs) in each CNV present in \code{cnvs.gr} to decide whether a CNV can be filtered out or not.
#' It returns an S3 object with 3 elments: \code{cnvs}, \code{variantsForEachCNV} and \code{filterParameters}. See return section for further details.
#'
#' A CNV deletion can be filtered out if there is at least \code{ht.deletions.threshold}% of heterozygous variants in the CNV.
#' A CNV duplication can be filtered out if the \code{score} is >= \code{dup.threshold.score} after computing all heterozygous variants falling in the CNV.
#'
#' If a CNV can be filtered out, then the value TRUE is set in the \code{filter} column of the \code{cnvs} element.
#'
#' @param cnvs.gr \code{GRanges} containing CNVs to be filtered out. Use \code{loadCNVcalls} to load them.
#' @param vcfs List of \code{GRanges} containing all variants (SNV/indel) obtaining with the \code{loadVCFs} function.
#' @param expected.ht.mean Expected heterozygous SNV/indel allele frequency (defaults to 50)
#' @param expected.dup.ht.mean1 Expected heterozygous SNV/indel allele frequency when the variant IS NOT in the same allele than the CNV duplication call. (defaults to 33.3)
#' @param expected.dup.ht.mean2 Expected heterozygous SNV/indel allele frequency when the variant IS in the same allele than the CNV duplication call. (defaults to 66.6)
#' @param sigmoid.c1 Sigmoid c1 parameter. (defaults to 2)
#' @param sigmoid.c2.vector Vector containing sigmoid c2 parameters for the six sigmoids functions. (defaults to c(28, 38.3, 44.7, 55.3, 61.3, 71.3))
#' @param dup.threshold.score Limit value to decide if a CNV duplication can be filtered out or not. A CNV duplication can be filtered out if the total score computed from heterozygous variants in the CNV is equal or greater than \code{dup.threshold.score}.  (defaults to 0.5)
#' @param ht.deletions.threshold Minimum percentage of heterozygous variants falling in a CNV deletion to filter that CNV. (defaults to 30)
#' @param margin.pct Variants in the CNV but close to the ends of the CNV will be ignored. \code{margin.pct} defines the percentage
#'  of CNV length, located at each CNV limit, where variants will be ignored. For example, for a CNV chr1:1000-2000 and
#'  a \code{margin.pct} value of 10, variants within chr1:1000-1100 and chr1:1900-2000 will be ignored.
#' @param verbose Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A S3 object with 3 elements:
#' \itemize{
#' \item \code{cnvs}: \code{GRanges} with the input CNVs and the meta-columns added during the call:
#'   \itemize{
#'    \item \code{cnv.id}: CNV id
#'    \item \code{filter}: Set to TRUE if the CNV can be filtered out
#'    \item \code{n.total.variants}: Number of variants in the CNV
#'    \item \code{n.hm.variants}: Number of homozygous variants. They do not give any evidenced for confirming or discarding the CNV.
#'    \item \code{n.ht.discard.CNV}: For a CNV duplication, number of heterozygous variants in that discard the CNV (those with a positive score)
#'    \item \code{n.ht.confirm.CNV}: For a CNV duplication, number of heterozygous variants that confirm the CNV (those with a negative score)
#'    \item \code{ht.pct}: Percentage of heterozygous variants for deletion CNVs
#'    \item \code{score}: total score when computing all the variants scores
#'   }
#' \item \code{variantsForEachCNV}: named list where each name correspond to a CNV id and the value is a \code{data.frame} with all variants falling in that CNV
#' \item \code{filterParameters}: input parameters used for filtering
#'}
#'
#' @examples
#' # Load CNVs data
#' cnvs.file <- system.file("extdata", "DECoN.CNVcalls.csv", package = "CNVfilteR", mustWork = TRUE)
#' cnvs.gr <- loadCNVcalls(cnvs.file = cnvs.file, chr.column = "Chromosome", start.column = "Start", end.column = "End", cnv.column = "CNV.type", sample.column = "Sample")
#'
#' # Load VCFs data
#' vcf.files <- c(system.file("extdata", "variants.sample1.vcf.gz", package = "CNVfilteR", mustWork = TRUE),
#'                system.file("extdata", "variants.sample2.vcf.gz", package = "CNVfilteR", mustWork = TRUE))
#' vcfs <- loadVCFs(vcf.files, cnvs.gr = cnvs.gr)
#'
#' # Filter CNVs
#' results <- filterCNVs(cnvs.gr, vcfs)
#'
#' # Check CNVs that can be filtered out
#' as.data.frame(results$cnvs[results$cnvs$filter == TRUE])
#'
#'
#' @import assertthat
#' @importFrom IRanges subsetByOverlaps IRanges
#' @importFrom GenomicRanges mcols GRanges
#' @importFrom methods is
#' @importFrom regioneR toGRanges
#' @export filterCNVs
#'
filterCNVs <- function(cnvs.gr, vcfs, expected.ht.mean = 50, expected.dup.ht.mean1 = 33.3, expected.dup.ht.mean2 = 66.6,
                       sigmoid.c1 = 2, sigmoid.c2.vector = c(28, 38.3, 44.7, 55.3, 61.3, 71.3),
                       dup.threshold.score = 0.5, ht.deletions.threshold = 30, verbose = FALSE,
                       margin.pct = 0) {

  # Check input
  assertthat::assert_that(methods::is(cnvs.gr, "GRanges"))
  assertthat::assert_that(is.list(vcfs))
  for (v in vcfs){
    if (!methods::is(v, "GRanges"))
      stop("vcfs parameter should be a list of GRanges")
  }
  assertthat::assert_that(assertthat::is.number(expected.ht.mean))
  assertthat::assert_that(assertthat::is.number(expected.dup.ht.mean1))
  assertthat::assert_that(assertthat::is.number(expected.dup.ht.mean2))
  assertthat::assert_that(assertthat::is.number(sigmoid.c1))
  assertthat::assert_that(is.numeric(sigmoid.c2.vector) && length(sigmoid.c2.vector) == 6)
  assertthat::assert_that(assertthat::is.number(dup.threshold.score))
  assertthat::assert_that(assertthat::is.number(ht.deletions.threshold))
  assertthat::assert_that(is.logical(verbose))

  # intersections between two sigmoids curves. For later use
  sigmoid.int1 <- mean(c(expected.dup.ht.mean1, expected.ht.mean))
  sigmoid.int2 <- mean(c(expected.dup.ht.mean2, expected.ht.mean))

  variantsForEachCNV <- list()

  # Filter CNVs depending on variant calls.
  cnvs.df <- as.data.frame(cnvs.gr)  # to speed up the process
  cnvs.df$cnv.id <- ''
  cnvs.df$filter <- ''
  cnvs.df$n.total.variants <- ''
  cnvs.df$n.hm.variants <- ''
  cnvs.df$n.ht.discard.CNV <- ''
  cnvs.df$n.ht.confirm.CNV <- ''
  cnvs.df$ht.pct <- ''
  cnvs.df$score <- ''
  nCNVs <- nrow(cnvs.df)
  for (i in seq_len(nCNVs)){

    cnvId <- toString(i)
    cnvs.df[i, "cnv.id"] <- cnvId

    if (cnvs.df[i, "sample"] %in% names(vcfs)){
      variants <- vcfs[[cnvs.df[i, "sample"]]]
      matchingVariants <- IRanges::subsetByOverlaps(variants, cnvs.gr[i],  type = "any")

      # Filter variants by pct lateral margin
      lateralMargin <- margin.pct / 100.0 * cnvs.df[i, "width"]
      regions.to.include <- GenomicRanges::GRanges(seqnames = cnvs.df[i, "seqnames"],
                                    ranges = IRanges::IRanges(cnvs.df[i, "start"] + lateralMargin, cnvs.df[i, "end"] - lateralMargin))
      matchingVariants <- IRanges::subsetByOverlaps(matchingVariants, regions.to.include,  type = "any")
      matchingVariants.df <- as.data.frame(matchingVariants)  # to speed up the process

      # count number of ht and hm
      htMatchingVariants.indexes <- which(matchingVariants.df$type == "ht")
      htMatchingVariants <- matchingVariants.df[htMatchingVariants.indexes,]
      nHT <- nrow(htMatchingVariants)
      nHM <- nrow(matchingVariants.df[matchingVariants.df$type == "hm",])
      cnvs.df[i, "n.total.variants"] <- nrow(matchingVariants.df)
      cnvs.df[i, "n.hm.variants"] <- nHM

      # if CNV is deletion & ht / (ht + hm)  > ht.deletions.threshold > discard CNV
      if (cnvs.df[i, "cnv"] == "deletion" && nHT > 0 && (nHT / (nHT + nHM) > (ht.deletions.threshold / 100.0) )) {
        cnvs.df[i, "filter"] <- TRUE
        if (verbose){
          message(paste0("CNV deletion at ", toString(cnvs.gr[i]) ," for sample ", cnvs.df[i, "sample"], " can be filtered"))
        }
        cnvs.df[i, "n.ht.discard.CNV"] <- nHT
        cnvs.df[i, "ht.pct"] <- nHT / nrow(matchingVariants.df) * 100.0
      }

      # Add score column for those variants mathching a duplication CNV
      if (cnvs.df[i, "cnv"] == "duplication" & nrow(matchingVariants.df) > 0){
        matchingVariants.df$score <- 0
      }

      # if CNV is duplication  & exists ht
      if (cnvs.df[i, "cnv"] == "duplication" & nHT > 0) {

        total.score <- 0
        discardCNV <- 0
        confirmCNV <- 0
        for (j in seq_len(nHT)){

          # Score variant
          score <- getVariantScore(
            freq = htMatchingVariants[j, "alt.freq"],
            expected.ht.mean = expected.ht.mean,
            expected.dup.ht.mean1 = expected.dup.ht.mean1,
            expected.dup.ht.mean2 = expected.dup.ht.mean2,
            sigmoid.int1 = sigmoid.int1, sigmoid.int2 = sigmoid.int2,
            sigmoid.c1 = sigmoid.c1, sigmoid.c2.vector = sigmoid.c2.vector)

          if(score > 0){
            discardCNV <- discardCNV + 1
          } else if(score < 0) {
            confirmCNV <- confirmCNV + 1
          }

          # compute total score and save score
          total.score <- total.score + score
          matchingVariants.df[htMatchingVariants.indexes[j], "score"] <- score
        }

        # Discard CNV depending on dup.threshold.score
        if (total.score >= dup.threshold.score){
          cnvs.df[i, "filter"] <- TRUE
          if (verbose){
            message(paste0("CNV duplication at ", toString(cnvs.gr[i]) ,
                           " for sample ", cnvs.df[i, "sample"],
                           " can be filtered, score: ", round(total.score, 4)))
          }
        }

        cnvs.df[i, "score"] <- total.score
        cnvs.df[i, "n.ht.confirm.CNV"] <- confirmCNV
        cnvs.df[i, "n.ht.discard.CNV"] <- discardCNV
      }

      if (nrow(matchingVariants.df) > 0)
        variantsForEachCNV[[cnvId]] <- matchingVariants.df
    }
  }

  # Calculate % of CNVs to be filtered
  nFiltered <- nrow(cnvs.df[cnvs.df$filter == TRUE,])
  pct <- round(nFiltered/nCNVs*100, 2)
  if (verbose){
    message(paste0(nFiltered, " of ", nCNVs, " (", pct ,"%) CNVs can be filtered"))
  }

  # Calculate % with matching variants to be filtered
  nCNVsWithMatchingVariants <- nrow(cnvs.df[cnvs.df$n.total.variants > 0 ,])
  pct <- round(nFiltered/nCNVsWithMatchingVariants*100, 2)
  if (verbose){
    message(paste0(nFiltered, " of ", nCNVsWithMatchingVariants, " (", pct ,"%) CNVs with overlapping SNVs can be filtered"))
  }

  # convert again data.frame to GRanges
  cnvs.gr.final <- regioneR::toGRanges(cnvs.df)
  GenomicRanges::mcols(cnvs.gr.final) <- GenomicRanges::mcols(cnvs.gr.final)[,  -which(names(GenomicRanges::mcols(cnvs.gr.final)) %in% c("strand","width"))]

  # Create S3 object and return
  filterParameters <- list(expected.ht.mean = expected.ht.mean, expected.dup.ht.mean1 = expected.dup.ht.mean1,
                           expected.dup.ht.mean2 = expected.dup.ht.mean2, sigmoid.c1 = sigmoid.c1, sigmoid.c2.vector = sigmoid.c2.vector,
                           dup.threshold.score = dup.threshold.score, ht.deletions.threshold = ht.deletions.threshold)
  resultsObject <- list(cnvs = cnvs.gr.final, variantsForEachCNV = variantsForEachCNV, filterParameters = filterParameters)
  class(resultsObject) <- "CNVfilteR_results"

  return(resultsObject)
}
