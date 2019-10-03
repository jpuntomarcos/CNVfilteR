#' filterCNVs
#'
#' @description
#' Identifies those copy number calls that can be filtered
#'
#' @details
#' Checks all the variants (SNV and optionally INDELs) matching each CNV present in \code{cnvs.gr} to decide whether a CNV can be filtered or not.
#' It returs a S3 object with 3 elments: \code{cnvs}, \code{variantsForEachCNV} and \code{filterParameters}. See return section for further details.
#'
#' A CNV deletion can be filtered if there is at least \code{ht.deletions.threshold}% of heterozygous variants matching the CNV.
#' A CNV duplication can be filtered if the \code{score} is >= \code{dup.threshold.score} after computing all heterozygous variants matching the CNV.
#'
#' If a CNV can be filtered, then the value TRUE is set in the \code{filter} column of the \code{cnvs} element.
#'
#' @param cnvs.gr \code{GRanges} containing CNVs to be filtered. Use \code{loadCNVcalls} to load them.
#' @param vcfs List of \code{GRanges} containing all variants (SNV/indel) obtaining with the \code{loadVCFs} function.
#' @param expected.ht.mean Expected heterozygous SNV/indel depth. (defaults to 50)
#' @param expected.dup.ht.mean1 Expected heterozygous SNV/indel depth when the variant IS NOT in the same allele than the CNV duplication call. (defaults to 33)
#' @param expected.dup.ht.mean2 Expected heterozygous SNV/indel depth when the variant IS in the same allele than the CNV duplication call. (defaults to 66)
#' @param sigmoid.c1 Sigmoid c1 parameter. (defaults to 2)
#' @param sigmoid.c2.vector Vector containing sigmoid c2 parameters for the six sigmoids functions. (defaults to c(28, 39, 44, 56, 60, 72))
#' @param dup.threshold.score Limit value to decide if a CNV duplication can be filtered or not. A CNV duplication can be filtered if the total score computed from heterozygous variants matching the CNV is equal or greater than \code{dup.threshold.score}.  (defaults to 0.5)
#' @param ht.deletions.threshold Minimum percentage of heterozygous variants matching a CNV deletion to filter that CNV. (defaults to 15)
#' @param verbose Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A S3 object with 3 elements:
#'  - \code{cnvs}: \code{GRanges} with the input CNVs and the meta-columns added during the call:
#'    - \code{filter}: Set to TRUE if the CNV can be filtered
#'    - \code{nVariants}: Number of variants matching the CNV
#'    - \code{nInFavor}: For a CNV duplication, number of variants in favor of filtering (those with a positive score)
#'    - \code{nAgainst}: For a CNV duplication, number of variants against filtering (those with a negative score)
#'    - \code{score}: total score when computing all the variants scores
#'    - \code{id}: CNV id
#'  - \code{variantsForEachCNV}: named list where each name correspond to a CNV id and the value is a \code{data.frame} with all variants matching that CNV
#'  - \code{filterParameters}: input parameters used for filtering
#'
#'
#' @examples
#' # Load CNVs data
#' cnvs.file <- system.file("extdata", "DECoN.CNVcalls.csv", package = "CNVfilteR", mustWork = TRUE)
#' cnvs.gr <- loadCNVcalls(path = cnvs.file, chr.column = "Chromosome", start.column = "Start", end.column = "End", cnv.column = "CNV.type", sample.column = "Sample")
#'
#' # Load VCFs data
#' vcf.paths <- c(system.file("extdata", "variants.sample1.vcf.gz", package = "CNVfilteR", mustWork = TRUE),
#'                system.file("extdata", "variants.sample2.vcf.gz", package = "CNVfilteR", mustWork = TRUE))
#' vcfs <- loadVCFs(vcf.paths, cnvs.gr = cnvs.gr)
#'
#' # Filter CNVs
#' results <- filterCNVs(cnvs.gr, vcfs)
#'
#' # Check CNVs that can be filtered
#' as.data.frame(results$cnvs[results$cnvs$filter == TRUE])
#'
#'
#' @import assertthat
#' @importFrom pracma sigmoid
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges mcols
#' @export filterCNVs
#'
filterCNVs <- function(cnvs.gr, vcfs, expected.ht.mean = 50, expected.dup.ht.mean1 = 33, expected.dup.ht.mean2 = 66,
                       sigmoid.c1 = 2, sigmoid.c2.vector = c(28, 39, 44, 56, 60, 72),
                       dup.threshold.score = 0.5, ht.deletions.threshold = 15, verbose = FALSE) {

  # Check input
  # assert_that(is.data.frame(cnvs.df))
  # if (identical(names(cnvs.df), CNV_DF_COLUMNS) | identical(cnvs.df, CNV_DF_COLUMNS_WITH_SAMPLE))
  #   stop("Expected columns for data.frame are 'sample'(optional) 'chr' 'start' 'end' 'cnv'")
  # assert_that(is.character(vcfs))


  # intersections between two sigmoids curves. For later use
  sigmoid.int1 <- mean(c(expected.dup.ht.mean1, expected.ht.mean))
  sigmoid.int2 <- mean(c(expected.dup.ht.mean2, expected.ht.mean))

  variantsForEachCNV <- list()

  # Filter CNVs depending on variant calls.
  cnvs.gr$filter <- ''
  cnvs.gr$nVariants <- ''
  cnvs.gr$nInFavor <- ''
  cnvs.gr$nAgainst <- ''
  cnvs.gr$score <- ''
  cnvs.gr$id <- ''
  nCNVs <- length(cnvs.gr)
  for (i in seq_len(nCNVs)){

    cnvId <- toString(i)
    cnvs.gr[i]$id <- cnvId

    if (cnvs.gr[i]$sample %in% names(vcfs)){
      variants <- vcfs[[cnvs.gr[i]$sample]]
      matchingVariants <- subsetByOverlaps(variants, cnvs.gr[i],  type = "any")

      # count number of ht and hm
      htMatchingVariants <- matchingVariants[matchingVariants$type == "ht",]
      nHT <- nrow(mcols(htMatchingVariants))
      nHM <- nrow(mcols(matchingVariants[matchingVariants$type == "hm",]))
      cnvs.gr[i]$nVariants <- length(matchingVariants)

      # if CNV is deletion & ht / (ht + hm)  > ht.deletions.threshold > discard CNV
      if (cnvs.gr[i]$cnv == "deletion" && nHT > 0 && (nHT / (nHT + nHM) > (ht.deletions.threshold / 100.0) )) {
        cnvs.gr[i]$filter <- TRUE
        if (verbose){
          message(paste0("Discarded CNV deletion at ", toString(cnvs.gr[i]) ," for sample ", cnvs.gr[i]$sample))
        }
        cnvs.gr[i]$nInFavor <- nHT
      }

      # Add score column for those variants mathching a duplication CNV
      if (cnvs.gr[i]$cnv == "duplication" & length(matchingVariants) > 0){
        matchingVariants$score <- 0
      }

      # if CNV is duplication  & exists ht
      if (cnvs.gr[i]$cnv == "duplication" & nHT > 0) {

        total.score <- 0
        in.favor <- 0
        against <- 0
        for (j in seq_len(length(htMatchingVariants))){

          v <- htMatchingVariants[j]
          if (v$alt.freq <= expected.dup.ht.mean1){
            against <- against + 1
            c2 <- sigmoid.c2.vector[1]
            score <- - sigmoid(v$alt.freq , a = sigmoid.c1, b = c2)  # extreme left sigmoid gives evidence of dup suspicius ht > reduces score
          } else if (v$alt.freq <= sigmoid.int1) {
            against <- against + 1
            c2 <- sigmoid.c2.vector[2]
            score <- - sigmoid(v$alt.freq + (c2 - v$alt.freq)*2, a = sigmoid.c1, b = c2)  # left sigmoid gives evidence of dup suspicius ht > reduces score
          } else if (v$alt.freq <= expected.ht.mean){
            in.favor <- in.favor + 1
            c2 <- sigmoid.c2.vector[3]
            score <- sigmoid(v$alt.freq , a = sigmoid.c1, b = c2)  # central left sigmoid removes evidence of dup suspicius ht > increases score
          } else if (v$alt.freq <= sigmoid.int2){
            in.favor <- in.favor + 1
            c2 <- sigmoid.c2.vector[4]
            score <- sigmoid(v$alt.freq + (c2 - v$alt.freq)*2 , a = sigmoid.c1, b = c2)  # central right sigmoid removes evidence of dup suspicius ht > increases score
          } else if (v$alt.freq <= expected.dup.ht.mean2){
            against <- against + 1
            c2 <- sigmoid.c2.vector[5]
            score <- - sigmoid(v$alt.freq , a = sigmoid.c1, b = c2)  # right sigmoid gives evidence of dup suspicius ht > reduces score
          } else {
            against <- against + 1
            c2 <- sigmoid.c2.vector[6]
            score <- - sigmoid(v$alt.freq + (c2 - v$alt.freq)*2 , a = sigmoid.c1, b = c2)  # extreme right sigmoid gives evidence of dup suspicius ht > reduces score
          }

          # compute total score and save score
          total.score <- total.score + score
          matchingVariants[names(v)]$score <- score

        }

        # Discard CNV depending on dup.threshold.score
        if (total.score >= dup.threshold.score){
          cnvs.gr[i]$filter <- TRUE
          if (verbose){
            message(paste0("Discarded CNV duplication at ", toString(cnvs.gr[i]) ," for sample ", cnvs.gr[i]$sample,  ", score: ", round(total.score, 4)))
            message(matchingVariants)
          }
        } else {
          if (verbose){
            message(paste0("CNV duplication NOT discarded at ", toString(cnvs.gr[i]) ," for sample ", cnvs.gr[i]$sample, ", score: ", round(total.score, 4)))
            message(matchingVariants)
          }
        }

        cnvs.gr[i]$score <- total.score
        cnvs.gr[i]$nAgainst <- against
        cnvs.gr[i]$nInFavor <- in.favor
      }


      if (length(matchingVariants) > 0)
        variantsForEachCNV[[cnvId]] <- as.data.frame(matchingVariants)

    }

  }

  # Calculate % of CNVs to be filtered
  nFiltered <- length(cnvs.gr[cnvs.gr$filter == TRUE,])
  pct <- round(nFiltered/nCNVs*100, 2)
  message(paste0(nFiltered, " of ", nCNVs, " (", pct ,"%) CNVs can be filtered"))

  # Calculate % with matching variants to be filtered
  nCNVsWithMatchingVariants <- length(cnvs.gr[cnvs.gr$nVariants > 0 ,])
  pct <- round(nFiltered/nCNVsWithMatchingVariants*100, 2)
  message(paste0(nFiltered, " of ", nCNVsWithMatchingVariants, " (", pct ,"%) CNVs with overlapping SNVs can be filtered"))

  # Create S3 object and return
  filterParameters <- list(expected.ht.mean = expected.ht.mean, expected.dup.ht.mean1 = expected.dup.ht.mean1,
                           expected.dup.ht.mean2 = expected.dup.ht.mean2, sigmoid.c1 = sigmoid.c1, sigmoid.c2.vector = sigmoid.c2.vector,
                           dup.threshold.score = dup.threshold.score, ht.deletions.threshold = ht.deletions.threshold)
  resultsObject <- list(cnvs = cnvs.gr, variantsForEachCNV = variantsForEachCNV, filterParameters = filterParameters)
  class(resultsObject) <- "CNVfilteR_results"

  return(resultsObject)
}

