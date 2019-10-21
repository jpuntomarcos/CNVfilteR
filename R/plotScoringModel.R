#' plotVariantsForCNV
#'
#' @description
#' Plots scoring model used for CNV duplications
#'
#' @param expected.ht.mean Expected heterozygous SNV/indel allele frequency
#' @param expected.dup.ht.mean1 Expected heterozygous SNV/indel allele frequency when the variant IS NOT in the same allele than the CNV duplication call
#' @param expected.dup.ht.mean2 Expected heterozygous SNV/indel allele frequency when the variant IS in the same allele than the CNV duplication call
#' @param sigmoid.c1 Sigmoid c1 parameter
#' @param sigmoid.c2.vector Vector containing sigmoid c2 parameters for the six sigmoids functions
#'
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
#' # Plot scoring model for duplication CNVs
#' p <- results$filterParameters
#' plotScoringModel(expected.ht.mean = p$expected.ht.mean, expected.dup.ht.mean1 = p$expected.dup.ht.mean1,
#'                   expected.dup.ht.mean2 = p$expected.dup.ht.mean2, sigmoid.c1 = p$sigmoid.c1, sigmoid.c2.vector = p$sigmoid.c2.vector)
#'
#' @return
#' nothing
#'
#'
#' @importFrom graphics legend plot lines points text abline
#' @import assertthat
#' @export plotScoringModel
#'
#'
plotScoringModel <- function(expected.ht.mean, expected.dup.ht.mean1,
                             expected.dup.ht.mean2, sigmoid.c1,
                             sigmoid.c2.vector) {


  # Check input
  assertthat::assert_that(assertthat::is.number(expected.ht.mean))
  assertthat::assert_that(assertthat::is.number(expected.dup.ht.mean1))
  assertthat::assert_that(assertthat::is.number(expected.dup.ht.mean2))
  assertthat::assert_that(assertthat::is.number(sigmoid.c1))
  assertthat::assert_that(is.numeric(sigmoid.c2.vector) && length(sigmoid.c2.vector) == 6)

  # intersections between two sigmoids curves
  sigmoid.int1 <- mean(c(expected.dup.ht.mean1, expected.ht.mean))
  sigmoid.int2 <- mean(c(expected.dup.ht.mean2, expected.ht.mean))


  # plot first part of model, (score <= 0, those that confirm CNV duplication)
  x <- seq(20, sigmoid.int1, length.out = 1000)
  y <- lapply(x, getVariantScore,  expected.ht.mean, expected.dup.ht.mean1, expected.dup.ht.mean2,
              sigmoid.c1, sigmoid.c2.vector, sigmoid.int1, sigmoid.int2)

  graphics::plot(x, y, type = "l", col = CONFIRM_COLOR,
                 xlim=c(20, 80), ylim=c(-1,1),  bty="L", lwd=2,
                 xlab = "variant allele frequency", ylab = "score",
                 main = "Scoring model for duplication CNVs",
                 panel.first = c(graphics::abline(h = 0, col = "gray"),
                                 graphics::abline(h = c(1, 0.5, -0.5, -1), lty = 2, col = "lightgray")))

  # plot second part of model, (score >= 0, those that discard CNV duplication)
  x <- seq(sigmoid.int1 + 0.00001, sigmoid.int2, length.out = 1000)
  y <- lapply(x, getVariantScore,  expected.ht.mean, expected.dup.ht.mean1, expected.dup.ht.mean2,
              sigmoid.c1, sigmoid.c2.vector, sigmoid.int1, sigmoid.int2)
  graphics::lines(x, y, type = "l", col = DISCARD_COLOR, xlab = "", ylab = "",
                  lwd=2)

  # plot third part of model, (score <= 0, those that confirm CNV duplication)
  x <- seq(sigmoid.int2 + 0.00001, 80, length.out = 1000)
  y <- lapply(x, getVariantScore,  expected.ht.mean, expected.dup.ht.mean1, expected.dup.ht.mean2,
              sigmoid.c1, sigmoid.c2.vector, sigmoid.int1, sigmoid.int2)
  graphics::lines(x, y, type = "l", col = CONFIRM_COLOR, xlab = "", ylab = "",
                  lwd=2)


  # plot expected.dup.ht.mean1-expected.ht.mean-expected.dup.ht.mean2 marks
  graphics::points((c(expected.dup.ht.mean1, expected.dup.ht.mean2)), c(-1,-1), pch=20, cex=0.7)
  graphics::text(c(expected.dup.ht.mean1, expected.dup.ht.mean2), c(-1.04,-1.04), offset = 0,
       labels=c("expected.dup.ht.mean1", "expected.dup.ht.mean2"),
       cex= 0.7, col = "#666666")
  graphics::points((c(expected.ht.mean)), c(1), pch=20, cex=0.7)
  graphics::text(c(expected.ht.mean), c(0.97), xpd=TRUE,
         labels=c("expected.ht.mean"),
         cex= 0.7, pos=3, col = "#666666")

  # plot legend
  graphics::legend("topright",
                   legend=c("Evidence for discarding a dup. CNV",
                            "Evidence for confirming a dup. CNV"),
                   border = "white", bty = "o",
                   fill = c(DISCARD_COLOR, CONFIRM_COLOR),
                   cex = 0.8, bg = "white",  box.col = "gray")
}

