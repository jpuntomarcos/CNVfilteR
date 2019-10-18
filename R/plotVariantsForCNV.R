#' plotVariantsForCNV
#'
#' @description
#' Plots a CNV with all the variants matching it
#'
#' @param cnvfilter.results S3 object returned by \code{filterCNVs} function
#' @param cnv.id CNV id for which to plot variants
#' @param points.cex Points cex (size). (Defaults to 1)
#' @param points.pch Points pch (symbol). (Defaults to 19)
#' @param legend.x.pos Legend x position. (Defaults to 0.08)
#' @param legend.y.pos Legend y position. (Defaults to 0.25)
#' @param cnv.zoom.margin If TRUE, the zoom leaves an small margin at both sides of the CNV. False otherwise. (Defaults to TRUE)
#'
#' @return invisibly returns a \code{karyoplot} object
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
#' # Check CNVs that can be filtered
#' as.data.frame(results$cnvs[results$cnvs$filter == TRUE])
#'
#' # Plot one of them
#' plotVariantsForCNV(results, "3")
#'
#'
#' @importFrom CopyNumberPlots plotCopyNumberCalls
#' @import karyoploteR
#' @importFrom graphics legend
#' @importFrom GenomicRanges width start end
#' @importFrom methods is
#' @import assertthat
#' @export plotVariantsForCNV
#'
plotVariantsForCNV <- function(cnvfilter.results, cnv.id, points.cex = 1,
                               points.pch = 19,
                               legend.x.pos = 0.08, legend.y.pos = 0.25,
                               cnv.zoom.margin = TRUE) {

  # Check input
  assertthat::assert_that(methods::is(cnvfilter.results, "CNVfilteR_results"))
  assertthat::assert_that(assertthat::is.string(cnv.id))
  assertthat::assert_that(assertthat::is.number(points.cex))
  assertthat::assert_that(assertthat::is.number(points.pch))
  assertthat::assert_that(is.logical(cnv.zoom.margin))

  # Get filter params values for later use
  params <- cnvfilter.results$filterParameters

  # Load variants for desired cnv
  vars <- cnvfilter.results$variantsForEachCNV[[cnv.id]]

  # Move freq values to the 0-1 interval
  vars$alt.freq <- vars$alt.freq / 100.0

  # Add cn column (required by CopyNumberPlots)
  cnvs.gr <- auxAddCNcolumn(cnvfilter.results$cnvs)
  cnv <- cnvs.gr[cnvs.gr$cnv.id == cnv.id,]

  # Build main title
  title <-  paste0(cnv$cnv, " at ", toString(cnv))
  if (cnv$filter == "TRUE"){
    if (cnv$cnv == "deletion") {
      title <- paste0(title, ", it can be filtered")
    } else if (cnv$cnv == "duplication") {
      title <- paste0(title, ", it can be filtered with a score of ", round(as.numeric(cnv$score), 4))
    }
  }

  # Prepare vars
  r0 = 0.1; r1=0.95;
  duphtmean1 <- params$expected.dup.ht.mean1 / 100.0
  duphtmean2 <- params$expected.dup.ht.mean2 / 100.0
  htmean <- params$expected.ht.mean / 100.0
  zoomGR <- cnv
  if (cnv.zoom.margin) {
    GenomicRanges::start(zoomGR) <- GenomicRanges::start(zoomGR) - GenomicRanges::width(cnv) * 0.1
    GenomicRanges::end(zoomGR) <- GenomicRanges::end(zoomGR) + GenomicRanges::width(cnv) * 0.1
  }

  ## Plot ##

  # Plot different elements
  pp <- karyoploteR::getDefaultPlotParams(4)
  pp$leftmargin <- 0.15
  pp$ideogramheight <- 2
  pp$topmargin <- 40
  kp <- karyoploteR::plotKaryotype(zoom = zoomGR, plot.type = 4,
                                   plot.params = pp, main = title)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = GenomicRanges::width(cnv) / 5.0,
                                cex = 0.7, digits = 5, add.units = TRUE)
  CopyNumberPlots::plotCopyNumberCalls(kp, cnv, cn.colors = CNV_COLORS, r0=0,
                                       r1=0.03, labels = "CNV")
  karyoploteR::kpAddLabels(kp, r0 = r0, r1 = r1, labels = c("variant allele frequency"),
                           srt = 90, pos = 3, cex = 0.8, label.margin = 0.1)

  # Add guide lines depending on CNV type
  if (cnv$cnv == "deletion") {
    karyoploteR::kpAxis(kp, ymin=0, ymax=1, r0=r0, r1=r1, col="gray50", cex=0.8,
                        labels = c("0", htmean, "1"), tick.pos = c(0, htmean, 1))
    karyoploteR::kpAbline(kp, h=c(htmean), col="gray80", ymin=0, ymax=1, r0=r0,
                          r1=r1, lty=2, lwd = 1)
  } else if (cnv$cnv == "duplication"){
    karyoploteR::kpAxis(kp, ymin=0, ymax=1, r0=r0, r1=r1, col="gray50", cex=0.8,
                        labels = c("0", duphtmean1, htmean, duphtmean2, "1"),
                        tick.pos = c(0, duphtmean1, htmean, duphtmean2, 1))
    karyoploteR::kpAbline(kp, h=c(duphtmean1, htmean, duphtmean2), col="gray80",
                          ymin=0, ymax=1, r0=r0, r1=r1, lty=2, lwd = 1)
  }


  # Draw legend
  graphics::legend(x=legend.x.pos, y=legend.y.pos, legend=c("Discard CNV", "Confirm CNV", "Neutral", "CNV deletion", "CNV duplication"),
        pch = c(points.pch, points.pch, points.pch, NA, NA), pt.cex=1, bg = "white",  box.col = "gray",
        col = c(DISCARD_COLOR, CONFIRM_COLOR, NEUTRAL_COLOR, NA, NA),
        border = "white",
        fill = c(NA, NA, NA, CNV_COLORS[2], CNV_COLORS[4]), ncol=2, cex = 0.8,
        bty="l")

  # Add points depending o CNV type
  if (cnv$cnv == "deletion") {
    v_sub <- vars[vars$type=="ht",]
    karyoploteR::kpPoints(kp, chr = as.character(v_sub$seqnames), x=v_sub$start, y=v_sub$alt.freq, col=DISCARD_COLOR, r0 = r0, r1=r1, cex=points.cex, pch=points.pch)
    v_sub <- vars[vars$type=="hm",]
    karyoploteR::kpPoints(kp, chr = as.character(v_sub$seqnames), x=v_sub$start, y=v_sub$alt.freq, col=NEUTRAL_COLOR, r0 = r0, r1=r1, cex=points.cex, pch=points.pch)
  } else if (cnv$cnv == "duplication") {
    v_sub <- vars[vars$score == 0,]
    karyoploteR::kpPoints(kp, chr = as.character(v_sub$seqnames), x=v_sub$start, y=v_sub$alt.freq, col=NEUTRAL_COLOR, r0 = r0, r1=r1, cex=points.cex, pch=points.pch)
    v_sub <- vars[vars$score > 0,]
    karyoploteR::kpPoints(kp, chr = as.character(v_sub$seqnames), x=v_sub$start, y=v_sub$alt.freq, col=DISCARD_COLOR, r0 = r0, r1=r1, cex=points.cex, pch=points.pch)
    v_sub <- vars[vars$score < 0,]
    karyoploteR::kpPoints(kp, chr = as.character(v_sub$seqnames), x=v_sub$start, y=v_sub$alt.freq, col=CONFIRM_COLOR, r0 = r0, r1=r1, cex=points.cex, pch=points.pch)
  }


  return(invisible(kp))
}

