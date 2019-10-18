#' getVariantScore
#'
#' @description
#' Returns score for a given allele frequency
#'
#' @details
#' Returns a value between -1 and 1. If the allele frequency increases the
#' evidence of discarding a CNV, then the score is positive. If the allele
#' frequency decreases the evidence for discarding a CNV, the score is negative.
#'
#' The model is based on the fuzzy logic and the score is calculated
#' using sigmoids. See the vignette to get more details.
#'
#' @param freq Variant allele frequency
#' @param expected.ht.mean Expected heterozygous SNV/indel allele frequency
#' @param expected.dup.ht.mean1 Expected heterozygous SNV/indel allele frequency when the variant IS NOT in the same allele than the CNV duplication call
#' @param expected.dup.ht.mean2 Expected heterozygous SNV/indel allele frequency when the variant IS in the same allele than the CNV duplication call
#' @param sigmoid.c1 Sigmoid c1 parameter
#' @param sigmoid.c2.vector Vector containing sigmoid c2 parameters for the six sigmoids functions
#' @param sigmoid.int1 Sigmoid int 1
#' @param sigmoid.int2 Sigmoid int 2
#'
#' @return variant score in the [-1, 1] range
#'
#' @importFrom pracma sigmoid
#'
getVariantScore <- function(freq, expected.ht.mean, expected.dup.ht.mean1, expected.dup.ht.mean2,
                            sigmoid.c1, sigmoid.c2.vector, sigmoid.int1, sigmoid.int2) {

  if (freq <= expected.dup.ht.mean1){
    c2 <- sigmoid.c2.vector[1]
    score <- - pracma::sigmoid(freq , a = sigmoid.c1, b = c2)  # extreme left sigmoid gives evidence of dup suspicius ht > reduces score
  } else if (freq <= sigmoid.int1) {
    c2 <- sigmoid.c2.vector[2]
    score <- - pracma::sigmoid(freq + (c2 - freq)*2, a = sigmoid.c1, b = c2)  # left sigmoid gives evidence of dup suspicius ht > reduces score
  } else if (freq <= expected.ht.mean){
    c2 <- sigmoid.c2.vector[3]
    score <- pracma::sigmoid(freq , a = sigmoid.c1, b = c2)  # central left sigmoid removes evidence of dup suspicius ht > increases score
  } else if (freq <= sigmoid.int2){
    c2 <- sigmoid.c2.vector[4]
    score <- pracma::sigmoid(freq + (c2 - freq)*2 , a = sigmoid.c1, b = c2)  # central right sigmoid removes evidence of dup suspicius ht > increases score
  } else if (freq <= expected.dup.ht.mean2){
    c2 <- sigmoid.c2.vector[5]
    score <- - pracma::sigmoid(freq , a = sigmoid.c1, b = c2)  # right sigmoid gives evidence of dup suspicius ht > reduces score
  } else {
    c2 <- sigmoid.c2.vector[6]
    score <- - pracma::sigmoid(freq + (c2 - freq)*2 , a = sigmoid.c1, b = c2)  # extreme right sigmoid gives evidence of dup suspicius ht > reduces score
  }

  return(score)
}
