#' Compute p-value using the Measure Concentration Algorithm
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients.
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#' @param accuracy an error bound on approximating a multinomial distribution
#' @param error.type an error type on approximating a multinomial distribution ("average", "upper", "lower")
#'
#' @return \code{jaccard.test.mca} returns a list consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @importFrom stats rbinom pchisq rnorm runif
#' @importFrom IsoSpecR IsoSpecify
#' @export jaccard.test.mca
#' @useDynLib jaccard
#' 
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test.mca(x,y)
jaccard.test.mca <- function(x, y, px = NULL, py = NULL, accuracy = 1e-05, error.type = "average", verbose = TRUE) {
  if (length(x) != length(y)) stop("Length mismatch")
  m <- length(x)
  null.p<-FALSE
  if (is.null(px) | is.null(py)) {
    px <- mean(x)
    py <- mean(y)
    null.p <- TRUE
  }
   
  x <- as.logical(x)
  y <- as.logical(y)
  expectation <- (px * py)/(px + py - px * py)
  j.obs <- sum(x & y)/sum(x | y) - expectation
  
    
  if (null.p) {
    pvalue = jaccard_mca_rcpp_known_p(px,py,m,j.obs,accuracy)
  } else {
    pvalue = jaccard_mca_rcpp(px,py,m,j.obs,accuracy)
  }
  epsilon = 1- accuracy  
  pvalue <- switch(error.type, lower = pvalue, average = pvalue/epsilon, 
                   upper = pvalue + 1 - epsilon)

  return(
    list(
      statistics = j.obs,
      pvalue = pvalue,
      expectation = expectation,
      accuracy = accuracy,
      error.type = error.type
      )
  )
}
