#' Compute p-value using the exact solution
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients.
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#'
#' @return \code{jaccard.test.exact} returns a list consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @importFrom stats rbinom pchisq rnorm runif
#' @importFrom combinat dmnom xsimplex nsimplex
#' @import magrittr
#' @export jaccard.test.exact
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test.exact(x,y)
jaccard.test.exact <- function(x, y, px = NULL, py = NULL, verbose = TRUE) {
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
  
    
  if(null.p) {
    pvalue = jaccard_mca_rcpp_known_p(px,py,m,j.obs,0)
  } else {
    pvalue = jaccard_mca_rcpp(px,py,m,j.obs,0)
  }

  return(
    list(
      statistics = j.obs,
      pvalue = pvalue,
      expectation = expectation
    )
  )
}