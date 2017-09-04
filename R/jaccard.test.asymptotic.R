#' Compute p-value using an asymptotic approximation
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients.
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#'
#' @return \code{jaccard.test.asymptotic} returns a list consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @importFrom stats rbinom pchisq rnorm runif
#' @export jaccard.test.asymptotic
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test.asymptotic(x,y)
jaccard.test.asymptotic <- function(x, y, px = NULL, py = NULL, verbose = TRUE) {
  # length of fingerprints
  if(length(x) != length(y)) stop("Length mismatch")
  m <- length(x)

  # probabilities of ones
  if(is.null(px) | is.null(py)){
    px <- mean(x)
    py <- mean(y)
  }

  #converting x,y to bool
  x <- as.logical(x)
  y <- as.logical(y)

  expectation <- (px*py)/(px+py-px*py)
  j <- sum(x&y)/sum(x|y) - expectation

  q <- c(px*py,px+py-2*px*py)
  qq <- q[1]+q[2]
  sigma <- q[1]*q[2]*(1-q[2])/(qq^3)
  norm <- sqrt(m)*(j)/sqrt(sigma)
  pvalue <- pchisq(norm*norm,1,lower.tail = F)

  return(
    list(
      statistics = j,
      pvalue = pvalue,
      expectation = expectation
    )
  )
}