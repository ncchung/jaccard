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
  # length of fingerprints
  if(length(x) != length(y)) stop("Length mismatch")
  m <- length(x)
  null.p <- FALSE
  # probabilities of ones
  if(is.null(px) | is.null(py)){
    px <- max(mean(x),mean(y))
    py <- min(mean(x),mean(y))
    null.p <- TRUE
  }

  #converting x,y to bool
  x <- as.logical(x)
  y <- as.logical(y)

  expectation <- (px*py)/(px+py-px*py)
  j.obs <- sum(x&y)/sum(x|y) - expectation

  if(verbose == TRUE & nsimplex(3,m)>10^5) print("The exact version  could take a long time")
  result<- t(xsimplex(3,m)) %>% as.data.frame
  names(result)=c('N1','N2','N3')
  result$prob <- apply(result,1,dmnom,prob=c(px*py,px+py-2*px*py,(1-px)*(1-py)))

  # computing Jaccard for sample
  j.sample <- data.frame(prob=result$prob, 
                         index = result$N1/(result$N2+result$N1+((result$N2+result$N1)==0))+((result$N2+result$N1)==0)*expectation,
                         px = ((result$N2+result$N1)==0)*px+((result$N2+result$N1)>=0)*pmin(result$N1+result$N2,(2*result$N1+result$N2)*py/(px+py))/m,
                         py = ((result$N2+result$N1)==0)*py+((result$N2+result$N1)>=0)*pmax(result$N1,(2*result$N1+result$N2)*px/(px+py))/m
                  )
  epsilon <- kahanSum(j.sample$prob)

  # centering
  if(null.p){
    j.sample$index <- j.sample$index - j.sample$px*j.sample$py/(j.sample$px+j.sample$py-j.sample$px*j.sample$py)
  } else {
    j.sample$index <- j.sample$index - expectation
  }

  #computing p-value
  ind  <- which(abs(j.sample$index)>= abs(j.obs))
  pvalue <- kahanSum(j.sample$prob[ind])

  return(
    list(
      statistics = j.obs,
      pvalue = pvalue,
      expectation = expectation
    )
  )
}