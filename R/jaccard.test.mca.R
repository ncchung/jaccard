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
jaccard.test.mca <- function(x, y, px = NULL, py = NULL, verbose = TRUE, accuracy = 1e-05, error.type = "average") {
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

  # Defining auxilary variables used by IsoSpec
  inputData <- data.frame(
    element = c('C','C','C'),
    isotope = c('N1','N2','N3'),
    mass  	= c(1.0, 1.0,1.0),
    abundance = c(px*py, px+py-2*px*py, 1+px*py-px-py)
  )
  result <- IsoSpecify( molecule=c(C=m), stopCondition=1-accuracy, isotopes=inputData,showCounts=T, trim=T )
  result <- as.data.frame(result)
  #computing jaccard index for sample
  j.sample <- data.frame(prob=exp(result$logProb), 
                         index = result$N1/(result$N2+result$N1+((result$N2+result$N1)==0))+((result$N2+result$N1)==0)*expectation,
                         px = ((result$N2+result$N1)==0)*px+((result$N2+result$N1)>=0)*pmin(result$N1+result$N2,(2*result$N1+result$N2)*py/(px+py))/m,
                         py = ((result$N2+result$N1)==0)*py+((result$N2+result$N1)>=0)*pmax(result$N1,(2*result$N1+result$N2)*px/(px+py))/m
              )
  epsilon <- kahanSum(j.sample$prob)

  # centering jaccard index
  if(null.p){
    j.sample$index <- j.sample$index - j.sample$px*j.sample$py/(j.sample$px+j.sample$py-j.sample$px*j.sample$py)
  } else {
    j.sample$index <- j.sample$index - expectation
  }

  #computing p-value
  ind  <- which(abs(j.sample$index)>= abs(j.obs))
  pvalue <- kahanSum(j.sample$prob[ind])
  pvalue <- switch (error.type,
                    "lower" = pvalue,
                    "average" = pvalue/epsilon,
                    "upper" = pvalue+1-epsilon
                    )

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
