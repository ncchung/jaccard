#' Compute p-value using the EC-BLAST method
#'
#' In the EC-BLAST paper, Rahman et al. (2014) provide the following description:
#' The mean (μ) and s.d. (σ) of the similarity scores are used to
#' define the z score, z = (Tw – μ)/σ. For the purpose of calculating
#' the P value, only hits with T > 0 are considered. The P value w
#' is derived from the z score using an extreme value distribution 
#' P = 1 – exp(−e−zπ/√(6) − Γ′ (1)), where the Euler-Mascheroni constant Γ′ (1) ≈ 0.577215665.
#'
#' @param j a numeric vector of observed Jaccard coefficients (uncentered)
#' @return \code{jaccard.rahman} returns a numeric vector of p-values
#'
#' @references Rahman, Cuesta, Furnham, Holliday, and Thornton (2014) EC-BLAST: a tool to automatically search and compare enzyme reactions. Nature Methods, 11(2) \url{http://www.nature.com/nmeth/journal/v11/n2/full/nmeth.2803.html}
#'
#' @export jaccard.rahman
jaccard.rahman <- function(j) {
	mu = mean(j)
	sd = sd(j)
	z = (j - mu) / sd
	p = 1 - exp(-exp(-z*pi/sqrt(6)-0.577215665))
	return(p)
}