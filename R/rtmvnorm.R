


#' Random truncated multivaraite normal generator
#'
#' @param m Mean vector
#' @param V Variance-covariance matrix
#' @param w current values
#' @importFrom stats rnorm
#' @details This generates a random multivariate truncated normal random variable using 
#' using the approach of Li and Ghosh (2015). The first element is unbounded and the remaining 
#' elements are bounded below by 0.
#' @references Li Y. and Ghosh S.K. (2015) "Efficient Sampling Methods for Truncated Multivariate Normal and Student-t Distributions Subject to Linear Inequality Constraints." Journal of Statistical Theory and Practice. 9(4) 712-732.
#'
#' @export


rtmvnorm <- function(m,V,w){
  return(drop(rtmvnormc(m=m, V=V, w=w, n=length(m))))
}

