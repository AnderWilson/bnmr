

#' Random univariate truncated normal generator
#'
#' @param l lower bound
#' @param u upper bound
#' @return A single draw from a univariate standard normal truncated between l and u.
#' @details The truncated standard univaraite random variable is generated with the mixed sampler of Li and Ghosh (2015).
#' @references Li Y. and Ghosh S.K. (2015) "Efficient Sampling Methods for Truncated Multivariate Normal and Student-t Distributions Subject to Linear Inequality Constraints." Journal of Statistical Theory and Practice. 9(4) 712-732.
#' @examples 
#' rtnorm(l=0)  # bounded below at 0
#' rtnorm(u=0)  # bounded above at 0
#' rtnorm(l=0,u=1)  # bounded between 0 and 1
#'
#' @export


rtnorm <- function(l=-Inf,u=Inf){
  return(rtnormc(l,u))
}


