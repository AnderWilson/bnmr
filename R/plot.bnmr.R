

#' Plot method for class bnmr
#'
#' @param x An onject of class bnmr.
#' @param ... Arguments passed to plot function.
#'
#' @importFrom graphics lines plot
#' @export

plot.bnmr <- function(x, ...){
  
  plot(y~x, data=data.frame(x=x$x, y=x$y), las=1, ...)
  lines(fitted~x, data=x$fitted )
  lines(lower~x, data=x$fitted, lty=2 )
  lines(upper~x, data=x$fitted, lty=2 )
  
}

