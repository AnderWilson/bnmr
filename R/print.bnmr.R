

#' Print method for class bnmr
#'
#' @param x An onject of class bnmr.
#' @param ... Not used.
#'
#' @export

print.bnmr <- function(x, ...){
  
  cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Model:\n",x$modelstatement,"\n", sep = "")
  
}



