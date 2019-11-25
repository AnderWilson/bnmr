

#' Print method for class summary.bnmr
#'
#' @param x An onject of class summary.bnmr.
#' @param ... Not used.
#'
#' @export

print.summary.bnmr <- function(x, ...){
  
  cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Model:\n",x$modelstatement,"\n", sep = "")
  
  cat("Parameters:\n")
  print(round(x$parameters,3))
  cat("\n", sep = "")
  
  cat("Posterior summary:\n",
    "  Pr(  flat  | - ) = ",round(x$pr_flat,3),
  "\n  Pr( linear | - ) = ",round(x$pr_linear,3),"\n\n"
  ,sep="")
  
  cat("MCMC Specifications:\n",
      "  Iterations: ",x$specs$niter,"\n",
      "  Burn-in   : ",x$specs$nburn,"\n",
      "  Thinning  : ",x$specs$nthin,"\n", sep="")
}



