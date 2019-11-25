

#' Summary method for class bnmr
#'
#' @param object An onject of class bnmr.
#' @param ... Not used.
#' 
#' @importFrom stats sd quantile
#'
#' @export

summary.bnmr <- function(object, ...){
  
  
  parameters <- data.frame(
    mean=colMeans(object$parms[round(object$specs$nburn/object$specs$nthin+1):nrow(object$parms),]),
    sd = apply(object$parms[round(object$specs$nburn/object$specs$nthin+1):nrow(object$parms),],2,sd),
    q2.5 = apply(object$parms[round(object$specs$nburn/object$specs$nthin+1):nrow(object$parms),],2,quantile, 0.025),
    q97.5 = apply(object$parms[round(object$specs$nburn/object$specs$nthin+1):nrow(object$parms),],2,quantile, 0.025)
  )
  row.names(parameters) <- colnames(object$parms)
  
  pr_flat = mean(rowMeans(object$beta[round(object$specs$nburn/object$specs$nthin+1):nrow(object$beta),-1])==0)
  pr_linear = mean(apply(object$beta[round(object$specs$nburn/object$specs$nthin+1):nrow(object$beta),-1],1,function(a) length(unique(a)))==1)
    
  out <- list(parameters=parameters, pr_flat=pr_flat, pr_linear=pr_linear, call=object$call, modelstatement=object$modelstatement, specs=object$specs)
  class(out) <- "summary.bnmr"
  return(out)
  
}



