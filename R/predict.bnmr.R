

#' Predict method for class bnmr
#'
#' @param object An onject of class bnmr.
#' @param newdata A data.frame with newdata.
#' @param interval Interval type: none, confidence, or prediction.
#' @param deriv Derivative of function to be predicted.
#' @param totalmass The total mass to be destributed under the first derivative. This is only allowed when deriv=1.
#' @param ... Not used.
#'
#' @export

predict.bnmr <- function(object, 
                         newdata, 
                         interval = c("confidence", "prediction", "none"),
                         deriv=0, 
                         totalmass, 
                         ...){
  
 interval <- match.arg(interval)

 # make deriv feasible
 deriv <- min(max(0,round(deriv)),object$specs$M-1)
 
 if(deriv>0 & interval=="prediction"){
   message("Warning: prediction intervals are only available with deriv=0.")
 }
 
 
  
 if(missing(newdata)){
   
   if(deriv==0){# old data with no derivative
     if(interval=="confidence"){
       
       if(!missing(totalmass)) message("totalmass ignored with deriv=0.")
       return(object$fitted)
       
     }else if(interval=="prediction"){
       
       if(!missing(totalmass)) message("totalmass ignored with deriv=0.")
       
       f <- object$basis%*%object$gamma[,-c(1:round(object$specs$nburn/object$specs$nthin))]
       e <- scale(matrix(rnorm(prod(dim(f))),nrow(f),ncol(f)), 
                  center=FALSE, 
                  scale=1/object$parms$sigma[-c(1:round(object$specs$nburn/object$specs$nthin))])
       yhat <- f+e
       return(
         data.frame(
           fitted=object$fitted$fitted,
           lower=apply(yhat,1,quantile,0.025),
           upper=apply(yhat,1,quantile,0.975),
           x=object$fitted$x))
       
     }else{
       
       if(!missing(totalmass)) message("totalmass ignored with deriv=0.")
       
       return(object$fitted[,c("fitted","x")])
       
     }
   }else if(deriv==1){# old data with first derivative
     if(interval=="confidence"){
       
       if(missing(totalmass)){
         
         return(object$derivative)
         
       }else{
         
         x <- object$x
         xs <- (x-min(object$x))/diff(range(object$x))
         
         scl <- (object$gamma[object$specs$M+1,-c(1:round(object$specs$nburn/object$specs$nthin))]-object$gamma[1,-c(1:round(object$specs$nburn/object$specs$nthin))])/totalmass
         
         B <- buildBPbasis(xs, M=object$specs$M, deriv=1)
         f <- (B$B %*% B$A) %*% object$gamma[,-c(1:round(object$specs$nburn/object$specs$nthin))]
         f <- scale(f,center=FALSE,scale=scl)
         
         out <- data.frame(fitted=rowMeans(f),
                           lower=NA,
                           upper=NA,
                           x=x)
         
         out$lower=apply(f,1,quantile,0.025, na.rm=TRUE)
         out$upper=apply(f,1,quantile,0.975, na.rm=TRUE)
         
         return(out)
       }
       
     }else{
       
       return(object$derivative[,c("fitted","x")])
       
     }
   }else{# old data with >first derivative
     
     if(!missing(totalmass)) message("totalmass ignored with deriv>1.")
     
     x <- object$x
     xs <- (x-min(object$x))/diff(range(object$x))
     
     
     
     B <- buildBPbasis(xs, M=object$specs$M, deriv=deriv)
     f <- (B$B %*% B$A) %*% object$gamma[,-c(1:round(object$specs$nburn/object$specs$nthin))]
     

     out <- data.frame(fitted=rowMeans(f),
                       lower=NA,
                       upper=NA,
                       x=x)
     
     if(interval=="confidence"){
       
       out$lower=apply(f,1,quantile,0.025, na.rm=TRUE)
       out$upper=apply(f,1,quantile,0.975, na.rm=TRUE)
       
       return(out)
       
     }else{
       
       return(object$derivative[,c("fitted","x")])
       
     }
   }
   
 }else{# new data
   
  
   
   # fit with new data
   x <- newdata
   xs <- (x-min(object$x))/diff(range(object$x))
   inrange <- which(x>min(object$x) & x<max(object$x))
   B <- buildBPbasis(xs[inrange], M=object$specs$M, deriv=deriv)
   f <- (B$B %*% B$A) %*% object$gamma[,-c(1:round(object$specs$nburn/object$specs$nthin))]
   
   if(!missing(totalmass)){
     # scale
     scl <- (object$gamma[object$specs$M+1,-c(1:round(object$specs$nburn/object$specs$nthin))]-object$gamma[1,-c(1:round(object$specs$nburn/object$specs$nthin))])/totalmass
     f <- scale(f,center=FALSE,scale=scl)
   }
   
   
   out <- data.frame(fitted=NA,
                     lower=NA,
                     upper=NA,
                     x=newdata)
   out$fitted[inrange] <- rowMeans(f)
   
   if(interval=="confidence"){
     
     out$lower[inrange]=apply(f,1,quantile,0.025, na.rm=TRUE)
     out$upper[inrange]=apply(f,1,quantile,0.975, na.rm=TRUE)
     
     return(out)
     
   }else if(interval=="prediction"){
     
     e <- scale(matrix(rnorm(prod(dim(f))),nrow(f),ncol(f)), 
             center=FALSE, 
             scale=1/object$parms$sigma[-c(1:round(object$specs$nburn/object$specs$nthin))])
     yhat <- f+e 
     
     out$lower=apply(yhat,1,quantile,0.025, na.rm=TRUE)
     out$upper=apply(yhat,1,quantile,0.975, na.rm=TRUE)
     
     return(out)
     
   }else{
     
     return(out[,c("fitted","x")])
     
   }
 }
  
}



