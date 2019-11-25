#' Build Bernstein Polynomial Basis, Derivatives, and Transformation functions
#'
#' @param t time or the predictor
#' @param M the number of basis functions
#' @param deriv The derivative order (defalt is 0)
#' @importFrom stats dbinom
#'
#' @export

buildBPbasis <- function(t,M, deriv=0){

  numbasis <- M+1-deriv
  
  inbounds <- which(t>=0 & t<=1)
  
  B<- matrix(NA,length(t),numbasis)
  for(k in 1:numbasis){
    B[inbounds,k]<- dbinom(k-1,numbasis-1,t[inbounds])
  }
  
  A <- diag(M+1)
  if(deriv>0){

    for(i in M:(M-deriv+1)){
      B <- B*M

      A0<- diag(-1,i,i+1)
      for (j in 1:i){
        A0[j,j+1] = 1
      }
      A <- A0 %*% A
    }

  }

  return(list(B=B,A=A))
}
