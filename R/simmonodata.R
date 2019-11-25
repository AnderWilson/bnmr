
#' Simulate monotone data
#'
#' @param n the number of observations
#' @param scenario Four scenaries y=f(x) + e to generate data from, where f(x) is flat, linear, flatnon, piecewise, and wavy
#' @param errorsd The error variance
#' @return dataframe with columns for x, f(x), and y.
#' @importFrom stats rnorm runif
#' @examples A <- simmonodata(n=100,scenario="flat")
#'
#' @export

simmonodata<-function(n=100,scenario="flat", errorsd=0.25){

  x<-runif(n,0,1)
  x <- x[order(x)]
  grideval <- data.frame(f=rep(0,100),fprime=rep(0,100),x=seq(min(x),max(x),length=100))
  
  if(scenario=="flat"){
    
    f <- rep(1,n)
    fprime <- rep(0,n)
    
    grideval$f <- rep(1,100)
    grideval$fprime <- rep(0,100)
    
  }else  if(scenario=="linear"){
    
    f <- x
    fprime <- rep(1,n)
    
    grideval$f <- grideval$x
    grideval$fprime <- rep(1,100)
    
  }else if(scenario=="flatnon"){
    
    fprime <-f <- rep(0,n)
    f[x>0.5] <- (2*(x[x>0.5]-0.5) )^2
    fprime[x>0.5] <- 4*(2*(x[x>0.5]-0.5))
    
    grideval$f[grideval$x>0.5] <- (2*(grideval$x[grideval$x>0.5]-0.5) )^2
    grideval$fprime[grideval$x>0.5] <- 4*(2*(grideval$x[grideval$x>0.5]-0.5))
    
  }else if(scenario=="piecewise"){
    
    fprime <-f <- rep(0,n)
    f[x>.4] <- pmin(2*x[x>.4]-.8,.6)
    fprime[x>.4 & x<.7] <- 2
    
    
    grideval$f[grideval$x>.4] <- pmin(2*grideval$x[grideval$x>.4]-.8,.6)
    grideval$fprime[grideval$x>.4 & grideval$x<.7] <- 2
    
  }else if(scenario=="wavy"){
    
    f <- sin(3*pi*x)/(3*pi)+x
    fprime <- cos(3*pi*x) + 1
    
    grideval$f <- sin(3*pi*grideval$x)/(3*pi)+grideval$x
    grideval$fprime <- cos(3*pi*grideval$x) + 1
    
  }

  y <- f + rnorm(n,0,errorsd)

  return(list(x=x,f=f,y=y, fprime=fprime, grideval=grideval))
}


