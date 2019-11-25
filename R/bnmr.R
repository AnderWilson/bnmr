

#' Bayesin Nonparametric Monotone Regression
#'
#' @param y vector of outcomes
#' @param x vector of predictors
#' @param M order of the Bernstein polynomial basis to be used.
#' @param niter number of iterations including burnin.
#' @param nburn number of iteractions discarded as burnin.
#' @param nthin thining number for posterior sample.
#' @param priors a list of priors. The list may contain elements alpha_a, alpha_b, pi_a, pi_b, mu0, and phi0 all scalars. The prios assigned will be alpha~Gamma(alpha_a,alpha_b), pi~beta(pi_a,pi_b), mu~N(mu0,phi0^2), theta_0~Normal(0,intsd^2).  The defaults are: alpha_a=1, alpha_b=1, pi_a=1, pi_b=1, mu0=0.5, phi0=0.25, intsd=2.
#' @importFrom stats rbeta var runif rgamma qnorm dnorm rnorm pnorm quantile
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export



# y <- dat$y; x <- dat$x; M=20;niter=500;nburn= round(niter/2); nthin=1; set.seed(1234)
# alpha_a <- 2; alpha_b <- 1; pi_a <- 10; pi_b <- 1


bnmr <- function(y,x,M,niter,nburn=round(niter/2), nthin=1, priors){
  
  
  n <- length(y)
  ys <- (y-mean(y))/sd(y)
  xs <- (x-min(x))/diff(range(x))
  
  # default order of polynomial
  if(missing(M)) M <- round(n/2)
  
  #setup basis
  B <- buildBPbasis(xs,M,deriv=0)$B
  A <- rbind(0,buildBPbasis(xs,M,deriv=1)$A)
  A[1,1] <- 1
  BA <- B %*% solve(A)
  
  if(missing(priors)){
    # prior on alpha
    alpha_a <- 1
    alpha_b <- 1
    # prior on pi
    pi_a <- 1
    pi_b <- 1
    mu0 <- 0.5
    phi0 <- 0.25
    intsd <- 2
  }else{
    if(is.null(priors$alpha_a)){
      alpha_a <- 1
    }else{
      alpha_a <- max(priors$alpha_a[1],1)
    }
    if(is.null(priors$alpha_b)){
      alpha_b <- 1
    }else{
      alpha_b <- max(priors$alpha_b[1],1)
    }
    if(is.null(priors$pi_a)){
      pi_a <- 1
    }else{
      pi_a <- max(priors$pi_a[1],1/2)
    }
    if(is.null(priors$pi_b)){
      pi_b <- 1
    }else{
      pi_b <- max(priors$pi_b[1],1/2)
    }
    if(is.null(priors$mu0)){
      mu0 <- 0.5
    }else{
      mu0 <- priors$mu0[1]
    }
    if(is.null(priors$phi0)){
      phi0 <- 0.25
    }else{
      phi0 <- priors$phi0[1]
    }
    if(is.null(priors$intvar)){
      intsd <- 2
    }else{
      intsd <- priors$intsd[1]
    }
    
  }
  

  modelstatement <- paste0(
  "  Y_i     ~  N[basis(x_i,order=",M,")*beta,sigma^2]\n",
  "  beta   =  A * theta\n",
  "  theta_0  ~  Normal(0,",intsd,"^2)  (intercept)\n",
  "  theta_j  ~  pi*DiracDelta_0 + (1-pi)*DP(alpha P_0)  (j=1,...,M)\n",
  "  P_0     =  N(",mu0,",",phi0,"^2)\n",
  "  alpha   ~  gamma(",alpha_a,",",alpha_b,")\n",
  "  pi      ~  beta(",pi_a,",",pi_b,")\n")
  cat("Model:\n")
  cat(modelstatement)

  
  # MCMC
  up <- mcmc(ys, BA, M, pi_a, pi_b, alpha_a, alpha_b,n,mu0,1/phi0^2,1/intsd^2,nthin,round(niter/nthin))
  beta_keep <- up[[1]]
  hyper_keep <- up[[2]]


  colnames(hyper_keep) <- c("alpha","sigma", "phi")
  
  
  # format and rescale
  hyper_keep[,"sigma"] <- sd(y)/sqrt(hyper_keep[,"sigma"])
  hyper_keep[,"phi"] <- 1/sqrt(hyper_keep[,"phi"])
  beta_keep <- beta_keep*sd(y)
  beta_keep[,1] <- beta_keep[,1] + mean(y)
  
  # estiamte function
  beta_post <- beta_keep[-c(1:round(nburn/nthin)),]
  gamma <- solve(A) %*% t(beta_keep[-c(1:round(nburn/nthin)),])
  f <- B%*%gamma
  fitted <- data.frame(fitted=rowMeans(f),
                       lower=apply(f,1,quantile, 0.025),
                       upper=apply(f,1,quantile, 0.975),
                       x=x)
  
  # estiamte first derivative
  fprime <- M*buildBPbasis(xs,M-1,deriv=0)$B%*%t(beta_keep[-c(1:round(nburn/nthin)),])[-1,]
  derivative <- data.frame(fitted=rowMeans(fprime),
                           lower=apply(fprime,1,quantile, 0.025),
                           upper=apply(fprime,1,quantile, 0.975),
                           x=x)
  
  out <- list(fitted=fitted, 
              derivative=derivative, 
              f_all=f,
              fprime_all=fprime,
              basis=B, 
              A=A, 
              beta=beta_keep, 
              gamma=solve(A) %*% t(beta_keep), 
              parms=as.data.frame(hyper_keep),
              specs=list(priors=list(mu0=mu0, phi0=phi0, pi_a=pi_a, pi_b=pi_b, alpha_a=alpha_a, alpha_b=alpha_b), nburn=nburn,niter=niter,nthin=nthin, M=M),
              call = match.call(),
              modelstatement=modelstatement,
              y=y,x=x)
  class(out) <- "bnmr"
  
  return(out)
}


