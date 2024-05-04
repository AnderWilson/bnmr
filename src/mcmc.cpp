#include "RcppArmadillo.h"
using namespace Rcpp;



// [[Rcpp::interfaces(r, cpp)]]





double rtnorm1(double a) {
  
  double x = 0.0;
  double y;
  
  if(a<0){
    // normal rejection sampling
    while(x==0){
      y = R::rnorm(0,1);
      if(y>a){
        x=y;
      }
    }
  }else if(a<0.2570){
    // half-normal rejection sampling
    while(x==0){
      y = fabs(R::rnorm(0,1));
      if(y>a){
        x=y;
      }
    }
  }else{
    // one-sided translated-exponential rejection sampling
    while(x==0){
      double lambdastar = (a + sqrt(a*a+4.0))/2.0;
      y = R::rexp(1)/lambdastar + a;
      if(R::runif(0,1) < exp(-0.5*y*y + lambdastar*y - 0.5*lambdastar + log(lambdastar))){
        x=y;
      }
    }
  }
  
  return x;
}


double rtnorm2(double a, double b) {
  
  double x = 0.0;
  double y;
  
  if(b > a + sqrt(2*M_PI)){
    // normal rejection sampling
    while(x==0){
      y = R::rnorm(0,1);
      if(y>a & y<b){
        x=y;
      }
    }
  }else{
    // uniform rejection sampling
    while(x==0){
      y = R::runif(a,b);
      if(R::runif(0,1)<exp(-y*y/2)){
        x=y;
      }
    }
  }
  
  return x;
}


double rtnorm3(double a, double b) {
  
  double x = 0.0;
  double y;
  
  if(a<0.2570){
    if(b>a + sqrt(M_PI/2)*exp(a*a/2)){
      // half-normal rejection sampling
      while(x==0){
        y = fabs(R::rnorm(0,1));
        if(y>a & y<b){
          x=y;
        }
      }
    }else{
      // uniform rejection sampling
      while(x==0){
        y = R::runif(a,b);
        if(R::runif(0,1)<exp((a*a-y*y)/2)){
          x=y;
        }
      }
    }
  }else{
    if(b>a + 2/(a + sqrt(a*a + 4))*exp((a*a-a*sqrt(a*a+4))/4 + 0.5)){
      // two-sided translated-exponential rejection sampling
      while(x==0){
        double lambdastar = (a + sqrt(a*a+4))/2;
        y =  a-log(R::runif(exp((a-b)*lambdastar),1))/lambdastar;
        if(R::runif(0,1) < exp(-y*y/2 + lambdastar*y - lambdastar/2 + log(lambdastar))){
          x=y;
        }
      }
    }else{
      // uniform rejection sampling
      while(x==0){
        y = R::runif(a,b);
        if(R::runif(0,1)<exp((a*a-y*y)/2)){
          x=y;
        }
      }
    }
  }
  
  return x;
}


// [[Rcpp::export]]
double rtnormc(double a, double b) {
  
  double x=0;
  
  if(!std::isfinite(a) && !std::isfinite(b)){
    // Not truncated
    x=R::rnorm(0,1);
  }else if(std::isfinite(a) && !std::isfinite(b)){
    // CASE 1
    x= rtnorm1(a);
  }else if(std::isfinite(a) && std::isfinite(b) && a < 0 && b > 0){
    // CASE 2
    x=rtnorm2(a,b);
  }else if(std::isfinite(a) && std::isfinite(b) && a >= 0){
    // CASE 3
    x=rtnorm3(a,b);
  }else if(!std::isfinite(a) && std::isfinite(b)){
    // CASE 4
    x=-rtnorm1(-b);
  }else if(std::isfinite(a) && std::isfinite(b) && a < 0){
    // CASE 5
    x=-rtnorm3(-b,-a);
  }
  
  return x;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rtmvnormc(arma::vec m, arma::mat V, arma::vec w, int n) {
  
  arma::mat R = chol(V, "lower");
  arma::vec x = solve(trimatl(R), (w - m));
  arma::vec a = -m;
  a(0) = -arma::datum::inf;
  int i = 0;
  int j = 0;
  double amin = -arma::datum::inf;
  double amax =  arma::datum::inf;
  
  for(i=0; i<n; i++){
    amin = -arma::datum::inf;
    amax =  arma::datum::inf;
    for(j=i; j<n; j++){
      x(i) = 0;
      if(R(j,i)>0){
        amin = std::max(amin, as_scalar( (a(j)-R.row(j)*x) / R.submat(j,i,j,i) ) );
      }else if(R(j,i)<0){
        amax = std::min(amax, as_scalar( (a(j)-R.row(j)*x) / R.submat(j,i,j,i) ) );   
      }
    }
    x(i) = rtnormc(amin,amax);
  }
  
  w = R*x + m;
  
  return w;
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double updatealpha(double alpha, arma::vec nc, double alpha_a, double alpha_b, int M) {
  
  if(nc.n_elem==1){
    // no non-null groups so pull alpha from prior
    alpha = R::rgamma(alpha_a,alpha_b);
  }else{
    // some non-null groups so use fill conditional
    double eta = R::rbeta(alpha+1,(M-nc(0)));
    double pieta = (alpha_a+nc.n_elem-2)/((M-nc(0))*(alpha_b-log(eta))) / ((alpha_a+nc.n_elem-2)/((M-nc(0))*(alpha_b-log(eta)))+1);
    if(R::runif(0,1)<pieta){
      alpha = R::rgamma(alpha_a+nc.n_elem-1 , 1/(alpha_b-log(eta)) );
    }else{
      alpha = R::rgamma(alpha_a+nc.n_elem-2 , 1/(alpha_b-log(eta)) );
    }
  }
  return alpha;
  
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List betatothetaandS(arma::vec beta) {
  int no0 = 1-1*any(beta == 0.0);
  arma::vec theta;
  if(no0==1){
    // add an empty location for zeros
    arma::vec beta2 = beta;
    beta2(0) = 0;
    theta = arma::unique(beta2.subvec(0,beta.n_elem-1));  
  }else{
    // there are zeros so we don't need to add it.
    theta = arma::unique(beta.subvec(1,beta.n_elem-1));  
  }
  arma::vec nc(theta.n_elem, arma::fill::zeros);
  arma::vec S(beta.n_elem-1, arma::fill::zeros);
  arma::uvec mtch;
  for(int i=no0; i<theta.n_elem; i++ ){
    // start at no0 to leave a zero group when there otherwise aren't any zeros
    mtch = find(beta==theta(i));
    S.elem(mtch-1).fill(i);
    nc(i) = mtch.n_elem;
  }
  
  
  return  List::create(Named("theta") = arma::vectorise(theta),
                       Named("nc")=nc,
                       Named("S")=S);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double updatebeta(arma::vec y, arma::mat BA, arma::vec beta, 
                  arma::vec nc,
                  arma::vec theta,
                  int M,
                  double pi_a,
                  double pi_b,
                  double alpha,
                  double sig2inv,
                  double phi2inv,
                  double mu,
                  int j,
                  int n) {
  
  
  int k;
  int i;
  
  // residuals without column j
  arma::vec r = y-BA*beta;
  r += BA.col(j) * beta(j);
  
  // remove jth coef from nc count
  nc(find(theta==beta(j))) -= 1;
  arma::vec loglik(nc.n_elem+1, arma::fill::zeros);
  
  double v_prime = 1/(phi2inv + sig2inv*arma::sum(arma::square(BA.col(j))));
  // double m_prime = v_prime*(phi2inv*mu + sig2inv*arma::sum(r));
  double m_prime = v_prime*(phi2inv*mu + sig2inv*as_scalar(r.t()*BA.col(j)));
  
  
  loglik(0) = log(nc(0) + pi_a) - log(M -1 + pi_a + pi_b) ;
  for(k=1; k<nc.n_elem; k++){
    loglik(k) = log(M-nc(0) + pi_b) + log(nc(k)) 
    - log(M -1 + pi_a + pi_b) - log(M-nc(0)  + alpha);
  }
  loglik(nc.n_elem) = log(M-nc(0) + pi_b) + log(alpha) 
    - log(M -1 + pi_a + pi_b) - log(M-nc(0)  + alpha)
    + R::pnorm(0,m_prime,sqrt(v_prime),0, 1)
    - R::pnorm(0,mu,sqrt(1/phi2inv), 0, 1)
    - .5*log(phi2inv) + (m_prime*m_prime)/(2*v_prime) - (phi2inv*mu*mu)/2 + log(v_prime)/2;
    
    for(i=0; i<n; i++){
      loglik(0) += R::dnorm(r(i), 0, sqrt(1/sig2inv), 1);
      for(k=1; k<nc.n_elem; k++){
        loglik(k) += R::dnorm(r(i),  BA(i,j)*theta(k), sqrt(1/sig2inv), 1);
      }
      loglik(nc.n_elem) += R::dnorm(r(i),0,sqrt(1/sig2inv), 1);
    }
    
    arma::uvec Snew = find(R::runif(0,1) < cumsum(exp(loglik - (log(arma::sum(exp(loglik - arma::max(loglik)))) + arma::max(loglik)))),1,"first");
    double theta_new;
    
    if(Snew(0)<theta.n_elem){
      theta_new = theta(Snew(0));
    }else{
      double vj = 1/(sig2inv*as_scalar( (BA.col(j)).t()*BA.col(j) ) + phi2inv);
      double mj = vj * (sig2inv*as_scalar((BA.col(j)).t()*r) + phi2inv*mu);
      theta_new = rtnormc(-mj/sqrt(vj),arma::datum::inf)*sqrt(vj) + mj;
    }
    return theta_new;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mcmc(arma::vec y, 
              arma::mat BA, 
              int M,
              double pi_a,
              double pi_b,
              double alpha_a, 
              double alpha_b,
              int n,
              double mu0,
              double phi2inv,
              double intvarinv,
              int thin,
              int loops) {
  
  
  double sig2inv = 1;
  double alpha=1;
  
  arma::vec beta;
  beta.zeros(M+1);
  for(int i=1; i<M+1; i++){
    beta(i) = R::runif(0,1)/10;
  }
  
  arma::vec inttheta=beta;
  Rcpp::List Stheta = betatothetaandS(beta);
  arma::vec theta = Stheta["theta"];
  arma::vec S = Stheta["S"];
  arma::mat trans;
  arma::vec res;
  arma::mat BAtrans;
  arma::mat BABA;
  arma::mat v;
  arma::mat m;
  int st;
  double intercept;
  
  
  arma::mat beta_keep(loops,beta.n_elem);
  arma::mat hyper_keep(loops,3);  
    
  for(int outr = 0; outr < loops; outr ++){
  for(int tn = 0; tn < thin; tn ++){

  // update precision
  res = y-BA * beta;
  sig2inv = R::rgamma(0.01 + n/2, 1/( 0.01+ as_scalar( res.t()*res )/2));

  // update prior on alpha
  // assume alpha~gamma(1,1), k=length(nc)-1
  alpha = updatealpha(alpha, Stheta["nc"], alpha_a, alpha_b, M);
      
  // update betas
  st = ceil(R::runif(0,M));
  for(int j=st; j<M+1; j++){
    beta(j) = updatebeta(y, BA, beta, Stheta["nc"], Stheta["theta"],
                     M, pi_a, pi_b, alpha, sig2inv, phi2inv, mu0, j, n);
    Stheta = betatothetaandS(beta);
  }
  for(int j=1; j<st; j++){
    beta(j) = updatebeta(y, BA, beta, Stheta["nc"], Stheta["theta"],
         M, pi_a, pi_b, alpha, sig2inv, phi2inv, mu0, j, n);
    Stheta = betatothetaandS(beta);
  }

  // make transformation matrix (beta to theta)
  arma::vec theta = Stheta["theta"];
  arma::vec S = Stheta["S"];
  trans.zeros(M+1, theta.n_elem);
  trans(0,0) = 1;
  for(int i=1; i<theta.n_elem; i++){
    trans.elem(find(S==i) + i*(M+1) +1 ) += 1;
  }
  
  // transform basis from beta to theta
  BAtrans = BA*trans;
  
  // add prior variance
  BABA = sig2inv*(BAtrans.t()*BAtrans);
  BABA.diag() += phi2inv;
  BABA(0,0) -= (phi2inv - intvarinv);
  BABA(0,0) += 0.01;
  
  v =inv_sympd(BABA);
  m =  sig2inv*(BAtrans.t()*y);
  if(theta.n_elem>1){
    m.rows(1,theta.n_elem-1) += phi2inv*mu0;
  }
  m = v * m;
  
  // update inttheta
  intercept = inttheta(0);
  inttheta = theta;
  inttheta(0) = intercept;
  
  // update theta including intercept and then phi2inv
  if(theta.n_elem>1){
    inttheta.subvec(1,inttheta.n_elem-1) = theta.subvec(1,inttheta.n_elem-1);
    inttheta = rtmvnormc(m,v,inttheta,inttheta.n_elem );
  }else{
    inttheta = R::rnorm(0,1) * sqrt(v(0,0)) + m(0,0);
  }
  
  // make beta
  beta = trans*inttheta;
  
  }
  
  
  beta_keep.row(outr) = beta.t();
  hyper_keep(outr,0) = alpha;
  hyper_keep(outr,1) = sig2inv;
  hyper_keep(outr,2) = phi2inv;
  }
  
  return  List::create(Named("beta_keep") = beta_keep,
                       Named("hyper_keep")=hyper_keep);
}

