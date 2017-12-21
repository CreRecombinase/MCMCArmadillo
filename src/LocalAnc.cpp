// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends("RcppArmadillo")]]

const double log2pi = std::log(2.0 * M_PI);


// [[Rcpp::export]]
arma::mat em_with_zero_mean_c(arma::mat y,
                              int maxit){
  //EM for empirical covariance matrix when y has missing values
  int orig_p = y.n_cols;
  arma::vec vars = arma::zeros<arma::vec>(orig_p);
  for (int i=0; i < orig_p; ++i){
    arma::vec ycol = y.col(i);
    arma::uvec finiteind = find_finite(ycol);
    arma::vec yy = ycol(finiteind);
    vars(i) = sum((yy-mean(yy))%(yy-mean(yy)));
  }
  arma::uvec valid_ind = find(vars>1e-6);
  y = y.cols(valid_ind);
  int p = y.n_cols;
  int n = y.n_rows;
  arma::rowvec mu = arma::zeros<arma::rowvec>(p);
  arma::mat y_imputed = y;
  for (int j = 0; j < p; ++j){
    arma::uvec colind = arma::zeros<arma::uvec>(1);
    colind(0) = j;
    arma::uvec nawhere = find_nonfinite(y_imputed.col(j));
    arma::uvec nonnawhere = find_finite(y_imputed.col(j));
    arma::vec tempcolmean = mean(y_imputed(nonnawhere, colind), 0);
    y_imputed(nawhere, colind).fill(tempcolmean(0));
  }
  arma::mat oldSigma = y_imputed.t() * y_imputed / n;
  arma::mat Sigma = oldSigma;
  double diff = 1;
  int it = 1;
  while (diff>0.001 && it < maxit){
    arma::mat bias = arma::zeros<arma::mat>(p,p);
    for (int i=0; i<n; ++i){
      arma::rowvec tempdat = y.row(i);
      arma::uvec ind = find_finite(tempdat);
      arma::uvec nind = find_nonfinite(tempdat);
      if (0 < ind.size() && ind.size() < p){
        //MAKE THIS PART FASTER
        bias(nind, nind) += Sigma(nind, nind) - Sigma(nind, ind) * (Sigma(ind, ind).i()) * Sigma(ind, nind);
        arma::uvec rowind = arma::zeros<arma::uvec>(1);
        rowind(0) = i;
        arma::mat yvec = y(rowind, ind);
        //MAKE THIS PART FASTER
        y_imputed(rowind, nind) = (Sigma(nind, ind)*(Sigma(ind, ind).i())*y(rowind, ind).t()).t();
      }
    }
    Sigma = (y_imputed.t() * y_imputed + bias)/n;
    arma::mat diffmat = (Sigma-oldSigma);
    arma::mat diffsq = diffmat%diffmat;
    diff = accu(diffsq);
    oldSigma = Sigma;
    it = it + 1;
  }
  arma::mat finalSigma = arma::zeros<arma::mat>(orig_p, orig_p);
  finalSigma.submat(valid_ind, valid_ind.t()) = Sigma;
  return finalSigma;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n,
                      arma::vec mu,
                      arma::mat Sigma) {
  //returns random multivariate normal vectors with mean mu and covariance Sigma
  //input : integer n for the number of vectors you'd like to draw
  //      : vector mu for the mean
  //      : matrix Sigma for the covariance - needs to be psd
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dmvnrm_arma(arma::rowvec x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  //returns the density of a multivariate normal vector
  //input : a rowvector x whose density you'd like to know
  //      : a rowvector mean for the mean of mvn
  //      : a matrix Sigma for covariance, needs to be psd
  //      : a boolean logd, true if you like the log density
  int xdim = x.n_cols;
  if(xdim==0){return 0;}
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  arma::vec z = rooti*arma::trans(x-mean);
  out = constants - 0.5*arma::sum(z%z)+rootisum;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double get_sigmabeta_from_h_c(double h,
                              arma::vec gam,
                              arma::mat Sigma,
                              arma::vec X,
                              int T){
  //convert h to sigmabeta conditioning on gamma and Sigma
  int n = X.size();
  arma::vec ds = Sigma.diag();
  double num = h * sum(ds);
  arma::uvec ind = find(gam == 1);
  double denom = (1-h)*sum(ds(ind)) * sum(X%X)/n;
  return num/denom;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_h_from_sigmabeta_c(arma::vec X, double sigmabeta,
                              arma::mat Sigma, arma::vec gam,
                              int n, int T){
  //converts sigmabeta to h conditioning on gamma and Sigma
  arma::uvec ind = find(gam==1);
  arma::vec ds = Sigma.diag();
  double num = sum(X%X)/n * sum(ds(ind)) * sigmabeta;
  double denom = num + sum(ds);
  return num/denom;
}


// [[Rcpp::export]]
arma::vec get_target_c(arma::vec X, arma::mat Y, double sigmabeta,
                       arma::mat Sigma, arma::vec gam, arma::vec beta){
  //get the target likelihood circumventing the missing value issue
  int T = Y.n_cols;
  int n = Y.n_rows;
  double L = 0;
  double B = 0;
  double G = 0;
  for (int i=0; i < n; ++i){
    arma::uvec naind = find_finite(Y.row(i).t());
    if(naind.size()>0){
      arma::uvec rowind = arma::zeros<arma::uvec>(1);
      rowind(0) = i;
      arma::rowvec Ytemp = Y(rowind, naind.t());
      L = L + dmvnrm_arma(Ytemp,
                          X(i)*beta(naind).t(),
                          Sigma(naind,naind.t()),
                          true);
    }
  }
  arma::uvec ind = find(gam==1);
  int s = ind.size();
  if(s>0){
    arma::vec ds = Sigma.diag();
    for (int j=0; j<s; ++j){
      int newind = ind(j);
      B = B + R::dnorm(beta(newind), 0, sqrt(sigmabeta*ds(newind)), true);
    }
  }
  G = log(std::tgamma(s+1)*std::tgamma(T-s+1)/std::tgamma(T+2));
  arma::vec out = arma::zeros<arma::vec>(3);
  out(0) = L;
  out(1) = B;
  out(2) = G;
  return out;
}

// [[Rcpp::export]]
int sample_index(int size, NumericVector prob = NumericVector::create()){
  //sample one number from 1:size
  arma::vec sequence = arma::linspace<arma::vec>(1, size, size);
  arma::vec out = Rcpp::RcppArmadillo::sample(sequence, size, false, prob);
  return out(0);
}


// [[Rcpp::export]]
Rcpp::List update_gamma_c(arma::vec X, arma::mat Y, arma::vec gam){
  //update gamma once
  int changeind = 0;
  arma::vec newgam = gam;
  int T = gam.size();
  arma::vec prob = arma::zeros<arma::vec>(2);
  prob(0) = 0.5; prob(1) = 0.5;
  arma::uvec ind0 = find(gam==0);
  arma::uvec ind1 = find(gam==1);
  int s = ind1.size();
  NumericVector prob2 = wrap(prob);
  int cas = sample_index(2, prob2);
  if(s==0){
    cas = 1;
  }else if(s==T){
    cas = 2;
  }
  if (cas==1){
    arma::rowvec marcor = arma::zeros<arma::rowvec>(ind0.size());
    for (int t=0; t<ind0.size(); ++t){
      arma::vec temp = Y.col(ind0(t));
      arma::uvec tempind = find_finite(temp);
      marcor(t) = abs(sum(temp(tempind)%X(tempind)))/tempind.size();
    }
    int add = 1;
    if(s<(T-1)){
      NumericVector marcor2 = wrap(marcor);
      add = sample_index(ind0.size(), marcor2);
    }
    newgam(ind0(add-1)) = 1;
    changeind = ind0(add-1);
  }else if(cas==2){
    int remove = sample_index(s);
    remove = ind1(remove-1);
    newgam(remove) = 0;
    changeind = remove;
  }
  return(
    Rcpp::List::create(
      Rcpp::Named("gam") = newgam,
      Rcpp::Named("changeind") = changeind)
  );
}

// [[Rcpp::export]]
arma::vec betagam_accept_c(arma::vec X,
                           arma::mat Y,
                           double sigmabeta1,
                           arma::mat inputSigma,
                           double Vbeta,
                           arma::vec gam1,
                           arma::vec beta1,
                           arma::vec gam2,
                           arma::vec beta2,
                           int changeind,
                           int change){
  //compute the target likelihood and the proposal ratio
  //to decide if you should accept the proposed beta and gamma
  double newtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam2,beta2));
  double oldtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam1,beta1));
  double proposal_ratio = R::dnorm(beta1(changeind)-beta2(changeind),0,sqrt(Vbeta),true);
  int T = gam1.size();
  int s1 = sum(gam1==1);
  int s2 = sum(gam2==1);
  arma::rowvec marcor = arma::zeros<arma::rowvec>(T);
  for (int t=0; t<T; ++t){
    arma::vec temp = Y.col(t);
    arma::uvec tempind = find_finite(temp);
    marcor(t) = abs(sum(temp(tempind)%X(tempind)))/tempind.size();
  }
  if(change==1){
    arma::uvec ind1 = find(gam1==0);
    double temp1 = marcor(changeind)/sum(marcor(ind1));
    proposal_ratio = -log(temp1)-log(s2)-proposal_ratio;
  }else{
    arma::uvec ind2 = find(gam2==0);
    double temp2 = marcor(changeind)/sum(marcor(ind2));
    proposal_ratio = log(temp2)+log(s1)+proposal_ratio;
  }
  double final_ratio = newtarget-oldtarget+proposal_ratio;
  arma::vec out = arma::zeros<arma::vec>(4);
  out(0) = final_ratio;
  out(1) = newtarget;
  out(2) = oldtarget;
  out(3) = proposal_ratio;
  return(out);
}

// [[Rcpp::export]]
Rcpp::List update_betagam_c(arma::vec X,
                            arma::mat Y,
                            arma::vec gam1,
                            arma::vec beta1,
                            arma::mat Sigma,
                            double sigmabeta,
                            double Vbeta,
                            int bgiter){
  //update and beta and gamma 'bgiter' times
  for (int i=1; i<bgiter; ++i){
    Rcpp::List temp = update_gamma_c(X,Y,gam1);
    arma::vec gam2 = as<arma::vec>(temp["gam"]);
    arma::vec beta2 = beta1 % gam2;
    arma::uvec ind = find(gam2==1);
    beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
    int changeind = temp["changeind"];
    int change = gam2(changeind);
    arma::vec A = betagam_accept_c(X,Y,sigmabeta,
                                   Sigma,Vbeta,
                                   gam1,beta1,
                                   gam2,beta2,
                                   changeind,change);
    NumericVector check2 = runif(1);
    double check = check2(0);
    if(exp(A(0))>check){
      gam1 = gam2; beta1 = beta2;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("gam") = gam1,
    Rcpp::Named("beta") = beta1
  );
}



// [[Rcpp::export]]
Rcpp::List update_h_c(double initialh, int hiter, arma::vec gam, arma::vec beta,
                      arma::mat Sig, arma::vec X, int T){
  double h1 = initialh;
  double sigbeta1 = get_sigmabeta_from_h_c(initialh, gam, Sig, X, T);
  arma::vec lik = arma::zeros<arma::vec>(hiter);
  arma::vec ds = Sig.diag();
  for (int i=1; i<hiter; ++i){
    double h2 = h1;
    NumericVector(rr) = runif(1, -0.1, 0.1);
    double r = rr(0);
    h2 = h2 + r;
    if(h2<0){h2 = abs(h2);}
    if(h2>1){h2 = 2-h2;}
    arma::uvec ind = find(gam==1);
    double sigmabeta1 = get_sigmabeta_from_h_c(h1, gam, Sig, X, T);
    double sigmabeta2 = get_sigmabeta_from_h_c(h2, gam, Sig, X, T);
    double lik1 = 0; double lik2 = 0;
    for (int j=0; j < ind.size(); ++j){
      int newind = ind(j);
      lik1 = lik1 + R::dnorm(beta(newind), 0, sqrt(sigmabeta1*ds(newind)), true);
      lik2 = lik2 + R::dnorm(beta(newind), 0, sqrt(sigmabeta2*ds(newind)), true);
    }
    double acceptanceprob = exp(lik2-lik1);
    arma::vec ee = runif(1);
    double e = ee(0);
    if(e<acceptanceprob){
      h1 = h2; sigbeta1 = sigmabeta2;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("h") = h1,
    Rcpp::Named("sigbeta") = sigbeta1
  );
}



// [[Rcpp::export]]
arma::cube rinvwish_c(int n, int v, arma::mat S){
  //draw a matrix from inverse wishart distribution with parameters S and v
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = chol(inv_sympd(S), "lower");
  arma::cube sims(p, p, n, arma::fill::zeros);
  for(int j = 0; j < n; j++){
    arma::mat A(p,p, arma::fill::zeros);
    for(int i = 0; i < p; i++){
      int df = v - (i + 1) + 1; //zero-indexing
      A(i,i) = sqrt(R::rchisq(df));
    }
    for(int row = 1; row < p; row++){
      for(int col = 0; col < row; col++){
        A(row, col) = R::rnorm(0,1);
      }
    }
    arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
    sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}

// [[Rcpp::export]]
arma::mat update_Sigma_c(int n, int nu, arma::vec X, arma::vec beta, arma::mat Phi, arma::mat Y){
  int T = Y.n_cols;
  arma::mat X2 = arma::zeros<arma::mat>(n, 1);
  X2.col(0) = X;
  arma::mat beta2 = arma::zeros<arma::mat>(T,1);
  beta2.col(0) = beta;
  arma::mat r = Y - X2 * beta2.t();
  arma::mat emp = em_with_zero_mean_c(r,100);
  arma::cube res = rinvwish_c(1, n+nu, emp*n + Phi*nu);
  return res.slice(0);
}

// [[Rcpp::export]]
Rcpp::List update_gamma_sw_c(arma::vec X,
                             arma::mat Y,
                             arma::vec gam,
                             arma::rowvec marcor){
  int changeind = 0;
  //flip marcor
  arma::rowvec marcor2 = (max(marcor) + min(marcor)) - marcor;
  arma::vec newgam = gam;
  int T = gam.size();
  arma::vec prob = arma::zeros<arma::vec>(2);
  prob(0) = 0.5; prob(1) = 0.5;
  arma::uvec ind0 = find(gam==0);
  arma::uvec ind1 = find(gam==1);
  int s = ind1.size();
  NumericVector prob2 = wrap(prob);
  int cas = sample_index(2, prob2);
  if(s==0){
    cas = 1;
  }else if(s==T){
    cas = 2;
  }
  if (cas==1){
    int add = 1;
    if(s<(T-1)){
      arma::vec mc = marcor(ind0);
      NumericVector mc2 = wrap(mc);
      add = sample_index(ind0.size(), mc2);
    }
    newgam(ind0(add-1)) = 1;
    changeind = ind0(add-1);
  }
  if(cas==2){
    int remove=1;
    if(s > 1){
      arma::vec mc = marcor2(ind1);
      NumericVector mc2 = wrap(mc);
      remove = sample_index(ind1.size(), mc2);
    }
    newgam(ind1(remove-1)) = 0;
    changeind = ind1(remove-1);
  }
  return(
    Rcpp::List::create(
      Rcpp::Named("gam") = newgam,
      Rcpp::Named("changeind") = changeind)
  );
}

// [[Rcpp::export]]
arma::vec betagam_accept_sw_c(arma::vec X,
                              arma::mat Y,
                              double sigmabeta1,
                              arma::mat inputSigma,
                              double Vbeta,
                              arma::vec gam1,
                              arma::vec beta1,
                              arma::vec gam2,
                              arma::vec beta2,
                              int changeind,
                              int change){
  double newtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam2,beta2));
  double oldtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam1,beta1));
  double proposal_ratio = R::dnorm(beta1(changeind)-beta2(changeind),0,sqrt(Vbeta),true);
  int T = gam1.size();
  arma::rowvec marcor = arma::zeros<arma::rowvec>(T);
  for (int t=0; t<T; ++t){
    arma::vec temp = Y.col(t);
    arma::uvec tempind = find_finite(temp);
    marcor(t) = abs(sum(temp(tempind)%X(tempind)))/tempind.size();
  }
  arma::rowvec marcor2 = min(marcor)+max(marcor)-marcor;
  if(change==1){
    arma::uvec ind1 = find(gam1==0);
    arma::uvec ind2 = find(gam2==1);
    double tempadd = marcor(changeind)/sum(marcor(ind1));
    double tempremove = marcor2(changeind)/sum(marcor2(ind2));
    proposal_ratio = -log(tempadd)+log(tempremove)-proposal_ratio;
  }else{
    arma::uvec ind1 = find(gam1==1);
    arma::uvec ind2 = find(gam2==0);
    double tempadd = marcor(changeind)/sum(marcor(ind2));
    double tempremove = marcor2(changeind) / sum(marcor2(ind1));
    proposal_ratio = log(tempadd)-log(tempremove)+proposal_ratio;
  }
  double final_ratio = newtarget-oldtarget+proposal_ratio;
  arma::vec out = arma::zeros<arma::vec>(4);
  out(0) = final_ratio;
  out(1) = newtarget;
  out(2) = oldtarget;
  out(3) = proposal_ratio;
  return(out);
}

// [[Rcpp::export]]
Rcpp::List update_betagam_sw_c(arma::vec X,
                               arma::mat Y,
                               arma::vec gam1,
                               arma::vec beta1,
                               arma::mat Sigma,
                               arma::rowvec marcor,
                               double sigmabeta,
                               double Vbeta,
                               int bgiter,
                               int smallworlditer){
  int T = gam1.size();
  arma::mat outgamma = arma::zeros<arma::mat>(T,bgiter);
  arma::mat outbeta = arma::zeros<arma::mat>(T,bgiter);
  outgamma.col(0) = gam1;
  outbeta.col(0) = beta1;
  arma::vec tar = arma::zeros<arma::vec>(bgiter);
  for (int i=1; i<bgiter; ++i){
    Rcpp::List temp = update_gamma_sw_c(X,Y,outgamma.col(i-1), marcor);
    arma::vec gam1 = outgamma.col(i-1);
    arma::vec beta1 = outbeta.col(i-1);
    //small world proposal
    if(i%10==0){
      double proposal_ratio = 0;
      arma::vec betatemp1 = beta1;
      arma::vec gamtemp1 = gam1;
      arma::vec betatemp2 = betatemp1;
      arma::vec gamtemp2 = gamtemp1;
      for (int j=0; j < smallworlditer; ++j){
        Rcpp::List temp = update_gamma_sw_c(X,Y,gamtemp1, marcor);
        arma::vec gamtemp2 = as<arma::vec>(temp["gam"]);
        arma::vec betatemp2 = betatemp1 % gamtemp2;
        arma::uvec ind = find(gamtemp2==1);
        betatemp2(ind) = betatemp1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
        int changeind = temp["changeind"];
        int change = gamtemp2(changeind);
        double proposaliter = R::dnorm(betatemp1(changeind)-betatemp2(changeind),
                                       0,sqrt(Vbeta), true);
        arma::rowvec marcor2 = -marcor + max(marcor) + 0.01;
        if(change==1){
          arma::uvec ind1 = find(gamtemp1==0);
          arma::uvec ind2 = find(gamtemp2==1);
          double tempadd = marcor(changeind)/sum(marcor(ind1));
          double tempremove = marcor2(changeind)/sum(marcor2(ind2));
          proposaliter = -log(tempadd)-log(tempremove)-proposaliter;
        }else{
          arma::uvec ind1 = find(gamtemp1==1);
          arma::uvec ind2 = find(gamtemp2==0);
          double tempadd = marcor(changeind)/sum(marcor(ind2));
          double tempremove = marcor2(changeind) / sum(marcor2(ind1));
          proposaliter = log(tempadd)+log(tempremove)+proposaliter;
        }
        proposal_ratio = proposal_ratio + proposaliter;
        gamtemp1 = gamtemp2; betatemp1 = betatemp2;
      }
      arma::vec gam2 = gamtemp2;
      arma::vec beta2 = betatemp2;
      double newtarget = sum(get_target_c(X,Y,sigmabeta,Sigma,gam2,beta2));
      double oldtarget = sum(get_target_c(X,Y,sigmabeta,Sigma,gam1,beta1));
      double A = newtarget-oldtarget + proposal_ratio;
      arma::vec check2 = runif(1,0,1); double check = check2(0);
      if(exp(A) > check){
        tar(i) = newtarget;
        outgamma.col(i)= gam2;
        outbeta.col(i) = beta2;
      }else{
        tar(i) = oldtarget;
        outgamma.col(i) = gam1;
        outbeta.col(i) = beta1;
      }
    }else{
      Rcpp::List temp = update_gamma_sw_c(X,Y,gam1, marcor);
      arma::vec gam2 = as<arma::vec>(temp["gam"]);
      arma::vec beta2 = beta1 % gam2;
      arma::uvec ind = find(gam2==1);
      beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
      int changeind = temp["changeind"];
      int change = gam2(changeind);
      arma::vec A = betagam_accept_c(X,Y,sigmabeta,
                                     Sigma,Vbeta,
                                     gam1,beta1,
                                     gam2,beta2,
                                     changeind,change);
      NumericVector check2 = runif(1);
      double check = check2(0);
      if(exp(A(0))>check){
        tar(i) = A(1);
        outgamma.col(i) = gam2; outbeta.col(i) = beta2;
      }else{
        tar(i) = A(2);
        outgamma.col(i) = gam1; outbeta.col(i) = beta1;
      }
    }
  }
  
  arma::vec outgamma2 = outgamma.col(bgiter-1);
  arma::vec outbeta2 = outbeta.col(bgiter-1);
  return Rcpp::List::create(
    Rcpp::Named("gam")= outgamma2,
    Rcpp::Named("beta") = outbeta2,
    Rcpp::Named("tar") = tar
  );
}



// [[Rcpp::export]]
Rcpp::List doMCMC_c(arma::vec X,
                    arma::mat Y,
                    int n,
                    int T,
                    arma::mat Phi,
                    int nu,
                    arma::vec initialbeta,
                    arma::vec initialgamma,
                    arma::mat initialSigma,
                    double initialsigmabeta,
                    arma::rowvec marcor,
                    double Vbeta,
                    int niter,
                    int bgiter,
                    int hiter,
                    int switer){
  //empty arrays to save values
  arma::mat outbeta = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb = arma::zeros<arma::vec>(niter);
  arma::vec outh = arma::zeros<arma::vec>(niter);
  arma::mat tar = arma::zeros<arma::mat>(3, niter);
  //initialize
  outbeta.col(0) = initialbeta;
  outgam.col(0) = initialgamma;
  outSigma.slice(0) = initialSigma;
  outsb(0) = initialsigmabeta;
  tar.col(0) = arma::zeros<arma::vec>(3);
  outh(0) = get_h_from_sigmabeta_c(X,outsb(0),outSigma.slice(0),
       outgam.col(0), n, T);
  for (int i=1; i<niter; ++i){
    arma::vec gam1    = outgam.col(i-1);
    arma::vec beta1   = outbeta.col(i-1);
    arma::mat Sigma1  = outSigma.slice(i-1);
    double sigmabeta1 = outsb(i-1);
    double h1         = outh(i-1);
    Rcpp::List bg = update_betagam_sw_c(X,Y,gam1,beta1,Sigma1,
                                        abs(marcor),sigmabeta1,Vbeta,bgiter,switer);
    arma::vec gam2  = as<arma::vec>(bg["gam"]);
    arma::vec beta2 = as<arma::vec>(bg["beta"]);
    arma::mat Sigma2 = update_Sigma_c(n,nu,X,beta2,Phi,Y);
    Rcpp::List hsig = update_h_c(h1, hiter, gam2, beta2, Sigma2, X, T);
    outh(i) = hsig["h"];
    outsb(i) = hsig["sigbeta"];
    outgam.col(i) = gam2;
    outbeta.col(i) = beta2;
    outSigma.slice(i) = Sigma2;
    if(!arma::is_finite(outsb(i))){
      outsb(i) = 1000;
    }
    tar.col(i) = get_target_c(X,Y,outsb(i), Sigma2,gam2, beta2);
    cout << i << "\n";
  }
  return Rcpp::List::create(
    Rcpp::Named("gam") = wrap(outgam.t()),
    Rcpp::Named("beta") = wrap(outbeta.t()),
    Rcpp::Named("sigbeta") = wrap(outsb),
    Rcpp::Named("Sigma") = wrap(outSigma),
    Rcpp::Named("tar") = wrap(tar.t())
  );
  
}


// [[Rcpp::export]]
Rcpp::List run2chains_c(arma::vec X,
                        arma::mat Y,
                        Rcpp::List initial_chain1,
                        Rcpp::List initial_chain2,
                        arma::mat Phi,
                        int niter = 1000,
                        int bgiter = 500,
                        int hiter = 50,
                        int switer = 50,
                        int burnin = 5){
  //initialize if not user-defined
  int T = Y.n_cols;
  int n = Y.n_rows;
  int nu = T+5;
  
  //marginal correlation
  arma::rowvec marcor = arma::zeros<arma::rowvec>(T);
  for (int t=0; t<T; ++t){
    arma::vec ycol = Y.col(t);
    arma::uvec non_null = find_finite(ycol);
    marcor(t) = sum(ycol(non_null)%X(non_null))/non_null.size();
  }
  //initialize Vbeta
  double Vbeta = sum(marcor%marcor) * 0.01;
  
  arma::mat outbeta1 = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam1 = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma1 = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb1 = arma::zeros<arma::vec>(niter);
  arma::vec outh1 = arma::zeros<arma::vec>(niter);
  arma::mat tar1 = arma::zeros<arma::mat>(3, niter);
  
  arma::mat outbeta2 = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam2 = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma2 = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb2 = arma::zeros<arma::vec>(niter);
  arma::vec outh2 = arma::zeros<arma::vec>(niter);
  arma::mat tar2 = arma::zeros<arma::mat>(3, niter);
  
  
  outbeta1.col(0)    = as<arma::vec>(initial_chain1["beta"]);
  outgam1.col(0)     = as<arma::vec>(initial_chain1["gamma"]);
  outSigma1.slice(0) = as<arma::mat>(initial_chain1["Sigma"]);
  outsb1(0)          = initial_chain1["sigmabeta"];
  
  outbeta2.col(0)    = as<arma::vec>(initial_chain2["beta"]);
  outgam2.col(0)     = as<arma::vec>(initial_chain2["gamma"]);
  outSigma2.slice(0) = as<arma::mat>(initial_chain2["Sigma"]);
  outsb2(0)          = initial_chain2["sigmabeta"];
  
  for (int i=1; i<niter; ++i){
    //chain 1 update
    Rcpp::List bg = update_betagam_sw_c(X,
                                        Y,
                                        outgam1.col(i-1),
                                        outbeta1.col(i-1),
                                        outSigma1.slice(i-1),
                                        abs(marcor),
                                        outsb1(i-1),
                                        Vbeta,
                                        bgiter,
                                        switer);
    outgam1.col(i)  = as<arma::vec>(bg["gam"]);
    outbeta1.col(i) = as<arma::vec>(bg["beta"]);
    outSigma1.slice(i) = update_Sigma_c(n,nu,X,outbeta1.col(i),Phi,Y);
    Rcpp::List hsig = update_h_c(outh1[i-1],
                                 hiter,
                                 outgam1.col(i),
                                 outbeta1.col(i),
                                 outSigma1.slice(i),
                                 X,
                                 T);
    outh1(i) = hsig["h"];
    outsb1(i) = hsig["sigbeta"];
    if(!arma::is_finite(outsb1(i))){
      outsb1(i) = 1000;
    }
    tar1.col(i) = get_target_c(X,
             Y,
             outsb1(i),
             outSigma1.slice(i),
             outgam1.col(i),
             outbeta1.col(i));
    
    //chain 2 update
    bg = update_betagam_sw_c(X,
                             Y,
                             outgam2.col(i-1),
                             outbeta2.col(i-1),
                             outSigma2.slice(i-1),
                             abs(marcor),
                             outsb2(i-1),
                             Vbeta,
                             bgiter,
                             switer);
    outgam2.col(i)  = as<arma::vec>(bg["gam"]);
    outbeta2.col(i) = as<arma::vec>(bg["beta"]);
    outSigma2.slice(i) = update_Sigma_c(n,nu,X,outbeta1.col(i),Phi,Y);
    hsig = update_h_c(outh1[i-1],
                      hiter,
                      outgam2.col(i),
                      outbeta2.col(i),
                      outSigma2.slice(i),
                      X,
                      T);
    outh2(i) = hsig["h"];
    outsb2(i) = hsig["sigbeta"];
    if(!arma::is_finite(outsb2(i))){
      outsb2(i) = 1000;
    }
    tar2.col(i) = get_target_c(X,
             Y,
             outsb2(i),
             outSigma2.slice(i),
             outgam2.col(i),
             outbeta2.col(i));
    
    //convergence criterion
    if(i>2*burnin && i%5==0){
      arma::vec rowmean1 = mean(outgam1.cols(burnin,i), 1);
      arma::vec rowmean2 = mean(outgam2.cols(burnin,i), 1);
      if(all(rowmean1<0.5) & all(rowmean2<0.5)){
        cout<< "both chains selected no variables - converged!";
        outgam1 = outgam1.cols(0,i);   outgam2 = outgam2.cols(0,i);
        outbeta1 = outbeta1.cols(0,i); outbeta2 = outbeta2.cols(0,i);
        outSigma1 = outSigma1.slices(0,i); outSigma2 = outSigma2.slices(0,i);
        outh1 = outh1(0,i); outh2 = outh2(0,i);
        outsb1 = outsb1(0,i); outsb2 = outsb2(0,i);
        break;
      }else{
        arma::uvec est1 = find(rowmean1 > 0.5);
        arma::uvec est2 = find(rowmean2 > 0.5);
        if(est1.size()==est2.size() && all(est1==est2)){
          arma::mat tmpbeta1 = outbeta1.cols(burnin,i);
          arma::mat tmpbeta2 = outbeta2.cols(burnin,i);
          arma::mat tmpgam1 = outgam1.cols(burnin,i);
          arma::mat tmpgam2 = outgam2.cols(burnin,i);
          double diff = 0;
          for (int k = 0; k < est1.size(); ++k){
            int kk = est1(k);
            arma::uvec ones = find(tmpgam1.row(kk)==1);
            arma::rowvec tmptmpbeta1 = tmpbeta1.row(kk);
            double beta1 = mean(tmptmpbeta1(ones));
            arma::uvec ones2 = find(tmpgam2.row(kk)==1);
            arma::rowvec tmptmpbeta2 = tmpbeta2.row(kk);
            double beta2 = mean(tmptmpbeta2(ones));
            diff = diff + (beta1-beta2)*(beta1-beta2);
          }
          if(diff/est1.size() < 1e-2){
            cout<< "beta difference is small between the two chains - converged!\n";
            outSigma1.shed_slices(i+1, niter-1);  outSigma2.shed_slices(i+1,niter-1);
            outgam1.shed_cols(i+1,niter-1);       outgam2.shed_cols(i+1,niter-1);
            outbeta1.shed_cols(i+1, niter-1);     outbeta2.shed_cols(i+1, niter-1);
            outh1.shed_rows(i+1,niter-1);         outh2.shed_rows(i+1, niter-1);
            outsb1.shed_rows(i+1,niter-1);        outsb2.shed_rows(i+1,niter-1);
            break;
          }
        }
      }
    }
    cout << i << "\n";
  }
  return Rcpp::List::create(
    Rcpp::Named("chain1") = Rcpp::List::create(
      Rcpp::Named("gamma") = outgam1.t(),
      Rcpp::Named("beta") = outbeta1.t(),
      Rcpp::Named("Sigma") = outSigma1,
      Rcpp::Named("sigmabeta") = outsb1,
      Rcpp::Named("h") = outh1
    ),
    Rcpp::Named("chain2") = Rcpp::List::create(
      Rcpp::Named("gamma") = outgam2.t(),
      Rcpp::Named("beta") = outbeta2.t(),
      Rcpp::Named("Sigma") = outSigma2,
      Rcpp::Named("sigmabeta") = outsb2,
      Rcpp::Named("h") = outh2
    )
  );
}







