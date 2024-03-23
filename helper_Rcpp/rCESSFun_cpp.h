#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]

// Calculate relative conditional effective sample size
double rCESSFun_cpp(arma::vec W, arma::vec logL, double num, double phi){
  
  arma::vec logw = num * logL;
  // logw.print("logw:");
  double logmax = max(logw);
  arma::vec w = exp(logw - logmax);
  
  double rt = pow((sum(w.t()*W)), 2) / (sum(W.t()*pow(w, 2))) - phi;
  
  return rt;
}
