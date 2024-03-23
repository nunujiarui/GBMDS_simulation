#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "reference_d_x_cpp.h"
using namespace Rcpp;
// [[Rcpp::export]]

// Calculate the reference ratio in log scale

double logReferenceRatio_cpp(arma::mat new_x, arma::mat current_x,
                             arma::mat prev_x){
  
  double sum1 = sum(reference_d_x_cpp(new_x, prev_x));
  double sum2 = sum(reference_d_x_cpp(current_x, prev_x));
  
  double ratio = sum1 - sum2;
  
  return (ratio);
  
}
