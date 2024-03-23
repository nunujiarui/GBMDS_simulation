#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
// Calculate the SSR value

double SSRFun_cpp(arma::mat d_mat, arma::mat delta_mat){
  
  d_mat.diag().fill(0);
  delta_mat.diag().fill(0);
  
  arma::mat A = pow(d_mat - delta_mat, 2);
  
  double SSR = 0.5*(accu(A));
  
  return SSR;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# SSRFun <- function(d.mat, delta.mat){
#   
#   # d.mat[!lower.tri(d.mat)] <- 0
#   # delta.mat[!lower.tri(delta.mat)] <- 0
#   diag(d.mat) <- 0
#   diag(delta.mat) <- 0
#   
#   ssr <- 0.5*sum((d.mat-delta.mat)^2)
#   
#   return(SSR = ssr)
#   
# }
# 
# erin = matrix(data = c(2,4,5,4,3,1,5,1,6), nrow = 3)
# #erin
# 
# rex <- matrix(data = c(3,2,1,2,5,4,1,4,2), nrow = 3)
# #rex
# 
# SSRFun(d.mat = erin, delta.mat = rex)
# SSRFun_cpp(d_mat = erin, delta_mat = rex)
*/
