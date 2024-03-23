#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 The three functions below, cholesky_tridiagonal, forward_algorithm, 
 and backward_algorithm were copied from the source code of package 
 stochvol by Georg Kastner that was shared under licence GPL (>= 2). 
 They were modified by the author of the current note on 30 August 2021
 and are shared under the GPL (>= 2) licence as well.
 */
List cholesky_tridiagonal(const vec& omega_diag, const double& omega_offdiag) {
  const int T = omega_diag.n_elem - 1;
  vec chol_diag(T+1);
  vec chol_offdiag(T+1);
  chol_diag[0] = std::sqrt(omega_diag[0]);
  for (int j = 1; j < T+1; j++) {
    chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
    chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
  }
  return List::create(_["chol_diag"]=chol_diag, _["chol_offdiag"]=chol_offdiag);
}

vec forward_algorithm(const vec& chol_diag, const vec& chol_offdiag, const vec& covector) {
  const int T = chol_diag.n_elem - 1;
  vec htmp(T+1);
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < T+1; j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
  return htmp;
}

vec backward_algorithm(const vec& chol_diag, const vec& chol_offdiag, const vec& htmp) {
  const int T = chol_diag.size() - 1;
  vec h(T+1);
  h[T] = htmp[T] / chol_diag[T];
  for (int j = T-1; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
  }
  return h;
}

// [[Rcpp::export]]
mat rmvnorm_arma_stochvol(int n, mat precision, vec location){
  
  int     T                 = precision.n_rows;
  vec     precision_diag    = precision.diag();
  double  precision_offdiag = precision(1,0);
  
  List  precision_chol  = cholesky_tridiagonal(precision_diag, precision_offdiag);
  vec   aa              = forward_algorithm(precision_chol["chol_diag"],
                                            precision_chol["chol_offdiag"],
                                                          location);
  mat draw(T, n);
  vec epsilon;
  for (int i=0; i<n; i++){
    epsilon     = rnorm(T);
    draw.col(i) = backward_algorithm(precision_chol["chol_diag"], 
             precision_chol["chol_offdiag"],
                           aa + epsilon);
  }
  
  return draw;
}


/*** R

# library(mgcv)
# 
# set.seed(12345)
# n           = 10
# T           = 250
# s           = rgamma(1, shape=10, scale=10)
# precision   = rgamma(1, shape=10, scale=10) * diag(T) + 2 * s * diag(T)
# sdiag(precision, 1)  = -s
# sdiag(precision, -1) = -s
# location    = as.matrix(rnorm(T))
# 
# covariance  = solve(precision)
# 
# rmvnorm_r = function(n, precision, location){
#   covariance = solve(precision)
#   draw = mvtnorm::rmvnorm(n, mean = location, sigma = covariance)
#   #draw        = backsolve(t(precision.L), matrix(rep(a, n), ncol=n) + epsilon)
#   
#   return(t(draw))
# }
# rmvnorm_r(1, precision, location) %>% View()
# 
# rmvnorm_bandchol = function(n, precision, location){
#   T           = dim(precision)[1]
#   precision.L = t(bandchol(precision))
#   
#   epsilon     = matrix(rnorm(n * T), ncol=n)
#   a           = forwardsolve(precision.L, location)
#   draw        = backsolve(t(precision.L), matrix(rep(a, n), ncol=n) + epsilon)
#   
#   return(draw)
# }
# 
# rmvnorm_bandchol(1, precision, location) %>% View()
# rmvnorm_arma_stochvol(1, precision, location) %>% View()
# 
# library(microbenchmark)
# 
# microbenchmark(cpp.stochvol = rmvnorm_arma_stochvol(n, precision, location),
#                r.solve1     = rmvnorm_r(n, precision, location),
#                r.solve2     = rmvnorm_bandchol(n, precision, location),
#                #check  = "equal",
#                setup  = set.seed(123)
# )


*/
