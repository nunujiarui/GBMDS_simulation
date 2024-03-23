#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include "dmvnrm_arma_fast.h"

using namespace Rcpp;
// using namespace arma;

// [[Rcpp::export]]

// Calculate the density of x
arma::vec reference_d_x_cpp(arma::mat x, arma::mat prev_x){
  
  int n = x.n_rows;
  int p = x.n_cols;
 
  arma::vec v(p, arma::fill::value(0.01));
  arma::mat D_mat = diagmat(v);
  
  arma::vec d_x_vec;
  arma::vec d_x(n);

  for (int i = 0; i < n; i++){
  
    arma::mat x_mat = x.row(i);
    //arma::rowvec prev_x_vec = prev_x.row(i);
    
    d_x_vec = dmvnrm_arma_fast(x_mat, prev_x.row(i), D_mat, true);
    //d_x_vec.print("d_x_vec:");
    
    d_x[i] = as_scalar(d_x_vec);
    //d_x.print("d_x:");
  }
  
  // d_x = vectorise(d_x);
  // D_mat.print("sd:");
  
  return(d_x);
  
}


/*** R

reference.d.x <- function(x, prev.x){
  
  n <- nrow(x)
  p <- ncol(x)
  reference.x.sd <- diag(rep(0.01, p))
  x.list <- lapply(seq_len(nrow(x)), function(i) x[i,])
  prev.x.list <- lapply(seq_len(nrow(prev.x)), function(i) prev.x[i,])
  
  reference.x.sd.list <- rep(list(reference.x.sd), n)
  
  d.x <- mapply(mvtnorm::dmvnorm, x.list, prev.x.list, reference.x.sd.list, 
                log = TRUE)
  
  return(d.x)
  
}

# reference.d.x(x = matrix(data = c(1.5,0.3,-2,1.1), nrow = 2),
#               prev.x = matrix(data = c(1.2,-0.1,-1.3,0.9), nrow = 2))
# # -26.232707  -7.232707
# 
# reference.d.x(x = matrix(data = c(1.5,0.3,-0.8,-2,1.1, 0.5), nrow = 3),
#               prev.x = matrix(data = c(1.2,-0.1,-0.5,-1.3,0.9, 0.7), nrow = 3))
# 
# #  -26.232707  -7.232707  -3.732707

reference_d_x_cpp(x = matrix(data = c(1.5,0.3,-2,1.1), nrow = 2),
              prev_x = matrix(data = c(1.2,-0.1,-1.3,0.9), nrow = 2))
# -26.232707  -7.232707

test = reference_d_x_cpp(x = matrix(data = c(1.5,0.3,-0.8,-2,1.1, 0.5), nrow = 3),
                  prev_x = matrix(data = c(1.2,-0.1,-0.5,-1.3,0.9, 0.7), nrow = 3))
test
class(test)

*/
