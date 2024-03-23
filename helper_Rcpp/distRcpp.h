#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace Rcpp; 
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]

Rcpp::NumericMatrix distRcpp(Rcpp::NumericMatrix X){
  Rcpp::NumericMatrix XX(X);
  //Rcpp::Rcout << X << endl; 
  //Rcpp::Rcout << XX << endl;    
  int p = XX.nrow();
  Rcpp::NumericMatrix DIST(p,p);
  for(int i=0; i< p; i++){
    for(int j=0; j< p; j++){
      if(i!=j){
        // Euclidean distance
        DIST(i,j)=sqrt(sum((XX.row(i)-XX.row(j))*(XX.row(i)-XX.row(j))));
        // Cosine distance
        //DIST(i,j)= 1-sum(XX.row(i)*XX.row(j))/(sqrt(sum(pow(XX.row(i),2)))*sqrt(sum(pow(XX.row(j),2))));
      }
    }    
  }
  return DIST;
}   

/*** R

# set.seed(123)
# X <- rmvnorm(n=5, mean = rep(0, 3), sigma = diag(3))
# 
# distRcpp(X)
# 
# 1 - philentropy::distance(X, method = "cosine", mute.message = TRUE)

*/
