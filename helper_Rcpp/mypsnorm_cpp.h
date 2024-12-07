#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

double mypsnorm_cpp(double q, double mean, double sd, double psi){
  
  // Standardize:
  double q_t = (q-mean)/sd;
  double pi_val = atan2(0, -1);
  double m1 = 2/sqrt(2*pi_val);
  double mu = m1 * (psi - 1/psi);
  double sigma = sqrt((1-pow(m1,2))*(pow(psi,2)+1/pow(psi,2)) + 
                      2*pow(m1,2) - 1);
  double z = q_t*sigma + mu;
  
  // Compute
  double sig = (z >= 0) ? 1 : -1;
  double Psi = pow(psi, sig);
  double g = 2  / (psi + 1/psi);
  
  double temp = (z >= 0) ? 1 : 0;
  double Probability = temp - sig * g * Psi * R::pnorm(-abs(z)/Psi,0,1,true,false);
  
  return(Probability);
}

/*** R

# psnorm(0, -2,2, 0) # 0.8377841

*/
