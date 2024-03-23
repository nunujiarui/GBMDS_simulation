#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector dinvgamma_cpp(Rcpp::NumericVector x,
                              double shape,
                              double rate,
                              bool log) {
  Rcpp::NumericVector log_f = Rcpp::dgamma(1.0/x, shape, rate, true) - 2 * Rcpp::log(x);
  if (log) {
    return log_f;
  } else{
    return Rcpp::exp(log_f);
  }
}

/*** R

# dinvgamma_cpp(x = 329.752, #diag(para.result.l$lambda)[1],
#           shape = 2.5, rate = 0.000248196, #1/4029.07, 
#           log = TRUE)
# 
# dinvgamma_cpp(x = diag(para.result.l$lambda)[2],
#               shape = 2.5, rate = 1/726.039, log = TRUE)

*/
