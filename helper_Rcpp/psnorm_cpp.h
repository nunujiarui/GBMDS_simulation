#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
SEXP psnorm_cpp(double q, double mean, double sd, double psi){
  
  // Obtaining namespace of Matrix package
  Environment pkg = Environment::namespace_env("fGarch");
  
  // Picking up Matrix() function from Matrix package
  Function f = pkg["psnorm"];
  
  // Executing Matrix( m, sparse = TRIE )
  //return f( q, psi, Named("mean", -2), Named("sd", 2));
  SEXP res = f(Rcpp::_["q"] = q, Rcpp::_["mean"] = mean, 
               Rcpp::_["sd"] = sd, Rcpp::_["xi"] = psi);
  return(res);
}



/*** R

# psnorm_cpp(0, -2, 2, 1.5)
# 0.8377841

*/
