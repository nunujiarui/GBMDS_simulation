cppFunction("
NumericVector rnormRcpp(int n, double mean, double sd) {
  NumericVector result = Rcpp::rnorm(n, mean, sd);
  return result;
}")

cppFunction("
NumericVector rinvgammaRcpp(int n, double shape, double scale) {
  NumericVector result = 1.0/Rcpp::rgamma(n, shape, 1.0/scale);
  return result;
}")

