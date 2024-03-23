#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include "distRcpp.h"
// #include "SSRFun_cpp.h"
// #include "dinvgamma_cpp.h"
// #include "dmvnrm_arma_fast.h"

using namespace Rcpp;
// [[Rcpp::export]]

// dproposal function
double dproposalFun_cpp(arma::mat dist_mat,
                        Rcpp::List para_result_l, Rcpp::List para_result_r,
                        String metric, Rcpp::List hyperparList){

  double n_obj = dist_mat.n_rows;
  int m = n_obj * (n_obj - 1) / 2;
  
  // extract hyperparameters
  double a = hyperparList["a"];
  double b = hyperparList["b"];
  double alpha = hyperparList["alpha"];
  arma::vec beta = hyperparList["beta"];
  double constant_multiple = hyperparList["constant_multiple"];
  
  // extract parameter values
  arma::mat x_l = para_result_l["x"];
  double sigma2_l = para_result_l["sigma2"];
  arma::mat lambda_l = para_result_l["lambda"];
  arma::mat x_r = para_result_r["x"];
  double sigma2_r = para_result_r["sigma2"];
  arma::mat lambda_r = para_result_r["lambda"];
  int p = x_l.n_cols;
  
  // some transformations
  arma::colvec lambda_l_diag = lambda_l.diag();
  Rcpp::NumericVector lambda_l_diag_rcpp = NumericVector(lambda_l_diag.begin(),
                                                    lambda_l_diag.end());

  // calculate delta matrix and d matrix
  arma::mat d_mat = dist_mat;
  Rcpp::NumericMatrix delta_mat_l_rcpp = distRcpp(wrap(x_l));
  arma::mat delta_mat_l = Rcpp::as<Mat<double>>(delta_mat_l_rcpp);
  Rcpp::NumericMatrix delta_mat_r_rcpp = distRcpp(wrap(x_r));
  arma::mat delta_mat_r = Rcpp::as<Mat<double>>(delta_mat_r_rcpp);

  // calculate SSR
  // double SSR_l = SSRFun_cpp(d_mat, delta_mat_l);
  double SSR_r = SSRFun_cpp(d_mat, delta_mat_r);

  // lambda_j
  // The full conditional posterior distrbuiton of lambda_j is the inverse Gamma distribution
  // lambda_j ~ IG(alpha + n/2 , beta_j + s_j/2)
  // where s_j/n is the sample variance of the jth coordinates of x_i's
  arma::vec lambda_density_mat_diag(p);
  for (int j = 0; j < p; j++){
    // calculate shape parameter
    double shape_lambda = alpha + n_obj / 2;
    // calculate scale parameter
    double sj_n = var(x_r.col(j));
    double scale_lambda = beta[j] + sj_n * n_obj / 2;
    // Rcpp::NumericVector lambda_l_rcpp = as<NumericVector>(wrap(lambda_l_diag[j]));
    // Rcout << "x: " << lambda_l_diag_rcpp;
    // Rcout << "shape_lambda: " << shape_lambda;
    // Rcout << "scale_lambda: " << 1.0/scale_lambda;
    Rcpp::NumericVector temp = dinvgamma_cpp(lambda_l_diag_rcpp,
                                             shape_lambda,
                                             1.0/scale_lambda, true);
    // Rcout << "temp: " << temp;
    arma::vec temp_arma = as<arma::vec>(wrap(temp));
    //temp_arma.print();
    lambda_density_mat_diag[j] = temp_arma[j];
  }
  
  double lambda_density = sum(lambda_density_mat_diag);
  // Rcout << "lambda_density" << lambda_density;
  
  // x_i
  // A normal proposal density is used in the random walk Metropolis algorithm for
  // generation of x_i, i = 1, ... , n
  // Choose the variance of the normal proposal density to be a constant
  // multiple of sigma2/(n-1)
  arma::vec x_density_vec_temp;
  arma::vec x_density_vec(n_obj);
  arma::vec x_var(p, arma::fill::value(constant_multiple * sigma2_r / (n_obj - 1)));
  arma::mat x_sigma = diagmat(x_var);
  for (int i = 0; i < n_obj; i++){
    x_density_vec_temp = dmvnrm_arma_fast(x_l.row(i), x_r.row(i),
                                          x_sigma, true);
    // x_density_vec_temp.print();
    x_density_vec[i] = as_scalar(x_density_vec_temp);
    // x_density_vec.print();
  }
  double x_density = sum(x_density_vec);
  
  // sigma2
  // A normal proposal density is used in the random walk Metropolis algorithm for
  // generation of sigma2
  // Choose the variance of the normal proposal density to be proportional to the
  // variance of IG(m/2+a, SSR/2+b)
  double sigma2_sd = sqrt(constant_multiple*pow((SSR_r/2+b),2)/(pow((m/2+a-1),2)*(m/2+a-2)));
  double sigma2_density = R::dnorm(log(sigma2_l), sigma2_r, sigma2_sd, true);

  // Rcout << "lambda_density: " << lambda_density;
  // Rcout << "x_density: " << x_density;
  // Rcout << "sigma2_density: " << sigma2_density;
  
  double output = lambda_density + x_density + sigma2_density;

  return(output);

}



/*** R

# dproposalFun_cpp(dist_mat,
#                  para_result_l = para.result.l, para_result_r = para.result.r,
#                  metric, hyperparList)



*/
