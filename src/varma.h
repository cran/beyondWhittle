#ifndef __VARMA_INCLUDED__
#define __VARMA_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat var_forecast(arma::mat x, arma::mat ar, arma::mat noise);
arma::cx_cube varma_transfer2psd(Rcpp::ComplexVector transfer_ar_,
                                 Rcpp::ComplexVector transfer_ma_,
                                 arma::cx_mat sigma);
arma::cx_cube transfer_polynomial(Rcpp::NumericVector lambda, arma::mat coef);
arma::mat epsilon_var(arma::mat zt, arma::mat ar);
double sldmvnorm(arma::mat z_t, arma::mat Sigma);

#endif
