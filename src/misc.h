#ifndef __MISC_INCLUDED__
#define __MISC_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericMatrix acvMatrix(Rcpp::NumericVector acv);
arma::mat acvBlockMatrix(arma::mat acv);
double acceptanceRate(Rcpp::NumericVector trace);

#endif
