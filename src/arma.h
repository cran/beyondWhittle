#ifndef __ARMA_INCLUDED__
#define __ARMA_INCLUDED__

#include <Rcpp.h>

Rcpp::NumericMatrix pacf2AR(Rcpp::NumericVector pacf);
Rcpp::NumericVector genEpsARMAC(Rcpp::NumericVector zt, Rcpp::NumericVector ar, 
                                Rcpp::NumericVector ma);
Rcpp::NumericVector psd_arma(Rcpp::NumericVector freq, Rcpp::NumericVector ar, 
                             Rcpp::NumericVector ma, double sigma2);

#endif
