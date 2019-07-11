#ifndef __GIBBS_NP_HELP_INCLUDED__
#define __GIBBS_NP_HELP_INCLUDED__

#include <Rcpp.h>

Rcpp::NumericVector pFromV(Rcpp::NumericVector v);
Rcpp::NumericVector vFromP(Rcpp::NumericVector p, const double eps);
Rcpp::NumericVector mixtureWeight(Rcpp::NumericVector p, 
                                  Rcpp::NumericVector w, unsigned k);
Rcpp::NumericVector densityMixture(Rcpp::NumericVector weights, 
                                   Rcpp::NumericMatrix densities);
Rcpp::NumericVector unrollPsd(Rcpp::NumericVector qPsd, unsigned n);

#endif
