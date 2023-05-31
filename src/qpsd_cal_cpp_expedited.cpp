/*
 * This file contains a function for dealing with the matrices in
 * evaluating the time-varying spectral density function.
 */

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace arma ;
using namespace Rcpp ;

// [[Rcpp::export]]
arma::rowvec qpsd_cal_cpp_expedited(arma::mat& basis1,
                                    arma::mat& basis2,
                                    arma::colvec& p,
                                    arma::uvec& pexpend,
                                    arma::uvec& selector1cpp,
                                    arma::uvec& selector2cpp) {
  return(arma::sum((p.cols(pexpend) % basis1.rows(selector1cpp)) % basis2.rows(selector2cpp), 0));
  }
