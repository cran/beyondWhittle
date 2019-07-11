#ifndef __GIBBS_VNP_HELP_INCLUDED__
#define __GIBBS_VNP_HELP_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "a_gamma_process_prior.h"
#include "bernstein_gamma_psd.h"

arma::cube realValuedPsd(Rcpp::ComplexVector f_);
arma::cx_cube complexValuedPsd(Rcpp::NumericVector f_);

double llike_whittle(const arma::cx_mat& FZ, const arma::cx_cube& f);
double lprior_bernsteinGammaPsd(const bernsteinGammaPsd& f, 
                                const AGammaProcessPrior& ap,
                                double k_theta);
double lpost_bernsteinGammaWhittle(const arma::cx_mat& FZ,
                                   const bernsteinGammaPsd& f,
                                   const AGammaProcessPrior& ap,
                                   double k_theta);
arma::cx_mat mdft_cpp(const arma::mat& x);

arma::cx_cube get_f_matrix(arma::mat U_phi, arma::vec r, arma::vec Z,
                           int k, Rcpp::List dbList);
void mean_center(arma::mat& x);

#endif
