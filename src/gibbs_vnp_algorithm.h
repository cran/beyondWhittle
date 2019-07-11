#ifndef __GIBBS_VNP_ALGORITHM_INCLUDED__
#define __GIBBS_VNP_ALGORITHM_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List gibbs_multivariate_nuisance_cpp(arma::mat data,
                                           arma::cx_mat FZ,
                                           Rcpp::NumericVector eps_r,
                                           Rcpp::NumericVector eps_Z,
                                           Rcpp::NumericVector eps_U,
                                           int k_0,
                                           arma::vec r_0,
                                           arma::vec Z_0,
                                           arma::mat U_phi_0,
                                           arma::vec phi_def,
                                           double eta,
                                           double omega,
                                           arma::cx_mat Sigma,
                                           int Ntotal,
                                           int print_interval,
                                           double numerical_thresh,
                                           bool verbose,
                                           int L,
                                           double k_theta,
                                           Rcpp::List dbList);

#endif
