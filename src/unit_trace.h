#ifndef __UNIT_TRACE_INCLUDED__
#define __UNIT_TRACE_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::cx_mat unit_trace_U_from_phi_cpp(const arma::vec& phi);
arma::cx_cube get_U_cpp(const arma::mat& u_phi);
Rcpp::NumericVector unit_trace_x_from_phi(Rcpp::NumericVector phi);
arma::cx_mat unit_trace_L_from_x(arma::vec x);
double unit_trace_jacobian_log_determinant(const arma::vec& phi);
Rcpp::NumericVector unit_trace_p(unsigned d);
Rcpp::NumericVector unit_trace_q(unsigned d);

#endif
