#ifndef __MATRIX_CUBE_INCLUDED__
#define __MATRIX_CUBE_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::cx_cube cx_cube_from_ComplexVector(Rcpp::ComplexVector x);
arma::cube cube_from_NumericVector(Rcpp::NumericVector x);
bool hasEigenValueSmallerZero(arma::cx_mat A, double TOL);

#endif
