#include <cmath>
#include "matrix_cube.h"

using namespace Rcpp;

//' I/O: Only used within Rcpp
//' 
//' This workaround for parsing cubes was neccessary at development time
//' because there was a (presumable) bug in RcppArmadillo that sometimes
//' caused the parsing of arma::cx_cube objects to fail, such that the function
//' received an un-initialized object instead of the parsed one.
//' 
//' The workaround parses an Rcpp vector instead, and manually
//' copies the data in an arma::cx_cube object.
//' Besides being redundant, it also makes the code less readable and it is
//' hoped that this workaround can be removed in future revisions.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube cx_cube_from_ComplexVector(ComplexVector x) {
  const IntegerVector dim_x = x.attr("dim");
  arma::cx_vec x_vec(x);
  return arma::cx_cube(x_vec.begin(), dim_x[0], dim_x[1],
                       dim_x[2], true); // re-allocate memory
}

//' I/O: Only used within Rcpp
//' Note: Same workaround as \code{cx_cube_from_ComplexVector}
//' @keywords internal
// [[Rcpp::export]]
arma::cube cube_from_NumericVector(NumericVector x) {
  const IntegerVector dim_x = x.attr("dim");
  arma::vec x_vec(x);
  return arma::cube(x_vec.begin(), dim_x[0], dim_x[1],
                       dim_x[2], true); // re-allocate memory
}

//' Does a matrix have an eigenvalue smaller than 0?
//' @keywords internal
// [[Rcpp::export]]
bool hasEigenValueSmallerZero(arma::cx_mat A, double TOL=0.0) {
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, A);
  bool smallerZero = false;
  for (unsigned j=0; j < eigval.size(); ++j) {
    if (eigval(j).real() < TOL) {
      smallerZero = true;
    }
  }
  return smallerZero;
}
