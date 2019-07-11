#include "misc.h"
#include "matrix_cube.h"
using namespace Rcpp;
using namespace std;

//' Build an n times n Toeplitz matrix from the 
//' autocovariance values gamma(0),...,gamma(n-1)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix acvMatrix(NumericVector acv) {
  unsigned p = acv.size();
  NumericMatrix m(p, p);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m(i,j) = acv[index];
    }
  }
  return(m);
}

//' Build an nd times nd Block Toeplitz matrix from the
//' (d times d) autocovariances gamma(0),...,gamma(n-1)
//' @keywords internal
// [[Rcpp::export]]
arma::mat acvBlockMatrix(arma::mat acv) {
  const unsigned d = acv.n_rows;
  const unsigned p = acv.n_cols / d;
  arma::mat m(p*d, p*d);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m.submat(i*d,j*d,(i+1)*d-1,(j+1)*d-1) = acv.submat(0,
               index*d, d-1, (index+1)*d-1);
    }
  }
  return(m);
}

//' Computing acceptance rate based on trace
//' Note: Only use for traces from continous distributions!
//' @keywords internal
// [[Rcpp::export]]
double acceptanceRate(NumericVector trace) {
  unsigned rejections = 0;
  for (unsigned i=1; i < trace.length(); ++i) {
    rejections += (trace[i]==trace[i-1]);
  }
  double rejectionRate = (double)rejections / (double)trace.length();
  return 1 - rejectionRate;
}
