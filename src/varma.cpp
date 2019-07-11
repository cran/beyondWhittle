#include <cmath>
#include "varma.h"
#include "matrix_cube.h"

using namespace Rcpp;
using namespace std;

//' Get VARMA PSD from transfer polynomials 
//' Helping function for \code{psd_varma}
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube varma_transfer2psd(ComplexVector transfer_ar_,
                                 ComplexVector transfer_ma_,
                                 arma::cx_mat sigma) {
  const arma::cx_cube transfer_ar = cx_cube_from_ComplexVector(transfer_ar_);
  const arma::cx_cube transfer_ma = cx_cube_from_ComplexVector(transfer_ma_);
  const unsigned d = transfer_ar.n_rows;
  const unsigned N = transfer_ar.n_slices;
  const double pi = std::acos(-1.0);
  arma::cx_cube res(d,d,N);
  for (unsigned j=0; j<N; ++j) {
    const arma::cx_mat transfer_ar_inv = arma::inv(transfer_ar.slice(j));
    res.slice(j) = transfer_ar_inv * transfer_ma.slice(j) * sigma *
      transfer_ma.slice(j).t() * transfer_ar_inv.t() / 2.0 / pi;
  }
  return res;
}

//' VARMA transfer polynomials
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube transfer_polynomial(NumericVector lambda, arma::mat coef) {
  const unsigned d = coef.n_rows;
  const unsigned p = coef.n_cols / d;
  const unsigned N = lambda.size();
  const arma::cx_mat eye(d,d,arma::fill::eye);
  arma::cx_cube res(d,d,N);
  for (unsigned l=0; l < N; ++l) {
    res.slice(l) = eye;
    for (unsigned j=0; j < p; ++j) {
      res.slice(l) += coef.submat(0,j*d,d-1,(j+1)*d-1) *
        std::polar<double>(1.0, -lambda[l]*(double)(j+1));
    }
  }
  return res;
}

//' epsilon process (residuals) of VAR model
//' @keywords internal
// [[Rcpp::export]]
arma::mat epsilon_var(arma::mat zt, arma::mat ar) {
  const unsigned d = zt.n_cols;
  const unsigned n = zt.n_rows;
  const unsigned p = ar.n_cols / d;
  arma::mat res(n-p, d, arma::fill::zeros);
  for (unsigned t=p; t < n; ++t) {
    res.row(t-p) = zt.row(t);
    for (unsigned tt=1; tt<=p; ++tt) {
      res.row(t-p) -= zt.row(t-tt) *
        ar.submat(0,(tt-1)*d,d-1,tt*d-1).t();
    }
  }
  return res;
}

//' sum of multivariate normal log densities
//' with mean 0 and covariance Sigma, unnormalized
//' @keywords internal
// [[Rcpp::export]]
double sldmvnorm(arma::mat z_t, arma::mat Sigma) {
  double res(0.0);
  double log_det_val;
  double log_det_sign;
  arma::log_det(log_det_val,log_det_sign,Sigma);
  const arma::mat Sigma_inv = arma::inv(Sigma);
  for (unsigned j=0; j < z_t.n_rows; ++j) {
    const arma::vec z(z_t.row(j).st());
    arma::mat zSz = trans(z) * Sigma_inv * z;
    res += 0.5 * (-log_det_val - zSz(0,0));
  }
  return(res);
}
