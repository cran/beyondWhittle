#include "gibbs_vnp_help.h"
#include "matrix_cube.h"

using namespace std;
using namespace Rcpp;


//' Store imaginary parts above and real parts below the diagonal
//' @keywords internal
// [[Rcpp::export]]
arma::cube realValuedPsd(ComplexVector f_) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  // convert to {f11, Re(f12), Im(f12), f22}-representation:
  arma::cube res(f.n_rows, f.n_cols, f.n_slices, arma::fill::zeros);
  for (unsigned j=0; j < f.n_slices; ++j) {
    for (unsigned r=0; r < f.n_rows; ++r) {
      for (unsigned s=0; s < f.n_cols; ++s) {
        if (r <= s) {
          res(r,s,j) = f(r,s,j).real();
        } else {
          res(r,s,j) = f(r,s,j).imag();
        }
      }
    }
  }
  return res;
}

//' Inverse function to realValuedPsd
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube complexValuedPsd(NumericVector f_) {
  const arma::cube f = cube_from_NumericVector(f_);
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices, arma::fill::zeros);
  for (unsigned j=0; j < f.n_slices; ++j) {
    for (unsigned r=0; r < f.n_rows; ++r) {
      for (unsigned s=0; s < f.n_cols; ++s) {
        if (r < s) {
          res(r,s,j) = arma::cx_double(f(r,s,j),-f(s,r,j));
        } else if (r==s) {
          res(r,s,j) = arma::cx_double(f(r,s,j),0.0);
        } else {
          res(r,s,j) = arma::cx_double(f(s,r,j),f(r,s,j));
        }
      }
    }
  }
  return res;
}

// log Whittle likelihood (multivariate), unnormalized
// Note: omega=0 and omega=pi excluded!
double llike_whittle(const arma::cx_mat& FZ, const arma::cx_cube& f) {
  const int N = FZ.n_rows;
  double res(0.0);
  for (unsigned j=1; j < N-1; ++j) {
    const arma::cx_mat Sigma(2.0 * M_PI * f.slice(j));
    const arma::cx_vec z(FZ.row(j).st()); // transpose without hermitian
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,Sigma);
    const double log_determ = log_det_val.real();
    const arma::cx_mat zSz = arma::trans(z) * arma::inv(Sigma) * z;
    res += (-log_determ - zSz(0,0).real());
  }
  return res;
}

// Log prior for Bernstein Gamma Psd (multivariate)
double lprior_bernsteinGammaPsd(const bernsteinGammaPsd& f, 
                                const AGammaProcessPrior& ap,
                                double k_theta) {
  const int k = f.get_k();
  const double klk = (double)k * std::log((double)k);
  return -k_theta * klk + ap.lprior(f.get_Z(), f.get_r(), f.get_U());
}

// Log posterior of Bernstein Gamma Psd under Whittle (multivariate)
double lpost_bernsteinGammaWhittle(const arma::cx_mat& FZ,
                                   const bernsteinGammaPsd& f,
                                   const AGammaProcessPrior& ap,
                                   double k_theta) {
  return lprior_bernsteinGammaPsd(f, ap, k_theta) + 
    llike_whittle(FZ, f.eval());
}

// Complex FFT, CPP-internal only
// Note: Functionality equivalent to mdft(..., real=F)
arma::cx_mat mdft_cpp(const arma::mat& x) {
  const int n = x.n_rows;
  const int d = x.n_cols;
  arma::mat x_t(n,d);
  x_t.row(0) = x.row(n-1);
  for (unsigned j=1; j<n; ++j) {
    x_t.row(j) = x.row(j-1);
  }
  const arma::cx_mat fourier(arma::fft(x_t) / std::sqrt((double)n));
  const int N((n%2) ? (n+1)/2 : n/2+1);
  return fourier.rows(0, N-1);
}

//' Construct psd mixture
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube get_f_matrix(arma::mat U_phi, arma::vec r, arma::vec Z,
                           int k, Rcpp::List dbList) {
  return bernsteinGammaPsd(U_phi, r, Z, &dbList, k).eval();
}

// (Sample) mean center observations
void mean_center(arma::mat& x) {
  const arma::rowvec mu(arma::mean(x));
  for (int j=0; j<x.n_rows; ++j) {
    x.row(j) -= mu;
  }
}
