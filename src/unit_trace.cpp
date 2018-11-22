/*
 * This file contains methods for Hpd matrices with unit trace.
 * It is based on the hyperspherical coordinates approach from 
 * Mittelbach et al (2012), see also Section 3.4.1 in Meier (2018).
 * 
 * Note: needed for vnp algorithm in gibbs_vnp_algorithm.cpp
 * 
 * M. Mittelbach, B. Matthiesen and E. A. Jorswieck (2012).
 * "Sampling uniformly from the set of positive definite matrices with 
 * trace constraint". IEEE Transactions on Signal Processing, 60(5):2167–2179.
 * 
 * A. Meier (2018). "A Matrix Gamma Process and Applications to Bayesian 
 * Analysis of Multivariate Time Series". PhD thesis, OvGU Magdeburg.
 */

#include "unit_trace.h"
#include <math.h>
using namespace std;
using namespace Rcpp;

// Get U from phi
// Note: redundant to unit_trace_U_from_phi (see unit_trace.R),
// may be merged later
arma::cx_mat unit_trace_U_from_phi_cpp(const arma::vec& phi) {
  // x from phi
  const unsigned N = phi.size();
  double phiProd = 1.0;
  arma::vec x(N+1); // careful: un-initialized
  for (unsigned j=0; j<N; ++j) {
    x(j) = std::cos(phi(j)) * phiProd;
    phiProd *= std::sin(phi(j));
  }
  x(N) = phiProd;
  // L from x
  const unsigned d = std::sqrt(x.size());
  arma::cx_mat L(d, d, arma::fill::zeros);
  unsigned k=0;
  for (int i=0; i<d; ++i) {
    for (int j=0; j<i; ++j) {
      L(i,j) = arma::cx_double(x(k), -x(k+1));
      k += 2;
    }
    L(i,i) = x(k);
    ++k;
  }
  return L * L.t();
}

//' Get U from phi, vectorized, cpp internal only
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube get_U_cpp(const arma::mat& u_phi) {
  const int phi_dim = u_phi.n_rows;
  const int L=u_phi.n_cols;
  const int d=std::sqrt(phi_dim+1);
  arma::cx_cube res(d,d,L);
  for (int j=0; j<L; ++j) {
    res.slice(j) = unit_trace_U_from_phi_cpp(u_phi.col(j));
  }
  return res;
}

//' Get x from phi, see (62) in Mittelbach et al.
//' @keywords internal
// [[Rcpp::export]]
NumericVector unit_trace_x_from_phi(NumericVector phi) {
  // length of phi: d^2 - 1
  // length of x:   d^2
  const unsigned N = phi.size();
  double phiProd = 1.0;
  NumericVector x(N+1);
  for (unsigned j=0; j<N; ++j) {
    x(j) = std::cos(phi(j)) * phiProd;
    phiProd *= std::sin(phi(j));
  }
  x(N) = phiProd;
  return x;
}

//' Get L (lower triangular Cholesky) from x
//' Called U^* in Mittelbach et al, see (60) there
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat unit_trace_L_from_x(arma::vec x) {
  const unsigned d = std::sqrt(x.size());
  arma::cx_mat res(d, d, arma::fill::zeros);
  unsigned k=0;
  for (int i=0; i<d; ++i) {
    for (int j=0; j<i; ++j) {
      res(i,j) = arma::cx_double(x(k), -x(k+1));
      k += 2;
    }
    res(i,i) = x(k);
    ++k;
  }
  return res;
}

// log determinant of Jacobian, see (3.22) in Meier (2018)
double unit_trace_jacobian_log_determinant(const arma::vec& phi) {
  const int N = phi.size();
  const unsigned n = std::sqrt(N+1); // actually d
  double res = 0.0;
  int i = 1; // take care to index phi by l-1
  for (int l=1; l <= N; ++l) {
    if (l == i*i) {
      // cosine part
      const int p_l = 2*(n-i)+1;
      res += (double)p_l * std::log(std::abs(std::cos(phi[l-1])));
      ++i;
    }
    // sine part
    const int i_2 = i-1; // take care not to use i (but i_2) for sine part
    const int m = l - (i_2)*(i_2);
    const int kappa_l = n-i_2-1;
    const int lambda_l = (i_2-1)*n + 1 + m;
    const int q_l = n*n + kappa_l*n - lambda_l;
    res += (double)q_l * std::log(std::abs(std::sin(phi[l-1])));
  }
  return res;
}

//' Get p vector, see (67) in Mittelbach et al.
//' @keywords internal
// [[Rcpp::export]]
NumericVector unit_trace_p(unsigned d) {
  const unsigned N=d*d-1;
  NumericVector res(N);
  int i=1;
  for (unsigned l=1; l<=N; ++l) {
    if (l==i*i) {
      res[l-1] = 2*(d-i)+1;
      ++i;
    } else {
      res[l-1] = 0;
    }
  }
  return res;
}

//' Get q vector, see (68) in Mittelbach et al.
//' @keywords internal
// [[Rcpp::export]]
NumericVector unit_trace_q(unsigned d) {
  const unsigned N=d*d-1;
  NumericVector res(N);
  int i=1;
  for (unsigned l=1; l<=N; ++l) {
    if (l == i*i) {
      ++i;
    }
    const int i_2 = i-1; // take care not to use i (but i_2) for sine part
    const int m = l - (i_2)*(i_2);
    const int kappa_l = d-i_2-1;
    const int lambda_l = (i_2-1)*d + 1 + m;
    res[l-1] = d*d + kappa_l*d - lambda_l;
  }
  return res;
}
