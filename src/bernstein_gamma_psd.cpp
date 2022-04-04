/*
 * This file contains a class for a Bernstein-Hpd-Gamma spectral density
 * matrix function in dimension dxd (d>1), see (5.2) in Meier (2018).
 * 
 * Note: needed for vnp algorithm in gibbs_vnp_algorithm.cpp
 * 
 * It is parametrized like in (5.8) in Meier (2018) as follows:
 *   U_phi:  (d*d-1)xL-matrix, spherical coordinate representation 
 *           of unit trace matrices U_1,...,U_L
 *   r:      L-vector, radial parts of Gamma process
 *   Z:      L-vector of mass allocation points 
 *           (here in [0,1] instead of [0,pi])
 *   dbList: list of precomputed Bernstein polynomials
 *   k:      integer, Bernstein polynomial degree
 * 
 * A. Meier (2018). "A Matrix Gamma Process and Applications to Bayesian 
 * Analysis of Multivariate Time Series". PhD thesis, OvGU Magdeburg.
 */

#include "bernstein_gamma_psd.h"

// constructor
bernsteinGammaPsd::bernsteinGammaPsd(const arma::mat& U_phi, 
                                     const arma::vec& r, 
                                     const arma::vec& Z, 
                                     const Rcpp::List* dbList,
                                     int k)  : r(r), 
                                     Z(Z), dbList(dbList), k(k) {
  U = get_U_cpp(U_phi);
  W = get_W(U, r, Z, k);
  update_f();
}
  
// get mixture weights W_1,...,W_k
arma::cx_cube bernsteinGammaPsd::get_W(const arma::cx_cube& U,
                                       const arma::vec& r,
                                       const arma::vec& Z,
                                       int k) const {
  const int d=U.n_cols;
  const int L=U.n_slices;
  // w
  arma::cx_cube w(d,d,L); // careful: no fill
  for (int j=0; j<L; ++j) {
    w.slice(j) = r(j)*U.slice(j);
  }
  // W
  arma::cx_cube W_ = arma::cx_cube(d, d, k, arma::fill::zeros);
  for (int j = 0; j < k; ++j) {
    for (int l = 0; l < L; ++l) {
      if (((double)j/k < Z[l]) && (Z[l] <= (double)(j+1)/k)) {
        W_.slice(j) += w.slice(l);
      }
    }
  }
  return W_;
}

// Add W * b_{j,k} to f (b_{j,k} being the j'th basis polynomial of degree k).
// Local recomputation, as benefitial in MH-within-Gibbs-steps.
void bernsteinGammaPsd::update_f_by_increment(const arma::cx_mat& W_incr,
                           int j) {

  const arma::mat db = dbList->at(k-1);
  const int n = db.n_cols;
  for (int i=0; i < n; ++i) {
    f.slice(i) += W_incr * db(j,i);
  }
}

// Full recomputation of internal f according
// to current parameter values
void bernsteinGammaPsd::update_f() {
  const arma::mat db = dbList->at(k-1);
  const int n = db.n_cols;
  const int d = U.n_cols;
  f = arma::cx_cube(d, d, n, arma::fill::zeros); 
  for (int j=0; j < k; ++j) {
    for (int i=0; i < n; ++i) {
      f.slice(i) += W.slice(j) * db(j,i);
    }
  }
}

// Which k-segment does x in [0,1] belong to?
int bernsteinGammaPsd::get_j(double x) const {
  const int res = std::ceil((double)k * x)-1;
  // assert (res >= 0);
  return res;
}

// getter
int bernsteinGammaPsd::get_k() const {
  return k;
}

// getter
const arma::cx_cube& bernsteinGammaPsd::eval() const {
  return f;
}

// wiggle k, and keep internal representation consistent
void bernsteinGammaPsd::update_k(int new_k) {
  k = new_k;
  W = get_W(U, r, Z, k);
  update_f();
}

// wiggle r, and keep internal representation consistent
void bernsteinGammaPsd::update_r(int l, double r_l) {
  const int j=get_j(Z(l)); // index of affected W
  const arma::cx_mat W_incr = (r_l-r(l))*U.slice(l);
  W.slice(j) += W_incr;
  r(l)=r_l;
  update_f_by_increment(W_incr, j);
}

// wiggle U_phi, and keep internal representation consistent
void bernsteinGammaPsd::update_U_phi(int l, const arma::vec& U_phi_l) {
  const int j=get_j(Z(l)); // index of affected W
  const arma::cx_mat U_l = unit_trace_U_from_phi_cpp(U_phi_l);
  const arma::cx_mat W_incr = r(l)*(U_l-U.slice(l));
  W.slice(j) += W_incr;
  U.slice(l) = U_l;
  update_f_by_increment(W_incr, j);
}

// wiggle Z, and keep internal representation consistent
void bernsteinGammaPsd::update_Z(int l, double Z_l) {
  const int j_old = get_j(Z(l));
  const int j_new = get_j(Z_l);
  if (j_old != j_new) {
    const arma::cx_mat W_incr = r(l)*U.slice(l);
    W.slice(j_old)-=W_incr;
    W.slice(j_new)+=W_incr;
    update_f_by_increment(-W_incr, j_old);
    update_f_by_increment(W_incr, j_new);
  }
  Z(l)=Z_l;
}

// getter
const arma::cx_cube& bernsteinGammaPsd::get_U() const {
  return U;
}

// getter
const arma::vec& bernsteinGammaPsd::get_r() const {
  return r;
}

// getter
const arma::vec& bernsteinGammaPsd::get_Z() const {
  return Z;
}
