/*
 * This file contains a class for the homogeneous AGamma process prior.
 * The process is defined on [0,1] and the stationary increments are
 * distributed according to the AGamma distribution, with
 * parameters eta, omega and Sigma.
 * 
 * Note: needed for vnp algorithm in gibbs_vnp_algorithm.cpp
 * 
 * For mathematical details on Hpd Gamma processes, 
 * and its numerical representation, see Chapter 3 in Meier (2018).
 * For details on the AGamma distribution and the AGamma process, see
 * Section 2.4 and Section 3.4.2 in Meier (2018).
 * 
 * A. Meier (2018). "A Matrix Gamma Process and Applications to Bayesian 
 * Analysis of Multivariate Time Series". PhD thesis, OvGU Magdeburg.
 */

#include <algorithm> // sort
#include "a_gamma_process_prior.h"
#include <boost/math/special_functions/expint.hpp>

// constructor
AGammaProcessPrior::AGammaProcessPrior(double eta, double omega, 
                                       arma::cx_mat Sigma) : eta(eta), 
                                       omega(omega), d(Sigma.n_cols), 
                                       dd((double)Sigma.n_cols), Sigma(Sigma),
                                       Sigma_inv(inv_sympd(Sigma)){}

// log alpha function -- see (3.26) in Meier (2018)
double AGammaProcessPrior::lalpha(const arma::cx_mat& U) const {
  return -dd * eta * std::log(beta(U)) + (eta-dd) * arma::log_det(U).real();
}

// beta function -- see Lemma 2.7 in Meier (2018)
double AGammaProcessPrior::beta(const arma::cx_mat& U) const {
  return arma::trace(Sigma_inv * U).real();
}

// Exponential integral function
double AGammaProcessPrior::e1(double x) {
  return -boost::math::expint(-x);
}

// tail of Gamma(a,b) measure
double AGammaProcessPrior::eab(double x, double a, double b) {
  return a * e1(b*x);
}

// log prior density, see Lemma 3.13 in Meier (2018)
double AGammaProcessPrior::lprior(const arma::vec& x, const arma::vec& r, 
                                  const arma::cx_cube& U) const {
  const int L = x.size();
  std::vector<double> w(L);
  double res=0.0;
  for (int j=0; j<L; ++j) {
    res += lalpha(U.slice(j)); // alpha(x,U) part
    const double beta_j = beta(U.slice(j));
    w[j] = eab(r(j), omega, beta_j);
    res -= (beta_j * r(j) + std::log(r(j))); // Jacobian
  }
  std::sort(w.begin(), w.end());
  res -= w[0];
  for (int j=1; j<L; ++j) {
    res -= (w[j]-w[j-1]); // radial parts
  }
  return res;
}
