/*
 * Several helping functions for the np algorithm, as implemented
 * in gibbs_np.R
 */

#include "gibbs_np_help.h"
using namespace Rcpp;

//' Get  p from v in Stick Breaking DP representation
//' @keywords internal
// [[Rcpp::export]]
NumericVector pFromV(NumericVector v) {
  unsigned L = v.size();
  NumericVector p(L + 1);
  double currentProduct = 1.0;
  double pSum = 0.0;
  for (unsigned l = 0; l < L; ++l) {
    p[l + 1] = currentProduct * v[l];
    currentProduct *= (1.0 - v[l]);
    pSum += p[l + 1];
  }
  p[0] = std::max(1.0 - pSum, 0.0); // account for numerical instabilities
  return p;
}

//' Get v from p (DP inverse stick breaking)
//' Note: p is assumed to have length L, i.e. it does NOT contain p_0
//' @keywords internal
// [[Rcpp::export]]
NumericVector vFromP(NumericVector p, const double eps=1e-8) {
  unsigned L = p.size();
  NumericVector v(L);
  double currentProduct = 1.0;
  for (unsigned l = 0; l < L; ++l) {
    v[l] = std::min(std::max(p[l] / currentProduct, eps),1.0-eps); // numerical stability
    //v[l] = p[l] / currentProduct;
    currentProduct *= (1.0 - v[l]);
  }
  return v;
}

//' Get mixture weights of Bernstein-Dirchlet-Mixtures 
//' @keywords internal
// [[Rcpp::export]]
NumericVector mixtureWeight(NumericVector p, NumericVector w, unsigned k) {
  typedef std::pair<double, double> wpType;
  std::vector<wpType> wp;
  for (unsigned l = 0; l < p.size(); ++l) {
    wp.push_back(wpType(w[l], p[l]));
  }
  std::sort(wp.begin(), wp.end());
  NumericVector weight(k);
  unsigned l = 0;
  for (unsigned j = 1; j <= k; ++j) {
    weight[j-1] = 0;
    double wMax = j / (double)k;
    while (l < wp.size() && wp[l].first <= wMax) {
      weight[j-1] += wp[l].second;
      l += 1;
    }
  }
  return weight;
}

//' Construct a density mixture from mixture weights and density functions.
//' @keywords internal
// [[Rcpp::export]]
NumericVector densityMixture(NumericVector weights, NumericMatrix densities) {
  if (weights.size() != densities.nrow()) {
    return(NumericVector());
  }
  const unsigned n = densities.ncol();
  NumericVector res(n);
  for (unsigned omega = 0; omega < n; ++omega) {
    res[omega] = 0.0;
  }
  for (unsigned j = 0; j < weights.size(); ++j) {
    for (unsigned omega = 0; omega < n; ++omega) {
      res[omega] += weights[j] * densities(j, omega);
    }
  }
  return(res);
}

//' Redundantly roll out a PSD from length N=floor(n/2) to length n
//' @keywords internal
// [[Rcpp::export]]
NumericVector unrollPsd(NumericVector qPsd, unsigned n) {
  NumericVector q(n);
  q[0] = qPsd[0];
  const unsigned N = (n-1)/2;
  for (unsigned i = 1; i <= N; ++i) {
    const unsigned j = 2 * i - 1;
    q[j] = qPsd[i];
    q[j+1] = qPsd[i];
  }
  if (!(n % 2)) {
    q[n-1] = qPsd[qPsd.size() - 1];
  }
  return(q);
}
