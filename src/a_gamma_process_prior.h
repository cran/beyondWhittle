#ifndef __AGAMMA_PROCESS_PRIOR_INCLUDED__
#define __AGAMMA_PROCESS_PRIOR_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>


class AGammaProcessPrior {
  
private:
  
  const double eta;
  const double omega;
  const unsigned d;
  const double dd;
  const arma::cx_mat Sigma;
  const arma::cx_mat Sigma_inv;
  double lalpha(const arma::cx_mat& U) const;
  double beta(const arma::cx_mat& U) const;
  static double e1(double x);
  static double eab(double x, double a, double b);
  
public:
  
  AGammaProcessPrior(double eta, double omega, arma::cx_mat Sigma);
  double lprior(const arma::vec& x, const arma::vec& r, 
                const arma::cx_cube& U) const;
  
};

#endif
