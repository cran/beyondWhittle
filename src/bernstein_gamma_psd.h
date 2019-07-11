#ifndef __BERNSTEIN_GAMMA_PSD_INCLUDED__
#define __BERNSTEIN_GAMMA_PSD_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "unit_trace.h"

//' multivariate psd class
class bernsteinGammaPsd {
  
private:
  
  arma::cx_cube f;
  arma::cx_cube W;
  arma::cx_cube U;
  arma::vec r;
  arma::vec Z;
  const Rcpp::List* dbList; 
  int k;
  
  arma::cx_cube get_W(const arma::cx_cube& U, const arma::vec& r,
                      const arma::vec& Z, int k) const;
  void update_f_by_increment(const arma::cx_mat& W_incr, int j);
  void update_f();
  int get_j(double x) const;
  
public:
  
  int get_k() const;
  const arma::cx_cube& get_U() const;
  const arma::vec& get_r() const;
  const arma::vec& get_Z() const;
  bernsteinGammaPsd(const arma::mat& U_phi, const arma::vec& r, 
                    const arma::vec& Z, const Rcpp::List* dbList,
                    int k);
  const arma::cx_cube& eval() const;
  void update_k(int new_k);
  void update_r(int l, double r_l);
  void update_U_phi(int l, const arma::vec& U_phi_l);
  void update_Z(int l, double Z_l);
};

#endif
