/*
 * Metropolis-within-Gibbs sampler for the spectral density matrix
 * under Whittle's likelihood and a Bernstein-Hpd-Gamma prior.
 * The algorithm is elaborated in Section 5.2.2 in Meier (2018).
 * 
 * The C++ core is called from the gibbs_vnp.R file, see also
 * the documentation of the gibbs_vnp method for further details.
 * 
 * A. Meier (2018). "A Matrix Gamma Process and Applications to Bayesian 
 * Analysis of Multivariate Time Series". PhD thesis, OvGU Magdeburg.
 */

#include <math.h> // M_PI
#include "gibbs_vnp_algorithm.h"
#include "gibbs_vnp_help.h"
#include "beyondWhittle_types.h"
#include "unit_trace.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

//' Gibbs sampler in Cpp
//' @keywords internal
// [[Rcpp::export]]
List gibbs_multivariate_nuisance_cpp(arma::mat data,         // data
                                     IntegerVector NA_pos,   // position of missing values
                                     arma::cx_mat FZ,        // Fourier transformed data
                                     NumericVector eps_r, // proposal variances (fixed here)
                                     NumericVector eps_Z, //  "
                                     NumericVector eps_U, //  "
                                     int k_0,          // starting values
                                     arma::vec r_0,    //  "
                                     arma::vec Z_0,    //  "
                                     arma::mat U_phi_0,//  "
                                     arma::vec phi_def,  // domain of definition of phi_U's
                                     double eta,         // AGamma process parameter
                                     double omega,       //  "
                                     arma::cx_mat Sigma, //  "
                                     int Ntotal,              // MCMC parameters
                                     int print_interval,      //  "
                                     double numerical_thresh, //  "
                                     bool verbose,            //  "
                                     int L,           // AGamma series truncation
                                     double k_theta,  // prior parameter for k
                                     List dbList) {   // precomputed Bernstein polynomials
  
  const double SQRT_PI = std::sqrt(M_PI);
  const int NUM_NA = NA_pos.size();
  
  const int n=data.n_rows;
  const int d=data.n_cols;
  const int kmax=dbList.length();
  const int phi_dim=d*d-1;
  
  // initialize storage
  IntegerVector k(Ntotal);
  arma::mat r(L, Ntotal);
  arma::mat Z(L, Ntotal);
  arma::cube U_phi(phi_dim, L, Ntotal);
  NumericVector lpostTrace(Ntotal);
  arma::mat missingValuesTrace(Ntotal, NUM_NA);
  double f_store = -1.0;
  int num_except = 0;

  // starting values
  k(0) = k_0;
  Z.col(0) = Z_0;
  r.col(0) = r_0;
  U_phi.slice(0) = U_phi_0;
  bernsteinGammaPsd psd_current(U_phi.slice(0), 
                                r.col(0), 
                                Z.col(0), 
                                &dbList, 
                                k(0));
  const AGammaProcessPrior a_prior(eta, omega, Sigma);
  f_store = lpost_bernsteinGammaWhittle(FZ, psd_current, a_prior, k_theta);
  lpostTrace(0) = f_store;
  for (int j=0; j<NUM_NA; ++j) {
    missingValuesTrace(0,j) = data[NA_pos[j]];
  }
  
  // MH Loop
  for (int i=0; i<Ntotal-1; ++i) {
    
    // Step 1: draw k
    const int k_old = k(i);
    int k_star = k_old + std::round(Rcpp::rt(1,1)[0]);
    while (k_star < 3 || k_star > kmax) {
      k_star = k_old + std::round(Rcpp::rt(1,1)[0]);
    }
    bernsteinGammaPsd psd_k_star(psd_current);
    psd_k_star.update_k(k_star);
    // MH-step
    bool k_accepted;
    const double f_k = f_store;
    double f_k_star;
    try {
      f_k_star = lpost_bernsteinGammaWhittle(FZ, psd_k_star, a_prior, k_theta);   
      const double alpha_k = std::min(0.0, f_k_star - f_k);
      k_accepted = (std::log(arma::randu()) < alpha_k);
    } catch (const std::runtime_error& e) {
      k_accepted = false;
      ++num_except;
    }
    if (k_accepted) {
      // accept
      k(i+1) = k_star;
      f_store = f_k_star;
      psd_current = psd_k_star;
    } else {
      // reject
      k(i+1) = k(i);
    }
    
    // Step 2: draw Z's
    arma::vec Z_old = Z.col(i);
    for (int l=0; l<L; ++l) {
      arma::vec Z_star = Z_old;
      Z_star(l) += Rcpp::runif(1, -1.0 * eps_Z(l), eps_Z(l))[0];
      if (Z_star(l)<0) Z_star(l)+=1;
      if (Z_star(l)>1) Z_star(l)-=1;
      //
      bernsteinGammaPsd psd_Z_star(psd_current);
      psd_Z_star.update_Z(l, Z_star(l));
      // MH-step
      bool Z_accepted;
      const double f_Z = f_store;
      double f_Z_star;
      try {
        f_Z_star = lpost_bernsteinGammaWhittle(FZ, psd_Z_star, a_prior, k_theta);
        const double alpha_Z = std::min(0.0, f_Z_star - f_Z);
        Z_accepted = (std::log(arma::randu()) < alpha_Z);
      } catch (const std::runtime_error& e) {
        Z_accepted = false;
        ++num_except;
      }
      if (Z_accepted) {
        // accept
        Z(l,i+1) = Z_star(l);
        f_store = f_Z_star;
        psd_current = psd_Z_star;
        Z_old = Z_star;
      } else {
        // reject
        Z(l,i+1) = Z(l,i);
      }
    }
    
    // Step 3: Draw r's
    arma::vec r_old = r.col(i);
    for (int l=0; l<L; ++l) {
      arma::vec r_star = r_old;
      r_star(l) = Rcpp::rlnorm(1, std::log(r_old(l)), std::sqrt(eps_r(l)))[0];
      // MH-step
      bool r_accepted;
      const double f_r = f_store;
      double f_r_star;
      bernsteinGammaPsd psd_r_star(psd_current);
      if (r_star(l) < numerical_thresh) {
        r_accepted = false;
      } else {
        psd_r_star.update_r(l, r_star(l));
        try {
          f_r_star = lpost_bernsteinGammaWhittle(FZ, psd_r_star, a_prior, k_theta);
          const double alpha_r = std::min(0.0, f_r_star + std::log(r_star(l)) - 
                                          f_r - std::log(r_old(l)));
          r_accepted = (std::log(arma::randu()) < alpha_r);
        } catch (const std::runtime_error& e) {
          r_accepted = false;
          ++num_except;
        }
      }
      if (r_accepted) {
        // accept
        r(l,i+1) = r_star(l);
        f_store = f_r_star;
        psd_current = psd_r_star;
        r_old = r_star;
      } else {
        // reject
        r(l,i+1) = r(l,i);
      }
    }
    
    // Step 4: Draw U's
    arma::mat phi_old = U_phi.slice(i);
    for (int l=0; l<L; ++l) {
      arma::mat phi_star = phi_old;
      NumericVector eps_star = Rcpp::runif(phi_dim, -1.0*eps_U(l), 1.0*eps_U(l));
      for (int j=0; j<phi_dim; ++j) {
        eps_star(j) *= phi_def(j);
        phi_star(j,l) += eps_star(j);
        if (phi_star(j,l)<0) phi_star(j,l)+=phi_def(j);
        if (phi_star(j,l)>phi_def(j)) phi_star(j,l)-=phi_def(j);
      }
      bernsteinGammaPsd psd_U_star(psd_current);
      psd_U_star.update_U_phi(l, phi_star.col(l));
      bool U_accepted;
      const double f_U = f_store;
      double f_U_star;
      try {
        f_U_star = lpost_bernsteinGammaWhittle(FZ, psd_U_star, a_prior, k_theta);
        const double alpha_U = std::min(0.0, f_U_star + unit_trace_jacobian_log_determinant(phi_star.col(l)) - 
                                        f_U - unit_trace_jacobian_log_determinant(phi_old.col(l)));
        U_accepted = (std::log(arma::randu()) < alpha_U);
      } catch (const std::runtime_error& e) {
        U_accepted = false;
        ++num_except;
      }
      if (U_accepted) {
        // accept
        for (int j=0; j<phi_dim; ++j) {
          U_phi(j,l,i+1) = phi_star(j,l);
        }
        psd_current = psd_U_star;
        f_store = f_U_star;
        phi_old = phi_star;
      } else {
        // reject
        for (int j=0; j<phi_dim; ++j) {
          U_phi(j,l,i+1) = U_phi(j,l,i);
        }
      }
    }
    
    // Step 5: missing values
    if (NUM_NA > 0) {
      const arma::cx_mat f_0 = psd_current.eval().slice(0);
      const double sigma_NA_prop = 2.0 * SQRT_PI * std::sqrt(arma::norm(f_0, 2));
      const NumericVector innov_prob(sigma_NA_prop * Rcpp::rt(NUM_NA,4.0));
      for (int j=0; j<NUM_NA; ++j) {
        const int jpos = NA_pos[j];
        arma::mat data_star(data);
        data_star[jpos] += innov_prob[j];
        mean_center(data_star);
        const arma::cx_mat FZ_star(mdft_cpp(data_star));
        bool NA_accepted;
        const double f_NA(f_store);
        double f_NA_star;
        try {
          f_NA_star = lpost_bernsteinGammaWhittle(FZ_star, psd_current, a_prior, k_theta);
          const double alpha_NA = std::min(0.0, f_NA_star - f_NA);
          NA_accepted = (std::log(arma::randu()) < alpha_NA);
        } catch (const std::runtime_error& e) {
          NA_accepted = false;
          ++num_except;
        }
        if (NA_accepted) {
          // accept
          data = data_star;
          FZ = FZ_star;
          f_store = f_NA_star;
          missingValuesTrace(i+1,j) = data_star[jpos];
        } else {
          // reject
          missingValuesTrace(i+1,j) = missingValuesTrace(i,j);
        }
      }
    }
    lpostTrace(i+1) = f_store;
  }
  
  List res;
  res["k"]=k; // trace of k (Bernstein polynomial degree)
  res["r"]=r; // trace of r (radial part of AGamma process)
  res["Z"]=Z; // trace or Z (mass allocation points of AGamma process)
  res["U_phi"]=U_phi; // trace of U_phi (spherical coordinate representation of unit trace matrices of AGamma process)
  res["lpostTrace"]=lpostTrace; // log posterior trace (for convergence diagnostics)
  res["num_except"]=num_except; // number of runtime exceptions caught due to numerical issues
  res["missingValuesTrace"]=missingValuesTrace; // trace of missing values (if any)
  return res;
}
