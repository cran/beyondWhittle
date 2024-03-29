# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' C++ function for computing AR coefficients from PACF.
#' See Section III in Barndorff-Nielsen and Schou (1973)
#' @references O. Barndorff-Nielsen and G. Schou (1973)
#' \emph{On the Parametrization of Autoregressive Models by Partial Autocorrelations}
#' Journal of Multivariate Analysis (3),408-419
#' <doi:10.1016/0047-259X(73)90030-4>
#' @keywords internal
pacf2AR <- function(pacf) {
    .Call(`_beyondWhittle_pacf2AR`, pacf)
}

#' Get epsilon process (i.e. model residuals) for ARMA(p,q)
#' @keywords internal
genEpsARMAC <- function(zt, ar, ma) {
    .Call(`_beyondWhittle_genEpsARMAC`, zt, ar, ma)
}

#' ARMA(p,q) spectral density function
#' 
#' Evaluate the ARMA(p,q) spectral density at some frequencies freq in [0,pi),
#' Note that no test for model stationarity is performed.
#' @details See section 4.4 in the referenced book
#' @param freq numeric vector of frequencies to evaluate the psd, 0 <= freq < pi
#' @param ar autoregressive coefficients of ARMA model (use numeric(0) for empty AR part)
#' @param ma moving average coefficients of ARMA model (use numeric(0) for empty MA part)
#' @param sigma2 the model innovation variance
#' @return numeric vector of the (real-valued) spectral density values
#' @references P. J. Brockwell and R. Davis (1996)
#' \emph{Time Series: Theory and Methods (Second Edition)}
#' @export
psd_arma <- function(freq, ar, ma, sigma2 = 1.0) {
    .Call(`_beyondWhittle_psd_arma`, freq, ar, ma, sigma2)
}

#' Get  p from v in Stick Breaking DP representation
#' @keywords internal
pFromV <- function(v) {
    .Call(`_beyondWhittle_pFromV`, v)
}

#' Get v from p (DP inverse stick breaking)
#' Note: p is assumed to have length L, i.e. it does NOT contain p_0
#' @keywords internal
vFromP <- function(p, eps = 1e-8) {
    .Call(`_beyondWhittle_vFromP`, p, eps)
}

#' Get mixture weights of Bernstein-Dirchlet-Mixtures 
#' @keywords internal
mixtureWeight <- function(p, w, k) {
    .Call(`_beyondWhittle_mixtureWeight`, p, w, k)
}

#' Construct a density mixture from mixture weights and density functions.
#' @keywords internal
densityMixture <- function(weights, densities) {
    .Call(`_beyondWhittle_densityMixture`, weights, densities)
}

#' Redundantly roll out a PSD from length N=floor(n/2) to length n
#' @keywords internal
unrollPsd <- function(qPsd, n) {
    .Call(`_beyondWhittle_unrollPsd`, qPsd, n)
}

#' Gibbs sampler in Cpp
#' @keywords internal
gibbs_multivariate_nuisance_cpp <- function(data, NA_pos, FZ, eps_r, eps_Z, eps_U, k_0, r_0, Z_0, U_phi_0, phi_def, eta, omega, Sigma, Ntotal, print_interval, numerical_thresh, verbose, L, k_theta, dbList) {
    .Call(`_beyondWhittle_gibbs_multivariate_nuisance_cpp`, data, NA_pos, FZ, eps_r, eps_Z, eps_U, k_0, r_0, Z_0, U_phi_0, phi_def, eta, omega, Sigma, Ntotal, print_interval, numerical_thresh, verbose, L, k_theta, dbList)
}

#' Store imaginary parts above and real parts below the diagonal
#' @keywords internal
realValuedPsd <- function(f_) {
    .Call(`_beyondWhittle_realValuedPsd`, f_)
}

#' Inverse function to realValuedPsd
#' @keywords internal
complexValuedPsd <- function(f_) {
    .Call(`_beyondWhittle_complexValuedPsd`, f_)
}

#' Construct psd mixture
#' @keywords internal
get_f_matrix <- function(U_phi, r, Z, k, dbList) {
    .Call(`_beyondWhittle_get_f_matrix`, U_phi, r, Z, k, dbList)
}

#' I/O: Only used within Rcpp
#' 
#' This workaround for parsing cubes was neccessary at development time
#' because there was a (presumable) bug in RcppArmadillo that sometimes
#' caused the parsing of arma::cx_cube objects to fail, such that the function
#' received an un-initialized object instead of the parsed one.
#' 
#' The workaround parses an Rcpp vector instead, and manually
#' copies the data in an arma::cx_cube object.
#' Besides being redundant, it also makes the code less readable and it is
#' hoped that this workaround can be removed in future revisions.
#' 
#' @keywords internal
cx_cube_from_ComplexVector <- function(x) {
    .Call(`_beyondWhittle_cx_cube_from_ComplexVector`, x)
}

#' I/O: Only used within Rcpp
#' Note: Same workaround as \code{cx_cube_from_ComplexVector}
#' @keywords internal
cube_from_NumericVector <- function(x) {
    .Call(`_beyondWhittle_cube_from_NumericVector`, x)
}

#' Does a matrix have an eigenvalue smaller than 0?
#' @keywords internal
hasEigenValueSmallerZero <- function(A, TOL = 0.0) {
    .Call(`_beyondWhittle_hasEigenValueSmallerZero`, A, TOL)
}

#' Build an n times n Toeplitz matrix from the 
#' autocovariance values gamma(0),...,gamma(n-1)
#' @keywords internal
acvMatrix <- function(acv) {
    .Call(`_beyondWhittle_acvMatrix`, acv)
}

#' Build an nd times nd Block Toeplitz matrix from the
#' (d times d) autocovariances gamma(0),...,gamma(n-1)
#' @keywords internal
acvBlockMatrix <- function(acv) {
    .Call(`_beyondWhittle_acvBlockMatrix`, acv)
}

#' Computing acceptance rate based on trace
#' Note: Only use for traces from continous distributions!
#' @keywords internal
acceptanceRate <- function(trace) {
    .Call(`_beyondWhittle_acceptanceRate`, trace)
}

qpsd_cal_cpp_expedited <- function(basis1, basis2, p, pexpend, selector1cpp, selector2cpp) {
    .Call(`_beyondWhittle_qpsd_cal_cpp_expedited`, basis1, basis2, p, pexpend, selector1cpp, selector2cpp)
}

#' Get U from phi, vectorized, cpp internal only
#' @keywords internal
get_U_cpp <- function(u_phi) {
    .Call(`_beyondWhittle_get_U_cpp`, u_phi)
}

#' Get x from phi, see (62) in Mittelbach et al.
#' @keywords internal
unit_trace_x_from_phi <- function(phi) {
    .Call(`_beyondWhittle_unit_trace_x_from_phi`, phi)
}

#' Get L (lower triangular Cholesky) from x
#' Called U^* in Mittelbach et al, see (60) there
#' @keywords internal
unit_trace_L_from_x <- function(x) {
    .Call(`_beyondWhittle_unit_trace_L_from_x`, x)
}

#' Get p vector, see (67) in Mittelbach et al.
#' @keywords internal
unit_trace_p <- function(d) {
    .Call(`_beyondWhittle_unit_trace_p`, d)
}

#' Get q vector, see (68) in Mittelbach et al.
#' @keywords internal
unit_trace_q <- function(d) {
    .Call(`_beyondWhittle_unit_trace_q`, d)
}

#' Get VARMA PSD from transfer polynomials 
#' Helping function for \code{psd_varma}
#' @keywords internal
varma_transfer2psd <- function(transfer_ar_, transfer_ma_, sigma) {
    .Call(`_beyondWhittle_varma_transfer2psd`, transfer_ar_, transfer_ma_, sigma)
}

#' VARMA transfer polynomials
#' @keywords internal
transfer_polynomial <- function(lambda, coef) {
    .Call(`_beyondWhittle_transfer_polynomial`, lambda, coef)
}

#' epsilon process (residuals) of VAR model
#' @keywords internal
epsilon_var <- function(zt, ar) {
    .Call(`_beyondWhittle_epsilon_var`, zt, ar)
}

#' sum of multivariate normal log densities
#' with mean 0 and covariance Sigma, unnormalized
#' @keywords internal
sldmvnorm <- function(z_t, Sigma) {
    .Call(`_beyondWhittle_sldmvnorm`, z_t, Sigma)
}

