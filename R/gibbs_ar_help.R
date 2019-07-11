#' Time domain AR(p) likelihood for nuisance/noise time series
#' @keywords internal
llike_AR <- function(x_noise, rho, sigma2, full_lik) {
  n <- length(x_noise)
  p <- length(rho)
  a <- pacf_to_ar(rho) 
  eps <- genEpsARMAC(x_noise, a, numeric(0))
  cll <- sum(dnorm(eps, mean=0, sd=sqrt(sigma2), log=T))
  if (full_lik && p>0) {
    zt_p <- head(x_noise, p)
    gamma_p <- ltsa::tacvfARMA(phi=a, maxLag=p-1, sigma2=sigma2)
    Gamma_p <- acvMatrix(gamma_p)
    Gamma_p_inv <- solve(Gamma_p)
    mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
    # marginal log likelihood (for x_1,...,x_p):
    mll <- -log(2*pi)*p/2 -log(det(Gamma_p)) / 2 + mll_unnormalized
  } else {
    mll <- 0
  }
  cll+mll
}

#' Log prior for PACF (~Beta) and sigma2 (~InverseGamma), unnormalized
#' @keywords internal
lprior_AR <- function(rho, rho.alpha, rho.beta, sigma2, sigma2.alpha, sigma2.beta) {
  sum(dbeta((rho+1)/2, rho.alpha, rho.beta, log=T)) -
    (sigma2.alpha + 1) * log(sigma2) - sigma2.beta / sigma2
}

#' Log Posterior = Log Prior + (conditional) Log Likelihood
#' @keywords internal
lpost_AR <- function(x_noise, 
                     rho, rho.alpha, rho.beta,
                     sigma2, sigma2.alpha, sigma2.beta,
                     full_lik) { 
  llike_AR(x_noise, rho, sigma2, full_lik) + 
    lprior_AR(rho, rho.alpha, rho.beta,
              sigma2, sigma2.alpha, sigma2.beta)
}