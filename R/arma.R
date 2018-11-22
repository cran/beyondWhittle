#' Convert partial autocorrelation coefficients to AR coefficients.
#' @details See Section 2 in Kirch et al (2018) or Section III in Barndorff-Nielsen and Schou (1973) for further details
#' @param pacf numeric vector of partial autocorrelations in (-1,1)
#' @return numeric vector of autoregressive model coefficients
#' @references C. Kirch et al
#' Supplemental material of
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' Bayesian Analysis
#' <doi:10.1214/18-BA1126SUPP>
#' @references O. Barndorff-Nielsen and G. Schou
#' On the parametrization of autoregressive models by partial autocorrelations
#' Journal of Multivariate Analysis (3),408-419
#' <doi:10.1016/0047-259X(73)90030-4>
#' @seealso \link[stats]{acf2AR}, \link[stats]{ARMAacf}
#' @export
pacf_to_ar <- function(pacf) {
  p <- length(pacf)
  if (p==0) {
    return(numeric(0))
  }
  if (p==1) {
    return(pacf)
  }
  if (p > 1) {
    return(pacf2AR(pacf)[p,])
  }
}

#' Negative ARMA(p, q) log likelihood
#' @keywords internal
arma_conditional <- function(zt, ar, ma, nll_fun, full_lik, ...) {
  if (!is.null(nll_fun)) {
    stop("nll_fun parameter deprecated")
  }
  if (!(is.null(ma) || length(ma)==0)) {
    stop("MA component not supported yet in likelihood")
  }
  sigma2 <- 1
  n <- length(zt)
  p <- length(ar)
  eps <- genEpsARMAC(zt, ar, numeric(0))
  # conditional log
  cll <- sum(dnorm(eps, mean=0, sd=sqrt(sigma2), log=T))
  if (full_lik) {
    zt_p <- head(zt, p)
    gamma_p <- ltsa::tacvfARMA(phi=ar, maxLag=p-1, sigma2=sigma2)
    Gamma_p <- acvMatrix(gamma_p)
    Gamma_p_inv <- solve(Gamma_p)
    mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
    # marginal log likelihood (for x_1,...,x_p):
    mll <- -log(2*pi)*p/2 -log(det(Gamma_p)) / 2 + mll_unnormalized
  } else {
    mll <- 0
  }
  -(cll+mll) # return negative log likelihood
}