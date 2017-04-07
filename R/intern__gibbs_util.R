# library(Rcpp)
# library(forecast)
# library(MASS)
# library(compiler)

#' Fourier frequencies, rescaled on the unit interval
#' @param n integer length in consideration
#' @return numeric vector of the rescaled Fourier frequencies
#' @export
omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

#' Convert partial autocorrelation coefficients to AR coefficients.
#' @param pacf numeric vector of partial autocorrelations in (-1,1)
#' @return numeric vector of autoregressive model coefficients
#' @export
pacfToAR <- function(pacf) {
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

#' Help function to compute the mean.
#' @keywords internal
fast_mean <- function(x) {
  sum(x) / length(x)
}

#' Compute F_n X_n with the real-valued Fourier matrix F_n
#' @keywords internal
fast_ft <- compiler::cmpfun(function(x) {
  # Function computes FZ (i.e. fast Fourier transformed data)
  # Outputs coefficients in correct order and rescaled
  n <- length(x)
  sqrt2 <- sqrt(2)
  sqrtn <- sqrt(n)
  # Cyclically shift so last observation becomes first
  x <- c(x[n], x[-n])  # Important since fft() uses 0:(n-1) but we use 1:n
  # FFT
  fourier <- fft(x)
  # Extract non-redundant real and imaginary coefficients in correct order and rescale
  FZ <- rep(NA, n)
  FZ[1] <- Re(fourier[1]) # First coefficient
  if (n %% 2) {
    N <- (n-1)/2
    FZ[2*(1:N)] <- sqrt2 * Re(fourier[2:(N+1)]) # Real coefficients
    FZ[2*(1:N)+1] <- sqrt2 * Im(fourier[2:(N+1)]) # Imaginary coefficients
  } else {
    FZ[n] <- Re(fourier[n / 2 + 1]) # Last coefficient
    FZ[2 * 1:(n / 2 - 1)] <- sqrt2 * Re(fourier[2:(n / 2)]) # Real coefficients
    FZ[2 * 1:(n / 2 - 1) + 1] <- sqrt2 * Im(fourier[2:(n / 2)]) # Imaginary coefficients
  }
  return(FZ / sqrtn)
})

#' Compute F_n^t X_n with the real-valued Fourier matrix F_n
#' @keywords internal
fast_ift <- compiler::cmpfun(function(x) {
  # Function computes inverse Fourier transform
  # Can be used for finding FCFZ
  n <- length(x)
  sqrtn <- sqrt(n)
  sqrtn2 <- sqrt(n / 2)
  # Construct complex vector
  CFZ <- rep(NA, n)
  CFZ[1] <- x[1] * sqrtn
  if (n %% 2) {
    N <- (n-1)/2
    CFZ[2:(N+1)] <- (x[2 * (1:N)] + 1i * x[2 * (1:N)+1] ) * sqrtn2
    CFZ[(N+2):n] <- rev(Conj(CFZ[2:(N+1)])) # Include complex complex conjugates
  } else {
    CFZ[n / 2 + 1] <- x[n] * sqrtn
    CFZ[2:(n / 2)] <- (x[2 * (1:(n / 2 - 1))] + x[2 * (1:(n / 2 - 1)) + 1] * 1i) * sqrtn2
    CFZ[(n / 2 + 2):n] <- rev(Conj(CFZ[2:(n / 2)])) # Include complex complex conjugates
  }
  # Inverse FFT (normalised)
  FCFZ <- fft(CFZ, inverse = TRUE) / n
  # Cyclically shift
  FCFZ <- c(FCFZ[-1], FCFZ[1])
  return(Re(FCFZ))
})

#' Density of t-distribution in terms of excess kurtosis
#' @keywords internal
dtex.kurt <- function(x, ex.kurt) {
  nu <- 6 / ex.kurt + 4
  dt(x, nu)
}

#' Help function: Uniform maximu
#' @keywords internal
uniformmax <- function(sample) {
  max(abs(sample - median(sample)) / mad(sample), na.rm=T)
}

#' Help function: Fuller Logarithm
#' @keywords internal
logfuller<-function(x, xi = 0.001){
  log(x + xi) - xi / (x + xi)
}

#' Help function: Rosenthal Logarithm
#' @keywords internal
logrosenthal <- function(x, inverse=F) {
  if (!inverse) return(sign(x) * log(1 + abs(x)))
  else return(sign(x) * (exp(abs(x)) - 1))
}

#' Help function: Standard error for kurtosis
#' @keywords internal
se_kurt <- function(n) {
  return(2 * sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3))) *
           sqrt((n ^ 2 - 1) / ((n - 3) * (n + 5))))
}

#' Negative log likelihood of iid standard normal observations [unit variance]
#' @keywords internal
nll_norm <- function(epsilon_t, ...) {
  m <- length(epsilon_t)
  cll <-  1 / 2 * (sum(epsilon_t ^ 2) + m * log(2*pi))
  return(cll)
}

#' unnormalized negative log likelihood of iid standard normal observations [unit variance]
#' @keywords internal
nll_norm_unnormalized <- function(epsilon_t, ...) {
  cll <-  1 / 2 * sum(epsilon_t ^ 2)
  return(cll)
}

#' negative log likelihood of iid t observations with given excess kurtosis [unit variance]
#' @keywords internal
nll_t <- function(epsilon_t, ex.kurt, ...) {
  nu <- (6 / ex.kurt) + 4
  sigma2 <- nu / (nu - 2)
  cll <- -sum(dt(epsilon_t * sqrt(sigma2),df=nu,log=T)) - log(sigma2)/2
  return(cll)
}

#' Help function for Generalized Gaussian
#' @keywords internal
generalizedGaussian.alpha <- function(beta) { # for unit variance
  lalpha <- 1 / 2 * (lgamma(1/beta) - lgamma(3/beta))
  alpha <- exp(lalpha)
  return(alpha)
}
#' Help function for Generalized Gaussian
#' @keywords internal
generalizedGaussian.kurtosis <- function(beta) {
  log.kurt <- lgamma(5 / beta) + lgamma(1 / beta) - 2*lgamma(3 / beta)
  ex.kurt <- exp(log.kurt) - 3
  return(ex.kurt)
}
#' Help function for Generalized Gaussian
#' @keywords internal
l_generalizedGaussian <- function(x, beta) {
  alpha <- generalizedGaussian.alpha(beta)
  llik <- log(beta) - log(2*alpha) - lgamma(1/beta) - (abs(x) / alpha)^beta
  return(llik)
}
#' Help function for Generalized Gaussian
#' @keywords internal
nll_generalizedGaussian <- function(epsilon_t, beta) {
  cll <- -sum(l_generalizedGaussian(epsilon_t, beta))
  return(cll)
}

#' Help function for visualizing
#' @keywords internal
plotPsdEstimate <- function(mcmc, lambda, psd.true) {
  plot(lambda, psd.true, type = "l", ylim = c(0, 2*max(psd.true)))
  lines(lambda, mcmc$fpsd.s, col = 2, lty = 2)
  lines(lambda, mcmc$fpsd.s05, col = 2, lty = 2)  # PW CI
  lines(lambda, mcmc$fpsd.s95, col = 2, lty = 2)  # PW CI
  lines(lambda, mcmc$log.confupper, col = 4, lty = 2)  # Uniform CI
  lines(lambda, mcmc$log.conflower, col = 4, lty = 2)  # Uniform CI
  title(mcmc$algorithm)
}
#' Help function for visualizing
#' @keywords internal
plotMCMC <- function(mcmc, lambda, psd.true, ex.kurt.true=ex.kurt.true) {
  plotPsdEstimate(mcmc, lambda, psd.true)
  plot.ts(mcmc$tau)
  plot.ts(mcmc$k)
  if (!(any(is.na(mcmc$ex.kurt)))) {
    plot.ts(mcmc$ex.kurt)
    abline(h = mean(mcmc$ex.kurt), col = 2, lty = 2)
    abline(h = median(mcmc$ex.kurt), col = 3, lty = 2)
    abline(h = ex.kurt.true, col = 4, lty = 2)
  } else {
    plot.ts(0)
  }
  if (!(any(is.na(mcmc$rho)))) { # TODO Adjust if rho is matrix
    plot.ts(as.numeric(mcmc$rho))
  } else {
    plot.ts(0)
  }
}
#' Help function for I/O
#' @keywords internal
filenameMCMC <- function(n, repN, ar.ex, ma.ex, ex.kurt.true, thin) {
  ar.string <- ifelse(length(ar.ex) > 0, paste0("_ar(", paste0(ar.ex, collapse="-"), ")"), "")
  ma.string <- ifelse(length(ma.ex) > 0, paste0("_ma(", paste0(ma.ex, collapse="-"), ")"), "")
  filename <- paste0("n", n, "_repN", repN, ar.string, ma.string, "_kurt", ex.kurt.true, "_thin", thin)
}
#' Help function for I/O
#' @keywords internal
reduceMemoryStorageMCMC <- function(mcmc) { # Discard memory intensive traces
    ret <- mcmc
    ret$fpsd.sample <- NULL

    ret$pdgrm <- NULL
    ret$V <- NULL
    ret$W <- NULL
    ret$df <- NULL

    ret$cy <- NULL
    ret$mixtureTrace <- NULL
    ret$phi_mu <- NULL
    ret$phi_sigma2 <- NULL

    return(ret)
  }
