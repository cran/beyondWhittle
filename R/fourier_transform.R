#' Fourier frequencies
#' 
#' Fourier frequencies on [0,pi], as defined by 2*pi*j/n for j=0,...,floor(n/2).
#' @param n integer
#' @return numeric vector of length floor(n/2)+1
#' @export
fourier_freq <- function(n) {
  pi*omegaFreq(n)
}

#' Fourier frequencies rescaled on the unit interval
#' @keywords internal
omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

#' Fast Fourier Transform
#' @details If \code{real}: computes F_n X_n with the real-valued Fourier 
#' transformation matrix F_n (see Section 2.1 in Kirch et al (2018)).
#' If \code{!real}: computes the complex-valued Fourier coefficients 
#' (see (4.5) in Meier (2018)).
#' @keywords internal
fast_ft <- compiler::cmpfun(function(x, real=T) {
  n <- length(x)
  sqrt2 <- sqrt(2)
  sqrtn <- sqrt(n)
  # Cyclically shift so last observation becomes first
  # Important since fft() uses 0:(n-1) but we use 1:n
  x <- c(x[n], x[-n])  
  fourier <- fft(x)
  if (real) {
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
  } else {
    N <- ifelse(n %% 2, (n+1)/2, n/2+1)
    FZ <- fourier[1:N]
  }
  return(FZ / sqrtn)
})

#' Fast Inverse Fourier Transform
#' @details inverse function of \code{fast_ft}
#' @keywords internal
fast_ift <- compiler::cmpfun(function(x, real=T, TOL=1e-15) {
  if (real) {
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
  } else {
    N <- length(x)
    n_is_even <- (abs(Im(x[N])) < TOL)
    n <- ifelse(n_is_even, 2*(N-1), 2*N-1)
    CFZ <- c(x, rev(Conj(x[-c(1,N)]))) * sqrt(n)
  }
  # Inverse FFT (normalised)
  Z <- fft(CFZ, inverse = TRUE) / n
  # Cyclically shift
  Z <- c(Z[-1], Z[1])
  if (real) {
    Z <- Re(Z)
  }
  return(Z)
})

#' Multivariate discrete (fast) Fourier Transform
#' @keywords internal
mdft <- function(Z, real=F) {
  FZ <- apply(Z, 2, fast_ft, real=real)
  return(FZ)
}

#' Multivariate inverse discrete (fast) Fourier Transform
#' @keywords internal
midft <- function(FZ, real=F) {
  Z <- apply(FZ, 2, fast_ift, real=real)
  return(Z)
}

#' Compute Periodgram matrix from (complex-valued) Fourier coefficients
#' @details see (4.7) in Meier (2018)
#' @keywords internal
mpdgrm <- function(FZ) {
  N <- nrow(FZ)
  d <- ncol(FZ)
  res <- array(data=NA, dim=c(d,d,N))
  for (j in 1:N) {
    res[,,j] <- FZ[j,] %*% Adj(FZ[j,])
  }
  res <- res / 2 / pi
  return(res)
}