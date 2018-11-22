##
## This file contains methods for Hpd matrices with unit trace.
## It is based on the hyperspherical coordinates approach from 
## Mittelbach et al (2012), see also Section 3.4.1 in Meier (2018).
## 
## M. Mittelbach, B. Matthiesen and E. A. Jorswieck (2012).
## "Sampling uniformly from the set of positive definite matrices with 
## trace constraint". IEEE Transactions on Signal Processing, 60(5):2167-2179.
##
## A. Meier (2018). "A Matrix Gamma Process and Applications to Bayesian 
## Analysis of Multivariate Time Series". PhD thesis, OvGU Magdeburg.
##

#' Get U (Hpd with unit trace) matrix from 
#' phi (hyperspherical coordinates) vector.
#' @keywords internal
unit_trace_U_from_phi <- function(phi) {
  x <- unit_trace_x_from_phi(phi)
  L <- unit_trace_L_from_x(x)
  return(L %*% Adj(L))
}

#' Range intervals I_l, see (63) in Mittelbach et al.
#' @keywords internal
unit_trace_I_l <- function(l) {
  stopifnot(length(l)==1)
  if (is_quadratic(l)) {
    I_l <- c(0,pi/2)
  } else {
    I_l <- c(0,pi)
  }
  I_l
}

#' Draw uniformly from Hpd matrices with unit trace
#' @keywords internal
unit_trace_runif <- function(n, d, verbose=F) {
  N <- d*d-1
  phi_res <- matrix(NA, nrow=N, ncol=n)
  U_res <- array(NA, dim=c(d,d,n))
  for (j in 1:n) {
    if (verbose) print(j)
    tmp <- unit_trace_runif_single(d)
    phi_res[,j] <- tmp$phi
    U_res[,,j] <- tmp$U
  }
  list(phi=phi_res, U=U_res)
}

#' Obtain one uniform draw from d times d Hpd matrices with unit trace
#' See Algorithm 2 in Mittelbach et al. (adjusted to complex case)
#' @keywords internal
unit_trace_runif_single <- function(d) {
  N <- d*d-1
  phi_res <- rep(NA, N)
  p <- unit_trace_p(d)
  q <- unit_trace_q(d)
  log_c <- unit_trace_log_c(p,q)
  log_d <- unit_trace_log_d(p,q)
  mu <- unit_trace_mu(p,q)
  sigma2 <- unit_trace_sigma2(p,q)
  log_nu <- unit_trace_nu(sigma2, log_c, log_d)
  for (l in 1:N) {
    accepted <- F
    while (!accepted) {
      sigma_l <- sqrt(sigma2[l])
      phi_star <- rnorm(1, mu[l], sigma_l)
      alpha <- unit_trace_log_f_l(phi_star,
                                p,
                                q,
                                log_c,
                                l) -
        dnorm(phi_star, mu[l], sigma_l, log=T) -
        log_nu[l]
      accepted <- (log(runif(1,0,1)) < alpha) 
    }
    phi_res[l] <- phi_star
  }
  U <- unit_trace_U_from_phi(phi_res)
  list(phi=phi_res, U=U)
}

#' Get mu vector, see (36) in Mittelbach et al.
#' Helping function for \code{unit_trace_runif}
#' @keywords internal
unit_trace_mu <- function(p, q) { 
  N <- length(p); stopifnot(length(q)==N); stopifnot(N>1)
  res <- rep(pi/2, N)
  i <- 1
  while (i*i < N) {
    l <- i*i
    res[l] <- atan(sqrt(q[l]/p[l]))
    i <- i+1
  }
  res
}

#' Get sigma2 vector, see (70) in Mittelbach et al.
#' Helping function for \code{unit_trace_runif}
#' @keywords internal
unit_trace_sigma2 <- function(p, q) {
  N <- length(p); stopifnot(length(q)==N); stopifnot(N>1)
  1 / (sqrt(p) + sqrt(q))^2
}

#' Get log(nu) vector, see (38) in Mittelbach et al.
#' Helping function for \code{unit_trace_runif}
#' @keywords internal
unit_trace_nu <- function(sigma2_vec, log_c_vec, log_d_vec) { 
  log(2*pi*sigma2_vec) / 2 + log_c_vec - log_d_vec
}

#' Get log(c) vector, see (70) in Mittelbach et al.
#' Helping function for \code{unit_trace_runif}
#' @keywords internal
unit_trace_log_c <- function(p, q) {
  N <- length(p); stopifnot(length(q)==N); stopifnot(N>1)
  res <- rep(NA, N)
  i <- 1
  for (l in 1:N) {
    if (l == i*i) {
      res[l] <- log(2) + lgamma( (p[l]+1)/2 + (q[l]+1)/2 ) - 
        lgamma( (p[l]+1)/2 ) - lgamma( (q[l]+1)/2 )
      i <- i+1
    } else {
      res[l] <- -log(pi)/2 + lgamma( (q[l]+1)/2 + 1/2 ) - 
        lgamma( (q[l]+1)/2 )
    }
  }
  res
}

#' Get log(d) vector, see (39) in Mittelbach et al, adjusted to complex case
#' Helping function for \code{unit_trace_runif}
#' @keywords internal
unit_trace_log_d <- function(p, q) { 
  N <- length(p); stopifnot(length(q)==N && N > 1)
  res <- rep(0,N)
  i <- 1
  while (i*i < N) {
    l <- i*i
    res[l] <- p[l]/2 * log(1+q[l]/p[l]) + q[l]/2 * log(1+p[l]/q[l])
    i <- i + 1
  }
  res
}

#' Is l quadratic?
#' @keywords internal
is_quadratic <- function(l, thresh=1e-15) {
  # a bit hacky
  sl <- sqrt(l)
  sl-as.integer(sl) < thresh
}

#' Get log(f_l), see (66) in Mittelbach et al.
#' Helping function for \code{unit_trace_runif}
#' @keywords internal
unit_trace_log_f_l <- function(phi, p, q, log_c, l) { 
  N <- length(p); stopifnot(length(q) == N && length(log_c) == N && N > 1)
  stopifnot(l >= 1 && l <= N) 
  I_l <- unit_trace_I_l(l)
  if (phi <= I_l[1] || phi >= I_l[2]) {
    lf_l <- -Inf
  }
  else {
    lf_l <- log_c[l] 
    if (p[l] != 0) {
      lf_l <- lf_l + p[l] * log(cos(phi))
    }
    if (q[l] != 0) {
      lf_l <- lf_l + q[l] * log(sin(phi))
    }
  }
  lf_l
}