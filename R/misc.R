#' mean center a numerical vector
#' @keywords internal
center <- function(x, ...) {
  return(x - mean(x, ...))
}

#' Help function to print MCMC state
#' @keywords internal
print_mcmc_state <- function(i, Ntotal, tim, tim0) {
  elapsed <- base::round((tim-tim0)[[3]], digits=2)
  txt <- paste0("MCMC iteration ", i, "/", Ntotal, " (CPU time: ", elapsed, "s)")
  print(txt)
}

#' Help function to print debugging messages
#' @keywords internal
print_warn <- function(msg) {
  print(msg)
  warning(msg)
}

#' Help function to compute the mean.
#' @keywords internal
fast_mean <- function(x) {
  sum(x) / length(x)
}

#' Uniform maximum, as needed for uniform credible intervals
#' @details see Section 4.1 in Kirch et al (2018)
#' @keywords internal
uniformmax <- function(sample) {
  tmp <- abs(sample - median(sample)) / mad(sample)
  if (all(is.na(tmp))) {
    # case of no deviation
    return(-Inf)
  } else {
    max(tmp, na.rm=T)
  }
}

#' Uniform credible intervals in matrix-valued case
#' @details see (6.5) in Meier (2018)
#' @keywords internal
uci_matrix <- function(fpsd.sample, alpha, uniform_among_components=F) {
  d <- dim(fpsd.sample)[1]
  N <- dim(fpsd.sample)[3]
  fpsd.uci05 <- fpsd.uci95 <- array(NA, dim=c(d, d, N))
  if (uniform_among_components) {
    # Use the same C_\alpha^* value for all components
    for (i in 1:d) {
      fpsd.sample[i,i,,] <- log(fpsd.sample[i,i,,])
    }
    fpsd.s <- apply(fpsd.sample, c(1,2,3), median)
    fpsd.mad <- apply(fpsd.sample, c(1,2,3), mad)
    fpsd.help <- uniformmax_multi(fpsd.sample)
    Cvalue <- quantile(fpsd.help, 1-alpha)
    fpsd.uci05 <- fpsd.s - Cvalue * fpsd.mad
    fpsd.uci95 <- fpsd.s + Cvalue * fpsd.mad
    for (i in 1:d) {
      fpsd.uci05[i,i,] <- exp(fpsd.uci05[i,i,])
      fpsd.uci95[i,i,] <- exp(fpsd.uci95[i,i,])
    }
  } else {
    # Use individual C_\alpha^* among each component
    for (i in 1:d) {
      for (j in 1:d) {
        uci_tmp <- uci_help(fpsd.sample[i,j,,], alpha, log=(i==j))
        fpsd.uci05[i,j,] <- uci_tmp$conflower
        fpsd.uci95[i,j,] <- uci_tmp$confupper
      }
    }
  }
  return(list(fpsd.uci05=fpsd.uci05, 
              fpsd.uci95=fpsd.uci95))
}

#' Helping function for \code{uci_matrix}
#' @keywords internal
uci_help <- function(fpsd.sample, alpha, log=F) {
  if (log) {
    fpsd.sample <- log(fpsd.sample) #logfuller(fpsd.sample)
  }
  fpsd.s <- apply(fpsd.sample, 1, median)
  fpsd.mad <- apply(fpsd.sample, 1, mad)
  fpsd.help <- apply(fpsd.sample, 1, uniformmax)
  Cvalue <- quantile(fpsd.help, 1-alpha)
  conflower <- fpsd.s - Cvalue * fpsd.mad
  confupper <- fpsd.s + Cvalue * fpsd.mad
  if (log) {
    conflower <- exp(conflower)
    confupper <- exp(confupper)
  }
  return(list(conflower=conflower,
              confupper=confupper))
}

#' Helping function for \code{uci_matrix}
#' @keywords internal
uniformmax_multi <- function(mSample) {
  d <- dim(mSample)[1]
  N <- dim(mSample)[3]
  N_sample <- dim(mSample)[4]
  C_help <- array(NA, dim=c(d,d,N,N_sample))
  for (j in 1:N) {
    for (r in 1:d) {
      for (s in 1:d) {
        C_help[r,s,j,] <- uniformmax_help(mSample[r,s,j,])
      }
    }
  }
  apply(C_help, 4, max, na.rm=T)
}

#' Helping function for \code{uci_matrix}
#' @keywords internal
uniformmax_help <- function(sample) {
  abs(sample - median(sample)) / mad(sample)
}

#' Fuller Logarithm
#' @details see Fuller (1996), page 496
#' @references W. Fuller (1996)
#' \emph{Introduction to Statistical Time Series}
#' Wiley Series in Probability and Statistics
#' @keywords internal
logfuller<-function(x, xi = 1e-8){
  log(x + xi) - xi / (x + xi)
}

#' Negative log likelihood of iid standard normal observations [unit variance]
#' Note: deprecated
#' @keywords internal
nll_norm <- function(epsilon_t, ...) {
  m <- length(epsilon_t)
  cll <-  1 / 2 * (sum(epsilon_t ^ 2) + m * log(2*pi))
  return(cll)
}

#' adjoint of complex matrix
#' @keywords internal
Adj <- function(m) {
  return(t(Conj(m)))
}

#' Check if a matrix is Hermitian positive definite
#' @keywords internal
is_hpd <- function(A, tol=1e-15) {
  (A==Adj(A)) && (!hasEigenValueSmallerZero(A, tol))
}

#' Check if a matrix is symmetric positive definite
#' @keywords internal
is_spd <- function(A, tol=1e-5) {
  (A==t(A)) && (!hasEigenValueSmallerZero(A, tol))
}

#' Get string representation for missing values position from vector index
#' @keywords internal
missingValues_str_help <- function(ind, n) {
  res <- c()
  for (j in ind) {
    t <- (j-1) %% n + 1 # time point
    d <- j %/% n + 1 # dimension
    res <- c(res, paste0(t, ",", d))
  }
  res
}
