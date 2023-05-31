#' Construct coarsened Bernstein polynomial basis of degree l on omega
#' @param omega numeric vector in [0,1] of evaluation points
#' @param l positive integer for the degree
#' @details See Appendix E.1 in Ghosal/Van der Vaart, Fundamentals, (2017)
#' @keywords internal
coarsened_bernstein <- function(omega, l) {
  res <- matrix(NA, nrow=l, ncol=length(omega))
  for (i in 1:l) {
    res[i,] <- coarsened_bernstein_i(omega, l, i)
  }
  res
}

#' Helping function for \code{coarsened_bernstein}
#' @keywords internal
coarsened_bernstein_i <- function(omega, l, i) {
  k <- l^2
  b_tmp <- 0 * omega
  for (j in ((i-1)*l+1):(i*l)) {
    b_tmp <- b_tmp + dbeta(omega, j, k+1-j)
  }
  b_tmp <- b_tmp / l
  b_tmp
}

#' Construct Bernstein polynomial basis of degree k on omega
#' @param omega numeric vector in [0,1] of evaluation points
#' @param k positive integer for the degree
#' @param coarsened bool flag indicating whether coarsened or standard Bernstein polynomials are used
#' @keywords internal
betaBasis_k <- function(omega, k, coarsened) {
  if (coarsened) {
    basis <- coarsened_bernstein(omega, k)
  } else {
    N <- length(omega)
    basis <- matrix(dbeta(omega,                                                 
                          rep(1:k, each = N),
                          rep(k:1, each = N)),
                    ncol = N,
                    byrow = TRUE)
  }
  basis
}

#' Construct Bernstein polynomial basises of degree up to kmax on omega
#' @param n positive integer determining the number of the (equidistant) evaluation points in [0,1]
#' @param kmax positive integer for the largest degree
#' @param bernstein_l,bernstein_r left and right truncation
#' @param coarsened bool flag indicating whether coarsened or standard Bernstein polynomials are used
#' @param verbose debugging parameter
#' @return A list of length kmax, where the k-th list element is a matrix containing the polynomial basis of degree k
#' @keywords internal
dbList <- function(n, kmax, bernstein_l=0, bernstein_r=1, coarsened=F, verbose=F) {
  db.list <- vector("list", kmax)
  omega <- omegaFreq(n); NN <- length(omega)
  omega_for_dblist <- seq(bernstein_l, bernstein_r, length.out=NN)
  if (verbose) {
    if (coarsened) {
      cat("Using coarsened Bernstein polynomials on (", bernstein_l , ",", bernstein_r, ")\n")
    } else {
      cat("Using standard Bernstein polynomials on (", bernstein_l , ",", bernstein_r, ")\n")
    }
  }
  if (verbose) cat("Precomputing polynomial basis functions ")
  for (kk in 1:kmax) {
    db.list[[kk]] <- betaBasis_k(omega_for_dblist, kk, coarsened)
    if (verbose && !(kk%%10)) cat(".")
  }
  if (verbose) cat(" done!\n")
  return(db.list)
}