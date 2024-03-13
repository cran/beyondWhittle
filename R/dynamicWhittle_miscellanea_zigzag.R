
#' Calculate the moving Fourier transform ordinates
#'
#' @param x A numeric vector containing time series.
#' @param m A positive integer indicating the size of window.
#' @param thinning_factor Selected from "1, 2, 3".
#'
#' @return A list containing the moving Fourier transform and corresponding time grid.
#' @export
#'
#' @references Y. Tang et al. (2023)
#' \emph{Bayesian nonparametric spectral analysis of locally stationary processes}
#' ArXiv preprint
#' <arXiv:2303.11561>
#' 
#' @examples
#' set.seed(1); x <- rnorm(1500)
#' local_moving_FT_zigzag(x, 50, 1)
local_moving_FT_zigzag <- function(x, m, thinning_factor){

  if(!(thinning_factor %in% 1:3)){return("thinning factor should be one of 1, 2 and 3")}

  n <- length(x) - 2*m

  if(thinning_factor == 1){grid_points <- (1:n)} else {

  index <- (floor(((1:n) - 1)/m) %% thinning_factor == 0) | (floor(((1:n) - 1)/m) %in% c(floor(n/m) - 1, floor(n/m)))

  grid_points <- (1:n)[index]


  }

  y <- integer(length(grid_points))

  for (t in grid_points) {
    xnew <- x[(t-m):(t+m) + m]
    lambda <- 2*pi*(1+(t-1)%%m)/(2*m+1)
    y[grid_points == t] <- (1/sqrt(2*m+1)) * sum(xnew * exp(-1i * (0:(2*m)) * lambda))
  }

  return(list(time_grid = grid_points, MF = y))
}


#' Calculation of log prior
#' @importFrom stats dbeta
#' @importFrom stats plogis
#' @importFrom stats dlogis
#' @keywords internal
lprior_dw <- function(tilde.v, tilde.w1, tilde.w2, k1, k2, tau,
                      M, g0.alpha, g0.beta, k1.theta, k2.theta, tau.alpha, tau.beta){

  lp_tilde.v <- sum(dbeta(plogis(tilde.v), 1, M, log = T) + dlogis(tilde.v, log = T) ) # log prior for V's - beta(1, M)
  lp_tilde.w1 <- sum(dbeta(plogis(tilde.w1), g0.alpha, g0.beta, log = T) + dlogis(tilde.w1, log = T) )# log prior for W1's - beta(a, b)
  lp_tilde.w2 <- sum(dbeta(plogis(tilde.w2), g0.alpha, g0.beta, log = T) + dlogis(tilde.w2, log = T) )# log prior for W1's - beta(a, b)
  lp_k1 <- -k1.theta * k1 * log(k1) # log prior for k1
  lp_k2 <- -k2.theta * k2 * log(k2) # log prior for k2
  lp_tau <- -(tau.alpha + 1) * log(tau) - tau.beta / tau # log prior for tau (Inverse Gamma)


  lp <- lp_tilde.v + lp_tilde.w1 + lp_tilde.w2 + lp_k1 + lp_k2 + lp_tau

  return(lp)

}

#' Evaluation of normalized time-varying spectral density function (based on posterior samples)
#' @keywords internal
qpsd_dw <- function(v, w1, w2, k1, k2, beta_basis_1_k, beta_basis_2_k) {

  p <- pFromV(v)

  selector1 <- findInterval(w1, (1:k1)/k1, left.open = T) + 1 # amount to ceiling(k1 * w1) but safer
  selector2 <- findInterval(w2, (1:k2)/k2, left.open = T) + 1 # amount to ceiling(k2 * w2) but safer

  psd <- crossprod(p * beta_basis_2_k[selector2, ],  beta_basis_1_k[selector1, ])
  return(psd)
}

#' Evaluation of normalized time-varying spectral density function (for MCMC algorithm)
#' @importFrom stats plogis
#' @keywords internal
qpsd_dw.tilde_zigzag_cpp_expedited <- function(tilde.v, tilde.w1, tilde.w2, k1, k2, beta_basis_1_k, beta_basis_2_k) {

  p <- pFromV(plogis(tilde.v))

  m <- ncol(beta_basis_2_k)

  selector1_cpp <- findInterval(plogis(tilde.w1), (1:k1)/k1, left.open = T)
  selector2_cpp <- findInterval(plogis(tilde.w2), (1:k2)/k2, left.open = T)

  psd <- qpsd_cal_cpp_expedited(beta_basis_1_k,
                                beta_basis_2_k,
                                p,
                                integer(ncol(beta_basis_1_k)),
                                selector1_cpp,
                                selector2_cpp)


  return(psd)
}


#' Calculating log likelihood
#' @keywords internal
llike_dw <- function(FZ, norm_psd, tau) {

  # Bivariate time-varying spectral density (defined on [0, 1]^2)
  f <- tau * norm_psd

  # Whittle log-likelihood
  ll <- -sum(log(f + 1e-100) + (FZ + 1e-100) / (f + 1e-100) ) # the added 1e-100 is for numeric stability

  return(ll)
}


#' Construct Bernstein polynomial basis of degree k on omega
#' @param omega numeric vector in \eqn{[0,1]} of evaluation points
#' @param k positive integer for the degree
#' @param bernstein_l left boundary of the dilated Bernstein polynomials
#' @param bernstein_r right boundary of the dilated Bernstein polynomials
#' @return A matrix
#' @examples
#' \dontrun{
#' omega <- seq(0, 1, by = 0.01)
#' betaBasis_k_dw(omega, 100, 0.1, 0.9)
#' }
#' @keywords internal
betaBasis_k_dw <- function(omega, k, bernstein_l = 0, bernstein_r = 1) {
  N <- length(omega)
  basis <- matrix(
    dbeta(bernstein_l + omega * (bernstein_r - bernstein_l), rep(1:k, each = N), rep(k:1, each = N))
    * (bernstein_r - bernstein_l)
    / (pbeta(bernstein_r, rep(1:k, each = N), rep(k:1, each = N)) - pbeta(bernstein_l, rep(1:k, each = N), rep(k:1, each = N))),
    ncol = N,
    byrow = TRUE)

  basis
}


#' Construct Bernstein polynomial bases of degree up to kmax on omega
#' @param omega numeric vector in \eqn{[0,1]} of evaluation points
#' @param kmax positive integer for the largest degree
#' @param bernstein_l,bernstein_r left and right truncation related to the dilation
#' @return A list of length kmax, where the k-th list element is a matrix containing the polynomial basis of degree k
#' @examples
#' \dontrun{
#' omega <- seq(0, 1, by = 0.01)
#' dbList_dw_Bern(omega, 100, 0.1, 0.9)
#' }
#' @keywords internal
dbList_dw_Bern <- function(omega, kmax, bernstein_l = 0, bernstein_r = 1) {
  db.list <- vector("list", kmax)
  cat("Precomputing beta basis functions ")
  for (kk in 1:kmax) {
    db.list[[kk]] <- betaBasis_k_dw(omega, kk, bernstein_l, bernstein_r)
    if (!(kk%%10)) cat(".")
  }
  cat(" done!\n")
  return(db.list)
}


#' Construct Bernstein polynomial bases of degree up to kmax on omega for frequency parameter lambda
#' @param omega numeric vector in \eqn{[0,1]} of evaluation points
#' @param kmax positive integer for the largest degree
#' @param bernstein_l,bernstein_r left and right truncation related to the dilation
#' @return A list of length kmax, where the k-th list element is a matrix containing the polynomial basis of degree k
#' @examples
#' \dontrun{
#' omega <- time_grid <- seq(0, 1, by = 0.01)
#' dbList_dw_Bern_for_lambda(omega, 100, 0.1, 0.9)
#' }
#' @keywords internal
dbList_dw_Bern_for_lambda <- function(omega, kmax, bernstein_l = 0, bernstein_r = 1, m, time_grid) {
  db.list <- vector("list", kmax)
  cat("Precomputing beta basis functions ")
  for (kk in 1:kmax) {
    x <- betaBasis_k_dw(omega, kk, bernstein_l, bernstein_r)
    db.list[[kk]] <- x[ ,(time_grid - 1)%%m + 1, drop = FALSE]
    if (!(kk%%10)) cat(".")
  }
  cat(" done!\n")
  return(db.list)
}

