#' Simulate from a Multivariate Normal Distribution
#' 
#' Produces one or more samples from the specified multivariate normal distribution.
#' @details This is a simple wrapper function based on \link[MASS]{mvrnorm},
#' to be used within \link[beyondWhittle]{sim_varma}
#' @param n sample size
#' @param d dimensionality
#' @param mu mean vector
#' @param Sigma covariance matrix
#' @param ... further arguments to be parsed to 
#' @return If n=1 a vector of length d, otherwise an n by d matrix with one sample in each row.
#' @export
rmvnorm <- function(n, d, mu=rep(0,d), Sigma=diag(d), ...) {
  MASS::mvrnorm(n, mu, Sigma, ...)
}


#' Simulate from a VARMA model
#' 
#' Simulate from a Vector Autoregressive Moving Average (VARMA) model.
#' Note that no test for model stationarity is performed.
#' @param model A list with component \code{ar} and/or \code{ma} giving the VAR and VMA 
#' coefficients respectively. An empty list gives an VARMA(0, 0) model, that is white noise.
#' @param n sample size
#' @param d positive integer for the dimensionality
#' @param rand.gen random vector generator, function of type rand.gen(n, d, ...)
#' @param burnin length of burnin period (initial samples that are discarded)
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @seealso \link[stats]{arima.sim} to simulate from univariate ARMA models
#' @return If n=1 a vector of length d, otherwise an n by d matrix with one sample in each row.
#' @examples 
#' \dontrun{
#' # Example: Draw from bivariate normal VAR(2) model
#' ar <- rbind(c(.5, 0, 0, 0), c(0, -.3, 0, -.5))
#' Sigma <- matrix(data=c(1, .9, .9, 1), nrow=2, ncol=2)
#' x <- sim_varma(n=256, d=2, model=list(ar=ar))
#' plot.ts(x)
#' }
#' @export
sim_varma <- function(model, n, d, rand.gen=rmvnorm,
                      burnin=1e4, ...) {
  if (!is.list(model)) {
    stop("'model' must be list")
  }
  if (is.null(model$ar)) {
    ar <- matrix(nrow=d, ncol=0)
  } else {
    if (d==1) {
      stopifnot(is.vector(model$ar))
      ar <- matrix(model$ar, nrow=1, ncol=length(model$ar))
    } else {
      stopifnot(nrow(model$ar) == d && !(ncol(model$ar) %% d))
      ar <- model$ar
    }
  }
  if (is.null(model$ma)) {
    ma <- matrix(nrow=d, ncol=0)
  } else {
    if (d==1) {
      stopifnot(is.vector(model$ma))
    } else {
      stopifnot(nrow(model$ma) == d && !(ncol(model$ma) %% d))
      ma <- model$ma
    }
  }
  p <- ncol(ar) / d
  q <- ncol(ma) / d
  stopifnot(burnin >= max(p,q))
  X_sim <- epsilon_sim <- matrix(NA, nrow=n+burnin, ncol=d)
  if (max(p,q) > 0) {
    X_sim[1:max(p,q),] <- rand.gen(n=max(p,q), d, ...)
    epsilon_sim[1:max(p,q),] <- rand.gen(n=max(p,q), d, ...)
  }
  for (t in (max(p,q)+1):nrow(X_sim)) {
    epsilon_sim[t,] <- rand.gen(n=1, d, ...)
    X_sim[t,] <- epsilon_sim[t,]
    for (j in seq_len(p)) {
      X_sim[t,] <- X_sim[t,] + t(ar[,((j-1)*d+1):(j*d)] %*% 
                                   matrix(X_sim[t-j,], nrow=d, ncol=1))
    }
    for (j in seq_len(q)) {
      X_sim[t,] <- X_sim[t,] + t(ma[,((j-1)*d+1):(j*d)] %*% 
                                   matrix(epsilon_sim[t-j,], nrow=d, ncol=1))
    }
  }
  if (burnin > 0) {
    return(X_sim[-(1:burnin),])
  } else {
    return(X_sim)
  }
}

#' Convert vector parametrization (beta) to matrix-parametrization (phi),
#' the latter as e.g. used in MTS::VAR()$ar
#' @param beta coefficient vector, of dimension K*d*d
#' @param K positive integer, vector dimensionality
#' @param p nonnegarive integer, VAR order
#' @return K times K*p coefficient matrix
#' @keywords internal
phiFromBeta_normalInverseWishart <- function(beta, K, p) {
  return(t(matrix(data=beta, nrow=K*p, ncol=K)))
}

#' VAR regressor matrix, see Section 2.2.3 in Koop and Korobilis (2010)
#' @keywords internal
VAR_regressor_matrix <- function(data, var.order) {
  K <- ncol(data)
  n <- nrow(data)
  p <- var.order
  if (p>0) {
    ZZ <- NULL 
    for (t in (p+1):n) { 
      z_mt <- c(t(data[(t-1):(t-p),]))
      for (m in 1:K) {
        ZZ <- rbind(ZZ, c(
          rep(0, (m-1)*K*p),
          z_mt,
          rep(0, (K-m)*K*p)))
      }
    }
  } else {
    ZZ <- matrix(0, n*K, 0)
  }
  ZZ
}

#' VARMA(p,q) spectral density function
#' 
#' Evaluate the VARMA(p,q) spectral density at some frequencies freq in [0,pi).
#' Note that no test for model stationarity is performed.
#' @details See section 11.5 in the referenced book
#' @param freq numeric vector of frequencies to evaluate the psd, 0 <= freq < pi
#' @param ar autoregressive coeffient matrix (d times p*d) of VARMA model, defaults to empty VAR component
#' @param ma moving average coeffient matrix (d times p*d) of VARMA model, defaults to empty VAR component
#' @param Sigma positive definite innovation covariance matrix (d times d)
#' @references P. J. Brockwell and R. Davis (1996)
#' \emph{Time Series: Theory and Methods (Second Edition)}
#' @return an array containing the values of the varma psd matrix at freq
#' @export
psd_varma <- function(freq, 
                      ar=matrix(nrow=nrow(Sigma),ncol=0), 
                      ma=matrix(nrow=nrow(Sigma),ncol=0), 
                      Sigma) {
  psd_varma_help(freq, ar=ar, ma=ma, sigma=Sigma)$psd
}

#' helping function for psd_varma
#' @keywords internal
psd_varma_help <- function(freq, 
                           ar=matrix(nrow=nrow(sigma),ncol=0), 
                           ma=matrix(nrow=nrow(sigma),ncol=0), 
                           sigma) {
  d <- nrow(sigma)
  N <- length(freq)
  stopifnot(nrow(ar)==d && !(ncol(ar)%%d))
  stopifnot(nrow(ma)==d && !(ncol(ma)%%d))
  transfer_ar <- transfer_polynomial(freq, -ar) # note the minus
  transfer_ma <- transfer_polynomial(freq, ma)
  psd <- varma_transfer2psd(transfer_ar, transfer_ma, sigma)
  return(list(psd=psd,
              transfer_ar=transfer_ar,
              transfer_ma=transfer_ma))
}

#' VAR(p) likelihood
#' @keywords internal
llike_var <- function(zt, ar, sigma, full_lik) {
  if (full_lik) {
    llike_var_full(zt, ar, sigma)
  } else {
    llike_var_partial(zt, ar, sigma)
  }
}

#' VAR(p) partial likelihood (unnormalized)
#' Note: Fine for fixed p, but not suited for model comparison
#' @keywords internal
llike_var_partial <- function(zt, ar, sigma) {
  epsilon_t <- epsilon_var(zt, ar)
  ll <- sldmvnorm(epsilon_t, sigma)
  return(ll)
}

#' VAR(p) full likelihood
#' @keywords internal
llike_var_full <- function(zt, ar, sigma) {
  n <- nrow(zt)
  d <- ncol(zt)
  p <- ncol(ar) / nrow(ar)
  epsilon_t <- epsilon_var(zt, ar)
  cll <- sldmvnorm(epsilon_t, sigma)
  if (p>0) {
    zt_p <- c(t(zt[1:p,]))
    gamma_p <- VARMAcov_muted(Phi=ar, Sigma=sigma, lag=p-1)$autocov[,(1:(d*p))]
    Gamma_p <- acvBlockMatrix(gamma_p)
    Gamma_p_inv <- solve(Gamma_p)
    mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
    mll <- -log(det(Gamma_p)) / 2 + mll_unnormalized
  } else {
    mll <- 0
  }
  return(cll + mll)
}

#' Psi-weight calculation for a VARMA model.
#' NOTE: This is an exact copy of the MTS::PSIwgt function
#' (only with the plot functionality removed, as not needed).
#' This has to be done because the MTS package has been removed
#' from CRAN in April 2022.
#' @keywords internal
PSIwgt <- function(Phi=NULL,Theta=NULL,lag=12,plot=TRUE,output=FALSE){
  ### Compute the psi-weight matrices of a VARMA(p,q) model,
  #### Phi=[phi1, phi2, ..., phip]
  #### Theta=[theta1,theta2,...,thetaq]
  #### Sigma= residual covariance matrix
  q=0; p=0; k=0
  if(length(Theta) > 0){
    k=dim(Theta)[1]
    k1=dim(Theta)[2]
    q=floor(k1/k)
  }
  #
  if(length(Phi) > 0){
    k=dim(Phi)[1]
    k1=dim(Phi)[2]
    p=floor(k1/k)
  }
  #
  if(k < 1) k=1
  PSI=diag(k); WGT=c(PSI)
  #
  for (il in 1:lag){
    ilk=il*k
    tmp=matrix(0,k,k)
    if((q > 0) && (il <= q))tmp=-Theta[,(ilk-k+1):ilk]
    if(p > 0){
      iend=min(il,p)
      for (j in 1:iend){
        jdx=(il-j)
        kdx=j*k
        tmp=tmp+Phi[,(kdx-k+1):kdx]%*%PSI[,(jdx*k+1):(jdx*k+k)]
      }
      ## end  of p > 0.
    }
    PSI=cbind(PSI,tmp)
    WGT=cbind(WGT,c(tmp))
    ### end of il-loop
  }
  ## print the output if needed
  if(output){
    for (i in 1:lag){
      cat("Lag: ",i," psi-matrix","\n")
      ist=i*k
      print(round(PSI[,(ist+1):(ist+k)],5))
    }
    ## end print
  }
  ## plots the psi-weights
  if(plot){
    # Note: plotting functionality has been removed 
    # for compatibility reasons.
  }
  PSIwgt <- list(psi.weight=PSI,irf=WGT)
}

#' This is a nearly exact copy of the MTS::VARMAcov function, where 
#' the output commands at the end are removed.
#' This has to be done because the function is called repeatedly
#' within the MCMC algorithm.
#' For future versions of the package, a better solution is intended.
#' @keywords internal
VARMAcov_muted <- function (Phi = NULL, Theta = NULL, 
                            Sigma = NULL, lag = 12, trun = 120) { 
  m1 = PSIwgt(Phi = Phi, Theta = Theta, lag = trun, plot = FALSE)
  Psi = m1$psi.weight
  nc = dim(Psi)[2]
  k = dim(Psi)[1]
  if (is.null(Sigma)) {
    wk = Psi
  }
  else {
    wk = NULL
    for (i in 0:trun) {
      ist = i * k
      wk = cbind(wk, Psi[, (ist + 1):(ist + k)] %*% Sigma)
    }
  }
  Gam0 = wk %*% t(Psi)
  SE = diag(1/sqrt(diag(Gam0)))
  covmtx = Gam0
  cormtx = SE %*% Gam0 %*% SE
  for (i in 1:lag) {
    ist = i * k
    Gami = wk[, (ist + 1):nc] %*% t(Psi[, 1:(nc - ist)])
    covmtx = cbind(covmtx, Gami)
    cormtx = cbind(cormtx, SE %*% Gami %*% SE)
  }
  VARMAcov <- list(autocov = covmtx, ccm = cormtx)
}