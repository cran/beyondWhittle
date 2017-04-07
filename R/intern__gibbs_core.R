#library(Rcpp)
#library(ltsa)

#' Function that calculates correction matrix C_n ^ (-1/2) for ARMA model
#' @keywords internal
Cn <- function(n, g) {

  # n: Sample size
  # f: Estimated PSD using Choudhuri et al. (2004) method
  # ar: Estimated AR parameters
  # ma: Estimated MA parameters
  # sigma2: Estimated sigma^2

  # g = f / f_param

  # Angular frequencies [0, pi]
  pin <- 2 * pi / n
  lambda <- pin * (1:(n / 2 + 1) - 1)
  #lambda <- 2 * pi * (1:(n / 2 + 1) - 1) / n

  # Diagonals of correction matrix
  # Note the hard-coded psd_arma function here, and power of -0.5.
  cc <- rep(NA, n)
  cc[1] <- g[1]
  cc[n] <- g[n]
  cc[2 * 1:(n / 2 - 1)] <- cc[2 * 1:(n / 2 - 1) + 1] <- g[2 * 1:(n / 2 - 1) + 1]

  cc <- 1 / sqrt(cc)  # Scale correction matrix

  return(cc)  # Don't need a diagonal matrix with FFT implementation

}
# Compute a PSD in the Bernstein-Dirichlet parametrization.
qpsd <- function(omega, v, w, k, db.list, epsilon=1e-20) {
  p <- pFromV(v)
  weight <- mixtureWeight(p, w, k)
  psd <- densityMixture(weight, db.list[[k]])
  psd <- pmax(psd, epsilon)
  return(list(psd = psd,
              weight = weight,
              p=p))  # *** Do we need to output weight? ***
}

#' Deprecated, not in use.
#' @keywords internal
arma_partial <- function(zt, ar, ma, nll_fun, ...) {

  stop("function deprecated (need to invlove sigma2)")
    epsilon_t <- genEpsARMAC(zt, ar, ma)

    # Conditional log likelihood for ARMA(p, q)
    cll <- nll_fun(epsilon_t, ...)

    return(cll)
}

#' Function for computing the ARMA(p, q) conditional likelihood
#' Models with AR component: zt are from t = p + 1 to m
#' Models with MA component: eps_0 = ... = eps_{q+1} = 0
#' The input zt should be FCFZ, i.e., data corrected in freq. domain then IFTd
#' FCFZ should be FCFZ[-c(1, n)]; i.e., have removed elements 1 and n beforehand
#' No sigma2 in here anymore!
#' ...: Further arguments to be passed to likelihood function
#' @keywords internal
arma_conditional <- function(zt, ar, ma, nll_fun, sigma2, ...) {
  stopifnot(sigma2 == 1)

  # Conditional Log Likelihood cll
  epsilon_t <- genEpsARMAC(zt, ar, ma)
  cll <- -nll_fun(epsilon_t, ...)

  # Marginal Log Likelihood mll (Gaussian!)
  p <- length(ar)
  if (p == 0) {
    mll <- 0
  } else {
    zt_p <- zt[1:p] #head(zt, p)
    gamma_p <- ltsa::tacvfARMA(phi=ar, maxLag=p-1, sigma2=sigma2)
    Gamma_p <- acvMatrix(gamma_p)
    Gamma_p_inv <- solve(Gamma_p)
    mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
    mll <-  -log(det(Gamma_p)) / 2 + mll_unnormalized -log(2*pi)*p/2
  }

  return(-cll-mll) # return negative log likelihood
}

#' Log corrected parametric likelihood
#' @keywords internal
llike <- function(omega, FZ, ar, ma, v, w, k, tau, corrected, toggle.q, pdgrm, db.list, dist, nll_fun, ex.kurt, beta, f.alpha) {
  
  # Calculates Whittle or corrected log-likelihood (assume n even) for Gaussian errors
  n <- length(FZ)
  
  # Un-normalised PSD (defined on [0, 1])
  q.psd <- qpsd(omega, v, w, k, db.list)$psd
  q <- unrollPsd(q.psd, n)
  
  # Normalised PSD (defined on [0, pi])
  f <- tau * q
  #f_param <- psd_arma(pi*omega, ar, ma,sigma2.ex)
  
  # do not consider the following frequencies in the likelihood,
  # as they correspond to the mean (or alternating mean)
  if (n %% 2) {
    excludedFrequecies <- c(1,n)
  } else {
    excludedFrequecies <- 1
  }
  
  # Corrected log-likelihood
  if (corrected == TRUE) {
    
    #     if (dist=="student" && ex.kurt != 0) {
    #       df <- 6 / ex.kurt + 4
    #       sigma2 <- df / (df - 2) # sigma2 of likelihood (!=1 only possible for t data, have to correct it here)
    #     } else {
    #       sigma2 <- 1
    #     }
    sigma2 <- 1
    
    if (toggle.q) {
      C <- sqrt(unrollPsd(psd_arma(pi*omega,ar,ma,1),n)^(1-f.alpha) / f) 
    } else {
      C <- sqrt(unrollPsd(psd_arma(pi*omega,ar,ma,1),n) / f) #Cn(n, f)
    }
    
    # Input for ARMA parametric likelihood - Inverse FT
    FCFZ <- fast_ift(C * FZ)
    
    # Calculate ARMA parametric conditional likelihood
    if (dist == "generalized") {
      ex.kurt <- beta # QUICK HACK: Abuse ex.kurt to parse beta to nll_fun
    }
    p.arma <- arma_conditional(FCFZ, ar, ma, nll_fun, sigma2, ex.kurt)
    
    # Corrected Gaussian log-likelihood 
    # llike <- sum(log(C[-c(1,n)])) - p.arma   # Note: The minus sign here.
    llike <- sum(log(C[-excludedFrequecies])) - p.arma   # Note: The minus sign here.
    
  }
  
  # Whittle log-likelihood
  if (corrected == FALSE) {
    llike <- -sum(log(f[-excludedFrequecies] * 2 * pi) + pdgrm[-excludedFrequecies] / (f[-excludedFrequecies] * 2 * pi))
    llike <- llike / 2
  } 
  
  return(llike)
  
}

#' Log corrected posterior
#' @keywords internal
lpost <- function(omega, FZ, ar, ma, v, w, k, tau,
                  M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                  corrected, toggle.q, pdgrm, db.list,
                  dist, nll_fun, ex.kurt, kurt.lambda, beta, beta.alpha, beta.beta, f.alpha, rho, rho.alpha, rho.beta) {

  # Unnormalised log posterior
  ll <- llike(omega, FZ, ar, ma, v, w, k, tau, corrected, toggle.q, pdgrm, db.list, dist, nll_fun, ex.kurt, beta, f.alpha)
  lp <- lprior(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta, dist, ex.kurt, kurt.lambda, beta, beta.alpha, beta.beta, f.alpha, corrected, rho, rho.alpha, rho.beta)

  if (is.na(ll)) {
    ll_params <- list(omega=omega, 
                      FZ=FZ, 
                      ar=ar, 
                      ma=ma, 
                      v=v, 
                      w=w, 
                      k=k,
                      tau=tau, 
                      corrected=corrected, 
                      toggle.q=toggle.q, 
                      pdgrm=pdgrm, 
                      dist=dist, 
                      nll_fun=nll_fun, 
                      ex.kurt=ex.kurt, 
                      beta=beta, 
                      f.alpha=f.alpha)
    TXT <- paste("Likelihood evaluated as NA. Parameters: ", ll_params)
    stop(TXT)
  }
  if (is.na(lp)) {
    lp_params <- list(v=v, 
                      w=w, 
                      k=k, 
                      tau=tau, 
                      M=M, 
                      g0.alpha=g0.alpha, 
                      g0.beta=g0.beta, 
                      k.theta=k.theta, 
                      tau.alpha=tau.alpha, 
                      tau.beta=tau.beta, 
                      dist=dist, 
                      ex.kurt=ex.kurt, 
                      kurt.lambda=kurt.lambda, 
                      beta=beta, 
                      beta.alpha=beta.alpha, 
                      beta.beta=beta.beta, 
                      f.alpha=f.alpha, 
                      corrected=corrected, 
                      rho=rho, 
                      rho.alpha=rho.alpha, 
                      rho.beta=rho.beta)
    TXT <- paste("Prior evaluated as NA. Parameters: ", lp_params)
    stop(TXT)
  }
  return(ll+lp)
}

#' Log prior of Bernstein-Dirichlet mixture and parametric working model -- all unnormalized
#' @keywords internal
lprior <- function(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                   dist, ex.kurt, kurt.lambda, beta, beta.alpha, beta.beta, f.alpha, corrected,
                   rho, rho.alpha, rho.beta) {

  # log joint prior - all unnormalised
  # Hyperparameters are M, g0.a, g0.b, k.theta, tau.alpha, tau.beta
  # NOTE: Flat prior on f.alpha

  logprior <- (M - 1) * sum(log(1 - v)) +  # log prior for V's - beta(1, M)
    sum((g0.alpha - 1) * log(w) + (g0.beta - 1) * log(1 - w)) -  # log prior for Z's - beta(a, b)
    k.theta * k * log(k) -   # log prior for k
    (tau.alpha + 1) * log(tau) - tau.beta / tau # log prior for tau (Inverse Gamma)

  if (dist=="student") {
    logprior <- logprior - kurt.lambda * ex.kurt # log prior for ex.kurt (Exponential)
  }
  if (dist=="generalized") {
    logprior <- logprior + dgamma(beta, beta.alpha, beta.beta, log=T) # log prior for beta (gamma)
  }

  # Beta prior on PACF
  if (corrected && !is.null(rho)) {
    logprior <- logprior + sum(dbeta((rho+1)/2, rho.alpha, rho.beta, log=T))
  }

  return(logprior)

}
