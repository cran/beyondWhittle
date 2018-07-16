#' Compute a PSD in the Bernstein-Dirichlet parametrization.
#' @keywords internal
qpsd <- function(omega, v, w, k, beta_basis_k, epsilon=1e-20) {
  p <- pFromV(v)
  weight <- mixtureWeight(p, w, k)
  psd <- densityMixture(weight, beta_basis_k)
  psd <- pmax(psd, epsilon)
  return(list(psd = psd,
              weight = weight,
              p=p))  # *** Do we need to output weight? ***
}

#' Function for computing the ARMA(p, q) conditional likelihood
#' Models with AR component: zt are from t = p + 1 to m
#' Models with MA component: eps_0 = ... = eps_{q+1} = 0
#' The input zt should be FCFZ, i.e., data corrected in freq. domain then IFTd
#' FCFZ should be FCFZ[-c(1, n)]; i.e., have removed elements 1 and n beforehand
#' No sigma2 in here anymore!
#' ...: Further arguments to be passed to likelihood function
#' @keywords internal
arma_conditional <- function(zt, ar, ma, nll_fun, full_lik, ...) {
  # sigma2 <- 1
  # 
  # # Conditional Log Likelihood cll
  # epsilon_t <- genEpsARMAC(zt, ar, ma)
  # cll <- -nll_fun(epsilon_t, ...)
  # 
  # # Marginal Log Likelihood mll (Gaussian!)
  # p <- length(ar)
  # #if (p == 0) {
  #   mll <- 0
  # #} else {
  # #  zt_p <- zt[1:p] #head(zt, p)
  # #  gamma_p <- tacvfARMA(phi=ar, maxLag=p-1, sigma2=sigma2)
  # #  Gamma_p <- acvMatrix(gamma_p)
  # #  Gamma_p_inv <- solve(Gamma_p)
  # #  mll_unnormalized <- as.numeric(-1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p)
  # #  mll <-  -log(det(Gamma_p)) / 2 + mll_unnormalized -log(2*pi)*p/2
  # #}
  # 
  # return(-cll-mll) # return negative log likelihood
  
  if (!(is.null(ma) || length(ma)==0)) {
    stop("MA component not supported yet in likelihood")
  }
  sigma2 <- 1
  n <- length(zt)
  p <- length(ar)
  eps <- genEpsARC(zt, ar)
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

#' Log corrected parametric likelihood
#' @keywords internal
llike <- function(omega, FZ, ar, ma, v, w, k, tau, corrected, toggle.q, 
                  pdgrm, beta_basis_k, nll_fun, f.alpha, excludeBoundary, full_lik) {
  
  # Calculates Whittle or corrected log-likelihood (assume n even) for Gaussian errors
  n <- length(FZ)
  
  if (!(n %% 2)) {
    boundaryFrequecies <- c(1,n)
  } else {
    boundaryFrequecies <- 1
  }
  
  # Un-normalised PSD (defined on [0, 1])
  q.psd <- qpsd(omega, v, w, k, beta_basis_k)$psd
  q <- unrollPsd(q.psd, n)
  
  # Normalised PSD (defined on [0, pi])
  f <- tau * q
  #f_param <- psd_arma(pi*omega, ar, ma,sigma2.ex)
  
  
  # Corrected log-likelihood
  if (corrected) {
    
    if (toggle.q) {
      C <- sqrt(unrollPsd(psd_arma(pi*omega,ar,ma,1),n)^(1-f.alpha) / f) 
    } else {
      C <- sqrt(unrollPsd(psd_arma(pi*omega,ar,ma,1),n) / f) #Cn(n, f)
    }
    
    if (excludeBoundary) {
      C[boundaryFrequecies] <- 0
    }
    
    # Input for ARMA parametric likelihood - Inverse FT
    FCFZ <- fast_ift(C * FZ)
    
    p.arma <- arma_conditional(FCFZ, ar, ma, nll_fun, full_lik)
    
    # Corrected Gaussian log-likelihood 
    # llike <- sum(log(C[-c(1,n)])) - p.arma   # Note: The minus sign here.
    
    if (excludeBoundary) {
      llike <- sum(log(C[-boundaryFrequecies])) - p.arma   # Note: The minus sign here.
    } else {
      llike <- sum(log(C)) - p.arma   # Note: The minus sign here.
    }
    
    
  }
  
  # Whittle log-likelihood
  if (corrected == FALSE) {
    
    pdgrm_scaling <- c(pi, rep(2*pi, n-1))
    if (!(n%%2)) pdgrm_scaling[n] <- pi
    
    if (excludeBoundary) {
      llike <- -sum(log(f[-boundaryFrequecies] * 2 * pi) + pdgrm[-boundaryFrequecies] / (f[-boundaryFrequecies] * 2*pi))
    } else {
      llike <- -sum(log(f * pdgrm_scaling) + pdgrm / (f * pdgrm_scaling))
    }
    llike <- llike / 2
  } 
  
  return(llike)
  
}

#' Log corrected posterior
#' @keywords internal
lpost <- function(omega, FZ, ar, ma, v, w, k, tau, 
                  M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                  corrected, toggle.q, pdgrm, beta_basis_k,
                  nll_fun, f.alpha, rho, rho.alpha, rho.beta, 
                  excludeBoundary, full_lik) {

  # Unnormalised log posterior
  ll <- llike(omega, FZ, ar, ma, v, w, k, tau, corrected, toggle.q, pdgrm, beta_basis_k, nll_fun, f.alpha, excludeBoundary, full_lik)
  lp <- lprior(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta, f.alpha, corrected, rho, rho.alpha, rho.beta)

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
                      nll_fun=nll_fun, 
                      f.alpha=f.alpha,
                      excludeBoundary=excludeBoundary,
                      full_lik=full_lik)
    TXT <- paste("Likelihood evaluated as NA. Parameters: ", ll_params)
    print(TXT)
    save(list="ll_params", file="beyondWhittle_llike_debug")
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
                      f.alpha=f.alpha, 
                      corrected=corrected, 
                      rho=rho, 
                      rho.alpha=rho.alpha, 
                      rho.beta=rho.beta)
    TXT <- paste("Prior evaluated as NA. Parameters: ", lp_params)
    print(TXT)
    save(list="lp_params", file="beyondWhittle_lprior_debug")
    stop(TXT)
  }
  return(ll+lp)
}

#' Log prior of Bernstein-Dirichlet mixture and parametric working model -- all unnormalized
#' @keywords internal
lprior <- function(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta, 
                   f.alpha, corrected, 
                   rho, rho.alpha, rho.beta) {
  
  # log joint prior - all unnormalised
  # Hyperparameters are M, g0.a, g0.b, k.theta, tau.alpha, tau.beta
  # NOTE: Flat prior on f.alpha
  
  logprior <- (M - 1) * sum(log(1 - v)) +  # log prior for V's - beta(1, M)
    sum((g0.alpha - 1) * log(w) + (g0.beta - 1) * log(1 - w)) -  # log prior for Z's - beta(a, b)
    k.theta * k * log(k) -   # log prior for k  
    (tau.alpha + 1) * log(tau) - tau.beta / tau # log prior for tau (Inverse Gamma)
  
  # Beta prior on PACF
  if (corrected && !is.null(rho)) {
    logprior <- logprior + sum(dbeta((rho+1)/2, rho.alpha, rho.beta, log=T))
  }
  
  return(logprior)
  
}
