#' Generate a random samples from a Dirichlet distribution
#' @param alpha numeric vector of positive concentration parameter
#' @return a vector of the same length as alpha
#' @keywords internal
my_rdirichlet <- function(alpha) {
  ga <- rgamma(length(alpha), alpha, rep(1,length(alpha)))
  ga / sum(ga)
}

#' Log determinant of stick breaking transformation V -> p
#' @details See Section 3 in Choudhuri et al (2004)
#' @keywords internal
logDet_stickBreaking <- function(v) {
  L <- length(v)
  sum(log(1-v)*((L-1):0))
}

#' Compute normalized PSD in the Bernstein-Dirichlet parametrization.
#' @details See (5) in Choudhuri et al (2004)
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


#' Log corrected parametric AR likelihood (Gaussian)
#' @details See (5) in Kirch et al (2018)
#' @keywords internal
llike <- function(omega, FZ, ar, ma, v, w, k, tau, corrected, toggle.q, 
                  pdgrm, beta_basis_k, nll_fun, f.alpha, excludeBoundary, full_lik) {
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

#' Log posterior = log prior + log corrected parametric likelihood
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
#' @details See Section 3 in Kirch et al (2018).
#' Hyperparameters are M, g0.a, g0.b, k.theta, tau.alpha, tau.beta.
#' Note: Flat prior on f.alpha.
#' @keywords internal
lprior <- function(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta, 
                   f.alpha, corrected, 
                   rho, rho.alpha, rho.beta) {
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
