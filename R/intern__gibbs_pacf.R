#library(Rcpp)
#library(ltsa) # tacvfARMA
# source("gibbs_util.R")

#' Negative log likelihood values for scree-type plots
#'
#' (Approximate) negative maximum log-likelihood for for different autoregressive orders to produce scree-type plots.
#' @details By default, the maximum likelihood is approximated by the Yule-Walker method, due to numerical stabililty and computational speed. Further details can be found in the simulation study section in the referenced paper.
#' @param data numeric vector of data
#' @param order.max maximum autoregressive order to consider
#' @param method character string giving the method used to fit the model, to be forwarded to \code{stats::\link{ar}}
#' @return a data frame containing the autoregressive orders \code{p} and the corresponding negative log likelihood values \code{nll}
#' @references C. Kirch et al. (2017)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' <arXiv:1701.04846>
#' @examples 
#' \dontrun{
#' 
#' ###
#' ### Interactive visual inspection for the sunspot data
#' ###
#' 
#' data <- sqrt(as.numeric(sunspot.year))
#' data <- data <- data - mean(data)
#' 
#' screeType <- ar_screeType(data, order.max=15)
#' 
#' # Determine the autoregressive order by an interactive visual inspection of the scree-type plot
#' plot(x=screeType$p, y=screeType$nll, type="b")
#' p_ind <- identify(x=screeType$p, y=screeType$nll, n=1, labels=screeType$p)
#' print(screeType$p[p_ind])
#' }
#' @export
ar_screeType <- function(data, order.max, method="yw") {
  stopifnot(order.max > 0)
  if (abs(mean(data)) > 1e-4) {
    data <- data - mean(data)
    warning("Data has been mean centered for your convenience")
  }
  p_vals <- 0:order.max
  nll_val <- rep(NA, length(p_vals))
  for (p in p_vals) {
    if (p==0) {
      nll_val[p+1] <- -sum(dnorm(data, mean=0, sd=sd(data), log=T))
    }
    if (p > 0) {
      ar_p <- ar(data, aic=F, order.max=p)
      if (length(ar_p$ar) != p) {
        stop(paste0("Something went wrong with the AR estimation for p=", p))
      }
      nll_val[p+1] <- -ar_lik(data, ar_p$ar, mean=0, sd=sqrt(ar_p$var.pred), log=T)
    }
  }
  return(data.frame(p=p_vals, nll=nll_val))
}

#' Full likelihood of an autoregressive time series model with i.i.d. normal innovations
#' @param x numeric vector of data
#' @param ar vector of ar parameters
#' @param mean the innovation mean
#' @param sd the innovation standard deviation
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @return numeric value for the likelihood or log-likelihood
#' @export
ar_lik <- function(x, ar, mean=0, sd=1, log=F) {
  stopifnot(length(x) > 0)
  stopifnot(length(mean) == 1)
  stopifnot(length(sd) == 1)
  epsilon_t <- genEpsARMAC(x, ar=ar, ma=numeric(0))
  cll <- sum(dnorm(epsilon_t, mean=mean, sd=sd, log=T)) # conditional part
  p <- length(ar)
  if (length(p) > 0) {
    zt_p <- head(x, p)
    gamma_p <- ltsa::tacvfARMA(phi=ar, maxLag=p-1, sigma2=sqrt(sd))
    Gamma_p <- acvMatrix(gamma_p)
    Gamma_p_inv <- solve(Gamma_p)
    mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
    mll <- -log(2*pi)*p/2 -log(det(Gamma_p)) / 2 + mll_unnormalized # marginal part
  } else {
    mll <- 0
  }
  res <- cll + mll
  if (!log) {
    res <- exp(res)
  }
  return(res)
}

#' Conditional Gaussian likelihood (good enough for sampling steps)
#' @keywords internal
llike_pacf <- function(data, psi, sigma2) {
  ar <- pacf2AR(psi)[length(psi),]
  epsilon_t <- genEpsARMAC(data, ar=ar, ma=numeric(0))
  cll <- sum(dnorm(epsilon_t, mean=0, sd=sqrt(sigma2), log=T))
  return(cll)
}
#' Full Gaussian likelihood (includes marginals; only needed for model comparison)
#' @keywords internal
llike_full <- function(data, psi, sigma2) {
  ar <- pacf2AR(psi)[length(psi),]
  epsilon_t <- genEpsARMAC(data, ar=ar, ma=numeric(0))
  cll <- sum(dnorm(epsilon_t, mean=0, sd=sqrt(sigma2), log=T))
  p <- length(ar)
  zt_p <- head(data, p)
  gamma_p <- ltsa::tacvfARMA(phi=ar, maxLag=p-1, sigma2=sigma2)
  Gamma_p <- acvMatrix(gamma_p)
  Gamma_p_inv <- solve(Gamma_p)
  mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
  mll <- -log(2*pi)*p/2 -log(det(Gamma_p)) / 2 + mll_unnormalized
  return(cll + mll)
}
#' Log-Prior of PACF parametrization
#' @keywords internal
lprior_pacf <- function(psi, sigma2, sigma2.alpha, sigma2.beta, psi.alpha, psi.beta) {
  logprior <- (-sigma2.alpha - 1) * log(sigma2) - sigma2.beta / sigma2 +
    sum(dbeta((psi+1)/2, psi.alpha, psi.beta, log=T))
  return(logprior)
}
#' Log-posterior
#' @keywords internal
lpost_pacf <- function(data, psi, sigma2, sigma2.alpha, sigma2.beta, psi.alpha, psi.beta) {
  lprior_pacf(psi, sigma2, sigma2.alpha, sigma2.beta, psi.alpha, psi.beta) +
    llike_pacf(data, psi, sigma2)
}
#' Gibbs sampler of AR model with PACF parametrization
#' @keywords internal
gibbs_pacf <- function(data,
                       Ntotal,
                       burnin,
                       Nadaptive=burnin,
                       adaption.batchSize=50,
                       adaption.targetAcceptanceRate=0.44,
                       ar.order,
                       sigma2.alpha=0.001,
                       sigma2.beta=0.001,
                       mu.prop=rep(0,ar.order),
                       var.prop=rep(1/length(data),ar.order), # starting value according to asymptotic theory of ML estimate
                       psi.alpha=rep(1,ar.order),  # defaults to uniform prior on psi
                       psi.beta=rep(1,ar.order)) { # defaults to uniform prior on psi
  stopifnot(ar.order >= 0)
  stopifnot(length(mu.prop) == ar.order)
  stopifnot(length(var.prop) == ar.order)
  stopifnot(length(psi.alpha) == ar.order)
  stopifnot(length(psi.beta) == ar.order)

  n <- length(data)

  if (ar.order == 0) {
    # Dummy AR(1) return model with zero coefficient (for computational simplicity)
    psi <- rep(0, Ntotal - burnin)
    sigma2 <- 1 / rgamma(Ntotal - burnin, sigma2.alpha + n / 2, sigma2.beta + sum(data^2) / 2)
    deviance <- rep(NA, Ntotal - burnin)
    for (isample in seq_len(Ntotal - burnin)) {
      deviance[isample] <- -2*llike_full(data, 0, sigma2[isample]) # include marginal likelihood
    }
    DIC.FIT <- mean(deviance)
    # DIC.ENP <- mean(deviance) + 2*llike_full(data, 0, mean(sigma2)) # effective number of parameters
    DIC.ENP <- var(deviance) / 2
    DIC <- list(DIC=DIC.FIT+DIC.ENP, DIC.ENP=DIC.ENP)
    return(list(psi=psi,
                sigma2=sigma2,
                ar.order=ar.order,
                DIC=DIC))
  }

  psi <- matrix(NA, nrow=ar.order, ncol=Ntotal)
  sigma2 <- rep(NA, Ntotal)
  deviance <- rep(NA, Ntotal)

  sigma2[1] <- var(data)
  psi[, 1] <- mu.prop
  deviance[1] <- -2*llike_full(data, psi[, 1], sigma2[1])

  sd.prop <- sqrt(var.prop)
  lsd.prop <- log(sd.prop)
  for (i in 1:(Ntotal-1)) {

    if (i%%200 == 0) {
      cat("iteration ", i, "\n")
    }
    
    ## Adaption step: Adjust propsal variance
    if ((i < Nadaptive) && (i > 1) && (i %% adaption.batchSize == 1)) {
      adaption.delta <- min(0.1, 1/(i^(1/3))) # c.f. Rosenthal
      batch.psi <- psi[, (i-adaption.batchSize):(i-1)]
      if (class(batch.psi)=="numeric") { # one psi param
        batch.psi.acceptanceRate <- acceptanceRate(batch.psi)
      } else { # several psi params
        stopifnot(class(batch.psi)=="matrix")
        batch.psi.acceptanceRate <- apply(batch.psi, 1, acceptanceRate)
      }
      lsd.prop <- lsd.prop + ((batch.psi.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      sd.prop <- exp(lsd.prop)
    }

    ## First, sample PACF with RWMWG
    psi.old <- psi[, i]
    for (p in 1:ar.order) {
      psi.star <- psi.old
      psi.star[p] <- psi.old[p] + rnorm(1, 0, sd.prop[p])

      ## Accept // Reject
      if (abs(psi.star[p]) < 1) {
        f.psi.old <- lpost_pacf(data,
                                psi.old,
                                sigma2[i],
                                sigma2.alpha,
                                sigma2.beta,
                                psi.alpha,
                                psi.beta)
        f.psi.star <- lpost_pacf(data,
                           psi.star,
                           sigma2[i],
                           sigma2.alpha,
                           sigma2.beta,
                           psi.alpha,
                           psi.beta)
        alpha.psi <- min(0, f.psi.star -
                           f.psi.old) # Note: Symmetric (normal) density in MH step
        if (log(runif(1,0,1)) < alpha.psi) {
          psi[p, i+1] <- psi.star[p] # Accept proposal
          psi.old[p] <- psi.star[p] # "old" for next Gibbs step
        } else {
          psi[p, i+1] <- psi.old[p] # Reject and use previous
        }
      } else {
        psi[p, i+1] <- psi.old[p] # Reject and use previous
      }
    }

    ## Now, sample variance sigma2
    epsilon_t <- genEpsARMAC(data, ar=pacf2AR(psi[,i+1])[ar.order,], ma=numeric(0))
    sigma2[i+1] <- 1 / rgamma(1, sigma2.alpha + (n - ar.order) / 2, sigma2.beta + sum(epsilon_t^2) / 2)
#     gamma_p <- tacvfARMA(phi=pacf2AR(psi[,i+1])[ar.order,], maxLag=p-1, sigma2=1) # Note the unit variance
#     Gamma_p <- acvMatrix(gamma_p)
#     Gamma_p_inv <- solve(Gamma_p)
#     zt_p <- head(data, p)
#     sigma2[i+1] <- 1 / rgamma(1, sigma2.alpha + n / 2,
#                               sigma2.beta + (t(zt_p) %*% Gamma_p_inv %*% zt_p + sum(epsilon_t^2)) / 2)
    deviance[i+1] <- -2*llike_full(data, psi[, i+1], sigma2[i+1]) # include marginal likelihood
  }

  ## DIC
  deviance <- deviance[-(1:burnin)]
  if (ar.order == 1) {
    psi.mean <- mean(psi)
  }
  if (ar.order > 1) {
    psi.mean <- apply(psi, 1, mean)
  }
  DIC.FIT <- mean(deviance)
  # DIC.ENP <- mean(deviance) + 2*llike_full(data, psi.mean, mean(sigma2)) # effective number of parameters
  DIC.ENP <- var(deviance) / 2
  DIC <- list(DIC=DIC.FIT+DIC.ENP, DIC.ENP=DIC.ENP)

  list(psi=psi[,-(1:burnin)],
       sigma2=sigma2[-(1:burnin)],
       ar.order=ar.order,
       DIC=DIC)
}
