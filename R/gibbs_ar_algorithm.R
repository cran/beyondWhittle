#' Gibbs sampler for Bayesian AR model in PACF parametrization,
#' including support for TS to be a nuisance parameter
#' @importFrom stats rbinom
#' @keywords internal
gibbs_AR_nuisance_intern <- function(data, mcmc_params, prior_params, model_params, full_lik=F) {
  # MCMC parameters
  stopifnot(!is.null(mcmc_params$Ntotal)); stopifnot(mcmc_params$Ntotal>0)
  Ntotal <- mcmc_params$Ntotal
  stopifnot(!is.null(mcmc_params$burnin)); stopifnot(mcmc_params$burnin>=0 && mcmc_params$burnin<Ntotal)
  burnin <- mcmc_params$burnin
  stopifnot(!is.null(mcmc_params$thin)); stopifnot(mcmc_params$thin>=1)
  thin <- mcmc_params$thin
  stopifnot(!is.null(mcmc_params$print_interval)); stopifnot(mcmc_params$print_interval>0)
  print_interval <- mcmc_params$print_interval
  
  # Adaption parameters
  stopifnot(!is.null(mcmc_params$Nadaptive)); stopifnot(mcmc_params$Nadaptive>=0)
  Nadaptive <- mcmc_params$Nadaptive
  stopifnot(!is.null(mcmc_params$adaption.batchSize)); stopifnot(mcmc_params$adaption.batchSize>0)
  adaption.batchSize <- mcmc_params$adaption.batchSize
  stopifnot(!is.null(mcmc_params$adaption.targetAcceptanceRate)); stopifnot(mcmc_params$adaption.targetAcceptanceRate>0 && mcmc_params$adaption.targetAcceptanceRate<1)
  adaption.targetAcceptanceRate <- mcmc_params$adaption.targetAcceptanceRate
  
  # Prior parameters
  stopifnot(!is.null(prior_params$ar.order))
  stopifnot(prior_params$ar.order >= 0)
  ar.order <- prior_params$ar.order
  stopifnot(!is.null(prior_params$rho.alpha) && !is.null(prior_params$rho.beta))
  if (ar.order > 0) {
    stopifnot(prior_params$rho.alpha>0 && prior_params$rho.beta>0)    
  }
  stopifnot(length(prior_params$rho.alpha)==prior_params$ar.order && length(prior_params$rho.beta)==prior_params$ar.order)
  rho.alpha <- prior_params$rho.alpha
  rho.beta <- prior_params$rho.beta
  stopifnot(!is.null(prior_params$sigma2.alpha) && !is.null(prior_params$sigma2.beta)); stopifnot(prior_params$sigma2.alpha > 0 && prior_params$sigma2.beta > 0)
  sigma2.alpha <- prior_params$sigma2.alpha
  sigma2.beta <- prior_params$sigma2.beta
  
  # Model paramaters
  stopifnot(!is.null(model_params$theta_dim)); stopifnot(model_params$theta_dim >= 0)
  theta_dim <- model_params$theta_dim
  stopifnot(!is.null(model_params$get_noise)); stopifnot(class(model_params$get_noise)=="function")
  get_noise <- model_params$get_noise
  stopifnot(!is.null(model_params$initialize_theta)); stopifnot(class(model_params$initialize_theta)=="function")
  initialize_theta <- model_params$initialize_theta
  stopifnot(!is.null(model_params$lprior_theta)); stopifnot(class(model_params$lprior_theta)=="function")
  lprior_theta <- model_params$lprior_theta
  stopifnot(!is.null(model_params$propose_next_theta)); stopifnot(class(model_params$propose_next_theta)=="function")
  propose_next_theta <- model_params$propose_next_theta
  # Note: model_params$excludeBoundary not really needed -- the AR approach
  # does not use the frequency domain representation of the data
  # --> Only used to decide on mean centering the data or not
  stopifnot(!is.null(model_params$excludeBoundary))
  excludeBoundary <- model_params$excludeBoundary
  n <- length(data)

  # handle missing values: Initialize
  NA_pos <- which(is.na(data))
  num_NA <- length(NA_pos)
  has_NA <- (num_NA > 0)
  missingValues_trace <- matrix(NA, nrow=Ntotal, ncol=num_NA)
  missingValues_trace[1,] <- rep(0, num_NA)
  data[NA_pos] <- missingValues_trace[1,]
  if (has_NA) {
    cat("NOTE: Missing values at the following time positions:", 
        NA_pos, "treating them as random.\n")
  }
  
  # Initialize storage
  omega <- omegaFreq(n)
  lambda <- pi*omega
  sigma2_trace <- lpostTrace <- deviance <- rep(NA, Ntotal)
  rho_trace <- matrix(NA, nrow=ar.order, ncol=Ntotal)
  theta_trace <- matrix(NA, nrow=theta_dim, ncol=Ntotal)
  theta_trace[,1] <- initialize_theta(data, NULL)
  if (ar.order > 0) {
    ar_start <- ar(data, order.max=ar.order, aic=F)
    sigma2_trace[1] <- ar_start$var.pred
    rho_trace[,1] <- ARMAacf(ar=ar_start$ar, pacf=T)
  } else {
    sigma2_trace[1] <- var(data)
  }
  lsd.prop.rho <- rep(-log(n), ar.order)/2
  var.prop.rho <- exp(2*lsd.prop.rho)
  
  # Metropolis-within-Gibbs sampler
  tim0 <- proc.time()
  for (j in 2:Ntotal) {
    
    if (j%%print_interval == 0) {
      tim <- proc.time()
      print_mcmc_state(j, Ntotal, tim, tim0)
    }
    
    ##
    ## 00: Initialize and adaptive MCMC
    ##
    if (j==2) {
      # Marginal posterior of time series (without prior of non-nuisance part)
      noise <- get_noise(data, theta_trace[,j-1])  
      a <- pacf_to_ar(rho_trace[,j-1]) #pacf2AR(rho_trace[,j-1])[ar.order,]
      f_param <- unrollPsd(psd_arma(lambda, a, numeric(0), sigma2_trace[j-1]), n)
      f.store <- lpost_AR(noise, 
                          rho_trace[,j-1], rho.alpha, rho.beta,
                          sigma2_trace[j-1], sigma2.alpha, sigma2.beta,
                          full_lik)
      lpostTrace[j-1] <- f.store + lprior_theta(theta_trace[,j-1])
      deviance[j-1] <- -2*llike_AR(noise, rho_trace[,j-1], sigma2_trace[j-1], full_lik=T)
    }
    if ((ar.order > 0) && (j < Nadaptive) && (j > 1) 
        && (j %% adaption.batchSize == 1)) {
      batch <- (j-adaption.batchSize):(j-1)
      adaption.delta <- min(0.1, 1/(j^(1/3))) # c.f. Rosenthal
      batch.rho <- rho_trace[, batch, drop=F]
      stopifnot(class(batch.rho)=="matrix")
      batch.rho.acceptanceRate <- apply(batch.rho, 1, acceptanceRate) 
      lsd.prop.rho <- lsd.prop.rho + ((batch.rho.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      var.prop.rho <- exp(2*lsd.prop.rho)
    }
    
    ##
    ## 01: Draw rho_i's (double step with theta)
    ##
    if (ar.order > 0) {
      rho_prev <- rho_trace[,j-1]
      theta_prev <- theta_trace[,j-1]
      for (pp in 1:ar.order) {
        rejected <- F
        rho_star <- rho_prev
        if (rbinom(1,1,0.75)) {
          rho_star[pp] <- rho_trace[pp,j-1] + rnorm(1, 0, sqrt(var.prop.rho[pp]))
          if (abs(rho_star[pp]) >= 1) rejected <- T
        } else {
          # fallback proposal
          rho_star[pp] <- runif(1,-1,1)
        }
        if (!rejected) {
          a_star <- pacf2AR(rho_star)[ar.order,]
          f_param_star <- unrollPsd(psd_arma(lambda, a_star, numeric(0), sigma2_trace[j-1]), n)
          theta_star <- propose_next_theta(data=data, f=f_param_star, previous_theta=theta_prev, NULL)
          noise_star <- get_noise(data, theta_star)  
          f.rho.star <- lpost_AR(noise_star, 
                                 rho_star, rho.alpha, rho.beta,
                                 sigma2_trace[j-1], sigma2.alpha, sigma2.beta,
                                 full_lik)
          f.rho <- f.store
          alpha.rho <- min(0, f.rho.star + lprior_theta(theta_star) -
                             f.rho - lprior_theta(theta_prev))
          rejected <- (log(runif(1,0,1)) >= alpha.rho)
        }
        if (!rejected) {
          # accept
          rho_trace[pp,j] <- rho_prev[pp] <- rho_star[pp]
          theta_trace[,j] <- theta_prev <- theta_star
          noise <- noise_star
          a <- a_star
          f_param <- f_param_star
        } else {
          # Reject and use previous
          rho_trace[pp,j] <- rho_prev[pp]
          theta_trace[,j] <- theta_prev
        }
      }      
    }
    
    ##
    ## 02: Draw theta
    ##
    if (theta_dim > 0) {
      theta_star <- propose_next_theta(data=data, f=f_param, previous_theta=theta_prev, NULL)
      noise_star <- get_noise(data, theta_star)  
      f.theta.star <- lpost_AR(noise_star, 
                             rho_trace[,j], rho.alpha, rho.beta,
                             sigma2_trace[j-1], sigma2.alpha, sigma2.beta,
                             full_lik)
      f.theta <- f.store
      # Note: theta_prev defined, from previous double step (since ar.order>0)
      alpha.theta <- min(0, f.theta.star + lprior_theta(theta_star) -
                           f.theta - lprior_theta(theta_prev)) 
      rejected <- (log(runif(1,0,1)) >= alpha.theta)
      if (!rejected) {
        # accept
        theta_trace[,j] <- theta_prev <- theta_star
        noise <- noise_star
      } else {
        # Reject and use previous
        theta_trace[,j] <- theta_prev
      }
    }
    
    ##
    ## 03: sigma2 from conjugate
    ##
    eps <- genEpsARMAC(noise, a, numeric(0)) 
    sigma2_trace[j] <- 1 / rgamma(1, shape=sigma2.alpha + (n-1)/2,
                                  rate=sigma2.beta+sum(eps^2)/2)
    f_param <- unrollPsd(psd_arma(lambda, a, numeric(0), sigma2_trace[j]), n)
    f.store <- lpost_AR(noise, 
                        rho_trace[,j], rho.alpha, rho.beta,
                        sigma2_trace[j], sigma2.alpha, sigma2.beta,
                        full_lik)
    lpostTrace[j] <- f.store + lprior_theta(theta_trace[,j])
    deviance[j] <- -2*llike_AR(noise, rho_trace[,j], sigma2_trace[j], full_lik=T)
    
    ##
    ## handle missing values: draw with MH
    ##
    if (has_NA) {
      sigma2_NA_prop <- 4 * pi * f_param[1]
      for (iii in seq_len(length(NA_pos))) {
        i <- NA_pos[iii]
        missingValues_j_star <- data[i] + sqrt(sigma2_NA_prop) * rt(1, 4)
        data_star <- data; data_star[i] <- missingValues_j_star
        data_star <- data_star - mean(data_star)
        noise_star <- get_noise(data_star, theta_trace[,j])  # noise = data - signal
        f.NA <- f.store
        f.NA.star <- lpost_AR(noise_star, 
                              rho_trace[,j], rho.alpha, rho.beta,
                              sigma2_trace[j], sigma2.alpha, sigma2.beta,
                              full_lik)
        alphaNA <- min(0, f.NA.star - f.NA)
        if (log(runif(1,0,1)) < alphaNA) {
          missingValues_trace[j,iii] <- missingValues_j_star
          data <- data_star
          noise <- noise_star
          f.store <- f.NA.star
        } else {
          missingValues_trace[j,iii] <- missingValues_trace[j-1,iii]
        }
      }
    }
  }
  
  ##
  ## Post processing
  ##
  keep <- seq(from=burnin+1, to=Ntotal, by=thin)
  theta_trace <- theta_trace[, keep, drop=F]
  sigma2_trace <- sigma2_trace[keep]
  rho_trace <- rho_trace[, keep, drop=F]
  lpostTrace <- lpostTrace[keep]
  deviance <- deviance[keep]
  missingValues_trace <- missingValues_trace[keep,,drop=F]
  colnames(missingValues_trace) <- NA_pos
  
  # DIC
  DIC.FIT <- mean(deviance)
  DIC.ENP <- var(deviance) / 2
  DIC <- list(DIC=DIC.FIT+DIC.ENP, DIC.ENP=DIC.ENP)
  
  # Construct spectral estimates and credible regions
  N_MCMC_IT <- length(keep)
  N <- length(lambda)
  fpsd_store <- log_fpsd_store <- array(data=NA, dim=c(N_MCMC_IT, N))
  for (i in 1:N_MCMC_IT) {
    if (ar.order == 0) {
      ar_store <- numeric(0)
    }
    if (ar.order == 1) {
      ar_store <- rho_trace[,i]
    }
    if (ar.order > 1) {
      ar_store <- pacf2AR(rho_trace[,i])[ar.order,]
    }
    sigma2_store <- sigma2_trace[i]
    fpsd_store[i,] <- psd_arma(lambda,ar_store,numeric(0),sigma2_store)
    log_fpsd_store[i,] <- logfuller(fpsd_store[i,])
  }
  fpsd.s <- apply(fpsd_store, 2, median)
  fpsd.mean <- apply(fpsd_store, 2, mean)
  fpsd.s05 <- apply(fpsd_store, 2, quantile, 0.05)
  fpsd.s95 <- apply(fpsd_store, 2, quantile, 0.95)
  log.fpsd.s <- apply(log_fpsd_store, 2, median)
  log.fpsd.mad <- apply(log_fpsd_store, 2, mad)
  log.fpsd.help <- apply(log_fpsd_store, 2, uniformmax)
  log.Cvalue <- quantile(log.fpsd.help, 0.9)
  log.conflower <- exp(log.fpsd.s - log.Cvalue * log.fpsd.mad)
  log.confupper <- exp(log.fpsd.s + log.Cvalue * log.fpsd.mad)
  
  list(psi=rho_trace,
       sigma2=sigma2_trace,
       theta=theta_trace,
       ar.order=ar.order,
       DIC=DIC,
       fpsd.s=fpsd.s,
       fpsd.s05=fpsd.s05,
       fpsd.s95=fpsd.s95,
       log.conflower=log.conflower,
       log.confupper=log.confupper,
       missingValues_trace=missingValues_trace,
       lpostTrace=lpostTrace)
}
