#' Gibbs sampling algorithm for VAR model
#' @importFrom utils tail
#' @importFrom stats rWishart
#' @keywords internal
gibbs_VAR_nuisance_intern <- function(data,
                               mcmc_params,
                               prior_params,
                               model_params,
                               full_lik=F) {
  K <- ncol(data)
  n <- nrow(data)
  
  # MCMC parameters
  stopifnot(!is.null(mcmc_params$Ntotal)); stopifnot(mcmc_params$Ntotal>0)
  Ntotal <- mcmc_params$Ntotal
  stopifnot(!is.null(mcmc_params$burnin)); stopifnot(mcmc_params$burnin>=0 && mcmc_params$burnin<Ntotal)
  burnin <- mcmc_params$burnin
  stopifnot(!is.null(mcmc_params$thin)); stopifnot(mcmc_params$thin>=1)
  thin <- mcmc_params$thin
  stopifnot(!is.null(mcmc_params$print_interval)); stopifnot(mcmc_params$print_interval>0)
  print_interval <- mcmc_params$print_interval
  # Note: No adaptive MCMC employed here
  
  # Prior parameters
  stopifnot(!is.null(prior_params$var.order))
  p <- prior_params$var.order
  stopifnot(!is.null(prior_params$beta_prior))
  stopifnot(length(prior_params$beta_prior)==K*K*p)
  beta_prior <- prior_params$beta_prior
  stopifnot(!is.null(prior_params$V_prior))
  stopifnot(class(prior_params$V_prior)=="matrix")
  stopifnot(ncol(prior_params$V_prior)==K*K*p && nrow(prior_params$V_prior)==K*K*p)
  if (p>0) stopifnot(is_spd(prior_params$V_prior))
  V_prior <- prior_params$V_prior
  stopifnot(!is.null(prior_params$S_prior))
  stopifnot(class(prior_params$S_prior)=="matrix")
  stopifnot(ncol(prior_params$S_prior)==K && nrow(prior_params$S_prior)==K)
  stopifnot(is_spd(prior_params$S_prior))
  S_prior <- prior_params$S_prior
  stopifnot(!is.null(prior_params$nu_prior))
  stopifnot(prior_params$nu_prior >= 0) # Note that prior is improper for 0 <= nu <= K-1
  nu_prior <- prior_params$nu_prior

  # Model paramaters
  stopifnot(!is.null(model_params$theta_dim)); stopifnot(model_params$theta_dim >= 0)
  theta_dim <- model_params$theta_dim
  stopifnot(!is.null(model_params$get_noise)); stopifnot(class(model_params$get_noise)=="function")
  get_noise <- model_params$get_noise
  stopifnot(!is.null(model_params$get_data)); stopifnot(class(model_params$get_data)=="function")
  get_data <- model_params$get_data
  stopifnot(!is.null(model_params$initialize_theta)); stopifnot(class(model_params$initialize_theta)=="function")
  initialize_theta <- model_params$initialize_theta
  stopifnot(!is.null(model_params$lprior_theta)); stopifnot(class(model_params$lprior_theta)=="function")
  lprior_theta <- model_params$lprior_theta
  stopifnot(!is.null(model_params$propose_next_theta)); stopifnot(class(model_params$propose_next_theta)=="function")
  propose_next_theta <- model_params$propose_next_theta

  TT <- n-p 
  if (p>0) {
    V_prior_inv <- solve(V_prior)
  } else {
    V_prior_inv <- matrix(0,0,0)
  }
  S_prior_inv <- solve(S_prior)
  
  # handle missing values: Initialize
  NA_tim <- which(apply(data, 1, function(z) any(is.na(z))))
  NA_pos <- which(is.na(data))
  num_NA <- length(NA_pos)
  has_NA <- (num_NA > 0)
  missingValues_trace <- matrix(NA, nrow=Ntotal, ncol=num_NA)
  missingValues_trace[1,] <- rep(0, num_NA)
  data[NA_pos] <- missingValues_trace[1,]
  if (has_NA) {
    cat("NOTE: Missing values at the following time points:", 
        NA_tim, "treating them as random.\n")
  }
  
  # Initiate values
  lambda <- pi*omegaFreq(n)
  N <- length(lambda)
  Sigma_trace <- array(NA, dim=c(K,K,Ntotal))
  Sigma_inv_trace <- array(NA, dim=c(K,K,Ntotal)) # redundant, kept for convenience
  beta_trace <- array(NA, dim=c(K*K*p, Ntotal))
  phi_trace <- array(NA, dim=c(K,K*p,Ntotal)) # redundant, kept for convenience
  # real valued {f11, Re(f12), Im(f12), f22}-representation here:
  psd_trace <- array(NA, dim=c(K,K,N,Ntotal)) # redundant, kept for convenience 
  theta_trace <- matrix(NA, nrow=theta_dim, ncol=Ntotal)
  theta_trace[,1] <- initialize_theta(data, NULL)
  noise <- get_noise(data, theta_trace[,1])
  Y_mat <- apply(noise, 2, tail, n-p) 
  Y_vec <- c(t(Y_mat))
  ZZ <- VAR_regressor_matrix(noise, p)
  
  lpost_trace <- rep(NA, Ntotal)
  if (p>0) {
    a1 <- ar(noise,order.max=p,aic=F,method="ols")
    Sigma_trace[,,1] <- a1$var.pred
    beta_trace[,1] <- c(a1$ar)
  } else {
    Sigma_trace[,,1] <- var(noise)
  }
  Sigma_inv_trace[,,1] <- solve(Sigma_trace[,,1])
  phi_trace[,,1] <- phiFromBeta_normalInverseWishart(beta_trace[,1], K, p)
  f_param <- psd_varma(lambda, ar=phi_trace[,,1], Sigma=Sigma_trace[,,1])
  # Internally, we store the PSDs in a non-redundant real-valued representation,
  # to facilitate the postprocessing (computation of quantiles)
  psd_trace[,,,1] <- realValuedPsd(f_param)
  lpost_trace[1] <- llike_var(noise, phi_trace[,,1], Sigma_trace[,,1], full_lik) +
    lprior_theta(theta_trace[,1])
  
  # Main MCMC loop
  tim0 <- proc.time()
  for (isample in 1:(Ntotal-1)) {
    if ((isample+1)%%print_interval == 0) {
      tim <- proc.time()
      print_mcmc_state(isample+1, Ntotal, tim, tim0)
    }

    # 1) Full conditional of beta
    if (p>0) {
      V_post_inv <- V_prior_inv
      beta_post <- V_prior_inv %*% beta_prior
      for (i in 1:TT) {
        Z_t <- ZZ[((i-1)*K+1):(i*K),]
        Y_t <- Y_vec[((i-1)*K+1):(i*K)]
        V_post_inv <- V_post_inv + t(Z_t) %*% Sigma_inv_trace[,,isample] %*% Z_t
        beta_post <- beta_post + t(Z_t) %*% Sigma_inv_trace[,,isample] %*% Y_t
      }
      V_post <- solve(V_post_inv)
      beta_post <- V_post %*% beta_post
      beta_trace[,isample+1] <- MASS::mvrnorm(1, mu=beta_post, Sigma=V_post)
      phi_trace[,,isample+1] <- phiFromBeta_normalInverseWishart(beta_trace[,isample+1], K, p)
    }
    f_param <- psd_varma(lambda, ar=phi_trace[,,isample+1], Sigma=Sigma_trace[,,isample])
    
    # 2) Draw theta (MH)
    if (theta_dim > 0) {
      theta_prev <- theta_trace[,isample]
      theta_prop <- propose_next_theta(data=data, f=f_param, previous_theta=theta_prev, NULL)
      theta_star <- theta_prop$theta_star
      
      noise_star <- get_noise(data, theta_star)
      f.theta_star <- llike_var(noise_star, 
                                phi_trace[,,isample+1], 
                                Sigma_trace[,,isample],
                                full_lik) +
        lprior_theta(theta_star)
      f.theta <- llike_var(noise, 
                           phi_trace[,,isample+1], 
                           Sigma_trace[,,isample],
                           full_lik) +
        lprior_theta(theta_trace[,isample])
      alpha.theta <- min(0, f.theta_star - f.theta + 
                           theta_prop$lprop_previous_theta - 
                           theta_prop$lprop_theta_star)
      if (log(runif(1,0,1)) < alpha.theta) {
        # accept
        theta_trace[,isample+1] <- theta_star
        noise <- noise_star
        Y_mat <- apply(noise, 2, tail, n-p) 
        Y_vec <- c(t(Y_mat))
        ZZ <- VAR_regressor_matrix(noise, p)
      } else {
        # reject and use previous
        theta_trace[,isample+1] <- theta_trace[,isample]
      }
    }
    
    # 3) Full conditional of Sigma_inv
    nu_post <- TT + nu_prior
    S_post <- S_prior
    for (i in 1:TT) {
      Z_t <- ZZ[((i-1)*K+1):(i*K),]
      Y_t <- Y_vec[((i-1)*K+1):(i*K)]
      ymZb <- Y_t - Z_t %*% beta_trace[,isample+1]
      S_post <- S_post + ymZb %*% t(ymZb)
    }
    S_post_inv <- solve(S_post)
    Sigma_inv_trace[,,isample+1] <- rWishart(n=1, df=nu_post, Sigma=S_post_inv)[,,1]
    Sigma_tmp <- solve(Sigma_inv_trace[,,isample+1])
    Sigma_trace[,,isample+1] <- Sigma_tmp

    # Update traces
    psd_trace[,,,isample+1] <- realValuedPsd(f_param)
    stopifnot(!any(is.na(psd_trace[,,,isample+1])))
    lpost_trace[isample+1] <- llike_var(noise, phi_trace[,,isample+1], 
                                        Sigma_trace[,,isample+1], full_lik) +
      lprior_theta(theta_trace[,isample+1])
    
    ##
    ## handle missing values: draw with MH
    ##
    if (has_NA) {
      sigma2_NA_prop <- 4 * pi * norm(f_param[,,1], type="2")
      for (iii in seq_len(length(NA_pos))) {
        i <- NA_pos[iii]
        missingValues_j_star <- data[i] + sqrt(sigma2_NA_prop) * rt(1, 4)
        data_star <- data; data_star[i] <- missingValues_j_star
        data_star <- data_star - mean(data_star)
        noise_star <- get_noise(data_star, theta_trace[,isample+1])  # noise = data - signal
        f.NA <- lpost_trace[isample+1]
        f.NA.star <- llike_var(noise_star, phi_trace[,,isample+1], 
                               Sigma_trace[,,isample+1], full_lik) +
          lprior_theta(theta_trace[,isample+1])
        alphaNA <- min(0, f.NA.star - f.NA)
        if (log(runif(1,0,1)) < alphaNA) {
          missingValues_trace[isample+1,iii] <- missingValues_j_star
          data <- data_star
          noise <- noise_star
          lpost_trace[isample+1] <- f.NA.star
        } else {
          missingValues_trace[isample+1,iii] <- missingValues_trace[isample,iii]
        }
      }
    }
  }
  
  # Remove burnin-period and return traces
  keep <- seq(from=burnin+1, to=Ntotal, by=thin)
  fpsd.s <- apply(psd_trace[,,,keep], c(1,2,3), median)
  fpsd.mean <- apply(psd_trace[,,,keep], c(1,2,3), mean)
  fpsd.s05 <- apply(psd_trace[,,,keep], c(1,2,3), quantile, 0.05)
  fpsd.s95 <- apply(psd_trace[,,,keep], c(1,2,3), quantile, 0.95)
  missingValues_trace <- missingValues_trace[keep,,drop=F]
  colnames(missingValues_trace) <- missingValues_str_help(NA_pos, n)
  
  # Uniform credible intervals
  alpha_uci <- 0.1
  uci_tmp <- uci_matrix(fpsd.sample=psd_trace[,,,keep], 
                        alpha=alpha_uci)
  fpsd.uci05 <- uci_tmp$fpsd.uci05
  fpsd.uci95 <- uci_tmp$fpsd.uci95
  uuci_tmp <- uci_matrix(fpsd.sample=psd_trace[,,,keep], 
                         alpha=alpha_uci, 
                         uniform_among_components=T)
  fpsd.uuci05 <- uuci_tmp$fpsd.uci05
  fpsd.uuci95 <- uuci_tmp$fpsd.uci95
  
  return(list(beta=beta_trace[,keep],
              phi=phi_trace[,,keep],
              Sigma=Sigma_trace[,,keep],
              Sigma_inv=Sigma_inv_trace[,,keep],
              fpsd.s=fpsd.s,
              fpsd.mean=fpsd.mean,
              fpsd.s05=fpsd.s05,
              fpsd.s95=fpsd.s95,
              fpsd.uci05=fpsd.uuci05,
              fpsd.uci95=fpsd.uuci95,
              missingValues_trace=missingValues_trace,
              lpost=lpost_trace, 
              p=p,
              theta=theta_trace[,keep]))
}
