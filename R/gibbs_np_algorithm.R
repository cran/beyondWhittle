#' Gibbs sampler for corrected parametric likelihood + Bernstein-Dirichlet mixture,
#' including possibility of using time series as mere nuisance parameter
#' @importFrom stats ARMAacf
#' @importFrom stats rbinom
#' @keywords internal
gibbs_nuisance <- function(data, 
                           mcmc_params, 
                           corrected, 
                           prior_params,
                           model_params,
                           full_lik=NULL) # full_lik only in use for NP 
{

  # Note: "Assuming symmetric theta proposals in MH steps"

  # MCMC parameters
  stopifnot(!is.null(mcmc_params$Ntotal)); stopifnot(mcmc_params$Ntotal>0)
  Ntotal <- mcmc_params$Ntotal
  stopifnot(!is.null(mcmc_params$burnin)); stopifnot(mcmc_params$burnin>=0 && mcmc_params$burnin<Ntotal)
  burnin <- mcmc_params$burnin
  stopifnot(!is.null(mcmc_params$thin)); stopifnot(mcmc_params$thin>=1)
  thin <- mcmc_params$thin
  stopifnot(!is.null(mcmc_params$print_interval)); stopifnot(mcmc_params$print_interval>0)
  print_interval <- mcmc_params$print_interval
  stopifnot(!is.null(mcmc_params$numerical_thresh)); stopifnot(mcmc_params$numerical_thresh>0)
  NUMERICAL_THRESH <- mcmc_params$numerical_thresh
  #
  # Adaption mcmc parameters: See corrected && toggle branch
  #
  
  # PRIOR PAREMETERS
  stopifnot(!is.null(prior_params$M)); stopifnot(prior_params$M > 0)
  M <- prior_params$M
  stopifnot(!is.null(prior_params$g0.alpha) && !is.null(prior_params$g0.beta)); stopifnot(prior_params$g0.alpha > 0 && prior_params$g0.beta > 0)
  g0.alpha <- prior_params$g0.alpha
  g0.beta <- prior_params$g0.beta
  stopifnot(!is.null(prior_params$k.theta)); stopifnot(prior_params$k.theta > 0)
  k.theta <- prior_params$k.theta
  stopifnot(!is.null(prior_params$tau.alpha) && !is.null(prior_params$tau.beta)); stopifnot(prior_params$tau.alpha > 0 && prior_params$tau.beta > 0)
  tau.alpha <- prior_params$tau.alpha
  tau.beta <- prior_params$tau.beta
  stopifnot(!is.null(prior_params$kmax)); stopifnot(prior_params$kmax > 0)
  kmax <- prior_params$kmax
  stopifnot(!is.null(prior_params$bernstein_l) && !is.null(prior_params$bernstein_r)); stopifnot(prior_params$bernstein_l >= 0 && prior_params$bernstein_r <= 1)
  bernstein_l <- prior_params$bernstein_l
  bernstein_r <- prior_params$bernstein_r
  stopifnot(!is.null(prior_params$bernstein_coars))
  bernstein_coars <- prior_params$bernstein_coars
  stopifnot(!is.null(prior_params$L)); stopifnot(prior_params$L > 0)
  L <- prior_params$L
  if (corrected) {
    stopifnot(!is.null(prior_params$toggle))
    toggle <- prior_params$toggle
    stopifnot(!is.null(prior_params$prior.q))
    prior.q <- prior_params$prior.q
    stopifnot(!is.null(prior_params$ar.order)); stopifnot(prior_params$ar.order>0)
    ar.order <- prior_params$ar.order
    stopifnot(!is.null(prior_params$alpha.toggle))
    alpha.toggle <- prior_params$alpha.toggle
    if (toggle) {
      AR.fit <- NULL
      MA.fit <- NULL
      stopifnot(!is.null(prior_params$rho.alpha) && !is.null(prior_params$rho.beta)); stopifnot(prior_params$rho.alpha>0 && prior_params$rho.beta>0)
      stopifnot(length(prior_params$rho.alpha)==prior_params$ar.order && length(prior_params$rho.beta)==prior_params$ar.order)
      rho.alpha <- prior_params$rho.alpha
      rho.beta <- prior_params$rho.beta
      #
      # Adaption mcmc parameters
      #
      stopifnot(!is.null(mcmc_params$Nadaptive)); stopifnot(mcmc_params$Nadaptive>0)
      Nadaptive <- mcmc_params$Nadaptive
      stopifnot(!is.null(mcmc_params$adaption.batchSize)); stopifnot(mcmc_params$adaption.batchSize>0)
      adaption.batchSize <- mcmc_params$adaption.batchSize
      stopifnot(!is.null(mcmc_params$adaption.targetAcceptanceRate)); stopifnot(mcmc_params$adaption.targetAcceptanceRate>0 && mcmc_params$adaption.targetAcceptanceRate<1)
      adaption.targetAcceptanceRate <- mcmc_params$adaption.targetAcceptanceRate
      #
      #
      #
    } else {
      stopifnot(!is.null(prior_params$AR.fit) && !is.null(prior_params$MA.fit))
      AR.fit <- prior_params$AR.fit
      MA.fit <- prior_params$MA.fit
      rho.alpha <- NULL
      rho.beta <- NULL
    }
    
  } else {
    toggle <- alpha.toggle <- prior.q <- FALSE
    f.alpha <- NULL
    rho <- NULL
  }
  
  n <- length(data)
  
  # handle missing values: Initialize
  NA_pos <- which(is.na(data))
  num_NA <- length(NA_pos)
  has_NA <- (num_NA > 0)
  missingValues_trace <- matrix(NA, nrow=Ntotal, ncol=num_NA)
  missingValues_trace[1,] <- rep(0,num_NA)
  data[NA_pos] <- missingValues_trace[1,]
  if (has_NA) {
    cat("Note: Missing values at the following time positions:", 
        NA_pos, "treating them as random.\n")
  }
  
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
  stopifnot(!is.null(model_params$excludeBoundary))
  excludeBoundary <- model_params$excludeBoundary # (!is.null(model_params$excludeBoundary) && model_params$excludeBoundary)
  if (!(n %% 2)) {
    boundaryFrequecies <- c(1,n)
  } else {
    boundaryFrequecies <- 1
  }

  # Basic error messages - TO COMPLETE
  if (burnin > Ntotal) {
    stop("Burn-in parameter is larger than number of iterations")
  }
  
  if (corrected) {
    if (toggle) {
      stopifnot(ar.order > 0)
      stopifnot(ar.order==length(rho.alpha))
      stopifnot(ar.order==length(rho.beta))
    } else {
      stopifnot(!is.null(AR.fit) && !is.null(MA.fit))
      ar.order <- length(AR.fit)
    }
  } else {
    stopifnot(!toggle)
    stopifnot(!prior.q)
  }
  
  # Frequencies on unit interval
  twon <- 2 / n
  omega <- twon * (1:(n / 2 + 1) - 1)  #2 * (1:(n / 2 + 1) - 1) / n
  NN <- length(omega)
  omega_for_db_list <- seq(bernstein_l, bernstein_r, length.out=NN)
  
  # Angular frequencies on [0, pi]
  lambda <- pi * omega  #2 * pi * (1:(n / 2 + 1) - 1) / n

  if (is.finite(kmax)) {
    db.list <- dbList(n, kmax, bernstein_l, bernstein_r, 
                      bernstein_coars, verbose=T)
  }

  # Open objects for storage
  lpostTrace <- rep(NA, Ntotal)
  tau <- rep(NA, Ntotal)
  V <- matrix(NA, nrow = L, ncol = Ntotal)
  W <- matrix(NA, nrow = L + 1, ncol = Ntotal)  # Called Z in Choudhuri
  k <- rep(NA, Ntotal)
  theta <- matrix(NA, nrow=theta_dim, ncol=Ntotal)
  
  # Starting values
  tau[1] <- var(data) / (2 * pi)
  if (is.infinite(kmax)) {
    k[1] <- sample(100:1000,1)
    beta_basis_k <- betaBasis_k(omega_for_db_list, k[1], bernstein_coars)
  } else {
    k[1] <- round(kmax / 2)
    beta_basis_k <- db.list[[k[1]]]
  }
  #p_1 <- rep
  V[, 1] <- vFromP(rep(1/(L+1),L)) #rbeta(L, 1, M)
  W[, 1] <- seq(from=1/(2*k[1]), to=1-1/(2*k[1]), length.out=L+1) # rbeta(L + 1, g0.alpha, g0.beta)  # g0.alpha = g0.beta = 1 gives U[0,1]
  
  while (min(qpsd(omega,V[, 1],W[,1],k[1],beta_basis_k)$psd) <= NUMERICAL_THRESH) {
    V[, 1] <- rbeta(L, 1, M)
    W[, 1] <- rbeta(L + 1, g0.alpha, g0.beta)  # g0.alpha = g0.beta = 1 gives U[0,1]  
  }
  
  if (corrected) {
    if (toggle) {
      rho <- matrix(NA, nrow=ar.order, ncol=Ntotal) # trace of PACFs
      rho[,1] <- ARMAacf(ar(data, order.max=ar.order, aic=F)$ar, pacf=T) # Start at Yule Walker estimate
      AR <- pacf2AR(rho[,1])[ar.order,] # current AR(p) parameters
      MA <- numeric(0)
      lsd.prop.rho <- rep(-log(n), ar.order)/2
      var.prop.rho <- exp(2*lsd.prop.rho)
    } else {
      rho <- NULL
      AR <- AR.fit
      MA <- MA.fit
    }
  } else {
    rho <- NULL
    AR <- NULL
    MA <- NULL
  }
  theta[,1] <- initialize_theta(data)
  
  eps <- seq(1, L + 1) / (seq(1, L + 1) + 2 * sqrt(n))  # Metropolis proposal parameters for V and W.
  # First one only used for W0.
  
  nll_fun <- NULL # deprecated
  
  #####
  
  if (prior.q) {
    if (alpha.toggle) {
      f.alpha <- rep(NA, Ntotal)
      f.alpha[1] <- 1/2
    } else {
      f.alpha <- rep(1, Ntotal)
    }
  } else {
    stopifnot(!alpha.toggle)
    f.alpha <- NULL
  }
  
  # Metropolis-within-Gibbs sampler
  tim0 <- proc.time()
  for (i in 1:(Ntotal-1)) {
    
    if ((i+1)%%print_interval == 0) {
      tim <- proc.time()
      print_mcmc_state(i+1, Ntotal, tim, tim0)
    }
    
    # Adaption step: Adjust propsal variance 
    if (corrected && toggle) {
      if ((i < Nadaptive) && (i > 1) && (i %% adaption.batchSize == 1)) {
        batch <- (i-adaption.batchSize):(i-1)
        adaption.delta <- min(0.1, 1/(i^(1/3))) # c.f. Rosenthal
        batch.rho <- rho[, batch]
        if (class(batch.rho)=="numeric") { # one rho param
          batch.rho.acceptanceRate <- acceptanceRate(batch.rho)
        } else { # several rho params
          stopifnot(class(batch.rho)=="matrix")
          batch.rho.acceptanceRate <- apply(batch.rho, 1, acceptanceRate) 
        }
        lsd.prop.rho <- lsd.prop.rho + ((batch.rho.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
        var.prop.rho <- exp(2*lsd.prop.rho)
      }
    }
    
    noise <- get_noise(data, theta[,i])  # noise = data - signal
    FZ <- fast_ft(noise)  # Frequency domain
    pdgrm <- abs(FZ) ^ 2
    pdgrm_scaling <- c(pi, rep(2*pi, n-1))
    if (!(n%%2)) pdgrm_scaling[n] <- pi
    
    
    # STEP 1: Sample Bernstein polynomial degree k
    
    k.star <- round(rt(1, 1)) + k[i]  # Cauchy distribution discretized
    while (k.star < 1 || k.star > kmax) {  # A bit hacky
      k.star <- round(rt(1, 1)) + k[i]
    }
    if (is.infinite(kmax)) {
      beta_basis_k_star <- betaBasis_k(omega_for_db_list, k.star, bernstein_coars)
    } else {
      beta_basis_k_star <- db.list[[k.star]]
    }
    
    f.k.star <- lpost(omega,
                      FZ,
                      AR,
                      MA,
                      V[,i],
                      W[, i],
                      k.star,  # Proposed value of k
                      tau[i],
                      M,
                      g0.alpha,
                      g0.beta,
                      k.theta,
                      tau.alpha,
                      tau.beta,
                      corrected,
                      prior.q,
                      pdgrm,
                      beta_basis_k_star, # Corresponding beta basis functions
                      nll_fun,
                      f.alpha[i],
                      rho[,i],
                      rho.alpha,
                      rho.beta,
                      excludeBoundary,
                      full_lik)
    
    # log posterior of previous iteration
    f.k <- lpost(omega,
                 FZ,
                 AR,
                 MA,
                 V[, i],
                 W[, i],
                 k[i], # Old value of k
                 tau[i],
                 M,
                 g0.alpha,
                 g0.beta,
                 k.theta,
                 tau.alpha,
                 tau.beta,
                 corrected,
                 prior.q,
                 pdgrm,
                 beta_basis_k,
                 nll_fun,
                 f.alpha[i],
                 rho[,i],
                 rho.alpha,
                 rho.beta,
                 excludeBoundary,
                 full_lik)      
    
    if (i==1) {
      lpostTrace[1] <- f.k + lprior_theta(theta[,1])
    }
    
    #####
    # Accept/reject
    
    alpha1 <- min(0, f.k.star - f.k)
    
    if (log(runif(1, 0, 1)) < alpha1 && min(qpsd(omega,V[,i],W[,i],k.star,beta_basis_k_star)$psd)>NUMERICAL_THRESH) {
      k[i + 1] <- k.star  # Accept k.star
      f.store <- f.k.star
      beta_basis_k <- beta_basis_k_star
    } else {
      k[i + 1] <- k[i]  # Reject and use previous
      f.store <- f.k
    }

    
    # Step 2: Metropolis-within-Gibbs step for V (EXPENSIVE)
    for (l in 1:L) {
      
      V.star <- V.old <- V[, i]
      if (l > 1) {
        for (il in 1:(l - 1)) {
          V.star[il] <- V.old[il] <- V[il, i + 1]
        }
      }
      
      # Uniform proposal (V[,i] - eps, V[,i] + eps) on (0,1)
      V.star[l] <- runif(1, V.star[l] - eps[l], V.star[l] + eps[l])
      V.star[l][V.star[l] > 1] <- V.star[l] - 1  # Puts in [0, 1]
      V.star[l][V.star[l] < 0] <- V.star[l] + 1  # Puts in [0, 1]
      
      # log posterior for proposal
      f.V.star <- lpost(omega,
                        FZ,
                        AR,
                        MA,
                        V.star,  # V.star here
                        W[, i],
                        k[i + 1], # i + 1 here since we've already sampled it
                        tau[i],
                        M,
                        g0.alpha,
                        g0.beta,
                        k.theta,
                        tau.alpha,
                        tau.beta,
                        corrected,
                        prior.q,
                        pdgrm,
                        beta_basis_k,
                        nll_fun,
                        f.alpha[i],
                        rho[,i],
                        rho.alpha,
                        rho.beta,
                        excludeBoundary,
                        full_lik)
      
      # log posterior of previous iteration
      f.V <- f.store      
      
      # Accept/reject
      alpha2 <- min(0, f.V.star - f.V)  # log acceptance ratio
      if (log(runif(1, 0, 1)) < alpha2 && min(qpsd(omega,V.star,W[,i],k[i+1],beta_basis_k)$psd)>NUMERICAL_THRESH) {
        V[l, i + 1] <- V.star[l]  # Accept V.star
        f.store <- f.V.star
      } else {
        V[l, i + 1] <- V[l, i]  # Reject and use previous
      }
      
    }  # END: Step 2.    

    # Step 3: Metropolis-within-Gibbs step for W (EXPENSIVE)
    for (l in 1:(L + 1)) {   
      
      W.star <- W.old <- W[, i]
      if (l > 1) {
        for (il in 1:(l - 1)) {
          W.star[il] <- W.old[il] <- W[il, i + 1]
        }
      }
      
      # Uniform proposal from (W[,i] - eps, W[,i] + eps) on (0,1)
      W.star[l] <- runif(1, W.star[l] - eps[l], W.star[l] + eps[l])
      W.star[l][W.star[l] > 1] <- W.star[l] - 1  # Puts in [0, 1]
      W.star[l][W.star[l] < 0] <- W.star[l] + 1  # Puts in [0, 1]
      
      if (k[i+1] > 2 && W.star[l] <= 2/k[i+1]) {
        ##
        ## Proposal landed in left interval --> Also propose new weights
        ##
        which_left <- which(W.star <= 2/k[i+1])
        p_current <- pFromV(V[, i + 1])
        p_left <- p_current[which_left]
        weight_left <- sum(p_left)
        p_left_star <- weight_left * my_rdirichlet(rep(1,length(p_left))) # propose uniform on simplex
        p_star <- p_current
        p_star[which_left] <- p_left_star
        v_star <- vFromP(p_star[-1])
        # log posterior for proposal
        f.W.star <- lpost(omega,
                          FZ,
                          AR,
                          MA,
                          v_star,  # Adjusted in left interval here
                          W.star,  # W.star here
                          k[i + 1], # i + 1 here since already sampled
                          tau[i],
                          M,
                          g0.alpha,
                          g0.beta,
                          k.theta,
                          tau.alpha,
                          tau.beta,
                          corrected,
                          prior.q,
                          pdgrm,
                          beta_basis_k,
                          nll_fun,
                          f.alpha[i],
                          rho[,i],
                          rho.alpha,
                          rho.beta,
                          excludeBoundary,
                          full_lik)
        
        # log posterior for previous iteration
        f.W <- f.store   
        
        # Accept/reject
        alpha3 <- min(0, f.W.star - logDet_stickBreaking(v_star) - 
                        f.W + logDet_stickBreaking(V[,i+1]))  # log acceptance ratio

        if(log(runif(1, 0, 1)) < alpha3 && min(qpsd(omega,v_star,W.star,k[i+1],beta_basis_k)$psd)>NUMERICAL_THRESH) {
          W[l, i + 1] <- W.star[l]  # Accept W.star AND V.star
          V[,i+1] <- v_star
          f.store <- f.W.star
        } else {
          W[l, i + 1] <- W[l, i]  # Reject and use previous
        }
      } else {
        ##
        ## Proposal not in left interval --> MH as usual
        ##
        # log posterior for proposal
        f.W.star <- lpost(omega,
                          FZ,
                          AR,
                          MA,
                          V[, i + 1],  # i + 1 here since already sampled
                          W.star,  # W.star here
                          k[i + 1], # i + 1 here since already sampled
                          tau[i],
                          M,
                          g0.alpha,
                          g0.beta,
                          k.theta,
                          tau.alpha,
                          tau.beta,
                          corrected,
                          prior.q,
                          pdgrm,
                          beta_basis_k,
                          nll_fun,
                          f.alpha[i],
                          rho[,i],
                          rho.alpha,
                          rho.beta,
                          excludeBoundary,
                          full_lik)
        
        # log posterior for previous iteration
        f.W <- f.store   
        
        # Accept/reject
        alpha3 <- min(0, f.W.star - f.W)  # log acceptance ratio
        if(log(runif(1, 0, 1)) < alpha3 && min(qpsd(omega,V[, i + 1],W.star,k[i+1],beta_basis_k)$psd)>NUMERICAL_THRESH) {
          W[l, i + 1] <- W.star[l]  # Accept W.star
          f.store <- f.W.star
        } else {
          W[l, i + 1] <- W[l, i]  # Reject and use previous
        }
      }
    }  # END: Step 3.
    
    
    # Step 4: Directly sample tau from conjugate Inverse-Gamma density
    q.psd <- qpsd(omega, V[, i + 1], W[, i + 1], k[i + 1], beta_basis_k)$psd
    q <- unrollPsd(q.psd, n)
    
    if (corrected == FALSE) {  # For Whittle likelihood

      if (excludeBoundary) {
        tau[i + 1] <- 1 / rgamma(1, tau.alpha + (n-length(boundaryFrequecies)) / 2,
                                 tau.beta + sum(pdgrm[-boundaryFrequecies] / q[-boundaryFrequecies] / (2*pi)) / 2)
      } else {
        tau[i + 1] <- 1 / rgamma(1, tau.alpha + (n-length(boundaryFrequecies)) / 2,
                                 tau.beta + sum(pdgrm / q / pdgrm_scaling) / 2)
      }

    }  # CHECK: Should this be 2pi here or in pdgrm?
    
    if (corrected) {  # For Corrected likelihood
      
      if (prior.q) {
        B <- 1 / sqrt(q / (unrollPsd(psd_arma(lambda,AR,MA,1),n))^(1-f.alpha[i]))
      } else {
        B <- 1 / sqrt(q / unrollPsd(psd_arma(lambda,AR,MA,1),n))
      }

      if (excludeBoundary) {
        B[boundaryFrequecies] <- 0
      }
      
      # Input for ARMA parametric conditional likelihood - Inverse FFT
      FBFZ <- fast_ift(B * FZ)

      # Calculate ARMA parametric conditional likelihood
      p.arma <- 1 / 2 * sum(genEpsARMAC(FBFZ, AR, MA)^2)
      
      tau[i + 1] <- 1 / rgamma(1, tau.alpha + (n-length(boundaryFrequecies)) / 2, tau.beta + p.arma)
      
      f.store <- lpost(omega,
                       FZ,
                       AR,
                       MA,
                       V[, i + 1],  
                       W[, i + 1],  
                       k[i + 1], 
                       tau[i + 1],
                       M,
                       g0.alpha,
                       g0.beta,
                       k.theta,
                       tau.alpha,
                       tau.beta,
                       corrected,
                       prior.q,
                       pdgrm,
                       beta_basis_k,
                       nll_fun,
                       f.alpha[i],
                       rho[,i],
                       rho.alpha,
                       rho.beta,
                       excludeBoundary,
                       full_lik)
      
      #####
      # SAMPLE PACF, rho, here:
      #####
      if (toggle) {
        ##
        ## Double Gibbs step for theta and rho_i's
        ##
        rho_prev <- rho[, i]
        theta_prev <- theta[, i]
        for (pp in 1:ar.order) {
          rejected <- F
          rho_star <- rho_prev
          if (rbinom(1,1,0.75)) {
            # Proposal with learned scaling
            rho_star[pp] <- rho[pp,i] + rnorm(1, 0, sqrt(var.prop.rho[pp]))
            if (abs(rho_star[pp])>=1) rejected <- T

          } else {
            # default fallback
            rho_star[pp] <- runif(1,-1,1)
          }
          if (!rejected) {
            AR_star <- pacf2AR(rho_star)[ar.order,]
            f_param_star_tmp <- psd_arma(lambda, AR_star, MA, 1)
            if (min(f_param_star_tmp) < NUMERICAL_THRESH ||
                max(f_param_star_tmp) > 1/NUMERICAL_THRESH) {
              rejected <- T
            }
          }
          if (!rejected) {
            f_param_star <- unrollPsd(f_param_star_tmp, n)^f.alpha[i]
            f_for_theta_star <- tau[i + 1] * q * f_param_star
            theta_star <- propose_next_theta(data=data, f=f_for_theta_star, previous_theta=theta_prev, NULL)
            noise_star <- get_noise(data, theta_star)  # noise = data - signal
            FZ_star <- fast_ft(noise_star)  # Frequency domain
            pdgrm_star <- abs(FZ_star)^2
            
            f.rho.star <- lpost(omega,
                                FZ_star,
                                AR_star, 
                                MA,
                                V[, i + 1],
                                W[, i + 1],
                                k[i + 1],
                                tau[i + 1],
                                M,
                                g0.alpha,
                                g0.beta,
                                k.theta,
                                tau.alpha,
                                tau.beta,
                                corrected,
                                prior.q,
                                pdgrm_star,
                                beta_basis_k,
                                nll_fun,
                                f.alpha[i],
                                rho_star,
                                rho.alpha,
                                rho.beta,
                                excludeBoundary,
                                full_lik)
            f.rho <- f.store
            alpha_rho <- min(0, f.rho.star + lprior_theta(theta_star) - 
                               f.rho - lprior_theta(theta_prev))
            rejected <- (log(runif(1,0,1)) >= alpha_rho)
          } 
          if (!rejected) {
            # accept proposal for rho and theta
            rho[pp,i+1] <- rho_prev[pp] <- rho_star[pp]
            AR <- AR_star
            theta[,i+1] <- theta_prev <- theta_star
            noise <- noise_star
            FZ <- FZ_star
            pdgrm <- pdgrm_star
            f.store <- f.rho
          } else {
            # reject and use previous
            rho[pp,i+1] <- rho_prev[pp]
            theta[,i+1] <- theta_prev
          }
        }
      } # END: Toggle
      
      
      ### Adaption step: SAMPLE f.alpha parameter
      if (prior.q && alpha.toggle) {
        sd.alpha <- 0.1 # TODO This needs tuning
        f.alpha.star <- f.alpha[i] + rnorm(1, 0, sd.alpha)
        f.alpha.star <- max(f.alpha.star, 0)
        f.alpha.star <- min(f.alpha.star, 1)
        ff.alpha.star <- lpost(omega,
                               FZ,
                               AR,
                               MA,
                               V[, i + 1],
                               W[, i + 1],
                               k[i + 1],
                               tau[i + 1],
                               M,
                               g0.alpha,
                               g0.beta,
                               k.theta,
                               tau.alpha,
                               tau.beta,
                               corrected,
                               prior.q,
                               pdgrm,
                               beta_basis_k,
                               nll_fun,
                               f.alpha.star,
                               rho[,i+1],
                               rho.alpha,
                               rho.beta,
                               excludeBoundary,
                               full_lik)  
        ff.alpha <- f.store
        # Accept/reject
        alphaAlpha <- min(0, ff.alpha.star - ff.alpha +
                            dnorm(f.alpha[i],f.alpha.star,sd.alpha,log=TRUE) -
                            dnorm(f.alpha.star,f.alpha[i],sd.alpha,log=TRUE))
        # Note that rho proposals with abs(rho) >= 1 are rejected
        if(log(runif(1, 0, 1)) < alphaAlpha) {
          f.alpha[i + 1] <- f.alpha.star
          f.store <- ff.alpha.star
        } else {
          f.alpha[i + 1] <- f.alpha[i]  # Reject and use previous
        }
      } else {
        f.alpha[i + 1] <- f.alpha[i]
      }
    }  # END: Corrected

    
    ##
    ## NOTE: Sample parameter of interest
    ##
        
    # MH Step for theta
    if (theta_dim > 0) {
      #if (!toggle) {
      if (corrected) {
        f_param <- unrollPsd(psd_arma(lambda, AR, MA, 1), n)^f.alpha[i + 1]
        f_for_theta <- tau[i + 1] * q * f_param
      } else {
        f_for_theta <- tau[i + 1] * q
      }
      if (corrected && toggle) {
        stopifnot(!is.na(theta[,i+1]))
        previous_theta <- theta[,i+1] # Note we already sampled theta in double step
      } else {
        stopifnot(is.na(theta[,i+1]))
        previous_theta <- theta[,i]
      }
      
      theta_star <- propose_next_theta(data=data, f=f_for_theta, previous_theta=previous_theta, NULL)
      noise_star <- get_noise(data, theta_star)  # noise = data - signal
      FZ_star <- fast_ft(noise_star)  # Frequency domain
      pdgrm_star <- abs(FZ_star) ^ 2
      f.theta.star <- lpost(omega,
                            FZ_star,
                            AR, 
                            MA,
                            V[, i + 1],
                            W[, i + 1],
                            k[i + 1],
                            tau[i + 1],
                            M,
                            g0.alpha,
                            g0.beta,
                            k.theta,
                            tau.alpha,
                            tau.beta,
                            corrected,
                            prior.q,
                            pdgrm_star,
                            beta_basis_k,
                            nll_fun,
                            f.alpha[i+1],
                            rho[,i+1],
                            rho.alpha,
                            rho.beta,
                            excludeBoundary,
                            full_lik)
      f.theta <- f.store
      # Accept/reject
      # NOTE: Assuming symmetric proposal for theta
      alphaTheta <- min(0, f.theta.star + lprior_theta(theta_star) - 
                          f.theta - lprior_theta(theta[,i]))
      if (log(runif(1,0,1)) < alphaTheta) {
        theta[,i+1] <- theta_star
        f.store <- f.theta.star
      } else {
        theta[,i+1] <- theta[,i]
      }
      #}
    }
    
    ##
    ## handle missing values: draw with MH
    ##
    if (has_NA) {
      # Note: q is up to date
      if (corrected && prior.q) {
        f_param_current <- unrollPsd(psd_arma(lambda, AR, MA, 1), n)^f.alpha[i+1]
        f_current <- tau[i+1] * q * f_param_current 
      }
      else {
        f_current <- tau[i+1] * q
      }
      sigma2_NA_prop <- 4 * pi * f_current[1] 
      for (jjj in seq_len(length(NA_pos))) {
        j <- NA_pos[jjj]
        missingValues_j_star <- data[j] + sqrt(sigma2_NA_prop) * rt(1, 4)
        data_star <- data; data_star[j] <- missingValues_j_star
        data_star <- data_star - mean(data_star)
        noise_star <- get_noise(data_star, theta[,i+1])  # noise = data - signal
        FZ_star <- fast_ft(noise_star)  # Frequency domain
        pdgrm_star <- abs(FZ_star) ^ 2
        f.NA <- f.store
        f.NA.star <- lpost(omega,
                           FZ_star,
                           AR, 
                           MA,
                           V[, i + 1],
                           W[, i + 1],
                           k[i + 1],
                           tau[i + 1],
                           M,
                           g0.alpha,
                           g0.beta,
                           k.theta,
                           tau.alpha,
                           tau.beta,
                           corrected,
                           prior.q,
                           pdgrm_star,
                           beta_basis_k,
                           nll_fun,
                           f.alpha[i+1],
                           rho[,i+1],
                           rho.alpha,
                           rho.beta,
                           excludeBoundary,
                           full_lik)
        alphaNA <- min(0, f.NA.star - f.NA)
        if (log(runif(1,0,1)) < alphaNA) {
          missingValues_trace[i+1,jjj] <- missingValues_j_star
          data <- data_star
          f.store <- f.NA.star
        } else {
          missingValues_trace[i+1,jjj] <- missingValues_trace[i,jjj]
        }
      }
    }
    ##
    ##
    ##

    lpostTrace[i+1] <- f.store + lprior_theta(theta[,i+1])
  }  # END: MCMC loop
  
  # ADD WARNING IF A LOT OF SAMPLES FROM KMAX - MAY NEED TO BE LARGER
  if (is.finite(kmax) && min(abs(k-kmax), na.rm=T) < 5) {
    print_warn("Reached kmax during sampling. Consider using a higher kmax value")
  }
  
  # Which iterations to keep
  keep <- seq(burnin + 1, Ntotal, by = thin)
  lpostTrace <- lpostTrace[keep]
  k <- k[keep]
  tau <- tau[keep]
  V <- V[, keep]
  W <- W[, keep]
  if (toggle) {
    rho <- matrix(rho[, keep], nrow=ar.order)
  }
  missingValues_trace <- missingValues_trace[keep,,drop=F]
  colnames(missingValues_trace) <- NA_pos
  
  if (prior.q) {
    f.alpha <- f.alpha[keep]
  }
  theta <- theta[,keep,drop=F]
  
  fpsd.sample <- log.fpsd.sample <- matrix(NA, nrow = length(omega), ncol = length(keep))
  acov.sample <- rep(NA, length(keep))
  acor.sample <- rep(NA, length(keep))
  deviance.sample <- rep(NA, length(keep)) # ameier
  
  for (isample in 1:length(keep)) {
    if (corrected) {
      if (toggle) {
        AR <- pacf2AR(rho[,isample])[ar.order,] # TODO Store traces of AR(p) parameters, too?
        MA <- numeric(0)
      } else {
        AR <- AR.fit
        MA <- MA.fit
      }
    }
    
    if (is.infinite(kmax)) {
      beta_basis_k <- betaBasis_k(omega_for_db_list, k[isample], bernstein_coars)
    } else {
      beta_basis_k <- db.list[[k[isample]]]
    }
    
    qstore <- qpsd(omega, V[, isample], W[, isample], k[isample], beta_basis_k)
    
    if (corrected && prior.q) {
      f_param <- psd_arma(lambda, AR, MA, 1)^f.alpha[isample]
      fpsd.sample[, isample] <- tau[isample] * qstore$psd * f_param  # Multiplied by f_param
    }
    else {
      fpsd.sample[, isample] <- tau[isample] * qstore$psd
    }
    
    # Create transformed version
    log.fpsd.sample[, isample] <- logfuller(fpsd.sample[, isample])

    # Calculate autocovariance/autocorrelation *by monte carlo integration*
    weight <- qstore$weight  # Weights
    jstar <- sample(1:k[isample], size = 1000, prob = weight, replace = TRUE)  # Why 1000?
    u <- rbeta(1000, jstar, k[isample] - jstar + 1)  
    if (corrected && prior.q) {
      f_param_sample <- psd_arma(u * pi, AR, MA, 1)^f.alpha[isample]
      acov.sample[isample] <- 2 * pi * tau[isample] * 
        fast_mean(cos(u * pi) * f_param_sample) # lag 1
      avar <- 2 * pi * tau[isample] * 
        fast_mean(f_param_sample) # lag 0
      acor.sample[isample] <- acov.sample[isample] / avar
    } else {
      avar <- 2 * pi * tau[isample]
      acor.sample[isample] <- fast_mean(cos(u * pi))
      acov.sample[isample] <- avar * acor.sample[isample]
      deviance.sample <- NA
    }
  }
  
  # PSD reconstructed (with 90% CI)
  fpsd.s <- apply(fpsd.sample, 1, median)
  fpsd.s05 <- apply(fpsd.sample, 1, quantile, probs = 0.05)
  fpsd.s95 <- apply(fpsd.sample, 1, quantile, probs = 0.95)
  fpsd.mean <- apply(fpsd.sample, 1, mean)
  
  fpsd.mad <- apply(fpsd.sample, 1, mad)
  fpsd.help <- apply(fpsd.sample, 1, uniformmax)
  Cvalue <- quantile(fpsd.help, 0.9)
  confupper <- fpsd.s + Cvalue * fpsd.mad
  conflower <- fpsd.s - Cvalue * fpsd.mad
  
  # Transformed versions of these for uniform CI construction
  log.fpsd.s <- apply(log.fpsd.sample, 1, median)
  log.fpsd.mad <- apply(log.fpsd.sample, 1, mad)
  log.fpsd.help <- apply(log.fpsd.sample, 1, uniformmax)
  log.Cvalue <- quantile(log.fpsd.help, 0.9)
  log.confupper <- exp(log.fpsd.s + log.Cvalue * log.fpsd.mad)
  log.conflower <- exp(log.fpsd.s - log.Cvalue * log.fpsd.mad)
  
  # Autocovariance/autocorellation: Monte Carlo estimate (posterior mean).
  acov.mean <- mean(acov.sample)
  acov.median <- median(acov.sample)
  acov.sd = sd(acov.sample)
  acov.90ci <- quantile(acov.sample, probs = c(0.05, 0.95))  # This is also the uniform CI?
  
  acor.mean <- mean(acor.sample)
  acor.median <- median(acor.sample)
  acor.sd <- sd(acor.sample)
  acor.90ci <- quantile(acor.sample, probs = c(0.05, 0.95))
  
  # Autocovariance/autocorellation: As given by posterior median psd.
  # TODO: MC integration here, too?
  acov <- sum(fpsd.s * cos(lambda)) * 4 * pi / n
  avar <- sum(fpsd.s) * 4 * pi / n
  acor <- acov / avar
  
  f.mean <- unrollPsd(fpsd.mean, n)
  if (corrected) {
    if (toggle) {
      if (ar.order > 1) {
        ar.deviance <- apply(apply(rho, 2, function(e) {pacf2AR(e)[ar.order,]}), 1, mean)
      } else {
        ar.deviance <- mean(apply(rho, 2, function(e) {pacf2AR(e)[ar.order,]}))
      }
    } else {
      if (ar.order > 0) {
        ar.deviance <- pacf2AR(AR.fit)[ar.order,]
      } else {
        ar.deviance <- numeric(0)
      }
    }
    if (prior.q) {
      C <- sqrt(unrollPsd(psd_arma(lambda,ar.deviance,numeric(0),1),n)/f.mean)
    } else {
      C <- 1 / sqrt(f.mean)
    }
    C[1] <- C[n] <- 0
    FCFZ <- fast_ift(C * FZ)
  }
  
  return(list(fpsd.s = fpsd.s,
              fpsd.s05 = fpsd.s05,     
              fpsd.s95 = fpsd.s95,
              fpsd.mean = fpsd.mean,
              acor = acor,
              acov = acov,
              acor.mean = acor.mean,
              acor.median = acor.median,
              acor.sd = acor.sd,
              acor.90ci = acor.90ci,
              acov.mean = acov.mean,
              acov.median = acov.median,
              acov.sd = acov.sd,
              acov.90ci = acov.90ci,
              pdgrm = pdgrm,
              k = k,
              tau = tau,     
              V = V,   
              W = W,
              rho = rho,
              data = data,
              fpsd.mad = fpsd.mad,
              fpsd.help = fpsd.help,
              confupper = confupper,
              conflower = conflower,
              log.fpsd.s = log.fpsd.s,
              log.fpsd.mad = log.fpsd.mad,
              log.fpsd.help = log.fpsd.mad,
              log.confupper = log.confupper,
              log.conflower = log.conflower,
              f.alpha=f.alpha,
              lpostTrace=lpostTrace,
              theta=theta,
              missingValues_trace=missingValues_trace))
}
