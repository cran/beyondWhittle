#' Gibbs sampler for corrected parametric likelihood + Bernstein-Dirichlet mixture,
#' including possibility of using time series as mere nuisance parameter
#' Multivariate case
#' @importFrom stats dlnorm integrate rexp rlnorm
#' @keywords internal
gibbs_multivariate_nuisance <- function(data,
                                        mcmc_params,
                                        corrected,
                                        prior_params,
                                        model_params) {
  
  stopifnot(!corrected)
  
  if (ncol(data) == 1) {
    stop("Use 1D algorithm for univariate data")
  }
  d <- ncol(data)
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
  stopifnot(!is.null(mcmc_params$numerical_thresh)); stopifnot(mcmc_params$numerical_thresh>0)
  NUMERICAL_THRESH <- mcmc_params$numerical_thresh # >= 1e-12 recommended
  stopifnot(!is.null(mcmc_params$verbose))
  verbose <- mcmc_params$verbose
  
  # Adaptive MCMC parameters
  stopifnot(!is.null(mcmc_params$Nadaptive)); stopifnot(mcmc_params$Nadaptive>0)
  Nadaptive <- mcmc_params$Nadaptive
  stopifnot(!is.null(mcmc_params$adaption.batchSize)); stopifnot(mcmc_params$adaption.batchSize>0)
  adaption.batchSize <- mcmc_params$adaption.batchSize
  stopifnot(!is.null(mcmc_params$adaption.targetAcceptanceRate)); stopifnot(mcmc_params$adaption.targetAcceptanceRate>0 && mcmc_params$adaption.targetAcceptanceRate<1)
  adaption.targetAcceptanceRate <- mcmc_params$adaption.targetAcceptanceRate
  
  # AGAMMA PRIOR PAREMETERS
  stopifnot(!is.null(prior_params$eta)); stopifnot(prior_params$eta > d-1) # eta fixed
  eta <- prior_params$eta
  stopifnot(!is.null(prior_params$omega)); stopifnot(prior_params$omega > 0)
  omega_prior <- prior_params$omega
  stopifnot(!is.null(prior_params$Sigma)); stopifnot(class(prior_params$Sigma)=="matrix");
  Sigma <- prior_params$Sigma
  stopifnot(!is.null(prior_params$k.theta)); stopifnot(prior_params$k.theta > 0)
  k.theta <- prior_params$k.theta
  stopifnot(!is.null(prior_params$kmax)); stopifnot(prior_params$kmax > 0)
  kmax <- prior_params$kmax
  stopifnot(!is.null(prior_params$bernstein_l) && !is.null(prior_params$bernstein_r)); stopifnot(prior_params$bernstein_l >= 0 && prior_params$bernstein_r <= 1)
  bernstein_l <- prior_params$bernstein_l
  bernstein_r <- prior_params$bernstein_r
  stopifnot(!is.null(prior_params$coarsened))
  coarsened <- prior_params$coarsened
  stopifnot(!is.null(prior_params$L)); stopifnot(prior_params$L > 0)
  L <- prior_params$L

  # handle missing values: Initialize
  NA_tim <- which(apply(data, 1, function(z) any(is.na(z))))
  NA_pos <- which(is.na(data))
  num_NA <- length(NA_pos)
  has_NA <- (num_NA > 0)
  missingValues_trace <- matrix(NA, nrow=Ntotal, ncol=num_NA)
  missingValues_trace[1,] <- rep(0,num_NA)
  data[NA_pos] <- missingValues_trace[1,]
  if (has_NA) {
    cat("NOTE: Missing values at the following time points:",
        NA_tim, "treating them as random.\n")
  }
    
  omega_freq <- omegaFreq(n)
  N <- length(omega_freq)
  lambda <- pi * omega_freq
  if (is.finite(kmax)) {
    db.list <- dbList(n, kmax, bernstein_l, bernstein_r, 
                      coarsened, verbose=T) 
  } else {
    stop("kmax=Inf not supported yet")
  }
  
  # initialize storage
  lpostTrace <- rep(NA, Ntotal)
  U__phi <- array(NA, dim=c(d*d-1, L, Ntotal)) # Mittelbach et al
  r <- matrix(NA, nrow=L, ncol=Ntotal)
  Z <- matrix(NA, nrow=L, ncol=Ntotal)
  k <- rep(NA, Ntotal)
  num_except <- 0
  
  # starting values
  k_0 <- round(kmax / 2)
  Z_0 <- seq(1/(2*k_0), 1-1/(2*k_0), length.out=L)
  r_0 <- rep(1/L, L)
  U_phi_0 <- unit_trace_runif(L, d)$phi
  U__IND_ALL <- 1:(d*d-1)
  U__IND_PIHALF <- (1:(d-1))^2
  U__IND_PI <- setdiff(U__IND_ALL, U__IND_PIHALF)
  U__SCALING <- U__IND_ALL
  U__SCALING[U__IND_PI] <- pi
  U__SCALING[U__IND_PIHALF] <- pi / 2
  
  # proposal variances
  eps_r <- eps_Z <- eps_U <- seq(1, L) / (seq(1, L) + 2 * sqrt(n))
  lsd_r <- log(eps_r) / 2
  lsd_Z <- log(eps_Z) / 2
  lsd_U <- log(eps_U) / 2
  
  FZ <- mdft(data)
  
  ##
  ##  MH-within Gibbs sampler: C++ in batches
  ##
  N_BATCH <- min(adaption.batchSize, 50) # cpp batch size
  while (print_interval %% N_BATCH) print_interval <- print_interval + 1
  tim0 <- proc.time()
  i_s <- 1:(Ntotal/N_BATCH)
  for (i in i_s) {
    
    # last batch may be larger than N_BATCH
    if (i==tail(i_s,1) && Ntotal%%N_BATCH>0) {
      N_MCMC <- N_BATCH + Ntotal%%N_BATCH
      CURRENT_BATCH <- ((i-1)*N_BATCH+1):(i*N_BATCH + Ntotal%%N_BATCH) 
    } else {
      N_MCMC <- N_BATCH
      CURRENT_BATCH <- ((i-1)*N_BATCH+1):(i*N_BATCH)
    }
    
    # adaption of proposal variances
    i_mcmc <- (i-1)*N_BATCH+1
    if ((i_mcmc < Nadaptive) && (i_mcmc > 1) && (i_mcmc %% adaption.batchSize == 1)) {
      batch <- (i_mcmc-adaption.batchSize):(i_mcmc-1)
      adaption.delta <- min(0.05, 1/(i_mcmc^(1/2))) # c.f. Rosenthal
      ### r
      batch.r <- r[, batch]
      batch.r.acceptanceRate <- apply(batch.r, 1, acceptanceRate) 
      lsd_r <- lsd_r + ((batch.r.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      eps_r <- exp(2*lsd_r)
      ### U (first component suffices)
      batch.U <- U__phi[1,,batch]
      batch.U.acceptanceRate <- apply(batch.U, 1, acceptanceRate)
      lsd_U <- lsd_U + ((batch.U.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      eps_U <- exp(2*lsd_U)
    }
    
    # call C++
    mcmc_i <- gibbs_multivariate_nuisance_cpp(data=data,
                                              # note the -1 because cpp is 0-indexed
                                              NA_pos=NA_pos-1, 
                                              FZ=FZ,
                                              eps_r=eps_r,
                                              eps_Z=eps_Z,
                                              eps_U=eps_U,
                                              k_0=k_0,
                                              r_0=r_0,
                                              Z_0=Z_0,
                                              U_phi_0=U_phi_0,
                                              phi_def=U__SCALING,
                                              eta=eta,
                                              omega=omega_prior,
                                              Sigma=Sigma,
                                              Ntotal=N_MCMC,
                                              print_interval=print_interval,
                                              numerical_thresh=NUMERICAL_THRESH,
                                              verbose=F,
                                              L=L,
                                              k_theta=k.theta,
                                              dbList=db.list)
    
    # collect results
    k[CURRENT_BATCH] <- mcmc_i$k
    Z[,CURRENT_BATCH] <- mcmc_i$Z
    r[,CURRENT_BATCH] <- mcmc_i$r
    U__phi[,,CURRENT_BATCH] <- mcmc_i$U_phi
    lpostTrace[CURRENT_BATCH] <- mcmc_i$lpostTrace
    num_except <- num_except + mcmc_i$num_except
    missingValues_trace[CURRENT_BATCH,] = mcmc_i$missingValuesTrace
    
    # starting values for next iteration
    k_0 <- mcmc_i$k[N_MCMC]
    Z_0 <- mcmc_i$Z[,N_MCMC]
    r_0 <- mcmc_i$r[,N_MCMC]
    U_phi_0 <- mcmc_i$U_phi[,,N_MCMC]
    data[NA_pos] <- mcmc_i$missingValuesTrace[N_MCMC,]
    if (verbose || !(((i+1)*N_MCMC)%%print_interval)) {
      tim <- proc.time()
      print_mcmc_state((i+1)*N_MCMC, Ntotal, tim, tim0)
    }
  }
  if (is.finite(kmax) && min(abs(k-kmax), na.rm=T) < 5) {
    print_warn("Reached kmax during sampling. Consider using a higher kmax value")
  }
  if (num_except > 0) {
    print_warn(paste0(num_except, " MH proposals discarded due to numerical issues"))
  }
  
  ##
  ## Post processing
  ##
  cat("Postprocessing...")
  keep <- seq(burnin+1, Ntotal, by=thin)
  k <- k[keep]
  r <- r[,keep]
  Z <- Z[,keep]
  U__phi <- U__phi[,,keep]
  U <- array(NA, dim=c(d,d,L,length(keep)))
  missingValues_trace <- missingValues_trace[keep,,drop=F]
  colnames(missingValues_trace) <- missingValues_str_help(NA_pos, n)
  data[NA_pos] <- NA
  
  fpsd.sample <- array(NA, dim=c(d, d, N, length(keep)))
  for (isample in 1:length(keep)) {
    U[,,,isample] <- get_U_cpp(U__phi[,,isample])
    f_sample <- get_f_matrix(U__phi[,,isample], r[,isample], Z[,isample], k[isample], db.list)
    # Internally, we store the PSDs in a non-redundant real-valued representation,
    # to facilitate the postprocessing (computation of quantiles)
    fpsd.sample[,,,isample] <- realValuedPsd(f_sample)
  }
  
  fpsd.s <- apply(fpsd.sample, c(1,2,3), median)
  fpsd.mean <- apply(fpsd.sample, c(1,2,3), mean)
  fpsd.s05 <- apply(fpsd.sample, c(1,2,3), quantile, 0.05)
  fpsd.s95 <- apply(fpsd.sample, c(1,2,3), quantile, 0.95)
  
  alpha_uci <- 0.1 # same as in 1D
  uci_tmp <- uci_matrix(fpsd.sample, alpha=alpha_uci)
  fpsd.uci05 <- uci_tmp$fpsd.uci05
  fpsd.uci95 <- uci_tmp$fpsd.uci95
  uuci_tmp <- uci_matrix(fpsd.sample, alpha=alpha_uci, uniform_among_components=T)
  fpsd.uuci05 <- uuci_tmp$fpsd.uci05
  fpsd.uuci95 <- uuci_tmp$fpsd.uci95
  
  print("done!")
  
  ##
  ## Return stuff
  ##
  return(list(data=data,
              k=k,
              r=r,
              Z=Z,
              U=U,
              lpostTrace=lpostTrace,
              fpsd.s=fpsd.s,
              fpsd.mean=fpsd.mean,
              fpsd.s05=fpsd.s05,
              fpsd.s95=fpsd.s95,
              fpsd.uuci05=fpsd.uuci05,
              fpsd.uuci95=fpsd.uuci95,
              missingValues_trace=missingValues_trace))
}
