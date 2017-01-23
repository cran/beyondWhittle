#' Gibbs sampler for corrected parametric likelihood + Bernstein-Dirichlet mixture
#' @keywords internal
gibbs_toggle_ar <- function(data,
                         Ntotal,
                         burnin,
                         thin = 1,
                         M = 1,
                         g0.alpha = 1,
                         g0.beta = 1,
                         k.theta = 0.01,
                         tau.alpha = 0.001,
                         tau.beta = 0.001,
                         kmax = 100,
                         L = max(20, length(data) ^ (1 / 3)),
                         corrected = TRUE,
                         prior.q=FALSE,  # use iff corrected: prior.q for (Bernstein) prior on q, !prior.q for prior on f
                         alpha.toggle=F, # use iff prior.q: prior on f/f_param^alpha (alpha==1 if alpha.toggle=F, alpha==0 if prior.q=F)
                         toggle=FALSE,
                         Nadaptive = 0, # use iff toggle: number of adaptive iterations for rho proposal variance
                         adaption.batchSize = 50, # use iff toggle: size of batches to compute acceptance for adaption
                         adaption.targetAcceptanceRate = 0.44, # use iff toggle
                         mu.prop = NULL,  # use iff toggle
                         var.prop.rho = NULL, #  " "
                         AR.fit=NULL,     # use iff !toggle
                         MA.fit=NULL,     #  " "
                         dist = "normal",
                         kurt.toggle = TRUE,     # use iff dist=="student"
                         kurt.lambda = 0.1,      # use iff kurt.toggle: exponential prior parameter
                         kurt.thresh = 6/96,     #   " "              : threshold towards normal likelihood
                         kurt.star.alpha = NULL, #   " "              : proposal parameter
                         kurt.star.beta = NULL,  #   " "              : proposal parameter
                         kurt.fit=NULL,          # use iff !kurt.toggle
                         beta.toggle = F,        # use iff dist=="generalized"
                         beta.alpha = 2,                      # prior parameter for beta (generalized gaussian)
                         beta.beta = beta.alpha / 2,          #   " "
                         beta.fit=NULL,
                         rho.alpha=rep(1, length(mu.prop)), # use iff prior && toggle
                         rho.beta=rep(1, length(mu.prop)))  # use iff prior && toggle

{
  if (!corrected) {
    stopifnot(dist=="normal") # TODO Remove corrected parameter
  }

  if (dist=="generalized" && !beta.toggle) {
    stopifnot(!is.null(beta.fit))
  }

  if (dist != "generalized") {
    beta <- NULL # ameier HACK
    beta.alpha <- NULL # ameier HACK
    beta.beta <- NULL # ameier HACK
  }


  n <- length(data)

  # Tolerance for mean centering
  tol <- 1e-4

  # Mean center - necessary for our implementation of FFT
  if (abs(mean(data)) > tol) {
    data <- data - mean(data)
    warning("Data has been mean centered for your convenience")
  }

  # Basic error messages - TO COMPLETE
  if (burnin > Ntotal) {
    stop("Burn-in parameter is larger than number of iterations")
  }

  if (corrected) {
     if (toggle) {
      # We expect not null here.
      # Note: To parse empty parameter sets, use numeric(0) instead of NULL.
      stopifnot(!is.null(mu.prop) && !is.null(var.prop.rho))
      ar.order <- length(mu.prop)
      stopifnot(ar.order > 0)
      stopifnot(ar.order==length(mu.prop))
      stopifnot(ar.order==length(var.prop.rho))
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

  if (dist=="student") {
    if (kurt.toggle) {
      stopifnot(!is.null(kurt.star.alpha) && !is.null(kurt.star.beta))
    } else {
      stopifnot(!is.null(kurt.fit))
      stopifnot(kurt.fit > 0) # need positive excess kurtosis for student distribution.
    }
  }

  # FFT data to frequency domain
  FZ <- fast_ft(data)
  FZ[1] <- FZ[n] <- 0

  # Periodogram
  pdgrm <- abs(FZ) ^ 2  # Note the length is n here. Note the missing 1 / sqrt(2*pi).

  # Frequencies on unit interval
  twon <- 2 / n
  omega <- twon * (1:(n / 2 + 1) - 1)  #2 * (1:(n / 2 + 1) - 1) / n

  # Angular frequencies on [0, pi]
  lambda <- pi * omega  #2 * pi * (1:(n / 2 + 1) - 1) / n

  #####
  # Store list for beta mixtures - WORK IN PROGRESS
  # TO DO: A more elegant solution to the storage problem
  # Can come across huge memory/storage problems with large kmax and n
  db.names <- paste("db", 1:kmax, sep = "")
  db.list <- vector("list", kmax)
  names(db.list) <- db.names

  for (kk in 1:kmax) {
    db.list[[kk]] <- matrix(dbeta(omega,
                                  rep(1:kk, each = length(omega)),
                                  rep(kk:1, each = length(omega))),
                            ncol = length(omega),
                            byrow = TRUE)
  }
  #####

  # Open objects for storage
  lpostTrace <- rep(NA, Ntotal)
  tau <- rep(NA, Ntotal)
  V <- matrix(NA, nrow = L, ncol = Ntotal)
  W <- matrix(NA, nrow = L + 1, ncol = Ntotal)  # Called Z in Choudhuri
  k <- rep(NA, Ntotal)

  # Starting values
  tau[1] <- var(data) / (2 * pi)
  V[, 1] <- rbeta(L, 1, M)
  W[, 1] <- rbeta(L + 1, g0.alpha, g0.beta)  # g0.alpha = g0.beta = 1 gives U[0,1]
  k[1] <- round(kmax / 2)
  if (corrected) {
    if (toggle) {
      rho <- matrix(NA, nrow=ar.order, ncol=Ntotal) # trace of PACFs
      rho[,1] <- mu.prop
      AR <- pacf2AR(rho[,1])[ar.order,] # current AR(p) parameters
      MA <- numeric(0)
      lsd.prop.rho <- log(var.prop.rho) / 2
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

  eps <- seq(1, L + 1) / (seq(1, L + 1) + 2 * sqrt(n))  # Metropolis proposal parameters for V and W.
  # First one only used for W0.

  if (dist == "student") {
    if (kurt.toggle) {
      ex.kurt <- rep(NA, Ntotal)
      ex.kurt[1] <- kurt.star.alpha * kurt.star.beta # mean as starting value
      if (ex.kurt[1] >= kurt.thresh) {
        nll_fun <- nll_t
      } else {
        nll_fun <- nll_norm
      }
    } else {
      ex.kurt <- rep(kurt.fit, Ntotal)
      nll_fun <- nll_t
    }
  }
  if (dist == "normal") {
    sigma2 <- rep(1, Ntotal)
    ex.kurt <- rep(0, Ntotal) # Note: Important to set this consistently -- to keep the prior consistent!
    df <- rep(Inf, Ntotal)
    nll_fun <- nll_norm
  }
  ## ameier
  if (dist=="generalized") {
    df <- rep(NA, Ntotal)
    sigma2 <- rep(1, Ntotal)
    nll_fun <- nll_generalizedGaussian
    if (beta.toggle) {
      ex.kurt <- rep(NA, Ntotal)
      beta <- rep(NA, Ntotal)
      beta[1] <- 2 # always start with normal dist
      ex.kurt[1] <- generalizedGaussian.kurtosis(beta[1])
      beta.thresh.lower <- 0.1 # make sure proposals remain numericcaly tractable
      beta.thresh.upper <- 500
      var.prop.beta <- 0.25
      lsd.prop.beta <- log(var.prop.beta) / 2
    } else {
      ex.kurt <- rep(generalizedGaussian.kurtosis(beta.fit), Ntotal)
      beta <- rep(beta.fit, Ntotal)
    }
  }
  if (dist != "normal") {
    var.prop.tau <- 0.01 # initial guess
    lsd.prop.tau <- log(var.prop.tau) / 2
  }
  ###

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
  for (i in 1:(Ntotal-1)) {

    if (i%%200 == 0) {
      cat("iteration ", i, "\n")
    }

    # Adaption step: Adjust propsal variance
    if ((i < Nadaptive) && (i > 1) && (i %% adaption.batchSize == 1)) {
      batch <- (i-adaption.batchSize):(i-1)
      adaption.delta <- min(0.1, 1/(i^(1/3))) # c.f. Rosenthal
      ### rho
      if (corrected && toggle) {
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
      ### tau
      if (dist != "normal") {
        batch.tau <- tau[batch]
        batch.tau.acceptanceRate <- acceptanceRate(batch.tau)
        lsd.prop.tau <- lsd.prop.tau + ((batch.tau.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
        var.prop.tau <- exp(2*lsd.prop.tau)
      }
      ### beta
      if (dist == "generalized" && beta.toggle) {
        batch.beta <- beta[batch]
        batch.beta.acceptanceRate <- acceptanceRate(batch.beta)
        lsd.prop.beta <- lsd.prop.beta + ((batch.beta.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
        var.prop.beta <- exp(2*lsd.prop.beta)
      }
    }

    #####
    # Step 1: Sample from full conditional of k: WORK IN PROGRESS
    # Metropolis proposal for k
    #k.star <- round(rt(1, 1, k[i]))  # Cauchy distribution discretized
    #while (k.star < 1 || k.star > kmax) {  # A bit hacky
    #  k.star <- round(rt(1, 1, k[i]))
    #}
    k.star <- round(rt(1, 1)) + k[i]  # Cauchy distribution discretized
    while (k.star < 1 || k.star > kmax) {  # A bit hacky
      k.star <- round(rt(1, 1)) + k[i]
    }
    #####

    # log posterior for proposal
    f.k.star <- lpost(omega,
                      FZ,
                      AR,
                      MA,
                      V[, i],
                      W[, i],
                      k.star,  # Sampled value of k
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
                      db.list,
                      dist,
                      nll_fun,
                      ex.kurt[i],
                      kurt.lambda,
                      beta[i],
                      beta.alpha,
                      beta.beta,
                      f.alpha[i],
                      rho[,i],
                      rho.alpha,
                      rho.beta)

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
                 db.list,
                 dist,
                 nll_fun,
                 ex.kurt[i],
                 kurt.lambda,
                 beta[i],
                 beta.alpha,
                 beta.beta,
                 f.alpha[i],
                 rho[,i],
                 rho.alpha,
                 rho.beta)

    if (i==1) {
      lpostTrace[1] <- f.k
    }

    #####
    # Accept/reject
    alpha1 <- min(0, f.k.star - f.k)  # log acceptance ratio
    if (log(runif(1, 0, 1)) < alpha1) {
      k[i + 1] <- k.star  # Accept k.star
      f.store <- f.k.star
    } else {
      k[i + 1] <- k[i]  # Reject and use previous
      f.store <- f.k
    }
    #####

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
                        db.list,
                        dist,
                        nll_fun,
                        ex.kurt[i],
                        kurt.lambda,
                        beta[i],
                        beta.alpha,
                        beta.beta,
                        f.alpha[i],
                        rho[,i],
                        rho.alpha,
                        rho.beta)

      # log posterior of previous iteration
      f.V <- f.store

      # Accept/reject
      alpha2 <- min(0, f.V.star - f.V)  # log acceptance ratio
      if (log(runif(1, 0, 1)) < alpha2) {
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
                        db.list,
                        dist,
                        nll_fun,
                        ex.kurt[i],
                        kurt.lambda,
                        beta[i],
                        beta.alpha,
                        beta.beta,
                        f.alpha[i],
                        rho[,i],
                        rho.alpha,
                        rho.beta)

      # log posterior for previous iteration
      f.W <- f.store

      # Accept/reject
      alpha3 <- min(0, f.W.star - f.W)  # log acceptance ratio
      if(log(runif(1, 0, 1)) < alpha3) {
        W[l, i + 1] <- W.star[l]  # Accept W.star
        f.store <- f.W.star
      } else {
        W[l, i + 1] <- W[l, i]  # Reject and use previous
      }

    }  # END: Step 3.


    # Step 4: Directly sample tau from conjugate Inverse-Gamma density
    q.psd <- qpsd(omega, V[, i + 1], W[, i + 1], k[i + 1], db.list)$psd
    q <- rep(NA, n)
    q[1] <- q.psd[1]
    q[n] <- q.psd[length(q.psd)]
    q[2 * 1:(n / 2 - 1)] <- q[2 * 1:(n / 2 - 1) + 1] <- q.psd[1:(n / 2 - 1) + 1]

    if (corrected == FALSE) {  # For Whittle likelihood

      # Note the (n - 2) here - we remove the first and last terms
      tau[i + 1] <- 1 / rgamma(1, tau.alpha + (n - 2) / 2,
                               tau.beta + sum(pdgrm[2:(n - 1)] / q[2:(n - 1)]) / (2 * pi) / 2)

    }  # CHECK: Should this be 2pi here or in pdgrm?

    if (corrected == TRUE) {  # For Corrected likelihood

      if (dist == "normal") {  # Conjugacy for Normal likelihood!

        # Correction matrix C_n^-0.5: Notice this is for q, not f.
        # Let C = tau * B.  B is not normalised.  Needed for full conditional.
        #B <- Cn(n, q, AR.fit, MA.fit, sigma2.fit)  # Remove first and last for B
        # The q here is actually q / f_param
        if (prior.q) {
          B <- Cn(n, q / (unrollPsd(psd_arma(lambda,AR,MA,1),n))^(1-f.alpha[i]))
        } else {
          B <- Cn(n, q / unrollPsd(psd_arma(lambda,AR,MA,1),n))
        }

        B[1] <- B[n] <- 0

        # Input for ARMA parametric conditional likelihood - Inverse FFT
        FBFZ <- fast_ift(B * FZ)

        # Calculate ARMA parametric conditional likelihood
        #p.arma <- arma_conditional(FBFZ, AR.fit, MA.fit, sigma2.fit)
        #p.arma <- arma_conditional(FBFZ, AR, MA, nll_norm_unnormalized, NULL)
        p.arma <- 1 / 2 * sum(genEpsARMAC(FBFZ, AR, MA)^2)

        # Note the (n - 2) here - we remove the first and last terms
        # Subtract p for ARMA since only p + 1 to n terms in conditional likelihood
        #tau[i + 1] <- 1 / rgamma(1, tau.alpha + n / 2, tau.beta + p.arma)  # Note the plus.
        tau[i + 1] <- 1 / rgamma(1, tau.alpha + (n - 2) / 2, tau.beta + p.arma)

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
                         db.list,
                         dist,
                         nll_fun,
                         ex.kurt[i],
                         kurt.lambda,
                         beta[i],
                         beta.alpha,
                         beta.beta,
                         f.alpha[i],
                         rho[,i],
                         rho.alpha,
                         rho.beta)
      }


      #####
      # NOTE: tau is NOT CONJUGATE anymore.  Use MH step!
      #####

      if (dist != "normal") {

        # Proposal for tau
        sd.prop.tau <- sqrt(var.prop.tau)
        tau.star <- tau[i] + rnorm(1, mean = 0, sd = sd.prop.tau)

        if (tau.star < 0) {
          tau[i + 1] <- tau[i]
        } else {
          # log posterior for proposal
          f.tau.star <- lpost(omega,
                              FZ,
                              AR,
                              MA,
                              V[, i + 1],
                              W[, i + 1],
                              k[i + 1],
                              tau.star,
                              M,
                              g0.alpha,
                              g0.beta,
                              k.theta,
                              tau.alpha,
                              tau.beta,
                              corrected,
                              prior.q,
                              pdgrm,
                              db.list,
                              dist,
                              nll_fun,
                              ex.kurt[i],
                              kurt.lambda,
                              beta[i],
                              beta.alpha,
                              beta.beta,
                              f.alpha[i],
                              rho[,i],
                              rho.alpha,
                              rho.beta)


          # log posterior for previous iteration
          f.tau <- f.store

          # Accept/reject
          # ameier: MH-correction for random walk proposal
          alpha4 <- min(0,
                        f.tau.star -
                          f.tau +
                          dnorm(tau[i], mean=tau.star, sd=sd.prop.tau, log=T) -
                          dnorm(tau.star, mean=tau[i], sd=sd.prop.tau, log=T))  # log acceptance ratio
          if(log(runif(1, 0, 1)) < alpha4) {
            tau[i + 1] <- tau.star  # Accept
            f.store <- f.tau.star
          } else {
            tau[i + 1] <- tau[i]
          }
        }



      }
      if (dist == "student") {
        # Step 5: Sample excess-kurtosis
        if (kurt.toggle) {
          #################################
          # UNDER MAINTENANCE: excess-kurtosis proposal
          #################################
          ex.kurt.star <-  rgamma(1, shape=kurt.star.alpha, scale=kurt.star.beta)
          while (ex.kurt.star < 0 || ex.kurt.star > 20) {
            ex.kurt.star <- rgamma(1, shape=kurt.star.alpha, scale=kurt.star.beta)
          }
          # Thresholding: Use normal likelihood for small ex.kurt values.
          if (ex.kurt.star >= kurt.thresh) {
            nll_fun.star <- nll_t
          } else {
            nll_fun.star <- nll_norm
          }

          # log posterior for proposal
          f.kurt.star <- lpost(omega,
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
                               db.list,
                               dist,
                               nll_fun.star,
                               ex.kurt.star,
                               kurt.lambda,
                               beta[i],
                               beta.alpha,
                               beta.beta,
                               f.alpha[i],
                               rho[,i],
                               rho.alpha,
                               rho.beta)

          # log posterior for previous iteration
          f.kurt <- f.store

          # Accept/reject
          alpha5 <- min(0, f.kurt.star - f.kurt +
                          dgamma(ex.kurt[i],shape=kurt.star.alpha, scale=kurt.star.beta,log=TRUE) -
                          dgamma(ex.kurt.star,shape=kurt.star.alpha, scale=kurt.star.beta,log=TRUE))  # log acceptance ratio
          if(log(runif(1, 0, 1)) < alpha5) {
            ex.kurt[i + 1] <- ex.kurt.star
            nll_fun <- nll_fun.star
            f.store <- f.kurt.star
          } else {
            ex.kurt[i + 1] <- ex.kurt[i]
          }
        }  # END: toggle excess kurtosis
      }  # END: Student-t stuff

      if (dist == "generalized" && beta.toggle) { # toggle beta
        sd.prop.beta <- sqrt(var.prop.beta)
        beta.star <- exp(log(beta[i]) + rnorm(1, 0, sd.prop.beta)) # TODO This needs tuning!
        #beta.star <- beta[i] + rnorm(1, 0, sd.prop.beta) # TODO This needs tuning!
        if (beta.star <= beta.thresh.lower || beta.star >= beta.thresh.upper) {
          beta[i+1] <- beta[i]  # Reject and use previous
        } else {
          # log posterior for proposal
          f.beta.star <- lpost(omega,
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
                               db.list,
                               dist,
                               nll_fun,
                               ex.kurt[i],
                               kurt.lambda,
                               beta.star,
                               beta.alpha,
                               beta.beta,
                               f.alpha[i],
                               rho[,i],
                               rho.alpha,
                               rho.beta)

          # log posterior for previous iteration
          f.beta <- f.store


          # Accept/reject
          alpha99 <- min(0, f.beta.star -
                           f.beta +
                           dnorm(log(beta[i]),log(beta.star),sd.prop.beta,log=TRUE) -
                           dnorm(log(beta.star),log(beta[i]),sd.prop.beta,log=TRUE))
          # Note that rho proposals with abs(rho) >= 1 are rejected
          if(log(runif(1, 0, 1)) < alpha99) {
            beta[i+1] <- beta.star
            f.store <- f.beta.star
          } else {
            beta[i+1] <- beta[i]  # Reject and use previous
          }
        }
        ex.kurt[i+1] <- generalizedGaussian.kurtosis(beta[i+1])
      }

      #####
      # SAMPLE PACF, rho, here:
      #####
      if (toggle) {

        for (pp in 1:ar.order) {
          rho.star <- rho.old <- rho[, i]
          if (pp > 1) {
            for (ip in 1:(pp-1)) {
              rho.star[ip] <- rho.old[ip] <- rho[ip, i+1]
            }
          }

          rho.star[pp] <- rho.old[pp] + rnorm(1, 0, sqrt(var.prop.rho[pp]))

          if (abs(rho.star[pp]) >= 1) {
            rho.star[pp] <- rho.old[pp] # TODO TUNE
          }

          AR.star <- pacf2AR(rho.star)[ar.order,]
          AR.old <- pacf2AR(rho.old)[ar.order,]

          # log posterior for proposal
          f.rho.star <- lpost(omega,
                              FZ,
                              AR.star, # AR.star here
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
                              db.list,
                              dist,
                              nll_fun,
                              ex.kurt[i + 1],
                              kurt.lambda,
                              beta[i+1],
                              beta.alpha,
                              beta.beta,
                              f.alpha[i],
                              rho.star,
                              rho.alpha,
                              rho.beta)

          f.rho <- f.store

          # Accept/reject
          alpha6 <- min(0, f.rho.star - f.rho +
                          dnorm(rho.old[pp],rho.star[pp],sqrt(var.prop.rho[pp]),log=TRUE) -
                          dnorm(rho.star[pp],rho.old[pp],sqrt(var.prop.rho[pp]),log=TRUE))
          # Note that rho proposals with abs(rho) >= 1 are rejected
          if(abs(rho.star[pp]) < 1 && log(runif(1, 0, 1)) < alpha6) {
            rho[pp,i + 1] <- rho.star[pp]
            AR <- AR.star
            f.store <- f.rho.star
          } else {
            rho[pp,i + 1] <- rho[pp,i]  # Reject and use previous
          }
        }
        #####
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
                               db.list,
                               dist,
                               nll_fun,
                               ex.kurt[i + 1],
                               kurt.lambda,
                               beta[i],
                               beta.alpha,
                               beta.beta,
                               f.alpha.star,
                               rho[,i+1],
                               rho.alpha,
                               rho.beta)
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
    lpostTrace[i+1] <- f.store
  }  # END: MCMC loop

  # ADD WARNING IF A LOT OF SAMPLES FROM KMAX - MAY NEED TO BE LARGER

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
  beta <- beta[keep]
  ex.kurt <- ex.kurt[keep]

  if (dist == "student") {
    df <- 6 / ex.kurt + 4
    sigma2 <- df / (df - 2)
    sigma2[is.na(sigma2)] <- 1 # correction for kurtosis = 0
  }
  if (dist == "normal") {
    df <- df[keep]
    sigma2 <- sigma2[keep]
  }
  if (dist == "generalized") {
    df <- df[keep]
    sigma2 <- sigma2[keep]
  }
  if (prior.q) {
    f.alpha <- f.alpha[keep]
  }

  fpsd.sample <- log.fpsd.sample <- log.fpsd.sample2 <- matrix(NA, nrow = length(omega), ncol = length(keep))
  acov.sample <- rep(NA, length(keep))
  acor.sample <- rep(NA, length(keep))
  deviance.sample <- rep(NA, length(keep)) # ameier

  for (isample in 1:length(keep)) {
    if (toggle) {
      AR <- pacf2AR(rho[,isample])[ar.order,] # TODO Store traces of AR(p) parameters, too?
      MA <- numeric(0)
    } else {
      AR <- AR.fit
      MA <- MA.fit
    }

    qstore <- qpsd(omega, V[, isample], W[, isample], k[isample], db.list)

    if (prior.q) {
      f_param <- psd_arma(lambda, AR, MA, 1)^f.alpha[isample]
      fpsd.sample[, isample] <- tau[isample] * qstore$psd * f_param  # Multiplied by f_param
    }
    else {
      fpsd.sample[, isample] <- tau[isample] * qstore$psd
    }

    # Create transformed version
    log.fpsd.sample[, isample] <- logfuller(fpsd.sample[, isample])
    log.fpsd.sample2[, isample] <- logrosenthal(fpsd.sample[, isample])

    #     # Calculate autocovariance/autocorrelation *numerically*
    #     acov.sample[isample] <- sum(fpsd.sample[, isample] * cos(lambda)) * 4 * pi / n  # Use original fpsd.sample, not transformed
    #     avar <- sum(fpsd.sample[, isample]) * 4 * pi / n
    #     acor.sample[isample] <- acov.sample[isample] / avar

    # Calculate autocovariance/autocorrelation *by monte carlo integration*
    weight <- qstore$weight  # Weights
    jstar <- sample(1:k[isample], size = 1000, prob = weight, replace = TRUE)  # Why 1000?
    u <- rbeta(1000, jstar, k[isample] - jstar + 1)
    if (prior.q) {
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

  # Transformed versions of these for uniform CI construction [ALTERNATIVE VERSION WITH ROSENTHAL LOGARITHM]
  log.fpsd.s2 <- apply(log.fpsd.sample2, 1, median)
  log.fpsd.mad2 <- apply(log.fpsd.sample2, 1, mad)
  log.fpsd.help2 <- apply(log.fpsd.sample2, 1, uniformmax)
  log.Cvalue2 <- quantile(log.fpsd.help2, 0.9)
  log.confupper2 <- logrosenthal(log.fpsd.s2 + log.Cvalue2 * log.fpsd.mad2, inverse=T)
  log.conflower2 <- logrosenthal(log.fpsd.s2 - log.Cvalue2 * log.fpsd.mad2, inverse=T)

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
    p.arma <- NA #arma_conditional(FCFZ, ar.deviance, numeric(0), nll_fun, ex.kurt)
    deviance.mean <- NA #-2 * (sum(log(C[-c(1, n)])) - p.arma)
  } else {
    deviance.mean <- NA # sum(log(f.mean[2:(n - 1)] * 2 * pi) + pdgrm[2:(n - 1)] / (f.mean[2:(n - 1)] * 2 * pi))
  }
  DIC <- NA # 2 * mean(deviance.sample) - deviance.mean # ameier

  # Some textual description of algorithm type.
  if (corrected) {
    algorithm <- paste0("C_p=", ar.order)
    if (prior.q) {
      algorithm <- paste0(algorithm, "_q-prior")
    } else {
      algorithm <- paste0(algorithm, "_f-prior")
    }
    if (toggle) {
      algorithm <- paste0(algorithm, "_toggle")
    } else {
      AR.string <- paste0("(", paste0(AR.fit, collapse=","), ")")
      MA.string <- paste0("(", paste0(MA.fit, collapse=","), ")")
      algorithm <- paste0(algorithm, "_AR=", AR.string, "_MA=", MA.string)
    }
    if (dist=="normal") {
      algorithm <- paste0(algorithm, "_normal")
    }
    if (dist=="student") {
      algorithm <- paste0(algorithm, "_student")
      if (!kurt.toggle) {
        algorithm <- paste0(algorithm, "_ex.kurt=", kurt.fit)
      }
    }
    if (dist=="generalized") {
      algorithm <- paste0(algorithm, "_gen")
      if (!kurt.toggle) {
        algorithm <- paste0(algorithm, "_beta=", beta.fit)
      }
    }
  } else {
    algorithm <- "Whittle"
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
              DIC = DIC,
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
              log.fpsd.s2 = log.fpsd.s2,
              log.fpsd.mad2 = log.fpsd.mad2,
              log.fpsd.help2 = log.fpsd.mad2,
              log.confupper2 = log.confupper2,
              log.conflower2 = log.conflower2,
              ex.kurt = ex.kurt,
              df = df,
              sigma2 = sigma2,
              beta=beta,
              algorithm=algorithm,
              f.alpha=f.alpha,
              lpostTrace=lpostTrace))
}  # Close function
