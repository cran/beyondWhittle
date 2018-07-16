#' Time series model X_t=e_t, E[e_t]=0, without additional parameters
#' Not public, since encapsulated in an extra function, for easier usage
#' @keywords internal
psd_dummy_model <- function() {
  theta_dim <- 0
  excludeBoundary <- T
  get_noise <- function(data, theta, ...) {
    # mean centered version
    data - mean(data)
  }
  propose_next_theta <- function(data, f, previous_theta, ...) {
    # dummy
    numeric(0)
  }
  initialize_theta <- function(data, ...) {
    # dummy
    numeric(0)
  }
  lprior_theta <- function(theta, ...) {
    # dummy
    0
  }
  model_params <- list(theta_dim=theta_dim,
                       get_noise=get_noise,
                       propose_next_theta=propose_next_theta,
                       lprior_theta=lprior_theta,
                       initialize_theta=initialize_theta,
                       excludeBoundary=excludeBoundary)
  return(model_params)
}


#' Normal mean model, with nuisance time series
#' 
#' This class represents the mean model X_t=mu+e_t
#' with mu~N(mu.mu0,mu.sd^2) and e_t being a nuisance parameter time series
#' @details The returned object of this function is intended for usage within
#' \link[beyondWhittle]{gibbs_AR_nuisance}, \link[beyondWhittle]{gibbs_NP_nuisance}
#' and \link[beyondWhittle]{gibbs_NPC_nuisance}.
#' The method \code{propose_next_theta} is optimized to be close to the
#' marginal posterior of mu in the model.
#' The proposal scaling can be controlled with the parameter \code{prop.scaling},
#' where larger values yield a broader (smaller values yield narrower) 
#' proposal distribution.
#' @param mu.mu0,mu.sd A priori mean and standard deviation of mu
#' @param prop.scaling Scaling parameter for generating Metropolis-Hastings
#' proposals of parameter of interest theta=mu
#' @return S3 \code{nuisanceModel} object representing the model parameter theta=mu
#' of interest, containing the following fields:
#'   \item{theta_dim}{Dimension of parameter of interest (here: \code{theta_dim=1})}
#'   \item{excludeBoundary}{Logical; Should the outermost Fourier frequencies be 
#'   ignored in the frequency domain representation? (here: \code{excludeBoundary=F})}
#'   \item{get_noise}{Function taking the two arguments \code{data,theta} 
#'   to compute the nuisance/noise time series e_t from data and parameter 
#'   theta of interest. (here: e_t=data-theta)}
#'   \item{propose_next_theta}{Function taking the parameters \code{data} 
#'   (Numeric vector of input data), \code{f} (Numeric Vector of current 
#'   spectral density at the Fourier frequencies within the Gibbs sampling algorithm) 
#'   and \code{previous_theta} (Previously sampled value of mu) and 
#'   returning a new proposal value for mu}
#'   \item{initialize_theta}{Function taking the Numeric Vector \code{data} of input
#'   data as argument to generate an initial value for mu to start an MCMC
#'   algorithm (here: \code{mean(data)})}
#'   \item{lprior_theta}{Function; Log density of prior of mu (here: log density of N(mu.mu0,mu.sd^2))}
#' @keywords internal
nuisanceModel_mean <- function(mu.mu0=0, # = mu_0, a priori mean of mu
                               mu.sd=1e4, # = sigma_0, a priori standard deviation of mu
                               prop.scaling=1 # proposal scaling parameter
) {
  theta_dim <- 1
  excludeBoundary <- F
  get_noise <- function(data, theta, ...) {
    return(data - theta)
  }
  propose_next_theta <- function(data, f, previous_theta, ...) {
    n <- length(data)
    sd_prop <- prop.scaling * max(5*sqrt(pi*f[1]/n), 2/sqrt(n))
    rnorm(1, mean(data), sd_prop)
  }
  initialize_theta <- function(data, ...) {
    mean(data)
  }
  lprior_theta <- function(theta, ...) {
    dnorm(theta, mu.mu0, mu.sd, log=T)
  }
  model_params <- structure(list(theta_dim=theta_dim,
                                 get_noise=get_noise,
                                 propose_next_theta=propose_next_theta,
                                 lprior_theta=lprior_theta,
                                 initialize_theta=initialize_theta,
                                 excludeBoundary=excludeBoundary), 
                            class="nuisanceModel")
  return(model_params)
}

#' Normal linear trend model, with nuisance time series
#' 
#' This class represents the linear trend model X_t=bt + mu + e_t, t=1,...,n,
#' with intercept mu~N(mu.mu0,mu.sd^2) and slope b~N(b.mu1,b.sd^2)
#' and e_t being a nuisance parameter time series.
#' @details The parameters mu and b are assumed to be a priori independent.
#' The returned object of this function is intended for usage within
#' \link[beyondWhittle]{gibbs_AR_nuisance}, \link[beyondWhittle]{gibbs_NP_nuisance}
#' and \link[beyondWhittle]{gibbs_NPC_nuisance}.
#' The method \code{propose_next_theta} is optimized to be close to the
#' marginal joint posterior of (mu,b) in the model.
#' The proposal scaling can be controlled with the parameter \code{prop.scaling},
#' where larger values yield a broader (smaller values yield narrower) 
#' proposal distribution.
#' @param mu.mu0,mu.sd A priori mean and standard deviation of mu
#' @param b.mu1,b.sd A priori mean and standard deviation of b
#' @param prop.scaling Scaling parameter for generating Metropolis-Hastings
#' proposals of parameter of interest theta=c(mu,b)
#' @return S3 \code{nuisanceModel} object representing the model parameter theta=c(mu,b)
#' of interest, containing the following fields:
#'   \item{theta_dim}{Dimension of parameter of interest (here: \code{theta_dim=2})}
#'   \item{excludeBoundary}{Logical; Should the outermost Fourier frequencies be 
#'   ignored in the frequency domain representation? (here: \code{excludeBoundary=F})}
#'   \item{get_noise}{Function taking the two arguments \code{data,theta} 
#'   to compute the nuisance/noise time series e_t from data and parameter 
#'   theta of interest. (here: e_t=data-mu-b*(1:n))}
#'   \item{propose_next_theta}{Function taking the parameters \code{data} 
#'   (Numeric vector of input data), \code{f} (Numeric Vector of current 
#'   spectral density at the Fourier frequencies within the Gibbs sampling algorithm) 
#'   and \code{previous_theta} (Previously sampled value of c(mu,b)) and 
#'   returning a new proposal value for c(mu,b)}
#'   \item{initialize_theta}{Function taking the Numeric Vector \code{data} of input
#'   data as argument to generate an initial value for c(mu,b) to start an MCMC
#'   algorithm}
#'   \item{lprior_theta}{Function; Log density of prior of theta}
#' @importFrom MASS ginv
#' @importFrom MASS mvrnorm
#' @importFrom stats lm
#' @keywords internal
nuisanceModel_linearTrend <- function(mu.mu0=0, mu.sd=1e4, b.mu1=0, b.sd=1, 
                                           prop.scaling=1) {
  #
  # theta[1] = mu
  # theta[2] = b
  #
  theta_dim <- 2
  excludeBoundary <- F
  get_noise <- function(data, theta, ...) {
    n <- length(data)
    Yt <- 1:n
    data - theta[2]*Yt - theta[1]
  }
  propose_next_theta <- function(data, f, previous_theta, ...) {
    n <- length(data)
    N <- length(f)
    f_scaling <- c(pi, rep(2*pi, N-1))
    if (!(n%%2)) f_scaling[N] <- pi
    yf <- fast_ft(data)
    invVar <- 1/(f*f_scaling)
    Yt <- 1:n
    Xf = cbind(fast_ft(rep(1, n)), fast_ft(Yt))
    dummy <- matrix(NA, nrow = ncol(Xf), ncol = nrow(Xf))
    for (jj in 1:nrow(dummy)) {
      dummy[jj, ] <- invVar * Xf[, jj]
    }
    Sigma.inv <- dummy %*% Xf  # Xf is the FTd design matrix
    Sigma <- ginv(Sigma.inv)
    mu <- Sigma %*% dummy %*% yf  # yf is original FTd data (i.e., signal + noise)
    mvrnorm(1, mu, prop.scaling*2*Sigma)  # Propose intercept and slope from conjugate-ish distribution
  }
  initialize_theta <- function(data, ...) {
    n <- length(data)
    Yt <- 1:n
    as.numeric(lm(data ~ Yt)$coef)
  }
  lprior_theta <- function(theta, ...) {
    dnorm(theta[1], mu.mu0, mu.sd, log=T) + dnorm(theta[2], b.mu1, b.sd, log=T)
  }
  model_params <- list(theta_dim=theta_dim,
                       get_noise=get_noise,
                       propose_next_theta=propose_next_theta,
                       lprior_theta=lprior_theta,
                       initialize_theta=initialize_theta,
                       excludeBoundary=excludeBoundary)
  return(model_params)
}