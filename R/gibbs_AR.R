#' Gibbs sampler for an autoregressive model with PACF parametrization.
#'
#' Obtain samples of the posterior of a Bayesian autoregressive model of fixed order.
#' @details Partial Autocorrelation Structure (PACF, uniform prior) and the residual variance sigma2 (inverse gamma prior) is used as model parametrization.
#' The DIC is computed with two times the posterior variance of the deviance as effective number of parameters, see (7.10) in the referenced book by Gelman et al.
#' Further details can be found in the simulation study section in the referenced paper by C. Kirch et al.
#' @param data numeric vector
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param ar.order order of the autoregressive model (integer >= 0)
#' @param sigma2.alpha,sigma2.beta prior parameters for sigma2 (inverse gamma)
#' @return list containing the following fields:
#'
#'    \item{rho}{matrix containing traces of the PACF parameters (if p>0)}
#'    \item{sigma2}{trace of sigma2}
#'    \item{DIC}{a list containing the numeric value \code{DIC} of the Deviance Information Criterion (DIC) and the effective number of parameters \code{ENP}}
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{uniform credibility interval}
#' @references C. Kirch et al. (2017)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' <arXiv:1701.04846>
#' @references A. Gelman et al. (2013)
#' \emph{Bayesian Data Analysis, Third Edition}
#' @examples 
#' \dontrun{
#' 
#' ##
#' ## Example 1: Fit an AR(p) model to sunspot data:
#' ##
#' 
#' # Use this variable to set the autoregressive model order
#' ar.order <- 2
#' 
#' data <- sqrt(as.numeric(sunspot.year))
#' data <- data - mean(data)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_AR(data=data, Ntotal=20000, burnin=8000, ar.order=ar.order)
#' 
#' # Plot spectral estimates on log-scale (excluding the zero frequency).
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts(log(pdgrm[-1]), col="gray", 
#'   main=paste0("Sunspot AR results on logarithmic scale for p=", ar.order))
#' lines(log(mcmc$psd.median[-1]))
#' lines(log(mcmc$psd.p05[-1]),lty=2)
#' lines(log(mcmc$psd.p95[-1]),lty=2)
#' lines(log(mcmc$psd.u05[-1]),lty=3)
#' lines(log(mcmc$psd.u95[-1]),lty=3)
#' legend(x="topright", legend=c("periodogram", "pointwise median", 
#'   "pointwise CI", "uniform CI"), lty=c(1,1,2,3), col=c("gray", 1, 1, 1))
#' 
#' 
#' ##
#' ## Example 2: Fit an AR(p) model to high-peaked AR(1) data
#' ##
#' 
#' # Use this variable to set the autoregressive model order
#' ar.order <- 1
#'
#' n <- 256
#' data <- arima.sim(n=n, model=list(ar=0.95)) 
#' data <- data - mean(data)
#' psd_true <- psd_arma(pi*omegaFreq(n), ar=0.95, ma=numeric(0), sigma2=1)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_AR(data=data, Ntotal=20000, burnin=8000, ar.order=1)
#' 
#' # Plot spectral estimates
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts(pdgrm[-1], col="gray",
#'   main=paste0("AR(1) data AR results for p=", ar.order))
#' lines(mcmc$psd.median[-1])
#' lines(mcmc$psd.p05[-1],lty=2)
#' lines(mcmc$psd.p95[-1],lty=2)
#' lines(mcmc$psd.u05[-1],lty=3)
#' lines(mcmc$psd.u95[-1],lty=3)
#' lines(psd_true[-1],col=2)
#' legend(x="topright", legend=c("periodogram", "true psd", 
#'   "pointwise median", "pointwise CI", "uniform CI"), lty=c(1,1,1,2,3), 
#'   col=c("gray", "red", 1, 1, 1))
#' 
#' # Compute the Integrated Absolute Error (IAE) of posterior median
#' cat("IAE=", mean(abs(mcmc$psd.median-psd_true)[-1]) , sep="")
#' }
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline lines plot title
#' @importFrom stats ar dbeta dgamma dnorm dt fft mad median plot.ts quantile rbeta rgamma rnorm rt runif sd var
#' @importFrom utils head
#' @useDynLib beyondWhittle, .registration = TRUE
#' @export
gibbs_AR <- function(data,
                     Ntotal,
                     burnin,
                     ar.order,
                     sigma2.alpha=0.001,
                     sigma2.beta=0.001) {
  
  stopifnot(ar.order > 0)
  
  # Tolerance for mean centering
  tol <- 1e-4
  
  # Mean center - necessary for our implementation of FFT
  if (abs(mean(data)) > tol) {
    data <- data - mean(data)
    warning("Data has been mean centered for your convenience")
  }
  
  mcmc_ar <- gibbs_pacf(data=data,
                        Ntotal=Ntotal,
                        burnin=burnin,
                        Nadaptive=burnin,
                        adaption.batchSize=50,
                        adaption.targetAcceptanceRate=0.44,
                        ar.order=ar.order,
                        sigma2.alpha=sigma2.alpha,
                        sigma2.beta=sigma2.beta,
                        mu.prop=rep(0,ar.order),
                        var.prop=rep(1/length(data),ar.order),
                        psi.alpha=rep(1,ar.order),
                        psi.beta=rep(1,ar.order))
  
  ##
  ## Construct spectral estimates and credible regions
  ##
  N_MCMC_IT <- Ntotal-burnin
  n <- length(data)
  lambda <- pi*omegaFreq(n)
  N <- length(lambda)
  fpsd_store <- log_fpsd_store <- array(data=NA, dim=c(N_MCMC_IT, N))
  for (i in 1:N_MCMC_IT) {
    if (ar.order == 1) {
      ar_store <- mcmc_ar$psi[i]
    }
    if (ar.order > 1) {
      ar_store <- pacfToAR(mcmc_ar$psi[,i])
    }
    sigma2_store <- mcmc_ar$sigma2[i]
    fpsd_store[i,] <- psd_arma(lambda,ar_store,numeric(0),sigma2_store)
    log_fpsd_store[i,] <- logfuller(fpsd_store[i,])
  }
  psd.median <- apply(fpsd_store, 2, median)
  psd.mean <- apply(fpsd_store, 2, mean)
  
  psd.p05 <- apply(fpsd_store, 2, quantile, 0.05)
  psd.p95 <- apply(fpsd_store, 2, quantile, 0.95)
  
  log.fpsd.s <- apply(log_fpsd_store, 2, median)
  log.fpsd.mad <- apply(log_fpsd_store, 2, mad)
  log.fpsd.help <- apply(log_fpsd_store, 2, uniformmax)
  log.Cvalue <- quantile(log.fpsd.help, 0.9)
  psd.u05 <- exp(log.fpsd.s - log.Cvalue * log.fpsd.mad)
  psd.u95 <- exp(log.fpsd.s + log.Cvalue * log.fpsd.mad)
  
  return(structure(list(rho=mcmc_ar$psi,
                        sigma2=mcmc_ar$sigma2,
                        DIC=mcmc_ar$DIC,
                        psd.median=psd.median,
                        psd.mean=psd.mean,
                        psd.p05=psd.p05,
                        psd.p95=psd.p95,
                        psd.u05=psd.u05,
                        psd.u95=psd.u95),
                   class="gibbs_AR"))
}
