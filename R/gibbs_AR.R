#' Gibbs sampler for an autoregressive model with PACF parametrization.
#'
#' Obtain samples of the posterior of a Bayesian autoregressive model of fixed order.
#' @details Partial Autocorrelation Structure (PACF, uniform prior) and the residual variance sigma2 (inverse gamma prior) is used as model parametrization.
#' The DIC is computed with two times the posterior variance of the deviance as effective number of parameters, see (7.10) in the referenced book by Gelman et al.
#' Further details can be found in the simulation study section in the referenced paper by C. Kirch et al.
#' @param data numeric vector; NA values are interpreted as missing values and treated as random
#' @param ar.order order of the autoregressive model (integer >= 0)
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param print_interval Number of iterations, after which a status is printed to console
#' @param numerical_thresh Lower (numerical pointwise) bound for the spectral density
#' @param adaption.N total number of iterations, in which the proposal variances (of rho) are adapted
#' @param adaption.batchSize batch size of proposal adaption for the rho_i's (PACF)
#' @param adaption.tar target acceptance rate for the rho_i's (PACF)
#' @param full_lik logical; if TRUE, the full likelihood for all observations is used; if FALSE, the partial likelihood for the last n-p observations
#' @param rho.alpha,rho.beta prior parameters for the rho_i's: 2*(rho-0.5)~Beta(rho.alpha,rho.beta), default is Uniform(-1,1)
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
#' mcmc <- gibbs_AR(data=data, Ntotal=20000, burnin=8000, thin=4, ar.order=ar.order)
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
#' mcmc <- gibbs_AR(data=data, Ntotal=20000, burnin=8000, thin=4, ar.order=1)
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
                     ar.order,
                     Ntotal,
                     burnin,
                     thin,
                     print_interval=500,
                     numerical_thresh=1e-7,
                     adaption.N=burnin,
                     adaption.batchSize=50,
                     adaption.tar=.44,
                     full_lik=F,
                     rho.alpha=rep(1,ar.order),
                     rho.beta=rep(1,ar.order),
                     sigma2.alpha=0.001,
                     sigma2.beta=0.001) {
  
  stopifnot(ar.order > 0)
  
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval, # Note
                      numerical_thresh=numerical_thresh,
                      Nadaptive=adaption.N,
                      adaption.batchSize=adaption.batchSize,
                      adaption.targetAcceptanceRate=adaption.tar)
  prior_params <- list(ar.order=ar.order,
                    rho.alpha=rho.alpha,
                    rho.beta=rho.beta,
                    sigma2.alpha=sigma2.alpha,
                    sigma2.beta=sigma2.beta)
  model_params <- psd_dummy_model()
  
  mcmc_ar <- gibbs_AR_nuisance_intern(data=data,
                               mcmc_params=mcmc_params,
                               prior_params=prior_params,
                               model_params=model_params,
                               full_lik=full_lik)

  return(structure(list(rho=mcmc_ar$psi,
                        sigma2=mcmc_ar$sigma2,
                        DIC=mcmc_ar$DIC,
                        psd.median=mcmc_ar$fpsd.s,
                        psd.mean=mcmc_ar$fpsd.mean,
                        psd.p05=mcmc_ar$fpsd.s05,
                        psd.p95=mcmc_ar$fpsd.s95,
                        psd.u05=mcmc_ar$log.conflower,
                        psd.u95=mcmc_ar$log.confupper,
                        missing_values=mcmc_ar$missingValues_trace),
                   class="gibbs_AR"))
}
