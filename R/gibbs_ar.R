#' Gibbs sampler for an autoregressive model with PACF parametrization.
#'
#' Obtain samples of the posterior of a Bayesian autoregressive model of fixed order.
#' @details Partial Autocorrelation Structure (PACF, uniform prior) and the residual variance sigma2 (inverse gamma prior) is used as model parametrization.
#' The DIC is computed with two times the posterior variance of the deviance as effective number of parameters, see (7.10) in the referenced book by Gelman et al.
#' Further details can be found in the simulation study section in the referenced paper by C. Kirch et al.
#' For more information on the PACF parametrization, see the referenced paper by Barndorff-Nielsen and Schou.
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
#'    \item{lpost}{trace of log posterior}
#' @references C. Kirch et al. (2018)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' Bayesian Analysis
#' <doi:10.1214/18-BA1126>
#' @references A. Gelman et al. (2013)
#' \emph{Bayesian Data Analysis, Third Edition}
#' @references O. Barndorff-Nielsen and G. Schou
#' On the parametrization of autoregressive models by partial autocorrelations
#' Journal of Multivariate Analysis (3),408-419
#' <doi:10.1016/0047-259X(73)90030-4>
#' @examples 
#' \dontrun{
#' 
#' ##
#' ## Example 1: Fit an AR(p) model to sunspot data:
#' ##
#' 
#' # Use this variable to set the AR model order
#' p <- 2
#' 
#' data <- sqrt(as.numeric(sunspot.year))
#' data <- data - mean(data)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_ar(data=data, ar.order=p, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Plot spectral estimate, credible regions and periodogram on log-scale
#' plot(mcmc, log=T)
#' 
#' 
#' ##
#' ## Example 2: Fit an AR(p) model to high-peaked AR(1) data
#' ##
#' 
#' # Use this variable to set the AR model order
#' p <- 1
#'
#' n <- 256
#' data <- arima.sim(n=n, model=list(ar=0.95)) 
#' data <- data - mean(data)
#' omega <- fourier_freq(n)
#' psd_true <- psd_arma(omega, ar=0.95, ma=numeric(0), sigma2=1)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_ar(data=data, ar.order=p, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Compare estimate with true function (green)
#' plot(mcmc, log=F, pdgrm=F, credib="uniform")
#' lines(x=omega, y=psd_true, col=3, lwd=2)
#' 
#' # Compute the Integrated Absolute Error (IAE) of posterior median
#' cat("IAE=", mean(abs(mcmc$psd.median-psd_true)[-1]) , sep="")
#' }
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline lines plot title
#' @importFrom stats is.ts ar dbeta dgamma dnorm dt fft mad median plot.ts quantile rbeta rgamma rnorm rt runif sd var
#' @importFrom utils head
#' @useDynLib beyondWhittle, .registration = TRUE
#' @export
gibbs_ar <- function(data,
                     ar.order,
                     Ntotal,
                     burnin,
                     thin=1,
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
  
  stopifnot(ar.order >= 0)
  
  if (!is.ts(data) && (!is.vector(data) || !is.numeric(data))) {
    stop("'data' must be a numeric vector or time series")
  }
  
  if (abs(mean(data,na.rm=T)) > 1e-4) {
    data <- data - mean(data,na.rm=T)
    warning("Data has been mean centered")
  }
  
  cl <- match.call()
  
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
  
  # Call internal MCMC algorithm
  mcmc_ar <- gibbs_AR_nuisance_intern(data=data,
                               mcmc_params=mcmc_params,
                               prior_params=prior_params,
                               model_params=model_params,
                               full_lik=full_lik)

  return(structure(list(call=cl,
                        data=data,
                        rho=mcmc_ar$psi,
                        sigma2=mcmc_ar$sigma2,
                        DIC=mcmc_ar$DIC,
                        psd.median=mcmc_ar$fpsd.s,
                        psd.mean=mcmc_ar$fpsd.mean,
                        psd.p05=mcmc_ar$fpsd.s05,
                        psd.p95=mcmc_ar$fpsd.s95,
                        psd.u05=mcmc_ar$log.conflower,
                        psd.u95=mcmc_ar$log.confupper,
                        missing_values=mcmc_ar$missingValues_trace,
                        lpost=mcmc_ar$lpostTrace,
                        algo="gibbs_ar"),
                   class="gibbs_psd"))
}
