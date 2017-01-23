#' Gibbs sampler for Bayesian semiparametric inference with the corrected AR likelihood
#'
#' Obtain samples of the posterior of the corrected autoregressive likelihood in conjunction with a Bernstein-Dirichlet prior on the correction.
#' @details Partial Autocorrelation Structure (PACF, uniform prior) and the residual variance sigma2 (inverse gamma prior) is used as model parametrization.
#' A Bernstein-Dirichlet prior for c_eta with base measure Beta(g0.alpha, g0.beta) is used.
#' Further details can be found in the simulation study section in the referenced paper.
#' @param data numeric vector
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param ar.order order of the autoregressive model (integer > 0)
#' @param eta logical variable indicating whether the model confidence eta 
#' should be included in the inference (eta=T) or fixed to 1 (eta=F)
#' @param M DP base measure constant (> 0)
#' @param g0.alpha,g0.beta parameters of Beta base measure of DP
#' @param k.theta prior parameter for polynomial degree k (propto exp(-k.theta*k*log(k)))
#' @param tau.alpha,tau.beta prior parameters for tau (inverse gamma)
#' @param kmax upper bound for polynomial degree of Bernstein-Dirichlet mixture
#' @param L truncation parameter of DP in stick breaking representation
#' @return list containing the following fields:
#'
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{uniform credibility interval}
#'    \item{k,tau,V,W}{posterior traces of nonparametric correction}
#'    \item{rho}{posterior trace of model AR parameters (PACF parametrization)}
#'    \item{eta}{posterior trace of model confidence eta}
#' @references C. Kirch, R. Meyer et al.
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' <arXiv:1701.04846>
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
#' data <- sqrt(as.numeric(sunspot.year[-length(sunspot.year)]))
#' data <- data - mean(data)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_NPC(data=data, Ntotal=50000, burnin=30000, thin=4,
#'   ar.order=ar.order)
#' 
#' # Plot spectral estimates on log-scale (excluding the zero frequency).
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts(log(pdgrm[-1]), col="gray", 
#'   main=paste0("Sunspot NPC results on logarithmic scale for p=", ar.order))
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
#' mcmc <- gibbs_NPC(data=data, Ntotal=50000, burnin=30000, thin=4,
#'   ar.order=ar.order)
#' 
#' # Plot spectral estimates
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts(pdgrm[-1], col="gray",
#'   main=paste0("AR(1) data NPC results for p=", ar.order))
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
#' @useDynLib beyondWhittle
#' @export
gibbs_NPC <- function(data,
                      Ntotal,
                      burnin,
                      thin,
                      ar.order,
                      eta=T,
                      M=1,
                      g0.alpha=1,
                      g0.beta=1,
                      k.theta=0.01,
                      tau.alpha=0.001,
                      tau.beta=0.001,
                      kmax = 500,
                      L = max(20, length(data) ^ (1 / 3))) {
  mcmc_NPC <- gibbs_toggle_ar(data=data,
                              Ntotal=Ntotal,
                              burnin=burnin,
                              thin=thin,
                              M=M,
                              g0.alpha=g0.alpha,
                              g0.beta=g0.beta,
                              k.theta=k.theta,
                              tau.alpha=tau.alpha,
                              tau.beta=tau.beta,
                              kmax=kmax,
                              L=L,
                              corrected=T,
                              prior.q=T,
                              alpha.toggle=eta,
                              toggle=T,
                              Nadaptive=burnin,
                              adaption.batchSize = 50,
                              adaption.targetAcceptanceRate = 0.44,
                              mu.prop=rep(0,ar.order),
                              var.prop.rho=rep(1/length(data),ar.order),
                              dist = "normal",
                              rho.alpha=rep(1, ar.order),
                              rho.beta=rep(1, ar.order))
  return(structure(list(psd.median=mcmc_NPC$fpsd.s,
                        psd.p05=mcmc_NPC$fpsd.s05,
                        psd.p95=mcmc_NPC$fpsd.s95,
                        psd.mean=mcmc_NPC$fpsd.mean,
                        psd.u05=mcmc_NPC$log.conflower,
                        psd.u95=mcmc_NPC$log.confupper,
                        k=mcmc_NPC$k,
                        tau=mcmc_NPC$tau,
                        V=mcmc_NPC$V,
                        W=mcmc_NPC$W,
                        rho=mcmc_NPC$rho,
                        eta=mcmc_NPC$f.alpha),
                   class="gibbs_NPC"))
}
