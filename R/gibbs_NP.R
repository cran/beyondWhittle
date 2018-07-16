#' Gibbs sampler for Bayesian nonparametric inference with Whittle likelihood
#'
#' Obtain samples of the posterior of the Whittle likelihood in conjunction with a Bernstein-Dirichlet prior on the spectral density.
#' @details Further details can be found in the simulation study section in the references papers.
#' @param data numeric vector; NA values are interpreted as missing values and treated as random
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param print_interval Number of iterations, after which a status is printed to console
#' @param numerical_thresh Lower (numerical pointwise) bound for the spectral density
#' @param M DP base measure constant (> 0)
#' @param g0.alpha,g0.beta parameters of Beta base measure of DP
#' @param k.theta prior parameter for polynomial degree k (propto exp(-k.theta*k*log(k)))
#' @param tau.alpha,tau.beta prior parameters for tau (inverse gamma)
#' @param kmax upper bound for polynomial degree of Bernstein-Dirichlet mixture (can be set to Inf, algorithm is faster with kmax<Inf due to pre-computation of basis functions, but values 500<kmax<Inf are very memory intensive)
#' @param trunc_l,trunc_r left and right truncation of Bernstein polynomial basis functions, 0<=trunc_l<trunc_r<=1
#' @param coars flag indicating whether coarsened or default bernstein polynomials are used (see Appendix E.1 in Ghosal and van der Vaart 2017)
#' @param L truncation parameter of DP in stick breaking representation
#' @return list containing the following fields:
#'
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{uniform credibility interval}
#'    \item{k,tau,V,W}{posterior traces of PSD parameters}
#' @references C. Kirch et al. (2017)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' <arXiv:1701.04846>
#' @references N. Choudhuri et al. (2004)
#' \emph{Bayesian Estimation of the Spectral Density of a Time Series} <DOI:10.1198/016214504000000557>
#' @references S. Ghosal and A. van der Vaart (2017)
#' \emph{Fundamentals of Nonparametric Bayesian Inference} <DOI:10.1017/9781139029834>
#' @examples 
#' \dontrun{
#' 
#' ##
#' ## Example 1: Fit the NP model to sunspot data:
#' ##
#' 
#' data <- sqrt(as.numeric(sunspot.year))
#' data <- data - mean(data)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_NP(data=data, Ntotal=50000, burnin=20000, thin=4)
#' 
#' # Plot spectral estimates on log-scale (excluding the zero frequency).
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts(log(pdgrm[-1]), col="gray", 
#'   main=paste0("Sunspot NP results on logarithmic scale"))
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
#' ## Example 2: Fit the NP model to high-peaked AR(1) data
#' ##
#'
#' n <- 256
#' data <- arima.sim(n=n, model=list(ar=0.95)) 
#' data <- data - mean(data)
#' psd_true <- psd_arma(pi*omegaFreq(n), ar=0.95, ma=numeric(0), sigma2=1)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_NP(data=data, Ntotal=50000, burnin=20000, thin=4)
#' 
#' # Plot spectral estimates
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts(pdgrm[-1], col="gray",
#'   main=paste0("AR(1) data NP results"))
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
#' @useDynLib beyondWhittle, .registration = TRUE
#' @export
gibbs_NP <- function(data,
                     Ntotal,
                     burnin,
                     thin,
                     print_interval=500,
                     numerical_thresh=1e-7,
                     M=1,
                     g0.alpha=1,
                     g0.beta=1,
                     k.theta=0.01,
                     tau.alpha=0.001,
                     tau.beta=0.001,
                     kmax = 100*coars + 500*(!coars),
                     trunc_l = 0.1,
                     trunc_r = 0.9,
                     coars=F,
                     L = max(20, length(data) ^ (1 / 3))) {
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval,
                      numerical_thresh=numerical_thresh)
  prior_params <- list(M=M,
                       g0.alpha=g0.alpha,
                       g0.beta=g0.beta,
                       k.theta=k.theta,
                       tau.alpha=tau.alpha,
                       tau.beta=tau.beta,
                       kmax=kmax, 
                       bernstein_l=trunc_l, # Note
                       bernstein_r=trunc_r, # Note
                       bernstein_coars=coars,
                       L=L)
  model_params <- psd_dummy_model() # dummy model within nuisance context
  mcmc_NP <- gibbs_nuisance(data=data, 
                            mcmc_params=mcmc_params, 
                            corrected=F, 
                            prior_params=prior_params, 
                            model_params=model_params)
  return(structure(list(psd.median=mcmc_NP$fpsd.s,
                        psd.p05=mcmc_NP$fpsd.s05,
                        psd.p95=mcmc_NP$fpsd.s95,
                        psd.mean=mcmc_NP$fpsd.mean,
                        psd.u05=mcmc_NP$log.conflower,
                        psd.u95=mcmc_NP$log.confupper,
                        k=mcmc_NP$k,
                        tau=mcmc_NP$tau,
                        V=mcmc_NP$V,
                        W=mcmc_NP$W,
                        missing_values=mcmc_NP$missingValues_trace),
                   class="gibbs_NP"))
}
