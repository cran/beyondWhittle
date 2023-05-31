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
#'    \item{lpost}{trace of log posterior}
#' @references C. Kirch et al. (2018)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' Bayesian Analysis
#' <doi:10.1214/18-BA1126>
#' @references N. Choudhuri et al. (2004)
#' \emph{Bayesian Estimation of the Spectral Density of a Time Series}
#' JASA
#' <doi:10.1198/016214504000000557>
#' @references S. Ghosal and A. van der Vaart (2017)
#' \emph{Fundamentals of Nonparametric Bayesian Inference} <doi:10.1017/9781139029834>
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
#' mcmc <- gibbs_np(data=data, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Plot spectral estimate, credible regions and periodogram on log-scale
#' plot(mcmc, log=T)
#' 
#' 
#' ##
#' ## Example 2: Fit the NP model to high-peaked AR(1) data
#' ##
#'
#' n <- 256
#' data <- arima.sim(n=n, model=list(ar=0.95)) 
#' data <- data - mean(data)
#' omega <- fourier_freq(n)
#' psd_true <- psd_arma(omega, ar=0.95, ma=numeric(0), sigma2=1)
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_np(data=data, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Compare estimate with true function (green)
#' plot(mcmc, log=F, pdgrm=F, credib="uniform")
#' lines(x=omega, y=psd_true, col=3, lwd=2)
#' 
#' # Compute the Integrated Absolute Error (IAE) of posterior median
#' cat("IAE=", mean(abs(mcmc$psd.median-psd_true)[-1]) , sep="")
#' }
#' @importFrom Rcpp evalCpp
#' @importFrom stats is.ts
#' @useDynLib beyondWhittle, .registration = TRUE
#' @export
gibbs_np <- function(data,
                     Ntotal,
                     burnin,
                     thin=1,
                     print_interval=100,
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
  model_params <- psd_dummy_model() 
  
  # Call internal MCMC algorithm
  mcmc_NP <- gibbs_nuisance(data=data, 
                            mcmc_params=mcmc_params, 
                            corrected=F, 
                            prior_params=prior_params, 
                            model_params=model_params)
  return(structure(list(call=cl,
                        data=data,
                        psd.median=mcmc_NP$fpsd.s,
                        psd.p05=mcmc_NP$fpsd.s05,
                        psd.p95=mcmc_NP$fpsd.s95,
                        psd.mean=mcmc_NP$fpsd.mean,
                        psd.u05=mcmc_NP$log.conflower,
                        psd.u95=mcmc_NP$log.confupper,
                        k=mcmc_NP$k,
                        tau=mcmc_NP$tau,
                        V=mcmc_NP$V,
                        W=mcmc_NP$W,
                        missing_values=mcmc_NP$missingValues_trace,
                        lpost=mcmc_NP$lpostTrace,
                        algo="gibbs_np"),
                   class="gibbs_psd"))
}
