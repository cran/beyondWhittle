#' Gibbs sampler for Bayesian semiparametric inference with the corrected AR likelihood
#'
#' Obtain samples of the posterior of the corrected autoregressive likelihood in conjunction with a Bernstein-Dirichlet prior on the correction.
#' @details Partial Autocorrelation Structure (PACF, uniform prior) and the residual variance sigma2 (inverse gamma prior) is used as model parametrization.
#' A Bernstein-Dirichlet prior for c_eta with base measure Beta(g0.alpha, g0.beta) is used.
#' Further details can be found in the simulation study section in the referenced paper by Kirch et al.
#' For more information on the PACF parametrization, see the referenced paper by Barndorff-Nielsen and Schou.
#' @param data numeric vector; NA values are interpreted as missing values and treated as random
#' @param ar.order order of the autoregressive model (integer > 0)
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
#' @param eta logical variable indicating whether the model confidence eta 
#' should be included in the inference (eta=T) or fixed to 1 (eta=F)
#' @param M DP base measure constant (> 0)
#' @param g0.alpha,g0.beta parameters of Beta base measure of DP
#' @param k.theta prior parameter for polynomial degree k (propto exp(-k.theta*k*log(k)))
#' @param tau.alpha,tau.beta prior parameters for tau (inverse gamma)
#' @param trunc_l,trunc_r left and right truncation of Bernstein polynomial basis functions, 0<=trunc_l<trunc_r<=1
#' @param coars flag indicating whether coarsened or default bernstein polynomials are used (see Appendix E.1 in Ghosal and van der Vaart 2017)
#' @param kmax upper bound for polynomial degree of Bernstein-Dirichlet mixture (can be set to Inf, algorithm is faster with kmax<Inf due to pre-computation of basis functions, but values 500<kmax<Inf are very memory intensive)
#' @param L truncation parameter of DP in stick breaking representation
#' @return list containing the following fields:
#'
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{uniform credibility interval}
#'    \item{k,tau,V,W}{posterior traces of nonparametric correction}
#'    \item{rho}{posterior trace of model AR parameters (PACF parametrization)}
#'    \item{eta}{posterior trace of model confidence eta}
#'    \item{lpost}{trace of log posterior}
#' @references C. Kirch et al. (2018)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' Bayesian Analysis
#' <doi:10.1214/18-BA1126>
#' @references S. Ghosal and A. van der Vaart (2017)
#' \emph{Fundamentals of Nonparametric Bayesian Inference} <doi:10.1017/9781139029834>
#' @references O. Barndorff-Nielsen and G. Schou
#' On the parametrization of autoregressive models by partial autocorrelations
#' Journal of Multivariate Analysis (3),408-419
#' <doi:10.1016/0047-259X(73)90030-4>
#' @examples 
#' \dontrun{
#' 
#' ##
#' ## Example 1: Fit a nonparametrically corrected AR(p) model to sunspot data:
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
#' mcmc <- gibbs_npc(data=data, ar.order=p, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Plot spectral estimate, credible regions and periodogram on log-scale
#' plot(mcmc, log=T)
#' 
#' 
#' ##
#' ## Example 2: Fit a nonparametrically corrected AR(p) model to high-peaked AR(1) data
#' ##
#' 
#' # Use this variable to set the autoregressive model order
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
#' mcmc <- gibbs_npc(data=data, ar.order=p, Ntotal=10000, burnin=4000, thin=2)
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
gibbs_npc <- function(data,
                      ar.order,
                      Ntotal,
                      burnin,
                      thin=1,
                      print_interval=100,
                      numerical_thresh=1e-7,
                      adaption.N=burnin,
                      adaption.batchSize=50,
                      adaption.tar=.44,
                      full_lik=F,
                      rho.alpha=rep(1,ar.order),
                      rho.beta=rep(1,ar.order),
                      eta=T,
                      M=1,
                      g0.alpha=1,
                      g0.beta=1,
                      k.theta=0.01,
                      tau.alpha=0.001,
                      tau.beta=0.001,
                      trunc_l=0.1,
                      trunc_r=0.9,
                      coars=F,
                      kmax = 100*coars + 500*(!coars),
                      L = max(20, length(data) ^ (1 / 3))) {
  if (ar.order == 0) {
    stop("'ar.order' must be greater than 0. Use gibbs_np instead.")
  }
  
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
                      Nadaptive=adaption.N,
                      adaption.batchSize=adaption.batchSize,
                      adaption.targetAcceptanceRate=adaption.tar,
                      numerical_thresh=numerical_thresh) 
  prior_params <- list(M=M,
                       g0.alpha=g0.alpha,
                       g0.beta=g0.beta,
                       k.theta=k.theta,
                       tau.alpha=tau.alpha,
                       tau.beta=tau.beta,
                       kmax=kmax,
                       bernstein_l=trunc_l,
                       bernstein_r=trunc_r,
                       bernstein_coars=coars,
                       L=L,
                       toggle=T,
                       alpha.toggle=eta,
                       prior.q=T,
                       ar.order=ar.order,
                       rho.alpha=rho.alpha,
                       rho.beta=rho.beta)
  model_params <- psd_dummy_model() 
  
  # Call internal MCMC algorithm
  mcmc_NPC <- gibbs_nuisance(data=data, 
                             mcmc_params=mcmc_params, 
                             corrected=T, 
                             prior_params=prior_params, 
                             model_params=model_params,
                             full_lik=full_lik)
  return(structure(list(call=cl,
                        data=data,
                        psd.median=mcmc_NPC$fpsd.s,
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
                        eta=mcmc_NPC$f.alpha,
                        missing_values=mcmc_NPC$missingValues_trace,
                        lpost=mcmc_NPC$lpostTrace,
                        algo="gibbs_npc"),
                   class="gibbs_psd"))
}
