#' Gibbs sampler for multivaiate Bayesian nonparametric inference with Whittle likelihood
#'
#' Obtain samples of the posterior of the multivariate Whittle likelihood in conjunction with an Hpd AGamma process prior on the spectral density matrix.
#' @details A detailed description of the method can be found in Section 5 in Meier (2018).
#' @param data numeric matrix; NA values are interpreted as missing values and treated as random
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param print_interval Number of iterations, after which a status is printed to console
#' @param numerical_thresh Lower (numerical pointwise) bound for the eigenvalues of the spectral density
#' @param adaption.N total number of iterations, in which the proposal variances (of r and U) are adapted
#' @param adaption.batchSize batch size of proposal adaption
#' @param adaption.tar target acceptance rate for adapted parameters
#' @param eta AGamma process parameter, real number > ncol(data)-1
#' @param omega AGamma process parameter, positive constant
#' @param Sigma AGamma process parameter, Hpd matrix 
#' @param k.theta prior parameter for polynomial degree k (propto exp(-k.theta*k*log(k)))
#' @param kmax upper bound for polynomial degree of Bernstein-Dirichlet mixture (can be set to Inf, algorithm is faster with kmax<Inf due to pre-computation of basis functions, but values 500<kmax<Inf are very memory intensive)
#' @param trunc_l,trunc_r left and right truncation of Bernstein polynomial basis functions, 0<=trunc_l<trunc_r<=1
#' @param coars flag indicating whether coarsened or default bernstein polynomials are used (see Appendix E.1 in Ghosal and van der Vaart 2017)
#' @param L truncation parameter of Gamma process
#' @return list containing the following fields:
#'
#'    \item{r,x,U}{traces of the AGamma process parameters}
#'    \item{k}{posterior trace of polynomial degree}
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise, componentwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{uniform credibility interval, see (6.5) in Meier (2018)}
#'    \item{lpost}{trace of log posterior}
#' @references A. Meier (2018)
#' \emph{A Matrix Gamma Process and Applications to Bayesian Analysis of Multivariate Time Series}
#' PhD thesis, OvGU Magdeburg
#' <https://opendata.uni-halle.de//handle/1981185920/13470>
#' @examples 
#' \dontrun{
#' 
#' ##
#' ## Example: Fit multivariate NP model to SOI/Recruitment series:
#' ##
#' 
#' data <- cbind(as.numeric(astsa::soi-mean(astsa::soi)), 
#'               as.numeric(astsa::rec-mean(astsa::rec)) / 50)
#' data <- apply(data, 2, function(x) x-mean(x))
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_vnp(data=data, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Visualize results
#' plot(mcmc, log=T)
#' 
#' 
#' ##
#' ## Example 2: Fit multivariate NP model to VMA(1) data
#' ##
#' 
#' n <- 256
#' ma <- rbind(c(-0.75, 0.5), c(0.5, 0.75))
#' Sigma <- rbind(c(1, 0.5), c(0.5, 1))
#' data <- sim_varma(model=list(ma=ma), n=n, d=2)
#' data <- apply(data, 2, function(x) x-mean(x))
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_vnp(data=data, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Plot spectral estimate, credible regions and periodogram on log-scale
#' plot(mcmc, log=T)
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib beyondWhittle, .registration = TRUE
#' @export
gibbs_vnp <- function(data,
                      Ntotal,
                      burnin,
                      thin=1,
                      print_interval=100,
                      numerical_thresh=1e-7,
                      adaption.N=burnin,
                      adaption.batchSize=50,
                      adaption.tar=.44,
                      eta=ncol(data),
                      omega=ncol(data),
                      Sigma=1e4*diag(ncol(data)),
                      k.theta=0.01,
                      kmax = 100*coars + 500*(!coars),
                      trunc_l = 0.1,
                      trunc_r = 0.9,
                      coars=F,
                      L = max(20, length(data) ^ (1 / 3))) {
  
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("'data' must be numeric matrix with d columns and n rows")
  }
  
  d <- ncol(data)
  if (d<2) {
    stop("This function is not suited for univariate time series. Use gibbs_NP instead")
  }
  
  if (max(abs(apply(data,2,mean,na.rm=T))) > 1e-4) {
    data <- apply(data,2,center,na.rm=T)
    warning("Data has been mean centered")
  }
  
  if (eta <= d-1) {
    stop("eta must be a number greater than d-1")
  }
  if (omega <= 0) {
    stop("omega must be a positive number")
  }
  if (!(is.matrix(Sigma)) || (!is_hpd(Sigma)) || any(dim(Sigma)!=c(d,d))) {
    stop("Sigma must be a Hermitian positive definite d times d matrix")
  }
  
  cl <- match.call()
  
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval,
                      numerical_thresh=numerical_thresh,
                      verbose=F,
                      Nadaptive=adaption.N,
                      adaption.batchSize=adaption.batchSize,
                      adaption.targetAcceptanceRate=adaption.tar)
  prior_params <- list(eta=eta,
                       omega=omega,
                       Sigma=Sigma,
                       k.theta=k.theta,
                       kmax=kmax, 
                       bernstein_l=trunc_l, # Note
                       bernstein_r=trunc_r, # Note
                       coarsened=coars,
                       L=L)
  model_params <- psd_dummy_model() 
  
  # Call internal MCMC algorithm
  mcmc_VNP <- gibbs_multivariate_nuisance(data=data, 
                                           mcmc_params=mcmc_params, 
                                           corrected=F, 
                                           prior_params=prior_params, 
                                           model_params=model_params)
  #return(mcmc_VNP)
  return(structure(list(call=cl,
                        data=data,
                        psd.median=complexValuedPsd(mcmc_VNP$fpsd.s),
                        psd.p05=complexValuedPsd(mcmc_VNP$fpsd.s05),
                        psd.p95=complexValuedPsd(mcmc_VNP$fpsd.s95),
                        psd.mean=complexValuedPsd(mcmc_VNP$fpsd.mean),
                        psd.u05=complexValuedPsd(mcmc_VNP$fpsd.uuci05),
                        psd.u95=complexValuedPsd(mcmc_VNP$fpsd.uuci95),
                        k=mcmc_VNP$k,
                        r=mcmc_VNP$r,
                        x=mcmc_VNP$Z,
                        U=mcmc_VNP$U,
                        missing_values=mcmc_VNP$missingValues_trace,
                        lpost=mcmc_VNP$lpostTrace,
                        algo="gibbs_vnp"),
                   class="gibbs_psd"))
}
