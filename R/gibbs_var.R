#' Gibbs sampler for vector autoregressive model.
#'
#' Obtain samples of the posterior of a Bayesian VAR model of fixed order.
#' An independent Normal-Inverse-Wishart prior is employed.
#' @details See Section 2.2.3 in Koop and Korobilis (2010) or Section 6.2 in Meier (2018) for further details
#' @param data numeric matrix; NA values are interpreted as missing values and treated as random
#' @param ar.order order of the autoregressive model (integer >= 0)
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param print_interval Number of iterations, after which a status is printed to console
#' @param full_lik logical; if TRUE, the full likelihood for all observations is used; if FALSE, the partial likelihood for the last n-p observations
#' @param beta.mu prior mean of beta vector (normal)
#' @param beta.Sigma prior covariance matrix of beta vector
#' @param Sigma.S prior parameter for the innovation covariance matrix, symmetric positive definite matrix
#' @param Sigma.nu prior parameter for the innovation covariance matrix, nonnegative real number
#' @return list containing the following fields:
#'
#'    \item{beta}{matrix containing traces of the VAR parameter vector beta}
#'    \item{Sigma}{trace of innovation covariance Sigma}
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise, componentwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{uniform credibility interval, see (6.5) in Meier (2018)}
#'    \item{lpost}{trace of log posterior}
#' @references G. Koop and D. Korobilis (2010)
#' \emph{Bayesian Multivariate Time Series Methods for Empirical Macroeconomics}
#' Foundations and Trends in Econometrics
#' <doi:10.1561/0800000013>
#' @references A. Meier (2018)
#' \emph{A Matrix Gamma Process and Applications to Bayesian Analysis of Multivariate Time Series}
#' PhD thesis, OvGU Magdeburg
#' <https://opendata.uni-halle.de//handle/1981185920/13470>
#' @examples 
#' \dontrun{
#' 
#' ##
#' ## Example 1: Fit a VAR(p) model to SOI/Recruitment series:
#' ##
#' 
#' # Use this variable to set the VAR model order
#' p <- 5
#' 
#' data <- cbind(as.numeric(astsa::soi-mean(astsa::soi)), 
#'               as.numeric(astsa::rec-mean(astsa::rec)) / 50)
#' data <- apply(data, 2, function(x) x-mean(x))
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_var(data=data, ar.order=p, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Plot spectral estimate, credible regions and periodogram on log-scale
#' plot(mcmc, log=T)
#' 
#' 
#' 
#' ##
#' ## Example 2: Fit a VAR(p) model to VMA(1) data
#' ##
#' 
#' # Use this variable to set the VAR model order
#' p <- 5
#' 
#' n <- 256
#' ma <- rbind(c(-0.75, 0.5), c(0.5, 0.75))
#' Sigma <- rbind(c(1, 0.5), c(0.5, 1))
#' data <- sim_varma(model=list(ma=ma), n=n, d=2)
#' data <- apply(data, 2, function(x) x-mean(x))
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' mcmc <- gibbs_var(data=data, ar.order=p, Ntotal=10000, burnin=4000, thin=2)
#' 
#' # Plot spectral estimate, credible regions and periodogram on log-scale
#' plot(mcmc, log=T)
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib beyondWhittle, .registration = TRUE
#' @export
gibbs_var <- function(data,
                     ar.order,
                     Ntotal,
                     burnin,
                     thin=1,
                     print_interval=500,
                     full_lik=F,
                     beta.mu=rep(0,ar.order * ncol(data)^2),
                     beta.Sigma=1e4 * diag(ar.order * ncol(data)^2),
                     Sigma.S=1e-4 * diag(ncol(data)),
                     Sigma.nu=1e-4) {
  
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("'data' must be numeric matrix with d columns and n rows")
  }
  
  d <- ncol(data)
  if (d<2) {
    stop("This function is not suited for univariate time series. Use gibbs_AR instead")
  }

  if (max(abs(apply(data,2,mean,na.rm=T))) > 1e-4) {
    data <- apply(data,2,center,na.rm=T)
    warning("Data has been mean centered")
  }
  
  cl <- match.call()
  
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval)
  prior_params <- list(var.order=ar.order,
                    beta_prior=beta.mu,
                    V_prior=beta.Sigma,
                    S_prior=Sigma.S,
                    nu_prior=Sigma.nu)
  model_params <- psd_dummy_model()
  
  # Call internal MCMC algorithm
  mcmc_var <- gibbs_VAR_nuisance_intern(data=data,
                               mcmc_params=mcmc_params,
                               prior_params=prior_params,
                               model_params=model_params)

  return(structure(list(call=cl,
                        data=data,
                        beta=mcmc_var$beta,
                        Sigma=mcmc_var$Sigma,
                        psd.median=complexValuedPsd(mcmc_var$fpsd.s),
                        psd.mean=complexValuedPsd(mcmc_var$fpsd.mean),
                        psd.p05=complexValuedPsd(mcmc_var$fpsd.s05),
                        psd.p95=complexValuedPsd(mcmc_var$fpsd.s95),
                        psd.u05=complexValuedPsd(mcmc_var$fpsd.uci05),
                        psd.u95=complexValuedPsd(mcmc_var$fpsd.uci95),
                        missing_values=mcmc_var$missingValues_trace,
                        lpost=mcmc_var$lpost,
                        algo="gibbs_var"),
                   class="gibbs_psd"))
}
