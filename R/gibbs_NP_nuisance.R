#' Bayesian nonparametric inference in nuisance model with Whittle likelihood
#'
#' Obtain samples of the posterior of the Bayesian model X_t=g(theta)+e_t,
#' where the time series part e_t is modeled nonparametrically (by the Whittle 
#' likelihood in conjunction with a Bernstein-Dirichlet prior on the spectral density) 
#' and the link function g (depending on the parameter of interest theta) is provided by the user.
#' Examples for this framework include the simple mean model X_t=mu+e_t (theta=mu, g_t(theta)=theta) or
#' the linear trend model X_t=bt+mu+e_t (theta=c(mu,b), g(theta)=bt+mu) for t=1,...,n.
#' @details See \link[beyondWhittle]{gibbs_NP} for further details on the
#' nonparametric time series part. See \link[beyondWhittle]{nuisanceModel_mean} or
#' \link[beyondWhittle]{nuisanceModel_linearTrend} for two examplary nuisance models.
#' @param data numeric vector; NA values are interpreted as missing values and treated as random
#' @param nuisanceModel model of class \code{nuisanceModel}, see \link[beyondWhittle]{nuisanceModel_mean}
#' and \link[beyondWhittle]{nuisanceModel_linearTrend} for examples
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
#'    \item{theta}{posterior traces of the parameter of interest}
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
#' ## Example 1: Fit Linear trend model with nonparametric nuisance time series to temperature data
#' ##
#' 
#' data <- as.numeric(nhtemp)
#' n <- length(data)
#' 
#' # Initialize the linear trend model
#' nuisanceModel <- nuisanceModel_linearTrend()
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' Ntotal <- 50000; burnin <- 20000; thin <- 4; Neffective <- (Ntotal-burnin)/thin
#' mcmc <- gibbs_NP_nuisance(data, nuisanceModel, Ntotal=Ntotal, 
#'                           burnin=burnin, thin=thin)
#' 
#' # Reconstruct linear trend lines from posterior traces of theta=c(mu,b)
#' trend_lines <- array(NA, dim=c(n,Neffective))
#' for (j in 1:Neffective) {
#'   mu_j <- mcmc$theta[1,j]
#'   b_j <- mcmc$theta[2,j]
#'   trend_lines[,j] <- mu_j + b_j * (1:n)
#' }
#' trend.median <- apply(trend_lines, 1, median)
#' trend.p05 <- apply(trend_lines, 1, quantile, 0.05) # 90% CI
#' trend.p95 <- apply(trend_lines, 1, quantile, 0.95) # "
#' trend.p025 <- apply(trend_lines, 1, quantile, 0.025) # 95% CI
#' trend.p975 <- apply(trend_lines, 1, quantile, 0.975) # 
#' 
#' # Plot confidence bands for linear trend curve
#' par(mfcol=c(2,1),mar=c(4,4,2,2))
#' data_timePeriod <- start(nhtemp)[1]:end(nhtemp)[1]
#' plot(x=data_timePeriod, y=data, 
#'      main="New Haven temperature w/ NP linear trend estimate", 
#'      col="gray", type="l", xlab="Year", ylab="Avg temp. (deg. F)")
#' lines(x=data_timePeriod, y=trend.median)
#' lines(x=data_timePeriod, y=trend.p05, lty=2)
#' lines(x=data_timePeriod, y=trend.p95, lty=2)
#' lines(x=data_timePeriod, y=trend.p025, lty=3)
#' lines(x=data_timePeriod, y=trend.p975, lty=3)
#' legend(x="topleft", legend=c("data", "posterior median", 
#'                          "posterior 90% CI", "posterior 95% CI"), 
#'        lty=c(1,1,2,3), col=c("gray", 1, 1, 1), 
#'        cex=.75, x.intersp=.5, y.intersp=.5)
#' 
#' # Plot spectral estimate
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts((pdgrm[-1]), col="gray", 
#'         main="New Haven temperature NP spectral estimate")
#' lines((mcmc$psd.median[-1]))
#' lines((mcmc$psd.p05[-1]),lty=2)
#' lines((mcmc$psd.p95[-1]),lty=2)
#' lines((mcmc$psd.u05[-1]),lty=3)
#' lines((mcmc$psd.u95[-1]),lty=3)
#' legend(x="top", legend=c("periodogram", "pointwise median", 
#'                               "pointwise CI", "uniform CI"), 
#'        lty=c(1,1,2,3), col=c("gray", 1, 1, 1),
#'        cex=.75, x.intersp=.5, y.intersp=.5)
#' par(mfcol=c(1,1))
#' 
#' 
#' ##
#' ## Example 2: Fit mean (with nonparametric nuisance time series) model to AR(1) data
#' ##
#' 
#' n <- 256
#' mu_true <- 3
#' data <- arima.sim(n=n, model=list(ar=0.75)) 
#' data <- data + mu_true
#' psd_true <- psd_arma(pi*omegaFreq(n), ar=0.95, ma=numeric(0), sigma2=1)
#' 
#' # Initialize the mean model
#' nuisanceModel <- nuisanceModel_mean()
#' 
#' # If you run the example be aware that this may take several minutes
#' print("example may take some time to run")
#' Ntotal <- 50000; burnin <- 20000; thin <- 4; Neffective <- (Ntotal-burnin)/thin
#' mcmc <- gibbs_NP_nuisance(data, nuisanceModel, Ntotal=Ntotal, 
#'                           burnin=burnin, thin=thin)
#' 
#' # Plot posterior trace of mu
#' par(mfcol=c(2,1),mar=c(4,2,2,2))
#' plot.ts(mcmc$theta[1,], main="Posterior trace of mu", xlab="Iteration")
#' abline(h=mu_true, col=2, lty=2)
#' abline(h=mean(data), col=3, lty=2)
#' legend(x="topright", legend=c("true mean", "empirical mean of data"), lty=c(2,2), col=c(2,3))
#' 
#' # Plot spectral estimate
#' N <- length(mcmc$psd.median)
#' pdgrm <- (abs(fft(data))^2 / (2*pi*length(data)))[1:N]
#' plot.ts((pdgrm[-1]), col="gray", 
#'         main="AR(1) NP spectral estimate")
#' lines((mcmc$psd.median[-1]))
#' lines((mcmc$psd.p05[-1]),lty=2)
#' lines((mcmc$psd.p95[-1]),lty=2)
#' lines((mcmc$psd.u05[-1]),lty=3)
#' lines((mcmc$psd.u95[-1]),lty=3)
#' legend(x="top", legend=c("periodogram", "pointwise median", 
#'                          "pointwise CI", "uniform CI"), 
#'        lty=c(1,1,2,3), col=c("gray", 1, 1, 1))
#' par(mfcol=c(1,1))
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib beyondWhittle, .registration = TRUE
#' @keywords internal
gibbs_NP_nuisance <- function(data,
                              nuisanceModel,
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
  if (any(is.na(data))) {
    stop("Missing values in nuisance models not supported (yet)")
  }
  
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
  model_params <- nuisanceModel 
  mcmc_NP <- gibbs_nuisance(data=data, 
                            mcmc_params=mcmc_params, 
                            corrected=F, 
                            prior_params=prior_params, 
                            model_params=model_params)
  return(structure(list(theta=mcmc_NP$theta,
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
                        missing_values=mcmc_NP$missingValues_trace),
                   class="gibbs_NP_nuisance"))
}