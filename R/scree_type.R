#' Negative log AR likelihood values for scree-type plots
#'
#' (Approximate) negative maximum log-likelihood for for different autoregressive orders to produce scree-type plots.
#' @details By default, the maximum likelihood is approximated by the Yule-Walker method, due to numerical stabililty and computational speed. Further details can be found in the simulation study section in the referenced paper.
#' @param data numeric vector of data
#' @param order.max maximum autoregressive order to consider
#' @param method character string giving the method used to fit the model, to be forwarded to \code{stats::\link{ar}}
#' @return a data frame containing the autoregressive orders \code{p} and the corresponding negative log likelihood values \code{nll}
#' @references C. Kirch et al. (2018)
#' \emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
#' Bayesian Analysis
#' <doi:10.1214/18-BA1126>
#' @examples 
#' \dontrun{
#' 
#' ###
#' ### Interactive visual inspection for the sunspot data
#' ###
#' 
#' data <- sqrt(as.numeric(sunspot.year))
#' data <- data <- data - mean(data)
#' 
#' screeType <- scree_type_ar(data, order.max=15)
#' 
#' # Determine the autoregressive order by an interactive visual inspection of the scree-type plot
#' plot(x=screeType$p, y=screeType$nll, type="b")
#' p_ind <- identify(x=screeType$p, y=screeType$nll, n=1, labels=screeType$p)
#' print(screeType$p[p_ind])
#' }
#' @export
scree_type_ar <- function(data, order.max, method="yw") {
  stopifnot(order.max > 0)
  if (abs(mean(data)) > 1e-4) {
    data <- data - mean(data)
    warning("Data has been mean centered for your convenience")
  }
  p_vals <- 0:order.max
  nll_val <- rep(NA, length(p_vals))
  for (p in p_vals) {
    if (p==0) {
      nll_val[p+1] <- -sum(dnorm(data, mean=0, sd=sd(data), log=T))
    }
    if (p > 0) {
      ar_p <- ar(data, aic=F, order.max=p)
      if (length(ar_p$ar) != p) {
        stop(paste0("Something went wrong with the AR estimation for p=", p))
      }
      nll_val[p+1] <- -lik_ar(data, ar_p$ar, mean=0, sd=sqrt(ar_p$var.pred), log=T, full=T)
    }
  }
  return(data.frame(p=p_vals, nll=nll_val))
}


#' Likelihood of an autoregressive time series model with i.i.d. normal innovations
#' @param x numeric vector of data
#' @param ar vector of ar parameters
#' @param mean the innovation mean
#' @param sd the innovation standard deviation
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @param full logical; if TRUE, the full likelihood for all observations is computed; if FALSE, the partial likelihood for the last n-p observations
#' @return numeric value for the likelihood or log-likelihood
#' @keywords internal
lik_ar <- function(x, ar, mean=0, sd=1, log=F, full=T) {
  stopifnot(length(x) > 0)
  stopifnot(length(mean) == 1)
  stopifnot(length(sd) == 1)
  epsilon_t <- genEpsARMAC(x, ar=ar, ma=numeric(0))
  cll <- sum(dnorm(epsilon_t, mean=mean, sd=sd, log=T)) # conditional part
  p <- length(ar)
  if (length(p) > 0 && full) {
    zt_p <- head(x, p)
    gamma_p <- ltsa::tacvfARMA(phi=ar, maxLag=p-1, sigma2=sqrt(sd))
    Gamma_p <- acvMatrix(gamma_p)
    Gamma_p_inv <- solve(Gamma_p)
    mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
    mll <- -log(2*pi)*p/2 -log(det(Gamma_p)) / 2 + mll_unnormalized # marginal part
  } else {
    mll <- 0
  }
  res <- cll + mll
  if (!log) {
    res <- exp(res)
  }
  return(as.numeric(res))
}