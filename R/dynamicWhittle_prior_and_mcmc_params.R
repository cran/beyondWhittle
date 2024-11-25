#' Generate a list of parameter values in prior elicitation
#'
#' @param M DP base measure constant (> 0)
#' @param g0.alpha,g0.beta parameters of Beta base measure of DP
#' @param k1.theta prior parameter for polynomial corresponding to rescaled time (propto exp(-k1.theta*k1*log(k1)))
#' @param k2.theta prior parameter for polynomial corresponding to rescaled frequency (propto exp(-k2.theta*k2*log(k2)))
#' @param tau.alpha,tau.beta prior parameters for tau (inverse gamma)
#' @param k1max upper bound of the degrees of Bernstein polynomial
#' corresponding to rescaled time (for pre-computation of basis functions)
#' @param k2max upper bound of the degrees of Bernstein polynomial
#' corresponding to rescaled frequency (for pre-computation of basis functions)
#' @param bernstein1_l,bernstein1_r left and right truncation of Bernstein polynomial basis functions
#' for rescaled time, 0<=bernstein1_l<bernstein1_r<=1
#' @param bernstein2_l,bernstein2_r left and right truncation of Bernstein polynomial basis functions
#' for rescaled frequency, 0<=bernstein2_l<bernstein2_r<=1
#' @param L number of terms in fast DP approximation
#' @return A list of prior parameter values
#' @export
#'
bdp_dw_prior_params_gen <- function(M = 1,
                                    g0.alpha = 1,
                                    g0.beta = 1,
                                    k1.theta = 0.01,
                                    k2.theta = 0.01,
                                    tau.alpha = 0.001,
                                    tau.beta = 0.001,
                                    k1max = 100,
                                    k2max = 100,
                                    L = 10,
                                    bernstein1_l = 0.1,
                                    bernstein1_r = 0.9,
                                    bernstein2_l = 0.1,
                                    bernstein2_r = 0.9){

  # PRIOR PAREMETERS CHECK
  stopifnot(!is.null(M)); stopifnot(M > 0)

  stopifnot(!is.null(g0.alpha) && !is.null(g0.beta)); stopifnot(g0.alpha > 0 && g0.beta > 0)

  stopifnot(!is.null(k1.theta)); stopifnot(k1.theta > 0)

  stopifnot(!is.null(k2.theta)); stopifnot(k2.theta > 0)

  stopifnot(!is.null(tau.alpha) && !is.null(tau.beta)); stopifnot(tau.alpha > 0 && tau.beta > 0)

  stopifnot(!is.null(k1max)); stopifnot(k1max > 0)

  stopifnot(!is.null(k2max)); stopifnot(k2max > 0)

  stopifnot(!is.null(L)); stopifnot(L > 0)

  stopifnot(!is.null(bernstein1_l) && !is.null(bernstein1_r)); stopifnot(bernstein1_l >= 0 && bernstein1_r <= 1)

  stopifnot(!is.null(bernstein2_l) && !is.null(bernstein2_r)); stopifnot(bernstein2_l >= 0 && bernstein2_r <= 1)



  prior_params <- list(M = M,
                       g0.alpha = g0.alpha, g0.beta = g0.beta,
                       k1.theta = k1.theta,
                       k2.theta = k2.theta,
                       tau.alpha = tau.alpha, tau.beta = tau.beta,
                       k1max = k1max,
                       k2max = k2max,
                       L = L,
                       bernstein1_l = bernstein1_l, bernstein1_r = bernstein1_r,
                       bernstein2_l = bernstein2_l, bernstein2_r = bernstein2_r)

  return(prior_params)

}



#' Generate a list of values for MCMC algorithm
#'
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (for post-processing of the posterior sample)
#' @param adaptive.batchSize the batch size for the adaptive MCMC algorithm for sampling tau
#' @param adaptive.targetAcceptanceRate the target acceptance rate for the adaptive MCMC algorithm for sampling tau
#' @return A list of MCMC parameter values
#' @export
#'
bdp_dw_mcmc_params_gen <- function(Ntotal = 1.1e5,
                                   burnin = 6e4,
                                   thin = 10,
                                   adaptive.batchSize = 50,
                                   adaptive.targetAcceptanceRate = 0.44){


  # MCMC parameters check
  stopifnot(!is.null(Ntotal)); stopifnot(Ntotal>0)

  stopifnot(!is.null(burnin)); stopifnot(burnin>=0 && burnin<Ntotal)

  stopifnot(!is.null(thin)); stopifnot(thin>=1)

  # Adaptive MCMC parameters check
  stopifnot(!is.null(adaptive.batchSize)); stopifnot(adaptive.batchSize>0)

  stopifnot(!is.null(adaptive.targetAcceptanceRate)); stopifnot(adaptive.targetAcceptanceRate>0 && adaptive.targetAcceptanceRate<1)




  mcmc_params <- list(Ntotal = Ntotal,
                      burnin = burnin,
                      thin = thin,
                      adaption.batchSize = adaptive.batchSize,
                      adaption.targetAcceptanceRate = adaptive.targetAcceptanceRate)

  return(mcmc_params)

}

#' Calculating the estimated posterior mean, median and credible region (tv-PSD)
#' @importFrom stats qgamma
#'
#' @param post_sample the posterior sample generated by \link{bdp_dw_mcmc}.
#' @param rescaled_time,freq numeric vectors forming a rectangular grid on which the estimated tv-PSD is evaluated.
#' @param unif_CR a Boolean value (default FALSE) indicating whether to calculate the uniform credible region
#' rescaled_time must be in \eqn{[0,1]} and freq must be in \eqn{[0,\pi]}.
#' @return list containing the following fields:
#' 
#'    \item{tvpsd.mean,tvpsd.median}{posterior mean and pointwise posterior median (matrices of dimension length(rescaled_time) by length(freq))}
#'    \item{tvpsd.p05,tvpsd.p95}{90 percent pointwise credibility interval}
#'    \item{tvpsd.u05,tvpsd.u95}{90 percent uniform credibility interval if unif_CR = TRUE. Otherwise NA}
#' @export
#'
bdp_dw_est_post_stats <- function(post_sample, rescaled_time, freq, unif_CR = FALSE){

  stopifnot(!is.null(rescaled_time)); stopifnot(all(rescaled_time >= 0 & rescaled_time <= 1))

  stopifnot(!is.null(freq)); stopifnot(all(freq >= 0 & freq <= pi))
  
  rescaled_freq <- freq/pi

  k1max <- max(post_sample$k1)
  k2max <- max(post_sample$k2)

  bernstein1_l <- (post_sample$prior_params)$bernstein1_l
  bernstein1_r <- (post_sample$prior_params)$bernstein1_r
  bernstein2_l <- (post_sample$prior_params)$bernstein2_l
  bernstein2_r <- (post_sample$prior_params)$bernstein2_r
  
  L <- (post_sample$prior_params)$L
  M <- (post_sample$prior_params)$M

  db.list_1 <- dbList_dw_Bern(rescaled_time, k1max, bernstein1_l, bernstein1_r)
  db.list_2 <- dbList_dw_Bern(rescaled_freq, k2max, bernstein2_l, bernstein2_r)

  post_funval <- array(NA, dim = c(length(rescaled_time), length(rescaled_freq), length(post_sample$k1)))

  for (j in 1:(dim(post_funval)[3])) {
    
    k1 <- (post_sample$k1)[j]
    k2 <- (post_sample$k2)[j]
    tau <- (post_sample$tau)[j]
    E <- (post_sample$E)[,j]
    W1 <- (post_sample$W1)[,j]
    W2 <- (post_sample$W2)[,j]
    beta_basis_1_k <- db.list_1[[k1]]
    beta_basis_2_k <- db.list_2[[k2]]
    
    p <- qgamma(1 - cumsum(E[-(L+1)])/sum(E), 
                shape = M/L, rate = 1)
    unip <- p/sum(p)
    
    selector1 <- findInterval(W1, (1:k1)/k1, left.open = T) + 1 
    selector2 <- findInterval(W2, (1:k2)/k2, left.open = T) + 1 
    
    post_funval[ , , j] <- tau * t(crossprod(unip * beta_basis_2_k[selector2, ],  beta_basis_1_k[selector1, ]))

  }
  
  post_logfunval <- log(post_funval)

  funval_mean <- apply(post_funval, c(1,2), mean)
  
  logfunval_median <- apply(post_logfunval, c(1,2), median)
  
  pointwise_90 <- apply(post_funval, c(1,2), quantile, probs = c(0.05, 0.95))
  
  if(unif_CR){
    
    logfunval_mad <- apply(post_logfunval, c(1,2), mad)
  
    logfunval_CI_help <- apply(abs(post_logfunval - array(logfunval_median, dim = dim(post_funval)))/array(logfunval_mad, dim = dim(post_funval)),
                             3, max, na.rm = TRUE)
  
    logfunval_Cvalue <- quantile(logfunval_CI_help, 0.9)

    return(list(tvpsd.mean = funval_mean, 
                tvpsd.median = exp(logfunval_median), 
                tvpsd.p05 = pointwise_90[1, , ], 
                tvpsd.p95 = pointwise_90[2, , ],
                tvpsd.u05 = exp(logfunval_median - logfunval_Cvalue * logfunval_mad),
                tvpsd.u95 = exp(logfunval_median + logfunval_Cvalue * logfunval_mad)))
    
    
  } else{
    
    return(list(tvpsd.mean = funval_mean, 
                tvpsd.median = exp(logfunval_median), 
                tvpsd.p05 = pointwise_90[1, , ], 
                tvpsd.p95 = pointwise_90[2, , ],
                tvpsd.u05 = NA,
                tvpsd.u95 = NA))
    
  }
  
  

}


#' Estimating the Bayes factor of hypothesis "k1 = 1".
#'
#' @param post_sample the posterior sample generated by \link{bdp_dw_mcmc}.
#' @param precision a positive integer specifying the number of terms used in approximating
#' the normalizing constant of the prior probability mass function of k1. Default 1000.
#' @return The Savage-Dickey estimate of the Bayes factor and its theoretical upper bound. c.f. section 3.3 of Tang et al. (2023).
#' @references Tang et al. (2023)
#' \emph{Bayesian nonparametric spectral analysis of locally stationary processes}
#' ArXiv preprint
#' <arXiv:2303.11561>
#' @export
#'
bdp_dw_bayes_factor_k1 <- function(post_sample, precision = 1000){

  prior_params <- post_sample$prior_params

  k1.theta <- prior_params$k1.theta
  
  if(k1.theta == 0.01){prior_norm_const <- 27.280819186230847} else{
    
    prior_norm_const <- sum(exp(-k1.theta * (1:precision) * log(1:precision))) #  prior_prob <- 1/prior_norm_const
    
  }
  
  post_prob <- mean(post_sample$k1 == 1)

  bf <- post_prob * prior_norm_const

  return(matrix(c(bf, prior_norm_const), ncol = 1, dimnames = list(c("Bayes factor", "upper bound"), c(""))))


}

#' time-varying spectral density function of the tvARMA(1,2) processes for illustrations
#' @param rescaled_time,freq numeric vectors forming a rectangular grid on which the tv-PSD is evaluated.
#' @param dgp optional: the tv-ARMA models demonstrated in section 4.2 of Tang et al. (2023). 
#' Should be chosen from "LS1", "LS2" and "LS3". See section Details.
#' @param a1,b1,b2 If dgp is not supplied, these arguments can be used to specify customized tv-ARMA
#' process (up to order(1,2)). See Details.
#' rescaled_time must be in \eqn{[0,1]} and freq must be in \eqn{[0,\pi]}.
#' @return a matrix of dimension length(rescaled_time) by length(freq).
#' @details See \link{sim_tvarma12} for the precise definition of a tvARMA(1,2) process. The time-varying
#' spectral density function of this process is defined as
#' 
#' \ifelse{html}{
#' \out{<math>f(u,&lambda;) = (2&pi;)<sup>-1</sup>(1+b<sub>1</sub><sup>2</sup>(u)+b<sub>2</sub><sup>2</sup>(u)+2b<sub>1</sub>(u)(b<sub>2</sub>(u)+1)cos(&lambda;)+2b<sub>2</sub>(u)cos(2&lambda;))/(1+a<sub>1</sub><sup>2</sup>(u)-2a<sub>1</sub>(u)cos(&lambda;)), (u,&lambda;)&isin;[0,1]&times;[0,&pi;],
#' </math>}
#' }{\deqn{f(u,\lambda) = \frac{1}{2\pi}\frac{1 + b_1^2(u) + b_2^2(u) + 2b_1(u)(b_2(u)+1)\cos(\lambda) + 2b_2(u)\cos(2\lambda)}{1+a_1^2(u)-2a_1(u)\cos(\lambda)},\quad (u,\lambda)\in[0,1]\times[0,\pi],}}
#' where \eqn{u} is called rescaled time and \eqn{\lambda} is called frequency.
#' 
#' For dgp = "LS1", it is a tvMA(2) process (MA order is 2) with
#' 
#' \ifelse{html}{
#' \out{<math>a<sub>1</sub>(u)=0, b<sub>1</sub>(u)=1.122(1-1.178sin(&pi;/2 u)), 
#' b<sub>2</sub>(u)=-0.81. </math>}
#' }{\deqn{a_1(u)=0, b_1(u)= 1.122(1 - 1.178\sin(\pi/2 u)), b_2(u) = -0.81.}}
#' For dgp = "LS2", it is a tvMA(1) process (MA order is 1) with
#' 
#' \ifelse{html}{
#' \out{<math>a<sub>1</sub>(u)=0, b<sub>1</sub>(u)=1.1cos(1.5-cos(4&pi; u)), b<sub>2</sub>(u)=0. 
#' </math>}
#' }{\deqn{a_1(u)=0, b_1(u)= 1.1\cos\left(1.5 - \cos\left(4\pi u \right) \right), b_2(u) = 0.}}
#' For dgp = "LS3", it is a tvAR(1) process (MA order is 0) with
#' 
#' \ifelse{html}{
#' \out{<math>a<sub>1</sub>(u)=1.2u-0.6, b<sub>1</sub>(u)=0, b<sub>2</sub>(u)=0. </math>}
#' }{\deqn{a_1(u)=1.2u-0.6, b_1(u)= 0, b_2(u) = 0.}}
#' @examples 
#' \dontrun{
#' res_time <- seq(0, 1, by = 0.005); freq <- pi*seq(0, 1, by = 0.01)
#' true_tvPSD <- psd_tvarma12(rescaled_time = res_time, freq = freq, dgp = "LS2")
#' plot(true_tvPSD)
#' }
#' @references Tang et al. (2023)
#' \emph{Bayesian nonparametric spectral analysis of locally stationary processes}
#' ArXiv preprint
#' <arXiv:2303.11561>
#' @export
psd_tvarma12 <- function(rescaled_time, 
                         freq, 
                         dgp = NULL, 
                         a1 = function(u){rep(0, length(u))}, 
                         b1 = function(u){rep(0, length(u))}, 
                         b2 = function(u){rep(0, length(u))}){

  stopifnot(!is.null(rescaled_time)); stopifnot(all(rescaled_time >= 0 & rescaled_time <= 1))

  stopifnot(!is.null(freq)); stopifnot(all(freq >= 0 & freq <= pi))
  
  if(dgp == "LS1"){
    
    a1 <- function(u1){rep(0, length(u1))}
    b1 <- function(u1){1.122 * (1 - 1.718 * sin(pi/2 * u1))}
    b2 <- function(u1){rep(-0.81, length(u1))}
    
  } else if(dgp == "LS2"){
    
    a1 <- function(u1){rep(0, length(u1))}
    b1 <- function(u1){1.1*cos(1.5 - cos(4*pi*u1))}
    b2 <- function(u1){rep(0, length(u1))}
    
  } else if(dgp == "LS3"){
    
    a1 <- function(u1){1.2*u1 - 0.6}
    b1 <- function(u1){rep(0, length(u1))}
    b2 <- function(u1){rep(0, length(u1))}

  }
  
  stopifnot(is.function(a1)); stopifnot(is.function(b1)); stopifnot(is.function(b2))
  

  if(length(rescaled_time) == 1){
    
    ar_part <- 1 + a1(rescaled_time)^2 - 2 * a1(rescaled_time) * cos(freq)
    
    ma_part <- 1 + b1(rescaled_time)^2 + b2(rescaled_time)^2 + 2 * b1(rescaled_time) * (b2(rescaled_time) + 1) * cos(freq) + 2 * b2(rescaled_time) * cos(2*freq)
    
    psd <- 1/(2*pi) * ma_part/ar_part
    
  } else{
    
    ar_part <- 1 + a1(rescaled_time)^2 - 2 * a1(rescaled_time) %o% cos(freq)
    
    ma_part <- (1 + b1(rescaled_time)^2 + b2(rescaled_time)^2) + 2 * (b1(rescaled_time) * (b2(rescaled_time) + 1)) %o% cos(freq) + 2 * b2(rescaled_time) %o% cos(2*freq)
    
    psd <- 1/(2*pi) * ma_part/ar_part
    
  }
  
    result <- list(rescaled_time = rescaled_time,
                   frequency = freq,
                   tv_psd = psd)
    
    class(result) <- "bdp_dw_tv_psd"

    return(result)
    
}



#' simulate from the tvARMA(1,2) process for illustration
#' @importFrom stats rnorm rt rexp
#' @param len_d a positive integer indicating the length of the simulated process.
#' @param dgp optional: the tv-ARMA models demonstrated in section 4.2 of Tang et al. (2023). Should be chosen from "LS1", "LS2" and "LS3". See section Details.
#' @param ar_order,ma_order,a1,b1,b2 If dgp is not supplied, these arguments can be used to specify customized tv-ARMA
#' process (up to order(1,2)). See details.
#' @param innov_distribution optional: the distributions of innovation used in section 4.2.2 of Tang et al. (2023) . 
#' Should be chosen from "a", "b", "c". "a" denotes standard normal distribution, 
#' "b" indicates standardized Student-t distribution with degrees of freedom 4 and
#' "c" denotes standardized Pareto distribution with scale 1 and shape 4.
#' @param wn If innov_distribution is not specified, one may supply its own innovation sequence. Please make sure the length
#' of wn is at least the sum of len_d and the MA order of the process. If ma_order is specified, then MA order is exactly
#' ma_order. If dgp is specified, the MA order of "LS1", "LS2" and "LS3" can be found in section Details below.
#' @return a numeric vector of length len_d simulated from the given process.
#' @details This function simulates from the following time-varying Autoregressive Moving Average model with order (1,2):
#' 
#' \ifelse{html}{\out{<math>X<sub>t,T</sub> = a<sub>1</sub>(<sup>t</sup>&frasl;<sub>T</sub>) X<sub>t-1,T</sub> + w<sub>t</sub> + b<sub>1</sub>(<sup>t</sup>&frasl;<sub>T</sub>) w<sub>t-1</sub> + b<sub>2</sub>(<sup>t</sup>&frasl;<sub>T</sub>) w<sub>t-2</sub>, t=1,2,...,T,</math>}
#' }{\deqn{X_{t,T} = a_1(t/T)X_{t-1,T} + w_{t} + b_1(t/T) w_{t-1} + b_2(t/T) w_{t-2}, \quad t=1,2,\cdots,T,}}
#' where \eqn{T} is the length specified and \ifelse{html}{\out{{w<sub>t</sub>}}}{\eqn{\{w_t\}}} are 
#' a sequence of i.i.d. random variables with mean 0 and standard deviation 1.
#' 
#' For dgp = "LS1", it is a tvMA(2) process (MA order is 2) with
#' 
#' \ifelse{html}{
#' \out{<math>a<sub>1</sub>(u)=0, b<sub>1</sub>(u)=1.122(1-1.178sin(&pi;/2 u)), b<sub>2</sub>(u)=-0.81.
#' </math>}
#' }{\deqn{a_1(u)=0, b_1(u)= 1.122(1 - 1.178\sin(\pi/2 u)), b_2(u) = -0.81.}}
#' For dgp = "LS2", it is a tvMA(1) process (MA order is 1) with
#' 
#' \ifelse{html}{
#' \out{<math>a<sub>1</sub>(u)=0, b<sub>1</sub>(u)=1.1cos(1.5-cos(4&pi; u)), b<sub>2</sub>(u)=0.
#' </math>}
#' }{\deqn{a_1(u)=0, b_1(u)= 1.1\cos\left(1.5 - \cos\left(4\pi u \right) \right), b_2(u) = 0.}}
#' For dgp = "LS3", it is a tvAR(1) process (MA order is 0) with
#' 
#' \ifelse{html}{
#' \out{<math>a<sub>1</sub>(u)=1.2u-0.6, b<sub>1</sub>(u)=0, b<sub>2</sub>(u)=0.</math>}
#' }{\deqn{a_1(u)=1.2u-0.6, b_1(u)= 0, b_2(u) = 0.}}
#' 
#' @examples 
#' \dontrun{
#' sim_tvarma12(len_d = 1500, 
#' dgp = "LS2", 
#' innov_distribution = "a") # generate from LS2a
#' 
#' sim_tvarma12(len_d = 1500, 
#' dgp = "LS2", 
#' wn = rnorm(1502)) # again, generate from LS2a
#' 
#' sim_tvarma12(len_d = 1500, 
#' ar_order = 0, 
#' ma_order = 1, 
#' b1 = function(u){1.1*cos(1.5 - cos(4*pi*u))}, 
#' innov_distribution = "a") # again, generate from LS2a
#' 
#' sim_tvarma12(len_d = 1500, 
#' ar_order = 0, 
#' ma_order = 1, 
#' b1 = function(u){1.1*cos(1.5 - cos(4*pi*u))}, 
#' wn = rnorm(1502)) # again, generate from LS2a
#' }
#' @references Tang et al. (2023)
#' \emph{Bayesian nonparametric spectral analysis of locally stationary processes}
#' ArXiv preprint
#' <arXiv:2303.11561>
#' @export
sim_tvarma12 <- function(len_d, 
                         dgp = NULL, 
                         ar_order = 1, 
                         ma_order = 2, 
                         a1 = NULL, 
                         b1 = NULL, 
                         b2 = NULL, 
                         innov_distribution = NULL, 
                         wn = NULL){
  
  len_d <- as.integer(len_d)

  stopifnot(!is.null(len_d)); stopifnot(length(len_d) == 1 & len_d >= 3)
  
  if(is.null(dgp)){
    
    stopifnot(!is.null(ar_order) & !is.null(ma_order)); stopifnot(ar_order %in% 0:1); stopifnot(ma_order %in% 0:2)
    
    if(ar_order == 1){stopifnot(!is.null(a1)); stopifnot(is.function(a1))}
    if(ma_order == 1){stopifnot(!is.null(b1)); stopifnot(is.function(b1))}
    if(ma_order == 2){stopifnot(!is.null(b1)); stopifnot(is.function(b1)); stopifnot(!is.null(b2)); stopifnot(is.function(b2))}
    
    dgp <- "NA"
    
  } else {stopifnot(dgp %in% c("LS1", "LS2", "LS3"))}
  
  if(dgp == "LS1"){
    
    ar_order <- 0
    ma_order <- 2
    
    b1 <- function(u1){1.122 * (1 - 1.718 * sin(pi/2 * u1))}
    b2 <- function(u1){rep(-0.81, length(u1))}
    
  } else if(dgp == "LS2"){
    
    ar_order <- 0
    ma_order <- 1
    
    b1 <- function(u1){1.1*cos(1.5 - cos(4*pi*u1))}
    
  } else if(dgp == "LS3"){
    
    ar_order <- 1
    ma_order <- 0
    
    a1 <- function(u1){1.2*u1 - 0.6}
    
  }
  

  sim_data <- integer(len_d)
  
  stopifnot(!is.null(innov_distribution) | !is.null(wn))
  
  if(!is.null(innov_distribution)){
    
    stopifnot(innov_distribution %in% c("a", "b", "c"))
    
    if(innov_distribution == "a"){wn <- rnorm(len_d + ma_order, 0, 1)}
    if(innov_distribution == "b"){wn <- rt(len_d + ma_order, df = 3)/sqrt(3/(3-2))}
    if(innov_distribution == "c"){wn <- (1*exp(rexp(len_d + ma_order, 4)) - 4/3)/sqrt(2/9)}
    
  } else if(!is.null(wn)){
    
    stopifnot(is.numeric(wn) & is.vector(wn))
    
    stopifnot(length(wn) >= len_d + ma_order)
    
    wn <- (wn - mean(wn))/sd(wn)
    
  }

  
  
  
  if(ar_order == 0){
    
    if(ma_order == 0){sim_data <- wn[1:len_d]}
    
    if(ma_order == 1){
      
       for (t in 1:len_d) {
         sim_data[t] <- wn[t+1] +  b1(t/len_d) * wn[t]
    }
      
    }
    
    if(ma_order == 2){
      
      for (t in 1:len_d) {
        sim_data[t] <- wn[t+2] + b1(t/len_d) * wn[t+1] + b2(t/len_d) * wn[t] 
      }
      
    }
    
    
  } else{# ar_order == 1
    
    
    if(ma_order == 0){
      
      sim_data[1] <- rnorm(1)
      
      for (t in 2:len_d) {
        sim_data[t] <- a1(t/len_d) * sim_data[t-1] + wn[t]
      }
      
      
    }
    
    if(ma_order == 1){
      
      sim_data[1] <- rnorm(1)
      
      for (t in 2:len_d) {
        sim_data[t] <- a1(t/len_d) * sim_data[t-1] + wn[t+1] +  b1(t/len_d) * wn[t]
      }
      
    }
    
    if(ma_order == 2){
      
      sim_data[1] <- rnorm(1)
      
      for (t in 2:len_d) {
        sim_data[t] <- a1(t/len_d) * sim_data[t-1] + wn[t+2] + b1(t/len_d) * wn[t+1] + b2(t/len_d) * wn[t] 
      }
      
    }
    
  }
  
  

  return(sim_data)

}



#' BDP-DW method: performing posterior sampling and calculating statistics based on the posterior samples
#' 
#' @param data time series that needs to be analyzed
#' @param m window size needed to calculate moving periodogram.
#' @param likelihood_thinning the thinning factor of the dynamic Whittle likelihood.
#' @param monitor a Boolean value (default TRUE) indicating whether to display the real-time status
#' @param print_interval If monitor = TRUE, then this value indicates number of iterations after which a status is printed to console;
#' If monitor = FALSE, it does not have any effect
#' @param unif_CR a Boolean value (default FALSE) indicating whether to calculate the uniform credible region
#' @param rescaled_time,freq a set of grid lines in \eqn{[0,1]} and \eqn{[0,\pi]}, respectively, specifying where to evaluate the estimated tv-PSD
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (for post-processing of the posterior sample)
#' @param adaptive.batchSize the batch size for the adaptive MCMC algorithm for sampling tau
#' @param adaptive.targetAcceptanceRate the target acceptance rate for the adaptive MCMC algorithm for sampling tau
#' @param M DP base measure constant (> 0)
#' @param g0.alpha,g0.beta parameters of Beta base measure of DP
#' @param k1.theta prior parameter for polynomial corresponding to rescaled time (propto exp(-k1.theta*k1*log(k1)))
#' @param k2.theta prior parameter for polynomial corresponding to rescaled frequency (propto exp(-k2.theta*k2*log(k2)))
#' @param tau.alpha,tau.beta prior parameters for tau (inverse gamma)
#' @param k1max upper bound of the degrees of Bernstein polynomial
#' corresponding to rescaled time (for pre-computation of basis functions)
#' @param k2max upper bound of the degrees of Bernstein polynomial
#' corresponding to rescaled frequency (for pre-computation of basis functions)
#' @param L truncation parameter of DP in stick breaking representation
#' @param bernstein1_l,bernstein1_r left and right truncation of Bernstein polynomial basis functions
#' for rescaled time, 0<=bernstein1_l<bernstein1_r<=1
#' @param bernstein2_l,bernstein2_r left and right truncation of Bernstein polynomial basis functions
#' for rescaled frequency, 0<=bernstein2_l<bernstein2_r<=1
#' @return list containing the following fields:
#' 
#'    \item{k1,k2,tau,V,W1,W2}{posterior traces of PSD parameters}
#'    \item{lpost}{traces log posterior}
#'    \item{tim}{total run time}
#'    \item{bf_k1}{Savage-Dickey estimate of Bayes factor of hypothesis k1=1}
#'    \item{tvpsd.mean,tvpsd.median}{posterior mean and pointwise posterior median (matrices of dimension length(rescaled_time) by length(freq))}
#'    \item{tvpsd.p05,tvpsd.p95}{90 percent pointwise credibility interval}
#'    \item{tvpsd.u05,tvpsd.u95}{90 percent uniform credibility interval if unif_CR = TRUE. Otherwise NA}
#' @examples
#' \dontrun{
#'
#' ##
#' ## Example: Applying BDP-DW method to a multi-peaked heavy-tailed tvMA(1) process
#' ##
#'
#' # set seed
#' set.seed(2)
#' # set the length of time series
#' len_d <- 1500
#' # generate data from DGP LS2b defined in Section 4.2.2 of Tang et al. (2023). 
#' # see also ?sim_tvarma12
#' sim_data <- sim_tvarma12(len_d = 1500, dgp = "LS2", innov_distribution = "b")
#' set.seed(NULL)
#' # specify grid-points at which the tv-PSD is evaluated
#' res_time <- seq(0, 1, by = 0.005); freq <- pi * seq(0, 1, by = 0.01)
#' # calculate the true tv-PSD of DGP LS2b at the pre-specified grid
#' true_tvPSD <- psd_tvarma12(rescaled_time = res_time, freq = freq, dgp = "LS2")
#' # plot the true tv-PSD
#' # type ?plot.bdp_dw_tv_psd for more info
#' plot(true_tvPSD)
#'
#' # If you run the example be aware that this may take several minutes
#' print("This example may take some time to run")
#' result <- gibbs_bdp_dw(data = sim_data, 
#' m = 50, 
#' likelihood_thinning = 2, 
#' rescaled_time = res_time,
#'  freq = freq)
#'
#' # extract bayes factor and examine posterior summary
#' bayes_factor(result)
#' summary(result)
#' 
#' # compare estimate with true function
#' # type ?plot.bdp_dw_result for more info
#' par(mfrow = c(1,2))
#' 
#' plot(result, which = 1,
#' zlim = range(result$tvpsd.mean, true_tvPSD$tv_psd)
#' )
#' plot(true_tvPSD,
#' zlim = range(result$tvpsd.mean, true_tvPSD$tv_psd),
#' main = "true tv-PSD")
#' 
#' par(mfrow = c(1,1))
#'
#' }
#' @references Tang et al. (2023)
#' \emph{Bayesian nonparametric spectral analysis of locally stationary processes}
#' ArXiv preprint
#' <arXiv:2303.11561>
#' @export
gibbs_bdp_dw <- function(data,
                         m,
                         likelihood_thinning = 1,
                         monitor = TRUE,
                         print_interval = 100,
                         unif_CR = FALSE,
                         rescaled_time,
                         freq,
                         Ntotal = 1.1e5,
                         burnin = 6e4,
                         thin = 10,
                         adaptive.batchSize = 50,
                         adaptive.targetAcceptanceRate = 0.44,
                         M = 1,
                         g0.alpha = 1,
                         g0.beta = 1,
                         k1.theta = 0.01,
                         k2.theta = 0.01,
                         tau.alpha = 0.001,
                         tau.beta = 0.001,
                         k1max = 100,
                         k2max = 100,
                         L = 10,
                         bernstein1_l = 0.1,
                         bernstein1_r = 0.9,
                         bernstein2_l = 0.1,
                         bernstein2_r = 0.9){
  
  
  mcmc_params <- bdp_dw_mcmc_params_gen(Ntotal,
                                    burnin,
                                    thin,
                                    adaptive.batchSize,
                                    adaptive.targetAcceptanceRate)
  
  
  prior_params <- bdp_dw_prior_params_gen(M,
                                      g0.alpha,
                                      g0.beta,
                                      k1.theta,
                                      k2.theta,
                                      tau.alpha,
                                      tau.beta,
                                      k1max,
                                      k2max,
                                      L,
                                      bernstein1_l,
                                      bernstein1_r,
                                      bernstein2_l,
                                      bernstein2_r)
  
  cl <- match.call()
  
  post_sample <- bdp_dw_mcmc(data, m, likelihood_thinning, mcmc_params,
                         prior_params, monitor, print_interval)
  
  
  bf_k1 <- bdp_dw_bayes_factor_k1(post_sample)
  
  
  post_funval_stats <- bdp_dw_est_post_stats(post_sample = post_sample,
                                             rescaled_time = rescaled_time, 
                                             freq = freq, unif_CR = unif_CR)
  
  
  
  result <- list(call = cl, 
                 data = post_sample$data,
                 rescaled_time = rescaled_time,
                 frequency = freq,
                 k1 = post_sample$k1, 
                 k2 = post_sample$k2, 
                 tau = post_sample$tau,
                 E = post_sample$E, 
                 W1 = post_sample$W1, 
                 W2 = post_sample$W2, 
                 lpost = post_sample$lpost,
                 tim = post_sample$tim, 
                 bayes_factor_k1 = bf_k1, 
                 tvpsd.mean = post_funval_stats$tvpsd.mean,
                 tvpsd.median = post_funval_stats$tvpsd.median, 
                 tvpsd.p05 = post_funval_stats$tvpsd.p05,
                 tvpsd.p95 = post_funval_stats$tvpsd.p95,
                 unif_CR = unif_CR,
                 tvpsd.u05 = post_funval_stats$tvpsd.u05,
                 tvpsd.u95 = post_funval_stats$tvpsd.u95)
  
  class(result) <- "bdp_dw_result"
  
  return(result)
}




