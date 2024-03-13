#' Plot method for gibbs_psd class
#' @details Visualizes the spectral density estimate (pointwise posterior median), along with the periodogram
#' and credibility regions.
#' If the data has missing values, the periodogram is computed with a linearly
#' interpolated version of the data using \link[forecast]{na.interp}.
#' @method plot gibbs_psd
#' @param x an object of class gibbs_psd
#' @param pdgrm bool flag indicating whether periodogram is visualized or not
#' @param credib string indicating which credible regions are visualized. Possible values are "pointwise", "uniform", "both" and "none".
#' @param log logical value to determine if the individual spectra are visualized on a log scale
#' @param ... further arguments to be parsed to \code{plot.default}
#' @export
plot.gibbs_psd <- function(x, pdgrm=T, credib="both", log=T, ...) {
  mcmc <- x
  f_star <- mcmc$psd.median
  g_list <- list()
  lty <- col <- 1
  lwd <- 1.5
  # periodogram
  if (pdgrm) {
    if (!is.matrix(mcmc$data)) {
      # univariate
      if (any(is.na(mcmc$data))) {
        data <- as.numeric(forecast::na.interp(mcmc$data))
      } else {
        data <- mcmc$data
      }
      pdgrm <- abs(fast_ft(data, real=F))^2 / 2 / pi
    } else {
      # multivariate
      if (any(is.na(mcmc$data))) {
        data <- apply(mcmc$data, 2, 
                      function(z) as.numeric(forecast::na.interp(z)))
      } else {
        data <- mcmc$data
      }
      pdgrm <- mpdgrm(mdft(data))
    }
    g_list[[length(g_list)+1]] <- pdgrm
    lty <- c(lty, 1)
    col <- c(col, "gray")
    lwd <- c(lwd, 1)
  } 
  # pointwise intervals
  if (credib == "pointwise" || credib == "both") {
    g_list[[length(g_list)+1]] <- mcmc$psd.p05
    g_list[[length(g_list)+1]] <- mcmc$psd.p95
    lty <- c(lty, 2, 2)
    col <- c(col, 4, 4)
    lwd <- c(lwd, 1.5, 1.5)
  } 
  # uniform intervals
  if (credib == "uniform" || credib == "both") {
    g_list[[length(g_list)+1]] <- mcmc$psd.u05
    g_list[[length(g_list)+1]] <- mcmc$psd.u95
    lty <- c(lty, 3, 3)
    col <- c(col, 2, 2)
    lwd <- c(lwd, 1.5, 1.5)
  } 
  plotMPsd(f_star, g=g_list, lty=lty, col=col, lwd=lwd, log=log)
}

#' Convert psd vector to array
#' (compatibility: to use plotMPsd for univariate functions as well)
#' @keywords internal
psd_array <- function(f) {
  stopifnot(is.vector(f))
  return(array(f, dim=c(1,1,length(f))))
}

#' Visualization of multivariate PSDs
#' Used in \code{plot.gibbs_psd}
#' @importFrom graphics par
#' @keywords internal
plotMPsd <- function(f, # main psd to plot
                     g=NULL, # list of additional psds (CI's, grount truth, etc)
                     lty=rep(1, 1+length(g)), # line type
                     col=rep(1, 1+length(g)), # color
                     type=rep("l", 1+length(g)), # line or points (e.g. for pdgrm)?
                     pch=rep("0", 1+length(g)), # point type for points
                     lwd=rep(1, 1+length(g)),
                     log=F, # plot diagonals on log scale?
                     ylim.compound=T, # ylim according to c(f,g) instead of only g?
                     mar=c(4,4,4,3), 
                     mfrow=c(dim(f)[1],dim(f)[1]),
                     lambda_scaling=pi,
                     ylab_prefix="f",
                     xlab="Frequency",
                     excludeBoundary=T,
                     ...) {
  stopifnot(length(lty) == 1+length(g))
  stopifnot(length(col) == 1+length(g))
  stopifnot(length(type) == 1+length(g))
  if (!is.array(f)) {
    # dummy 1x1xN representation in univariate case
    f <- psd_array(f)
    g <- lapply(g, psd_array)
  }
  d <- dim(f)[1]
  N <- dim(f)[3]
  if (excludeBoundary) {
    lambda <- (1:(N-2))/(N-1)*lambda_scaling
    lambda_ind <- 2:(N-1)
  } else {
    lambda <- (0:(N-1))/(N-1)*lambda_scaling
    lambda_ind <- 1:N
  }
  if (d>1 && any(is.complex(f))) {
    f_plot <- realValuedPsd(f)
  } else {
    f_plot <- f
  }
  if (d>1 && any(sapply(g, is.complex))) {
    stopifnot(all(sapply(g, is.complex)))
    g_plot <- lapply(g, realValuedPsd)
  } else {
    g_plot <- g
  }
  par(mfrow=mfrow, mar=mar)
  for (i in 1:d) {
    for (j in 1:d) {
      ylab <- ylab_prefix
      if (d>1) ylab <- paste0(ylab, "_", i, j)
      if (i==j && log) {
        if (log) {
          ylab <- paste0("log(", ylab, ")")
        }
        if (ylim.compound) {
          ylim_min <- min(min(log(f_plot[i,j,lambda_ind])), 
                          as.numeric(sapply(g_plot, function(z) { 
                            if (is.null(z)) return(Inf) 
                            else return(min(log(z[i,j,lambda_ind])))})))
          ylim_max <- max(max(log(f_plot[i,j,lambda_ind])), 
                          as.numeric(sapply(g_plot, function(z) { 
                            if (is.null(z)) return(-Inf) 
                            else return(max(log(z[i,j,lambda_ind])))})))
          ylim=c(ylim_min, ylim_max)
        } else {
          ylim <- NULL
        }
        plot(x=lambda, 
             y=log(f_plot[i,j,lambda_ind]), 
             col=col[1], 
             lty=lty[1], 
             ylim=ylim, 
             type=type[1], 
             pch=pch[1],
             lwd=lwd[1],
             xlab=xlab, 
             ylab=ylab, 
             ...)
        for (ll in seq_len(length(g))) {
          lines(x=lambda, y=log((g_plot[[ll]])[i,j,lambda_ind]), col=col[1+ll], 
                lty=lty[1+ll], type=type[1+ll], pch=pch[1+ll], lwd=lwd[1+ll], ...)
        }
      } else {
        if (i<j) {
          ylab <- paste0("Re(", ylab, ")")
        } else {
          if (i==j) {
            ylab <- ylab
          } else {
            ylab <- paste0("Im(", ylab, ")")
          }
        }
        if (ylim.compound) {
          ylim_min <- min(min((f_plot[i,j,lambda_ind])), 
                          as.numeric(sapply(g_plot, function(z) { 
                            if (is.null(z)) return(Inf) 
                            else return(min((z[i,j,lambda_ind])))})))
          ylim_max <- max(max((f_plot[i,j,lambda_ind])), 
                          as.numeric(sapply(g_plot, function(z) { 
                            if (is.null(z)) return(-Inf) 
                            else return(max((z[i,j,lambda_ind])))})))
          ylim=c(ylim_min, ylim_max)
        } else {
          ylim <- NULL
        }
        plot(x=lambda, 
             y=f_plot[i,j,lambda_ind], 
             col=col[1], 
             lty=lty[1], 
             ylim=ylim, 
             type=type[1], 
             pch=pch[1], 
             lwd=lwd[1],
             xlab=xlab,
             ylab=ylab, 
             ...)
        for (ll in seq_len(length(g))) {
          lines(x=lambda, y=(g_plot[[ll]])[i,j,lambda_ind], col=col[1+ll], 
                lty=lty[1+ll], type=type[1+ll], pch=pch[1+ll], lwd=lwd[1+ll], ...)
        }
      }
    }
  }
  par(mfcol=c(1,1))
}

#' Helping function for print and summary (both are quite similar)
#' @param flag An indicator. flag=T indicates print, flag=F indicates summary
#' @keywords internal
print_summary_gibbs_psd_help <- function(mcmc, flag=T) {
  # title
  if (mcmc$algo=="gibbs_ar") {
    hstr <- "Bayesian autoregressive model with PACF parametrization"
    d <- 1
    n <- length(mcmc$data)
    N <- length(mcmc$sigma2)
    p <- dim(mcmc$rho)[1]
    if (flag) {
      print_estimates <- function() {
        if (p>0) {
          cat("Posterior median rho:\n")
          print(apply(mcmc$rho,1,median))        
        }
        cat("\nPosterior median sigma2:\n")
        print(median(mcmc$sigma2))
      }
    } else {
      print_estimates <- function() {
        if (p>0) {
          cat("Summary rho:\n")
          print(apply(mcmc$rho,1,summary))        
        }
        cat("\nSummary sigma2:\n")
        print(summary(mcmc$sigma2))
      }
    }
  }
  else if (mcmc$algo=="gibbs_np") {
    hstr <- "Bayesian nonparametric inference with Whittle likelihood"
    d <- 1
    n <- length(mcmc$data)
    N <- length(mcmc$tau)
    if (flag) {
      print_estimates <- function() {
        cat("Posterior median k:\n")
        print(median(mcmc$k))
        cat("\nPosterior median tau:\n")
        print(median(mcmc$tau))
        cat("\nPosterior median V:\n")
        print(apply(mcmc$V,1,median))  
        cat("\nPosterior median W:\n")
        print(apply(mcmc$W,1,median)) 
      }
    } else {
      print_estimates <- function() {
        cat("Summary k:\n")
        print(summary(mcmc$k))
        cat("\nSummary tau:\n")
        print(summary(mcmc$tau))
        cat("\nSummary V:\n")
        print(apply(mcmc$V,1,summary))  
        cat("\nSummary W:\n")
        print(apply(mcmc$W,1,summary)) 
      }
    }
  }
  else if (mcmc$algo=="gibbs_npc") {
    hstr <- "Bayesian semiparametric inference with corrected AR likelihood"
    d <- 1
    n <- length(mcmc$data)
    N <- length(mcmc$tau)
    if (flag) {
      print_estimates <- function() {
        cat("Posterior median k:\n")
        print(median(mcmc$k))
        cat("\nPosterior median tau:\n")
        print(median(mcmc$tau))
        cat("\nPosterior median V:\n")
        print(apply(mcmc$V,1,median))  
        cat("\nPosterior median W:\n")
        print(apply(mcmc$W,1,median)) 
        cat("\nPosterior median rho:\n")
        print(apply(mcmc$rho,1,median)) 
        cat("\nPosterior median eta:\n")
        print(median(mcmc$eta))
      }
    } else {
      print_estimates <- function() {
        cat("Summary k:\n")
        print(summary(mcmc$k))
        cat("\nSummary tau:\n")
        print(summary(mcmc$tau))
        cat("\nSummary V:\n")
        print(apply(mcmc$V,1,summary))  
        cat("\nSummary W:\n")
        print(apply(mcmc$W,1,summary)) 
        cat("\nSummary rho:\n")
        print(apply(mcmc$rho,1,summary)) 
        cat("\nSummary eta:\n")
        print(summary(mcmc$eta))
      }
    }
  }
  else if (mcmc$algo=="gibbs_var") {
    hstr <- "Bayesian vector autoregressive model"
    d <- ncol(mcmc$data)
    n <- nrow(mcmc$data)
    N <- dim(mcmc$beta)[2]
    p <- dim(mcmc$beta)[1]
    if (flag) {
      print_estimates <- function() {
        if (p>0) {
          cat("Posterior median beta:\n")
          print(apply(mcmc$beta,1,median))        
        }
        cat("\nPosterior median Sigma:\n")
        print(apply(mcmc$Sigma, c(1,2), median))
      }
    } else {
      print_estimates <- function() {
        if (p>0) {
          cat("Summary beta:\n")
          print(apply(mcmc$beta,1,summary))        
        }
        cat("\nSummary Sigma:\n")
        print(apply(mcmc$Sigma, c(1,2), summary))
      }
    }
  }
  else if (mcmc$algo=="gibbs_vnp") {
    hstr <- "Bayesian nonparametric multivariate inference with Whittle likelihood"
    d <- ncol(mcmc$data)
    n <- nrow(mcmc$data)
    N <- length(mcmc$k)
    if (flag) {
      print_estimates <- function() {
        cat("Posterior median k:\n")
        print(median(mcmc$k))
        cat("\nPosterior median r:\n")
        print(apply(mcmc$r, 1, median))
        cat("\nPosterior median x:\n")
        print(apply(mcmc$x, 1,median))  
        #cat("\nPosterior median U:\n")
        #print(apply(mcmc$U,c(1,2,3),median)) 
      }
    } else {
      print_estimates <- function() {
        cat("Summary k:\n")
        print(summary(mcmc$k))
        cat("\nSummary r:\n")
        print(apply(mcmc$r, 1, summary))
        cat("\nSummary x:\n")
        print(apply(mcmc$x, 1,summary))  
        #cat("\nPosterior median U:\n")
        #print(apply(mcmc$U,c(1,2,3),summary)) 
      }
    }
  }
  else stop("Unknown algorithm.")
  cat("\n", hstr, "\n", sep="")
  # call
  cat("\nCall:\n"); print(mcmc$call); cat("\n")
  print_estimates()
  if (any(is.na(mcmc$data))) {
    if (flag) {
      cat("\nPostrior median missing values:\n")
      print(apply(mcmc$missing_values, 2, median))      
    } else {
      cat("\nSummary missing values:\n")
      print(apply(mcmc$missing_values, 2, summary)) 
    }
  }
  cat("\nPosterior sample size: N=", N, "\n", sep="")
  cat("Number of observations: n=", n, "\n", sep="")
  cat("Dimension: d=", d, "\n", sep="")
}

#' Print method for gibbs_psd class
#' @param x object of class gibbs_psd
#' @param ... not in use
#' @export
print.gibbs_psd <- function(x, ...) {
  print_summary_gibbs_psd_help(x, T)
}

#' Summary method for gibbs_psd class
#' @param object object of class gibbs_psd
#' @param ... not in use
#' @export
summary.gibbs_psd <- function(object, ...) {
  print_summary_gibbs_psd_help(object, F)
}