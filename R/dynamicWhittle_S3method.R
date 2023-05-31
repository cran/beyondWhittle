
#' Plot method for bdp_dw_tv_psd class
#' @details Visualizes the spectral density function of time-varying 
#' @method plot bdp_dw_tv_psd
#' @param x an object of class bdp_dw_tv_psd
#' @param col choice of color, default hcl.colors(200, "Blue-Red 3").
#' @param ... further arguments to be parsed to \code{image.default}
#' @importFrom graphics image
#' @importFrom grDevices hcl.colors
#' @export
plot.bdp_dw_tv_psd <- function(x, col = hcl.colors(200, "Blue-Red 3"), ...){
  image(x = x$rescaled_time,
        xlab = "rescaled_time",
        y = x$frequency,
        ylab = "frequency",
        z = x$tv_psd,
        col = col,
        ...)
}


#' Print method for bdp_dw_result class
#' @param x object of class bdp_dw_result
#' @param ... not in use
#' @export
print.bdp_dw_result <- function(x, ...) {
  # title
  hstr <- "Bayesian nonparametric inference of time-varying spectrum with dynamic Whittle likelihood"
  n <- length(x$data)
  N <- length(x$tau)
  tim <- x$tim
  bf <- x$bayes_factor_k1
  
  print_estimates <- function() {
      cat("Posterior median k1:\n")
      print(median(x$k1))
      cat("Posterior median k2:\n")
      print(median(x$k2))
      cat("\nPosterior median tau:\n")
      print(median(x$tau))
      cat("\nPosterior median V:\n")
      print(apply(x$V,1,median))  
      cat("\nPosterior median W1:\n")
      print(apply(x$W1,1,median)) 
      cat("\nPosterior median W2:\n")
      print(apply(x$W2,1,median))
  }
  
  cat("\n", hstr, "\n", sep="")
  # call
  cat("\nCall:\n"); print(x$call); cat("\n")
  print_estimates()
  cat("\nEstimated Bayes factor of {k1=1}: BF=", bf[1,1], ". The theoretical upper bound is ", bf[2,1], "\n", sep="")
  cat("\nPosterior sample size: N=", N, "\n", sep="")
  cat("Number of observations: n=", n, "\n", sep="")
  cat("Run-time: ", unclass(tim), " ", units(tim), sep="")
 
}
  


#' Summary method for bdp_dw_result class
#' @param object object of class bdp_dw_result
#' @param ... not in use
#' @export
summary.bdp_dw_result <- function(object, ...) {
  # title
  hstr <- "Bayesian nonparametric inference of time-varying spectrum with dynamic Whittle likelihood"
  n <- length(object$data)
  N <- length(object$tau)
  tim <- object$tim
  bf <- object$bayes_factor_k1
  
  print_estimates <- function() {
    cat("Summary k1:\n")
    print(summary(object$k1))
    cat("Summary k2:\n")
    print(summary(object$k2))
    cat("\nSummary tau:\n")
    print(summary(object$tau))
    cat("\nSummary V:\n")
    print(apply(object$V,1,summary))  
    cat("\nSummary W1:\n")
    print(apply(object$W1,1,summary))
    cat("\nSummary W2:\n")
    print(apply(object$W2,1,summary))
  }
  
  cat("\n", hstr, "\n", sep="")
  # call
  cat("\nCall:\n"); print(object$call); cat("\n")
  print_estimates()
  cat("\nEstimated Bayes factor of {k1=1}: BF=", bf[1,1], ". The theoretical upper bound is ", bf[2,1], "\n", sep="")
  cat("\nPosterior sample size: N=", N, "\n", sep="")
  cat("Number of observations: n=", n, "\n", sep="")
  cat("Run-time: ", unclass(tim), " ", units(tim), sep="")
  
}



#' Plot method for bdp_dw_result class
#' @param x object of class bdp_dw_result
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:6. 1 indicates posterior mean,
#' 2 indicates posterior median, 3 for lower bound of pointwise 90 percent credible region, 4 for upper bound of 
#' pointewise 90 percent credible region, 5 indicates lower bound of uniform 90 percent credible region, 6 indicates
#' upper bound of uniform 90 percent credible region.
#' @param ask logical; if TRUE, the user is asked before each plot.
#' @param col choice of color, default hcl.colors(200, "Blue-Red 3").
#' @param ... other parameters to be passed through to \code{image.default}
#' @importFrom grDevices hcl.colors dev.flush dev.hold dev.interactive devAskNewPage
#' @export
plot.bdp_dw_result <- function(x, which = 1:4,
                           ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                           col = hcl.colors(200, "Blue-Red 3"), ...) {
  
  if (!inherits(x, "bdp_dw_result"))
    stop("use only with \"bdp_dw_result\" objects")
  if(!is.numeric(which) || any(which < 1) || any(which > 6))
    stop("'which' must be in 1:6")
  
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  if (show[1L]) {
    dev.hold()
    
    image(x = x$rescaled_time,
          xlab = "rescaled_time",
          y = x$frequency,
          ylab = "frequency",
          z = x$tvpsd.mean,
          main = "posterior mean",
          col = col,
          ...)
    
    dev.flush()
  }
  
  if (show[2L]) {
    dev.hold()
    
    image(x = x$rescaled_time,
          xlab = "rescaled_time",
          y = x$frequency,
          ylab = "frequency",
          z = x$tvpsd.median,
          main = "posterior median",
          col = col,
          ...)
    
    dev.flush()
  }
  
  if (show[3L]) {
    dev.hold()
    
    image(x = x$rescaled_time,
          xlab = "rescaled_time",
          y = x$frequency,
          ylab = "frequency",
          z = x$tvpsd.p05,
          main = "90 percent pointwise CR: lower bound",
          col = col,
          ...)
    
    dev.flush()
  }
  
  if (show[4L]) {
    dev.hold()
    
    image(x = x$rescaled_time,
          xlab = "rescaled_time",
          y = x$frequency,
          ylab = "frequency",
          z = x$tvpsd.p95,
          main = "90 percent pointwise CR: upper bound",
          col = col,
          ...)
    
    dev.flush()
  }
  
  
  
  if ((x$unif_CR) & show[5L]) {
    dev.hold()
    
    image(x = x$rescaled_time,
          xlab = "rescaled_time",
          y = x$frequency,
          ylab = "frequency",
          z = x$tvpsd.u05,
          main = "90 percent uniform CR: lower bound",
          col = col,
          ...)
    
    dev.flush()
  } else if (!(x$unif_CR) & show[5L]) {print("uniform credible region not calculated")}
  
  if ((x$unif_CR) & show[6L]) {
    dev.hold()
    
    image(x = x$rescaled_time,
          xlab = "rescaled_time",
          y = x$frequency,
          ylab = "frequency",
          z = x$tvpsd.u05,
          main = "90 percent uniform CR: upper bound",
          col = col,
          ...)
    
    dev.flush()
  } else if (!(x$unif_CR) & show[6L]) {print("uniform credible region not calculated")}
  
  invisible()
  
}

#' a generic method for bdp_dw_result class
#' @param obj object of class bdp_dw_result
#' @export
bayes_factor <- function(obj){UseMethod("bayes_factor")}


#' Extracting the Bayes factor of k1=1 from bdp_dw_result class
#' @param obj object of class bdp_dw_result
#' @return the estimated Bayes factor of k1=1
#' @export
bayes_factor.bdp_dw_result <- function(obj){obj$bayes_factor_k1}




