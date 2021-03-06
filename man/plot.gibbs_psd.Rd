% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psd_S3.R
\name{plot.gibbs_psd}
\alias{plot.gibbs_psd}
\title{Plot method for gibbs_psd class}
\usage{
\method{plot}{gibbs_psd}(x, pdgrm = T, credib = "both", log = T, ...)
}
\arguments{
\item{x}{an object of class gibbs_psd}

\item{pdgrm}{bool flag indicating whether periodogram is visualized or not}

\item{credib}{string indicating which credible regions are visualized. Possible values are "pointwise", "uniform", "both" and "none".}

\item{log}{logical value to determine if the individual spectra are visualized on a log scale}

\item{...}{further arguments to be parsed to \code{plot.default}}
}
\description{
Plot method for gibbs_psd class
}
\details{
Visualizes the spectral density estimate (pointwise posterior median), along with the periodogram
and credibility regions.
If the data has missing values, the periodogram is computed with a linearly
interpolated version of the data using \link[forecast]{na.interp}.
}
