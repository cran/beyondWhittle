\name{beyondWhittle-package}
\alias{beyondWhittle-package}
\alias{beyondWhittle}
\docType{package}
\title{
  Bayesian spectral inference for time series
}
\description{
  Bayesian parametric, nonparametric and semiparametric procedures for spectral density inference of univariate (locally) stationary time series and multivariate stationary time series
}
\details{The package contains several methods (parametric, nonparametric and semiparametric) for Bayesian spectral density inference.
The main algorithms to fit the models for univariate stationary time series are:
\itemize{
\item \link[beyondWhittle]{gibbs_ar}: Parametric, autoregressive (AR) model
\item \link[beyondWhittle]{gibbs_np}: Nonparametric model with Whittle's likelihood and Bernstein-Dirichlet prior from Choudhuri et al (2007)
\item \link[beyondWhittle]{gibbs_npc}: Semiparametric model with corrected AR likelihood and Bernstein-Dirichlet prior from Kirch et al (2018)
}
The package also contains the following models for multivariate stationary time series:
\itemize{
\item \link[beyondWhittle]{gibbs_var}: Parametric, vector autoregressive (VAR) model
\item \link[beyondWhittle]{gibbs_vnp}: Nonparametric model with Whittle's likelihood and Bernstein-Hpd-Gamma prior from Meier (2018)
}
The main function for univariate locally stationary time series is:
\itemize{
\item \link[beyondWhittle]{gibbs_bdp_dw}: Nonparametric model with BDP-DW approach from Tang et al (2023)
}
as well as some useful utility functions.
To get started, it is recommended to consider the examples and documentation of the functions listed above.
The work was supported by DFG grants KI 1443/3-1 and KI 1443/3-2.
}
\author{
Claudia Kirch, Renate Meyer, Matthew C. Edwards, Alexander Meier, Yifu Tang

Maintainer: Renate Meyer <renate.meyer@auckland.ac.nz>
}
\references{
N. Choudhuri, S. Ghosal and A. Roy (2004)
\emph{Bayesian estimation of the spectral density of a time series}
JASA
<doi:10.1198/016214504000000557>

C. Kirch, M. C. Edwards, A. Meier and R. Meyer (2018)
\emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood with a Focus on Bayesian Time Series Analysis}
Bayesian Analysis
<doi:10.1214/18-BA1126>

A. Meier (2018)
\emph{A matrix Gamma process and applications to Bayesian analysis of multivariate time series}
PhD thesis, OvGU Magdeburg
<doi:10.25673/13407>

Y. Tang, C. Kirch, J. E. Lee and R. Meyer (2023)
\emph{Bayesian nonparametric spectral analysis of locally stationary processes}
ArXiv preprint
<arXiv:2303.11561>

}
\keyword{ package }
