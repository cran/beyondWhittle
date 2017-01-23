\name{beyondWhittle-package}
\alias{beyondWhittle-package}
\alias{beyondWhittle}
\docType{package}
\title{
  Bayesian spectral inference for stationary time series
}
\description{
Implementations of a Bayesian parametric (autoregressive), a Bayesian nonparametric (Whittle likelihood with Bernstein-Dirichlet prior) and a Bayesian semiparametric (autoregressive likelihood with Bernstein-Dirichlet correction) procedure are provided
}
\details{
The work is based on the corrected parametric likelihood by Kirch, Meyer et al. 
It was supported by DFG grant KI 1443/3-1.
See the examples and the referenced paper for further details.
The function \code{gibbs_AR}, \code{gibbs_NP}, \code{gibbs_NPC} correspond to the procedures AR, NP and NPC outlined in the simulation study in the referenced paper.
The other functions are useful utility functions.
}
\author{
Claudia Kirch, Renate Meyer, Matthew C. Edwards, Alexander Meier

Maintainer: Alexander Meier <alexander.meier@ovgu.de>
}
\references{
C. Kirch et al. (2017)
\emph{Beyond Whittle: Nonparametric Correction of a Parametric Likelihood With a Focus on Bayesian Time Series Analysis}
<arXiv:1701.04846>
}
\keyword{ package }
