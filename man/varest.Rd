\name{varest}
\alias{varest}
\title{Variance estimate with sandwich formula based on the ES algorithm}
\description{Calculate the variance estimates using the sandwich formula based on the ES algorithm.}
\usage{
varest(Time, Status, X, Z, id, gamma, beta, kappa, gphi, gcor, bphi, bcor, 
       Lambda, w, model)
}
\arguments{
  \item{Time}{right censored data which is the follow up time.
}
  \item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.
}
  \item{X}{a matrix of covariates corresponding to the latency part.
}
  \item{Z}{a matrix of covariates corresponding to the incidence part.
}
  \item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
}
  \item{gamma}{the estimates for the incidence part.
}
  \item{beta}{the estimates for the latency part.
}
  \item{kappa}{the estimate of the shape parameter in the Weibull baseline hazard function when \code{model = "para"}.
}
  \item{gphi}{the estimate of the scale parameter \eqn{\phi_1} in the GEE for the incidence part.
}
  \item{gcor}{the estimate of the correlation parameter \eqn{\rho_1} in the GEE for the incidence part.
}
  \item{bphi}{the estimate of the scale parameter \eqn{\phi_2} in the GEE for the latency part.
}
  \item{bcor}{the estimate of the correlation parameter \eqn{\rho_2} in the GEE for the latency part.
}
  \item{Lambda}{the estimate of the cumulative baseline hazard function in the GEE for the latency part.
}
  \item{w}{conditional probability of a patient remains uncured.
}
  \item{model}{specifies your model, it can be \code{para} which represents the parametric PHMC model with two-parameter Weibull baseline survival function,  or \code{semi} which represents the semiparametric PHMC model.
}
}