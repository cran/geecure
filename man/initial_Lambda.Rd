\name{initial_Lambda}
\alias{initial_Lambda}
\title{Initial value of the cumulative baseline hazard function}
\description{Obtain the initial value of the cumulative baseline hazard function in the latency part through the GEE with the \code{independence} working correlation matrix.}
\usage{
initial_Lambda(Time, Status, X, Z, id, model, corstr)
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
  \item{model}{specifies your model, it can be \code{para} which represents the parametric PHMC model with two-parameter Weibull baseline survival function,  or \code{semi} which represents the semiparametric PHMC model.
}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence} and \code{exchangeable}.
}
}