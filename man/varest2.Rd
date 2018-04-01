\name{varest2}
\alias{varest2}
\title{Variance estimate with sandwich formula based on the EM algorithm}
\description{Calculate the variance estimates using the sandwich formula based on the EM algorithm.}
\usage{
varest2(Time, Status, X, Z, id, gamma, beta, bsurv, w)
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
  \item{bsurv}{the estimate of the baseline survival function for the latency part.
}
  \item{w}{conditional probability of a patient remains uncured.
}
}