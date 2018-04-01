\name{emes}
\alias{emes}
\title{Expectation-Maximization (EM) algorithm and Expectation-Solution (ES) algorithm}
\description{EM algorithm is based on Peng et al. (2007) and ES algorithm is based on Niu and Peng (2013). ES algorithm is an estension of the EM algorithm where the M-step of the EM algorithm is replaced by a step requiring the solution of a series of generalised estimating equations. Both algorithm are used for the analysis of survival cure data with potential correlation.}
\usage{emes(Time, Status, X, Z, id, corstr, stdz, esmax, eps)}
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
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence} and \code{exchangeable}.
}
  \item{stdz}{If it is TRUE, all the covariates in the \code{formula} and \code{cureform} are standardized. By default, \code{stdz = FALSE}.
}
  \item{esmax}{specifies the maximum iteration number. If the convergence criterion is not met, the ES iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last ES iteration. The default \code{esmax = 100}.
}
  \item{eps}{tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
}
}
