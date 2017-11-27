\name{geebt}
\alias{geebt}
\title{Generalized estimating equations for the latency part}
\description{Fit the PH model in the latency part with the GEE approach.}
\usage{geebt(Status, Lambda, X, beta, w, id, corstr)}
\arguments{
  \item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
}
  \item{Lambda}{initial cumulative baseline hazard function from the GEE with independence working corrlation matrix.
}
  \item{X}{a matrix of covariates corresponding to the latency part.
}
  \item{beta}{initial beta for the GEE for the latency part. We use 0 as the initial value. 
}
  \item{w}{conditional probability of a patient remains uncured at the mth iteration. We use Status as initial value.
}
  \item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence} and \code{exchangeable}.
}
}