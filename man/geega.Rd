\name{geega}
\alias{geega}
\title{Generalized estimating equations for the incidence part}
\description{Fit the logistic model in the incidence part with the GEE approach}
\usage{geega(w, Z, gamma, id, corstr)}
\arguments{
  \item{w}{conditional probability of a patient remains uncured at the mth iteration. We use Status as initial value.
}
  \item{Z}{a matrix of covariates corresponding to the incidence part.
}
  \item{gamma}{initial beta for the GEE for the latency part. We use 0 as the initial value.
}
  \item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence} and \code{exchangeable}.
}
}