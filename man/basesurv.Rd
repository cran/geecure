\name{basesurv}
\alias{basesurv}
\title{Estimation of the baseline survival function}
\description{The estimated baseline survival function based on the product-limit estimator (Kalbfleisch and Prentice, 2002), which is uesd to update the E-step in the ES algorithm.}
\usage{
basesurv(Time, Status, X, beta, Lambda, w, model)
}
\arguments{
  \item{Time}{right censored data which is the follow up time.
}
  \item{Status}{the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
}
  \item{X}{a matrix of covariates corresponding to the latency part.
}
  \item{beta}{initial beta from the GEE for the latency part.  
}
  \item{Lambda}{initial cumulative baseline hazard function from the GEE with independence working corrlation matrix.
}
  \item{w}{conditional probability of a patient remains uncured at the mth iteration. We use Status as initial value.
}
  \item{model}{specifies your model, it can be \code{para} which represents the parametric PHMC model with two-parameter Weibull baseline survival function,  or \code{semi} which represents the semiparametric PHMC model.
}
}
\references{
Kalbfleisch, J. D. and Prentice, R. L. (2002) \emph{The Statistical Analysis of Failure Time Data}. John Wiley & Sons, New York, 2nd edition.
}