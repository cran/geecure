\name{geecure2}
\alias{geecure2}
\title{Semiparametric marginal proportional hazards mixture cure model}
\description{Fit the semiparametric marginal proportional hazards mixture cure (PHMC) model for clustered failure time data. The function is based on the methods proposed by Peng et al. (2007) and Niu and Peng (2013). }
\usage{
geecure2(formula, cureform, data, id, corstr = c("independence", "exchangeable"), 
Var = TRUE, stdz = FALSE, boots = FALSE, nboot = 100, esmax = 100, eps = 1e-06)
}
\arguments{
  \item{formula}{a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. It is used to specify the covariate effects on the failure time of uncured subjects. See the documentation for survreg, Surv for details. The expression to the right of the "~" specifies the effect of covariates on the failure time of uncured patients.
}
  \item{cureform}{a formula expression, of the form \code{cureform ~ predictors}. It is used to specify the effects of covariates on the cure rate. A covariate may be used in both \code{formula} and \code{cureform}. 
}
  \item{data}{a data frame in which to interpret the variables named in the \code{formula} and the \code{cureform}.
}
  \item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence} and \code{exchangeable}.
}
  \item{Var}{If it is TRUE, the program returns Std.Error by the sandwich method. By default, \code{Var = TRUE}. 
}
  \item{stdz}{If it is TRUE, all the covariates in the \code{formula} and \code{cureform} are standardized. By default, \code{stdz = FALSE}.
}
  \item{boots}{If it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}.
}
  \item{nboot}{the number of bootstrap samples. The default is \code{nboot = 100}.
}
  \item{esmax}{specifies the maximum iteration number. If the convergence criterion is not met, the ES iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last ES iteration. The default \code{esmax = 100}.
}
  \item{eps}{tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
}
}
\details{The semiparametric marginal PHMC model is considered in this function. For cure rate, a logistic regression model is employed and the probability of being cured is given by \eqn{(1+\exp(\gamma Z))^{(-1)}}. For uncured subject, the failure time is modeled by the semiparametric proportional hazards model. A covariate can be used either in \code{formula} or in \code{cureform} or in both. When \code{corstr = independence}, the model parameters are estimated by the expectation-maximization (EM) algorithm and the standard error estimates are obtained from sandwich variance formula based on Peng et al. (2007). When \code{corstr = exchangeable}, \code{stdz = TRUE} and \code{boots = TRUE}, the model parameters are estimated by the expectation-solution (ES) algorithm and the standard error estimates are obtained from bootstrap variance formula based on and Niu et al. (2013).}
\value{An object of class \code{geecure2} is returned. It can be examined by \code{print.geecure2()}.}
\references{
Peng, Y., Taylor, J. M. G., and Yu, B. (2007) A marginal regression model for multivariate failure time data with a surviving fraction. \emph{Lifetime Data Analysis}, \bold{13}: 351-369.

Niu, Y. and Peng, Y. (2013) A semiparametric marginal mixture cure model for clustered survival data. \emph{Statistics in Medicine}, \bold{32}: 2364-2373.
}
\examples{
\dontshow{# for checks
           data(bmt)
           testdata <- bmt[(1:dim(bmt)[1])[bmt$Z9 == 2 | bmt$Z9 == 3], ] 
           geecuretest2 <- geecure2(Surv(T2, d3) ~ Z8, cureform = ~ Z8, data = testdata, 
           id = testdata$Z9, corstr= "independence")
}
\donttest{
# Example. Fit the marginal semiparametric PHMC model for the bmt data.
data(bmt)
geebmtind2 <- geecure2(Surv(T2, d3) ~ Z8, cureform = ~ Z8, data = bmt, id = bmt$Z9, 
              corstr= "independence")           
geebmtexch2 <- geecure2(Surv(T2, d3) ~ Z8, cureform = ~ Z8, data = bmt, id = bmt$Z9, 
               corstr= "exchangeable", stdz = TRUE, Var = FALSE)
}
}