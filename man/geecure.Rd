\name{geecure}
\alias{geecure}
\title{Marginal proportional hazards mixture cure model with generalzied estimating equations}
\description{Fit the marginal proportional hazards mixture cure (PHMC) model with the generalized estimating equations (GEE). GEE approach is generalized to the marginal PHMC model through the expectation-solution (ES) algorithm to account for the correlation among the cure statuses and the dependence among the failure times of uncured patients to improve the estimation efficiency.
}
\usage{
geecure(formula, cureform, data, id, model = c("para", "semi"),  
        corstr = c("independence", "exchangeable"), Var = TRUE, stdz = FALSE, 
        boots = FALSE, nboot = 100, esmax = 100, eps = 1e-06)
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
  \item{model}{specifies your model, it can be \code{para} which represents the parametric PHMC model with two-parameter Weibull baseline survival function,  or \code{semi} which represents the semiparametric PHMC model.
}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence} and \code{exchangeable}.
}
  \item{Var}{If it is TRUE, the program returns Std.Error by the sandwich method. By default, \code{Var = TRUE}. %If set to False, the program only returns estimators of coefficients. 
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
\details{The marginal PHMC model is considered in this function. For cure rate, a logistic regression model is employed and the probability of being cured is given by \eqn{(1+\exp(\gamma Z))^{(-1)}}. For uncured subject, the failure time is modeled by either the parametric proportional hazards model with Weibull baseline distributions or the semiparametric proportional hazards model. A covariate can be used either in \code{formula} or in \code{cureform} or in both. The model parameters are estimated by the expectation-solution (ES) algorithm and the standard error estimates are obtained from sandwich variance formula based on Niu and Peng (2014) and Niu et al. (2017).}
\value{An object of class \code{geecure} is returned. It can be examined by \code{print.geecure()}.}
\references{
Niu, Y. and Peng, Y. (2014) Marginal regression analysis of clustered failure time data with a cure fraction. \emph{Journal of Multivariate Analysis}, \bold{123}: 129-142.

Niu, Y., Song, L., Liu, Y, and Peng, Y. (2017) Modeling clustered long-term survivors using marginal mixture cure model. under revision.
}

\examples{
\dontshow{# for checks
           data(bmt)
           testdata <- bmt[(1:dim(bmt)[1])[bmt$Z9 == 2 | bmt$Z9 == 3], ] 
           geecuretest <- geecure(Surv(T2, d3) ~ factor(g), cureform = ~ factor(g), 
           data = testdata, model =  "semi", id = testdata$Z9, corstr= "independence")
}
\donttest{# Be patient, the following examples may take several minites on a faster computer.
# Example 1. Fit the marginal parametric PHMC model for the smoking cessation data.
data(smoking)
smoking$Time <- ifelse(smoking$Relapse == 0, smoking$Timept1, 
                       (smoking$Timept1 + smoking$Timept2)/2)

# plot KM survival curve of smoking cessation data
plot(survfit(Surv(Time, Relapse) ~ SexF + (SI.UC), data = smoking), 
      ylab = "Survival function", xlab = "Years", ylim = c(0.5, 1), 
      xlim = c(0, 6), lty = 1:4, col = 1:4)
legend(0.5, 0.63, c("SI/Male", "SI/Female", "UC/Female", "UC/Male"), cex = 1, 
        col = c(2, 4, 3, 1), lty = c(2, 4, 3, 1))

geesmokingind <- geecure(Surv(Time, Relapse) ~ SexF + Duration + SI.UC + F10Cigs + 
      SexF * SI.UC, cureform = ~ SexF + Duration + SI.UC + F10Cigs + SexF * SI.UC, 
      data = smoking,  model = "para", id = smoking$Zip, corstr = "independence")

geesmokingexch <- geecure(Surv(Time, Relapse) ~ SexF + Duration + SI.UC + F10Cigs + 
      SexF * SI.UC, cureform = ~ SexF + Duration + SI.UC + F10Cigs + SexF * SI.UC, 
      data = smoking,  model = "para", id = smoking$Zip, corstr = "exchangeable")

# Example 2. Fit the marginal semiparametric PHMC model for the bmt data.
data(bmt)
# plot KM survival curve of bmt data 
plot(survfit(Surv(T2, d3) ~ 1, data = bmt), xlab = "Days", ylab = "Survival Probability", 
     cex.lab = 1.7, cex.axis = 2, cex.main = 1.7, mark.time = TRUE)

geebmtind <- geecure(Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8, 
             data = bmt, model = "semi", id = bmt$Z9, corstr= "independence")           

geebmtexch <- geecure(Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8, 
              data = bmt, model = "semi", id = bmt$Z9, corstr= "exchangeable", 
              stdz = TRUE, boots = TRUE)
}
}