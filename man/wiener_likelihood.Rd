\name{wiener likelihood and criterion}
\alias{wiener_likelihood}
\alias{wiener_deviance}
\alias{wiener_aic}
\alias{wiener_bic}
\title{Wiener likelihood and criterion functions}
\description{
  \code{wiener_likelihood} computes the log-likelihood for given parameter
  values and data.
  \code{wiener_deviance} computes the deviance.
  \code{wiener_aic} computes the AIC.
  \code{wiener_bic} computes the BIC.
  These functions can be very useful in combination with the optim funcion,
  to estimate parameters (see example below).
}
\usage{
  wiener_likelihood(x, data)
  wiener_deviance(x, data)
  wiener_aic(x, data, loss=NULL)
  wiener_bic(x, data, loss=NULL)

}
\arguments{
  \item{x}{vector with the four parameter values: alpha, tau, beta, delta.}
  \item{data}{dataframe with data. Needs a reaction time column and a
  accuracy/response column.}  
  \item{loss}{Defaults to NULL, which means that the default computation is
  done by using -2*wiener_likelihood(x,data) in the formula. If not NULL,
  this can be a function to replace the default -2*wiener_likelihood(x,data)
  in the code and use a custom function instead.}
}
%\details{}
%\references{} 
%\note{}
%\author{}
\examples{
### Example 1: Parameter estimation
# generate random data
dat <- rwiener(100,2,.3,.5,0)

#compute likelihood
wiener_likelihood(c(2,.3,.5,0), dat)

# estimate parameters with optim
# onm <- optim(c(1,.1,.1,1),wiener_deviance,data=dat, method="Nelder-Mead")
# est <- optim(onm$par,wiener_deviance,data=dat, method="BFGS",hessian=TRUE)
# est$par # parameter estimates
#
# # the following code needs the MASS package 
# sqrt(diag(MASS::ginv(est$hessian))) # sd for parameters

### Example 2: Simple model comparison
# compare two models with deviance
wiener_deviance(c(3,.3,.5,0), dat)
wiener_deviance(c(3,.3,.5,0.5), dat)
# log-likelihood difference
wiener_likelihood(c(3,.3,.5,0), dat)-wiener_likelihood(c(3,.3,.5,0.5), dat)

### Example 3: Likelihood-Ratio-Test and Wald-Test
# # Suppose we have data from 2 conditions
# dat1 <- rwiener(100,2,.3,.5,-.5)
# dat2 <- rwiener(100,2,.3,.5,.5)
# onm1 <- optim(c(1,.1,.1,1),wiener_deviance,data=dat1, method="Nelder-Mead")
# est1 <- optim(onm1$par,wiener_deviance,data=dat1, method="BFGS",hessian=TRUE)
# wiener_likelihood(est1$par,dat1)+wiener_likelihood(est1$par,dat2) # combined loglike
# model_ll <- function(pars,delta,dat1,dat2) {
#   wiener_likelihood(pars,dat1)+
#     wiener_likelihood(c(pars[1:3],pars[4]+delta),dat2)
# }
# ## Likelihood ratio test
# # 0-model: delta=0; alt-model: delta=1
# model_ll(est1$par,1,dat1,dat2)
# # compute Likelihood ratio
# LR <- -2*model_ll(est1$par,0,dat1,dat2)+2*model_ll(est1$par,1,dat1,dat2)
# # compare with critical X^2(1) quantile, alpha=0.05
# LR > qchisq(0.95,1)
# # get p-value from X^2(1)
# pchisq(LR,1, lower.tail=FALSE)
# ## Wald-Test
# # estimate parameter delta and test for significance
# onm2 <- optim(c(1,.1,.1,1),wiener_deviance,data=dat2, method="Nelder-Mead")
# est2 <- optim(onm2$par,wiener_deviance,data=dat2, method="BFGS",hessian=TRUE)
# delta <- est2$par[4]-est1$par[4]
# # the following code needs the MASS package 
# est1.sd <- sqrt(diag(MASS::ginv(est1$hessian))) # sd for parameters
# WT <- (est1$par[4]-(est1$par[4]+delta))/est1.sd[4]
# # compare with critical quantile N(0,1), alpha=0.05
# abs(WT) > qnorm(0.975)
# # get p-value from N(0,1)
# pnorm(WT)

### Example 4: Custom AIC loss function
many_drifts <- function(x,datlist) {
  l = 0
  for (c in 1:length(datlist)) {
    l = l + wiener_deviance(x[c(1,2,3,c+3)],datlist[[c]])
  }
  return(l)
}
dat1 <- rwiener(n=100, alpha=2, tau=.3, beta=.5, delta=0.5)
dat2 <- rwiener(n=100, alpha=2, tau=.3, beta=.5, delta=1)
datlist <- list(dat1,dat2)
wiener_aic(x=c(2,.3,.5,.5,1), data=datlist, loss=many_drifts)

}
\keyword{wiener_likelihood}
\keyword{wiener_deviance}
\keyword{wiener_aic}
\keyword{wiener_bic}

\references{
Wabersich, D., & Vandekerckhove, J. (2014). The RWiener package: An R
package providing distribution functions for the Wiener diffusion model.
The R Journal, 6(1), 49-56.
}
