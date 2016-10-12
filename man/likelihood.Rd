\name{likelihood}
\alias{logLik.wdm}
\alias{deviance.wdm}
\alias{AIC.wdm}
\alias{BIC.wdm}
\title{Likelihood and criterion functions for wdm}
\description{
  \code{logLik.wdm} computes the log-likelihood.
  \code{deviance.wdm} computes the deviance.
  \code{AIC.wdm} computes the AIC.
  \code{BIC.wdm} computes the BIC.
}
\usage{
\method{logLik}{wdm}(object, \ldots)
\method{deviance}{wdm}(object, \ldots)
\method{AIC}{wdm}(object, \ldots)
\method{BIC}{wdm}(object, \ldots)
}
\arguments{
  \item{object}{a wdm object file or a list containing a \code{$par}
  vector with the model parameters, a \code{$data} data.frame with the data
  and optionally a \code{$loss} function.}
  \item{\ldots}{optional arguments}
}
\details{
  The \code{$par} vector with the (four) parameter values should be in the following
  order: alpha, tau, beta, delta.
  
  The \code{$data} data.frame with data needs a reaction time column and a
  accuracy/response column.

  The \code{$loss} function defaults to NULL, which means that the default
  computation is done by using the default formula.
  If not NULL, this can be a function to replace the default
  computation in the code and use a custom function
  instead. The custom function takes two arguments: the parameter vector
  and the data.frame with the data.

  These functions can be very useful in combination with the optim funcion,
  to estimate parameters manually. Check the examples below to see how to
  use the provided generic functions in a manual estimation routine.
}

\examples{
## generate random data
dat <- rwiener(100,3,.25,.5,0.8)

## fit wdm
wdm1 <- wdm(dat, alpha=3, tau=.25, beta=0.5)

## compute likelihood, AIC, BIC, deviance
logLik(wdm1)
AIC(wdm1)
BIC(wdm1)
deviance(wdm1)

## estimate parameters by calling optim manually
## first define necessary wrapper function
nll <- function(x, data) { 
  object <- wdm(data, alpha=x[1], tau=x[2], beta=x[3], delta=x[4])
  -logLik(object)
}
## call estimation routine
onm <- optim(c(1,.1,.1,1),nll,data=dat, method="Nelder-Mead")
est <- optim(onm$par,nll,data=dat, method="BFGS",hessian=TRUE)
est$par # parameter estimates
## use the obtained parameter estimates to create wdm object
wdm2 <- wdm(dat, alpha=est$par[1], tau=est$par[2], beta=est$par[3],
  delta=est$par[4])
## now the generic functions can be used again
logLik(wdm2)
}
\keyword{logLik.wdm}
\keyword{deviance.wdm}
\keyword{AIC.wdm}
\keyword{BIC.wdm}

\references{
Wabersich, D., & Vandekerckhove, J. (2014). The RWiener package: An R
package providing distribution functions for the Wiener diffusion model.
The R Journal, 6(1), 49-56.
}
