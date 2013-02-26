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
  wiener_likelihood(x, dat)
  wiener_deviance(x, dat)
  wiener_aic(x, dat)
  wiener_bic(x, dat)

}
\arguments{
  \item{x}{vector with the four parameter values.}
  \item{dat}{dataframe with data. Needs a reaction time column and a
  accuracy/response column.}  
}
%\details{}
%\references{} 
%navarro,fuss  
%blurton, kesselmaier, gondan
%\note{}
%\author{}
\examples{
# generate random data
dat <- rwiener(100,2,.3,.5,0)

#compute likelihood
wiener_likelihood(c(2,.3,.5,0), dat)

# estimate parameters with optim
onm <- optim(c(1,.1,.1,1),wiener_deviance,dat=dat, method="Nelder-Mead")
est <- optim(onm$par,wiener_deviance,dat=dat, method="BFGS",hessian=TRUE)
est$par # parameter estimates

# the following code needs the MASS package 
#library(MASS) # ginv - inverse of a matrix - needed
#sqrt(diag(ginv(est$hessian))) # sd for paramerets

}
\keyword{wiener_likelihood}
\keyword{wiener_deviance}
\keyword{wiener_aic}
\keyword{wiener_bic}
