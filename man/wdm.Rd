\name{wiener diffusion model fit}
\alias{wdm}
\alias{print.wdm}
\title{Wiener Diffusion model fit function}
\description{
  \code{wdm} creates parameter estimates for the four parameters
  of the wiener model.
}
\usage{
  wdm(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL,
    delta=NULL, xvar=NULL, xvar.par=NULL, start=NULL)
}
\arguments{
  \item{data}{is the data object containing data coming from a
  (hypothetical) Wiener Diffusion process. For further details on the data
  object, see \link{is.wiener}.}
  \item{yvar}{represents an optional vector, that can be used to define the
  names of the reaction time column. For further details on the data
  object, see \link{is.wiener}.}
  \item{alpha}{optional, can be used to fix the alpha parameter to the
  given value.}
  \item{tau}{optional, can be used to fix the tau parameter to the given
  value.}
  \item{beta}{optional, can be used to fix the beta parameter to the given
  value.}
  \item{delta}{optional, can be used to fix the delta parameter to the
  given value.}
  \item{xvar}{NULL (not yet fully implemented)}
  \item{xvar.par}{NULL (not yet implemented)}
  \item{start}{an optional vector with the four starting parameter values
  in the following order: alpha, tau, beta, delta.}  
}
\details{
  The \code{wdm} function calls an estimation routine, to estimate the
  model parameters. 
  
  If all but one parameters are fixed, a \code{"Brent
  (optim)"} type algorithm is used. For the estimation of more than one
  parameter, first a \code{"BFGS (optim)"} type algorithm is tried, if
  unsuccessful, a \code{"Newton type (nlm)"} algorithm is tried, if again
  unsuccessful, a \code{"Nelder-Mead (optim)"} algorithm is used.

  In case all parameters are set to fixed values, no estimation routine is
  called, but a \code{wdm} object will still be created.

  The returned \code{wdm} object is basically a list containing the
  parameter estimates in \code{$coefficients}. \code{$hessian} contains the
  numerically differentiated Hessian matrix (if available, else NULL).
  \code{$data} contains the data passed to the \code{wdm} function call.
  \code{$loglik} contains the log-likelihood value for the \code{wdm}
  object and its parameter estimates. \code{$estpar} contains a vector,
  that is \code{TRUE} if the respective parameter was estimated and
  \code{FALSE} if the respective parameter was set to a fixed value.
  Additional information is given in other list objects.
  }
%\references{} 
%\note{}
%\author{}
\examples{
\dontshow{set.seed(1)}%
## generate random data
d1 <- rwiener(100, 1,.3,.5,1)
d2 <- rwiener(100, 1,.3,.5,-1)
d1$group <- factor('g1', levels=c('g1', 'g2'))
d2$group <- factor('g2', levels=c('g1', 'g2'))
dat <- rbind(d1,d2)

## fit wdm
wdm1 <- wdm(dat)

## extract sd
sqrt(diag(MASS::ginv(wdm1$hessian)))

## fit wdm
wdm2 <- wdm(dat, beta=.5)
wdm3 <- wdm(dat, alpha=wdm1$coefficients[1], tau=wdm1$coefficients[2],
beta=wdm1$coefficients[3],
  xvar="group")

\dontshow{
if (!(round(sum(wdm1$coefficients),6) == 1.548682)) 
  stop("Example calculations changed despite fixed seed")
if (!(round(sum(wdm2$coefficients),6) == 1.58825)) 
  stop("Example calculations changed despite fixed seed")
if (!(round(sum(wdm3$coefficients),6) == 1.395183)) 
  stop("Example calculations changed despite fixed seed")
}
}
\keyword{wdm}
