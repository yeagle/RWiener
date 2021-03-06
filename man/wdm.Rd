\name{wdm}
\alias{wdm}
\alias{print.wdm}
\alias{vcov.wdm}
\alias{confint.wdm}
\alias{summary.wdm}
\alias{print.summary.wdm}
\title{Wiener diffusion model fit function}
\description{
  \code{wdm} creates parameter estimates for the four parameters
  of the Wiener model.
}
\usage{
  wdm(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL,
    delta=NULL, xvar=NULL, start=NULL, fixed=0)
  \method{vcov}{wdm}(object, \ldots, method="hessian")
}
\arguments{
  \item{data}{is the data object containing data coming from a
  (hypothetical) Wiener diffusion process. For further details on the data
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
  \item{xvar}{optional: group factor variable to estimate all unfixed
  parameters separate for the given groups.}
  \item{start}{an optional vector with the four starting parameter values
  in the following order: alpha, tau, beta, delta.}  
  \item{fixed}{a number indicating how many of the parameters are fixed
    (not free). This number will be subtracted from the number of free
    parameters. Defaults to 0.}  
  \item{method}{the method to use for estimating the covariance matrix of
    the parameter estimates. Options are \code{"opg"} for outer product of
  gradients or \code{"hessian"} to use the hessian matrix from the
  estimation routine. Defaults to \code{"hessian"}.}  
  \item{object}{a wdm object file or a list containing a \code{$par}
  vector with the model parameters, a \code{$data} data.frame with the data
  and optionally a \code{$loss} function.}
 \item{\dots}{arguments passed to methods.}
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

  The standard R functions \code{coef}, \code{vcov}, \code{confint},
  \code{summary} can be used with \code{wdm} objects.
  }
%\references{} 
%\note{}
%\author{}
\examples{
## generate random data
dat <- rbind(cbind(rwiener(100, 2,.3,.5,1), group=factor("A", c("A","B"))),
        cbind(rwiener(100,2,.3,.5,-1), group=factor("B", c("A", "B"))))

## fit wdm
wdm1 <- wdm(dat)

## extract parameters
coef(wdm1)

## further models
wdm2 <- wdm(dat, beta=.5)
wdm3 <- wdm(dat, alpha=wdm1$coefficients[1], tau=wdm1$coefficients[2],
            beta=wdm1$coefficients[3], xvar="group")
}
\keyword{wdm}
