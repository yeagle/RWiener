\name{Wiener Diffusion Model fit}
\alias{wdm}
\title{Wiener Diffusion model estimation function}
\description{
  \code{wdm} creates parameter estimates for the four parameters
  of the wiener model.
}
\usage{
  wdm(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL,
    delta=NULL, xvar=NULL, xvar.par=NULL, start=NULL)
}
\arguments{
  \item{data}{dataframe with data. Needs a reaction time column and a
  accuracy/response column.}  
  \item{start}{vector with the four starting parameter values: alpha, tau, beta, delta.}  
%  \item{fn}{A wiener function to be minimized (or maximized).}  
}
%\details{}
%\references{} 
%\note{}
%\author{}
\examples{
## generate example dataset, fix seed for reproducability
set.seed(1)

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
wdm3 <- wdm(dat, alpha=wdm1$par[1], tau=wdm1$par[2], beta=wdm1$par[3], xvar="group")
}
\keyword{wdm}
