\name{miscellaneous wiener model functions}
\alias{wdm}
\title{Wiener model parameter estimation function}
\description{
  \code{wdm} creates parameter estimates for the four parameters
  of the wiener model.
}
\usage{
  wdm(data, alpha=NULL, tau=NULL, beta=NULL, delta=NULL, start=NULL)
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
# generate random data
dat <- rwiener(100,2,.3,.5,0)

# optim
wdm(dat)

}
\keyword{estimate}
