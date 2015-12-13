\name{wiener estimate}
\alias{estimate}
\title{Wiener model parameter estimation function}
\description{
  \code{estimate} creates parameter estimates for the four parameters
  of the wiener model.
}
\usage{
  estimate(dat, par=c(1,.1,.1,1), fn=wiener_deviance)
}
\arguments{
  \item{dat}{dataframe with data. Needs a reaction time column and a
  accuracy/response column.}  
  \item{par}{vector with the four starting parameter values: alpha, tau, beta, delta.}  
  \item{fn}{A wiener function to be minimized (or maximized).}  
}
%\details{}
%\references{} 
%\note{}
%\author{}
\examples{
# generate random data
dat <- rwiener(100,2,.3,.5,0)

# optim
estimate(dat)

}
\keyword{estimate}
