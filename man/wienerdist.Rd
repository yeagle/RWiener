\name{wienerdist}
\alias{dwiener}
\alias{pwiener}
\alias{qwiener}
\alias{rwiener}
\title{Wiener process distribution functions}
\description{
  \code{dwiener} computes the wiener first passage time density.
  \code{pwiener} computes the CDF for the wiener first passage time density.
  \code{qwiener} computes the quantile for a given CDF value.
  \code{rwiener} generates random quantiles from a wiener process
  distribution.
}
\usage{
  dwiener(q, alpha,tau,beta,delta, resp="upper", give_log=FALSE)
  pwiener(q, alpha,tau,beta,delta, resp="upper")
  qwiener(p, alpha,tau,beta,delta, resp="upper")
  rwiener(n, alpha,tau,beta,delta)

}
\arguments{
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}  
  \item{n}{number of observations.}  
  \item{alpha}{boundary separation parameter.}
  \item{tau}{non-decision time parameter.}
  \item{beta}{bias parameter.}
  \item{delta}{drift rate parameter.}
  \item{resp}{respone: "upper", "lower", or "both"}
  \item{give_log}{function returns log, if this argument is TRUE}
}
%\details{}
%\references{} 
%navarro,fuss  
%blurton, kesselmaier, gondan
%\note{}
%\author{}
\examples{
## calculate density for reactiontime 1.45 and some parameters
dwiener(1.45, 2,0.3,0.5,0)

## calculate CDF for reactiontime 1.45 and some parameters
pwiener(1.45, 2,0.3,0.5,0)

## calculate density for reactiontime 1.45 and some parameters
qwiener(0.5, 2,0.3,0.5,0)

## generate one random value
rwiener(1, 2,0.3,0.5,0)
}
\keyword{dwiener}
\keyword{pwiener}
\keyword{qwiener}
\keyword{rwiener}
