\name{wienerdist}
\alias{dwiener}
\alias{pwiener}
\alias{qwiener}
\alias{rwiener}
\title{Wiener process distribution functions}
\description{
  \code{dwiener} computes the wiener first passage time density.
  \code{pwiener} computes the CDF for the wiener first passage time
  density.
  \code{qwiener} computes the quantile for a given CDF value.
  \code{rwiener} generates random quantiles from a wiener process
  distribution, based on the rejection based method.
  For all functions, the standard deviation of the diffusion process is
  fixed to 1.
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
\examples{
## calculate density for reactiontime 1.45, upper bound and some parameters
dwiener(1.45, 2,0.3,0.5,0)

## calculate CDF for reactiontime 1.45, upper bound and some parameters
pwiener(1.45, 2,0.3,0.5,0)

## calculate quantile for CDF value of 0.5, upper bound and some parameters
qwiener(0.5, 2,0.3,0.5,0)

## generate one random value
rwiener(1, 2,0.3,0.5,0)
}
\keyword{dwiener}
\keyword{pwiener}
\keyword{qwiener}
\keyword{rwiener}

\references{
Wabersich, D., & Vandekerckhove, J. (2014). The RWiener package: An R
package providing distribution functions for the Wiener diffusion model.
The R Journal, 6(1), 49-56.
}
