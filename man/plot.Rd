\name{plot}
\alias{plot.data.wiener}
\alias{wiener_plot}
\title{Wiener plot function}
\description{
  \code{plot} creates a density plot of correct and wrong responses
  for a given dataset.
}
\usage{
\method{plot}{data.wiener}(x, ...)
}
\arguments{
  \item{x}{data.wiener object, which is basically a data.frame with data.
  Needs a reaction time column and a accuracy/response column.}  
  \item{...}{Arguments to be passed to methods, such as graphical
  parameters.}  
}
%\details{}
%\references{} 
%\note{}
%\author{}
\examples{
## generate random data
dat <- rwiener(100,2,.3,.5,0)

## plot
plot(dat)
}
\keyword{plot.data.wiener}

\references{
Wabersich, D., & Vandekerckhove, J. (2014). The RWiener package: An R
package providing distribution functions for the Wiener diffusion model.
The R Journal, 6(1), 49-56.
}
