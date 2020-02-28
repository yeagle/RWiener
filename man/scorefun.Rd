\name{scorefun}
\alias{scorefun}
\alias{scorefun.wdm}
\alias{estfun}
\alias{estfun.wdm}
\title{Extract Empirical Estimating Functions}
\description{
 !EXPERIMENTAL FUNCTION!

 Generic function for extracting the empirical estimating functions
 of a fitted model. 

 !EXPERIMENTAL FUNCTION!
}
\usage{
scorefun(x, \dots)
estfun(x, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{\dots}{arguments passed to methods.}
}

\value{A matrix containing the empirical estimating functions.
Typically, this should be an \eqn{n \times k}{n x k} matrix corresponding
to \eqn{n} observations and \eqn{k} parameters. The columns should be named
as in \code{\link{coef}} or \code{\link{terms}}, respectively.

The estimating function (or score function) for a model - \code{scorefun} - is the derivative of the objective function
with respect to the parameter vector. The empirical estimating functions is
the evaluation of the estimating function at the observed data (\eqn{n} observations)
and the estimated parameters (of dimension \eqn{k}).

The \code{estfun} function is basically the score function, but with the
additional functionality to sum up the scores by the given covariable
\code{id} in the dataset.
}

\references{
  Zeileis A (2006), Object-Oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\examples{
## generate random data
dat <- rwiener(100,1,.2,.5,0.5)

## fit wdm
wdm1 <- wdm(dat)

## estimating function
scores <- scorefun(wdm1)

## print
head(scores)

## plot
par(mfrow=c(2,2))
plot(scores[,1]);plot(scores[,2]);plot(scores[,3]);plot(scores[,4])
}
