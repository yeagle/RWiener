\name{scorefun}
\alias{scorefun}
\alias{scorefun.wdm}
\alias{infofun}
\alias{infofun.wdm}
\title{Extract Empirical Estimating Functions}
\description{
 Generic function for extracting the empirical estimating functions
 of a fitted model.
}
\usage{
scorefun(x, \dots)
infofun(x, \dots)
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

The information function for a model - \code{infofun} - is the 2nd derivative of the objective
function with respect to the parameter vector. The empirical information functions is
the evaluation of the information function at the observed data (\eqn{n} observations)
and the estimated parameters (of dimension \eqn{k}).
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
infos <- infofun(wdm1)

## print
head(scores)
head(infos)
}
