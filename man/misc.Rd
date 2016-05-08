\name{miscellaneous wiener model functions}
\alias{is.wiener}
\alias{as.wiener}
\alias{reshape.numdata.wiener}
\alias{reshape}
\alias{reshape.default}
\alias{reshape.data.wiener}
\title{Miscellaneous Wiener Diffusion model functions}
\description{
  \code{is.wiener} bla.
}
\usage{
is.wiener(data)
as.wiener(data, yvar=c("q", "resp"))
\method{reshape}{numdata.wiener}(data, \ldots)
\method{reshape}{data.wiener}(data, \ldots)
}
\arguments{
  \item{data}{can be a data.wiener and/or data.frame with data (needs a reaction time column and a
  accuracy/response column). Further it can be a numdata.wiener and/or
  numeric with the data as single variable (lower bound reaction times are
  then represented as negative numbers, upper bound reaction times as
  positive numbers).}  
  \item{yvar}{represents an optional vector, that can be used to define the
  names of the reaction time column (first value) and the accuracy/response
  column (second value), if a data.wiener and/or data.frame is given as
  data.}
  \item{\ldots}{optional arguments:
  \code{yvar} (as described above) and \code{direction}: character string
  that can be used to define the desired format of the returned data.
  \code{"wide"} returns a numdata.wiener, \code{"long"} returns a
  data.wiener.}
}
\details{
  \code{data.wiener} and \code{numdata.wiener} are data objects that
  represent data coming from a Wiener Diffusion process. \code{data.wiener}
  uses a data.frame with 2 columns for the 2 response variables (\code{"q"}
  and \code{"resp"} by default). \code{numdata.wiener} emulates a single
  variable representation by using a vector, that stores the responses for
  the upper boundary as positive numbers and the responses for the lower
  boundary as negative numbers. 

  The \code{as.wiener} function can be used to create wiener data objects
  (\code{data.wiener} or \code{numdata.wiener}), that can be used by
  generic functions, e.g. \code{plot}.

  \code{is.wiener} checks if the given data is a wiener data object
  (\code{data.wiener} or \code{numdata.wiener}).

  \code{reshape.data.wiener} and \code{reshape.numdat.wiener} can be used
  to transform \code{data.wiener} objects to \code{numdata.wiener} objects
  and vice versa. The generic function \code{reshape}(data, \ldots) can
  be called for convenience.
  See \link[stats]{reshape} for default behavior of the reshape function.
  }
\examples{
## generate data
dat <- rwiener(100, 4,.35,.5,0.8)

## simple function calls
is.wiener(dat)
dat <- as.data.frame(dat)
dat <- as.wiener(dat)
y <- reshape(dat)
y
reshape(y)
}
\keyword{is.wiener}
\keyword{as.wiener}
\keyword{reshape.numdata.wiener}
\keyword{reshape.data.wiener}
