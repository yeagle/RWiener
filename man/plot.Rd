\name{wiener plot}
\alias{plot.wiener}
\alias{wiener_plot}
\title{Wiener plot function}
\description{
  \code{plot} creates a density plot of correct and wrong responses
  for a given dataset.
}
%\usage{
  %plot(x, ...)
%}
\arguments{
  \item{x}{dataframe with data. Needs a reaction time column and a
  accuracy/response column.}  
  \item{...}{Arguments to be passed to methods, such as graphical
  parameters.}  
}
%\details{}
%\references{} 
%\note{}
%\author{}
\examples{
# generate random data
dat <- rwiener(100,2,.3,.5,0)

# plot
plot(dat)

}
\keyword{plot.wiener}
