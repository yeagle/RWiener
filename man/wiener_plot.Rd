\name{wiener plot}
\alias{wiener_plot}
\title{Wiener plot function}
\description{
  \code{wiener_plot} creates a density plot of correct and wrong responses
  for a given dataset.
}
\usage{
  wiener_plot(dat)
}
\arguments{
  \item{dat}{dataframe with data. Needs a reaction time column and a
  accuracy/response column.}  
}
%\details{}
%\references{} 
%\note{}
%\author{}
\examples{
# generate random data
dat <- rwiener(100,2,.3,.5,0)

# plot
wiener_plot(dat)

}
\keyword{wiener_plot}
