\name{tests}
\alias{anova.wdm}
\alias{print.wlrt}
\alias{waldtest}
\alias{waldtest.wdm}
\alias{print.wwaldt}
\title{Wiener Diffusion model test functions}
\description{
  Calculates test scores and further information for \code{\link{wdm}}
  model objects.
}
\usage{
  \method{anova}{wdm}(object, \ldots, test="LRT")
}
\arguments{
  \item{object}{a wdm model object.}
  \item{test}{Statistical test to calculate, so far the only option is a
  Likelihood Ratio Test (LRT).}
 \item{\dots}{Further model objects or other arguments passed to methods.}
}
\details{
  The \code{anova.wdm} function calls the specified test and calculates the
  test statistics and other information for two or more nested
  \code{\link{wdm}} model objects.
}
%\references{} 
%\note{}
%\author{}
\examples{
# a random dataset
dat <- rbind(cbind(rwiener(100, 2,.3,.5,0), grp="A"),
cbind(rwiener(100,2,.3,.5,1), grp="B"))

# create nested wdm model objects (from specific to general)
wdm1 <- wdm(dat)
wdm2 <- wdm(dat,
alpha=coef(wdm1)[1], tau=coef(wdm1)[2], beta=coef(wdm1)[3],
xvar="grp")
wdm3 <- wdm(dat,
tau=coef(wdm1)[2],
xvar="grp")

# conduct LRT tests
anova1 <- anova(wdm1,wdm2,wdm3)
anova1

# extract exact p-value
anova1$pvalue
}
\keyword{anova.wdm}
