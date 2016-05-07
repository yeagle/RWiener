#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# wdm.R
#
# (c) 2016 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2016-05-07
# last mod 2016-05-07 20:29 DW
#

wdm <- function(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL, delta=NULL,
               xvar=NULL, xvar.par=NULL, start=NULL) {
  # save original function call
  cl <- match.call()

  # prepare passed arguments
  verifydata(data)
  if (is.numeric(data) & is.null(xvar)) data <- reshape.wiener(data, yvar=yvar)
  else if (length(yvar)==1) {
    cbind(reshape.wiener(data[,yvar]),data)
    yvar <- c("q", "resp")
  }
  fpar <- c("alpha"=unname(alpha), "tau"=unname(tau), 
    "beta"=unname(beta), "delta"=unname(delta))

  # estimate parameters
  if (!is.null(xvar)) {
    if(length(xvar)==1) {
      if(class(data[,xvar]) == "factor"){
        res <- list()
        res$par <- fpar
        for (l in levels(data[,xvar])) {
          est <- estfun(data[data[,xvar]==l,yvar], fpar, start)
          est$par <- est$par[!(names(est$par) %in% names(fpar))]
          names(est$par) <- paste(l,names(est$par), sep=":")
          res$par <- append(res$par, est$par)
          res$counts <- append(res$counts, est$counts)
          res$convergence <- append(res$convergence, est$convergence)
          res$message <- append(res$message, est$message)
          res$hessian <- append(res$hessian, est$hessian)
          res$loglik <- sum(res$loglik, est$loglik)
        }
      }
    }
  }
  else
    res <- estfun(data[,yvar], fpar, start)

  # prepare return object
  res$n <- length(data[,1])
  res$npar <- length(res$par)
  res$data <- data
  res$yvar <- yvar
  res$estpar <- c("alpha"=is.null(alpha), "tau"=is.null(tau),
                   "beta"=is.null(beta), "delta"=is.null(delta))
  res$call <- cl
  class(res) <- c("wdm")
  return(res)
}
