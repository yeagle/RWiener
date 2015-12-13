#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# mle.R
#
# (c) 2015 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2015-12-12
# last mod 2015-12-13 15:25 DW
#

estimate <- function(dat, par=c(1,.1,.1,1), fn=wiener_deviance) {
  est <- optim(par,fn,dat=dat, method="Nelder-Mead")
  pars <- est$par
  names(pars) <- c("alpha", "tau", "beta", "delta")
  return(pars)
}
