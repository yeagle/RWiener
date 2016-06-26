#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# tests.R
#
# (c) 2016 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2016-06-26
# last mod 2016-06-26 15:48 DW
#

wlrt <- function(wdm1, wdm2) {

  G2 <- 2*logLik(wdm2)-logLik(wdm1)
  Df <- wdm2$npar - wdm1$npar
  pvalue <- pchisq(G2, Df, lower.tail=FALSE)

  res <- list(
  G2 = G2, 
  Df = Df,
  pvalue = pvalue
  )

  return(res)
}

wscoret <- function(wdm1, wdm2) {
  res <- NULL
  return(res)
}

wwaldt <- function(wdm1, wdm2) {
  res <- NULL
  return(res)
}
