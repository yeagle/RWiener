#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# likelihood.R
#
# (c) 2015 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2015-12-12
# last mod 2015-12-12 12:31 DW
#

wiener_likelihood <- function(x, dat) {
 if (!check_wiener_pars(x[1],x[2],x[3],x[4])) {
    return(-Inf)
  }
  ll <- vector("double", length(dat[,1]))
  for (i in 1:length(dat[,1])) {
    ll[i] <- dwiener(as.double(dat[i,1]), x[1],x[2],x[3],x[4], 
                  resp=as.character(dat[i,2]), give_log=TRUE)
  }
  return(sum(ll))
}
  
wiener_deviance <- function(x, dat) {
  -2*wiener_likelihood(x,dat)
}

wiener_bic <- function(x, dat, loss=NULL) {
  if(is.null(loss)) {
    -2*wiener_likelihood(x,dat)+4*log(length(dat[,1]))
  }
  else {
    if(is.list(dat)) {
      loss(x,dat)+length(x)*log(length(dat[[1]][,1]))
    }
    else if (is.data.frame(dat)) {
      loss(x,dat)+length(x)*log(length(dat[,1]))
    }
    else {
      stop("don't know how to handle the dat object!")
    }

  }
}

wiener_aic <- function(x, dat, loss=NULL) {
  if(is.null(loss)) {
    -2*wiener_likelihood(x,dat)+4*2 
  }
  else {
    loss(x,dat)+length(x)*2 
  }
}
