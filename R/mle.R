estfun <- function(dat, fn=nlogLik.wiener, start=NULL) {

  if(!is.wiener(dat)) stop("supplied data is not of class dat.wiener!")
  if(!is.data.frame(dat)) dat <- reshape.wiener(dat)

  if (is.null(start))
  {
    start <- c(runif(1,0,2),min(dat$q)/2,runif(1,.1,.9),runif(1,-1,1))
  }

  onm <- optim(start,fn,dat=dat, method="Nelder-Mead")
  est <- optim(onm$par,fn,dat=dat, method="BFGS",hessian=TRUE)

  rval <- list(
    par = c("alpha"=est$par[1], "tau"=est$par[2], "beta"=est$par[3], "delta"=est$par[4]),
    value = est$value,
    counts = est$counts,
    convergence = est$convergence,
    message = est$message,
    hessian = est$hessian
  )

  return(rval)
}
