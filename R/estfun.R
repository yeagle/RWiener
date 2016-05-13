## empirical estimation function (score function) 
estfun <- function(x, ...) {
  UseMethod("estfun")
}

estfun.wdm <- function(x, ...) {
  y <- x$data[,x$yvar]

  alpha <- x$coefficients["alpha"]
  tau <- x$coefficients["tau"]
  beta <- x$coefficients["beta"]
  delta <- x$coefficients["delta"]

  n <- length(y[,1])
  res <- matrix(rep(NA,4*n), n,4)
  colnames(res) <- c("alpha", "tau", "beta", "delta")
  for (i in 1:n) {

    ## kappa and lambda calculations
    kst <- .Call(kappaST, y[i,1])
    klt <- .Call(kappaLT, y[i,1])
    lambda <- kst - klt
    if (lambda < 0) kappa <- kst
    else kappa <- klt

    if (y[i,2] == "lower") {
      res[i,1] <- .Call(sclalpha, y[i,1], alpha, tau, beta, delta, lambda, kappa)
      res[i,2] <- .Call(scltau, y[i,1], alpha, tau, beta, delta, lambda, kappa)
      res[i,3] <- .Call(sclbeta, y[i,1], alpha, tau, beta, delta, lambda, kappa)
      res[i,4] <- .Call(scldelta, y[i,1], alpha, tau, beta, delta, lambda, kappa)
    }
    else if (y[i,2] == "upper") {
      res[i,1] <- .Call(sclalpha, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
      res[i,2] <- .Call(scltau, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
      res[i,3] <- .Call(sclbeta, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
      res[i,4] <- .Call(scldelta, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
    }
  }

  ## aggregate scores by id
  if("id" %in% names(x$data)) {
    res <- cbind(res, id=x$data$id)
    res <- aggregate(. ~ id, sum, data=as.data.frame(res))[,-1]
  }

  return(res)
}
