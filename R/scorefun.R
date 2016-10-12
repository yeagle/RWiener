## estfun generic function
estfun <- function(x, ...) {
  UseMethod("estfun")
}

## define estfun for the wdm model object (needed by sctest function)
## same as scorefun, but aggregates by id (persons)
## empirical estimation function (score function) 
estfun.wdm <- function(x, ...) {
  res <- scorefun.wdm(x)
  if("id" %in% names(x$data)) {
    res <- cbind(res, id=x$data$id)
    res <- aggregate(. ~ id, sum, data=as.data.frame(res))[,-1]
  }
  return(res)
}

## empirical estimation function (score function) 
scorefun <- function(x, ...) {
  UseMethod("scorefun")
}

scorefun.wdm <- function(x, ...) {
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

  return(res)
}

## information function
infofun <- function(x, ...) {
  UseMethod("infofun")
}

infofun.wdm <- function(x, ...) {
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
      res[i,1] <- .Call(inlalpha, y[i,1], alpha, tau, beta, delta, lambda, kappa)
      res[i,2] <- .Call(inltau, y[i,1], alpha, tau, beta, delta, lambda, kappa)
      res[i,3] <- .Call(inlbeta, y[i,1], alpha, tau, beta, delta, lambda, kappa)
      res[i,4] <- .Call(inldelta, y[i,1], alpha, tau, beta, delta, lambda, kappa)
    }
    else if (y[i,2] == "upper") {
      res[i,1] <- .Call(inlalpha, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
      res[i,2] <- .Call(inltau, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
      res[i,3] <- .Call(inlbeta, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
      res[i,4] <- .Call(inldelta, y[i,1], alpha, tau, 1-beta, -delta, lambda, kappa)
    }
  }

  return(res)
}
