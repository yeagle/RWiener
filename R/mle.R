# internal function
esvec <- function(x, fpar) {
  if (!is.null(fpar)) {
    if ("alpha" %in% names(fpar))
      x[1] <- NA
    if ("tau" %in% names(fpar)) 
      x[2] <- NA
    if ("beta" %in% names(fpar)) 
      x[3] <- NA
    if ("delta" %in% names(fpar)) 
      x[4] <- NA

    x <- as.numeric(na.omit(x))
  }

  return(x)
}

# internal function
efn <- function(x, data, fpar=NULL) {
  object <- list(data=data)
  par <- numeric(4)

  if (is.null(fpar)) {
    par <- x
  }
  else {
    if("alpha" %in% names(fpar))
      par[1] <- fpar["alpha"]
    else {
      par[1] <- x[1]
      x <- x[-1]
    }
    if("tau" %in% names(fpar)) 
      par[2] <- fpar["tau"]
    else {
      par[2] <- x[1]
      x <- x[-1]
    }
    if("beta" %in% names(fpar)) 
      par[3] <- fpar["beta"]
    else {
      par[3] <- x[1]
      x <- x[-1]
    }
    if("delta" %in% names(fpar)) 
      par[4] <- fpar["delta"]
    else {
      par[4] <- x[1]
      x <- x[-1]
    }
  }
  object$par <- par

  res <- nlogLik.wdm(object)
  return(res)
}

# internal function
eparvec <- function(x, fpar=NULL) {
  res <- numeric(4)
  names(res) <- c("alpha", "tau", "beta", "delta")

  if (!is.null(fpar)) {
    if ("alpha" %in% names(fpar))
      res[1] <- NA
    if ("tau" %in% names(fpar)) 
      res[2] <- NA
    if ("beta" %in% names(fpar)) 
      res[3] <- NA
    if ("delta" %in% names(fpar)) 
      res[4] <- NA
  }

  res[!is.na(res)] <- x
  res[is.na(res)] <- fpar

  return(res)
}

# internal function
estfun <- function(data, fpar=NULL, start=NULL) {

  if (is.null(start))
  {
    start <- c(runif(1,1,2),min(data[,1])/3,runif(1,.2,.8),runif(1,-1,1))
  }

  start <- esvec(start, fpar)

  if (length(fpar)==3)
    onm <- optim(start,efn,data=data,fpar=fpar, method="Brent",
                 lower=-100, upper=100)
  else
    onm <- optim(start,efn,data=data,fpar=fpar, method="Nelder-Mead")
  est <- optim(onm$par,efn,data=data,fpar=fpar, method="BFGS",hessian=TRUE)

  par <- eparvec(est$par, fpar)

  res <- list(
    par = par,
    logLik = -est$value,
    counts = est$counts,
    convergence = est$convergence,
    message = est$message,
    hessian = est$hessian
  )

  return(res)
}
