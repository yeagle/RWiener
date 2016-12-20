## maximum likelihood estimation of wdm model parameters
wdm <- function(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL, delta=NULL,
               xvar=NULL, start=NULL, fixed=0) {
  # save original function call
  cl <- match.call()

  # prepare passed arguments
  verifydata(data)
  if (is.numeric(data) & is.null(xvar)) data <- revampwiener(data, yvar=yvar)
  else if (length(yvar)==1) {
    cbind(revampwiener(data[,yvar]), data)
    yvar <- c("q", "resp")
  }
  fpar <- c("alpha"=unname(alpha), "tau"=unname(tau), 
    "beta"=unname(beta), "delta"=unname(delta))

  # estimate parameters
  if (!is.null(xvar)) {
    if(length(xvar)==1) {
      if(class(data[,xvar]) == "factor"){
        res <- list()
        res$coefficients <- fpar
        for (l in levels(data[,xvar])) {
          est <- mle(data[data[,xvar]==l,yvar], fpar, start)
          est$coefficients <- est$coefficients[!(names(est$coefficients) %in% names(fpar))]
          names(est$coefficients) <- paste(l,names(est$coefficients), sep=":")
          res$coefficients <- append(res$coefficients, est$coefficients)
          res$counts <- append(res$counts, est$counts)
          res$algorithm <- append(res$algorithm, list(est$algorithm))
          res$convergence <- append(res$convergence, est$convergence)
          res$message <- append(res$message, est$message)
          res$hessian <- append(res$hessian, list(est$hessian))
          res$loglik <- sum(res$loglik, est$loglik)
        }
      } else stop("xvar has to be a factor")
    }
  }
  else
    res <- mle(data[,yvar], fpar, start)

  # prepare return object
  res$nobs <- length(data[,1])
  res$npar <- length(res$coefficients)-fixed
  res$data <- data
  res$yvar <- yvar
  res$estpar <- c("alpha"=is.null(alpha), "tau"=is.null(tau),
                   "beta"=is.null(beta), "delta"=is.null(delta))
  res$call <- cl
  class(res) <- c("wdm")
  return(res)
}

### mle function
## internal function
mle <- function(data, fpar=NULL, start=NULL) {

  if (is.null(start))
  {
    start <- c(runif(1,1,2),min(data[,1])/3,runif(1,.2,.8),runif(1,-1,1))
  }

  start <- esvec(start, fpar)

  if (length(start)==0) {
    est <- list(
      coefficients = fpar,
      data = data,
      convergence = NULL,
      hessian = NULL,
      algorithm = list(type="None (all parameters fixed)") )
      est$value <- -logLik.wdm(est)
  }
  else if (length(start)==1) {
    ## only one parameter: use 'Brent (optim)'
    est <- optim(start,efn,data=data,fpar=fpar, method="Brent",
      lower=-100, upper=100, hessian=TRUE, control=list(maxit=2000))
    est$algorithm <- list(type="Brent (optim)", counts=est$counts, message=est$message)
    est$coefficients <- est$par
  }
  else
  {
    ## first: try 'BFGS (optim)'
    est <- tryCatch(optim(start,efn,data=data,fpar=fpar,
      method="BFGS",hessian=TRUE, control=list(maxit=2000)), 
      error=function(e) NULL)
    if (!is.null(est)) {
      if (est$convergence == 0) {
        est$algorithm <- list(type="BFGS (optim)", counts=est$counts, message=est$message)
        est$coefficients <- est$par
      }
      else est <- NULL
    }
    if (is.null(est)) {
      ## second: try 'Nelder-Mead (optim)', then 'BFGS (optim)'
      opt <- optim(start,efn,data=data,fpar=fpar, method="Nelder-Mead", control=list(maxit=2000))
      if (!is.null(opt)) {
        est <- tryCatch(optim(opt$par,efn,data=data,fpar=fpar,
          method="BFGS",hessian=TRUE, control=list(maxit=2000)), 
          error=function(e) NULL)
        if (!is.null(est) && est$convergence == 0) {
            est$algorithm <- list(type="BFGS (optim) after Nelder-Mead", counts=est$counts, message=est$message)
            est$coefficients <- est$par
          }
        }
      }
    if (is.null(est)) {
      ## third: try 'Newton-type (nlm)'
      ## note: suppressWarnings used for nlm, as nlm inflates warning messages
      est <- tryCatch(suppressWarnings(nlm(efn,start,data=data,fpar=fpar,hessian=TRUE)),
        error=function(e) NULL)
      if (!is.null(est)) {
        if  (est$code < 3) {
          est$convergence <- est$code
          est$coefficients <- est$estimate; est$estimate <- NULL
          est$value <- est$minimum; est$minimum <- NULL
          est$counts <- c(iterations=est$iterations); est$iterations <- NULL
          est$algorithm <- list(type="Newton-type (nlm)",
            gradient=est$gradient, counts=est$counts, message=est$message); est$gradient <- NULL
        }
        else est <- NULL
      }
    }
    if (is.null(est)) {
      ## fourth: try 'Nelder-Mead (optim)'
      est <- optim(start,efn,data=data,fpar=fpar, method="Nelder-Mead", control=list(maxit=2000))
      est$algorithm <- list(type="Nelder-Mead", counts=est$counts, message=est$message)
      est$coefficients <- est$par
    }
  }
  est$message <- NULL; est$counts <- NULL

  par <- eparvec(est$coefficients, fpar)

  res <- list(
    coefficients = par,
    loglik = -est$value,
    convergence = est$convergence,
    hessian = est$hessian,
    algorithm = est$algorithm
  )

  return(res)
}

## internal function
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

## internal function
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
  object$coefficients <- par

  res <- nlogLik.wdm(object)
  return(res)
}

## internal function
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

## additional functions

vcov.wdm <- function(object, ..., method="hessian") {
  # opg-estimator (outer product of gradients)
  if (method=="opg")
    res <- solve(crossprod(scorefun(object)))
  else if (method=="hessian")
    if(is.list(object$hessian)) {
      res <- list()
      for (k in 1:length(object$hessian)) {
        res[[k]] <- solve(object$hessian[[k]])
      }
    }
    else
      res <- solve(object$hessian)
  else 
    stop("Wrong method specified")

  return(res)
}
