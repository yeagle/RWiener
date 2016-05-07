# internal function
verifydata <- function(data) {
  if (is.null(data)) stop("Error: missing values (no data supplied)")
  if (!is.wiener(data)) {
    if (!(is.data.frame(data) | is.numeric(data))) 
      stop("Error: supplied data in wrong format")
  }
}

is.wiener <- function(x) {
  res <- (inherits(x, "data.wiener") | inherits(x, "numdata.wiener"))
  return(res)
}

as.wiener <- function(data, yvar=c("q", "resp")) {
  if(is.data.frame(data) & (sum(as.numeric(colnames(data) == yvar))==2) )
  {
    class(data) <- c("data.wiener", "data.frame")
  }
  else if(is.numeric(data) | is.vector(data))
  {
    class(data) <- c("numdata.wiener", "numeric")
  }
  else stop("can only convert vectors (with + / - values for upper/lower bound) or data.frames (with 2 columns: 'q' and 'resp').")
  return(data)
}

reshape.wiener <- function(data, yvar=c("q", "resp"), direction="auto") {
  verifydata(data)

  if(is.data.frame(data) & (direction %in% c("wide", "auto")))
  {
    res <- data[,yvar[1]]
    for (i in 1:(length(data[,1])))
    {
      if(data[i,yvar[2]] == "upper") res[i] <- data[i,yvar[1]]
      else res[i] <- -data[i,yvar[1]]
    }
    class(res) <- c("numdata.wiener", "numeric")
  }
  else if ((is.vector(data) | is.numeric(data)) 
           & (direction %in% c("long", "auto"))) {
    res <- data.frame(as.numeric(abs(data)), factor((data>0), levels=c("TRUE", "FALSE"),
             labels=c("upper", "lower")))
    colnames(res) <- yvar[1:2]
    class(res) <- c("data.wiener", "data.frame")
  }
  else stop("Error: argument(s) not valid")

  return(res)
}

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
          res$logLik <- sum(res$logLik, est$logLik)
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
