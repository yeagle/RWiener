logLik.wdm <- function(object, ...) {
  data <- object$data
  x <- object$coefficients

  if (length(x) == 4) { 
    if (!verifypars(x[1],x[2],x[3],x[4])) {
      return(-Inf)
    }
    ll <- vector("double", length(data[,1]))
    for (i in 1:length(data[,1])) {
      ll[i] <- dwiener(as.double(data[i,1]), x[1],x[2],x[3],x[4], 
                    resp=as.character(data[i,2]), give_log=TRUE)
    }
    return(sum(ll))
  }
  else return(object$loglik)
}

## internal function
nlogLik.wdm <- function(object, ...) -logLik.wdm(object, ...)
  
deviance.wdm <- function(object, ...) {
  -2*logLik.wdm(object)
}

AIC.wdm <- function(object, ...) {
  if(is.null(object$loss)) {
    -2*logLik.wdm(object)+4*2 
  }
  else {
    data <- object$data
    x <- object$par
    object$loss(x,data)+length(x)*2 
  }
}

BIC.wdm <- function(object, ...) {
  if(is.null(object$loss)) {
    -2*logLik.wdm(object)+4*log(length(object$data[,1]))
  }
  else {
    data <- object$data
    x <- object$par
    if(is.list(data)) {
      object$loss(x,data)+length(x)*log(length(data[[1]][,1]))
    }
    else if (is.data.frame(data)) {
      object$loss(x,data)+length(x)*log(length(data[,1]))
    }
    else {
      stop("don't know how to handle the data object")
    }
  }
}
