logLik.wiener <- function(object, ...) {
  data <- object$data
  x <- object$par

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
  else return(object$logLik)
}

# internal function
nlogLik.wiener <- function(object, ...) -logLik.wiener(object, ...)
  
deviance.wiener <- function(object, ...) {
  -2*logLik.wiener(object)
}

AIC.wiener <- function(object, ...) {
  if(is.null(object$loss)) {
    -2*logLik.wiener(object)+4*2 
  }
  else {
    data <- object$data
    x <- object$par
    object$loss(x,data)+length(x)*2 
  }
}

BIC.wiener <- function(object, ...) {
  if(is.null(object$loss)) {
    -2*logLik.wiener(object)+4*log(length(object$data[,1]))
  }
  else {
    data <- object$data
    x <- object$par
    if(is.list(data)) {
      loss(x,data)+length(x)*log(length(data[[1]][,1]))
    }
    else if (is.data.frame(data)) {
      loss(x,data)+length(x)*log(length(data[,1]))
    }
    else {
      stop("don't know how to handle the data object!")
    }
  }
}
