is.wiener <- function(x) inherits(x, "data.wiener")

as.wiener <- function(data) {
  if(is.data.frame(data) & (sum(as.numeric(colnames(data) == c("q", "resp")))==2) )
  {
    class(data) <- c("data.wiener", class(data))
  }
  else if(is.numeric(data) | is.vector(data))
  {
    class(data) <- c("data.wiener", class(data))
  }
  else stop("can only convert vectors (with + / - values for upper/lower bound) or data.frames (with 2 columns: 'q' and 'resp').")
  return(data)
}

reshape.wiener <- function(data, direction="auto") {
  if(is.null(data)) stop("missing values (no data supplied)!")
  if(!is.wiener(data)) stop("supplied data is not of class data.wiener!")

  if(is.data.frame(data) & (direction %in% c("wide", "auto")))
  {
    rval <- data$q
    for (i in 1:(length(data[,1])))
    {
      if(data[i,]$resp == "upper") rval[i] <- data[i,]$q
      else rval[i] <- -data[i,]$q
    }
  }
  else if ((is.vector(data) | is.numeric(data)) 
           & (direction %in% c("long", "auto")))
    rval <- data.frame(q=abs(data),resp=factor((data>0), levels=c("TRUE", "FALSE"),
             labels=c("upper", "lower")))
  else stop("Error: argument(s) not valid")

  class(rval) <- c("data.wiener", class(rval))
  return(rval)
}

wm <- function(data, alpha=NULL, tau=NULL, beta=NULL, delta=NULL,
               start=NULL) {
  rval <- estfun(data, alpha, tau, beta, delta, start)
  rval$logLik <- -rval$value
  rval$value <- NULL
  rval$data <- data
  rval$estpar <- c("alpha"=is.null(alpha), "tau"=is.null(tau),
                   "beta"=is.null(beta), "delta"=is.null(delta))
  class(rval) <- c("wiener")
  return(rval)
}
