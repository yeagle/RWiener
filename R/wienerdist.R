# internal function
verifypars <- function(alpha,tau,beta,delta)
{
  if(!is.numeric(alpha) || !is.numeric(tau) || 
     !is.numeric(beta) || !is.numeric(delta)) {
    return(FALSE)
  }
  if(alpha > 0 & 
     tau > 0 &
     beta >= 0 & beta <= 1) return(TRUE)
  else return(FALSE)
}

dwiener <- function(q, alpha,tau,beta,delta, resp="upper", give_log=FALSE) 
{
  if (!verifypars(alpha,tau,beta,delta) ||
      !is.numeric(q) || !(is.character(resp) || is.factor(resp))) {
    stop("bad parameter values")
  }

  if (!(length(resp) == length(q))) {
    stop("arguments q and resp need to be of the same length")
  }

  if(class(resp) == "factor") {
    resp <- as.character(resp)
  }

  d <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0")

    if (resp[i] == "upper") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta, give_log)
    else if (resp[i] == "lower") 
      d[i] <- .Call(dwiener_c, -q[i], alpha,tau,beta,delta, give_log)
    else if (resp[i] == "both") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta, give_log) +
           .Call(dwiener_c, -q[i], alpha,tau,beta,delta, give_log)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(d[i])) d[i] <- 0
  }

  return(d)
}

pwiener <- function(q, alpha,tau,beta,delta, resp="upper")
{
  if (!verifypars(alpha,tau,beta,delta) ||
      !is.numeric(q) || !(is.character(resp) || is.factor(resp))) {
    stop("bad parameter values")
  }

  if (!(length(resp) == length(q))) {
    stop("arguments q and resp need to be of the same length")
  }

  if(class(resp) == "factor") {
    resp <- as.character(resp)
  }

  p <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0")

    if (resp[i] == "upper") 
      p[i] <- .Call(pwiener_c, q[i], alpha,tau,beta,delta)
    else if (resp[i] == "lower")
      p[i] <- .Call(pwiener_c, -q[i], alpha,tau,beta,delta)
    else if (resp[i] == "both")
      p[i] <- .Call(pwiener_full_c, q[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(p[i])) p[i] <- 0
  }

  return(p)
}

qwiener <- function(p, alpha,tau,beta,delta, resp="upper")
{
  if (!verifypars(alpha,tau,beta,delta) ||
      !is.numeric(p) || !(is.character(resp) || is.factor(resp))) {
    stop("bad parameter values")
  }

  if (!(length(resp) == length(p))) {
    stop("arguments p and resp need to be of the same length")
  }

  if(class(resp) == "factor") {
    resp <- as.character(resp)
  }


  q <- vector("double", length=length(q))
  for (i in 1:length(p)) {
    if (p[i]<0) stop("p must be > 0")

    if (resp[i] == "upper")
      q[i] <- .Call(qwiener_c, p[i], alpha,tau,beta,delta)
    else if (resp[i] == "lower")
      q[i] <- .Call(qwiener_c, -p[i], alpha,tau,beta,delta)
    else if (resp[i] == "both")
      q[i] <- .Call(qwiener_full_c, p[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(q[i])) p[i] <- 0
  }

  return(q)
}

rwiener <- function(n, alpha,tau,beta,delta)
{
  if (!verifypars(alpha,tau,beta,delta)) {
    stop("bad parameter values")
  }

  res <- data.frame(q=vector("double"),resp=factor(levels=c("upper", "lower")))

  for (i in 1:n) {
    r <- .Call(rwiener_c, alpha,tau,beta,delta)
    if (r >= 0) res[i,] <- c(r,"upper")
    else res[i,] <- c(abs(r),"lower")
  }

  res[,1] <- as.double(res[,1])
  class(res) <- c("data.wiener", class(res))
  return(res)
}
