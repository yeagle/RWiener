check_wiener_pars <- function(alpha,tau,beta,delta)
{
  if(!is.real(alpha) || !is.real(tau) || 
     !is.real(beta) || !is.real(delta)) {
    return(FALSE)
  }
  if(alpha > 0 & 
     tau > 0 &
     beta >= 0 & beta <= 1) return(TRUE)
  else return(FALSE)
}

dwiener <- function(q, alpha,tau,beta,delta, resp="upper") 
{
 ##d <- .Call(dwiener_c)
 #print("asdf")
 ##out <- .C("dwiener", d=as.double(5), package="RWiener")
 #d <- .Call(dwiener_c)
 #print("asdf2")
 ##return(out$d)
 #return(d);

  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.real(q) || !is.character(resp)) {
    stop("bad parameter values!")
  }

  d <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0!")

    if (resp == "upper") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta)
    else if (resp == "lower") 
      d[i] <- .Call(dwiener_c, -q[i], alpha,tau,beta,delta)
    else if (resp == "both") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta) +
           .Call(dwiener_c, -q[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(d[i])) d[i] <- 0
  }

  return(d)
}

pwiener <- function(q, alpha,tau,beta,delta, resp="upper")
{
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.real(q) || !is.character(resp)) {
    stop("bad parameter values!")
  }

  p <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0!")

    if (resp == "upper") 
      p[i] <- .Call(pwiener_c, q[i], alpha,tau,beta,delta)
    else if (resp == "lower")
      p[i] <- .Call(pwiener_c, -q[i], alpha,tau,beta,delta)
    else if (resp == "both")
      p[i] <- .Call(pwiener_full_c, q[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(p[i])) p[i] <- 0
  }

  return(p)
}

qwiener <- function(p, alpha,tau,beta,delta, resp="upper")
{
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.real(p) || !is.character(resp)) {
    stop("bad parameter values!")
  }

  q <- vector("double", length=length(q))
  for (i in 1:length(p)) {
    if (p[i]<0) stop("p must be > 0!")

    if (resp == "upper")
      q[i] <- .Call(qwiener_c, p[i], alpha,tau,beta,delta)
    else if (resp == "lower")
      q[i] <- .Call(qwiener_c, -p[i], alpha,tau,beta,delta)
    else if (resp == "both")
      q[i] <- .Call(qwiener_full_c, p[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(q[i])) p[i] <- 0
  }

  return(q)
}

rwiener <- function(n, alpha,tau,beta,delta)
{
  if (!check_wiener_pars(alpha,tau,beta,delta)) {
    stop("bad parameter values!")
  }

  rdat <- data.frame(q=vector(),resp=factor(levels=c("upper", "lower")))

  for (i in 1:n) {
    r <- .Call(rwiener_c, alpha,tau,beta,delta)
    if (r >= 0) rdat[i,] <- c(r,"upper")
    else rdat[i,] <- c(abs(r),"lower")
    
  }

  return(rdat)
}
