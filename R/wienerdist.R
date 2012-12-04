check_wiener_pars <- function(alpha,tau,beta,delta)
{
  if(alpha > 0 & 
     tau > 0 &
     beta >= 0 & beta <= 1) return(TRUE)
  else return(FALSE)
}

dwiener <- function(q, alpha,tau,beta,delta, resp="upper") 
{
  if (q<0) stop("q must be > 0!")
  if (!check_wiener_pars(alpha,tau,beta,delta)) stop("bad parameter
                                                      values")

  if (resp == "upper") 
    d <- .Call("dwiener", q, alpha,tau,beta,delta, PACKAGE="RWiener")
  else if (resp == "lower") 
    d <- .Call("dwiener", -q, alpha,tau,beta,delta, PACKAGE="RWiener")
  else if (resp == "both") 
    d <- .Call("dwiener", q, alpha,tau,beta,delta, PACKAGE="RWiener") +
         .Call("dwiener", -q, alpha,tau,beta,delta, PACKAGE="RWiener")
  else stop("resp must be either 'lower', 'upper' or 'both'")

  return(d)
}

pwiener <- function(q, alpha,tau,beta,delta, resp="upper")
{
  if (q<0) stop("q must be > 0!")
  if (!check_wiener_pars(alpha,tau,beta,delta)) stop("bad parameter
                                                      values")

  if (resp == "upper") 
    q <- .Call("pwiener", q, alpha,tau,beta,delta, PACKAGE="RWiener")
  else if (resp == "lower")
    q <- .Call("pwiener", -q, alpha,tau,beta,delta, PACKAGE="RWiener")
  else if (resp == "both")
    q <- .Call("pwiener_full", q, alpha,tau,beta,delta, PACKAGE="RWiener")
  else stop("resp must be either 'lower', 'upper' or 'both'")

  return(q)
}

qwiener <- function(p, alpha,tau,beta,delta, resp="upper")
{
  if (q<0) stop("q must be > 0!")
  if (!check_wiener_pars(alpha,tau,beta,delta)) stop("bad parameter
                                                      values")

  if (resp == "upper")
    p <- .Call("qwiener", p, alpha,tau,beta,delta, PACKAGE="RWiener")
  else if (resp == "lower")
    p <- .Call("qwiener", -p, alpha,tau,beta,delta, PACKAGE="RWiener")
  else if (resp == "both")
    p <- .Call("qwiener_full", p, alpha,tau,beta,delta, PACKAGE="RWiener")
  else stop("resp must be either 'lower', 'upper' or 'both'")

  return(p)
}

rwiener <- function(n, alpha,tau,beta,delta)
{
  if (!check_wiener_pars(alpha,tau,beta,delta)) stop("bad parameter
                                                      values")

  rdat <- data.frame(q=vector(),resp=factor(levels=c("upper", "lower")))

  for (i in 1:n) {
    r <- c(r, .Call("rwiener", alpha,tau,beta,delta, PACKAGE="RWiener"))
    if (r >= 0) rdat[i,] <- c(r,"upper")
    else rdat[i,] <- c(abs(r),"lower")
    
  }

  return(rdat)
}
