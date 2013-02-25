check_wiener_pars <- function(alpha,tau,beta,delta)
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
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.numeric(q) || !is.character(resp)) {
    stop("bad parameter values!")
  }

  d <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0!")

    if (resp == "upper") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta, give_log)
    else if (resp == "lower") 
      d[i] <- .Call(dwiener_c, -q[i], alpha,tau,beta,delta, give_log)
    else if (resp == "both") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta, give_log) +
           .Call(dwiener_c, -q[i], alpha,tau,beta,delta, give_log)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(d[i])) d[i] <- 0
  }

  return(d)
}

pwiener <- function(q, alpha,tau,beta,delta, resp="upper")
{
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.numeric(q) || !is.character(resp)) {
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
      !is.numeric(p) || !is.character(resp)) {
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

  rdat <- data.frame(q=vector("double"),resp=factor(levels=c("upper", "lower")))

  for (i in 1:n) {
    r <- .Call(rwiener_c, alpha,tau,beta,delta)
    if (r >= 0) rdat[i,] <- c(r,"upper")
    else rdat[i,] <- c(abs(r),"lower")
    
  }

  rdat[,1] < as.real(rdat[,1])
  return(rdat)
}

wiener_likelihood <- function(x, dat) {
 if (!check_wiener_pars(x[1],x[2],x[3],x[4])) {
    return(-Inf)
  }
  ll <- vector("double", length(dat[,1]))
  for (i in 1:length(dat[,1])) {
    ll[i] <- dwiener(as.real(dat[i,1]), x[1],x[2],x[3],x[4], 
                  resp=as.character(dat[i,2]), give_log=TRUE)
  }
  return(sum(ll))
}
  
wiener_deviance <- function(x, dat) {
  -2*wiener_likelihood(x,dat)
}
wiener_bic <- function(x, dat) {
  -2*wiener_likelihood(x,dat)+4*log(length(dat[,1]))
}
wiener_aic <- function(x, dat) {
  -2*wiener_likelihood(x,dat)+4*2 
}

