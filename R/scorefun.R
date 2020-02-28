## estfun generic function
estfun <- function(x, ...) {
  UseMethod("estfun")
}

## define estfun for the wdm model object (needed by sctest function)
## same as scorefun, but aggregates by id (persons)
## empirical estimation function (score function) 
estfun.wdm <- function(x, ...) {
  res <- scorefun.wdm(x)
  if("id" %in% names(x$data)) {
    res <- cbind(res, id=x$data$id)
    res <- aggregate(. ~ id, sum, data=as.data.frame(res))[,-1]
  }
  return(res)
}

## empirical estimation function (score function) 
scorefun <- function(x, ...) {
  UseMethod("scorefun")
}

scorefun.wdm <- function(x, ...) {
  y <- x$data[,x$yvar]

  alpha <- x$coefficients["alpha"]
  tau <- x$coefficients["tau"]
  beta <- x$coefficients["beta"]
  delta <- x$coefficients["delta"]

  n <- length(y[,1])
  res <- matrix(rep(NA,4*n), n,4)
  colnames(res) <- c("alpha", "tau", "beta", "delta")
  for (i in 1:n) {
    if (y[i,2] == "lower") {
      res[i,1] <- sclalpha(y[i,1], alpha, tau, beta, delta)
      res[i,2] <- scltau(y[i,1], alpha, tau, beta, delta)
      res[i,3] <- sclbeta(y[i,1], alpha, tau, beta, delta)
      res[i,4] <- scldelta(y[i,1], alpha, tau, beta, delta)
    }
    else if (y[i,2] == "upper") {
      res[i,1] <- sclalpha(y[i,1], alpha, tau, 1-beta, -delta)
      res[i,2] <- scltau(y[i,1], alpha, tau, 1-beta, -delta)
      res[i,3] <- sclbeta(y[i,1], alpha, tau, 1-beta, -delta)
      res[i,4] <- scldelta(y[i,1], alpha, tau, 1-beta, -delta)
    }
  }

  return(res)
}

pow <- function(x,y) x^y

kappaLT <- function(t, err=1e-10) {
  (sqrt(2)*sqrt(-log(pi*err*t)/t)/pi)
}

kappaST <- function(t, err=1e-10) {
  (sqrt(2)*sqrt(-t*log(2*sqrt(2)*sqrt(pi)*err*sqrt(t))) + 2)
}

difflogdl01alphaST <- function(t, alpha, beta, kappa) {
  res <- 0; res1<-0;res2<-0;res3<-0;res4<-0
  for (k in -ceiling((kappa-1)/2):floor((kappa-1)/2)) {
    res1 <- res1 + (-alpha*pow(beta, 3)*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t - 4*alpha*pow(beta, 2)*k*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t - 4*alpha*beta*pow(k, 2)*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t )
    res2 <- res2 + (-2*alpha*pow(beta, 2)*k*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t - 8*alpha*beta*pow(k, 2)*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t - 8*alpha*pow(k, 3)*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t)
    res3 <- res3 + (beta*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t))
    res4 <- res4 + (2*k*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t))
  }
  res <- (res1+res2)/(res3+res4) + 3/alpha
  return(res)
}

difflogdl01alphaLT <- function(t, alpha, beta, kappa) {
  res <- 0; res1<-0;res2<-0
  for (k in 1:ceiling(kappa)) {
    res1 <- res1 + (pow(pi, 2)*pow(k, 3)*t*exp(-0.5*pow(pi, 2)*pow(k, 2)*t/pow(alpha, 2))*sin(pi*beta*k)/pow(alpha, 3))
    res2 <- res2 + (k*exp(-0.5*pow(pi, 2)*pow(k, 2)*t/pow(alpha, 2))*sin(pi*beta*k))
  }
  res <- res1 / res2
  return(res)
}

difflogdl01alpha <- function(t, alpha, beta) {
  kst <- kappaST(t)
  klt <- kappaLT(t)
  wlam <- kst - klt
  if(wlam < 0) difflogdl01alphaST(t, alpha, beta, kst)
  else difflogdl01alphaLT(t, alpha, beta, klt)
}

sclalpha <- function(t, alpha, tau, beta, delta) {
  t <- t-tau
  res <- -beta*delta + difflogdl01alpha(t,alpha,beta) - 2/alpha
  return(res)
}

difflogdl01tauST <- function(t, tau, alpha, beta, kappa) {
  res <- 0; res1<-0;res2<-0;res3<-0;res4<-0
  for (k in -ceiling((kappa-1)/2):floor((kappa-1)/2)) {
    res1 <- res1 + (-2*pow(beta, 3)*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau))/pow(2*t - 2*tau, 2) - 8*pow(beta, 2)*k*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau))/pow(2*t - 2*tau, 2) - 8*beta*pow(k, 2)*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau))/pow(2*t - 2*tau, 2))
    res2 <- res2 + (-4*pow(beta, 2)*k*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau))/pow(2*t - 2*tau, 2) - 16*beta*pow(k, 2)*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau))/pow(2*t - 2*tau, 2) - 16*pow(k, 3)*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau))/pow(2*t - 2*tau, 2))
    res3 <- res3 + (beta*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau)))
    res4 <- res4 + (2*k*exp(-pow(beta, 2)/(2*t - 2*tau))*exp(-4*pow(k, 2)/(2*t - 2*tau))*exp(-4*beta*k/(2*t - 2*tau)))
  }
  res <- (res1+res2)/(res3+res4) + (1.5)/(t - tau)
  return(res)
}

difflogdl01tauLT <- function(t, tau, alpha, beta, kappa) {
  res <- 0; res1<-0;res2<-0
  for (k in 1:ceiling(kappa)) {
    res1 <- res1 + ((0.5)*pow(pi, 2)*pow(k, 3)*exp(-0.5*pow(pi, 2)*pow(k, 2)*t)*exp((0.5)*pow(pi, 2)*pow(k, 2)*tau)*sin(pi*beta*k))
    res2 <- res2 + (k*exp(-0.5*pow(pi, 2)*pow(k, 2)*t)*exp((0.5)*pow(pi, 2)*pow(k, 2)*tau)*sin(pi*beta*k))
  }
  res <- res1 / res2
  return(res)
}

difflogdl01tau <- function(t, tau, alpha, beta) {
  kst <- kappaST(t-tau)
  klt <- kappaLT(t-tau)
  wlam <- kst - klt
  if(wlam < 0) difflogdl01tauST(t, tau, alpha, beta, kst)
  else difflogdl01tauLT(t, tau, alpha, beta, klt)
}
scltau <- function(t, alpha, tau, beta, delta) {
  res <- delta^2/2 - 1/alpha^2 * difflogdl01tau(t, tau, alpha, beta)
  #res <- 0
  return(res)
}

difflogdl01betaST <- function(t, alpha, beta, kappa) {
  res <- 0; res1<-0;res2<-0;res3<-0;res4<-0
  for (k in -ceiling((kappa-1)/2):floor((kappa-1)/2)) {
    res1 <- res1 + (-2*pow(alpha, 2)*beta*k*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t - 4*pow(alpha, 2)*pow(k, 2)*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t)
    res2 <- res2 + (-pow(alpha, 2)*pow(beta, 2)*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t - 2*pow(alpha, 2)*beta*k*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t)/t + exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t))
    res3 <- res3 + (beta*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t))
    res4 <- res4 + (2*k*exp(-0.5*pow(alpha, 2)*pow(beta, 2)/t)*exp(-2*pow(alpha, 2)*pow(k, 2)/t)*exp(-2*pow(alpha, 2)*beta*k/t))

  }
  res <- (res1+res2)/(res3+res4)
  return(res)
}

difflogdl01betaLT <- function(t, alpha, beta, kappa) {
  res <- 0; res1<-0;res2<-0
  for (k in 1:ceiling(kappa)) {
    res1 <- res1 + (pi*pow(k, 2)*exp(-0.5*pow(pi, 2)*pow(k, 2)*t/pow(alpha, 2))*cos(pi*beta*k))
    res2 <- res2 + (k*exp(-0.5*pow(pi, 2)*pow(k, 2)*t/pow(alpha, 2))*sin(pi*beta*k))
  }
  res <- res1 / res2
  return(res)
}

difflogdl01beta <- function(t, alpha, beta) {
  kst <- kappaST(t)
  klt <- kappaLT(t)
  wlam <- kst - klt
  if(wlam < 0) difflogdl01betaST(t, alpha, beta, kst)
  else difflogdl01betaLT(t, alpha, beta, klt)
}

sclbeta <- function(t, alpha, tau, beta, delta) {
  t <- t-tau
  res <- -alpha*delta + difflogdl01beta(t, alpha, beta)
  return(res)
}

scldelta <- function(t, alpha, tau, beta, delta) {
  t <- t-tau
  res <- -alpha*beta - delta*t
  return(res)
}
