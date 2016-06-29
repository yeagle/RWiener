wlrt <- function(wdmspecific, wdmgeneral) {

  G2 <- 2 * ( logLik(wdmgeneral)-logLik(wdmspecific) )
  Df <- wdmgeneral$npar - wdmspecific$npar
  pvalue <- pchisq(G2, Df, lower.tail=FALSE)

  res <- list(
  G2 = G2, 
  Df = Df,
  pvalue = pvalue
  )

  return(res)
}

wscoret <- function(wdmh0) {
  sc <- scorefun(wdmh0)
  info <- infofun(wdmh0)
  n <- wdmh0$nobs

  W <- 1/n * sum( sc^2/info ) 

  Df <- 1
  pvalue <- pchisq(W, Df, lower.tail=FALSE)

  res <- list(
  W = W, 
  Df = Df,
  pvalue = pvalue
  )

  return(res)
}

wwaldt <- function(wdmh0, wdmh1) {
  res <- NULL
  return(res)
}
