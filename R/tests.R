anova.wdm <- function(object, ..., test="LRT") {
  cl <- match.call()

  if(test=="LRT") {
    for (i in 2:(length(cl)-1)) {
      if(!is.object(eval(cl[[i+1]]))) break
      else {
        wlrts <- eval(call("wlrt", cl[[i]], cl[[i+1]]))
        if (i == 2) res <- wlrts
        else {
          res$table <- rbind(res$table, wlrts$table[2,])
          res$models <- c(res$models, wlrts$models[2])
          res$G2 <- c(res$G2, wlrts$G2)
          res$Df <- c(res$Df, wlrts$Df)
          res$pvalue <- c(res$pvalue, wlrts$pvalue)
        }
      }
    } # end for loop

    rownames(res$table) <- 1:length(res$table[,1])
    res$call <- cl
    class(res) <- c(class(res), "anova")
    return(res)
  } # end if(test==LRT)
  else stop("Only LRT Test implemented yet")
}

## internal function
wlrt <- function(wdmspecific, wdmgeneral) {

  cl <- match.call()

  G2 <- 2 * ( logLik(wdmgeneral)-logLik(wdmspecific) )
  Df <- wdmgeneral$npar - wdmspecific$npar
  pvalue <- pchisq(G2, Df, lower.tail=FALSE)

  pvalue2 <- 1
  if (pvalue <= 0.10) pvalue2 <- 0.05
  if (pvalue <= 0.05) pvalue2 <- 0.05
  if (pvalue <= 0.01) pvalue2 <- 0.01
  if (pvalue <= 0.001) pvalue2 <- 0.001

  wtab <- rbind(as.double(c(length(coef(wdmspecific)), 
                  round(c(AIC(wdmspecific), BIC(wdmspecific), logLik(wdmspecific)), 4),
                  " ", " ", " ")),
                as.double(c(length(coef(wdmgeneral)), 
                  round(c(AIC(wdmgeneral), BIC(wdmgeneral), logLik(wdmgeneral)), 4),
                  round(c(G2, Df, pvalue2), 4))) )
  colnames(wtab) <- c("df", "AIC", "BIC", "logLik",
                      "LRT.G2", "LRT.df", "p-value <=")
  rownames(wtab) <- c(1,2) 
  models <- as.character(c(cl[[2]], cl[[3]]))

  res <- list(
  G2 = G2, 
  Df = Df,
  pvalue = pvalue,
  models = models,
  table = wtab
  )

  class(res) <- c(class(res), "wlrt")
  return(res)
}

print.wlrt <- function(x, ...) {
  cat("\n")
  for (i in 1:length(x$models)) {
    cat("Model ", i, ": ", x$models[i], sep="")
    cat("\n")
  }
  cat("\n")
  print(x$table)
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
