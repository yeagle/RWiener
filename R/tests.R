anova.wdm <- function(object, ..., test="LRT") {
  cl <- match.call()

  if(test=="LRT") {
    for (i in 2:(length(cl)-1)) {
      if(!is.object(eval(cl[[i+1]]))) break
      else {
        wdmspecific <- eval(cl[[i]])
        wdmgeneral <- eval(cl[[i+1]])
        
        G2 <- 2 * ( logLik(wdmgeneral)-logLik(wdmspecific) )
        Df <- wdmgeneral$npar - wdmspecific$npar
        pvalue <- pchisq(G2, Df, lower.tail=FALSE)

        pvalue2 <- 1
        if (pvalue <= 0.10) pvalue2 <- 0.05
        if (pvalue <= 0.05) pvalue2 <- 0.05
        if (pvalue <= 0.01) pvalue2 <- 0.01
        if (pvalue <= 0.001) pvalue2 <- 0.001

        wtab <- rbind(as.double(c(wdmspecific$npar, 
                        round(c(AIC(wdmspecific), BIC(wdmspecific), logLik(wdmspecific)), 4),
                        " ", " ", " ")),
                      as.double(c(wdmgeneral$npar, 
                        round(c(AIC(wdmgeneral), BIC(wdmgeneral), logLik(wdmgeneral)), 4),
                        round(c(G2, Df, pvalue2), 4))) )
        colnames(wtab) <- c("df", "AIC", "BIC", "logLik",
                            "LRT.G2", "LRT.df", "p-value <=")
        rownames(wtab) <- c(1,2) 
        models <- as.character(c(cl[[i]], cl[[i+1]]))

        wlrts <- list(
        G2 = G2, 
        Df = Df,
        pvalue = pvalue,
        models = models,
        table = wtab
        )

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
    class(res) <- c(class(res),"wlrt", "anova")
    return(res)
  } # end if(test==LRT)
  else stop("Only LRT Test implemented yet")
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

## define waldtest function to be generic
waldtest <- function(object, ...) UseMethod("waldtest")

## H1: wdm
## waldtest.wdm function
waldtest.wdm <- function(object, ..., theta="delta", theta0=0) {
  pars <- coef(object)
  vars <- diag(vcov(object)) 
  names(vars) <- names(pars)
  n <- nobs(object)

  # method 1: chi-squared distribution with df=1
  W2 <- ((pars[theta]-theta0)^2 / vars[theta])
  chisq <- pchisq(W2, 1, lower.tail=FALSE)
  # method 2: normal distribution
  W <- (pars[theta]-theta0)/(sqrt(vars[theta]))
  ND <- 2 * pnorm(abs(W), lower.tail=FALSE)

  res <- list(
  W2 = W2,
  W2.pvalue = chisq,
  W = W,
  pvalue = ND,
  test = list(theta=theta,theta0=theta0)
  )

  class(res) <- c(class(res), "wwaldt", "waldtest")
  return(res)
}

print.wwaldt <- function(x, ...) {
  wtab <- rbind(c(x$W2, x$W2.pvalue), 
                c(x$W, x$pvalue) )
  colnames(wtab) <- c("Test-Statistic", "pvalue")
  rownames(wtab) <- c("W2","W") 

  cat("\n")
  cat("Wiener Wald Test")
  cat("\n\n")
  cat("Null hypothesis: ", x$test$theta, "=", x$test$theta0 , sep="")
  cat("\n")
  #cat("W2 ~ Chisq, Df=1; p-value: upper tail (one-tailed)\n")
  cat("W ~ Normal Distribution; p-value: upper tail (two-tailed)\n")
  cat("\n")
  print(wtab[2,])
}
