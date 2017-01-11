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
        models <- as.character(c(cl[[i]], cl[[i+1]]))

        if (i == 2) {
          res <- list(
          models = models,
          modeldf = c(wdmspecific$npar, wdmgeneral$npar),
          model.AIC = c(AIC(wdmspecific), AIC(wdmgeneral)),
          model.BIC = c(BIC(wdmspecific), BIC(wdmgeneral)),
          model.loglik = c(logLik(wdmspecific), logLik(wdmgeneral)),
          G2 = c(NA, G2), 
          Df = c(NA, Df),
          pvalue = c(NA, pvalue)
          )
        }
        else {
          res$modeldf <- c(res$modeldf, wdmgeneral$npar)
          res$models <- c(res$models, models[2])
          res$model.AIC <- c(res$model.AIC, AIC(wdmgeneral))
          res$model.BIC <- c(res$model.BIC, BIC(wdmgeneral))
          res$model.loglik <- c(res$model.loglik, logLik(wdmgeneral))
          res$G2 <- c(res$G2, G2)
          res$Df <- c(res$Df, Df)
          res$pvalue <- c(res$pvalue, pvalue)
        }
      }
    } # end for loop
  } # end if(test==LRT)
  else stop("Only LRT Test implemented yet")

  out <- data.frame(df = res$modeldf, 
                    AIC = res$model.AIC, BIC = res$model.BIC, logLik = res$model.loglik,
                    Df = res$Df, LRT.G2 = res$G2,
                    pvalue = res$pvalue)
  dimnames(out) <- list(1:length(res$models), c("Model.df", "AIC", "BIC", "logLik",
                            "LRT.df", "LRT.G2", "p.value"))
  structure(out,
              heading = c("Model comparison Table with LRTs\n",
                          paste0("Model ", 1:length(res$models), ": ",
                                 res$models, collapse = "\n"),""),
  class = c("anova", "data.frame"))
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
  colnames(wtab) <- c("Test-Statistic", "p.value")
  rownames(wtab) <- c("W2","W") 

  cat("\n")
  cat("Wiener Wald Test")
  cat("\n\n")
  cat("Null hypothesis: ", x$test$theta, "=", x$test$theta0 , sep="")
  cat("\n")
  #cat("W2 ~ Chisq, Df=1; p-value: upper tail (one-tailed)\n")
  cat("W ~ Normal Distribution\n")
  cat("\n")
  print(wtab[2,])
}
