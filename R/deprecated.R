# document with all deprecated functions, kept for backwards compatibility

wiener_likelihood <- function(x, data) {
  #warning("wiener_likelihood is deprecated, use logLik.wiener instead")
  obj <- list(par=x,data=data)
  logLik.wiener(obj)
}

wiener_aic <- function(x, data, loss=NULL) {
  #warning("wiener_aic is deprecated, use AIC.wiener instead")
  obj <- list(par=x,data=data,loss=loss)
  AIC.wiener(obj)
}

wiener_bic <- function(x, data, loss=NULL) {
  #warning("wiener_bic is deprecated, use BIC.wiener instead")
  obj <- list(par=x,data=data,loss=loss)
  BIC.wiener(obj)
}

wiener_deviance <- function(x, data) {
  #warning("wiener_deviance is deprecated, use deviance.wiener instead")
  obj <- list(par=x,data=data)
  deviance.wiener(obj)
}

wiener_plot <- function(data) {
  #warning("wiener_plot is deprecated, use plot.data.wiener instead")
  plot.data.wiener(data)
}

