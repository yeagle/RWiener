is.wiener <- function(x) inherits(x, "dat.wiener")

as.wiener <- function(dat)
{
  if(is.data.frame(dat) & (sum(as.numeric(colnames(dat) == c("q", "resp")))==2) )
  {
    class(dat) <- c("dat.wiener", class(dat))
  }
  else if(is.numeric(dat) | is.vector(dat))
  {
    class(dat) <- c("dat.wiener", class(dat))
  }
  else stop("can only convert vectors (with + / - values for upper/lower bound) or data.frames (with 2 columns: 'q' and 'resp').")
  return(dat)
}

reshape.wiener <- function(dat)
{
  if(is.null(dat)) stop("missing values (no data supplied)!")
  if(!is.wiener(dat)) stop("supplied data is not of class dat.wiener!")

  if(is.data.frame(dat))
  {
    rval <- dat$q
    for (i in 1:(length(dat[,1])))
    {
      if(dat[i,]$resp == "upper") rval[i] <- dat[i,]$q
      else rval[i] <- -dat[i,]$q
    }
  }
  else
    rval <- data.frame(q=abs(dat),resp=factor((dat>0), levels=c("TRUE", "FALSE"),
             labels=c("upper", "lower")))

  rval <- as.wiener(rval)
  return(rval)
}
