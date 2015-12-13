#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# plot.R
#
# (c) 2015 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2015-12-12
# last mod 2015-12-13 14:38 DW
#

# Plot function by Rainer W. Alexandrowicz
plot.wiener <- function(x, ...)  {
  rt = as.double(x$q)                  # response time
  rc = as.numeric(x$resp)              # response cat: 1=up 2=lo
  dpos = try(density(rt[rc==1],from=0))  # density upper
  dneg = try(density(rt[rc==2],from=0))  # density lower
  maxt = max(pretty(max(rt)))            # overall max response time
  maxd = max(dpos$y,dneg$y)              # overall max density

  par(mar=c(0,5,0,0),mfcol=c(2,1),ask=FALSE)

  plot(dpos,xlim=c(0,maxt),ylim=c(0,maxd),las=2,lwd=2,col="green3",
      main="",ylab="",ask=FALSE)
      rug(rt[rc== 1],col="green3")
      mtext("Density of positive responses",side=2,line=4,cex=0.8)
  plot(dneg,xlim=c(0,maxt),ylim=c(maxd,0),las=2,lwd=2,col="red",
      main="",ylab="",ask=FALSE)
      mtext("Density of negative responses",side=2,line=4,cex=0.8)
      rug(rt[rc==2],col="red",side=3)
}

wiener_plot <- function(dat) {
  warning("wiener_plot is deprecated, use plot.wiener instead!")
  plot.wiener(dat)
}
