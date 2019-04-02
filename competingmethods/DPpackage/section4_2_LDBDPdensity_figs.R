library(DPpackage)
load('output/section4_2_LDBDPdensity.RData')

# MCMC summary
summary(out)
plot(out)

# jpeg('figs/LDBDPdensity.jpg',
#      res=300,width=6,height=6,units='in')
setEPS()
postscript('figs/LDBDPdensity.eps', width=6, height=6)
  par(mfrow=c(2,2))
  y <- alldat[, 1]
  
  # .76,.76
  ind <- 1
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    ylim=c(0,8),
    las=1,
    main=expression('x'[1]*'=.76, x'[2]*'=.76')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  themean <- (1 - min(y)) / (max(y) - min(y))
  thesd <- sqrt(.5)/(max(y) - min(y))
  curve(dnorm(x,themean,thesd),lwd=2,add=TRUE)
  
  # .9, .9
  ind <- 6
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    las=1,
    ylim=c(0,8),
    main=expression('x'[1]*'=.9, x'[2]*'=.9')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  themean <- (1 - min(y)) / (max(y) - min(y))
  thesd <- sqrt(.5)/(max(y) - min(y))
  curve(dnorm(x,themean,thesd),lwd=2,add=TRUE)
  
  # .1, .8
  ind <- 5 
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    las=1,
    main=expression('x'[1]*'=.1, x'[2]*'=.8')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  themean1 <- (1 - min(y)) / (max(y) - min(y))
  themean2 <- (5 - min(y)) / (max(y) - min(y))
  thesd <- 1/(max(y) - min(y))
  newdens <- function(x) .5*dnorm(x,themean1,thesd) + .5*dnorm(x,themean2,thesd)
  curve(newdens(x),lwd=2,add=TRUE)
  
  # .8, .1
  ind <- 4
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    ylim=c(0,3.5),
    las=1,
    main=expression('x'[1]*'=.8, x'[2]*'=.1')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  therange <- (max(y) - min(y))
  themean <- min(y)
  newgamma <- function(z) therange*dgamma(z*therange + themean,10,2)
  #curve(dgamma(x,10,2),col='red',add=TRUE)
  curve(newgamma(x),lwd=2,add=TRUE)
dev.off()