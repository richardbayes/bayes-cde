library(DPpackage)
load('output/section4_2_DPcdensity.RData')

# MCMC Summary
summary(out)
plot(out)

#jpeg('figs/DPcdensity.jpg',
#     res=300,width=6,height=6,units='in')
setEPS()
postscript('figs/DPcdensity.eps', width=8, height=3)
  par(mfrow=c(1,4))
  
  # .76,.76
  ind <- 3
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="",
    ylab="Density",
    ylim=c(0,.6),
    main=expression('x'[1]*'=.76, x'[2]*'=.76'),
    las=1
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  curve(dnorm(x,1,sqrt(.5)),add=TRUE,lwd=2)
  
  # .9, .9
  ind <- 10 
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="",
    ylab="",
    ylim=c(0,.6),
    las=1,
    main=expression('x'[1]*'=.9, x'[2]*'=.9')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  curve(dnorm(x,1,sqrt(.5)),lwd=2,add=TRUE)
  
  # .1, .8
  ind <- 8 
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="",
    ylab="",
    las=1,
    main=expression('x'[1]*'=.1, x'[2]*'=.8')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  newdens <- function(x) .5*dnorm(x,1,1) + .5*dnorm(x,5,1)
  curve(newdens(x),lwd=2,add=TRUE)
  
  # .8, .1
  ind <- 7
  plot(
    out$grid,
    out$densp.h[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="",
    ylab="",
    las=1,
    main=expression('x'[1]*'=.8, x'[2]*'=.1')
  )
  
  lines(out$grid,out$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densp.m[ind,],lwd=2,type="l",lty=2)
  curve(dgamma(x,10,2),lwd=2,add=TRUE)
  
  mtext("Y", 1, outer = TRUE, line = -2)
dev.off()