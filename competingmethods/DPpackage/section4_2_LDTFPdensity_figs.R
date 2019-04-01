library(DPpackage)
load('output/section4_2_LDTFPdensity.RData')

# MCMC results
summary(out)
plot(out)

# jpeg('figs/LDTFPdensity.jpg',
#      res=300,width=6,height=6,units='in')
setEPS()
postscript('figs/LDTFPdensityclean.eps', width=6, height=6)
  par(mfrow=c(2,2))
  # .76,.76
  ind <- 1
  plot(
    out$grid,
    out$densu[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    ylim=c(0,.8),
    las=1,
    main=expression('Posterior Density at x'[1]*'=.76, x'[2]*'=.76')
  )
  
  lines(out$grid,out$densl[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densm[ind,],lwd=2,type="l",lty=2)
  curve(dnorm(x,1,sqrt(.5)),lwd=2,add=TRUE)
  
  # .9, .9
  ind <- 6
  plot(
    out$grid,
    out$densu[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    las=1,
    ylim=c(0,.8),
    main=expression('Posterior Density at x'[1]*'=.9, x'[2]*'=.9')
  )
  
  lines(out$grid,out$densl[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densm[ind,],lwd=2,type="l",lty=2)
  curve(dnorm(x,1,sqrt(.5)),lwd=2,add=TRUE)
  
  # .1, .8
  ind <- 5
  plot(
    out$grid,
    out$densu[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    las=1,
    main=expression('Posterior Density at x'[1]*'=.1, x'[2]*'=.8')
  )
  
  lines(out$grid,out$densl[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densm[ind,],lwd=2,type="l",lty=2)
  newdens <- function(x) .5*dnorm(x,1,1) + .5*dnorm(x,5,1)
  curve(newdens(x),lwd=2,add=TRUE)
  
  # .8, .1
  ind <- 4
  plot(
    out$grid,
    out$densu[ind,],
    lwd=1,
    type="l",
    lty=2,
    xlab="Y",
    ylab="Density",
    ylim=c(0,.3),
    las=1,
    main=expression('Posterior Density at x'[1]*'=.8, x'[2]*'=.1')
  )
  lines(out$grid,out$densl[ind,],lwd=1,type="l",lty=2)
  lines(out$grid,out$densm[ind,],lwd=2,type="l",lty=2)
  curve(dgamma(x,10,2),lwd=2,add=TRUE)
dev.off()