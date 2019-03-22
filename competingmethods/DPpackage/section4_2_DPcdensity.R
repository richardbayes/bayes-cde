library(DPpackage)

alldat <- read.csv("../../data/section4_2.csv",header=FALSE)
dat <- alldat[,1:3]
#y <- dat[,1]
#X <- dat[,2:3]
dat <- dat[1:10, 1:2]
y <- dat[,1]
X <- dat[,2]
xpred <- dat[1, 2]

w = cov(dat)
wbar = apply(dat,2,mean)
prior <- list(a0=1,
              # b0=.01,
              b0 = 1,
              nu1=2 + ncol(dat),
              nu2 = 2 + ncol(dat),
              psiinv2=2*solve(w),
              tau1=6.01,
              tau2=2.01,
              m2=wbar,
              s2=.5*w)
mcmc <- list(nburn=10^1,nsave=10^1,ndisplay=1,nskip=1)
# xpred <- matrix(c(.75,.5,
#                   .5,.75,
#                   .76,.76,
#                   .74,.76,
#                   .76,.74,
#                   .74,.74,
#                   .8,.1,
#                   .1,.8,
#                   .5,.5,
#                   .9,.9),ncol=2,byrow=TRUE)

begin <- proc.time()
set.seed(29)
# With y
out2 <- DPcdensity(y=y, x = X,xpred = xpred,prior=prior,state = NULL,
                   mcmc=mcmc,status=TRUE,compute.band=TRUE, type.band = "PD")
end <- proc.time()
save.image('output/sim6_DPcdensity.RData')
# load('output/sim6_DPcdensity.RData')

if (FALSE) {
  # How well did it do?
  # Traceplots
  plot(out2$save.state[[1]][,1],type='l')
  
  plot(out2$save.state[[1]][,2],type='l')
  plot(out2$save.state[[1]][,3],type='l')
  plot(out2$save.state[[1]][,4],type='l')
  plot(out2$save.state[[1]][,5],type='l')
  plot(out2$save.state[[1]][,6],type='l')
  plot(out2$save.state[[1]][,7],type='l')
  plot(out2$save.state[[1]][,8],type='l')
  plot(out2$save.state[[1]][,9],type='l')
  plot(out2$save.state[[1]][,10],type='l')
  plot(out2$save.state[[1]][,11],type='l')
  plot(out2$save.state[[1]][,12],type='l')
  
  plot(out2)
  
  for(i in 1:nrow(xpred)){
    plot(out2$grid,out2$densp.h[i,],lwd=1,type="l",lty=2,
         xlab="values",ylab="density",main=paste0('x1=',xpred[i,1],
                                                  ', x2=',xpred[i,2]))
    
    lines(out2$grid,out2$densp.l[i,],lwd=1,type="l",lty=2)
    lines(out2$grid,out2$densp.m[i,],lwd=2,type="l",lty=1)
    #curve(dnorm(x,f(xpred[i,1]),dnorm(xpred[i,1],2,1)^2),
    #      col='red',add=TRUE)
  }
  
  
  
  #jpeg('../Paper/Biometrika20170718/figs/DPcdensity.jpg',
  #     res=300,width=6,height=6,units='in')
  setEPS()
  postscript('../Paper/Biometrika20170718/figs/DPcdensityclean.eps',
             width=6,height=6)
    par(mfrow=c(2,2))
    
    # .76,.76
    ind <- 3
    plot(out2$grid,out2$densp.h[ind,],lwd=1,type="l",lty=2,
         xlab="Y",ylab="Density",ylim=c(0,.6),
         main=expression('Posterior Density at x'[1]*'=.76, x'[2]*'=.76'),
         las=1)
    
    lines(out2$grid,out2$densp.l[ind,],lwd=1,type="l",lty=2)
    lines(out2$grid,out2$densp.m[ind,],lwd=2,type="l",lty=2)
    curve(dnorm(x,1,sqrt(.5)),add=TRUE,lwd=2)
    
    # .9, .9
    ind <- 10 
    plot(out2$grid,out2$densp.h[ind,],lwd=1,type="l",lty=2,
         xlab="Y",ylab="Density",
         ylim=c(0,.6),las=1,
         main=expression('Posterior Density at x'[1]*'=.9, x'[2]*'=.9'))
    
    lines(out2$grid,out2$densp.l[ind,],lwd=1,type="l",lty=2)
    lines(out2$grid,out2$densp.m[ind,],lwd=2,type="l",lty=2)
    curve(dnorm(x,1,sqrt(.5)),lwd=2,add=TRUE)
    
    # .1, .8
    ind <- 8 
    plot(out2$grid,out2$densp.h[ind,],lwd=1,type="l",lty=2,
         xlab="Y",ylab="Density",las=1,
         main=expression('Posterior Density at x'[1]*'=.1, x'[2]*'=.8'))
    
    lines(out2$grid,out2$densp.l[ind,],lwd=1,type="l",lty=2)
    lines(out2$grid,out2$densp.m[ind,],lwd=2,type="l",lty=2)
    newdens <- function(x) .5*dnorm(x,1,1) + .5*dnorm(x,5,1)
    curve(newdens(x),lwd=2,add=TRUE)
    
    # .8, .1
    ind <- 7
    plot(out2$grid,out2$densp.h[ind,],lwd=1,type="l",lty=2,
         xlab="Y",ylab="Density",las=1,
         main=expression('Posterior Density at x'[1]*'=.8, x'[2]*'=.1'))
    
    lines(out2$grid,out2$densp.l[ind,],lwd=1,type="l",lty=2)
    lines(out2$grid,out2$densp.m[ind,],lwd=2,type="l",lty=2)
    curve(dgamma(x,10,2),lwd=2,add=TRUE)
  dev.off()
  
  
  # .5, .5
  ind <- 9
  plot(out2$grid,out2$densp.h[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="density",ylim=c(0,.6),
       main=expression('Posterior Density at x'[1]*'=.76, x'[2]*'=.76'))
  
  lines(out2$grid,out2$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(out2$grid,out2$densp.m[ind,],lwd=2,type="l",lty=2)
  curve(dgamma(x,10,2),col='red',add=TRUE)
  curve(newdens(x),col='red',add=TRUE)
}






