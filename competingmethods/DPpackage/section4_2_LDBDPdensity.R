library(DPpackage)
library(splines)


alldat <- read.csv("../data/sim6.csv",header=FALSE)
datorig <- dat <- alldat[,1:3]
colnames(dat) <- c('y','x1','x2')
# Standardize y
dat$y <- (dat$y - min(dat$y))/(max(dat$y) - min(dat$y))
xpred <- matrix(c(.76,.76,
                  .74,.76,
                  .76,.74,
                  .8,.1,
                  .1,.8,
                  .9,.9),ncol=2,byrow=TRUE)
colnames(xpred) <- c("X1",'X2')
xpred <- cbind(intercept=rep(1,nrow(xpred)),xpred)
prior <- list(lambda=25,
              maxn=25,
              alpha=1,
              m0=rep(0,3), # should be the number of betas/coefficients
              S0=diag(3)*10^3,
              nu=5, # number of betas + 2
              psiinv=diag(3)*10^3)
mcmc <- list(nburn=10^4,nskip=2,ndisplay=1000,nsave=2*10^4)

begin <- proc.time()
  set.seed(10101)
  out <- LDBDPdensity(dat$y ~ dat$x1 + dat$x2,xpred,prior,mcmc,
               state=NULL,status=TRUE,compute.band=TRUE)
  # With Splines
  df <- 5
  xpred2 <- cbind(xpred[,1],bs(xpred[,2],df),bs(xpred[,3],df))
  prior2 <- list(lambda=25,
                maxn=25,
                alpha=1,
                m0=rep(0,11), # should be the number of betas/coefficients
                S0=diag(11)*10^3,
                nu=13, # number of betas + 2
                psiinv=diag(11)*10^3)
  set.seed(20202)
  out2 <- LDBDPdensity(dat$y ~ bs(dat$x1,df) + bs(dat$x2,df),xpred2,prior2,mcmc,
                              state=NULL,status=TRUE,compute.band=TRUE)
  
  xpred3 <- cbind(xpred,xpred[,2]*xpred[,3])
  prior3 <- list(lambda=25,
                 maxn=25,
                 alpha=1,
                 m0=rep(0,4), # should be the number of betas/coefficients
                 S0=diag(4)*10^3,
                 nu=6, # number of betas + 2
                 psiinv=diag(4)*10^3)
  set.seed(30303)
  out3 <- LDBDPdensity(dat$y ~ dat$x1 + dat$x2 + I(dat$x1*dat$x2),
                       xpred3,prior3,mcmc,
                       state=NULL,status=TRUE,compute.band=TRUE)
end <- proc.time()

#df <- 5
#out2 <- LDBDPdensity(dat$y ~ bs(dat$x1,df) + bs(dat$x2,df),xpred,prior,mcmc,
#                    state=NULL,status=TRUE,
#                    compute.band=TRUE)
# save.image('sim6_LDBDPdensity.RData')
# load('sim6_LDBDPdensity.RData')

outtmp <- out

# How well did it do?
# Traceplots
plot(outtmp$save.state[[1]][,1],type='l')
plot(outtmp$save.state[[1]][,2],type='l')
plot(outtmp$save.state[[1]][,3],type='l')
plot(outtmp$save.state[[1]][,4],type='l')
plot(outtmp$save.state[[1]][,5],type='l')
plot(outtmp$save.state[[1]][,6],type='l')
plot(outtmp$save.state[[1]][,7],type='l')
plot(outtmp$save.state[[1]][,8],type='l')
plot(outtmp$save.state[[1]][,9],type='l')
plot(outtmp$save.state[[1]][,10],type='l')
plot(outtmp$save.state[[1]][,11],type='l')
plot(outtmp$save.state[[1]][,12],type='l')

# plot(outtmp)

for(i in 1:nrow(xpred)){
  plot(outtmp$grid,outtmp$densp.h[i,],lwd=1,type="l",lty=2,
       xlab="values",ylab="density",main=paste0('x1=',xpred[i,1],
                                                ', x2=',xpred[i,2]))
  
  lines(outtmp$grid,outtmp$densp.l[i,],lwd=1,type="l",lty=2)
  lines(outtmp$grid,outtmp$densp.m[i,],lwd=2,type="l",lty=1)
  #curve(dnorm(x,f(xpred[i,1]),dnorm(xpred[i,1],2,1)^2),
  #      col='red',add=TRUE)
}


# jpeg('../Paper/Biometrika20170718/figs/LDBDPdensity.jpg',
#      res=300,width=6,height=6,units='in')
setEPS()
postscript('../Paper/Biometrika20170718/figs/LDBDPdensityclean.eps',
           width=6,height=6)
  par(mfrow=c(2,2))
  y <- datorig[,1]
  
  # .76,.76
  ind <- 1
  plot(outtmp$grid,outtmp$densp.h[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",ylim=c(0,8),las=1,
       main=expression('Posterior Density at x'[1]*'=.76, x'[2]*'=.76'))
  
  lines(outtmp$grid,outtmp$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(outtmp$grid,outtmp$densp.m[ind,],lwd=2,type="l",lty=2)
  themean <- (1 - min(y)) / (max(y) - min(y))
  thesd <- sqrt(.5)/(max(y) - min(y))
  curve(dnorm(x,themean,thesd),lwd=2,add=TRUE)
  
  # .9, .9
  ind <- 6
  plot(outtmp$grid,outtmp$densp.h[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",las=1,
       ylim=c(0,8),
       main=expression('Posterior Density at x'[1]*'=.9, x'[2]*'=.9'))
  
  lines(outtmp$grid,outtmp$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(outtmp$grid,outtmp$densp.m[ind,],lwd=2,type="l",lty=2)
  themean <- (1 - min(y)) / (max(y) - min(y))
  thesd <- sqrt(.5)/(max(y) - min(y))
  curve(dnorm(x,themean,thesd),lwd=2,add=TRUE)
  
  # .1, .8
  ind <- 5 
  plot(outtmp$grid,outtmp$densp.h[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",las=1,
       main=expression('Posterior Density at x'[1]*'=.1, x'[2]*'=.8'))
  
  lines(outtmp$grid,outtmp$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(outtmp$grid,outtmp$densp.m[ind,],lwd=2,type="l",lty=2)
  themean1 <- (1 - min(y)) / (max(y) - min(y))
  themean2 <- (5 - min(y)) / (max(y) - min(y))
  thesd <- 1/(max(y) - min(y))
  newdens <- function(x) .5*dnorm(x,themean1,thesd) + .5*dnorm(x,themean2,thesd)
  curve(newdens(x),lwd=2,add=TRUE)
  
  # .8, .1
  ind <- 4
  plot(outtmp$grid,outtmp$densp.h[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",ylim=c(0,3.5),las=1,
       main=expression('Posterior Density at x'[1]*'=.8, x'[2]*'=.1'))
  
  lines(outtmp$grid,outtmp$densp.l[ind,],lwd=1,type="l",lty=2)
  lines(outtmp$grid,outtmp$densp.m[ind,],lwd=2,type="l",lty=2)
  therange <- (max(y) - min(y))
  themean <- min(y)
  newgamma <- function(z) therange*dgamma(z*therange + themean,10,2)
  #curve(dgamma(x,10,2),col='red',add=TRUE)
  curve(newgamma(x),lwd=2,add=TRUE)
dev.off()






