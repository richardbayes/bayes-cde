library(DPpackage)
library(splines)

alldat <- read.csv("../data/sim6.csv",header=FALSE)
dat <- alldat[,1:3]
y <- dat[,1]
X <- dat[,2:3]
colnames(X) <- c("X1",'X2')

# Try using splines later...

X_design <- as.matrix(cbind(intercept=rep(1,nrow(X)),X))
xtf <- X_design
xpred <- matrix(c(.76,.76,
                  .74,.76,
                  .76,.74,
                  .8,.1,
                  .1,.8,
                  .9,.9),ncol=2,byrow=TRUE)
colnames(xpred) <- c("X1",'X2')
xpred <- cbind(intercept=rep(1,nrow(xpred)),xpred)
prediction <- list(xdenpred=xpred,
                   xtfdenpred=xpred,
                   xmedpred=xpred,
                   xtfmedpred=xpred,
                   quans=c(.025,.5,.975))
maxm <- round(log(nrow(X)/7.5,base=2)) # recommended in the paper
a0 <- 5; b0 <- 1; # used in the example in their paper
mub <- 0
Sb <- diag(ncol(X_design))*10^3
tau1 <- tau2 <- 10^(-4)
prior <- list(maxm=maxm,
              a0=a0,
              b0=b0,
              mub=mub,
              Sb=Sb,
              tau1=tau1,
              tau2=tau2)
mcmc <- list(nburn=2*10^4,nskip=2,nsave=2*10^4,ndisplay=1000)

begin <- proc.time()
  set.seed(5792)
  out <- LDTFPdensity(y,X_design,xtf,prediction,prior,mcmc,
               state=NULL,status=TRUE,compute.band=TRUE)
  # with splines to see if things improve
  df <- 5
  x1spline <- bs(X[,1],df)
  x2spline <- bs(X[,2],df)
  X_design2 <- as.matrix(cbind(intercept=rep(1,nrow(X)),x1spline,x2spline))
  xtf2 <- X_design2
  xpred2 <- cbind(xpred[,1],bs(xpred[,2],df),bs(xpred[,3],df))
  Sb <- diag(ncol(X_design2))*10^3
  prior2 <- list(maxm=maxm,
                a0=a0,
                b0=b0,
                mub=mub,
                Sb=Sb,
                tau1=tau1,
                tau2=tau2)
  prediction2 <- list(xdenpred=xpred2,
                     xtfdenpred=xpred2,
                     xmedpred=xpred2,
                     xtfmedpred=xpred2,
                     quans=c(.025,.5,.975))
  set.seed(10223)
  out2 <- LDTFPdensity(y,X_design2,xtf2,prediction2,prior2,mcmc,
                      state=NULL,status=TRUE,compute.band=TRUE)
  
  # With an interaction
  X_design3 <- as.matrix(cbind(X_design,X_design[,2]*X_design[,3]))
  colnames(X_design3)[4] <- 'X1*X2'
  xtf3 <- X_design3
  xpred3 <- cbind(xpred,xpred[,2]*xpred[,3])
  Sb <- diag(ncol(X_design3))*10^3
  prior3 <- list(maxm=maxm,
                 a0=a0,
                 b0=b0,
                 mub=mub,
                 Sb=Sb,
                 tau1=tau1,
                 tau2=tau2)
  prediction3 <- list(xdenpred=xpred3,
                      xtfdenpred=xpred3,
                      xmedpred=xpred3,
                      xtfmedpred=xpred3,
                      quans=c(.025,.5,.975))
  set.seed(3033033)
  out3 <- LDTFPdensity(y,X_design3,xtf3,prediction3,prior3,mcmc,
                       state=NULL,status=TRUE,compute.band=TRUE)
end <- proc.time()
#save.image('sim6_LDTFPdensity.RData')
#load('sim6_LDTFPdensity.RData')

# How well did it do?
# Traceplots
plot(out$save.state[[1]][,1],type='l')
plot(out$save.state[[1]][,2],type='l')
plot(out$save.state[[1]][,3],type='l')
plot(out$save.state[[1]][,4],type='l')
plot(out$save.state[[1]][,5],type='l')
plot(out$save.state[[1]][,6],type='l')
plot(out$save.state[[1]][,7],type='l')
plot(out$save.state[[1]][,8],type='l')
plot(out$save.state[[1]][,9],type='l')
plot(out$save.state[[1]][,10],type='l')
plot(out$save.state[[1]][,11],type='l')
plot(out$save.state[[1]][,12],type='l')

# plot(out)

tmpout <- out
  

# jpeg('../Paper/Biometrika20170718/figs/LDTFPdensity.jpg',
#      res=300,width=6,height=6,units='in')
setEPS()
postscript('../Paper/Biometrika20170718/figs/LDTFPdensityclean.eps',
     width=6,height=6)
  par(mfrow=c(2,2))
  # .76,.76
  ind <- 1
  plot(tmpout$grid,tmpout$densu[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",ylim=c(0,.8),las=1,
       main=expression('Posterior Density at x'[1]*'=.76, x'[2]*'=.76'))
  lines(tmpout$grid,tmpout$densl[ind,],lwd=1,type="l",lty=2)
  lines(tmpout$grid,tmpout$densm[ind,],lwd=2,type="l",lty=2)
  curve(dnorm(x,1,sqrt(.5)),lwd=2,add=TRUE)
  
  
  # .9, .9
  ind <- 6
  plot(tmpout$grid,tmpout$densu[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",las=1,
       ylim=c(0,.8),
       main=expression('Posterior Density at x'[1]*'=.9, x'[2]*'=.9'))
  
  lines(tmpout$grid,tmpout$densl[ind,],lwd=1,type="l",lty=2)
  lines(tmpout$grid,tmpout$densm[ind,],lwd=2,type="l",lty=2)
  curve(dnorm(x,1,sqrt(.5)),lwd=2,add=TRUE)
  
  # .1, .8
  ind <- 5
  plot(tmpout$grid,tmpout$densu[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",las=1,
       main=expression('Posterior Density at x'[1]*'=.1, x'[2]*'=.8'))
  
  lines(tmpout$grid,tmpout$densl[ind,],lwd=1,type="l",lty=2)
  lines(tmpout$grid,tmpout$densm[ind,],lwd=2,type="l",lty=2)
  newdens <- function(x) .5*dnorm(x,1,1) + .5*dnorm(x,5,1)
  curve(newdens(x),lwd=2,add=TRUE)
  
  # .8, .1
  ind <- 4
  plot(tmpout$grid,tmpout$densu[ind,],lwd=1,type="l",lty=2,
       xlab="Y",ylab="Density",ylim=c(0,.3),las=1,
       main=expression('Posterior Density at x'[1]*'=.8, x'[2]*'=.1'))
  lines(tmpout$grid,tmpout$densl[ind,],lwd=1,type="l",lty=2)
  lines(tmpout$grid,tmpout$densm[ind,],lwd=2,type="l",lty=2)
  curve(dgamma(x,10,2),lwd=2,add=TRUE)
dev.off()
  
  






