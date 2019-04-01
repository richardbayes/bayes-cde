library(DPpackage)

alldat <- read.csv("../../data/section4_2.csv",header=FALSE)
dat <- alldat[,1:3]
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
end <- proc.time()
 
save.image('output/section4_2_LDBDPdensity.RData')