library(DPpackage)

alldat <- read.csv("../../data/section4_2.csv",header=FALSE)
dat <- alldat[,1:3]
y <- dat[,1]
X <- dat[,2:3]

w = cov(dat)
wbar = apply(dat,2,mean)
prior <- list(a0=1,
              b0=.01,
              nu1=2 + ncol(dat),
              nu2 = 2 + ncol(dat),
              psiinv2=2*solve(w),
              tau1=6.01,
              tau2=2.01,
              m2=wbar,
              s2=.5*w)
mcmc <- list(nburn=10^4,nsave=10^4,ndisplay=1000,nskip=10)
xpred <- matrix(c(.75,.5,
                  .5,.75,
                  .76,.76,
                  .74,.76,
                  .76,.74,
                  .74,.74,
                  .8,.1,
                  .1,.8,
                  .5,.5,
                  .9,.9),ncol=2,byrow=TRUE)

begin <- proc.time()
set.seed(29)
out <- DPcdensity(y=y, x = X,xpred = xpred,prior=prior,state = NULL,
                   mcmc=mcmc,status=TRUE,compute.band=TRUE, type.band = "PD")
end <- proc.time()
save.image('output/section4_2_DPcdensity.RData')