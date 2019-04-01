library(DPpackage)

alldat <- read.csv("../../data/section4_2.csv",header=FALSE)
dat <- alldat[,1:3]
y <- dat[,1]
X <- dat[,2:3]
colnames(X) <- c("X1",'X2')

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
end <- proc.time()
save.image('output/section4_2_LDTFPdensity.RData')