library(fields)
library(ggplot2)
setwd("//sf1/Users/grad/richard/My Documents/Mallick/Density Regression/examples/")
f <- function(x1,x2){
  5*exp(-10*(x2-5)^2) + 2*x1
}

x1s <- seq(from=0,to=3,by=.1)
x2s <- seq(from=0,to=10,by=.1)

Z <- matrix(NA,length(x1s),length(x2s))
for(i in 1:length(x1s)){
  for(j in 1:length(x2s)){
    Z[i,j] <- f(x1s[i],x2s[j])
  }
}

image.plot(x1s,x2s,Z)  
contour(x1s,x2s,Z)

set.seed(259253)
n <- 1000
x1 <- runif(n,0,3)
x2 <- runif(n,0,10)

sum(x2 < 5.25 & x2 > 4.75)

y <- f(x1,x2) + rnorm(n,0,.1)

plot(x1,x2,col=y)

dat <- data.frame(x1=x1,x2=x2,y=y)

ggplot(dat,aes(x=x1,y=x2,colour=y)) + 
  geom_point(size=3)




