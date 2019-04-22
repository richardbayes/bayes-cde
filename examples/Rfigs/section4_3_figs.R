dat <- read.csv("data/section4_3_densities.csv", header = FALSE)
datx <- read.csv("data/section4_3_X.csv", header = FALSE)
colnames(dat) <- c(
  paste0("xt", 1:3),
  paste0("p", 1:3),
  "pq1l", "pq1u", "pq2l", "pq2u", "pq3l", "pq3u"
)
colnames(datx) <- c("x1", "x2", "prt")
X <- as.matrix(datx[, 1:2])
mu <- c(1, 5, 7.5)
Sigma <- matrix(c(1, .5, .1, .5, 1, .5, .1, .5, 1), 3, 3)
# center index on X
theinds <- c(20, 387, 4)
mu1 <- mu[1]
mu2 <- mu[2:3]
sigma11 <- Sigma[1, 1, drop = FALSE]
sigma12 <- Sigma[1, 2:3, drop = FALSE]
sigma21 <- t(sigma12)
sigma22 <- Sigma[2:3, 2:3, drop = FALSE]
get_mean <- function(x) {
  mu1 + sigma12 %*% solve(sigma22, x - mu2)
}
sd <- sigma11 - sigma12 %*% solve(sigma22, sigma21)


setEPS()
postscript("figs/section4_3.eps", width = 8, height = 3)
par(mfrow = c(1, 4))
plot(
  datx$x1,
  datx$x2,
  pch = c(1, 3, 6)[datx$prt],
  col = c("cornflowerblue", "darkolivegreen3", "firebrick")[datx$prt],
  xlab = expression(X[1]),
  ylab = expression(X[2]),
  las = 1
)
lwd1 <- 2
plot(
  dat$xt1,
  dat$p1,
  lty = 2,
  lwd = lwd1,
  type = "l",
  main = expression(paste(X[1] == 5.1, ", ", X[2] == 5.6)),
  ylim = c(0, .6),
  ylab = "Density",
  xlab = "Y",
  las = 1
)
lines(dat$xt1, dat$pq1l, lty = 2)
lines(dat$xt1, dat$pq1u, lty = 2)
curve(dnorm(x, get_mean(X[theinds[1], ]), sd), add = TRUE, lwd = lwd1)
plot(
  dat$xt2,
  dat$p2,
  lty = 2,
  lwd = lwd1,
  type = "l",
  main = expression(paste(X[1] == 4.3, ", ", X[2] == 6.8)),
  ylim = c(0, .6),
  ylab = "",
  xlab = "Y",
  las = 1
)
lines(dat$xt2, dat$pq2l, lty = 2)
lines(dat$xt2, dat$pq2u, lty = 2)
curve(dnorm(x, get_mean(X[theinds[2], ]), sd), add = TRUE, lwd = lwd1)
plot(
  dat$xt3,
  dat$p3,
  lty = 2,
  lwd = lwd1,
  type = "l",
  main = expression(paste(X[1] == 2.3, ", ", X[2] == 6.8)),
  ylim = c(0, .6),
  xlab = "Y",
  ylab = "",
  las = 1
)
lines(dat$xt3, dat$pq3l, lty = 2)
lines(dat$xt3, dat$pq3u, lty = 2)
curve(dnorm(x, get_mean(X[theinds[3], ]), sd), add = TRUE, lwd = lwd1)
dev.off()