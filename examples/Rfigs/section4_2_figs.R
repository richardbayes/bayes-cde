dat <- read.csv("data/section4_2.csv", header = FALSE)
colnames(dat) <- c(
  paste0("xt", 1:3),
  paste0("p", 1:3),
  "pq1l", "pq1u", "pq2l", "pq2u", "pq3l", "pq3u"
)

setEPS()
postscript("figs/section4_2dens.eps", width = 7, height = 4)
par(mfrow = c(1, 3))
lwd1 <- 2
plot(
  dat$xt1,
  dat$p1,
  lty = 2,
  lwd = lwd1,
  type = "l",
  main = expression(paste(X[1] > X[2], ", ", X[2] < .75)),
  ylim = c(0, .3),
  ylab = "Density",
  xlab = "",
  las = 1
)
lines(dat$xt1, dat$pq1l, lty = 2)
lines(dat$xt1, dat$pq1u, lty = 2)
curve(dgamma(x, 10, 2), add = TRUE, lwd = lwd1)
plot(
  dat$xt2,
  dat$p2,
  lty = 2,
  lwd = lwd1,
  type = "l",
  main = expression(paste(X[1] < X[2], ", ", X[1] < .75)),
  ylim = c(0, .25),
  ylab = "",
  xlab = "Y",
  las = 1
)
lines(dat$xt2, dat$pq2l, lty = 2)
lines(dat$xt2, dat$pq2u, lty = 2)
mixdens <- function(x) .5 * dnorm(x, 1, 1) + .5 * dnorm(x, 5, 1)
curve(mixdens(x), add = TRUE, lwd = lwd1)
plot(
  dat$xt3,
  dat$p3,
  lty = 2,
  lwd = lwd1,
  type = "l",
  main = expression(paste(X[1] > .75, ", ", X[2] > .75)),
  ylim = c(0, .7),
  xlab = "",
  ylab = "",
  las = 1
)
lines(dat$xt3, dat$pq3l, lty = 2)
lines(dat$xt3, dat$pq3u, lty = 2)
curve(dnorm(x, 1, sqrt(.5)), add = TRUE, lwd = lwd1)
dev.off()