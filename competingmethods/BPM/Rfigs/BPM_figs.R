dat <- read.csv("Rfigs/data/section4_2_BPM.csv", header = FALSE)
colnames(dat) <- c(
  "x",
  paste0("p", 1:4),
  "pq1l", "pq1u", "pq2l", "pq2u", "pq3l", "pq3u", "pq4l", "pq4u"
)

setEPS()
postscript("Rfigs/section4_2_BPM.eps", width = 8, height = 3)
  par(mfrow = c(1, 4))
  lwd1 <- 2
  plot(
    dat$x,
    dat$p3,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = expression(paste(X[1] == .76, ", ", X[2] == .76)),
    ylim = c(0, .6),
    ylab = "Density",
    xlab = "",
    las = 1
  )
  lines(dat$x, dat$pq3l, lty = 2)
  lines(dat$x, dat$pq3u, lty = 2)
  curve(dnorm(x, 1, sqrt(.5)), add = TRUE, lwd = lwd1)
  plot(
    dat$x,
    dat$p4,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = expression(paste(X[1] == .9, ", ", X[2] == .9)),
    ylim = c(0, .6),
    ylab = "",
    xlab = "",
    las = 1
  )
  lines(dat$x, dat$pq4l, lty = 2)
  lines(dat$x, dat$pq4u, lty = 2)
  curve(dnorm(x, 1, sqrt(.5)), add = TRUE, lwd = lwd1)
  plot(
    dat$x,
    dat$p1,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = expression(paste(X[1] == .1, ", ", X[2] == .8)),
    ylim = c(0, .25),
    ylab = "",
    xlab = "",
    las = 1
  )
  lines(dat$x, dat$pq1l, lty = 2)
  lines(dat$x, dat$pq1u, lty = 2)
  newdens <- function(x) .5 * dnorm(x, 1, 1) + .5 * dnorm(x, 5, 1)
  curve(newdens(x), add = TRUE, lwd = lwd1)
  plot(
    dat$x,
    dat$p2,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = expression(paste(X[1] == .8, ", ", X[2] == .1)),
    ylim = c(0, .3),
    ylab = "",
    xlab = "",
    las = 1
  )
  lines(dat$x, dat$pq2l, lty = 2)
  lines(dat$x, dat$pq2u, lty = 2)
  curve(dgamma(x, 10, 2), add = TRUE, lwd = lwd1)
  mtext("Y", 1, outer = TRUE, line = -2)
dev.off()