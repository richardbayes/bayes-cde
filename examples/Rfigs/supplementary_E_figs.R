dat <- read.csv("data/supplementary_E.csv", header = FALSE)
colnames(dat) <- c(
  paste0("xt", 1:3),
  paste0("p", 1:3),
  "pq1l", "pq1u", "pq2l", "pq2u", "pq3l", "pq3u"
)
setEPS()
postscript("figs/supplementary_E.eps", width = 8, height = 3)
  par(mfrow = c(1, 3))
  lwd1 <- 2
  plot(
    dat$xt1,
    dat$p1,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = "X < .25",
    ylim = c(0, 10),
    xlim = c(0, 1),
    ylab = "Density",
    xlab = "",
    las = 1
  )
  lines(dat$xt1, dat$pq1l, lty = 2)
  lines(dat$xt1, dat$pq1u, lty = 2)
  curve(dbeta(x, 30, 20), add = TRUE, lwd = lwd1)
  plot(
    dat$xt3,
    dat$p3,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = ".25 < X < .5",
    ylim = c(0, 10),
    xlim = c(0, 1),
    ylab = "",
    xlab = "Y",
    las = 1
  )
  lines(dat$xt3, dat$pq3l, lty = 2)
  lines(dat$xt3, dat$pq3u, lty = 2)
  curve(dbeta(x, 10, 30), add = TRUE, lwd = lwd1)
  plot(
    dat$xt2,
    dat$p2,
    lty = 2,
    lwd = lwd1,
    type = "l",
    main = "X > .75",
    ylim = c(0, 10),
    xlim = c(0, 1),
    ylab = "",
    xlab = "",
    las = 1
  )
  lines(dat$xt2, dat$pq2l, lty = 2)
  lines(dat$xt2, dat$pq2u, lty = 2)
  curve(dbeta(x, .5, .5), add = TRUE, lwd = lwd1)
dev.off()