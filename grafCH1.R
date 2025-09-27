rm(list = ls())
# Sample mean on the x axis
# arbitrarily chosen quantities to allow
# an effective representation
n <- 1000
sd <- 1
mu0 <- 0
mu1 <- 0.08
mi <- -0.15
ma <- 0.2
x <- seq(mi, ma, length.out = 10000)

# Values used to plot the densities under H_0 and H_1
y0 <- dnorm(x, mean = mu0, sd = sd / sqrt(n))
y1 <- dnorm(x, mean = mu1, sd = sd / sqrt(n))

# Critical threshold and probability of committing a first type error
alpha <- 0.2
ct <- mu0 + qnorm(1 - alpha) * (sd / sqrt(n))

# plot of the two densities (Figure 1.1)

plot(x, y0,
  type = "l", lwd = 1, col = "black",
  xlim = c(mi, ma), ylim = c(0, max(y0, y1)),
  xlab = expression(bar(x)), ylab = "probability density",
  main = bquote("Sample mean PDFs" ~ alpha == 0.2)
)

points(x, y1, type = "l", lwd = 1, col = "black")
abline(h = 0)
segments(ct, 0, ct, max(dnorm(ct, mean = mu1, sd = sd / sqrt(n)), dnorm(ct, mean = mu0, sd = sd / sqrt(n))), lwd = 1.5)
x1 <- x[x <= ct]
x2 <- x[x >= ct]

# Probability of committing a type II error (False Negative)
polygon(c(x1, rev(x1)),
  c(dnorm(x1, mu1, sd / sqrt(n)), rep(0, length(x1))),
  density = 20, angle = 45, col = "blue", border = NA
)

# Probability of committing a type I error (False positive)
polygon(c(x2, rev(x2)),
  c(dnorm(x2, mu0, sd / sqrt(n)), rep(0, length(x2))),
  density = 40, angle = -45, col = "red", border = NA
)

# Probability of making a true discovery (True positive)
polygon(c(x2, rev(x2)),
  c(dnorm(x2, mu1, sd / sqrt(n)), rep(0, length(x2))),
  density = 10, angle = 30, col = "green", border = NA
)

segments(ct, 0, ma + 1, 0, col = "purple", lwd = 4)
segments(mu0, 0, mu0, dnorm(mu0, mu0, sd / sqrt(n)), col = "black", lwd = 1, lty = "dashed")
segments(mu1, 0, mu1, dnorm(mu1, mu1, sd / sqrt(n)), col = "black", lwd = 1, lty = "dashed")

text(-0.04, 11, labels = expression(H[0]))
text(0.12, 11, labels = expression(H[1]))
legend(-0.15, 12,
  fill = c("blue", "green", "red"), bty = "n",
  legend = c(expression(beta), expression(pi), expression(alpha)),
  cex = 1, density = c(30, 30, 30), angle = c(45, 30, -45)
)

text(0.19, 0.8, labels = c("R"), col = "purple")
segments(mi - 1, 0, ct, 0, col = "dodgerblue", lwd = 4)
text(-0.11, 0.8, labels = c("A"), col = "dodgerblue")

# Power function (Figure 1.2)
alpha <- 0.25

ct <- mu0 + qnorm(1 - alpha) * (sd / sqrt(n))

power <- function(mu0, mu, alpha, sd, n) {
  p <- numeric(length(mu))
  for (i in 1:length(mu)) {
    p[i] <- 1 - pnorm(((mu0 - mu[i]) / (sd / sqrt(n))) + qnorm(1 - alpha))
  }
  return(p)
}
mm <- seq(-0.5, 0.5, length.out = 2000)
pp <- power(mu0, mu = mm, alpha, sd, n)
pp
min <- -0.1
max <- 0.15
plot(mm, pp,
  type = "l", col = "black", lwd = 1.5, xlim = c(min, max),
  main = "Power function", xlab = expression(mu),
  ylab = expression(pi * "(" * mu * ")")
)
abline(v = mu0, lty = "dashed")
abline(v = ct, lty = "dashed")
abline(h = 0)
abline(h = 1)
segments(0, 0, 0.3, 0, lwd = 3, col = "darkorange")
segments(-0.3, 0, 0, 0, lwd = 3, col = "lightblue")
text(0.12, 0.07, labels = c(expression(H[1])), col = "darkorange", cex = 1.2)
text(-0.087, 0.07, labels = c(expression(H[0])), col = "lightblue", cex = 1.2)
text(-0.01, 0.8, labels = c(expression(mu[0])))
text(-0.095, 0.29, labels = c(expression(alpha)), col = "red2", cex = 1.2)
text(ct + 0.03, 0.45, labels = expression(mu[0] + z[1 - alpha] * sigma / sqrt(n)))
eps <- 0.009
segments(min - eps, 0, min - eps, alpha, lwd = 3, col = "red2")
segments(min - eps, alpha, 0, alpha, lty = "dashed")
segments(min - eps, alpha, min - eps, 1, lwd = 3, col = "green")

# Power comparison (Figure 1.5)

alpha <- 0.1
ct <- mu0 + qnorm(1 - alpha) * (sd / sqrt(n))

power <- function(mu0, mu, alpha, sd, n) {
  p <- numeric(length(mu))
  for (i in 1:length(mu)) {
    p[i] <- 1 - pnorm(((mu0 - mu[i]) / (sd / sqrt(n))) + qnorm(1 - alpha))
  }
  return(p)
}
mm <- seq(-0.5, 0.5, length.out = 2000)
pp <- power(mu0, mu = mm, alpha, sd, n)
pp
min <- -0.1
max <- 0.15
plot(mm, pp,
  type = "l", col = "black", lwd = 1.5, xlim = c(min, max),
  main = "Power function comparison", xlab = expression(mu),
  ylab = expression(pi * "(" * mu * ")")
)
alpha1 <- 0.01
pp1 <- power(mu0, mu = mm, alpha1, sd, n)
points(mm, pp1, col = "red", type = "l", lwd = 1.5)
alpha2 <- 0.5
pp2 <- power(mu0, mu = mm, alpha2, sd, n)
points(mm, pp2, col = "turquoise", type = "l", lwd = 1.5)
legend(
  x = -0.10, y = 0.8,
  legend = c(
    expression(alpha == 0.5),
    expression(alpha == 0.1),
    expression(alpha == 0.01)
  ),
  col = c("turquoise", "black", "red"),
  lty = 1,
  lwd = 1.5,
  bty = "n",
  cex = 0.8
)
abline(v = mu0, lty = "dashed")
abline(h = 0)
abline(h = 1)
segments(0, 0, 0.3, 0, lwd = 3, col = "darkorange")
segments(-0.3, 0, 0, 0, lwd = 3, col = "lightblue")
text(0.12, 0.07, labels = c(expression(H[1])), col = "darkorange", cex = 1.2)
text(-0.087, 0.07, labels = c(expression(H[0])), col = "lightblue", cex = 1.2)
text(-0.01, 0.8, labels = c(expression(mu[0])))
