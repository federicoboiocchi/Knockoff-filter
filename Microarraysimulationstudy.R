# Author: Federico Boiocchi

# Microarray simulation study
# Power and FDR comparison using two-sample T tests (Naive,Bonferroni,Holm,Bhq)
rm(list = ls())
par(mfrow = c(1, 1))
par(pty = "m")
M <- 1000 # Montecarlo iterations
AA <- seq(0, 1.8, length.out = 100) # signal magnitude
len <- length(AA)
pwnaive <- numeric(len)
pwbonf <- numeric(len)
pwholm <- numeric(len)
pwbhq <- numeric(len)
fdrnaive <- numeric(len)
fdrbonf <- numeric(len)
fdrholm <- numeric(len)
fdrbhq <- numeric(len)
alpha <- 0.2 # FDR upper bound
n <- 50 # number of observations
p <- 100 # number of featurs
tp <- 15 # number of true significant features
fp <- p - tp
n1 <- 25 # dimension of the first group
n2 <- 25 # dimension of the second group
sd <- 1 # shared variance
# we are assuming unknown variance equal among groups
# (in simulation we fix a standard deviation of sd=1)

for (i in 1:len) {
  # cycle for different signal amplitudes
  A <- AA[i]
  pwnaive_iter <- numeric(M)
  pwbonf_iter <- numeric(M)
  pwholm_iter <- numeric(M)
  pwbhq_iter <- numeric(M)
  fdrnaive_iter <- numeric(M)
  fdrbonf_iter <- numeric(M)
  fdrholm_iter <- numeric(M)
  fdrbhq_iter <- numeric(M)

  for (m in 1:M) {
    # cycle for Montecarlo iterations 
    mu <- numeric(p)
    ind_tp <- sample(1:p, size = tp, replace = FALSE)
    ind_fp <- c(1:p)[-ind_tp]
    mu[ind_tp] <- A
    mu[ind_fp] <- 0
    pval <- numeric(p)

    # two sample T test with equal variance are done independently

    for (j in 1:p) {
      # cycle to test all p hypotheses using a Two sample T statistics
      x1 <- rnorm(n1, mean = 0, sd = sd)
      x2 <- rnorm(n2, mean = mu[j], sd = sd)
      x1bar <- mean(x1)
      x2bar <- mean(x2)
      var1 <- var(x1)
      var2 <- var(x2)
      varPool <- ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
      T_stat <- (x1bar - x2bar) / sqrt(varPool * (1 / n1 + 1 / n2))
      pval[j] <- 2 * (1 - pt(abs(T_stat), df = n1 + n2 - 2))
    }
    pval_asc <- sort(pval)
    indices <- c(1:p)
    # Naive
    ind_naive <- which(pval <= alpha)
    totp_naive <- length(ind_naive)
    tp_naive <- sum(ind_naive %in% ind_tp)
    fp_naive <- totp_naive - tp_naive
    # Bonferroni
    ind_bonf <- which(pval <= alpha / p)
    totp_bonf <- length(ind_bonf)
    tp_bonf <- sum(ind_bonf %in% ind_tp)
    fp_bonf <- totp_bonf - tp_bonf
    # Holm
    i0 <- min(which(pval_asc > alpha / (p - indices + 1)))
    ind_holm <- which(pval < alpha / (p - i0 + 1))
    totp_holm <- length(ind_holm)
    tp_holm <- sum(ind_holm %in% ind_tp)
    fp_holm <- totp_holm - tp_holm
    # BHq
    imax <- max(which(pval_asc <= (indices / p) * alpha))
    ind_bhq <- which(pval <= (imax / p) * alpha)
    totp_bhq <- length(ind_bhq)
    tp_bhq <- sum(ind_bhq %in% ind_tp)
    fp_bhq <- totp_bhq - tp_bhq

    pwnaive_iter[m] <- tp_naive / tp
    pwbonf_iter[m] <- tp_bonf / tp
    pwholm_iter[m] <- tp_holm / tp
    pwbhq_iter[m] <- tp_bhq / tp

    fdrnaive_iter[m] <- fp_naive / max(totp_naive, 1)
    fdrbonf_iter[m] <- fp_bonf / max(totp_bonf, 1)
    fdrholm_iter[m] <- fp_holm / max(totp_holm, 1)
    fdrbhq_iter[m] <- fp_bhq / max(totp_bhq, 1)
  }
  pwnaive[i] <- mean(pwnaive_iter)
  pwbonf[i] <- mean(pwbonf_iter)
  pwholm[i] <- mean(pwholm_iter)
  pwbhq[i] <- mean(pwbhq_iter)

  fdrnaive[i] <- mean(fdrnaive_iter)
  fdrbonf[i] <- mean(fdrbonf_iter)
  fdrholm[i] <- mean(fdrholm_iter)
  fdrbhq[i] <- mean(fdrbhq_iter)
}

# Power comparison (Figure 3.3)
plot(AA, pwbhq,
  type = "l", col = "blue", main = "Power Comparison", ylab = "Power",
  xlab = "Signal Magnitude", lwd = 1.25
)
points(AA, pwholm, type = "l", col = "red", lwd = 1.25)
points(AA, pwbonf, type = "l", col = "deepskyblue", lty = "dashed", lwd = 1.25)
points(AA, pwnaive, type = "l", col = "green", lwd = 1.25)
legend(
  x = 1.35, y = 0.8,
  legend = c("Naive", "BHq", "Holm", "Bonf."),
  col = c("green", "blue", "red", "deepskyblue"), # line colors
  lty = c(1, 1, 1, 2), # line types
  cex = 0.65,
  lwd = 1.25,
)


# FDR comparison (Figure 3.4)
plot(AA, fdrbhq,
  type = "l", ylim = c(-0.01, 0.9), col = "blue", main = "FDR comparison", xlab = "Signal magnitude",
  ylab = "FDR", lwd = 1.25
)
points(AA, fdrholm, , type = "l", col = "red", lwd = 1.25)
lines(AA, fdrbonf, lty = "dashed", col = "deepskyblue", lwd = 1.25)
points(AA, fdrnaive, type = "l", col = "green", lwd = 1.25)
abline(h = alpha, col = "black", lty = "dashed", lwd = 1.25)
legend(
  x = 1.38, y = 0.90,
  legend = c("Naive", "BHq", "Holm", "Bonf.", expression(alpha)),
  col = c("green", "blue", "red", "deepskyblue", "black"), # line colors
  lty = c(1, 1, 1, 2, 2), # line types
  cex = 0.65,
  lwd = 1.25,
)

