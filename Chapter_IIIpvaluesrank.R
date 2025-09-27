# MULTIPLE TESTING PROCEDURES COMPARISON
# p-value vs rank plot

rm(list = ls())
library(tidyverse)
set.seed(123)

m <- 1000 # number of total hypotheses tested
m0 <- 900 #  number of true null H0
m1 <- m - m0 # number of true alternatives H1

# tni = true null indexes
tni <- sample(1:m, size = m0, replace = FALSE)
# tai = true alternative indexes
tai <- c(1:m)[-tni]

# the order of the hypotheses tested is not important
# they are invariant with respect to the index j

# Simulation

# under true null H0j zj ~ N(0,1)
# under the alternative H1j zj ~ N(muj,1)
# the variance is fixed at 1
# We consider muj = mu (for all j in tai)

# zj is the test statistics we will use to compute the p-value

zj <- numeric(m)

# We will fill this vector zj drawing from a N(0,1) for the indexes tni
# and from a N(muj,1) for the indexes tai. In other words we skip the
# data generation phase and we directly generate test statistics.
# we are assuming a known variance

muj <- 3
sd <- 1
zj[tni] <- rnorm(m0, mean = 0, sd)
zj[tai] <- rnorm(m1, mean = muj, sd)

# two tailed test p-value
pval <- 2 * (1 - pnorm(abs(zj)))
alpha <- 0.05
pval_asc <- sort(pval, decreasing = FALSE)
q <- alpha

upper <- m # maximum index we want to represent on the x-axis
pval_naive <- pval_asc[1:upper]
ind <- 1:upper

indx <- numeric()
for (j in 1:upper) {
  indx[j] <- which(pval_asc[j] == pval)
}
indx
color <- numeric()

for (j in 1:upper) {
  if (indx[j] %in% tai) {
    color[j] <- 1
  } else {
    color[j] <- 0
  }
}

# Holm's threshold
thr_holm <- function(index) {
  q / (m - index + 1)
}
color <- factor(color)
levels(color) <- c("TRUE H0", "TRUE H1")
data <- data.frame(ind, pval_naive)
xmax <- 150

# p-values vs rank plot (Figure 3.1)
plot1 <- ggplot(data = data[1:xmax, ], mapping = aes(
  x = ind, y = pval_naive,
  color = color[1:xmax]
), xlab = "index j", ylab = "p-value") +
  theme_light() +
  geom_point(alpha = .4, size = 2.5) +
  scale_color_manual(values = c("red", "green")) +
  geom_abline(intercept = 0, slope = q / m, color = "blue") +
  geom_text(data = data[1, ], aes(x = 140, y = 0.009, label = "BHq"), color = "blue", size = 4) +
  geom_abline(intercept = q / m, slope = 0, colour = "deepskyblue") +
  geom_text(data = data[1, ], aes(x = 145, y = 0.002, label = "Bonf."), color = "deepskyblue", size = 4) +
  geom_abline(intercept = q, slope = 0, colour = "green") +
  geom_text(data = data[1, ], aes(x = 10, y = 0.052, label = "Naive"), color = "green", size = 4) +
  geom_function(fun = thr_holm, color = "red", lty = "dashed") +
  geom_text(data = data[1, ], aes(x = 120, y = 0.002, label = "Holm"), color = "red", size = 4) +
  labs(x = "index j", y = "p-values", color = NULL, title = "p-values vs rank: Procedures Comparison")
plot1

# same graph in logarithmic scale base = e

logpval <- log(pval_naive)
data <- data.frame(ind, logpval)

thr_naive <- function(index) {
  log(q)
}
thr_bonf <- function(index) {
  log(q / m)
}
thr_holm <- function(index) {
  log(q / (m - index + 1))
}
thr_BHq <- function(index) {
  log((q * index) / m)
}
# Naive
ind_naive <- which(pval <= alpha)
(TOTP_naive <- (length(ind_naive)))
# Bonferroni
ind_bonf <- which(pval <= alpha / m)
(TOTP_bonf <- (length(ind_bonf)))
padj_bonf <- p.adjust(pval, method = "bonferroni")
ind_bonf_adj <- which(padj_bonf <= alpha)
# Holm
pval_asc <- sort(pval)
indices <- 1:m
i0 <- min(which(pval_asc > alpha / (m - indices + 1)))
ind_holm <- which(pval < alpha / (m - i0 + 1))
TOTP_holm <- length(ind_holm)
padj_holm <- p.adjust(pval, method = "holm")
ind_holm_adj <- which(padj_holm <= alpha)
# BHq
imax <- max(which(pval_asc <= (indices / m) * q))
ind_bhq <- which(pval <= (imax / m) * q)
TOTP_BH <- length(ind_bhq)
padj_bhq <- p.adjust(pval, method = "BH")
ind_bhq_adj <- which(padj_bhq <= alpha)

xmax <- 150
ln_breaks <- log(c(1, 0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6))
ln_labels <- c("1", "0.1", "0.01", "0.001", "1e-4", "1e-5", "1e-6")

# p-values vs rank plot (log-scale) (Figure 3.2)
plot2 <- ggplot(data = data[1:xmax, ], mapping = aes(
  x = ind, y = logpval,
  color = color[1:xmax]
), xlab = "index j", ylab = "p-value") +
  theme_light() +
  geom_point(alpha = .4, size = 2.5) +
  scale_color_manual(values = c("red", "green")) +
  geom_function(fun = thr_BHq, color = "blue") +
  geom_text(data = data[1, ], aes(x = xmax - 10, y = log((q * xmax) / m) - 0.5, label = "BHq"), color = "blue", size = 4) +
  geom_function(fun = thr_bonf, colour = "deepskyblue") +
  geom_text(data = data[1, ], aes(x = xmax - 10, y = log(q / m) - 0.5, label = "Bonf."), color = "deepskyblue", size = 4) +
  geom_function(fun = thr_naive, colour = "green") +
  geom_text(data = data[1, ], aes(x = 10, y = log(q) + 0.5, label = "Naive"), color = "green", size = 4) +
  geom_function(fun = thr_holm, color = "red", lty = "dashed") +
  geom_text(data = data[1, ], aes(x = xmax - 10, y = log(q / (m - xmax + 1)) + 0.5, label = "Holm"), color = "red", size = 4) +
  geom_vline(xintercept = TOTP_BH, lty = "dotted", alpha = 0.25) +
  geom_vline(xintercept = TOTP_bonf, lty = "dotted", alpha = 0.25) +
  geom_vline(xintercept = TOTP_holm, lty = "dotted", alpha = 0.25) +
  geom_vline(xintercept = TOTP_naive, lty = "dotted", alpha = 0.25) +
  geom_text(data = data[1, ], aes(x = TOTP_bonf, y = -15, label = as.character(TOTP_bonf)), color = "red", size = 4) +
  geom_text(data = data[1, ], aes(x = TOTP_holm, y = -15, label = as.character(TOTP_holm)), color = "red", size = 4) +
  geom_text(data = data[1, ], aes(x = TOTP_BH, y = -15, label = as.character(TOTP_BH)), color = "red", size = 4) +
  geom_text(data = data[1, ], aes(x = TOTP_naive, y = -15, label = as.character(TOTP_naive)), color = "red", size = 4) +
  labs(x = "index j", y = "(p-value) log scale", color = NULL, title = "p-values vs rank: Procedures Comparison(log scale)") +
  scale_x_continuous(breaks = c(seq(0, xmax, by = 50)), minor_breaks = NULL, ) +
  scale_y_continuous(breaks = ln_breaks, labels = ln_labels)
plot2
