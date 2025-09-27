# Comparison of Power and FDR as functions of the feature correlation
# across six procedures: Naive, Bonferroni, Benjamini-Hochberg, Holm, Knockoff and Knockoff+
# (equicorrelated knockoffs and non-orthogonal design)
# design and error and response redrawn at each montecarlo iteration
# vector of betas fixed.

rm(list = ls())
set.seed(321)
par(mfrow = c(1, 1))
library(MASS) # to load the function mvrnorm
library(glmnet)
library(knockoff)

# artificial generation of the design matrix
# n >= 2p
n <- 300 # number of observations
p <- 150 # number of variables
tp <- 10 # number of true positives
A <- 3.5
q <- 0.2 # fdr upper bound
ind_tp <- c(1:tp)
ind_fp <- c((tp + 1):p)
M <- 800 # number of Montecarlo iterations
rho <- seq(0, 0.99, length.out = 30)
FDR_bhq <- numeric(length(rho))
PW_bhq <- numeric(length(rho))
FDR_knock <- numeric(length(rho))
PW_knock <- numeric(length(rho))
FDR_knock_plus <- numeric(length(rho))
PW_knock_plus <- numeric(length(rho))
FDR_naive <- numeric(length(rho))
PW_naive <- numeric(length(rho))
FDR_bonf <- numeric(length(rho))
PW_bonf <- numeric(length(rho))
FDR_holm <- numeric(length(rho))
PW_holm <- numeric(length(rho))
mu <- rep(0, p)
betas <- matrix(c(rep(A, tp), rep(0, p - tp)), ncol = 1)

for (i in 1:length(rho)) {
  FDP_iter_knock <- numeric(M)
  PW_iter_knock <- numeric(M)
  FDP_iter_bhq <- numeric(M)
  PW_iter_bhq <- numeric(M)
  FDP_iter_knock_plus <- numeric(M)
  PW_iter_knock_plus <- numeric(M)
  FDP_iter_naive <- numeric(M)
  PW_iter_naive <- numeric(M)
  FDP_iter_bonf <- numeric(M)
  PW_iter_bonf <- numeric(M)
  FDP_iter_holm <- numeric(M)
  PW_iter_holm <- numeric(M)
  Sigma <- matrix(data = rep(rho[i], p^2), p, p) + diag(1 - rho[i], p)

  for (m in 1:M) {
    X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
    Xc <- scale(X, center = T, scale = FALSE) # centered matrix
    Xcn <- apply(X, 2, function(x) x / sqrt(sum(x^2)))
    y <- Xcn %*% betas + rnorm(n)
    y <- y - mean(y)

    # equicorrelated knockoffs
    I <- diag(p) # identity matrix
    SIG <- crossprod(Xcn) # Gram matrix

    eig_min <- min(eigen(SIG)$values)
    s <- min(2 * eig_min, 1)
    diags <- diag(rep(s, p))

    # Inverse gram matrix
    invSIG <- solve(SIG)
    mat <- 2 * diags - diags %*% invSIG %*% diags + diag(1e-12, nrow = p)

    # Cholesky decomposition to find C
    C <- chol(mat)
    U <- Null(Xcn)[, 1:p] # U matrix orthogonal to the span of X
    Xtil <- Xcn %*% (I - invSIG %*% diags) + U %*% C # knockoff matrix
    Xtilc <- scale(Xtil, center = T, scale = FALSE) # centered matrix
    Xtilcn <- apply(Xtilc, 2, function(x) x / sqrt(sum(x^2)))
    Xtot <- as.matrix(cbind(Xcn, Xtilcn))

    nlambda <- 2 * p
    lambda_max <- max(abs(crossprod(X, y))) / n
    lambda_min <- lambda_max / 2000
    k <- (0:(nlambda - 1)) / nlambda
    lambda_val <- lambda_max * (lambda_min / lambda_max)^k

    # Lasso path to compute statistics Z_j and Z_j tilde
    fit <- glmnet(Xtot, y,
      alpha = 1,
      lambda = lambda_val,
      standardize = FALSE,
      standardize.response = FALSE
    )
    first_nz <- function(x) match(T, abs(x) > 0)
    first_nz_ind <- apply(fit$beta, 1, first_nz)
    sum(is.na(first_nz_ind))
    Z_j <- as.numeric(ifelse(is.na(first_nz_ind), 0, fit$lambda[first_nz_ind] * n))

    # compute the statistics W_j's
    W_j <- numeric(p)
    ind_orig <- 1:p
    W_j <- pmax(Z_j[ind_orig], Z_j[ind_orig + p]) * sign(Z_j[ind_orig] - Z_j[ind_orig + p])

    # Compute the data-dependent threshold for Knockoff and Knockoff+
    W <- unique(abs(W_j))
    W <- W[W != 0]
    FDP_plus <- numeric(length(W))
    FDP <- numeric(length(W))
    for (j in 1:length(W)) {
      t <- W[j]
      FDP_plus[j] <- (sum(W_j <= -t) + 1) / max(sum(W_j >= t), 1)
      FDP[j] <- sum(W_j <= -t) / max(sum(W_j >= t), 1)
    }

    j_plus <- which(FDP_plus <= q)
    j <- which(FDP <= q)

    if (length(j) == 0) {
      Th <- Inf
      ind_knock <- integer(0)
    } else {
      Th <- min(W[j])
      ind_knock <- which(W_j >= Th)
    }

    ind_knock <- which(W_j >= Th)
    TOTP_knockoff <- length(ind_knock)
    TP_knockoff <- sum(ind_knock %in% ind_tp)
    FP_knock <- sum(ind_knock %in% ind_fp)
    FDP_iter_knock[m] <- FP_knock / max(TOTP_knockoff, 1)
    PW_iter_knock[m] <- TP_knockoff / tp

    if (length(j_plus) == 0) {
      Th_plus <- Inf
      ind_knock_plus <- integer(0)
    } else {
      Th_plus <- min(W[j_plus])
      ind_knock_plus <- which(W_j >= Th_plus)
    }

    ind_knock_plus <- which(W_j >= Th_plus)
    TOTP_knockoff_plus <- length(ind_knock_plus)
    TP_knockoff_plus <- sum(ind_knock_plus %in% ind_tp)
    FP_knock_plus <- sum(ind_knock_plus %in% ind_fp)
    FDP_iter_knock_plus[m] <- FP_knock_plus / max(TOTP_knockoff_plus, 1)
    PW_iter_knock_plus[m] <- TP_knockoff_plus / tp

    # Benjamini-Hochberg
    mod <- lm(y ~ Xcn - 1) # no intercept
    pvalues <- coef(summary(mod))[, 4]
    cutoff <- max(c(0, which(sort(pvalues) <= q * (1:p) / p)))
    ind_bhq <- which(pvalues <= q * cutoff / p)
    TOTP_bhq <- length(ind_bhq)
    TP_bhq <- sum(ind_bhq %in% ind_tp)
    FP_bhq <- sum(ind_bhq %in% ind_fp)
    FDP_iter_bhq[m] <- FP_bhq / max(TOTP_bhq, 1)
    PW_iter_bhq[m] <- TP_bhq / tp

    # Naive
    alpha <- q
    ind_naive <- which(pvalues <= alpha)
    totp_naive <- length(ind_naive)
    tp_naive <- sum(ind_naive %in% ind_tp)
    fp_naive <- totp_naive - tp_naive
    FDP_iter_naive[m] <- fp_naive / max(totp_naive, 1)
    PW_iter_naive[m] <- tp_naive / tp

    # Bonferroni
    ind_bonf <- which(pvalues <= alpha / p)
    totp_bonf <- length(ind_bonf)
    tp_bonf <- sum(ind_bonf %in% ind_tp)
    fp_bonf <- totp_bonf - tp_bonf
    FDP_iter_bonf[m] <- fp_bonf / max(totp_bonf, 1)
    PW_iter_bonf[m] <- tp_bonf / tp

    # Holm
    indices <- c(1:p)
    i0 <- min(which(sort(pvalues) > alpha / (p - indices + 1)))
    ind_holm <- which(pvalues < alpha / (p - i0 + 1))
    totp_holm <- length(ind_holm)
    tp_holm <- sum(ind_holm %in% ind_tp)
    fp_holm <- totp_holm - tp_holm
    FDP_iter_holm[m] <- fp_holm / max(totp_holm, 1)
    PW_iter_holm[m] <- tp_holm / tp
  }
  FDR_knock[i] <- mean(FDP_iter_knock)
  PW_knock[i] <- mean(PW_iter_knock)

  FDR_bhq[i] <- mean(FDP_iter_bhq)
  PW_bhq[i] <- mean(PW_iter_bhq)

  FDR_knock_plus[i] <- mean(FDP_iter_knock_plus)
  PW_knock_plus[i] <- mean(PW_iter_knock_plus)

  FDR_naive[i] <- mean(FDP_iter_naive)
  PW_naive[i] <- mean(PW_iter_naive)

  FDR_bonf[i] <- mean(FDP_iter_bonf)
  PW_bonf[i] <- mean(PW_iter_bonf)

  FDR_holm[i] <- mean(FDP_iter_holm)
  PW_holm[i] <- mean(PW_iter_holm)
}

# Comparison of Powers as functions of feature correlation
# across six multiple testing approaches (Figure 3.8)
plot(rho, PW_knock * 100,
  type = "l", col = "darkorange",
  xlab = "Feature correlation", ylab = "Power(%)", main = "Power",
  ylim = c(-2, max(PW_naive * 100)), lwd = 1.2, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9
)
points(rho, PW_bhq * 100, type = "l", col = "blue", lwd = 1.2)
points(rho, PW_knock_plus * 100, type = "l", col = "purple", lwd = 1.2)
points(rho, PW_naive * 100, type = "l", col = "green", lwd = 1.2)
points(rho, PW_holm * 100, type = "l", col = "red", lwd = 1.2)
lines(rho, PW_bonf * 100, lty = "dashed", col = "deepskyblue", lwd = 1.2)
legend("topright",
  legend = c("Naive", "BHq", "Holm", "Bonf.", "Knockoff", "Knockoff+"),
  col = c(
    "green", "blue", "red", "deepskyblue",
    "darkorange", "purple"
  ),
  lty = c(1, 1, 1, 2, 1, 1),
  cex = 0.57,
  lwd = 1.1,
)

# Comparison of FDRs as functions of feature correlation
# across six multiple testing approaches (Figure 3.9)
plot(rho, FDR_knock * 100,
  ylim = c(0, 100),
  type = "l", col = "darkorange",
  xlab = "Feature correlation", ylab = "FDR(%)", main = "False Discovery Rate", lwd = 1.2, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9
)
abline(h = q * 100, lty = "dashed", lwd = 1.2, col = "grey50")
points(rho, FDR_bhq * 100, type = "l", col = "blue", lwd = 1.2)
points(rho, FDR_knock_plus * 100, type = "l", col = "purple", lwd = 1.2)
points(rho, FDR_naive * 100, type = "l", col = "green", lwd = 1.2)
points(rho, FDR_holm * 100, type = "l", col = "red", lwd = 1.2)
lines(rho, FDR_bonf * 100, lty = "dashed", col = "deepskyblue", lwd = 1.2)
legend(
  x = 0.75, y = 75,
  legend = c("Naive", "BHq", "Holm", "Bonf.", "Knockoff", "Knockoff+", expression(alpha)),
  col = c(
    "green", "blue", "red", "deepskyblue",
    "darkorange", "purple", "grey50"
  ),
  lty = c(1, 1, 1, 2, 1, 1, 2),
  cex = 0.57,
  lwd = 1.1,
)

