# Comparison of Power and FDR as functions of the sparsity level
# across six procedures: Naive, Bonferroni, Benjamini-Hochberg, Holm, Knockoff and Knockoff+
# (equicorrelated knockoffs and non-orthogonal design)
# (design matrix drawn only once, new noise at each montecarlo iteration)

rm(list = ls())
set.seed(321)
par(mfrow = c(1, 1))
library(MASS) # to load the function mvrnorm
library(glmnet) # to compute lasso penalization parameters

# artificial generation of the design matrix
# setting: n >= 2p
n <- 200 # number of observations
p <- 100 # number of variables
tp <- seq(1, 60, by = 1) # number of true positives
A <- 3.5
rho <- 0.4
q <- 0.2 # fdr upper bound
M <- 2000 # number of Montecarlo iterations
ltp <- length(tp)
FDR_bhq <- numeric(ltp)
PW_bhq <- numeric(ltp)
FDR_knock <- numeric(ltp)
PW_knock <- numeric(ltp)
FDR_knock_plus <- numeric(ltp)
PW_knock_plus <- numeric(ltp)
FDR_naive <- numeric(ltp)
PW_naive <- numeric(ltp)
FDR_bonf <- numeric(ltp)
PW_bonf <- numeric(ltp)
FDR_holm <- numeric(ltp)
PW_holm <- numeric(ltp)
Sigma <- matrix(data = rep(rho, p^2), p, p) + diag(1 - rho, p) # covariance matrix of the data
mu <- rep(0, p) # mean vector of the data
X <- mvrnorm(n = n, mu = mu, Sigma = Sigma) # design matrix
Xc <- scale(X, center = T, scale = FALSE) # centered matrix
Xcn <- apply(X, 2, function(x) x / sqrt(sum(x^2))) # normalized design

# equicorrelated knockoffs
I <- diag(p) # identity matrix
SIG <- crossprod(Xcn) # Gram matrix
eig_min <- min(abs(eigen(SIG)$values))
s <- min(2 * eig_min, 1) 
diags <- diag(rep(s, p))
invSIG <- solve(SIG) # inverse gram matrix
mat <- 2 * diags - diags %*% invSIG %*% diags + diag(1e-12, nrow = p)
C <- chol(mat) # cholesky decomposition to find C
U <- Null(Xcn)[, 1:p] # U matrix orthogonal to the span of X
Xtil <- Xcn %*% (I - invSIG %*% diags) + U %*% C # knockoff matrix
Xtilc <- scale(Xtil, center = T, scale = FALSE) # centered matrix
Xtilcn <- apply(Xtilc, 2, function(x) x / sqrt(sum(x^2)))
Xtot <- as.matrix(cbind(Xcn, Xtilcn))
nlambda <- 2 * p

for (i in 1:length(tp)) {
  ind_tp <- c(1:tp[i])
  ind_fp <- c((tp[i] + 1):p)
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
  betas <- matrix(c(rep(A, tp[i]), rep(0, p - tp[i])), ncol = 1)

  for (m in 1:M) {
    eps <- rnorm(n)
    y <- Xcn %*% betas + eps
    y <- y - mean(y)

    lambda_max <- max(abs(crossprod(X, y))) / n
    lambda_min <- lambda_max / 2000
    k <- (0:(nlambda - 1)) / nlambda
    lambda_val <- lambda_max * (lambda_min / lambda_max)^k

    # lasso path to compute statistics Z_j and Z_j tilde
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
    FDP <- numeric(length(W))
    FDP_plus <- numeric(length(W))
    for (j in 1:length(W)) {
      t <- W[j]
      FDP_plus[j] <- (sum(W_j <= -t) + 1) / max(sum(W_j >= t), 1)
      FDP[j] <- sum(W_j <= -t) / max(sum(W_j >= t), 1)
    }

    j <- which(FDP <= q)
    j_plus <- which(FDP_plus <= q)

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
    PW_iter_knock[m] <- TP_knockoff / tp[i]

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
    PW_iter_knock_plus[m] <- TP_knockoff_plus / tp[i]

    # Benjamini-Hochberg
    mod <- lm(y ~ Xcn - 1) # no intercept
    pvalues <- coef(summary(mod))[, 4]
    cutoff <- max(c(0, which(sort(pvalues) <= q * (1:p) / p)))
    ind_bhq <- which(pvalues <= q * cutoff / p)
    TOTP_bhq <- length(ind_bhq)
    TP_bhq <- sum(ind_bhq %in% ind_tp)
    FP_bhq <- sum(ind_bhq %in% ind_fp)
    FDP_iter_bhq[m] <- FP_bhq / max(TOTP_bhq, 1)
    PW_iter_bhq[m] <- TP_bhq / tp[i]

    # Naive
    alpha <- q
    ind_naive <- which(pvalues <= alpha)
    totp_naive <- length(ind_naive)
    tp_naive <- sum(ind_naive %in% ind_tp)
    fp_naive <- totp_naive - tp_naive
    FDP_iter_naive[m] <- fp_naive / max(totp_naive, 1)
    PW_iter_naive[m] <- tp_naive / tp[i]

    # Bonferroni
    ind_bonf <- which(pvalues <= alpha / p)
    totp_bonf <- length(ind_bonf)
    tp_bonf <- sum(ind_bonf %in% ind_tp)
    fp_bonf <- totp_bonf - tp_bonf
    FDP_iter_bonf[m] <- fp_bonf / max(totp_bonf, 1)
    PW_iter_bonf[m] <- tp_bonf / tp[i]

    # Holm
    indices <- c(1:p)
    i0 <- min(which(sort(pvalues) > alpha / (p - indices + 1)))
    ind_holm <- which(pvalues < alpha / (p - i0 + 1))
    totp_holm <- length(ind_holm)
    tp_holm <- sum(ind_holm %in% ind_tp)
    fp_holm <- totp_holm - tp_holm
    FDP_iter_holm[m] <- fp_holm / max(totp_holm, 1)
    PW_iter_holm[m] <- tp_holm / tp[i]
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

# Comparison of Powers as functions of sparsity level
# across six multiple testing approaches (Figure 3.6)
plot(tp, PW_knock * 100,
  type = "l", col = "darkorange",
  xlab = "Sparsity Level", ylab = "Power(%)", main = "Power",
  ylim = c(0, max(PW_naive * 100)), lwd = 1.2, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9
)
points(tp, PW_bhq * 100, type = "l", col = "blue", lwd = 1.2)
points(tp, PW_knock_plus * 100, type = "l", col = "purple", lwd = 1.2)
points(tp, PW_naive * 100, type = "l", col = "green", lwd = 1.2)
points(tp, PW_holm * 100, type = "l", col = "red", lwd = 1.2)
lines(tp, PW_bonf * 100, lty = "dashed", col = "deepskyblue", lwd = 1.2)
legend(
  x = 44.4, y = 46,
  legend = c("Naive", "BHq", "Holm", "Bonf.", "Knockoff", "Knockoff+"),
  col = c(
    "green", "blue", "red", "deepskyblue",
    "darkorange", "purple"
  ),
  lty = c(1, 1, 1, 2, 1, 1),
  cex = 0.6,
  lwd = 1.1,
)

# Comparison of FDRs as functions of sparsity level
# across six multiple testing approaches (Figure 3.7)
plot(tp, FDR_knock * 100,
  ylim = c(0, max(FDR_naive * 100)),
  type = "l", col = "darkorange",
  xlab = "Sparsity Level", ylab = "FDR(%)", main = "False Discovery Rate", lwd = 1.2,
  cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9
)
abline(h = q * 100, lty = "dashed", lwd = 1.2, col = "grey50")
points(tp, FDR_bhq * 100, type = "l", col = "blue", lwd = 1.2)
points(tp, FDR_knock_plus * 100, type = "l", col = "purple", lwd = 1.2)
points(tp, FDR_naive * 100, type = "l", col = "green", lwd = 1.2)
points(tp, FDR_holm * 100, type = "l", col = "red", lwd = 1.2)
lines(tp, FDR_bonf * 100, lty = "dashed", col = "deepskyblue", lwd = 1.2)
legend("topright",
  legend = c("Naive", "BHq", "Holm", "Bonf.", "Knockoff", "Knockoff+", expression(alpha)),
  col = c(
    "green", "blue", "red", "deepskyblue",
    "darkorange", "purple", "grey50"
  ),
  lty = c(1, 1, 1, 2, 1, 1, 2),
  cex = 0.65,
  lwd = 1.1,
)

