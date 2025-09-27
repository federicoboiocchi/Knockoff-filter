# KNOCKOFFS: CORE CODE

rm(list = ls())
library(MASS) # to load the function mvrnorm
library(glmnet)
library(knockoff)
# artificial generation of the design matrix

# Setting: n >= 2p
n <- 1500 # number of observations
p <- 80 # number of variables
rho <- 0.3 # correlation among variables
tp <- 12 # number of true positives
A <- 3 # signal amplitude
q <- 0.15 # fdr upper bound
ind_tp <- sample(x = c(1:p), size = tp, replace = FALSE)
ind_fp <- c(1:p)[-ind_tp]

mu <- rep(0, p)
Sigma <- matrix(data = rep(rho, p^2), nrow = p, ncol = p) + diag(rep(1 - rho, p))
X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
Xc <- scale(X, center = T, scale = FALSE) # centered matrix
Xcn <- apply(Xc, 2, function(x) x / sqrt(sum(x^2)))

# True positive variables

Xtp <- Xcn[, ind_tp]
z <- rnorm(n, mean = 0, sd = 1)

y <- A * rowSums(Xtp) + z
y <- (y - mean(y)) / sd(y)

# Equicorrelated knockoffs

I <- diag(p) # identity matrix
SIG <- crossprod(Xcn) # Gram Matrix
eig_min <- min(abs(eigen(SIG)$values))
s <- min(2 * eig_min, 1)
diags <- diag(rep(s, p))

# Inverse gram matrix
invSIG <- solve(SIG)
mat <- 2 * diags - diags %*% invSIG %*% diags + diag(1e-12, nrow = p)

# Cholesky decomposition to find C
C <- chol(mat)
# U matrix (orthonormal and orthogonal to the span of X)
U <- Null(Xcn)[, 1:p]

# apply(U,2,function(x) norm(x,type="2")) # columns are normalized
# (t(U)%*%Xcn) it can be verified it is a null matrix

# knockoff formula satisfying the two constraint on the correlation structure
Xtil <- Xcn %*% (I - invSIG %*% diags) + U %*% C
Xtilc <- scale(Xtil, center = T, scale = FALSE) # centered matrix
Xtilcn <- apply(Xtilc, 2, function(x) x / sqrt(sum(x^2)))
Xtot <- as.matrix(cbind(Xcn, Xtilcn))

nlambda <- 2 * p # number of penalization parameters to consider
lambda_max <- max(abs(crossprod(Xcn, y))) / n # maximum lambda at which all regression coefficients go to zero

lambda_min <- lambda_max / 2000 # (arbitrarily chosen)
k <- (0:(nlambda - 1)) / nlambda
lambda_val <- lambda_max * (lambda_min / lambda_max)^k 
# in this way we have a logarithmic spaced penalization parameter sequence
# between lambda_min and lambda_max, namely lambdas are more densely
# packed near zero and gradually become sparser as they grow.

# fitting the lasso path through glmnet to find the statistics Z_j and
# Z_j tilde.
fit <- glmnet(Xtot, y,
  alpha = 1,
  lambda = lambda_val,
  standardize = FALSE,
  standardize.response = FALSE,
  intercept = FALSE
)
first_nz <- function(x) match(T, abs(x) > 0)
first_nz_ind <- apply(fit$beta, 1, first_nz)
sum(is.na(first_nz_ind))
Z_j <- as.numeric(ifelse(is.na(first_nz_ind), 0, fit$lambda[first_nz_ind]) * n)

# compute the statistics W_j's
W_j <- numeric(p)
ind_orig <- 1:p
W_j <- pmax(Z_j[ind_orig], Z_j[ind_orig + p]) * sign(Z_j[ind_orig] - Z_j[ind_orig + p])

# Compute the data-dependent threshold Th
W <- unique(abs(W_j))
W <- W[W != 0]
FDP <- numeric(length(W))
for (j in 1:length(W)) {
  t <- W[j]
  FDP[j] <- (sum(W_j <= -t)) / max(sum(W_j >= t), 1)
}
j <- which(FDP <= q)
if (length(j) == 0) {
  Th <- Inf
  ind_knock <- integer(0)
} else {
  Th <- min(W[j])
  ind_knock <- which(W_j >= Th)
}

ind_knock <- which(W_j >= Th) # selected indices by the knockoff method
TOTP_knockoff <- length(ind_knock) # total number of discoveries
TP_knockoff <- sum(ind_knock %in% ind_tp) # number of true discoveries
FP_knock <- sum(ind_knock %in% ind_fp) # number of false discoveries
(FDP_knock <- FP_knock / (TOTP_knockoff + 1 / q))
(PW_iter <- TP_knockoff / tp)

# Knockoff pair plot (Figure 3.5)
lim <- max(Z_j) * 8 / 7
par(pty = "s")
plot(Z_j[ind_fp], Z_j[ind_fp + p],
  pch = 19, asp = 1,
  xlim = c(0, lim),
  ylim = c(0, lim),
  xlab = expression(Z[j]),
  ylab = expression(tilde(Z)[j]),
  cex = 0.8,
  main = expression(paste("Knockoff pairs: (", Z[j], ", ", tilde(Z)[j], ")"))
)
abline(a = 0, b = 1, lty = "dashed", col = "grey50", lwd = 1.5)
segments(x0 = Th, y0 = 0 - 5, x1 = Th, y1 = Th, col = "black", lwd = 1.5)
segments(x0 = 0 - 5, y0 = Th, x1 = Th, y1 = Th, lwd = 1.5)
points(Z_j[ind_tp], Z_j[ind_tp + p], col = "red", pch = 15, cex = 0.8)
xx <- c(-Th, -Th, Th, Th, lim * 3 / 2, lim * 3 / 2)
yy <- c(lim * 3 / 2, Th, Th, 0 - Th, 0 - Th, lim * 3 / 2)
polygon(xx, yy, border = NULL, col = rgb(0.5, 0.5, 0.5, alpha = 0.1))
legend("topright",
  inset = c(-0.58, 0),
  legend = c("Null features", "Non-null features"),
  col = c("black", "red"),
  pch = c(19, 15),
  xpd = NA
)

# Testing the Code programmed from scratch with the results of the actual
# knockoff package of Candes and Barber. We needed a
# customized version in order to adapt it to the construction shown in the article.
knock <- function(X) create.fixed(Xcn, method = c("equi"), randomize = F)
stats <- function(X, X_k, y, nlambda, standardize) {
  stat.lasso_lambdasmax(
    X = Xcn, X_k = Xtil, y = y, nlambda = 2 * p, ,
    standardize = FALSE
  )
}
result <- knockoff.filter(Xcn, y, knockoffs = knock, statistic = stats, fdr = q, offset = 0)
result # vector of indices selected by knockoff procedure implemente in the package knockoff
ind_knock # vector of indices selected by the code I've programmed

# result and ind_knock are pretty much always identical vectors.
# there could still be differences since they used a computationally more efficient
# way of constructing knockoffs (using SVD), while this code is less robust
# to errors.
Wstats <- stats(Xcn, Xtil, y, nlambda = 2 * p)
W_j
Wstats
# W_j, Wstats are the vector of test statistics pretty much always identical
# for the same reasons explained above.
thresh <- knockoff.threshold(Wstats, fdr = q, offset = 0)
thresh
Th # these are the final data dependent threshold
# pretty much always identical

# the following are the essential commands to navigate and inspect
# the functions of the Knockoff package of Candes (only those that were
# interesting for this work)
# and understand its structure
# getAnywhere(stat.lasso_lambdasmax)
# getAnywhere(stat.glmnet_lambdasmax)
# getAnywhere(lasso_max_lambda)
# getAnywhere(lasso_max_lambda_glmnet)
# getAnywhere(create.fixed)
# getAnywhere(create_equicorrelated)
# get("decompose", envir = asNamespace("knockoff"))
# getAnywhere(create.fixed)
