# Author: Federico Boiocchi June 2025

# Experiment on real data: HIV-1 resistance

# Our attenation is focused on HIV 1 resistance to Protease inhibitors, in particular
# to the antiviral NFV (Nelfinavir).
# link for the matrix predictors and responses (matrix of positions of the mutations and drug resistance measurements
# for 7 protease inhibitors) https://hivdb.stanford.edu/_wrapper/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt

# link for the TSM list containing relevant mutation of the HIV-1 protease regardless of the specific protein inhibitor used
# https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/NP_TSM/PI

rm(list = setdiff(ls(), c("PI_data", "PI_TSM")))

# position selected by treatement selected mutation (approximately the ground truth)
pos <- PI_TSM$V1
npos <- length(pos) # number of mutations selected
amm <- PI_TSM$V2

# Structure of PI_TSM: each row is a relevant position in the chain of
# HIV protease that has provably shown mutations. In the second column
# each cell is a list of mutation at a specific position in the Protease chain
# specified in the first column.

# The following is the code to create the full notations (like M46I etc ...)
# for the mutations in table PI_TSM (whose structure is not directly
# M46I, I309 etc...)

# wt is the wild type sequence (not-mutated of HIV-1) protease
# it has been taken from https://www.uniprot.org/uniprotkb/O90777/entry
wt <- unlist(strsplit("PQVTLWQRPIVTIKIGGQLKEALLDTGADDTVLEEMSLPGKWKPKMIGGIGGFIKVRQYDQVSIEICGHKAIGTVLIGPTPVNIIGRNLLTQLGCTLNF", split = ""))
names(wt) <- as.character(1:99)

amm <- strsplit(amm, " ")
amm_unlist <- unlist(amm)
namm <- length(amm_unlist)
wtmut <- list()
k <- 1
for (i in 1:npos) {
  position <- pos[i]
  amm_position <- unlist(amm[[i]])
  nmu <- length(amm_position)
  for (j in 1:nmu) {
    wtmut[[k]] <- paste(c(wt[position], position, amm_position[j]), collapse = "")
    k <- k + 1
  }
}
length(wtmut)
# wtmut contains the full list of mutations in TSM list for HIV-1 protease
# the list is composed by mutation with the notation such as M46I

data <- PI_data
# adjusting the column names
colnames(data) <- data[1, ] # naming the columns with the first row
data <- data[-1, ] # removing the first column of id's
n <- dim(data)[1] # the number of initial rows
# - are missing value that are structural
# the NA in the response drug resistance are important to take into account

summary(is.na(data)) # assessing the presence of NA in the drug resistance measurements
# We choose the drug Nelfinavir (among the Protease inhibitors)
# because it has the least amount of missing values.

# These are the indices of missing values in the response of drug resistance
ind_na_nfv <- which(is.na(data$NFV) == TRUE) # isolating the inidices
length(ind_na_nfv)

# The response is computed as a log-fold change, namely is the log
# base = 10 of the ratio between
# the concentration of drug to inhibit 50% of the replication of the virus
# when the HIV-1 Protease is mutated
# divided by the concentration of drug needed to reduce by 50% the replication of
# the virus when HIV 1 protease is wild type
y <- as.numeric(data$NFV[-ind_na_nfv])
# if y=1 the patient virus is resistant as wild type
# if y>1 patient virus is more resistant than wild type
# if y<1 patient virus is less resistant than wild type
# We also remove from the predictors the rows that have a missing value in the response
# of log-fold change of Nelfinavir
data <- data[-ind_na_nfv, ]
# We remove all responses (in this way dataX is the matrix of predictors)
dataX <- data[, -c(1:10)]
# structure of dataX: we have a matrix approximately 844x99
# where each columns represent a position. in this way we have for a fixed
# position/column all ammino acids mutated at that specific position
# obviously we could have several different mutations at the same position.
# We would like to have specific mutations like M46I instead of positions like P70
# as features. In dataX the value in cell (i,j) is either - if the the HIV-1
# Protease sample of the patient i doesn't show a mutation in position j
# or it is the first letter of the ammino acid that represents the mutation if the mut. is actually present.
# we would like to have only "-" or letters in dataX; We don't want any other symbol.
# therefore we replace "." with "-".
mut_list <- list()
n_X <- dim(dataX)[1]
p_X <- dim(dataX)[2]
anyNA(dataX)
for (j in 1:p_X) {
  for (i in 1:n_X) {
    if (dataX[i, j] == ".") {
      dataX[i, j] <- "-"
    }
  }
}
# In the following for cycles we create the list of all mutations in the standard notation
# that can be found in dataX. So mut_list is a list of characters such as
# M46I, I10L and so on.
k <- 1
for (j in 1:p_X) {
  for (i in 1:n_X) {
    if (dataX[i, j] != "-") {
      mut_list[k] <- paste(c(wt[j], j, dataX[i, j]), collapse = "")
      k <- k + 1
    }
  }
}
# Obviously we will have several mutations that will be exactly equal.
# therefore we take the vector of unique mutations umut.
numut <- length(unique(mut_list))
umut <- unique(mut_list)
# with the following code  we want to create a list of lists.
# more precisely we want a list for each HIV-1 sample including all mutations appearing for that
# sample. so we will have a list of n_X lists, where n_X list is the
# number of filtered patients. each list in the big list will have a different number
# of mutations since different samples have different mutations

mut_grouped_by_sample <- list() # list of list of mutations for each sample of HIV1 protease
for (i in 1:n_X) {
  mut_list_i <- list()
  k <- 1
  for (j in 1:p_X) {
    if (dataX[i, j] != "-") {
      mut_list_i[k] <- paste(c(wt[j], j, dataX[i, j]), collapse = "")
      k <- k + 1
    }
  }
  mut_grouped_by_sample[[i]] <- mut_list_i
}
str(mut_grouped_by_sample)
# now we are able to create a matrix having on the columns unique mutations in standard notations
# and on the rows the HIV-1 Protease samples. The entry in cell (i,j)
# will be either 0 or 1 depending whether sample i-th of HIV-1 Protease contains the mutation
# that labels column j-th or not.

X <- matrix(data = 0, nrow = n_X, ncol = length(umut))
colnames(X) <- umut
for (i in 1:n_X) {
  for (j in 1:length(umut)) {
    X[i, j] <- as.numeric(ifelse(umut[j] %in% unlist(mut_grouped_by_sample[[i]]), 1, 0))
  }
}

# We remove mutations that appears in less than 3 samples (namely we remove columns)
ind_rm <- which(apply(X, 2, sum) < 3)
length(ind_rm)
X <- X[, -ind_rm]

# we also remove duplicates in order to have a full rank matrix
X <- X[, which(!duplicated(t(X)) == T)]

dim(X)
# we standardize the response
y <- (y - mean(y)) / sd(y)
# we now apply the six methods

# Knockoff filter

library(MASS) # to load the function mvrnorm
library(glmnet)
library(knockoff)

n <- nrow(X) # number of observations
p <- ncol(X) # number of variables
q <- 0.2 # fdr upper bound
Xc <- scale(X, center = T, scale = FALSE) # centered matrix
Xcn <- apply(Xc, 2, function(x) x / sqrt(sum(x^2)))
Xcn <- as.matrix(Xcn)
attr(Xcn, "dimnames") <- NULL

# equicorrelated knockoffs
I <- diag(p) # identity matrix
SIG <- crossprod(Xcn) # Gram matrix
anyNA(SIG)
eig_min <- min(eigen(SIG)$values)
s <- min(2 * eig_min, 1)
diags <- diag(rep(s, p))

# inverse gram matrix
invSIG <- solve(SIG)
mat <- 2 * diags - diags %*% invSIG %*% diags + diag(1e-10, nrow = p)

# cholesky decomposition to find C
C <- chol(mat)
U <- Null(Xcn)[, 1:p]
Xtil <- Xcn %*% (I - invSIG %*% diags) + U %*% C
Xtilc <- scale(Xtil, center = T, scale = FALSE) # centered matrix
Xtilcn <- apply(Xtilc, 2, function(x) x / sqrt(sum(x^2)))
Xtot <- as.matrix(cbind(Xcn, Xtilcn))

nlambda <- 2 * p
lambda_max <- max(abs(crossprod(Xcn, y))) / n
nlambda <- 2 * p
lambda_min <- lambda_max / 2000
k <- (0:(nlambda - 1)) / nlambda
lambda_val <- lambda_max * (lambda_min / lambda_max)^k

fit <- glmnet(Xtot, y,
  alpha = 1,
  lambda = lambda_val,
  standardize = FALSE,
  standardize.response = FALSE, intercept = FALSE
)

first_nz <- function(x) match(T, abs(x) > 0)
first_nz_ind <- apply(fit$beta, 1, first_nz)
Z_j <- as.numeric(ifelse(is.na(first_nz_ind), 0, fit$lambda[first_nz_ind] * n))

# compute the statistics W_j's
W_j <- numeric(p)
ind_orig <- 1:p
W_j <- pmax(Z_j[ind_orig], Z_j[ind_orig + p]) * sign(Z_j[ind_orig] - Z_j[ind_orig + p])

# Compute the data-dependent threshold
W <- unique(abs(W_j))
W <- W[W != 0]
FDP <- numeric(length(W))
FDP_plus <- numeric(length(W))
for (j in 1:length(W)) {
  t <- W[j]
  FDP[j] <- (sum(W_j <= -t)) / max(sum(W_j >= t), 1)
  FDP_plus[j] <- (sum(W_j <= -t) + 1) / max(sum(W_j >= t), 1)
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

# ind_knock represent the indices of the columns that have been selected.
# However mutations doesn't correspond directly to columns
# so we need to find the mutations that have been selected
# selected mutations are knock_mut.
mutations <- colnames(X)
knock_mut <- mutations[ind_knock]
# at this point we consider the set of approximately true mutations deriving from the
# TSM list. even though we have for each position the type of mutations
# and we could compare standard notation mutation with the ones selected
# we decide only to compare positions selected and not the complete mutations
# since the list of TSM is approximated and besides that it is general for PIs not
# specific for Nelfinavir.

# we therefore extract the positions of mutations selected
pos_sel_knock <- numeric(length(knock_mut))
for (i in 1:length(knock_mut)) {
  pos_sel_knock[i] <- as.numeric(gsub("\\D", "", knock_mut[i]))
}

# we take unique values of positions selected and we look for
# how many unique positions selected are also in the TSM list of positions
tp_knock_sel <- sum(unique(pos_sel_knock) %in% pos)
tot_knock_sel <- length(unique(pos_sel_knock))
tot_in_TSM <- length(pos)
fp_knock_sel <- tot_knock_sel - tp_knock_sel
# then we repeat the same thing for the other multiple testing methods
if (length(j_plus) == 0) {
  Th_plus <- Inf
  ind_knock_plus <- integer(0)
} else {
  Th_plus <- min(W[j_plus])
  ind_knock_plus <- which(W_j >= Th_plus)
}

ind_knock_plus <- which(W_j >= Th_plus)
TOTP_knockoff_plus <- length(ind_knock_plus)
knock_plus_mut <- mutations[ind_knock_plus]
pos_sel_knock_plus <- numeric(length(knock_plus_mut))
for (i in 1:length(knock_plus_mut)) {
  pos_sel_knock_plus[i] <- as.numeric(gsub("\\D", "", knock_plus_mut[i]))
}
tp_knock_plus_sel <- sum(unique(pos_sel_knock_plus) %in% pos)
tot_knock_plus_sel <- length(unique(pos_sel_knock_plus))
fp_knock_plus_sel <- tot_knock_plus_sel - tp_knock_plus_sel

# knockoff graph code (Figure 3.13)

pos_tot <- numeric()
for (j in 1:p) {
  pos_tot[j] <- as.numeric(gsub("\\D","",mutations[j]))
}
index <- which(pos_tot%in%pos)
par(mfrow=c(1,1))
par(pty="s")
plot(Z_j[ind_orig],Z_j[ind_orig+p],pch=19,
     xlab = expression(Z[j]),
     ylab = expression(tilde(Z)[j]),cex=0.5,xlim=c(0,max(Z_j[ind_orig])),
     ylim=c(0,max(Z_j[ind_orig+p])),main=expression(paste("Knockoff pairs: (", Z[j], ", ", tilde(Z)[j], ")")))
abline(a=0,b=1,lty="dashed",col="grey20",lwd=1)
segments(x0=Th,y0=-Th,x1=Th,y1=Th,col="grey20",lwd=1,lty="dashed")
segments(x0=-Th,y0=Th,x1=Th,y1=Th,col="grey20",lwd=1,lty="dashed")
points(Z_j[index],Z_j[index+p],pch=19,col="red",cex=0.5)

lim<-max(Z_j)
xx1 <- c(Th,Th,6,lim*3/2,lim*3/2) 
yy1<-c(-Th-10,Th,6,lim*3/2,-Th-10)
xx2<-c(-Th-10,-Th-10,Th,6)
yy2<-c(lim*3/2,Th,Th,6)
xx3 <- c(-Th-10,-Th-10,Th,Th)
yy3 <- c(-Th-10,Th,Th,-Th-10)
polygon(xx1, yy1, border = NA, col = rgb(0.2, 0.6, 0.3, alpha = 0.33)) 
polygon(xx2, yy2, border = NA, col = rgb(0.8, 0.3, 0.2, alpha = 0.33))  
polygon(xx3, yy3, border = NA, col = rgb(0.0, 0.0, 1.0, alpha = 0.33))  # muted green
legend("topright",
       inset = c(-0.5, 0), 
       legend = c("Null mutations", "Non-null mutations"),
       col = c("black", "red"),
       pch = c(19, 19),
       xpd = NA,
       cex=0.6)
legend("bottomright",  
       inset = c(-0.52, 0.59),
       legend = c("Selected mut.", "Non selected mut.", "Ignored mut."),  # labels
       fill = c(rgb(0.2, 0.6, 0.3, alpha = 0.33),
                rgb(0.8, 0.3, 0.2, alpha = 0.33),
                rgb(0.0, 0.0, 1.0, alpha = 0.33)),
       xpd = NA,
       cex=0.6) 

# Benjamini-Hochberg
mod <- lm(y ~ Xcn - 1) # no intercept
pvalues <- as.numeric(coef(summary(mod))[, 4])
cutoff <- max(c(0, which(sort(pvalues) <= q * (1:p) / p)))
ind_bhq <- which(pvalues <= q * cutoff / p)
TOTP_bhq <- length(ind_bhq)
bhq_mut <- mutations[ind_bhq]
bhq_mut <- mutations[ind_bhq]
pos_sel_bhq <- numeric(length(bhq_mut))
for (i in 1:length(bhq_mut)) {
  pos_sel_bhq[i] <- as.numeric(gsub("\\D", "", bhq_mut[i]))
}
tp_bhq_sel <- sum(unique(pos_sel_bhq) %in% pos)
tot_bhq_sel <- length(unique(pos_sel_bhq))
fp_bhq_sel <- tot_bhq_sel - tp_bhq_sel

# Naive
alpha <- q
ind_naive <- which(pvalues <= alpha)
TOTP_naive <- length(ind_naive)
naive_mut <- mutations[ind_naive]
naive_mut <- mutations[ind_naive]
pos_sel_naive <- numeric(length(naive_mut))
for (i in 1:length(naive_mut)) {
  pos_sel_naive[i] <- as.numeric(gsub("\\D", "", naive_mut[i]))
}
tp_naive_sel <- sum(unique(pos_sel_naive) %in% pos)
tot_naive_sel <- length(unique(pos_sel_naive))
fp_naive_sel <- tot_naive_sel - tp_naive_sel

# Bonferroni
ind_bonf <- which(pvalues <= alpha / p)
bonf_mut <- mutations[ind_bonf]
pos_sel_bonf <- numeric(length(bonf_mut))
for (i in 1:length(bonf_mut)) {
  pos_sel_bonf[i] <- as.numeric(gsub("\\D", "", bonf_mut[i]))
}
tp_bonf_sel <- sum(unique(pos_sel_bonf) %in% pos)
tot_bonf_sel <- length(unique(pos_sel_bonf))
fp_bonf_sel <- tot_bonf_sel - tp_bonf_sel

# Holm
indices <- c(1:p)
i0 <- min(which(sort(pvalues) > alpha / (p - indices + 1)))
ind_holm <- which(pvalues < alpha / (p - i0 + 1))
TOTP_holm <- length(ind_holm)
holm_mut <- mutations[ind_holm]
pos_sel_holm <- numeric(length(holm_mut))
for (i in 1:length(holm_mut)) {
  pos_sel_holm[i] <- as.numeric(gsub("\\D", "", holm_mut[i]))
}
tp_holm_sel <- sum(unique(pos_sel_holm) %in% pos)
tot_holm_sel <- length(unique(pos_sel_holm))
fp_holm_sel <- tot_holm_sel - tp_holm_sel

# Bar plot to assess the perfomance of Multiple testing procedure
# on HIV data. The following is the code for figure (3.12)

library(tidyverse)
data_comp <- data.frame(
  Procedures = c("Naive", "Bonferroni", "Holm", "BHq", "Knockoff", "Knockoff+"),
  Not_in_TSM_list = c(fp_naive_sel, fp_bonf_sel, fp_holm_sel, fp_bhq_sel, fp_knock_sel, fp_knock_plus_sel),
  In_TSM_list = c(tp_naive_sel, tp_bonf_sel, tp_holm_sel, tp_bhq_sel, tp_knock_sel, tp_knock_plus_sel)
)

data_comp_long <- data_comp %>%
  gather(key = "Type", value = "Count", Not_in_TSM_list, In_TSM_list)

data_comp_long$Type <- recode(data_comp_long$Type,
  "Not_in_TSM_list" = "Not in TSM list",
  "In_TSM_list" = "In TSM list"
)

data_comp_long$Type <- factor(data_comp_long$Type, levels = c("Not in TSM list", "In TSM list"))
colors <- c("Not in TSM list" = "salmon", "In TSM list" = "blue")
ggplot(data_comp_long, aes(x = Procedures, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.7, color = "black") +
  scale_fill_manual(values = colors, name = NULL) +
  labs(
    y = "# HIV-1 protease positions selected",
    title = "Resistance to Nelfinavir"
  ) +
  geom_hline(yintercept = length(pos), lty = "dashed", color = "grey50") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# the grey line represents the number of positions selected in the TSM
# list. each bar represents the total amount of positions selected
# where the blue part are true positives, namely positions selected appearing in the TSM list, while
# the orange part are false positives, namely it represents for each approach the number of
# mutations selected that were not in the TSM list.


