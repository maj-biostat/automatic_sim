# automatic_functions.R

library(varapproxr)
library(mystan)
library(mvtnorm)
library(doParallel)
library(Matrix)


# Sim dat
sim_arm_dat <- function(n, p) {
  data.frame(x = 1:length(p), n = n, p = p, y = rbinom(length(p), n, p))
}

# Generate random treatment allocations
# according to mass weighted urn design
mass_weighted_urn_design <- function(
  target_alloc,
  sample_size,
  alpha = 4
) {
  arms <- length(target_alloc)
  prob_alloc <- target_alloc / sum(target_alloc)
  # Masses
  x <- matrix(0, sample_size + 1, arms)
  x[1, ] <- alpha * prob_alloc
  # Sample size
  n <- matrix(0, sample_size + 1, arms)
  # Random number
  y <- runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size + 1, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation Predictability
  g <- rep(0, sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)

  imbalance_cap <- sqrt(sum(((alpha - 1)*(1 - prob_alloc) + (arms - 1))^2))

  for(i in 2:(sample_size + 1)) {
    # Update allocation probabilities
    p[i - 1, ] <- pmax(alpha * prob_alloc - n[i - 1, ] + (i - 1)*prob_alloc, 0)
    p[i - 1, ] <- p[i - 1, ] / sum(p[i - 1, ])
    trt[i-1] <- findInterval(y[i - 1], c(0, cumsum(p[i - 1, ])))
    # Update sample sizes
    n[i, ] <- n[i - 1, ]
    n[i, trt[i-1]] <- n[i, trt[i-1]] + 1
    # Update urn masses
    x[i, trt[i-1]] <- x[i - 1, trt[i-1]] - 1 + prob_alloc[trt[i-1]]
    x[i, -trt[i-1]] <- x[i - 1, -trt[i-1]] + prob_alloc[-trt[i-1]]
    # Calculate imbalance
    d[i - 1] <- sqrt(sum((n[i, ] - (i - 1)*prob_alloc)^2))
    # Calculate allocation predictability
    g[i] <- d[i - 1] / alpha
  }
  return(list(
    max_imbalance_bound = imbalance_cap,
    imbalance = d,
    alloc_predict = g,
    rand_num = y,
    trt = trt,
    mass = x,
    sample_size = n,
    selection_prob = p))
}

# Square root transformed posterior probabilities
# p - probability best
brar1 <- function(p) {
  stopifnot(all(p >= 0))
  r <- sqrt(p)
  return(r / sum(r))
}

# Root transformed proportional to fractional sample size
# p - probability best
brar2 <- function(p, n, N) {
  stopifnot(all(p >= 0))
  r <- sqrt(p ^ (n / N))
  return(r / sum(r))
}

brar3 <- function(p, n, V, no_alloc_thres = 0.01, fix_ctrl = NULL, min_alloc = 0) {
  stopifnot(all(p >= 0))
  stopifnot(all(n > 0))
  m <- length(p)


  # Set min allocation and dropping
  w <- rep(min_alloc, m)
  w[which(p < no_alloc_thres)] <- 0
  w_rem <- 1 - sum(w)
  p[which(p < no_alloc_thres)] <- 0

  r <- sqrt(p * V / n)
  w <- w + w_rem * r / sum(r)

  if(is.null(fix_ctrl)) {
    return(w)
  } else {
    if(!(fix_ctrl > 0 & fix_ctrl < 1)) stop("fix_ctrl must be between 0 and 1.")
    # Re-distribute
    w[1] <- fix_ctrl
    w[-1] <- w[-1] / sum(w[-1]) * (1 - w[1])
    return(w)
  }
}

update_alloc <- function(w, zero_ind) {
  w[zero_ind] <- 0
  w[-zero_ind] <- w[-zero_ind] / sum(w[-zero_ind])
  return(w)
}


# Probability best cell mean based
# on posterior draws
# mat is a matrix of posterior draws
# lvls is the optional levels numbering
prob_best <- function(mat, lvls = NULL) {
  if(is.null(lvls)) lvls <- seq_len(ncol(mat))
  return(prop.table(table(factor(apply(mat, 1, which.max), levels = lvls))))
}

# mat is a matrix of posterior draws
# A   is a design matrix determining average group effects
prob_in_best <- function(mat, A, ...) {
  prob_best(mat %*% t(A))
}

# Rank each row and return probability of each rank
prob_rank <- function(mat) {
  k <- ncol(mat)
  apply(apply(apply(mat, 1, function(x) rank(-x)), 1, function(z) table(factor(z, levels = 1:k))), 2, prop.table)
}

# Probability that each column in mat
# is noninferior to colvec_best
# by the margin eps
prob_each_noninferior <- function(mat, colvec_best, eps) {
  dmat <- sweep(mat, 1, colvec_best, "-")
  return(apply(dmat, 2, function(x) mean(x > -eps)))
}

prob_all_noninferior <- function(mat, colvec_best, eps) {
  dmat <- sweep(mat, 1, colvec_best, "-")
  mean(apply(dmat, 1, function(x) all(x > -eps)))
}

prob_all_equivalent <- function(mat, eps) {
  a <- rep(1/ncol(mat), ncol(mat))
  mean(apply(sweep(mat, 1, drop(mat %*% a), `-`), 1, function(x) all(abs(x) <= eps)))
}

prob_superior <- function(mat, col, eps) {
  mean(apply(mat, 1, function(x) all(x[col] - x[-col] > eps)))
}

# Do pairwise comparison of posterior draws
# cell i,j gives probability that
# parameter j is greater than parameter i
# mat is a matrix of posterior draws
pairwise_comp <- function(mat) {
  out <- apply(mat, 2, function(x) colMeans(x > mat))
  return(out)
}

# Calculate pairwise difference matrix
# under transformation trans
pairwise_diff <- function(mat, trans = I) {
  trans <- match.fun(trans)
  pair_comp <- combn(ncol(mat), 2)
  trans_mat <- trans(mat)
  pair_mat <- apply(pair_comp, 2, function(x) trans_mat[,x[1]] - trans_mat[,x[2]])
  colnames(pair_mat) <- apply(pair_comp, 2, paste, collapse = "-")
  return(pair_mat)
}

# Testing for Euclidean equiality
# between all columns in mat with
# radius delta
prob_equi <- function(mat, delta, ...) {
  dmat <- pairwise_diff(mat, ...)
  mean(apply(dmat, 1, crossprod) < delta^2)
}

# Assumes that control is first column
# Assess probability better than control
control_comp <- function(mat) {
  out <- colMeans(mat[, -1] > mat[, 1])
  return(out)
}

# Use variational Bayes
est_model_approx <- function(X, y, n, draws = 1e4, ...) {
  mod <- vb_logistic_n(X, y, n, alg = "sj", ...)
  post_draws <- rmvnorm(draws, mod$mu[, 1], sigma = mod$Sigma)
  mu = as.numeric(mod$mu)
  return(list(post_draws = post_draws, mu = mu, sigma = mod$Sigma))
}

# Use Stan
est_model <- function(X, y, n, draws = 1e4, ...) {
  mod <- glm_binomial(X, y, n, iter = 2*draws/4, refresh = 0, ...)
  post_draws <- as.matrix(mod)[, 1:length(y)]
  p_best <- prop.table(table(factor(apply(post_draws, 1, which.max), levels = 1:ncol(post_draws))))
  mu <- apply(post_draws, 2, mean)
  sigma <- diag(var(post_draws))
  return(list(post_draws = post_draws, p_best = p_best, mu = mu, sigma = sigma))
}
