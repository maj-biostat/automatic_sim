source("r/automatic_sim_functions.R")

rho <- 1/3
gamma <- 1/4
sigsq <- 1.5^2
Sigma0 <- matrix(0, 13, 13)
diag(Sigma0) <- 1
Sigma0[rbind(
  c(2,3), c(2,4), c(3,4), c(3,2), c(4,2), c(4,3),
  c(5,6), c(5,7), c(6,7), c(6,5), c(7,5), c(7,6),
  c(8,9), c(8,10), c(9,10), c(9,8), c(10,8), c(10,9),
  c(11, 12), c(11, 13), c(12, 13), c(12, 11), c(13, 11), c(13, 12)
)] <- rho
Sigma0[rbind(
  c(2, 5), c(2, 8), c(2, 11), c(5,2), c(8,2), c(11,2),
  c(3, 6), c(3, 9), c(3, 12), c(6, 3), c(9, 3), c(12, 3),
  c(4, 7), c(4, 10), c(4, 13), c(7, 4), c(10, 4), c(13, 4),
  c(5, 8), c(5, 11), c(11, 5), c(8, 5),
  c(6, 9), c(6, 12), c(12, 6), c(9, 6),
  c(7, 10), c(7, 13), c(13, 7), c(10, 7),
  c(8, 11), c(11, 8), c(9, 12), c(12, 9), c(10, 13), c(13, 10)
)] <- gamma
Sigma0 <- sig * Sigma0

sim_scenarios <- list(
  null1 = list(
    p = rep(0.75, 13),
    equ_thres = 0.5,
    sup_thres = 0.85,
    drop_thres_best = 0.01,
    fix_ctrl = 0.05,
    Sigma0 = Sigma
  )
)



o <- simulate_trial_basic(
  nA = 3, nB = 4, n_max = 1e4, burn_in = 1000, n_update = 1000, sup_thres = 0.9, drop_thres_best = 0.001,
  allow_stopping = T, brar = 3, fix_ctrl = NULL, Sigma0 = Sigma0, eps = 0.2, equ_thres = 0.5)
o <- simulate_trial_basic(
  nA = 3, nB = 4, p = c(rep(0.75, 10), rep(0.85, 3)), n_max = 1e4, burn_in = 1000, n_update = 1000,
  sup_thres = 0.9, drop_thres_best = 0.001,
  allow_stopping = T, brar = 3, fix_ctrl = NULL, Sigma0 = Sigma0, eps = 0.2, equ_thres = 0.5)

fin_mu <- tail(o$mu, 1)[1, ]
fin_sig <- tail(o$sigma, 1)[1, ]
for(i in 1:5) {
  if(i == 1) {
    curve(dnorm(x, fin_mu[1], sqrt(fin_sig[1])), xlim = c(-1, 2), ylim = c(0, 10))
  }
  curve(dnorm(x, fin_mu[i], sqrt(fin_sig[i])), col = i, add = TRUE, n = 1e3)
}
dd <- plogis(rmvnorm(1e4, fin_mu, diag(fin_sig)))
for(i in 1:5) {
  if(i == 1) {
    plot(density(dd[, 1]), xlim = c(0.5, 1))
  }
  lines(density(dd[ ,i]), col = i)
}

m <- c(0.98,0.99,1,1,1.01)
S <- diag(0.001,5)
A <- cbind(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), rep(-1,4))
1 - pnorm(-0.1, m[1] - m[5], sqrt(S[1,1] + S[5,5]))
pmvnorm(-0.1, Inf, m[-5] - m[5], sigma = A %*% S %*% t(A))
D <- rmvnorm(1e4, m, S)
prob_all_noninferior(D[,-5], D[,5], 0.1)
prob_noninferior_to_best(D[, -5], D[, 5], 0.1)

which.max(fin_mu)
A <- matrix(c(-1, -1, -1, -1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 , 0, 0, 0, 0, 1), 4, 5)
M <- A %*% fin_mu
S <- A %*% diag(fin_sig) %*% t(A)
1 - pnorm(-0.1, fin_mu[5] - fin_mu[1], sqrt(fin_sig[1] + fin_sig[5]))
pmvnorm(rep(-0.1, 4), rep(Inf, 4), M[, 1], sigma = S)

p0 <- rep(0.75, 13)
stpar <- system.time(
  res0 <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01,
    allow_stopping = T, brar = 3, fix_ctrl = 0.05, Sigma0 = Sigma0, eps = 0.2))
summarise_scenario_sim(res0)
saveRDS(res0, "sim_out/sim2/null_effect_scenario1.rds")

p0 <- rep(0.75, 13)
stpar <- system.time(
  res0 <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.001,
    allow_stopping = T, brar = 3, fix_ctrl = 0.05, Sigma0 = Sigma0, eps = 0.2))
summarise_scenario_sim(res0)
saveRDS(res0, "sim_out/sim2/null_effect_scenario2.rds")

p <- c(rep(0.75, 12), rep(0.85, 1))
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.9, drop_thres_best = 0.01,
    allow_stopping = T, brar = 3, fix_ctrl = 1 / length(p), Sigma0 = Sigma0, eps = 0.2))
saveRDS(res, "sim_out/sim2/large_effect_one_sup1.rds")
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.9, drop_thres_best = 0.001,
    allow_stopping = T, brar = 3, fix_ctrl = 1 / length(p), Sigma0 = Sigma0, eps = 0.2))
saveRDS(res, "sim_out/sim2/large_effect_one_sup2.rds")

p <- c(rep(0.75, 10), rep(0.85, 3))
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.9, drop_thres_best = 0.01,
    allow_stopping = T, brar = 3, fix_ctrl = 1 / length(p), Sigma0 = Sigma0, eps = 0.2))
saveRDS(res, "sim_out/sim2/large_effect_three_sup1.rds")
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.9, drop_thres_best = 0.001,
    allow_stopping = T, brar = 3, fix_ctrl = 1 / length(p), Sigma0 = Sigma0, eps = 0.2))
saveRDS(res, "sim_out/sim2/large_effect_three_sup2.rds")

p <- c(rep(0.75, 8), rep(0.85, 5))
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.9, drop_thres_best = 0.01,
    allow_stopping = T, brar = 3, fix_ctrl = 1 / length(p), Sigma0 = Sigma0, eps = 0.2))
saveRDS(res, "sim_out/sim2/large_effect_five_sup1.rds")
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    equ_thres = 0.5, sup_thres = 0.9, drop_thres_best = 0.001,
    allow_stopping = T, brar = 3, fix_ctrl = 1 / length(p), Sigma0 = Sigma0, eps = 0.2))
saveRDS(res, "sim_out/sim2/large_effect_five_sup2.rds")
