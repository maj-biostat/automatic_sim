source("r/automatic_sim_noninferior.R")

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
Sigma0 <- sigsq * Sigma0

p0 <- rep(0.75, 5)
stpar <- system.time(
  res0 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 2, nB = 2, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.75, sup_thres = 0.95, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/5, Sigma0 = diag(2, 5))
  )
summarise_scenario_noninf(res0)
saveRDS(res0, "sim_out/noninf_5/null_effect_scenario1.rds")

stpar <- system.time(
  res0 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 2, nB = 2, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.75, sup_thres = 0.95, drop_thres_best = 0.001, brar = 3,
    fix_ctrl = 1/5, Sigma0 = diag(2, 5))
)
summarise_scenario_noninf(res0)
saveRDS(res0, "sim_out/noninf_5/null_effect_scenario2.rds")

p1 <- c(rep(0.75, 4), 0.85)
stpar <- system.time(
  res1 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 2, nB = 2, p = p1, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.75, sup_thres = 0.95, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/5, Sigma0 = diag(2, 5))
)
summarise_scenario_noninf(res1)
saveRDS(res1, "sim_out/noninf_5/one_effect_scenario1.rds")

stpar <- system.time(
  res1 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 2, nB = 2, p = p1, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.75, sup_thres = 0.95, drop_thres_best = 0.001, brar = 3,
    fix_ctrl = 1/5, Sigma0 = diag(2, 5))
)
summarise_scenario_noninf(res1)
saveRDS(res1, "sim_out/noninf_5/one_effect_scenario2.rds")


p3 <- c(rep(0.75, 2), rep(0.85, 3))
stpar <- system.time(
  res3 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 2, nB = 2, p = p3, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.75, sup_thres = 0.95, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/5, Sigma0 = diag(2, 5))
)
summarise_scenario_noninf(res3)
saveRDS(res3, "sim_out/noninf_5/three_effect_scenario1.rds")

stpar <- system.time(
  res3 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 2, nB = 2, p = p3, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.75, sup_thres = 0.95, drop_thres_best = 0.001, brar = 3,
    fix_ctrl = 1/5, Sigma0 = diag(2, 5))
)
summarise_scenario_noninf(res3)
saveRDS(res3, "sim_out/noninf_5/three_effect_scenario2.rds")



#===============#
# THIRTEEN ARMS #
#===============#

for(i in 1:13) {
  p <- c(rep(0.75, 13 - i + 1), rep(0.85, i - 1))
  res <- simulate_scenario_noninf_par(
    sims = 2, nA = 3, nB = 4, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
saveRDS(res, paste0("sim_out/noninf_13/", sprintf("%02d", i - 1), "_effect_scenario1.rds"))
}

p0 <- rep(0.75, 13)
stpar <- system.time(
  res0 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
)
summarise_scenario_noninf(res0)
saveRDS(res0, "sim_out/noninf_13/null_effect_scenario1.rds")

p1 <- c(rep(0.75, 12), 0.85)
stpar <- system.time(
  res1 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p1, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res1)
saveRDS(res1, "sim_out/noninf_13/one_effect_scenario1.rds")

p2 <- c(rep(0.75, 11), rep(0.85, 2))
stpar <- system.time(
  res2 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p2, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res2)
saveRDS(res1, "sim_out/noninf_13/two_effect_scenario1.rds")

p3 <- c(rep(0.75, 11), rep(0.85, 2))
stpar <- system.time(
  res2 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p2, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res2)
saveRDS(res1, "sim_out/noninf_13/two_effect_scenario1.rds")

stpar <- system.time(
  res0 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.001, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res0)
saveRDS(res0, "sim_out/noninf_13/null_effect_scenario2.rds")



stpar <- system.time(
  res1 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p1, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.001, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res1)
saveRDS(res1, "sim_out/noninf_13/one_effect_scenario2.rds")


p3 <- c(rep(0.75, 10), rep(0.85, 3))
stpar <- system.time(
  res3 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p3, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res3)
saveRDS(res3, "sim_out/noninf_13/three_effect_scenario1.rds")

stpar <- system.time(
  res3 <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p3, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.001, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(2, 13))
)
summarise_scenario_noninf(res3)
saveRDS(res3, "sim_out/noninf_13/three_effect_scenario2.rds")
