source("r/automatic_sim_noninferior.R")

for(i in 1:13) {
  p <- c(rep(0.75, 13 - i + 1), rep(0.85, i - 1))
  res <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
  saveRDS(res, paste0("sim_out/", sprintf("%02d", i - 1), "_effect_noninf_scenario1.rds"))
}

for(i in 1:13) {
  p <- c(rep(0.75, 13 - i + 1), rep(0.8, i - 1))
  res <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
  saveRDS(res, paste0("sim_out/", sprintf("%02d", i - 1), "_effect_noninf_scenario2.rds"))
}

for(i in 1:13) {
  p <- c(rep(0.8, 13 - i + 1), rep(0.9, i - 1))
  res <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
  saveRDS(res, paste0("sim_out/", sprintf("%02d", i - 1), "_effect_noninf_scenario3.rds"))
}

for(i in 1:13) {
  p <- c(rep(0.8, 13 - i + 1), rep(0.85, i - 1))
  res <- simulate_scenario_noninf_par(
    sims = 1000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
  saveRDS(res, paste0("sim_out/", sprintf("%02d", i - 1), "_effect_noninf_scenario4.rds"))
}
