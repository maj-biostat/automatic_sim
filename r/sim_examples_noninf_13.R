source("r/automatic_sim_noninferior.R")

for(i in 1:13) {
  p <- c(rep(0.75, 13 - i + 1), rep(0.85, i - 1))
  res <- simulate_scenario_noninf_par(
    sims = 2, nA = 3, nB = 4, p = p0, n_max = 1e4, burn_in = 500, n_update = 500,
    eps = 0.15, non_thres = 0.5, sup_thres = 0.85, drop_thres_best = 0.01, brar = 3,
    fix_ctrl = 1/13, Sigma0 = diag(10, 13))
  saveRDS(res, paste0("sim_out/", sprintf("%02d", i - 1), "_effect_scenario1.rds"))
}
