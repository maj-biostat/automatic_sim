source("r/automatic_sim_functions.R")

# null scenario
p <- rep(0.75, 13)
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    10000, nA = 3, nB = 4, p = p, n_max = 1e4, burn_in = 420, n_update = 500, eff_thres = 0.9, 
    allow_stopping = T, brar = 3, fix_ctrl = 0.05, Sigma0 = diag(1.5, 13)))
summarise_scenario_sim(res)
saveRDS(res, "sim_out/null_effect_scenario.rds")

# cell effect scenario
p <- rep(0.75, 13)
p[2] <- 0.85
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    10000, nA = 3, nB = 4, p = p, 
    n_max = 1e4, burn_in = 420, n_update = 500, eff_thres = 0.9,
    allow_stopping = T, brar = 3, fix_ctrl = 0.05, Sigma0 = diag(1.5, 13)))
summarise_scenario_sim(res)
saveRDS(res, "sim_out/cell_effect_scenario.rds")

# level effect scenario
p <- rep(0.75, 13)
p[2:4] <- 0.85
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    10000, nA = 3, nB = 4, p = p, 
    n_max = 1e4, burn_in = 420, n_update = 500, eff_thres = 0.9,
    allow_stopping = T, brar = 3, fix_ctrl = 0.05, Sigma0 = diag(1.5, 13)))
summarise_scenario_sim(res)
saveRDS(res, "sim_out/level_effect_scenario.rds")

# Control is best
p <- rep(0.75, 13)
p[1] <- 0.85
stpar <- system.time(
  res <- simulate_scenario_basic_par(
    10000, nA = 3, nB = 4, p = p, 
    n_max = 1e4, burn_in = 420, n_update = 500, eff_thres = 0.9,
    allow_stopping = T, brar = 3, fix_ctrl = 0.05, Sigma0 = diag(1.5, 13)))
summarise_scenario_sim(res)
saveRDS(res, "sim_out/control_best_scenario.rds")
