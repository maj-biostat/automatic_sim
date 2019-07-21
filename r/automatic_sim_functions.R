

# nA - number of levels in factor A
# nB - number of levels in factor B
# p  - cell response probability
# n_max - maximum sample size
# burn_in - timing of first interim
# n_update - timing of subsequent interims
# sup_thres - threshold for declaring superiority
# equ_thres - threshold for declaring non-inferiority relative to best
# drop_thres_best - probability for setting allocation to zero
simulate_trial_basic <- function(
  nA = 3,
  nB = 4,
  p = rep(0.75, nA * nB + 1),
  n_max = 1e4,
  burn_in = 1e3,
  n_update = 1e3,
  sup_thres = 0.95,
  equ_thres = 0.95,
  drop_thres_best = 0.01,
  drop_thres_ctrl = drop_thres_best,
  eps = 0.1,
  brar = 3,
  rand = 'mwu',
  fix_ctrl = NULL,
  allow_stopping = TRUE,
  allow_equiv_stop = TRUE,
  allow_noninf_stop = TRUE,
  min_alloc = 0,
  mu0 = rep(0, nA * nB + 1),
  Sigma0 = diag(1, nA * nB + 1),
  ...
) {
  # Argument checks
  if(length(p) != (1 + nA * nB)) stop("p size doesn't match nA + nB + 1")
  if(!(all(p > 0) & all(p < 1))) stop("all p must be between 0 and 1.")
  if(burn_in > n_max | n_update > n_max) stop("burn_in and n_update should be less than n_max.")
  if(!(sup_thres > 0 & sup_thres < 1)) stop("sup_thres must be between 0 and 1.")
  if(!(drop_thres_best >= 0 & drop_thres_best <= 1)) stop("drop_thres_best must be between 0 and 1.")
  if(!(drop_thres_ctrl >= 0 & drop_thres_ctrl <= 1)) stop("drop_thres_ctrl must be between 0 and 1.")
  if(!brar %in% c(0, 1,2,3)) {
    warning("brar was neither 0, 1, 2, or 3. Value was set to 1 by default")
    brar <- 1
  }
  if(!rand %in% c("mwu", "simple")) {
    warning("rand was neither 'mwu' nor 'simple'. Value set to 'mwu' by default.")
    rand <- "mwu"
  }

  # Design
  n_arms <- length(p)
  a <- c(0, rep(1:nA, nB))
  b <- c(0, rep(1:nB, each = nA))
  A <- cbind(0, Matrix::bdiag(replicate(nA, list(1)))[, rep(1:nA, nB)]/nB)
  B <- cbind(0, Matrix::bdiag(replicate(nB, list(t(rep(1, nA)))))/nA)
  Xdes <- diag(n_arms)

  # Setup interim sample sizes
  n_seq <- unique(c(seq(burn_in, n_max, by = n_update), n_max))
  n_interims <- length(n_seq)
  n_draw <- diff(c(0, n_seq))

  # Data storage
  alloc_prob <- matrix(0, n_interims + 1, n_arms)
  if(is.null(fix_ctrl)) {
    alloc_prob[1, ] <- rep(1/n_arms, n_arms)
  } else {
    alloc_prob[1, ] <- c(fix_ctrl, rep((1 - fix_ctrl) / (n_arms - 1), n_arms - 1))
  }

  n <- matrix(0, n_interims, n_arms)
  y <- matrix(0, n_interims, n_arms)
  cum_n <- rep(0, n_arms)
  cum_y <- rep(0, n_arms)

  p_best <- matrix(0, n_interims, n_arms)
  # p_best_A <- matrix(0, n_interims, nA)
  # p_best_B <- matrix(0, n_interims, nB)
  p_better_than_ctrl <- matrix(0, n_interims, n_arms - 1)
  p_noninferior_to_best <- matrix(0, n_interims, n_arms)
  p_superior <- rep(0, n_interims)
  p_noninferior <- rep(0, n_interims)
  p_equivalent <- rep(0, n_interims)

  i_best <- rep(0, n_interims)
  i_equ <- NULL
  i_sup <- NULL
  i_non <- NULL

  is_best <- matrix(0, n_interims, n_arms)
  is_best_A <- matrix(0, n_interims, nA)
  is_best_B <- matrix(0, n_interims, nB)
  is_worse_than_ctrl <- matrix(0, n_interims, n_arms - 1)
  is_noninferior_to_best <- matrix(0, n_interims, n_arms)
  is_competing <- matrix(1, n_interims + 1, n_arms)
  n_competing <- rep(0, n_interims)

  superior <- FALSE
  noninferior <- FALSE
  equivalent <- FALSE

  # Posterior parameters
  mu <- matrix(0, n_interims, n_arms)
  sigma <- matrix(0, n_interims, n_arms)

  for(i in 1:n_interims) {
    if(rand == "mwu") {
      mwu <- mass_weighted_urn_design(alloc_prob[i, ], n_draw[i], ...)
      n[i, ] <- mwu$sample_size[nrow(mwu$sample_size), ]
    } else if (rand == "simple") {
      alloc <- sample(1:n_arms, n_draw[i], alloc_prob[i, ], replace = T)
      n[i, ] <- as.numeric(table(factor(alloc, levels = 1:n_arms)))
    }

    y[i, ] <- rbinom(n_arms, n[i, ], p)
    cum_n <- cum_n + n[i, ]
    cum_y <- cum_y + y[i, ]

    mod_out <-
      est_model_approx(
        Xdes, cum_y, cum_n,
        mu0 = mu0, Sigma0 = Sigma0,
        mu_init = rep(0, ncol(Xdes)), Sigma_init = diag(ncol(Xdes))
      )

    mu[i, ] <- mod_out$mu
    sigma[i, ] <- diag(mod_out$sigma)

    # hdi[, , i] <- mapply(function(x,y)
    #   qnorm(c(0.25, 0.975), x, y), mu[i, ], sqrt(sigma[i, ]))
    #
    p_best[i, ] <- prob_best(mod_out$post_draws)
    i_best[i] <- which.max(mod_out$mu)
    # p_best_A[i, ] <- prob_in_best(mod_out$post_draws, A)
    # p_best_B[i, ] <- prob_in_best(mod_out$post_draws, B)
    p_better_than_ctrl[i, ] <- control_comp(mod_out$post_draws)
    # p_noninferior_to_best[i, -i_best[i]] <-
    #   prob_each_noninferior(mod_out$post_draws[, -i_best[i]],
    #                            mod_out$post_draws[, i_best[i]],
    #                            eps)
    # Best is noninferior to itself
    # p_noninferior_to_best[i, i_best[i]] <- 1

    is_best[i, ] <- p_best[i, ] > sup_thres
    # is_best_A[i, ] <- p_best_A[i, ] > sup_thres
    # is_best_B[i, ] <- p_best_B[i, ] > sup_thres
    is_worse_than_ctrl[i, ] <- p_better_than_ctrl[i, ] < drop_thres_ctrl
    # is_noninferior_to_best[i, ] <- p_noninferior_to_best[i, ] > equ_thres

    # Update allocations
    if(brar == 0) { # No RAR
      alloc_prob[i + 1, ] <- alloc_prob[i, ]
    } else if(brar == 1) {
      alloc_prob[i + 1, ] <- brar1(p_best[i, ])
    } else if (brar == 2) {
      alloc_prob[i + 1, ] <- brar2(p_best[i, ], n_draw[i], n_max)
    } else {
      alloc_prob[i + 1, ] <- brar3(p_best[i, ], cum_n, sigma[i, ], no_alloc_thres = drop_thres_best, fix_ctrl = fix_ctrl, min_alloc = 0)
      # Redistribute allocation if any worse than control
      if(any(is_worse_than_ctrl[i, ] == 1))
        update_alloc(alloc_prob[i, ], which(is_worse_than_ctrl[i, ] == 1) + 1)
    }

    p_superior[i] <- prob_superior(mod_out$post_draws, which.max(mod_out$mu), eps)

    # Are any competing for superiority?
    is_competing[i + 1, ] <- as.numeric(p_best[i, ] > drop_thres_best)
    # How many competing?
    n_competing[i] <- sum(is_competing[i + 1, ])
    if(n_competing[i] > 1) {
      p_noninferior[i] <-
        prob_all_noninferior(mod_out$post_draws[, which(is_competing[i + 1, ] == 1)][, -i_best[i], drop = F],
                             mod_out$post_draws[, i_best[i]],
                             eps)
      p_equivalent[i] <-
        prob_all_equivalent(mod_out$post_draws[, which(is_competing[i + 1, ] == 1)], eps)
    }

    superior <- any(is_best[i, ] == 1)
    noninferior <- p_noninferior[i] > equ_thres
    equivalent <- p_equivalent[i] > equ_thres

    if(superior)
      i_sup <- which(is_best[i, ] == 1)
    if(noninferior)
      i_non <- which(is_competing[i + 1, ] == 1)
    if(equivalent)
      i_equ <- which(is_competing[i + 1, ] == 1)

    if(allow_stopping & i < n_interims) {

      # Do stopping checks for superiority
      # and for equivalence of arms with non-zero allocation
      if(superior || (noninferior & allow_noninf_stop) || (equivalent & allow_equiv_stop)) {

        return(list(
          n = n[1:i, , drop = F],
          y = y[1:i, , drop = F],
          alloc_prob = alloc_prob[1:(i+1), , drop = F],
          p_best = p_best[1:i, , drop = F],
          # p_best_A = p_best_A[1:i, , drop = F],
          # p_best_B = p_best_B[1:i, , drop = F],
          p_better_than_ctrl = p_better_than_ctrl[1:i, , drop = F],
          # p_noninferior_to_best = p_noninferior_to_best[1:i, , drop = F],
          p_superior = p_superior[1:i],
          p_noninferior = p_noninferior[1:i],
          p_equivalent = p_equivalent[1:i],
          i_best = i_best[1:i],
          is_best = is_best[1:i, , drop = F],
          # is_best_A = is_best_A[1:i, , drop = F],
          # is_best_B = is_best_B[1:i, , drop = F],
          is_worse_than_ctrl = is_worse_than_ctrl[1:i, , drop = F],
          # is_noninferior_to_best = is_noninferior_to_best[1:i, , drop = F],
          is_competing = is_competing[1:(i+1), , drop = F],
          mu = mu[1:i, , drop = F],
          sigma = sigma[1:i, , drop = F],
          Sigma = mod_out$sigma,
          # hdi = hdi[, , 1:i, drop = F],
          p = p,
          stopped_early = 1,
          superior = superior,
          noninferior = noninferior,
          equivalent = equivalent,
          i_sup = i_sup, i_non = i_non, i_equ = i_equ
        ))
      }
    }


  }
  return(list(
    n = n,
    y = y,
    alloc_prob = alloc_prob,
    p_best = p_best,
    # p_best_A = p_best_A,
    # p_best_B = p_best_B,
    p_better_than_ctrl = p_better_than_ctrl,
    # p_noninferior_to_best = p_noninferior_to_best,
    p_superior = p_superior,
    p_noninferior = p_noninferior,
    p_equivalent = p_equivalent,
    i_best = i_best,
    is_best = is_best,
    # is_best_A = is_best_A,
    # is_best_B = is_best_B,
    is_worse_than_ctrl = is_worse_than_ctrl,
    # is_noninferior_to_best = is_noninferior_to_best,
    is_competing = is_competing,
    mu = mu,
    sigma = sigma,
    Sigma = mod_out$sigma,
    # hdi = hdi,
    p = p,
    stopped_early = 0,
    superior = superior,
    noninferior = noninferior,
    equivalent = equivalent,
    i_sup = i_sup, i_non = i_non, i_equ = i_equ
  ))
}


simulate_scenario_basic_par <- function(
  sims,
  ...
) {
  ex_fun <- c("simulate_trial_basic", "mass_weighted_urn_design",
              "brar1", "brar2", "brar3", "prob_each_noninferior", "prob_superior",
              "prob_best", "prob_in_best", "est_model_approx", "prob_all_noninferior",
              "control_comp", "update_alloc", "prob_all_equivalent")
  ex_pkg <- c("mvtnorm", "varapproxr", "Matrix")
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  res <- foreach(i = seq_len(sims), .export = ex_fun, .packages = ex_pkg) %dopar% {
    simulate_trial_basic(...)
  }
  return(res)
}



simulate_scenario_basic_par_chunks <- function(
  sims,
  ...
) {
  ex_fun <- c("simulate_trial_basic", "mass_weighted_urn_design",
              "brar1", "brar2", "brar3", "prob_best",
              "prob_in_best", "est_model_approx", "pairwise_comp",
              "prob_noninferior_to_best")
  ex_pkg <- c("mvtnorm", "varapproxr", "Matrix")
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  inds <- split(seq_len(sims), sort(rep_len(seq_len(ncores), sims)))

  res <- foreach(l = seq_along(inds), .export = ex_fun, .packages = ex_pkg) %dopar% {
    in_res <- lapply(inds[[l]], function(i) simulate_trial_basic(...))
    in_res
  }
  return(do.call(c, res))
}

simulate_scenario_basic <- function(
  sims,
  ...
) {
  res <- lapply(seq_len(sims), function(i) simulate_trial_basic(...))
  return(res)
}


summarise_scenario_sim <- function(sim) {
  the_best <- which(sim[[1]][["p"]] == max(sim[[1]][["p"]]))
  if(length(the_best) == 1) {
    prob_identify_correct <-
      mean(unlist(lapply(sim, function(x) {
        competing <- tail(x[["i_best"]], 1)
        identical(the_best, competing) && x[["superior"]]
      })))
  } else {
    prob_identify_correct <-
      mean(unlist(lapply(sim, function(x) {
        competing <- which(tail(x[["is_competing"]], 1) == 1)
        identical(the_best, competing) && x[["equivalent"]]
      })))
  }
  list(
    prob_identify_correct = prob_identify_correct,
    superior =
      mean(sapply(sim, function(x) x[["superior"]])),
    noninferior =
      mean(sapply(sim, function(x) x[["noninferior"]])),
    equivalent =
      mean(sapply(sim, function(x) x[["equivalent"]])),
    stopped_early =
      mean(sapply(sim, function(x) x[["stopped_early"]])),
    stopping_time =
      table(sapply(sim, function(x) nrow(x[["n"]]))),
    expected_ss =
      apply(do.call(rbind, lapply(sim, function(x) sum(x[["n"]]))), 2, mean),
    n_competing =
      table(do.call(rbind, lapply(sim, function(x) sum(tail(x[["is_competing"]], 1))))),
    prob_competing_at_end =
      apply(do.call(rbind, lapply(sim, function(x) tail(x[["is_competing"]], 1))), 2, mean),
    prob_equivalent_at_end =
      apply(do.call(rbind,lapply(sim[unlist(lapply(sim, function(x) x[["equivalent"]]))], function(x) tail(x[["is_competing"]],1))), 2, mean),
    prob_find_effect =
      sum(apply(do.call(rbind, lapply(sim, function(x) x[["is_best"]][nrow(x[["is_best"]]), ])), 2, mean)),
    prob_declare_best =
      apply(do.call(rbind, lapply(sim, function(x) x[["is_best"]][nrow(x[["is_best"]]), ])), 2, mean),
    avg_estimate =
      apply(do.call(rbind, lapply(sim, function(x) tail(x[["mu"]], 1)[1, ])), 2, mean),
    avg_estimate_trans = apply(do.call(rbind, lapply(sim, function(x) plogis(tail(x[["mu"]], 1)[1, ]))), 2, mean),
    avg_allocations = apply(do.call(rbind, lapply(sim, function(x) apply(x[["n"]], 2, sum))), 2, mean),
    avg_num_equiv = mean(unlist(lapply(sim[unlist(lapply(sim, function(x) x[["equivalent"]]))], function(x) length(x[["i_equ"]]))))
  )
}
