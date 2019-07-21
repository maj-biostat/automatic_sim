source("r/automatic_sim_misc_functions.R")

simulate_trial_noninf <- function(
  nA = 3,
  nB = 4,
  p = rep(0.75, nA * nB + 1),
  n_max = 1e4,
  burn_in = 1e3,
  n_update = 1e3,
  sup_thres = 0.95,
  non_thres = 0.95,
  drop_thres_best = 0.01,
  drop_thres_ctrl = 0.01,
  eps = 0.1,
  brar = 3,
  rand = 'mwu',
  fix_ctrl = NULL,
  mu0 = rep(0, nA * nB + 1),
  Sigma0 = diag(1, nA * nB + 1),
  ...
) {
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

  # Store posterior probabilities
  p_best <- matrix(0, n_interims, n_arms)
  p_better_than_ctrl <- matrix(0, n_interims, n_arms - 1)
  p_noninferior_to_best <- matrix(0, n_interims, n_arms)
  p_superior <- rep(0, n_interims)
  p_noninferior <- rep(0, n_interims)

  # Store flags
  i_sup <- NULL
  i_non <- NULL
  i_best <- rep(0, n_interims)
  is_best <- matrix(0, n_interims, n_arms)
  is_worse_than_ctrl <- matrix(0, n_interims, n_arms - 1)
  is_noninferior_to_best <- matrix(0, n_interims, n_arms)
  is_competing <- matrix(1, n_interims + 1, n_arms)
  n_competing <- rep(0, n_interims)

  # Stopping flags
  superior <- FALSE
  noninferior <- FALSE

  # Posterior parameters
  mu <- matrix(0, n_interims, n_arms)
  sigma <- matrix(0, n_interims, n_arms)

  # Loop through analyses
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
    p_best[i, ] <- prob_best(mod_out$post_draws)
    i_best[i] <- which.max(mu[i, ])
    p_better_than_ctrl[i, ] <- control_comp(mod_out$post_draws)
    is_best[i, ] <- p_best[i, ] > sup_thres
    is_worse_than_ctrl[i, ] <- p_better_than_ctrl[i, ] < drop_thres_ctrl

    # Update allocations
    if(brar == 0) { # No RAR
      alloc_prob[i + 1, ] <- alloc_prob[i, ]
    } else if(brar == 1) {
      alloc_prob[i + 1, ] <- brar1(p_best[i, ])
    } else if (brar == 2) {
      alloc_prob[i + 1, ] <- brar2(p_best[i, ], n_draw[i], n_max)
    } else {
      alloc_prob[i + 1, ] <- brar3(p_best[i, ], cum_n, sigma[i, ],
                                   no_alloc_thres = drop_thres_best,
                                   fix_ctrl = fix_ctrl, min_alloc = 0)
      # Redistribute allocation if any worse than control
      if(any(is_worse_than_ctrl[i, ] == 1))
        update_alloc(alloc_prob[i, ], which(is_worse_than_ctrl[i, ] == 1) + 1)
    }

    # Probability that the current best
    # arm is superior to all others by
    # at least eps
    p_superior[i] <- prob_superior(mod_out$post_draws, which.max(mod_out$mu), eps)

    # Are any competing for superiority?
    is_competing[i + 1, ] <- as.numeric(p_best[i, ] > drop_thres_best)
    # How many are competing?
    n_competing[i] <- sum(is_competing[i + 1, ])
    # If multiple are competing then
    # check if they are noninferior to
    # the current best
    if(n_competing[i] > 1) {
      p_noninferior[i] <-
        prob_all_noninferior(mod_out$post_draws[, which(is_competing[i + 1, ] == 1)][, -i_best[i], drop = F],
                             mod_out$post_draws[, i_best[i]],
                             eps)
    }

    superior <- any(is_best[i, ] == 1)
    noninferior <- p_noninferior[i] > non_thres

    # Which is superior or which are noninferior
    if(superior)
      i_sup <- which(is_best[i, ] == 1)
    if(noninferior)
      i_non <- which(is_competing[i + 1, ] == 1)

    if(i < n_interims) {

      # Do stopping checks for superiority
      # and for equivalence of arms with non-zero allocation
      if(superior || noninferior) {

        return(list(
          n = n[1:i, , drop = F],
          y = y[1:i, , drop = F],
          alloc_prob = alloc_prob[1:(i+1), , drop = F],
          p_best = p_best[1:i, , drop = F],
          p_better_than_ctrl = p_better_than_ctrl[1:i, , drop = F],
          p_superior = p_superior[1:i],
          p_noninferior = p_noninferior[1:i],
          i_best = i_best[1:i],
          is_best = is_best[1:i, , drop = F],
          is_worse_than_ctrl = is_worse_than_ctrl[1:i, , drop = F],
          is_competing = is_competing[1:(i+1), , drop = F],
          mu = mu[1:i, , drop = F],
          sigma = sigma[1:i, , drop = F],
          Sigma = mod_out$sigma,
          p = p,
          stopped_early = 1,
          superior = superior,
          noninferior = noninferior,
          i_sup = i_sup, i_non = i_non
        ))
      }
    }
  }
  return(list(
    n = n,
    y = y,
    alloc_prob = alloc_prob,
    p_best = p_best,
    p_better_than_ctrl = p_better_than_ctrl,
    p_superior = p_superior,
    p_noninferior = p_noninferior,
    i_best = i_best,
    is_best = is_best,
    is_worse_than_ctrl = is_worse_than_ctrl,
    is_competing = is_competing,
    mu = mu,
    sigma = sigma,
    Sigma = mod_out$sigma,
    p = p,
    stopped_early = 0,
    superior = superior,
    noninferior = noninferior,
    i_sup = i_sup, i_non = i_non
  ))
}


simulate_scenario_noninf_par <- function(
  sims,
  ...
) {
  ex_fun <- c("simulate_trial_noninf", "mass_weighted_urn_design",
              "brar1", "brar2", "brar3", "prob_each_noninferior", "prob_superior",
              "prob_best", "prob_in_best", "est_model_approx", "prob_all_noninferior",
              "control_comp", "update_alloc", "prob_all_equivalent")
  ex_pkg <- c("mvtnorm", "varapproxr", "Matrix")
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  res <- foreach(i = seq_len(sims), .export = ex_fun, .packages = ex_pkg) %dopar% simulate_trial_noninf(...)
  return(res)
}

summarise_scenario_noninf <- function(sim) {
  the_best <- which(sim[[1]][["p"]] == max(sim[[1]][["p"]]))
  any_noninf <- any(sapply(sim, function(x) x[["noninferior"]]))
  if(any_noninf) {
    avg_num_noninf <- mean(unlist(lapply(sim[unlist(lapply(sim, function(x) x[["noninferior"]]))], function(x) length(x[["i_non"]]))))
    prob_noninf_at_end <-
      apply(do.call(rbind,lapply(sim[unlist(lapply(sim, function(x) x[["noninferior"]]))], function(x) tail(x[["is_competing"]],1))), 2, mean)
  } else {
    prob_noninf_at_end <- NA
    avg_num_noninf <- NA
  }
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
        identical(the_best, competing) && x[["noninferior"]]
      })))
  }
  list(
    prob_identify_correct = prob_identify_correct,
    superior =
      mean(sapply(sim, function(x) x[["superior"]])),
    noninferior =
      mean(sapply(sim, function(x) x[["noninferior"]])),
    stopped_early =
      mean(sapply(sim, function(x) x[["stopped_early"]])),
    stopping_time =
      table(sapply(sim, function(x) nrow(x[["n"]]))),
    expected_ss =
      apply(do.call(rbind, lapply(sim, function(x) sum(x[["n"]]))), 2, mean),
    n_competing =
      table(do.call(rbind, lapply(sim, function(x) sum(tail(x[["is_competing"]], 1))))),
    prob_noninf_at_end = prob_noninf_at_end,
    prob_competing_at_end =
      apply(do.call(rbind, lapply(sim, function(x) tail(x[["is_competing"]], 1))), 2, mean),
    prob_find_effect =
      sum(apply(do.call(rbind, lapply(sim, function(x) x[["is_best"]][nrow(x[["is_best"]]), ])), 2, mean)),
    prob_declare_best =
      apply(do.call(rbind, lapply(sim, function(x) x[["is_best"]][nrow(x[["is_best"]]), ])), 2, mean),
    avg_estimate =
      apply(do.call(rbind, lapply(sim, function(x) tail(x[["mu"]], 1)[1, ])), 2, mean),
    avg_estimate_trans = apply(do.call(rbind, lapply(sim, function(x) plogis(tail(x[["mu"]], 1)[1, ]))), 2, mean),
    avg_allocations = apply(do.call(rbind, lapply(sim, function(x) apply(x[["n"]], 2, sum))), 2, mean),
    avg_num_noninf = avg_num_noninf
  )
}
