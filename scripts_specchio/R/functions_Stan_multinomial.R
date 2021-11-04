#' @title Computes seroprevalences
#'
#' @description for given variable and reference
#'
#' @param df serprevalance draws dataframe
#' @param var variable of interest
#' @param ref reference value
#'
#' @return return
computeSeroPrev = function(df, var, ref) {
  colnames(df)[colnames(df) == var] = "val"
  df$var = var
  df %>%
    group_by(sim, val, var) %>%
    summarize(p = weighted.mean(seropos, strata_pop)) %>%
    group_by(sim) %>%
    mutate(coef_val = ifelse(val == ref, NA,
                             ifelse(p > p[val == ref], 1, -1)
    ))
}

#' @title Get Geneva age categories
#'
#' @description Reads in GE population data by age and computes categories for given age cuts
#'
#' @param age_cuts Vector of ages to cut the population at
#'
#' @return tibble of Sex, ageCat, population for that strata, and % of the total
getGE_age_cats = function(age_cuts) {
  # GVA population from Office Cantonal STAT
  GE_pop_all = read_csv("data/processed/tidied_GEpopByAgeSexOn2020-12-31.csv",
                        comment = "#",
                        col_types = "cfd"
  )
  GE_pop_all[GE_pop_all$Age == ">= 105", "Age"] = "105"
  GE_pop_all = GE_pop_all %>% mutate(
    Age = as.integer(Age),
    age_cat = cut(Age, age_cuts, right = FALSE, include.lowest = TRUE)
  )

  ## find population-level proportions of age and sex for sero_dat
  GE_pop_cats = GE_pop_all %>%
    drop_na() %>%
    count(Sex, age_cat, wt = Total, name = "strata_pop") %>%
    mutate(percent = strata_pop / sum(strata_pop))
  return(GE_pop_cats)
}

#' @title Run basic Stan analysis
#' @description Runs the basic Stan seroprevalence model not accounting for household random effects
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param coef_eqn formula in character format of variables to model
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param n_cores number of cores to use for parallel computation
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#' @param seed Random seed
#' @param chains number of chains to use in Stan
#' @param iter number of total iterations per chain in Stan
#' @param warmup number of warmup iterations per chain in Stan
#' @param Stan_control List of control parameters to be passed on to Stan
#' @param session_ID Unique identifier for linking output files
#'
#' @return a list with parameter posteriors and results
run_analysis_stan_no_hh = function(model_script,
                                   dat,
                                   coef_eqn,
                                   pos_control,
                                   neg_control,
                                   control_tp,
                                   control_fp,
                                   n_cores = detectCores() - 2,
                                   sex_ref = "Female",
                                   age_ref = "[20,50)",
                                   seed,
                                   chains,
                                   iter,
                                   warmup,
                                   Stan_control,
                                   session_ID,
                                   ...) {
  stopifnot(nrow(dat) == sum(dat$UEP_result == 1) + sum(dat$UEP_result == 0))
  # Set model matrix
  X = model.matrix(as.formula(paste("~", coef_eqn)), data = dat)

  # name of the stan output file
  stan_sampling_filename = paste0(
    "output/results/", session_ID, "-stan_fit_",
    basename(model_script), "_", chains,
    "_", iter, "_", warmup, "_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".rds"
  )

  # Run stan model
  cat(paste0("Starting sampling session ", session_ID, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  stan_posterior = rstan::sampling(stan_model(model_script),
                                   data = list(
                                     N_survey = nrow(dat),
                                     p_vars = ncol(X),
                                     X = X,
                                     survey_pos = dat$UEP_result,
                                     J_studies = length(pos_control),
                                     N_pos_control = pos_control,
                                     control_tp = control_tp,
                                     N_neg_control = neg_control,
                                     control_fp = control_fp
                                   ),
                                   chains = chains,
                                   iter = iter,
                                   warmup = warmup,
                                   seed = seed,
                                   control = Stan_control
  )
  cat(paste0("Finished sampling at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\nSaving as RDS file ", stan_sampling_filename, "\n"))
  saveRDS(stan_posterior, stan_sampling_filename)
  cat(paste0("Finished saving RDS ", stan_sampling_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  return(list(stan_posterior = stan_posterior, model_matrix = X))
}

#' @title Run Stan analysis with household random effect
#' @description Runs the Stan seroprevalence model accounting for household random effects
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param coef_eqn formula in character format of variables to model
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param n_cores number of cores to use for parallel computation
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#' @param seed Random seed
#' @param chains number of chains to use in Stan
#' @param iter number of total iterations per chain in Stan
#' @param warmup number of warmup iterations per chain in Stan
#' @param Stan_control List of control parameters to be passed on to Stan
#' @param session_ID Unique identifier for linking output files
#'
#' @return a list with parameter posteriors and results
run_multinomial_with_hh = function(model_script,
                                   dat,
                                   coef_eqn,
                                   vaccinated,
                                   pos_control_S,
                                   neg_control_S,
                                   control_tp_S,
                                   control_fp_S,
                                   pos_control_N,
                                   neg_control_N,
                                   control_tp_N,
                                   control_fp_N,
                                   n_cores = detectCores() - 2,
                                   sex_ref = "Female",
                                   age_ref = "[20,50)",
                                   seed,
                                   chains,
                                   iter,
                                   warmup,
                                   Stan_control,
                                   session_ID,
                                   ...) {
  stopifnot((nrow(dat) == sum(dat$UEP_S_result == 1) + sum(dat$UEP_S_result == 0)) | (nrow(dat) == sum(dat$UEP_N_result == 1) + sum(dat$UEP_N_result == 0)))
  # Set model matrix
  X = model.matrix(as.formula(paste("~", coef_eqn)), data = dat)

  # Unique household ids from 1 to N
  u_household_ids = unique(dat$hh_id)
  dat$u_household_id = map_dbl(dat$hh_id, ~ which(u_household_ids == .))

  # name of the stan output file
  stan_sampling_filename = paste0(
    "output/results/", session_ID, "-stan_fit_",
    basename(model_script), "_", chains,
    "_", iter, "_", warmup, "_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".rds"
  )

  # Run stan model
  cat(paste0("Starting sampling session ", session_ID, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  stan_posterior = rstan::sampling(stan_model(model_script),
                                   data = list(
                                     N_survey = nrow(dat),
                                     H = length(u_household_ids),
                                     survey_pos = matrix(c(dat$UEP_S_result, dat$UEP_N_result), ncol=2),
                                     hh = dat$u_household_id,
                                     p_vars = ncol(X),
                                     X = X,
                                     vaccinated = vaccinated,
                                     p_vacc_response = .999, # don't set to 1 to avoid errors
                                     p_infect_prior = c(10, 1 ,1),
                                     J_studies_S = length(pos_control_S),
                                     J_studies_N = length(pos_control_N),
                                     N_pos_control_S = pos_control_S,
                                     control_tp_S = control_tp_S,
                                     N_neg_control_S = neg_control_S,
                                     control_fp_S = control_fp_S,
                                     N_pos_control_N = pos_control_N,
                                     control_tp_N = control_tp_N,
                                     N_neg_control_N = neg_control_N,
                                     control_fp_N = control_fp_N
                                   ),
                                   chains = chains,
                                   iter = iter,
                                   warmup = warmup,
                                   seed = seed,
                                   control = Stan_control
  )
  cat(paste0("Finished sampling at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\nSaving as RDS file ", stan_sampling_filename, "\n"))
  saveRDS(stan_posterior, stan_sampling_filename)
  cat(paste0("Finished saving RDS ", stan_sampling_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  return(list(stan_posterior = stan_posterior, model_matrix = X))
}


#' @title Run Stan analysis with household random effect and prob prediction
#' @description Runs the Stan seroprevalence model accounting for household random effects
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param coef_eqn formula in character format of variables to model
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param n_cores number of cores to use for parallel computation
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#' @param seed Random seed
#' @param chains number of chains to use in Stan
#' @param iter number of total iterations per chain in Stan
#' @param warmup number of warmup iterations per chain in Stan
#' @param Stan_control List of control parameters to be passed on to Stan
#' @param session_ID Unique identifier for linking output files
#'
#' @return a list with parameter posteriors and results
run_multinomial_with_hh_pred <-  function(model_script,
                                          dat,
                                          coef_eqn,
                                          vaccinated,
                                          pos_control_S,
                                          neg_control_S,
                                          control_tp_S,
                                          control_fp_S,
                                          pos_control_N,
                                          neg_control_N,
                                          control_tp_N,
                                          control_fp_N,
                                          n_cores = detectCores() - 2,
                                          sex_ref = "Female",
                                          age_ref = "[20,50)",
                                          seed,
                                          chains,
                                          iter,
                                          warmup,
                                          Stan_control,
                                          session_ID,
                                          pop_age_cats,
                                          ...) {
  stopifnot((nrow(dat) == sum(dat$UEP_S_result == 1) + sum(dat$UEP_S_result == 0)) | (nrow(dat) == sum(dat$UEP_N_result == 1) + sum(dat$UEP_N_result == 0)))
  # Set model matrix
  X = model.matrix(as.formula(paste("~", coef_eqn)), data = dat)

  # Get all unique combinations
  X_comb = pop_age_cats %>%
    select(age_cat, Sex) %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  # Unique household ids from 1 to N
  u_household_ids = unique(dat$hh_id)
  dat$u_household_id = map_dbl(dat$hh_id, ~ which(u_household_ids == .))

  # name of the stan output file
  stan_sampling_filename = paste0(
    "~/", session_ID, "-stan_fit_",
    basename(model_script), "_", chains,
    "_", iter, "_", warmup, "_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".rds"
  )

  # Run stan model
  cat(paste0("Starting sampling session ", session_ID, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  stan_model_pred <- stan_model(model_script, auto_write = T)

  H <- length(u_household_ids)
  logit <- function(x) log(x/(1-x))
  param_initializer <- function() {
    list(
      logit_sens_S = logit(runif(2, 0.8, 1)),
      logit_sens_N = logit(runif(2, 0.8, 1)),
      logit_spec_S = logit(runif(2, 0.8, 1)),
      logit_spec_N = logit(runif(2, 0.8, 1)),
      sigma_h = abs(rnorm(1, 1, .05)),
      sigma_vacc = abs(rnorm(1, 1, .05)),
      eta_h = rnorm(H, 0, .05),
      eta_h_vacc = rnorm(H, 0, .05),
      beta = rnorm(10, 0, .1),
      beta_vacc = rnorm(10, 0, .1),
      sigma_logit_spec_S = abs(rnorm(1, 1, .05)),
      sigma_logit_spec_N = abs(rnorm(1, 1, .05)),
      sigma_logit_sens_S = abs(rnorm(1, 1, .05)),
      sigma_logit_sens_N = abs(rnorm(1, 1, .05)),
      mu_logit_spec_S = logit(runif(1, .8, 1)),
      mu_logit_spec_N = logit(runif(1, .8, 1)),
      mu_logit_sens_S = logit(runif(1, .8, 1)),
      mu_logit_sens_N = logit(runif(1, .8, 1)),
      p_vacc_response = runif(1, .8, 1),
      p_infect_response = c(runif(1, .9, 1)) %>% c(., (1-.)/2, (1-.)/2)
    )
  }


  param_inits <- map(1:chains,~ param_initializer())

  stan_posterior = rstan::sampling(stan_model_pred,
                                   data = list(
                                     N_survey = nrow(dat),
                                     H = H,
                                     survey_pos = matrix(c(dat$UEP_S_result, dat$UEP_N_result), ncol=2),
                                     hh = dat$u_household_id,
                                     p_vars = ncol(X),
                                     X = X,
                                     X_comb = X_comb,
                                     N_comb = nrow(X_comb),
                                     vaccinated = vaccinated,
                                     p_infect_prior = c(10, 1 ,1),
                                     J_studies_S = length(pos_control_S),
                                     J_studies_N = length(pos_control_N),
                                     N_pos_control_S = pos_control_S,
                                     control_tp_S = control_tp_S,
                                     N_neg_control_S = neg_control_S,
                                     control_fp_S = control_fp_S,
                                     N_pos_control_N = pos_control_N,
                                     control_tp_N = control_tp_N,
                                     N_neg_control_N = neg_control_N,
                                     control_fp_N = control_fp_N
                                   ),
                                   chains = chains,
                                   iter = iter,
                                   warmup = warmup,
                                   seed = seed,
                                   control = Stan_control,
                                   init = param_inits,
                                   pars = c("p", "p_vacc", "p_vacc_inf", "ll_vacc", "P1", "P2", "P3", "Pv", "theta", "eta_h", "eta_h_vacc", "P_pospos", "P_negpos", "P_posneg", "P_negneg"),
                                   include=FALSE
  )

  cat(paste0("Finished sampling at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\nSaving as RDS file ", stan_sampling_filename, "\n"))
  saveRDS(stan_posterior, stan_sampling_filename)
  cat(paste0("Finished saving RDS ", stan_sampling_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  return(list(stan_posterior = stan_posterior, model_matrix = X))
}

#' @title Calculate seroprevalence probabilities for all strata
#' @description Calculates seroprevalence probabilities by taking the inv logit of the parameters %*% strata
#'
#' @param stan_fit Output from a Stan sampling function i.e. the posterior
#' @param pop_age_cats age categories
#' @param coef_eqn formula in character format of variables to model
#' @param n_cores number of cores to use for parallel computation
#' @param session_ID Unique identifier for linking output files
#'
#' @return List of population strata and their probabilities, extracted Stan parameters and the model matrix used in the probability calculation
calc_seroprev_probs = function(stan_fit,
                               pop_age_cats,
                               coef_eqn,
                               n_cores,
                               session_ID) {

  ## extract parameters
  beta = rstan::extract(stan_fit, pars = "beta")[[1]]

  pop_cat_mat = pop_age_cats %>%
    select(age_cat, Sex) %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  stopifnot(ncol(beta) == ncol(pop_cat_mat)) # otherwise the matrix multiplication ain't gonna work...

  # name of the seroprevalence probabilities output file
  seroprev_probs_filename = paste0(
    "output/results/",
    session_ID,
    "_integrated-probs_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  )

  cat(paste0("Starting to calculate seroprevalence probabilities at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  ## setup parallel code
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  pop_cats_probs = foreach(
    i = 1:nrow(pop_cat_mat),
    .combine = rbind,
    .inorder = F,
    .packages = c("tidyverse", "foreach"), .errorhandling = "pass"
  ) %dopar% {
    foreach(
      j = 1:nrow(beta),
      .combine = rbind,
      .inorder = T
    ) %do% {
      # Compute probability
      prob = plogis(beta[j, , drop = F] %*% t(pop_cat_mat[i, , drop = F]))

      tibble_row(
        age_cat = pop_age_cats$age_cat[i],
        Sex = pop_age_cats$Sex[i],
        strata_pop = pop_age_cats$strata_pop[i],
        seropos = prob,
        sim = j
      )
    }
  }
  parallel::stopCluster(cl)
  cat(paste0("Finished calculating seroprevalence probabilities at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
  saveRDS(pop_cats_probs, seroprev_probs_filename)

  return(list(pop_cats_probs = pop_cats_probs, params = rstan::extract(stan_fit), pop_cat_mat = pop_cat_mat))
}

#' @title Calculate integrated seroprevalence probabilities for all strata
#' @description Calculates seroprevalence probabilities by taking the integral of the inv logit of the quantile normal distribution with mean parameters %*% strata and household sigma
#'
#' @param stan_fit Output from a Stan sampling function i.e. the posterior
#' @param pop_age_cats age categories
#' @param coef_eqn formula in character format of variables to model
#' @param n_cores number of cores to use for parallel computation
#' @param session_ID Unique identifier for linking output files
#'
#' @return List of population strata and their probabilities, extracted Stan parameters and the model matrix used in the probability calculation
calc_integrated_seroprev_probs = function(stan_fit,
                                          pop_age_cats,
                                          coef_eqn,
                                          n_cores,
                                          session_ID) {

  ## extract parameters
  beta = rstan::extract(stan_fit, pars = "beta")[[1]]
  sigma = rstan::extract(stan_fit, pars = "sigma_h")[[1]]

  pop_cat_mat = pop_age_cats %>%
    select(age_cat, Sex) %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  stopifnot(ncol(beta) == ncol(pop_cat_mat)) # otherwise the matrix multiplication ain't gonna work...

  # name of the seroprevalence probabilities output file
  seroprev_probs_filename = paste0(
    "output/results/",
    session_ID,
    "_integrated-probs_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  )

  cat(paste0("Starting to calculate seroprevalence probabilities at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  ## setup parallel code
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  pop_cats_probs = foreach(
    i = 1:nrow(pop_cat_mat),
    .combine = rbind,
    .inorder = F,
    .packages = c("tidyverse", "foreach"), .errorhandling = "pass"
  ) %dopar% {
    foreach(
      j = 1:nrow(beta),
      .combine = rbind,
      .inorder = T
    ) %do% {
      # Compute probability integrating across household random effects
      prob = integrate(function(x) {
        plogis(qnorm(
          x, beta[j, , drop = F] %*% t(pop_cat_mat[i, , drop = F]),
          sigma[j]
        ))
      }, 0, 1, rel.tol = 1e-8, stop.on.error = FALSE)

      tibble_row(
        age_cat = pop_age_cats$age_cat[i],
        Sex = pop_age_cats$Sex[i],
        strata_pop = pop_age_cats$strata_pop[i],
        seropos = prob[[1]],
        sim = j,
        message = prob[[4]]
      )
    }
  }
  parallel::stopCluster(cl)
  cat(paste0("Finished calculating seroprevalence probabilities at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
  saveRDS(pop_cats_probs, seroprev_probs_filename)

  failures = pop_cats_probs[which(pop_cats_probs$message != "OK"), ]
  if (nrow(failures) > 0L) {
    cat("*********WARNING**********\n There were ", nrow(failures), " integration failures\n")
    print(failures)
  }

  return(list(pop_cats_probs = pop_cats_probs, params = rstan::extract(stan_fit), pop_cat_mat = pop_cat_mat, failures = failures))
}


#' @title Calculate integrated seroprevalence probabilities for all strata
#' @description Calculates seroprevalence probabilities by taking the integral of the inv logit of the quantile normal distribution with mean parameters %*% strata and household sigma
#'
#' @param stan_fit Output from a Stan sampling function i.e. the posterior
#' @param pop_age_cats age categories
#' @param coef_eqn formula in character format of variables to model
#' @param n_cores number of cores to use for parallel computation
#' @param session_ID Unique identifier for linking output files
#'
#' @return List of population strata and their probabilities, extracted Stan parameters and the model matrix used in the probability calculation
calc_integrated_seroprev_probs_pred = function(stan_fit,
                                               pop_age_cats,
                                               coef_eqn,
                                               n_cores,
                                               session_ID) {

  ## extract probabilities
  post_probs <- rstan::extract(stan_fit, pars = "post_probs")[[1]]
  post_prob_infect <- rstan::extract(stan_fit, pars = "p_est")[[1]]
  post_prob_vacc = rstan::extract(stan_fit, pars = "p_est_vacc")[[1]]
  post_prob_any <- rstan::extract(stan_fit, pars = "p_any")[[1]]

  pop_cat_mat = pop_age_cats %>%
    select(age_cat, Sex) %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  stopifnot(ncol(beta) == ncol(pop_cat_mat)) # otherwise the matrix multiplication ain't gonna work...

  # name of the seroprevalence probabilities output file
  seroprev_probs_filename = paste0(
    "~/",
    session_ID,
    "_integrated-probs_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  )

  pop_cats_resp_probs <- map_df(1:nrow(pop_cat_mat), function(i) {
    as_tibble(post_probs[, i, ]) %>%
      magrittr::set_colnames(c("S+N+", "S-N+", "S+N-", "S-N-")) %>%
      mutate(
        age_cat = pop_age_cats$age_cat[i],
        Sex = pop_age_cats$Sex[i],
        strata_pop = pop_age_cats$strata_pop[i],
        sim = row_number()
      )
  })

  pop_cats_infect_probs <- map_df(1:nrow(pop_cat_mat), function(i) {
      tibble(
        age_cat = pop_age_cats$age_cat[i],
        Sex = pop_age_cats$Sex[i],
        strata_pop = pop_age_cats$strata_pop[i],
        seropos = post_prob_infect[, i]
      ) %>%
      mutate(sim = row_number())
  })

  pop_cats_vacc_probs = map_df(1:nrow(pop_cat_mat), function(i) {
    tibble(
      age_cat = pop_age_cats$age_cat[i],
      Sex = pop_age_cats$Sex[i],
      strata_pop = pop_age_cats$strata_pop[i],
      seropos = post_prob_vacc[, i]
    ) %>%
      mutate(sim = row_number())
  })

  pop_cats_any_probs <- map_df(1:nrow(pop_cat_mat), function(i) {
    tibble(
      age_cat = pop_age_cats$age_cat[i],
      Sex = pop_age_cats$Sex[i],
      strata_pop = pop_age_cats$strata_pop[i],
      seropos = post_prob_any[, i]
    ) %>%
      mutate(sim = row_number())
  })

  all_sero_probs = list(pop_cats_resp_probs = pop_cats_resp_probs,
                        pop_cats_infect_probs = pop_cats_infect_probs,
                        pop_cats_vacc_probs = pop_cats_vacc_probs,
                        pop_cats_any_probs = pop_cats_any_probs,
                        # params = rstan::extract(stan_fit),
                        pop_cat_mat = pop_cat_mat)
  saveRDS(all_sero_probs, seroprev_probs_filename)

  return(all_sero_probs)
}

#' @title Compute population strata weighted probabilities
#' @description Calculates seroprevalence probabilities by weighting by strata population size
#'
#' @param pop_cats_probs Population strata and their probabilities
#'
#' @return Tibble of overall estimates, age specific estimates and sex specific estimates
compute_weighted_estimates = function(pop_cats_probs) {

  ## overall estimate
  overall_re = pop_cats_probs %>%
    mutate(
      var = "Overall",
      val = ""
    ) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, strata_pop)) %>%
    ungroup()

  ## find age specific probabilities in order to make relative risks
  age_re = pop_cats_probs %>%
    filter(Sex == sex_ref) %>%
    mutate(var = "Age") %>%
    rename(val = age_cat) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, strata_pop)) %>%
    ungroup()

  # sex-specific probabilities
  sex_re = pop_cats_probs %>%
    filter(age_cat == age_ref) %>%
    mutate(var = "Sex") %>%
    rename(val = Sex) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, strata_pop)) %>%
    ungroup()

  return(bind_rows(
    overall_re,
    sex_re,
    age_re
  ))
}

#' @title Make final table
#' @description Calculates seropositivity and seroprevalence for age, sex and overall strata
#'
#' @param pop_cats_probs Population strata and their probabilities
#' @param subset_estimates Tibble of overall estimates, age specific estimates and sex specific estimates
#' @param dat Input sero data
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#'
#' @return Formatted data frame of age, sex, overall, N, N pos and neg for S and N tests, % vaccinated, seroprevalence +/- credible interval, Bayesian p-value
make_final_table = function(pop_cats_probs, subset_estimates, input_dat, age_ref, sex_ref) {

  ## Compute age estimates
  age_seroprev = computeSeroPrev(pop_cats_probs, "age_cat", age_ref) %>%
    bind_rows(
      subset_estimates %>%
        filter(var == "Overall") %>%
        mutate(var = "age_cat", val = "all") %>%
        mutate(coef_val = as.numeric(NA))
    )

  age_counts = input_dat %>%
    rbind(input_dat %>% mutate(age_cat = "all")) %>%
    group_by(age_cat) %>%
    summarize(
      n = n(),
      pos_S = sum(UEP_S_result),
      neg_S = n() - pos_S,
      pos_N = sum(UEP_N_result),
      neg_N = n() - pos_N,
      vacc = sum(vaccinated)
    ) %>%
    mutate(age_cat = as.character(age_cat))

  age_res = left_join(age_seroprev, age_counts, by = c("val" = "age_cat")) %>%
    mutate(var = "Age")

  ## Compute sex estimates
  sex_seroprev = computeSeroPrev(pop_cats_probs, "Sex", sex_ref)

  sex_counts = input_dat %>%
    mutate(val = Sex) %>%
    group_by(val) %>%
    summarize(
      n = n(),
      pos_S = sum(UEP_S_result),
      neg_S = n() - pos_S,
      pos_N = sum(UEP_N_result),
      neg_N = n() - pos_N,
      vacc = sum(vaccinated)
    )

  sex_res = left_join(sex_seroprev, sex_counts)

  # Combine results for seropositive estimate
  seropos = bind_rows(age_res, sex_res) %>%
    group_by(var, val, n, pos_S, neg_S, pos_N, neg_N, vacc) %>%
    summarize(
      p_mean = mean(p),
      p_quantile.025 = quantile(p, probs = .025),
      p_quantile.975 = quantile(p, probs = .975),
      p = ifelse(is.na(mean(coef_val)), "--",
                 min(2 * c(mean(coef_val > 0), mean(coef_val < 0))) %>%
                   formatC(4, format = "f")
      )
    ) %>%
    ungroup() %>%
    mutate(
      `Seroprevalence (95% CI)` = paste0(
        (100 * p_mean) %>% formatC(1, format = "f"), " (",
        (100 * p_quantile.025) %>% formatC(1, format = "f"), "-",
        (100 * p_quantile.975) %>% formatC(1, format = "f"), ")"
      ),
      pos_S = paste0(pos_S, " (", formatC(100 * pos_S / n, 1, format = "f"), "%)"),
      neg_S = paste0(neg_S, " (", formatC(100 * neg_S / n, 1, format = "f"), "%)"),
      pos_N = paste0(pos_N, " (", formatC(100 * pos_N / n, 1, format = "f"), "%)"),
      neg_N = paste0(neg_N, " (", formatC(100 * neg_N / n, 1, format = "f"), "%)"),
      vacc = paste0(vacc, " (", formatC(100 * vacc / n, 1, format = "f"), "%)")
    ) %>%
    rename(
      Category = val, Obs = n, Vaccinated = vacc, `S Test positive` = pos_S, `S Test negative` = neg_S, `N Test positive` = pos_N, `N Test negative` = neg_N
    ) %>%
    mutate(Category = factor(Category, levels = c(str_sort(unique(levels(input_dat$age_cat)), numeric = TRUE), "Male", "Female", "all"))) %>%
    relocate(Vaccinated, .after="Obs") %>%
    arrange(Category) %>%
    select(-starts_with("p_"))

  # Formatted final table
  final_table = seropos %>% arrange(Category)

  return(final_table)
}

#' @title Calculate relative risks
#' @description Calculates seropositivity and relative risks for age and sex strata
#'
#' @param subset_estimates Tibble of overall estimates, age specific estimates and sex specific estimates
#' @param dat Input sero data
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#'
#' @return Formatted data frame of age, sex, N, N pos and neg for S and N tests, % vaccinated, relative risk +/- CI, Bayesian p-value
calc_rrisk = function(subset_estimates, input_dat, age_ref, sex_ref) {
  # Compute relative risk
  rrisk = subset_estimates %>%
    filter(var == "Age") %>%
    group_by(sim) %>%
    mutate(rr = ifelse(val == age_ref, NA, p / p[val == age_ref])) %>%
    ungroup() %>%
    left_join(
      input_dat %>%
        group_by(age_cat) %>%
        summarize(
          n = n(),
          pos_S = sum(UEP_S_result),
          neg_S = n() - pos_S,
          pos_N = sum(UEP_N_result),
          neg_N = n() - pos_N,
          vacc = sum(vaccinated)
        ) %>% mutate(age_cat = as.character(age_cat)),
      by = c("val" = "age_cat")
    ) %>%
    bind_rows(
      subset_estimates %>%
        filter(var == "Sex") %>%
        group_by(sim) %>%
        mutate(rr = ifelse(val == sex_ref, NA, p / p[val == sex_ref])) %>%
        ungroup() %>%
        mutate(val = ifelse(val == sex_ref, "Female", "Male")) %>%
        left_join(input_dat %>%
                    mutate(val = ifelse(Sex == sex_ref, "Female", "Male")) %>%
                    group_by(val) %>%
                    summarize(
                      n = n(),
                      pos_S = sum(UEP_S_result),
                      neg_S = n() - pos_S,
                      pos_N = sum(UEP_N_result),
                      neg_N = n() - pos_N,
                      vacc = sum(vaccinated)
                    ))
    ) %>%
    group_by(var, val, n, pos_S, neg_S, pos_N, neg_N, vacc) %>%
    summarize(
      `Relative risk (95% CI)` = ifelse(is.na(mean(rr)), "--",
                                        paste0(
                                          mean(rr, na.rm = T) %>%
                                            formatC(2, format = "f"),
                                          " (", quantile(rr, probs = .025, na.rm = T) %>%
                                            formatC(2, format = "f"), "-",
                                          quantile(rr, probs = .975, na.rm = T) %>%
                                            formatC(2, format = "f"), ")"
                                        )
      ),
      p = ifelse(is.na(mean(rr)), "--",
                 min(2 * c(
                   mean(rr > 1, na.rm = T),
                   mean(rr < 1, na.rm = T)
                 )) %>%
                   formatC(4, format = "f")
      )
    ) %>%
    ungroup() %>%
    mutate(
      pos_S = paste0(pos_S, " (", formatC(100 * pos_S / n, 1, format = "f"), "%)"),
      neg_S = paste0(neg_S, " (", formatC(100 * neg_S / n, 1, format = "f"), "%)"),
      pos_N = paste0(pos_N, " (", formatC(100 * pos_N / n, 1, format = "f"), "%)"),
      neg_N = paste0(neg_N, " (", formatC(100 * neg_N / n, 1, format = "f"), "%)"),
      vacc = paste0(vacc, " (", formatC(100 * vacc / n, 1, format = "f"), "%)")
    ) %>%
    rename(
      Category = val, Obs = n, Vaccinated = vacc, `S Test positive` = pos_S, `S Test negative` = neg_S, `N Test positive` = pos_N, `N Test negative` = neg_N
    ) %>%
    select(-var) %>%
    mutate(Category = factor(Category, levels = c(str_sort(unique(levels(input_dat$age_cat)), numeric = TRUE), "Male", "Female"))) %>%
    arrange(Category)

  return(rrisk)
}
