# Preamble ----------------------------------------------------------------
library(tidyverse)
library(rstan)
library(uuid) # for generating a new Universally Unique Identifier
library(parallel)
library(foreach)

source("scripts_specchio/R/functions_Stan_multinomial.R")

# use a session ID for similar filenames from same code run
session_ID = substr(uuid::UUIDgenerate(), 1, 8)

## Stan control settings that can be changed for testing
# (for production use 4 chains and at least 1500 iter, 250 warmup)
n_chains = 4
n_iter = 1500
n_warmup = 250
random_seed = 1543

with_households = TRUE
adults_Only = FALSE
input_dat = read_csv(here::here("data", "processed", "**** YOUR DATA HERE *****.csv"),
                     col_types = "cicccl")
# if (adults_Only) {
#   input_dat = input_dat %>%
#     filter(age >= 18)
# }

if (with_households) {
  Stan_script = "scripts_specchio/Stan/seroprev-hh-multinomial_pred_probvacc.stan"
} else {
  stop("Without households data selection not implemented yet")
}

## Stan settings, don't need to change
options(mc.cores = 8)
p_delta = 0.99
n_treedepth = 20
rstan_options(auto_write = TRUE)

# Define model --------------------------------------------------------------
model = "Sex + age_cat"

## Define desired age cuts
age_cuts = c(0, 6, 12, 18, 25, 35, 50, 65, 75, 105)
# if (adults_Only) {
#   age_cuts = c(18, 25, 35, 50, 65, 75, 105)
# }

age_ref = "[25,35)"

sex_ref = "Female"

ID_note = "_HH_multinomial"
# if (adults_Only) {
#   ID_note = paste0(ID_note, "_adultsOnly")
# }

session_ID = paste0(session_ID, ID_note)

# specchio data -------------------------------------------
input_dat = input_dat %>%
  mutate(
    age_cat = factor(cut(age, age_cuts, right = FALSE, include.lowest = TRUE)),
    Sex = fct_recode(sex, "Female" = "f", "Male" = "m"),
    UEP_S_result = fct_recode(
      S_interp,
      "0" = "neg", "1" = "pos"
    ) %>%
      as.character() %>%
      as.integer(),
    UEP_N_result = fct_recode(
      N_interp,
      "0" = "neg", "1" = "pos"
    ) %>%
      as.character() %>%
      as.integer(),
    vaccinated = as.integer(vaccinated)
  ) %>%
  group_by(hh_id) %>%
  mutate(
    hh_obs = n(),
    hh_inf_S = sum(UEP_S_result) - UEP_S_result,
    other_inf_S = hh_inf_S > 0,
    hh_inf_N = sum(UEP_N_result) - UEP_N_result,
    other_inf_N = hh_inf_N > 0,
  ) %>%
  ungroup() %>%
  droplevels()

# some fct_relevelling
input_dat = input_dat %>%
  mutate(
    Sex = fct_relevel(Sex, ref = "Female"),
    age_cat = fct_relevel(age_cat, ref = age_ref)
  ) %>%
  group_by(hh_id, hh_obs) %>%
  mutate(tot_inf_S = sum(UEP_S_result),
         tot_inf_N = sum(UEP_N_result)) %>%
  arrange(desc(hh_obs), tot_inf_S, hh_id)


# Control data ------------------------------------------------------------
# Roche S data from lab validation and Roche website
## number of positive controls from validation data
pos_control_S = c(172, 1423)
## number of negative controls
neg_control_S = c(185, 5991)
## number of true positives for cases
control_tp_S = c(159, 1406)
## number of false positives for controls (1-specificity)
control_fp_S = c(1, 1)

# Roche N data from lab validation and
# https://www.thelancet.com/action/showPdf?pii=S1473-3099%2820%2930634-4
# Table at bottom of pg 4 "Samples with complete data available taken ≥14 days post symptom onset"
## number of positive controls from validation data
pos_control_N = c(172, 561)
## number of negative controls
neg_control_N = c(185, 976)
## number of true positives for cases
control_tp_N = c(151, 543)
## number of false positives for controls (1-specificity)
control_fp_N = c(0, 2)

# Geneva population ages data ---------------------------------------------------------------
population_data = getGE_age_cats(age_cuts)

population_data = population_data %>%
  mutate(
    Sex = fct_relevel(Sex, ref = "Female"),
    age_cat = fct_relevel(age_cat, ref = age_ref)
  )

population_data = population_data %>%
  filter(age_cat %in% levels(input_dat$age_cat),
         Sex %in% levels(input_dat$Sex)) %>%
  droplevels()

# Where to save output --------------------------------------------------------------
if (!dir.exists("output/results")) {
  dir.create("output/results")
}

output_result_list_filename = paste0(
  "~/",
  session_ID,
  "_results-list_",
  format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
  ".rds"
)

# Run models --------------------------------------------------------------
if (with_households) {
  ## model the overall seroprevalence with household effect
  fun = run_multinomial_with_hh_pred
} else {
  stop("Without households model not implemented yet")
}

calc_seropos = fun(
  model_script = Stan_script,
  dat = input_dat,
  coef_eqn = model,
  vaccinated = input_dat %>% pull(vaccinated),
  pos_control_S = pos_control_S,
  neg_control_S = neg_control_S,
  control_tp_S = control_tp_S,
  control_fp_S = control_fp_S,
  pos_control_N = pos_control_N,
  neg_control_N = neg_control_N,
  control_tp_N = control_tp_N,
  control_fp_N = control_fp_N,
  n_cores = getOption("mc.cores"),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  Stan_control = list(
    adapt_delta = p_delta,
    max_treedepth = n_treedepth
  ),
  seed = random_seed,
  pop_age_cats = population_data,
  session_ID = session_ID,
  sex_ref = sex_ref,
  age_ref = age_ref
)

# Calculate seroprevalence probabilities --------------------------------------

if (with_households) {
  seroprev_probs = calc_integrated_seroprev_probs_pred(calc_seropos$stan_posterior,
                                                       coef_eqn = model,
                                                       pop_age_cats = population_data,
                                                       session_ID = session_ID
  )
} else {
  stop("Without households calculation of probs not implemented yet")
}

pop_cats_resp_probs = seroprev_probs$pop_cats_resp_probs
pop_cats_infect_probs = seroprev_probs$pop_cats_infect_probs
pop_cats_vacc_probs = seroprev_probs$pop_cats_vacc_probs
pop_cats_any_probs = seroprev_probs$pop_cats_any_probs

subset_estimates_pospos = compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S+N+`, age_cat, Sex, strata_pop, sim))
subset_estimates_negpos = compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S-N+`, age_cat, Sex, strata_pop, sim))
subset_estimates_posneg = compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S+N-`, age_cat, Sex, strata_pop, sim))
subset_estimates_negneg = compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S-N-`, age_cat, Sex, strata_pop, sim))

subset_estimates_infect = compute_weighted_estimates(pop_cats_infect_probs)
subset_estimates_vacc = compute_weighted_estimates(pop_cats_vacc_probs)
subset_estimates_any = compute_weighted_estimates(pop_cats_any_probs)

all_results_list = c(
  model_matrix = list(calc_seropos$model_matrix),
  seroprev_probs,
  subset_estimates_pospos = list(subset_estimates_pospos),
  subset_estimates_negpos = list(subset_estimates_negpos),
  subset_estimates_posneg = list(subset_estimates_posneg),
  subset_estimates_negneg = list(subset_estimates_negneg),
  subset_estimates_infect = list(subset_estimates_infect),
  subset_estimates_vacc = list(subset_estimates_vacc),
  subset_estimates_any = list(subset_estimates_any)
)
saveRDS(all_results_list, output_result_list_filename)
cat("\n-------\n", output_result_list_filename, " saved at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n\n", sep = "")

# Save final output table ----------------------------------------------------
final_table = full_join(
  make_final_table(pop_cats_any_probs, subset_estimates_any, input_dat, age_ref, sex_ref) %>%
    rename(Any_antibody_response=`Seroprevalence (95% CI)`) %>%
    select(-p),
  make_final_table(pop_cats_infect_probs, subset_estimates_infect, input_dat, age_ref, sex_ref) %>%
    rename(natural_infection=`Seroprevalence (95% CI)`) %>%
    relocate(p, .after = last_col()),
  by=c("var", "Category", "Obs", "Vaccinated", "S Test positive", "S Test negative", "N Test positive", "N Test negative")) %>%
  full_join(make_final_table(pop_cats_vacc_probs, subset_estimates_vacc, input_dat, age_ref, sex_ref) %>% select(-p),
                             by=c("var", "Category", "Obs", "Vaccinated", "S Test positive", "S Test negative", "N Test positive", "N Test negative")) %>%
              rename(estimated_vaccination=`Seroprevalence (95% CI)`) %>%
  relocate(estimated_vaccination, .after = "Vaccinated")

res_table_file = paste0("output/results/", session_ID, "_results_table_", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")
final_table %>%
  write_excel_csv(res_table_file)

cat("\n---- Done ----\n", res_table_file, " written at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n-------\n", sep = "")

rrisk_any = calc_rrisk(all_results_list$subset_estimates_any, input_dat, age_ref, sex_ref)
rrisk_any %>%
  select(Category, Obs, `Relative risk (95% CI)`, p) %>%
  write_excel_csv(here::here("output", "results", paste0(session_ID, "_rrisk_any_", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")))
rrisk_infect = calc_rrisk(all_results_list$subset_estimates_infect, input_dat, age_ref, sex_ref)
rrisk_infect %>%
  select(Category, Obs, `Relative risk (95% CI)`, p) %>%
  write_excel_csv(here::here("output", "results", paste0(session_ID, "_rrisk_infect_", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")))

# library(shinystan)
# launch_shinystan(my_stanfit)
