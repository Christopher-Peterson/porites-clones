# Run the mixed models w/ brms
suppressPackageStartupMessages({
  library(tidyverse)
  library(bridgesampling)
  library(brms)
  library(rlang)
  library(glue)
  library(loo)
})


if(!exists('argv')) argv = commandArgs(TRUE)

model_name = argv[1] %|% "fixed_full_allbiol"
seed = argv[2] |> as.integer() %|% 90923L
adapt_delta = argv[3] |> as.numeric() %|% .95
max_treedepth = argv[4] |> as.numeric() %|% 10

set.seed(seed)

# Datasets
datasets = list(
  growth = read_csv('data/2022_growth_nouniques_age_NEW.csv') |>
    mutate(len_scaled = len_ct |> scale() |> as.numeric()),
  chlor = 'data/2022_chl_nouniques_age_NEW.csv' |> read_csv() |>
    mutate(R_fin_scaled = R_fin |> scale() |> as.numeric())
    ) |> map(\(x) x |> mutate(loc = factor(loc, levels = c('home', 'away'))) |>
           select(-age) |> rename(age = age2))
data_key = tribble(~'dataset', ~'y',
                   'growth', 'len_scaled',
                   'chlor', 'R_fin_scaled') |>
  # Create directory structure for output
  mutate(direc = file.path('out', dataset, 'model_runs', model_name))

# Create the directories
data_key$direc |> walk(dir.create, recursive = TRUE, showWarnings = FALSE)

# Some functions for automated stan model exploration ####

### Functions for setting up the model
# Select model terms

source("R/model_configs.r")

# Compile the model
compile_model = \(model_parms, model_data = pick_data('len_scaled'), ...) {
  # browser()
  # ... extra args to brm
  fit = brm(model_parms$formula, model_data, prior = model_parms$prior,
            stanvar = model_parms$base_stanvar + model_parms$hp(),
            save_pars = save_pars(all = TRUE),
            # backend = 'cmdstanr', #file = model_parms$cmp_file,
            iter = 0, chains = 0, ...)
  list(cmp = fit, model_parms = model_parms, misc = list(...), hp = model_parms$hp)
}
# Fit the model; use use cmp_model$hp() to update the model parameters
# Use argument newdata to change the y
fit_model = \(cmp_model, new_hp, ...) { # Pass argument newdata to change y's
  # browser()
  update_hp = !is_missing(new_hp)
  if(update_hp) {
    hp = do.call(cmp_model$model_parms$hp, new_hp)
  } else hp = cmp_model$model_parms$hp()
  # browser()
  brms:::update.brmsfit(cmp_model$cmp,
                        stanvars = cmp_model$model_parms$base_stanvar + hp,
                        ..., save_pars = save_pars(all = TRUE))
}

# Helper function for picking names
name_subset = \(x, pat) x[names(x) |> str_subset(pat)]



# Fit and process model
run_model = \(y_var, comp_model, loo_iter = 2000, bridge_iter = 6000, cores = 4, ...) {
  model_data  = pick_data(y_var)

  control_list = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)

  model_fit_loo = fit_model(comp_model, newdata = model_data, iter = loo_iter, cores = cores, ..., control = control_list)
  loo_results = model_fit_loo |> loo(moment_match = TRUE, reloo = TRUE, cores = cores, ...)

  model_fit_bridge = fit_model(comp_model, newdata = model_data, iter = bridge_iter, cores = cores, ..., control = control_list)
  bridge_results = bridge_sampler(model_fit_bridge)

  out_dir = data_key |> filter(y == y_var) |> pull(direc)

  write_rds(model_fit_loo, file.path(out_dir, 'model_fit.rds'))
  write_rds(model_fit_bridge, file.path(out_dir, 'model_fit_large.rds'))
  write_rds(loo_results, file.path(out_dir, 'loo.rds'))
  write_rds(bridge_results, file.path(out_dir, 'marginal_likelihood.rds'))
  # model_summary = summary(model_fit_loo)
  # write_rds(model_summary, file.path(out_dir, 'model_smry.rds'))
  invisible()
}


### Compile the model ####

compiled_model = pick_model(model_name) |> compile_model()

### Run all models

data_key$y |> walk(run_model, comp_model = compiled_model)


