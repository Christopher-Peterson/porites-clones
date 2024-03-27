# Sensitivity analysis
suppressPackageStartupMessages({
  library(tidyverse)
  library(bridgesampling)
  library(brms)
  library(rlang)
  library(glue)
})

# argv = c("fixed_nage_nloc_nseither", "chlor","28391")
if(!exists('argv')) argv = commandArgs(TRUE)

model_name = argv[1] %|% "fixed_nage_nloc_nclone"
yvar = argv[2] %|% 'growth'
seed = argv[3] |> as.integer() %|% 90923L
adapt_delta = argv[4] |> as.numeric() %|% .95
max_treedepth = argv[5] |> as.numeric() %|% 10

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

compiled_model = pick_model(model_name) |> compile_model()




## Helper Functions ####
zero_na = \(x) {
  x[is.na(x)] <- 0
  x
}
add_sd = \(x, y) sqrt(x^2 + y^2)
pivot_parms = \(w) w |> pivot_longer(everything(), names_to = 'param', values_to = 'val')
get_ci = \(.x) {
  quants = quantile(.x, probs = c(.025, .05, .25, .5, .75, .95, .975))
  quant_names = names(quants) |> str_remove('%') |> paste0('q',n=_)
  list(tibble(!!!(quants |> set_names(quant_names))))
}
tidy_density_from_zero = \(x) {
  density(x, from = 0)[c('x', 'y')] |> as_tibble() |> list()
}
tidy_density_unbound = \(x) {
  density(x)[c('x', 'y')] |> as_tibble() |> list()
}

# Get the standard deviation parameters & calculate proportions
get_summary = \(df, from_zero = FALSE) {
  den_fun = ifelse(from_zero, tidy_density_from_zero, tidy_density_unbound)
  densities = df |> summarize(across(everything(), den_fun)) |>
    pivot_parms() |> unnest(val) |> rename(density = y)
  cis = df |>
    pivot_parms() |> unnest(val) |> group_by(param) |>
    summarize(ci = get_ci(val))
  list(ci = cis, density = densities)
}



### Parameter sweet all models ####

# Fit and process model

run_model = \(hp_list, y_var = 'len_scaled', comp_model = compiled_model, iter = 2000, cores = 4, ...) {
  model_data  = pick_data(y_var)
  control_list = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)

  model_fit = fit_model(comp_model, newdata = model_data, iter = iter,
                        new_hp = hp_list,
                        cores = cores, ..., control = control_list)

  model_fit
}

get_props = \(fit) {
  posterior::as_draws_df(fit) |>
    select(starts_with(c('prop', 'sd_subclone', 'sd_clone'))) |>
    select(-starts_with('prop_rand'))
}

# Function to figure out which parameters need to be swept

id_prop_parms = \(cmp) {
  cmp$model_parms$hp |> formals() |> names() |>
    str_subset('^hp_.+_a$') |> str_remove('hp_') |> str_remove('_a$')
  # Convieniently, the hp_rands weren't given a/b parameterizations
}

# Run sensitivty analyses on all valid props
sweep_props = \(cmp) {
  props_available = id_prop_parms(cmp)
  y_name = data_key |> filter(dataset == yvar) |> pull(y) # yvar is global command arg
  ab_df = bind_rows(
    tibble(a= 1, b = 1:10),
    tibble(a= 1:10,b = 1),
  ) |> distinct()

  results = map_dfr(props_available, \(parm) {
    prop_smry = ab_df |>
      rename_with(\(x) paste('hp', parm, x, sep = '_')) |>
      transpose() |> map(\(x) {
        run_model(x, y_name, compiled_model) |> get_props() |> get_summary()
      })
    model_summary = ab_df |>
      mutate(summary =  prop_smry, param = parm, model = model_name)
    model_summary
  })
  results
}

sensitivity_output = sweep_props(compiled_model)

## Run sensitivity ####

file.path('out/sensitivity', yvar) |> dir.create(showWarnings = FALSE, recursive = TRUE)
out_file = glue('out/sensitivity/{yvar}/{model_name}.rds')

write_rds(sensitivity_output, out_file)

#
# model_smry = read_rds(out_file)
#
# model_dat = model_smry |>
#   mutate(ci = map(summary, 'ci'), density = map(summary, 'density')) |>
#   select(-draws, -summary) |>
#   mutate(beta_med = qbeta(.5, hp_subclone_site_a, hp_subclone_site_b)) |>
#   arrange(beta_med) |>
#   mutate(x_rank = 1:n())
#
# theme_set(theme_classic())
#
# tst_plt = model_dat |> select(beta_med, x_rank, density, starts_with('hp')) |> unnest(density) |>
#   filter(param == 'prop_subclone_site') |>
#   mutate(theor_dens = dbeta(x, hp_subclone_site_a, hp_subclone_site_b)) |>
#   ggplot() +
#   geom_slab(aes(x = x_rank, y = x, thickness = density),
#             fill = NA, color = "red", linewidth = 1) +
#   geom_slab(aes(x = x_rank, y = x, thickness = theor_dens),
#             fill = NA, color = "black", linetype = 1, linewidth = .5) +
#   # xlab(+
#   ylab("Portion of variance due to Ramet:Site vs. Ramet") +
#   theme(axis.ticks.y = element_blank(),
#         axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
#   scale_x_continuous("Prior on Ramet:Site vs. Ramet (black line) vs Posterior (red line)",
#                      breaks = c(3, 11, 17), labels = c('Ramet Biased', 'Even Prior', "Ramet:Site Biased")) +
#   ylim(0, 1) + coord_flip()
#
# ggsave('figures/posterior_sens_test.png', tst_plt, dpi = 300, width = 5, height = 8)
#
# # model_dat |> select(beta_med, x_rank, density, starts_with('hp')) |> unnest(density) |>
# #   filter(param == 'prop_subclone_site') |>
# #   mutate(theor_dens = dbeta(x, hp_subclone_site_a, hp_subclone_site_b)) |>
# #   ggplot() +
# #   geom_slab(aes(y = x_rank, x = x, thickness = density),
# #             fill = NA, color = "red", linewidth = .7) +
# #   geom_slab(aes(y = x_rank, x = x, thickness = theor_dens),
# #             fill = NA, color = "black", linetype = 3, linewidth = .5) +
# #   ylab("Prior on Ramet:Site vs. Ramet (dashed red line)") +
# #   xlab("Portion of variance due to Ramet:Site vs. Ramet") +
# #   theme(axis.text.x = element_blank()) +
# #   xlim(0, 1)
# #
#
# model_dat |> ggplot(aes(x = x_rank)) +
#   geom_function()


