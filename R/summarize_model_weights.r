suppressPackageStartupMessages({
  library(tidyverse)
  library(memoise)
  library(brms)
  library(rlang)
  library(glue)
  library(rlang)
})
bf = brms::bf

if(!exists('argv')) argv = commandArgs(TRUE)
job_file = argv[1] %|% 'jobs/sensitivity.job'

source('R/model_configs.r')
#### Functions ####
# Get model weights
get_weights = function(yvar = c('growth', 'chlor')) {
  # browser()
  yvar = match.arg(yvar)
  direcs = file.path('out', yvar, "model_runs") |> dir(full.names = TRUE)

  loo_files = file.path(direcs, 'loo.rds') |> keep(file.exists)


  loo_names = loo_files |> dirname() |> basename()
  loo_results = map(loo_files, read_rds) |> set_names(loo_names)
  loo_weights = loo::loo_model_weights(loo_results)

  bridge_files = file.path(direcs, 'marginal_likelihood.rds') |> keep(file.exists)
  bridge_names = bridge_files |> dirname() |> basename()
  bridge_results = map(bridge_files, read_rds) |> set_names(bridge_names)

  bf_posterior_probs = do.call(post_prob, c(bridge_results |> unname(), list(model_names = loo_names)))

  bf_tibble = tibble(posterior_prob = unname(bf_posterior_probs), model_name = names(bf_posterior_probs))

  joint_weights = loo_weights |> as_tibble(rownames = 'model_name') |>
    rename(loo_weight = x) |>
    left_join(bf_tibble) |>
    arrange(desc(loo_weight), desc(posterior_prob))
  joint_weights |>  write_csv(file.path('out', yvar, 'joint_weights.csv'))
  invisible(joint_weights)
}
# Get formula from name
form_from_name = memoise::memoise(\(nm) {
  pick_model(nm) |> _$formula |> as.character() |> _[1]
})

#### Get Model Weights ####
joint_weights = c('growth', 'chlor') |>   map_dfr(\(x) get_weights(x) |> mutate(y = x))
# If joint_weight csv files already exist: comment out above and uncomment this:
# joint_weights = map_dfr(c('growth', 'chlor'), \(x) file.path('out', x, 'joint_weights.csv') |> read_csv() |> mutate(y = x))

all_models = joint_weights |> arrange(y, desc(posterior_prob)) |>
  select(-loo_weight) |> group_by(y)|>
  mutate(run_weight = cumsum(posterior_prob)) |>
  mutate(formula = map_chr(model_name, form_from_name) |> str_remove('y ~ '))


write_csv(all_models, 'out/model_weights_with_formula.csv')

top_models = all_models |>
  filter(run_weight <= .9) |> mutate(rank = 1:n()) |>
  mutate(keep = rank * (run_weight < 0.85) ) |>
  filter(rank <= (max(keep) + 1)) |> select(-keep)

# Create a job for sensitivity analysis ####
set.seed(142312598)
commands = top_models |> ungroup() |>
  select(model_name, y) |>
  mutate(seed = rdunif(n(), 65536)) |>
  mutate(cmd = glue('r-brms R/batch_sensitivity.r {model_name} {y} {seed}')) # go w/ default args for others

write_lines(commands$cmd, job_file)

# Format top models for table ####
formatted_top_models = top_models |>
  arrange(desc(y), run_weight) |>
  # select(-model_name) |>
  mutate(across(c(posterior_prob, run_weight), \(x) sprintf('%.3f', x))) |>
  select(Response = y, Model = formula, `Posterior Probability` = posterior_prob,
         `Cummulative Probability` = run_weight, model_name) |>
  mutate(Model = Model |>
           str_replace_all('subclone', 'Ramet') |>
           str_replace_all('clone', 'Genet') |>
           str_replace_all('trans_site', 'Site') |>
           str_replace('rope', "Rope") |>
           str_replace_all('loc', 'Location') |>
           str_replace_all('age', 'Age'),
         Response = recode(Response, growth = "Growth", chlor = "Chlorophyll"))

formatted_top_models |> select(-model_name) |> write_csv('top_models.csv')
formatted_top_models |> write_csv('formatted_top_models.csv')
