# To model average

suppressPackageStartupMessages({
  library(tidyverse)
  library(loo)
  library(brms)
  library(fs)
  library(matrixStats)
  library(posterior)
})


# Get the loo / bma weights

# Identify the number of samples for each model under each approach
# Set the seed & use the rng

# I guess I'll go w/ the model_fit_loo options, then.


seed = 99879823
joint_weights = map_dfr(c('growth', 'chlor'), \(x) file.path('out', x, 'joint_weights.csv') |> read_csv() |> mutate(y = x))


n_iter = 12000
joint_samples = joint_weights |> group_by(y) |>
  mutate(samples_bf = round(posterior_prob / max(posterior_prob) * n_iter)) |>
  select(y, model_name, samples_bf)

set.seed(seed)
sample_ids = joint_samples |>
  mutate(samples_bf_idx =  map(samples_bf, \(x) sample(1:n_iter, x)) )

write_rds(sample_ids, 'out/model_average_sample_ids.rds')
# sample_ids = read_rds('out/model_average_sample_ids.rds')
get_samples = \(y, model_name, samples_loo_idx, samples_bf_idx, ...) {
  # browser()
  fit = file.path('out', y, 'model_runs', model_name, 'model_fit_large.rds') |> read_rds()
  # browser()
  # I will just hope that this is in the right order; I can't see why it shouldn't be
  epreds =   posterior_epred(fit, re_formula = NA)
  fixed_sd = matrixStats::rowSds(epreds)

  fit_post = posterior::as_draws_df(fit) |> mutate(sigma_fixed = fixed_sd)
  bf_draws = fit_post |> filter(.draw %in% samples_bf_idx)
  # loo_draws = fit_post |> filter(.draw %in% samples_loo_idx)
  bf_draws
}

sampled_draws = sample_ids |> nest() |>
  mutate(bf_draws = map2( data, y, \(df, y) df |> mutate(y=y) |>  pmap(get_samples))) |>
  # mutate(loo_draws = map(draws, 'loo'), bf_draws = map(draws, 'bf')) |>
  select(y, bf_draws)

write_rds(sampled_draws, 'out/model_average_draws.rds')
sampled_draws = read_rds('out/model_average_draws.rds')

# Now, let's dummy these into fake posterior objects

pseudo_posterior = \(draw_list) {

  map(draw_list, \(l) l |> map_dfr(as_tibble) |>
    ungroup() |>
    mutate(.chain = 1, .iteration = 1:n(), .draw = .iteration) |>
    as_draws_df()
  )
}

merged_draws = sampled_draws |>
  unnest(bf_draws) |>
  reframe(bf_draws = map_dfr(bf_draws, as_tibble) ) |>
  nest(bf_draws = bf_draws) |>
  mutate(bf_draws = pseudo_posterior(bf_draws))
write_rds(merged_draws, 'out/model_average_ensemble_draws.rds')
# merged_draws = read_rds('out/model_average_ensemble_draws.rds')

# merged_draws = read_rds('out/model_average_ensemble_draws.rds')
# now I have the draws... what to do?

# I'm missing the yhats...


# Can't work directly w/ existing prop_variables; need to re-define them in terms of sigmas
# Get posterior densities and CI's

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

# Subset and manipulate specific parameters
sigma_parms = \(x) {
  x |> select(starts_with('sd_'), sigma, sigma_total, sigma_fixed) |> zero_na() |>
    mutate(sigma_total = add_sd(sigma_total, sigma_fixed)) |>
    mutate(
      sd_clone_total = add_sd(sd_clone__Intercept, `sd_clone:trans_site__Intercept`),
      sd_subclone_total = add_sd(sd_subclone__Intercept, `sd_subclone:trans_site__Intercept`),
      sd_biol_total = add_sd(sd_clone_total, sd_subclone_total),
      sd_resid = add_sd(sd_rope__Intercept, sigma) )|>
    mutate(across(-sigma_total, \(x) x^2 / sigma_total^2, .names = 'prop_{.col}'))
}
fixed_parms = \(x) x |> select(starts_with('b_')) |> zero_na()

model_summaries_sigma = merged_draws |>
  group_by(y) |>
  mutate(summaries = map(bf_draws, \(x) sigma_parms(x) |> get_summary(from_zero = TRUE))) |>
  ungroup() |>
  mutate(ci = map(summaries, 'ci'), density = map(summaries, 'density')) |>
  select(-bf_draws, -summaries)
write_rds(model_summaries_sigma, 'out/model_average_summaries.rds')


# Fixed effect summaries
fixed_effect_summaries = merged_draws |> # pivot_longer(-y, names_to = c('weight_set'), values_to = 'draw_df') |>
  group_by(y) |>
  mutate(summaries = map(bf_draws, \(x) fixed_parms(x) |> get_summary() )) |>
  ungroup() |>
  mutate(ci = map(summaries, 'ci'), density = map(summaries, 'density')) |>
  select(-bf_draws, -summaries)
write_rds(fixed_effect_summaries, 'out/fixed_effect_summaries.rds')




