# Plot model averages
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(loo)
  library(brms)
  library(glue)
  library(fs)
  library(ggdist)
  library(matrixStats)
  library(posterior)
})

# Utility functions ####
rename_variances = \(x) recode(x,
                               "sd_clone__Intercept" = 'Genet (w/o site)|sd',
                               "sd_clone:trans_site__Intercept" = 'Genet : Site|sd',
                               "sd_rope__Intercept" = 'Rope|sd',
                               "sd_subclone__Intercept" = 'Ramet (w/o site)|sd',
                               "sd_subclone:trans_site__Intercept" = 'Ramet:Site|sd',
                               "sigma" = 'Residual|sd',
                               "sigma_total" = 'Total S.D.|sd',
                               "sigma_fixed" = 'Fixed Effects (mostly site)|sd',
                               "sd_clone_total" = 'Genet (Total)|sd',
                               "sd_biol_total" = 'Genet + Ramet (Total)|sd',
                               "sd_resid" = 'Residual + Rope|sd',
                               "sd_subclone_total" = 'Ramet (Total)|sd',
                               "prop_sd_clone__Intercept" = 'Genet (w/o site)|prop',
                               "prop_sd_clone:trans_site__Intercept" =  'Genet:Site|prop',
                               "prop_sd_rope__Intercept" =  'Rope|prop',
                               "prop_sd_subclone__Intercept" =  'Ramet (w/o site)|prop',
                               "prop_sd_subclone:trans_site__Intercept" =  'Ramet:Site|prop',
                               "prop_sigma" =  'Residual|prop',
                               "prop_sigma_fixed" =  'Fixed Effects (mostly site)|prop',
                               "prop_sd_clone_total" =  'Genet (Total)|prop',
                               "prop_sd_subclone_total" =  'Ramet (Total)|prop',
                               "prop_sd_biol_total" =  'Genet+Ramet (Total)|prop',
                               "prop_sd_resid" =  'Residual+Rope|prop',
                               # These are for the terms instead
                               "age" = 'Age (f)',
                               "age:loc" = 'Age:Home/Away (f)',
                               "clone" = 'Genet',
                               "clone:trans_site" = 'Genet:Site' ,
                               "loc" = 'Home/Away (f)',
                               "rope" = 'Rope',
                               "subclone" = 'Ramet',
                               "subclone:trans_site" = 'Ramet:Site',
                               "trans_site" = 'Site (f)'
)

prep_sigma_data = \(data, resp = 'growth') {
  # browser()
  dens = data |> filter(y == resp) |>
    select( density) |>
    unnest(density) |>
    mutate(param = rename_variances(param)) |>
    filter(!str_detect(param, "sd$")) |>
    mutate(param = str_remove(param, fixed('|prop')))
  # browser()
  ci_tmp = data |> filter(y == resp) |> select( ci) |>
    unnest(ci) |> unnest(ci) |> # needs to be done twice for some reason
    mutate(param = rename_variances(param)) |>
    filter(!str_detect(param, "sd$")) |>
    mutate(param = str_remove(param, fixed('|prop'))) |>
    # Now let's rehsape for quantiles
    pivot_longer(starts_with('q'), names_to = 'quantile', values_to = 'val') |>
    mutate(CI = recode(quantile, q2.5 = "95%", q97.5 = "95%", q5 = "90%", q95 = "90%", q25 = "50%", q75 = "50%", q50 = "Median"),
           side = recode(quantile, q2.5 = "Lower", q97.5 = "Upper", q5 = "Lower", q95 = "Upper", q25 = "Lower", q75 = "Upper", q50 = "Median"),)

  medians = ci_tmp |> filter(CI == "Median") |> rename(x = val)
  ci = ci_tmp |> filter(CI != "Median") |> select(-quantile) |> pivot_wider(names_from = side, values_from = val)
  list(dens = dens, ci = ci, median = medians)
}

plot_sigmas = \(data, resp = 'growth') {
  .d = prep_sigma_data(data, resp)

  ggplot() + aes(y = param) +
    geom_slab(aes(x = x, thickness = density), normalize = 'xy', data = .d$dens, fill = grey(.8), color = grey(.4), linewidth = .5) +
    geom_interval(aes(xmin = Lower, xmax = Upper, linewidth = CI), orientation = 'horizontal', data = .d$ci) +
    scale_linewidth_manual('Credible Interval', values = c(1.75,1.25,.75))+
    geom_point(aes(x = x), data = .d$median, shape = "|", size = 3) +
    xlab("Portion of total variance explained") + ylab("Variance Component")
}

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

# Pull up the weight list ####
clean_summaries = read_rds('out/model_average_summaries.rds')
joint_weights = map_dfr(c('growth', 'chlor'), \(x) file.path('out', x, 'joint_weights.csv') |> read_csv() |> mutate(y = x))
source("R/model_configs.r")

model_terms = pick_model('all') |> # pick_model('all')
  map(\(x) x$formula$formula |> terms() |>  attr('term.labels') |> str_remove(fixed('1 | ')))
# Named list w/ terms for each mode.

models_per_term = tibble(model_name = names(model_terms)) |>
  mutate(terms = model_terms) |>
  # chop(terms) |>
  inner_join(joint_weights) |>
  unchop(terms) |>
  group_by(y, terms) |>
  summarize(posterior_prob = sum(posterior_prob)) |>
  mutate(terms = rename_variances(terms))
plot_model_terms = \(yvar = 'growth') {
  models_per_term |> filter(y == yvar) |>
    ggplot(aes(y = terms, x = posterior_prob)) +
    geom_col() +
    ylab('') +
    xlab("Model weight inclusion")
}
# Set the global themes
theme_set(theme_classic()) ; theme_update(
  axis.text.y = element_text(size = 10, color = 'black'),
  axis.text.x = element_text(size = 9.5),
  axis.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  plot.title = element_text(size = 11),
  plot.margin = margin(0,0,0,0),
  plot.background = element_blank(),
  plot.tag = element_text(size = 10, face = 'italic', hjust = 0, vjust = 1, color = grey(.2)),
  plot.tag.position = c(0,1)
)

# Make the sigma plots ####
growth_sigmas = plot_sigmas(clean_summaries, 'growth') + ggtitle("A")
growth_terms = plot_model_terms('growth') + ggtitle("C") # theme(plot.margin = margin(0,0, 0, -20))
chlor_sigmas = plot_sigmas(clean_summaries, 'chlor') + ggtitle("B") # ggtitle("Chlorophyll")
chlor_terms = plot_model_terms('chlor') + ggtitle("D") #+ theme(plot.margin = margin(0,0, 0, -20))

# joint_sigma_area = \(w_l, w_r, o = 0) {
#   c(
#     area(2, l = 1,       r = 1 + w_l),
#     area(2, l = 1+w_l-o, r = 1 + w_l-o+w_r),
#     area(1, l = 1,       r = 1 + w_l),
#     area(1, l = 1+w_l-o, r = 1 + w_l-o+w_r)
#   )
# }

joint_sigma_posterior = (growth_sigmas  + plot_spacer() + growth_terms + plot_spacer())  + chlor_sigmas + plot_spacer() +  chlor_terms + plot_spacer() +
  patchwork::plot_layout(ncol = 4, widths = c(5.5, -1.9,  3, 0.2), guides = 'collect') +
  # patchwork::plot_layout(design = joint_sigma_area(8, 4, -1), ncol = 2, guides = 'collect') +
  patchwork::plot_annotation(tag_levels = list(c("Growth", "", "Chlorophyll",  ""))) &
  theme(legend.position = 'bottom')

scl = 1; ggsave('figures/variance_components/joint_ensemble_sigmas.png',
                  joint_sigma_posterior, bg = 'white',
       dpi = 300, width = 7 * scl, height = 10.5/1.8 * scl, scale = scl)
scl = 1; ggsave('figures/variance_components/joint_ensemble_sigmas.pdf',
                joint_sigma_posterior, bg = 'white',
                dpi = 300, width = 7 * scl, height = 10.5/1.8 * scl, scale = scl)

# scl = 1.1; ggsave('figures/variance_components/joint_ensemble_sigmas.png',
#                   joint_sigma_posterior, bg = 'white',
#                   dpi = 300, width = 7 * scl, height = 10.5/2 * scl, scale = scl)
#

rename_fixed = \(x) recode(x,
                           "b_Intercept" = "Intercept (Site M1)",
                           "b_age" = "Age",
                           "b_age:locaway" = "Age : Location",
                           "b_locaway" = "Location (Away)",
                           "b_trans_siteM2" = "Transplant Site (M2 vs. M1)",
                           "b_trans_siteR" = "Transplant Site (R vs. M1)")

## Fixed Effect Plots ####

fixed_effect_summaries = read_rds('out/fixed_effect_summaries.rds')

prep_beta_data = \(data, resp = 'growth') {
  # browser()
  dens = data |> filter(y == resp) |>
    select( density) |>
    unnest(density) |>
    mutate(param = rename_fixed(param))
  # browser()
  ci_tmp = data |> filter(y == resp) |> select( ci) |>
    unnest(ci) |> unnest(ci) |> # needs to be done twice for some reason
    mutate(param = rename_fixed(param)) |>
    # Now let's rehsape for quantiles
    pivot_longer(starts_with('q'), names_to = 'quantile', values_to = 'val') |>
    mutate(CI = recode(quantile, q2.5 = "95%", q97.5 = "95%", q5 = "90%", q95 = "90%", q25 = "50%", q75 = "50%", q50 = "Median"),
           side = recode(quantile, q2.5 = "Lower", q97.5 = "Upper", q5 = "Lower", q95 = "Upper", q25 = "Lower", q75 = "Upper", q50 = "Median"),)

  medians = ci_tmp |> filter(CI == "Median") |> rename(x = val)
  ci = ci_tmp |> filter(CI != "Median") |> select(-quantile) |> pivot_wider(names_from = side, values_from = val)
  list(dens = dens, ci = ci, median = medians)
}

plot_betas = \(data, resp = 'growth') {
  .d = prep_beta_data(data, resp)

  ggplot() + aes(y = param) +
    geom_vline(xintercept = 0, linetype = 2, color = grey(.7)) +
    geom_slab(aes(x = x, thickness = density), normalize = 'xy', data = .d$dens, fill = grey(.8), color = grey(.4), linewidth = .15) +
    geom_interval(aes(xmin = Lower, xmax = Upper, linewidth = CI), orientation = 'horizontal', data = .d$ci) +
    scale_linewidth_manual('Credible Interval', values = c(1.75,1.25,.75))+
    geom_point(aes(x = x), data = .d$median, shape = "|", size = 3) +
    theme_classic() + xlab("Coefficient Estimate") + ylab("Parameter")
}

growth_betas = plot_betas(fixed_effect_summaries, 'growth') + ggtitle("Growth")
chlor_betas = plot_betas(fixed_effect_summaries, 'chlor') + ggtitle("Chlorophyll")

beta_plot = growth_betas / chlor_betas + patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = list(c("A", 'B'))) &
  theme(legend.position = 'bottom')
ggsave('figures/fixed_effects.png', beta_plot, width = 6, height = 6, dpi = 300)

# Make a table for assistance w/ writing

combined_table = bind_rows(fixed_effect_summaries, clean_summaries |> filter()) |>
  select(y, ci) |> unnest(ci) |> unnest(ci) |>
  mutate(param_name = param |> rename_fixed() |> rename_variances() |>
           str_remove(fixed('(fixed)')) |> str_remove(fixed('|prop')) |>
           str_remove(fixed('(w/o site)')) |> trimws(),
         terms = param_name |>
           str_remove(fixed('(Away)')) |>
           str_remove(fixed('(Total)')) |> trimws()) |>
  filter(str_detect(param, '^prop')|str_detect(param, '^b')) |>
  select(-param) |>
  left_join(models_per_term |>
              mutate(terms = terms |>
                       str_remove(fixed('(fixed)')) |>
                       str_replace('Home/Away', 'Location') |>
                       trimws()))
  # add in the percent inclusion
combined_table |>

  mutate(across(c(q50, q2.5, q97.5, posterior_prob), \(x) sprintf("%.1f%%", x*100))) |>
  # mutate(incl = sprintf(posterior_prob*100))
  mutate(inclusion = if_else(posterior_prob == 'NA%', '', paste0('; ', posterior_prob, ' inclusion')),
         txt = glue('{q50} [95% CI: {q2.5} - {q97.5}{inclusion}]')) |>
  select(y, param_name, txt) |> arrange(desc(y)) |>
  filter(y == 'growth') |> select(-y) |> as.data.frame()

