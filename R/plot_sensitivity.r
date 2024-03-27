# Plot prior sensitivity data
library(tidyverse)
library(ggdist)
theme_set(theme_classic())

chlor_data = dir('out/sensitivity/chlor', full.names = TRUE) |>
  map_dfr(read_rds) |>
  mutate(yvar = 'chlor')

growth_data = dir('out/sensitivity/growth', full.names = TRUE) |>
  map_dfr(read_rds) |>
  mutate(yvar = 'growth')


prep_sens_data = \(dat) {
  dat |> unnest_wider(summary) |>
    mutate(beta_med = qbeta(.5, a, b)) |>
    arrange(yvar, param, model, beta_med) |>
    group_by(yvar, param, model) |>
    mutate(x_rank = 1:n()) |>
    nest() |> ungroup()
}

dens_data = \(df, base_param) {
  # browser()
  df |> select(beta_med, x_rank, density, a, b, model) |>
    unnest(density) |>
    filter(param == paste('prop', base_param, sep = '_')) |>
    group_by(a, b, model) |>
    mutate(x_theor = seq(0, 1, length.out = n()),
           theor_dens = dbeta(x_theor, a, b))
}

sens_data = bind_rows(growth_data, chlor_data) |> prep_sens_data()
formatted_top_models = read_csv('formatted_top_models.csv')
all_models = read_csv('out/model_weights_with_formula.csv')
sens_data |> select(1:3) |> map(unique)

param_key = tribble(
  ~param, ~avar, ~bvar,
  'biol', 'Genet + Ramet', 'Incidental',
  'clone_site', 'Genet:Site', "Genet",
  "genet", "Genet", "Ramet",
  "subclone_site", "Ramet:Site", "Ramet"
)
sens_plot_data = sens_data |> left_join(param_key, 'param') |>
  rename(model_name = model) |>
  mutate(Response = recode(yvar, growth = 'Growth', chlor = 'Chlorophyll')) |>
  left_join(all_models, by = c('yvar' = 'y', 'model_name')) |>
  select(-model_name) |>
  arrange(rev(yvar), param, run_weight) |>
  mutate(formula = formula |>
           str_replace_all('subclone', 'Ramet') |>
           str_replace_all('clone', 'Genet') |>
           str_replace_all('trans_site', 'Site') |>
           str_replace('rope', "Rope") |>
           str_replace_all('loc', 'Location') |>
           str_replace_all('age', 'Age') |>
           str_wrap(50),
  # Fix formula format
         model = paste0(formula, "\n(", sprintf('%.1f%%', posterior_prob*100),
                        " posterior probability)" ) |>
                  fct_reorder(posterior_prob, .desc = TRUE)) |>
  select(param, Response, avar, bvar, model, data) |>
  unnest(data) |> nest(data = c(model, a, b, ci, density, beta_med, x_rank))


# Add in model inclusion proportion

# There's some duplication going on w/ the red lines...
plot_sens = \(data, Response, param, avar, bvar, ...) {
  # browser()
  data |> dens_data(base_param = param) |> ggplot() +
    geom_slab(aes(x = x_rank, y = x, thickness = density),
              fill = alpha('red', .5), color = NA) +
    geom_slab(aes(x = x_rank, y = x_theor, thickness = theor_dens),
              fill = NA, color = alpha("black", 0.8), linetype = 1, linewidth = .5) +
    geom_slab(aes(x = x_rank, y = x, thickness = density),
              fill = NA, color = "red", linewidth = .35) +
    ylab(glue("Partition of {bvar} vs. {avar} for {Response}")) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    scale_x_continuous(glue("Prior (black) and posterior (red) distributions\nof {avar} vs. {bvar} partition"),
                       breaks = c(3, 11, 17),
                       labels = paste(c(bvar, 'Even', avar),
                                      c("Biased", "Prior", "Biased"),
                                      sep = '\n')) +
    ylim(0, 1) + coord_flip() +
    facet_wrap(~model)
}

dir.create('figures/sensitivity')
sens_plots = sens_plot_data |> mutate(plot = pmap(cur_data(), plot_sens))

sens_plots |>
  mutate(filename = glue("figures/sensitivity/{Response}_{param}.png")) |>
  select(filename, plot) |>
  pwalk(ggsave, width = 12, height = 8, dpi = 300)


# Quick condensed version of sens_plot

condensed_data = sens_plot_data |>
  rename(focal_param = param) |>
  unnest(data) |> select(-ci) |>
  unnest(density) |>
  filter(paste0('prop_', focal_param) == param) |>
  filter(a == 10 | b == 10 | a == b) |>
  mutate(bvar = if_else(str_detect(avar, 'Site'), paste(bvar, '(Avg)'), bvar )) |>
  arrange(desc(a), b) |>
  mutate(contrast = paste(avar, 'vs.', bvar),
         prior = case_when(a == 10 ~ paste(avar, 'biased'),
                           b == 10 ~ paste(bvar, 'biased'),
                           a == b ~ 'Even') |> fct_inorder(),
         Response = as.factor(Response) |> fct_rev())
condensed_prior = condensed_data |>
  distinct(contrast, prior, Response, a, b) |>
  mutate(data = map2(a,b, \(A,B) tibble(x = seq(0, 1, length.out = 512)) |> mutate(density = dbeta(x, A, B) ) )) |>
  unnest(data)

condensed_sens_plot = ggplot(condensed_data, aes(x = x, y = prior)) +
  facet_grid(contrast ~ Response, scale = 'free_y', switch = 'y') +
  geom_slab(aes(thickness = density),
            fill = alpha('red', .1), color = NA) +
  geom_slab(aes(thickness = density), data = condensed_prior,
            fill = NA, color = alpha("black", 1), linetype = 1, linewidth = .5) +
  geom_slab(aes(thickness = density),
            fill = NA, color = alpha("red", .8), linewidth = .2) +
  geom_slab(aes(thickness = density), data = condensed_prior,
            fill = NA, color = alpha("black", 5), linetype = 1, linewidth = .5) +
  theme_classic() + theme(strip.background = element_blank(), strip.placement = 'outside') +
  ylab('') + xlab('Partition Proportion') + xlim(c(0, 1))
ggsave('figures/sensitivity/condensed.png', condensed_sens_plot, width = 6, height = 8, dpi = 300)
