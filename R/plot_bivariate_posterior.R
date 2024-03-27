# N_DRAWS = draw_details %>% unchop(.draw) %>% nrow()
# Get the credible intervals of effect sizes & the frequency of the factor's inclusion in models

all_model_draws = read_rds('out/model_average_ensemble_draws.rds')

add_var = \(x, y) {
  x0 = replace_na(x, 0)
  y0 = replace_na(y, 0)
  out = sqrt(x0^2 + y0^2)
  out[is.na(x) & is.na(y)] = NA
  out
}
growth_sd_draws = all_model_draws |> filter(y == 'growth') |> pull(bf_draws) |> _[[1]] |> as_tibble() |> select(starts_with('sd')) |>
  rename_with(\(x) str_remove(x, '__Intercept') |> str_remove('sd_')) |>
  mutate(clone_total = add_var(clone, `clone:trans_site`), subclone_total = add_var(subclone, `subclone:trans_site`))
chlor_sd_draws = all_model_draws |> filter(y == 'chlor') |> pull(bf_draws) |> _[[1]] |> as_tibble() |> select(starts_with('sd')) |>
  rename_with(\(x) str_remove(x, '__Intercept') |> str_remove('sd_')) |>
  mutate(clone_total = add_var(clone, `clone:trans_site`), subclone_total = add_var(subclone, `subclone:trans_site`))


make_bivariate_plot = \(data, xcol, ycol, xname, yname,
                        xsplit = c('split', 'sole'), ysplit = c('sole', 'split'),
                        var_cutoff = 0.5, split_lab_pos = var_cutoff, aspect = 3,
                        dens_n = 400, dens_breaks = 2^(2:8), marg_pad = 0.1,
                        split_text_size = 4,
                        ...) {
  # Shoudl x & y margins be split by the other variable's presence/absence?
  xsplit = match.arg(xsplit)
  ysplit = match.arg(ysplit)

  # Adjust aspect for changing var_cutoff; it's normalized for 1 right now
  aspect = aspect / var_cutoff

# browser()
  # Make the data
  plot_dat = data |> select(xvar = {{xcol}}, yvar = {{ycol}}) |>
    mutate(across(xvar:yvar, \(z) z^2)) |> # Convert to variances
    # Remove double-missing variables; replace single missings w/ 0
    filter(!(is.na(xvar) & is.na(yvar))) |>
    replace_na(list(xvar = 0, yvar = 0))

  density_2d_data = plot_dat |>
    # Density estimate
    with(MASS::kde2d(xvar, yvar, n = dens_n)) |>
    # convert to data frame
    with({
      # Reminder: x and y aren't the same as xvar and yvar
      n = length(x)
      tibble(x = rep(x, times = n),
             y = rep(y, each = n),
             z = as.numeric(z))
    }) |>
    mutate(alpha = if_else(z < 1, sqrt(z), 1))


  # marginal plots
  split_nm = \(x) paste(c("Without", "With"), x)
  make_sole_marginal = \(df, var, ...) {
    df |>
      select(v = {{var}}) |> filter(v != 0) |>
      ggplot() +  ggdist::stat_halfeye(#aes(side = side),
                                       scale = .8, n = 10000, side = 'topleft',
                                       slab_linewidth = 0.25, slab_color = grey(.4),
                                       slab_linetype = 1, slab_fill = alpha(grey(.6), .5))
  }
  make_split_marginal = \(df, var, split_by, split_name) {
    dat = df |> select(v = {{var}}, s = {{split_by}}) |>
      mutate(side = if_else(s == 0, 'top', 'bottom')) |>
      ggplot() +
      ggdist::stat_halfeye(aes(side = side), scale = .8, n = 10000,
                           slab_color = grey(.4), slab_linewidth = 0.25,
                           slab_linetype = 1, slab_fill = alpha(grey(.6), .5))+
      scale_side_mirrored(guide = 'none') + theme(legend.position = 'none')
  }
  x_marg_theme = \(ref_lim = 1) {
    # ref_lim is the lower limit for the y axis; It may be different for sole vs. split margins
    list(
      theme(axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(),
            axis.text.x = element_text(size = 11),
            axis.title.x = element_text(size = 5 * .pt)), # scale equivalent to geom_text s
      scale_x_continuous(expand = c(0,0)),
      scale_y_discrete(expand = expansion(add = .05)),
      coord_fixed(xlim = c(0, var_cutoff), ylim = c(ref_lim - marg_pad, ref_lim + 1 + marg_pad), ratio = 1/aspect, expand = FALSE)
  )  }
  y_marg_theme = \(ref_lim = 1) {
    list(
      theme(axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(),
            axis.text.y = element_text(size = 11),
            axis.title.y = element_text(size = 5 * .pt)),
      scale_y_continuous(expand = c(0,0)),
      scale_x_discrete(expand = expansion(add = .05)),
      coord_fixed(ylim = c(0, var_cutoff), xlim = c(ref_lim  - marg_pad, ref_lim + 1 + marg_pad), ratio = aspect, expand = FALSE)
    ) }
# browser()
  make_marginal = list(
    x = list(
      split = \(df, var, split_by, split_name, ..., lab_pos = 1.5 ) {

    make_split_marginal(df, enexpr(var), enexpr(split_by), split_name) +
      aes(x = v, y = side) +
      x_marg_theme(1) +
      annotate("text", x = lab_pos, y = 1, label = split_nm(split_name)[1], hjust = 1, vjust = -0.3, size = split_text_size) +
      annotate("text", x = lab_pos, y = 2, label = split_nm(split_name)[2], hjust = 1, vjust = 1.3, size = split_text_size)
  },
      sole  = \(df, var, ...) {
    make_sole_marginal(df, enexpr(var)) +  aes(x = v) + x_marg_theme(-1)
  }
    ),
    y = list(
      split = \(df, var, split_by, split_name, ..., lab_pos = 1.5 ) {

      make_split_marginal(df, enexpr(var), enexpr(split_by), split_name) +
        aes(y = v, x = side) +
        y_marg_theme(1) +
        annotate("text", y = lab_pos, x = 1,angle = 90, label = split_nm(split_name)[1], vjust = 1.3, hjust = 1, size = split_text_size) +
        annotate("text", y = lab_pos, x = 2, angle = 90, label = split_nm(split_name)[2], vjust = -0.3, hjust = 1, size = split_text_size)
    },
      sole = \(df, var, ...) {
        make_sole_marginal(df, enexpr(var)) +  aes(y = v) + y_marg_theme(-1)
    }
    )
  )
  ## 2d density plot

  make_density_plot = \(data) {
    data |>
      mutate(alpha = pmin(z, 1),
             z = pmax(.5, z)) |>
      ggplot(aes(x, y, fill = z,  alpha = alpha)) +
      geom_raster() +
      scale_fill_viridis_c("Posterior\nDensity", trans =  "log1p",
                           breaks = dens_breaks,
                           rescaler =~pmax(.x, 1) |>  scales::rescale(),
                           direction = -1, ) +
      scale_alpha_continuous(range = c(0, 1),
                             # rescaler = ~if_else(.x < 1, .x^2, 1),
                             guide = "none") +
      theme_classic() +
      # scale_fill_viridis_c(direction = -1) +
      coord_fixed(expand = FALSE, ylim = c(0, var_cutoff), xlim = c(0, var_cutoff)) +
      theme(legend.position = c(1, 1), legend.justification = c(1,1),
          axis.title = element_blank(),
          plot.margin = margin(),
          legend.box.background = element_rect(fill = "#00000000", color = "#00000000"),
          legend.background = element_rect(fill = "#00000000", color = "#00000000"),
          axis.text = element_blank())
  }

  xmarg_plot = make_marginal$x[[xsplit]](df = plot_dat, var = xvar, split_by = yvar, split_name = yname, lab_pos = split_lab_pos) +
    xlab(xname)
  ymarg_plot = make_marginal$y[[ysplit]](df = plot_dat, var = yvar, split_by = xvar, split_name = xname, lab_pos = split_lab_pos) +
    ylab(yname)
  joint_plot = make_density_plot(density_2d_data)
  # Render the plots

  out_plot = ymarg_plot + joint_plot +
    plot_spacer() + xmarg_plot +
    plot_layout(ncol = 2)
  out_plot
}
theme_set(theme_classic())

# Growth model

make_tags = \(letters = LETTERS[1:3]) {
  map(letters, \(x) c(x, '', '')) |> unlist() |> list()
  # we want c('A', '', ''), c('B', '', ''), ect
}


growth_ramet = make_bivariate_plot(growth_sd_draws, xcol = subclone, ycol = `subclone:trans_site`, xname = 'Ramet (Average)',
                    yname = 'Ramet : Site', aspect = 4, dens_breaks = 2^(1:8))
growth_genet = make_bivariate_plot(growth_sd_draws, xcol = clone, ycol = `clone:trans_site`, xname = 'Genet (Average)',
                    yname = 'Genet : Site', var_cutoff = .25, dens_n = 1000, dens_breaks = 4^(1:6))

growth_gen_ram = make_bivariate_plot(growth_sd_draws, xcol = clone_total, ycol = subclone_total, xname = 'Genet (Total)',
                    yname = 'Ramet (Total)',  ysplit = 'split', var_cutoff = .5, dens_n = 1000, dens_breaks = 4^(1:6))
# all_growth = (((growth_ramet) / (growth_genet) / (growth_gen_ram))  )
# ggsave('figures/bivariate_variances/growth_trio.png', all_growth, width = 5, height = 11, dpi = 300)



chlor_ramet = make_bivariate_plot(chlor_sd_draws, xcol = subclone, ycol = `subclone:trans_site`, xname = 'Ramet (Average)',
                    yname = 'Ramet : Site', var_cutoff = 0.3, aspect = 4, dens_breaks = 4^(1:5))
chlor_genet = make_bivariate_plot(chlor_sd_draws, xcol = clone, ycol = `clone:trans_site`, xname = 'Genet (Average)',
                    yname = 'Genet : Site', dens_n = 1000, var_cutoff = .35,
                    aspect = 4, dens_breaks = 4^(1:6))

chlor_gen_ram = make_bivariate_plot(chlor_sd_draws, xcol = clone_total, ycol = subclone_total, xname = 'Genet (Total)',
                    yname = 'Ramet (Total)',  ysplit = 'split',
                    dens_n = 1000, var_cutoff = .45, aspect = 4, dens_breaks = 2^(1:8))

# all_chlor = (chlor_ramet) / (chlor_genet) / (chlor_gen_ram)
# ggsave('figures/bivariate_variances/chlor_trio.png', all_chlor, width = 5, height = 12, dpi = 300)

plot_list = list(
  A = growth_gen_ram, B = growth_ramet,  C= growth_genet,
  D = chlor_gen_ram, E = chlor_ramet, F = chlor_genet
)
dsgn = "
AD
BE
CF"
# full_plot = (all_growth) | (all_chlor) + plot_annotation(tag_levels = make_tags(LETTERS[1:6]))

full_plot = wrap_plots(plot_list, ncol = 2, byrow = FALSE, design = dsgn) +
  plot_annotation(tag_levels = make_tags(LETTERS[1:6]))

ggsave('figures/bivariate_variances/full.png', full_plot, width = 9, height = 12, dpi = 300)


