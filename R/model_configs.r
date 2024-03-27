set_formula = \( ...) {
  formula_sets = list(
    # No rope, it's always in
    fixed_site = .~  trans_site+.,
    rand_site = .~  (1|trans_site)+.,
    clone = .~(1 | clone)+.,
    clone_site = .~(1  | clone:trans_site) +.,
    # clone_site = .~(1 + trans_site | clone) +.,
    subclone = .~(1 | subclone)+.,
    subclone_site = .~(1  | subclone:trans_site) +.,
    # subclone_site = .~(1 + trans_site | subclone) +.,
    age = .~ . + age,
    loc = .~ . + loc,
    age_loc = .~. + age:loc
  )
  # ... is a list of terms (besides rope)
  lhs = y ~ (1|rope)
  terms = rlang::ensyms(...) |> map_chr(as.character)
  # browser()
  terms = match.arg(terms, names(formula_sets), several.ok = TRUE)
  # browser()
  list(form = bf(reduce(formula_sets[terms], update.formula, .init = lhs)),
       terms = terms)
}

# Set up the model structure to match the formula
# Return object should be everything needed to run the model, plus an hp() function for setting hyperpriors
set_params = \(formula_set) {
  # define specifications from terms

  add_if = \(spec, terms, true, false = NULL) {
    if(all(terms %in% formula_set$terms)) c(spec, true) else c(spec, false)
  }
  # Determine model components vs. specified terms
  specs = 'always' |> add_if(c('clone', 'subclone'), 'genet') |>
    add_if('clone_site', c('clone_site')) |>
    add_if('subclone_site', 'subclone_site') |>
    add_if('clone', c('biol', 'clone')) |> add_if('subclone', c('biol', 'subclone')) |>
    # add_if('rand_site', 'rope_site', 'just_rope') |>
    add_if('fixed_site', 'fixed_site') |>
    add_if('age', 'age') |>
    add_if('loc', 'loc') |>
    add_if('age_loc', 'age_loc') |>
    unique()

  spec_subset = \(lst) {
    idx = specs[specs %in% names(lst)]
    lst[idx]
  }
  stanvar_list = list(
    always = stanvar(scode = 'real<lower=0> sigma_total;', block = 'parameters') +
      stanvar(scode = 'lprior += student_t_lpdf(sigma_total | hp_sigma_nu, 0, hp_sigma_scale) - student_t_lccdf(0 | hp_sigma_nu, 0, hp_sigma_scale);',
              block = 'tparameters', position = 'end') +
      stanvar(scode = 'simplex[n_rand] prop_rand;', block = 'parameters') +
      stanvar(scode = 'lprior += dirichlet_lpdf(prop_rand | hp_rand);',  block = 'tparameters', position = 'end'),
    biol = stanvar(scode = 'real<lower=0,upper=1> prop_biol;', block = 'parameters') +
      stanvar(scode = 'lprior += beta_lpdf(prop_biol | hp_biol_a, hp_biol_b);',  block = 'tparameters', position = 'end'),
    genet = stanvar(scode = 'real<lower=0,upper=1> prop_genet;', block = 'parameters') +
      stanvar(scode = 'lprior += beta_lpdf(prop_genet | hp_genet_a, hp_genet_b);',  block = 'tparameters', position = 'end'),
    clone_site = stanvar(scode = 'real<lower=0,upper=1> prop_clone_site;', block = 'parameters') +
      stanvar(scode = 'lprior += beta_lpdf(prop_clone_site | hp_clone_site_a, hp_clone_site_b);',  block = 'tparameters', position = 'end'),
    subclone_site = stanvar(scode = 'real<lower=0,upper=1> prop_subclone_site;', block = 'parameters') +
      stanvar(scode = 'lprior += beta_lpdf(prop_subclone_site | hp_subclone_site_a, hp_subclone_site_b);',  block = 'tparameters', position = 'end')
  )

  # We have a few cases:
  # Add prop_genet_site and prop_ramet_site

  # Itteratively build up stan code based on conditions
  add_component = \(x, component, cond, code) {
    if(isTRUE(cond)) x[[component]] = paste(x[[component]], code)
    x
  }

  opts = list(clone = '', subclone = '', biol = '', clone_site = '', subclone_site = '') |>
    add_component('biol', 'biol' %in% specs, 'sqrt(1-prop_biol) *') |>
    add_component('clone', 'genet' %in% specs, "* sqrt(prop_genet)") |>
    add_component('subclone', 'genet' %in% specs, "* sqrt(1-prop_genet)") |>
    add_component('clone_site', 'genet' %in% specs, "* sqrt(prop_genet)") |>
    add_component('subclone_site', 'genet' %in% specs, "* sqrt(1-prop_genet)") |>

    add_component('clone', 'clone_site' %in% specs, "* sqrt(1-prop_clone_site)") |>
    add_component('subclone', 'subclone_site' %in% specs, "* sqrt(1-prop_subclone_site)") |>
    add_component('clone_site', 'clone_site' %in% specs, "* sqrt(prop_clone_site)") |>
    add_component('subclone_site', 'subclone_site' %in% specs, "* sqrt(prop_subclone_site)")


  # browser()
  prior_list = list(
    always = set_prior('normal(0, hp_intercept)',  class = 'Intercept') +
      set_prior(glue('constant({opts$biol}sqrt(prop_rand[1])  * sigma_total)'), class = 'sigma') +
      set_prior(glue('constant({opts$biol}sqrt(prop_rand[2])  * sigma_total)'), class = 'sd', group = 'rope'),
    age = set_prior('normal(0, hp_age)',  class = 'b', coef = 'age'),
    loc = set_prior('normal(0, hp_loc)',  class = 'b', coef = 'locaway'),
    age_loc = set_prior('normal(0, hp_age_loc)',  class = 'b', coef = 'age:locaway'),
    fixed_site =  set_prior('normal(0, hp_fixed_site)',  class = 'b', coef = paste0('trans_site', c('M2', 'R'))),
    rand_site = set_prior(glue('constant({opts$biol}sqrt(prop_rand[3])  * sigma_total)'),
                          class = 'sd', group = 'trans_site'),
    clone =  set_prior(glue('constant(sqrt(prop_biol) {opts$clone} * sigma_total)'), class = 'sd', group = 'clone'),
    subclone = set_prior(glue('constant(sqrt(prop_biol) {opts$subclone} * sigma_total)'), class = 'sd', group = 'subclone'),
    clone_site =  set_prior(glue('constant(sqrt(prop_biol) {opts$clone_site} * sigma_total)'), class = 'sd', group = 'clone:trans_site'),
    subclone_site = set_prior(glue('constant(sqrt(prop_biol) {opts$subclone_site} * sigma_total)'), class = 'sd', group = 'subclone:trans_site')
  )
  # browser()
  n_rand = as.integer(2L + ('rand_site' %in% formula_set$terms))

  # Add log_lik in cmdstan to gc
  # rand_terms = glue::glue("r_{i}_1[J_{i}[n]] * Z_{i}_1[n]", i = 1:n_rand) |> glue::glue_collapse(sep = ' + ')
  # log_lik = stanvar(scode = glue::glue("vector[N] log_lik;
  #       {
  #       vector[N] mu_base = Intercept + Xc * b;
  #       for (n in 1:N) {
  #         real mu = mu_base[n]  +  <<rand_terms>>;
  #         log_lik[n] = normal_lpdf(Y[n] | mu, sigma);
  #
  #       }
  #        }", .open = "<<", .close = ">>"),
  #       block  = 'genquant')

  # I'm getting a weird error about sizes being vectors.  NOt sure what the deal w/ it is.


  stanvars = spec_subset(stanvar_list) |> reduce(`+`) #+ log_lik
  priors = spec_subset(prior_list) |> reduce(`+`)

  # Set default hyper priors
  hp_list = list(
    always = list(hp_intercept = 2, hp_sigma_nu = 7, hp_sigma_scale=1,
                  n_rand = n_rand, hp_rand = rep(1, n_rand)),
    biol = list(hp_biol_a = 1, hp_biol_b = 1),
    genet = list(hp_genet_a = 1, hp_genet_b = 1),
    age = list(hp_age = 1),
    loc = list(hp_loc = 1),
    age_loc = list(hp_age_loc = .5),
    clone_site = list(hp_clone_site_a = 1, hp_clone_site_b = 1),
    subclone_site = list(hp_subclone_site_a = 1, hp_subclone_site_b = 1),
    fixed_site = list(hp_fixed_site = 1.5)
  )
  hp_defaults = spec_subset(hp_list) |> flatten()
  # Create a function that will auto set them but allow for easy modification:
  hp_function = local({
    hp_selfnamed = names(hp_defaults) |> map(rlang::sym) |> set_names(names(hp_defaults))
    hp_self_call = rlang::call2(list, !!!hp_selfnamed)
    hp_body = rlang::expr(imap(!!hp_self_call, stanvar) |> reduce(`+`))
    rlang::new_function(args = hp_defaults,  body = hp_body)
  })
  model_name = paste0('model_', paste0(specs, collapse = '_'))
  list(formula = formula_set$form,
       prior = priors, base_stanvar = stanvars,
       model_name = model_name,
       hp = hp_function)
}


pick_data = \(name, data_list = datasets) {
  set = data_key |> filter(y == name) |> pull(dataset)
  data_list[[set]] |> mutate(y = !!ensym(name))
}
pick_model = \(name = model_name) {
  models = list(
    null_full_allbiol    = set_formula( subclone, clone, age, loc, age_loc, clone_site, subclone_site)   ,
    fixed_full_allbiol         = set_formula(fixed_site, subclone, clone, age, loc, age_loc, clone_site, subclone_site)  ,
    # no interaction
    null_nint_allbiol          = set_formula( subclone, clone, age, loc, clone_site, subclone_site) ,
    fixed_nint_allbiol         = set_formula(fixed_site, subclone, clone,age, loc, clone_site, subclone_site) ,
    # no age
    null_nage_allbiol          = set_formula( subclone, clone, loc, clone_site, subclone_site) ,
    fixed_nage_allbiol         = set_formula(fixed_site, subclone, clone,loc, clone_site, subclone_site) ,
    # no loc
    null_nloc_allbiol          = set_formula( subclone, clone, age,  clone_site, subclone_site) ,
    fixed_nloc_allbiol         = set_formula(fixed_site, subclone, clone,age,  clone_site, subclone_site) ,
    # no age or loc
    null_nage_nloc_allbiol    = set_formula( subclone, clone,  clone_site, subclone_site) ,
    fixed_nage_nloc_allbiol   = set_formula(fixed_site, subclone, clone, clone_site, subclone_site) ,

    # Less full versions; consolidated for ease of column editing

    # nsclone: No site:clone
    null_full_nsclone       = set_formula(  subclone, clone,  subclone_site, age, loc, age_loc     )   ,
    null_nint_nsclone       = set_formula(  subclone, clone,  subclone_site, age, loc     ) ,
    null_nage_nsclone       = set_formula(  subclone, clone,  subclone_site, loc     ) ,
    null_nloc_nsclone       = set_formula(  subclone, clone,  subclone_site, age     ) ,
    null_nage_nloc_nsclone  = set_formula(  subclone, clone,  subclone_site     ) ,
    fixed_nage_nsclone      = set_formula(fixed_site, subclone, clone,  subclone_site, age, loc, age_loc     )   ,
    fixed_nloc_nsclone      = set_formula(fixed_site, subclone, clone,  subclone_site, age, loc     ) ,
    fixed_full_nsclone      = set_formula(fixed_site, subclone, clone,  subclone_site, loc     ) ,
    fixed_nint_nsclone      = set_formula(fixed_site, subclone, clone,  subclone_site, age     ) ,
    fixed_nage_nloc_nsclone = set_formula(fixed_site, subclone, clone,  subclone_site     ) ,
    # nssub: No site:subclone
    null_full_nssub       = set_formula(  subclone, clone, clone_site, age, loc, age_loc     )   ,
    null_nint_nssub       = set_formula(  subclone, clone, clone_site, age, loc     ) ,
    null_nage_nssub       = set_formula(  subclone, clone, clone_site, loc     ) ,
    null_nloc_nssub       = set_formula(  subclone, clone, clone_site, age     ) ,
    null_nage_nloc_nssub  = set_formula(  subclone, clone, clone_site     ) ,
    fixed_nage_nssub      = set_formula(fixed_site, subclone, clone, clone_site, age, loc, age_loc     )   ,
    fixed_nloc_nssub      = set_formula(fixed_site, subclone, clone, clone_site, age, loc     ) ,
    fixed_full_nssub      = set_formula(fixed_site, subclone, clone, clone_site, loc     ) ,
    fixed_nint_nssub      = set_formula(fixed_site, subclone, clone, clone_site, age     ) ,
    fixed_nage_nloc_nssub = set_formula(fixed_site, subclone, clone, clone_site     ) ,
    # Nseither: No site:clone or site:subclone
    null_full_nseither       = set_formula(  subclone, clone, age, loc, age_loc     )   ,
    null_nint_nseither       = set_formula(  subclone, clone, age, loc     ) ,
    null_nage_nseither       = set_formula(  subclone, clone, loc     ) ,
    null_nloc_nseither       = set_formula(  subclone, clone, age     ) ,
    null_nage_nloc_nseither  = set_formula(  subclone, clone     ) ,
    fixed_nage_nseither      = set_formula(fixed_site, subclone, clone, age, loc, age_loc     )   ,
    fixed_nloc_nseither      = set_formula(fixed_site, subclone, clone, age, loc     ) ,
    fixed_full_nseither      = set_formula(fixed_site, subclone, clone, loc     ) ,
    fixed_nint_nseither      = set_formula(fixed_site, subclone, clone, age     ) ,
    fixed_nage_nloc_nseither = set_formula(fixed_site, subclone, clone     ) ,

    # nclone: No clone
    null_full_nclone       = set_formula(  subclone, subclone_site, age, loc, age_loc     )   ,
    null_nint_nclone       = set_formula(  subclone, subclone_site, age, loc     ) ,
    null_nage_nclone       = set_formula(  subclone, subclone_site, loc     ) ,
    null_nloc_nclone       = set_formula(  subclone, subclone_site, age     ) ,
    null_nage_nloc_nclone  = set_formula(  subclone, subclone_site     ) ,
    fixed_nage_nclone      = set_formula(fixed_site, subclone, subclone_site, age, loc, age_loc     )   ,
    fixed_nloc_nclone      = set_formula(fixed_site, subclone, subclone_site, age, loc     ) ,
    fixed_full_nclone      = set_formula(fixed_site, subclone, subclone_site, loc     ) ,
    fixed_nint_nclone      = set_formula(fixed_site, subclone, subclone_site, age     ) ,
    fixed_nage_nloc_nclone = set_formula(fixed_site, subclone, subclone_site     ) ,
    # nsub: No subclone
    null_full_nsub       = set_formula(   clone, clone_site, age, loc, age_loc     )   ,
    null_nint_nsub       = set_formula(   clone, clone_site, age, loc     ) ,
    null_nage_nsub       = set_formula(   clone, clone_site, loc     ) ,
    null_nloc_nsub       = set_formula(   clone, clone_site, age     ) ,
    null_nage_nloc_nsub  = set_formula(   clone, clone_site     ) ,
    fixed_nage_nsub      = set_formula(fixed_site,  clone, clone_site, age, loc, age_loc     )   ,
    fixed_nloc_nsub      = set_formula(fixed_site,  clone, clone_site, age, loc     ) ,
    fixed_full_nsub      = set_formula(fixed_site,  clone, clone_site, loc     ) ,
    fixed_nint_nsub      = set_formula(fixed_site,  clone, clone_site, age     ) ,
    fixed_nage_nloc_nsub = set_formula(fixed_site,  clone, clone_site     ) ,

    # nclone_nssub: No clone, no site:subclone
    null_full_nclone_nssub       = set_formula(  subclone, age, loc, age_loc     )   ,
    null_nint_nclone_nssub       = set_formula(  subclone, age, loc     ) ,
    null_nage_nclone_nssub       = set_formula(  subclone, loc     ) ,
    null_nloc_nclone_nssub       = set_formula(  subclone, age     ) ,
    null_nage_nloc_nclone_nssub  = set_formula(  subclone     ) ,
    fixed_nage_nclone_nssub      = set_formula(fixed_site, subclone, age, loc, age_loc     )   ,
    fixed_nloc_nclone_nssub      = set_formula(fixed_site, subclone, age, loc     ) ,
    fixed_full_nclone_nssub      = set_formula(fixed_site, subclone, loc     ) ,
    fixed_nint_nclone_nssub      = set_formula(fixed_site, subclone, age     ) ,
    fixed_nage_nloc_nclone_nssub = set_formula(fixed_site, subclone     ) ,
    # nsub_nsclone: No subclone, no site:clone
    null_full_nsub_nsclone       = set_formula(   clone, age, loc, age_loc     )   ,
    null_nint_nsub_nsclone       = set_formula(   clone, age, loc     ) ,
    null_nage_nsub_nsclone       = set_formula(   clone, loc     ) ,
    null_nloc_nsub_nsclone       = set_formula(   clone, age     ) ,
    null_nage_nloc_nsub_nsclone  = set_formula(   clone     ) ,
    fixed_nage_nsub_nsclone      = set_formula(fixed_site,  clone, age, loc, age_loc     )   ,
    fixed_nloc_nsub_nsclone      = set_formula(fixed_site,  clone, age, loc     ) ,
    fixed_full_nsub_nsclone      = set_formula(fixed_site,  clone, loc     ) ,
    fixed_nint_nsub_nsclone      = set_formula(fixed_site,  clone, age     ) ,
    fixed_nage_nloc_nsub_nsclone = set_formula(fixed_site,  clone     ) ,

    # nbiol: No clone or subclone
    null_full_nbiol       = set_formula( age, loc, age_loc     )   ,
    null_nint_nbiol       = set_formula( age, loc     ) ,
    null_nage_nbiol       = set_formula( loc     ) ,
    null_nloc_nbiol       = set_formula( age     ) ,
    # null_nage_nloc_nbiol  = set_formula(null_site     ) ,
    fixed_nage_nbiol      = set_formula(fixed_site, age, loc, age_loc     )   ,
    fixed_nloc_nbiol      = set_formula(fixed_site, age, loc     ) ,
    fixed_full_nbiol      = set_formula(fixed_site, loc     ) ,
    fixed_nint_nbiol      = set_formula(fixed_site, age     ) ,
    fixed_nage_nloc_nbiol = set_formula(fixed_site     )
  )
  if(name == 'all') return(models |> map(set_params))
  models[[name]] |> set_params()
}
