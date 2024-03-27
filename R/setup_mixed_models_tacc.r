suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(rlang)
  library(purrr)
  library(glue)
})

if(!exists('argv')) argv = commandArgs(TRUE)
job_file = argv[1] %|% 'jobs/run_mixed_models.job'
model_list = c("null_full_allbiol", "fixed_full_allbiol", "null_nint_allbiol", "fixed_nint_allbiol",
             "null_nage_allbiol", "fixed_nage_allbiol", "null_nloc_allbiol", "fixed_nloc_allbiol",
             "null_nage_nloc_allbiol", "fixed_nage_nloc_allbiol", "null_full_nsclone", "null_nint_nsclone",
             "null_nage_nsclone", "null_nloc_nsclone", "null_nage_nloc_nsclone", "fixed_nage_nsclone",
             "fixed_nloc_nsclone", "fixed_full_nsclone", "fixed_nint_nsclone", "fixed_nage_nloc_nsclone",
             "null_full_nssub", "null_nint_nssub", "null_nage_nssub", "null_nloc_nssub", "null_nage_nloc_nssub",
             "fixed_nage_nssub", "fixed_nloc_nssub", "fixed_full_nssub", "fixed_nint_nssub", "fixed_nage_nloc_nssub",
             "null_full_nseither", "null_nint_nseither", "null_nage_nseither", "null_nloc_nseither",
             "null_nage_nloc_nseither", "fixed_nage_nseither", "fixed_nloc_nseither", "fixed_full_nseither",
             "fixed_nint_nseither", "fixed_nage_nloc_nseither", "null_full_nclone", "null_nint_nclone",
             "null_nage_nclone", "null_nloc_nclone", "null_nage_nloc_nclone", "fixed_nage_nclone", "fixed_nloc_nclone",
             "fixed_full_nclone", "fixed_nint_nclone", "fixed_nage_nloc_nclone", "null_full_nsub", "null_nint_nsub",
             "null_nage_nsub", "null_nloc_nsub", "null_nage_nloc_nsub", "fixed_nage_nsub", "fixed_nloc_nsub",
             "fixed_full_nsub", "fixed_nint_nsub", "fixed_nage_nloc_nsub", "null_full_nclone_nssub",
             "null_nint_nclone_nssub", "null_nage_nclone_nssub", "null_nloc_nclone_nssub", "null_nage_nloc_nclone_nssub",
             "fixed_nage_nclone_nssub", "fixed_nloc_nclone_nssub", "fixed_full_nclone_nssub", "fixed_nint_nclone_nssub",
             "fixed_nage_nloc_nclone_nssub", "null_full_nsub_nsclone", "null_nint_nsub_nsclone", "null_nage_nsub_nsclone",
             "null_nloc_nsub_nsclone", "null_nage_nloc_nsub_nsclone", "fixed_nage_nsub_nsclone",
             "fixed_nloc_nsub_nsclone", "fixed_full_nsub_nsclone", "fixed_nint_nsub_nsclone",
             "fixed_nage_nloc_nsub_nsclone", "null_full_nbiol", "null_nint_nbiol", "null_nage_nbiol",
             "null_nloc_nbiol", "null_nage_nloc_nbiol", "fixed_nage_nbiol", "fixed_nloc_nbiol", "fixed_full_nbiol",
             "fixed_nint_nbiol", "fixed_nage_nloc_nbiol")
model_list_fixed = str_subset(model_list,  'fixed_')
model_list_null = str_subset(model_list, 'null_')
# Split them out because I first ran this on just the fixed, so this way they'll keep the same seeds & null will get new seeds.
full_model_list = c(model_list_fixed, model_list_null)

set.seed(23902)
seeds = rdunif(length(full_model_list), 1, 65536*1024)
cmds = glue::glue("r-brms R/batch_mixed_models.r {full_model_list} {seeds}")
dir.create('jobs', showWarnings = FALSE)
dir.create('logs', showWarnings = FALSE)
cmds |>
  write_lines(job_file)

