
# Bayesian Mixed Model analysis for “The role of intra-clonal variation in facilitating adaptation in a predominately asexual coral population”

These scripts were designed for use on the Lonestar 6 at the TExas
Advanced Computing Center (TACC); they should work on other HPC systems
with some minor adjustments.

## Setup

[Install the docker image](setup/docker) locally and on TACC. Clone the
git repository to TACC.

## Run Analyses

Run the mixed model analysis on TACC:

``` bash
sbatch slurm/mixed_models.slurm
```

After the mixed model job completes, run the sensitivity analysis

``` bash
sbatch slurm/sensitivity.slurm
```

## Create Plots

After mixed model and sensitivity jobs complete, open an interactive
session on TACC `(idev)` and run the following scripts:

``` bash
r-brms plot_model_averages.r # r-brms is Rscript wrapped in the docker image
r-brms plot_bivariate_posterior.r
r-brms plot_sensitivity.r
```

If you wish to do this locally, copy the files from `out` in TACC and
run these R scripts through the docker image’s version of RStudio.
