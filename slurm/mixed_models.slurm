#!/bin/bash

#-------------------------------------------------------
# SBATCH -J mixed_models  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 4
#SBATCH -n 45
#SBATCH -p development
#SBATCH -o logs/mixed_models.o
#SBATCH -e logs/mixed_models.e
#SBATCH -t 2:00:00

#------------------------------------------------------

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=$PWD/jobs/mixed_models.job

# Setup the model runs
r-brms R/setup_mixed_models_tacc.r $LAUNCHER_JOB_FILE
export LAUNCHER_WORKDIR=$PWD
# Run the models
$LAUNCHER_DIR/paramrun
# Combine the results & setup the sensitivty job
r-brms R/summarize_model_weights.r $PWD/jobs/sensitivity.job
# Run the model averaging
r-brms R/model_average.r
