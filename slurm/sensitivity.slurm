#!/bin/bash

#-------------------------------------------------------
#SBATCH -J sensitivity  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p development
#SBATCH -o logs/sens.o
#SBATCH -e logs/sens.e
#SBATCH -t 1:00:00

#------------------------------------------------------

ml launcher
module unload xalt

export LAUNCHER_JOB_FILE=$PWD/jobs/sensitivity.job
export LAUNCHER_WORKDIR=$PWD
$LAUNCHER_DIR/paramrun

