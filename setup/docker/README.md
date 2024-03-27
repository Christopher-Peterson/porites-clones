Docker Setup
================

# Using the docker image locally

To setup the docker image for the first time

``` bash
TAG=crpeters/r-brms:4.3.1
# Download the image
sudo docker pull $TAG

# Configure the image:
# Set the working directory as your current directory (the project root)
LOCAL_DIR=$PWD 
CONF_DIR=$LOCAL_DIR/setup/docker/prefs
# If you'd like to access the git repo from within the image, mount your ssh directory
# Otherwise, comment this line out
SSH_MOUNT="-v $HOME/.ssh:/home/rstudio/.ssh "
CMDSTAN_DIR="$HOME/cmdstan" # Local directory where cmdstan will be installed.
NAME=leaf_miner
PROJECT=leaf_miner
# uncomment to link stanCompose directory
STAN_COMPOSE="-v $LOCAL_DIR/../stanCompose:/home/rstudio/stanCompose "

USER_PASSWORD="1234" # if you wish to set a password for the rstudio account; otherwise, comment out the next line
SET_PASSWORD="-e PASSWORD=${USER_PASSWORD}"
PORT=8781 # port in web browser
# Run the container
sudo docker run -d -p 127.0.0.1:$PORT:8787 \
  -v $LOCAL_DIR:/home/rstudio/$PROJECT \
  -v $CMDSTAN_DIR:/home/rstudio/cmdstan \
  -v $CONF_DIR:/home/rstudio/.config/rstudio \
  $STAN_COMPOSE \ 
   $SSH_MOUNT \
   $SET_PASSWORD \
  -e DISABLE_AUTH=true -e USERID=$UID --name $NAME $TAG
```

To open Rstudio, go to [localhost:8781](localhost:8781) in a browser.

# Using the docker image on an HPC with Apptainer

## Installation on TACC

This image is designed to work on TACC (or other HPC systems,
potentially with some modifications) via Singularity/Apptainer.

Login to TACC and open an interactive session with `idev`. Then enter
the following:

``` bash
mkdir -p $WORK/apptainer/bin
mkdir -p $SCRATCH/apptainer/bin
# Add $Scratch/apptainer/bin to your $PATH in .bashrc

idev
cdw apptainer
module load tacc-apptainer 

apptainer pull -F docker://crpeters/r-leafminer:4.3.1 # CHANGE THIS
cp r-*.sif $SCRATCH/apptainer
```

Next, add `$WORK/singularity/bin` to your `$PATH` (usually in your
`.bashrc`). Upload the `singularity_mask.sh` script to
`$WORK/singularity` and use it to create a wrapper for Rscript that uses
the docker image:

``` bash
# still in $SCRATCH/apptainer
source apptainer_mask.sh
# Create the r-brms program, which is like Rscript but uses the image
apptainer_mask r-leafminer $SCRATCH/singularity/r-leafminer_4.3.1.sif Rscript
```

# Installing cmdstan

The first time you set up the image, you may need to install cmdstan. On
TACC, download it from github

``` bash
cmdstan_ver=2.33.1 # update this as needed
cmdstan_file="https://github.com/stan-dev/cmdstan/releases/download/v${cmdstan_ver}/cmdstan-${cmdstan_ver}.tar.gz"
mkdir $SCRATCH/cmdstan
cds cmdstan
wget $cmdstan_file 
# Extract the package
tar -xvzf cmdstan*tar.gz
cd cmdstan-${cmdstan_ver}
```

Create `make/local` with vim, nano, or another text editor, and enter
the following in it:

``` bash
# If you change anything in here, use make clean-all before rebuilding
# Enable threading
STAN_THREADS=true
# This example enables pedantic mode
STANCFLAGS+= --warn-pedantic
# Enable C++ compiler and linker optimization recommended by Stan developers.
STAN_CPP_OPTIMS=true
# Remove range checks from the model for faster runtime. 
STAN_NO_RANGE_CHECKS=true
# Adding other arbitrary C++ compiler flags
CXXFLAGS+= -funroll-loops -march=skylake-avx512
```

Now, open an `idev` terminal, open up the docker image in apptainer, and
build it.

``` bash
idev
ml tacc-apptainer
IMAGE=$SCRATCH/apptainer/r-brms_4.3.1.sif
# Open the image
apptainer shell -H /home $IMAGE

make build -j24 # build w/ 24 cores
```
