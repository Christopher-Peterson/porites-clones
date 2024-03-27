# Create binary wrappers for apptainer images
  # It will create one image that runs your command ($CMD) and another that launches the shell
# Before using, you should: 
  # mkdir $WORK/apptainer/bin -p
  # mkdir $SCRATCH/apptainer/bin -p
  # edit your .bashrc file to include $SCRATCH/apptainer/bin in your $PATH
  # apptainer pull your image to $WORK/apptainer and copy it to $SCRATCH/apptainer
# This should be saved in $WORK/singularity
# Source this in the install slurm script

# Be sure to add $WORK/apptainer/bin to $PATH-

function apptainer_mask {
  # Args: bin_name, image file, command = bin_name, home directory = /home (inside the image)
  local BIN_NAME=$1
  local SIF=$2 # Should contain whole path
  local CMD=${3:-$BIN_NAME}
  local HDIR=${4:-/home}
  local BIN=$WORK/apptainer/bin/$BIN_NAME
  local SHELL_BIN=${BIN}-shell
  local SCRATCH_BIN=$SCRATCH/apptainer/bin/$BIN_NAME
  local SCRATCH_SHELL=$SCRATCH/apptainer/bin/${BIN_NAME}-shell
  local WORK_SIF=${SIF/$SCRATCH/$WORK}
  # Setup the base of the file
  echo "#!/bin/bash
    module load tacc-apptainer #tacc-singularity 
    HDIR=$HDIR
    IMAGE=$SIF
    CMD=$CMD" > $BIN
  # It is probably best to NOT have this in here, and instead add a function that checks the image instead
  echo 'touch $IMAGE' >> $BIN
  # Add a check to see if the $SCRATCH version doesn't exist; if so, copy it from work
#   echo "[ ! -f $SIF ] && cp $WORK_SIF $SIF" >> $BIN # let's not do this for now
  cp $BIN $SHELL_BIN
  # Add the final line using single quotes, to avoid escaping
  echo 'apptainer exec -H $HDIR $IMAGE $CMD "$@" ' >> $BIN
  echo 'apptainer shell -H $HDIR $IMAGE' >> $SHELL_BIN
  # Make it executable
  chmod +x $BIN $SHELL_BIN
  cp $BIN $SCRATCH_BIN
  cp $SHELL_BIN $SCRATCH_SHELL
}