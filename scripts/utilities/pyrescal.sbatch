#!/bin/bash

# example of using sbatch to run multiple instances of example_pyrescal 

# ensure that RESCAL_SNOW_ROOT/data_runs/ exists
# make it if it doesn't
if [ ! -d "$RESCAL_SNOW_ROOT/data_runs/" ]; then
    mkdir "$RESCAL_SNOW_ROOT/data_runs/"
fi


#SBATCH --job-name r0
#SBATCH --output=log.txt
#SBATCH --qos normal
#SBATCH --mail-user=user@domain.com
#SBATCH --mail-type=ALL
#SBATCH --time=00:05:00
#SBATCH --ntasks 8
srun --wait=0 --cpus-per-task=2 --ntasks=8 example_pyrescal.py

