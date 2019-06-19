#!/bin/bash
#SBATCH --qos normal
#SBATCH --mail-user=yourusername@yourserver.com
#SBATCH --mail-type=ALL
#SBATCH --time=00:05:00
#SBATCH --ntasks 10


# This is an example script for submitting runs to a computing cluster
# Current options:
#   - mail-user : set to your email to receive notifications about run start/stop/crash
#   - time      : time for run (currently: 1 hour, 01:00:00)
#   - ntasks    : number of processors requested (expect 1 processor per ReSCAL instance)
# The syntax of these commands will almost certainly need to be modified for your computing resources
#  along with the number of processors per node, the run time for a given simulation, etc.

# Srun launches N instances of the "submit.sh" script,
#  which in turn will direct each of N processors to a separate directory
#  where each processor will run an instance of ReSCAL.
# Use the "--wait=0" (or equivalent for your system) flag 
#  to prevent one finished/failed run from stopping others.

echo "run_parallel_rescal.sbatch : starting to run ReSCAL in parallel"

srun --nodes=1 --cpus-per-task=1 --ntasks=10 --wait=0 submit.sh
