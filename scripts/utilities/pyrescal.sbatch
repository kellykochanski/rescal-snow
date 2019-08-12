#!/bin/bash
#SBATCH --job-name r0
#SBATCH --output=log.txt
#SBATCH --qos normal
#SBATCH --mail-user=defazio1@llnl.gov
#SBATCH --mail-type=ALL
#SBATCH --time=00:05:00
#SBATCH --ntasks 8
srun --wait=0 --cpus-per-task=2 --ntasks=8 example_pyrescal

