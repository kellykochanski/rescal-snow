#!/bin/bash

# Set up this run in the correct directory
#Configured for LL CPU not RC

output_root="../test_gaussian_parallel"
echo "Looking for jobs to run under ${output_root}"

output_dirs=(${output_root}/*)

output_dir="${output_dirs[${PMI_RANK}]}"
cd ${output_dir}
echo "Process ${PMI_RANK} running in  ${output_dir}"

# Need to be able to read and execute this directory
chmod u+rwx *

./run.run
