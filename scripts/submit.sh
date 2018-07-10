#!/bin/bash

# Set up this run in the correct directory


output_root="../new_runs3"
echo ${output_root}

output_dirs=(${output_root}/*)

output_dir="${output_dirs[${PMI_RANK}]}"
cd ${output_dir}
echo "Process ${PMI_RANK} in  ${output_dir}"

# Need to be able to read and execute this directory
chmod u+rwx *

./run.run
