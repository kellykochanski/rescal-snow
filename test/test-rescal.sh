#!/bin/bash

# Kelly Kochanski Sept 27 2019
# Test script for Rescal-snow
# This script:
# - builds Rescal-snow
# - Runs one example (scripts/snow_cone.run)
# - Compares one output png to test data

set -o errexit  #Exit if any command fails
set -o pipefail #Exit if any pipe fails
set -o nounset  #Exit if an unset variable is de-referenced

# Allow the script to be run from project root, e.g., test/test-utilities.sh or
# from the test directory.
test_dir="$(cd "$(dirname "$0")" && pwd)"

# Find number of available cores
if command -v ncpus > /dev/null 2>&1 ; then
    NCPU=$(nproc --all)
elif [[ "$OSTYPE" == "[Dd]arwin"* ]]; then
    NCPU=$(sysctl -n hw.logicalcpu)
else
    NCPU=8
fi

echo "Attempting to build Rescal-snow..."
cd "${test_dir%/}/.." || exit 1
rm -rf build || true # Don't error if not present
mkdir build
(cd build
cmake -Wdev -DCMAKE_BUILD_TYPE=Release ..
make -j ${NCPU})


# Test that Rescal-snow runs
echo "Attempting to run Rescal-snow using the snow_cone test case..."
echo "(This may take a few minutes)"
cd "${test_dir%/}/../scripts"
./snow_cone.run
diff SNO00003_t0.png ../test/SNO00003_t0.png

# Test the visualization example
# echo "Attempting visualization example..."
# cp ../test/viz_example_test.py .
# python3 viz_example_test.py
# rm viz_example_test.py
# diff out/ALTI00030_t0_image.png ../test/ALTI00030_t0_image.png

# Test rescal_utilities for setting up parallel runs
echo "Attempting rescal_utilities example..."
cd "${test_dir%/}/../scripts/utilities"
python3 param_space_exploration_example.py

cd "${test_dir%/}/../"
echo "Printing parallel directories:"
ls test_parallel_runs
echo "Printing directory contents:"
ls test_parallel_runs/tauMin0_lambdaI0.001
ls test_parallel_runs/tauMin0_lambdaI0.01
ls test_parallel_runs/tauMin1000_lambdaI0.001
ls test_parallel_runs/tauMin1000_lambdaI0.01
ls test_parallel_runs/tauMin100_lambdaI0.001
ls test_parallel_runs/tauMin100_lambdaI0.01
ls test_parallel_runs/tauMin200_lambdaI0.001
ls test_parallel_runs/tauMin200_lambdaI0.01
ls test_parallel_runs/tauMin300_lambdaI0.001
ls test_parallel_runs/tauMin300_lambdaI0.01
rm -rf test_parallel_runs

# Remove results of test
cd "${test_dir%/}/../scripts"
./clean -f
rm -rf out
cd ..

echo "----Testing complete.----"
