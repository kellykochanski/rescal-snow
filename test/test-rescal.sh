#!/bin/bash

# Kelly Kochanski Sept 27 2019
# Test script for Rescal-snow
# This script:
# - builds Rescal-snow
# - Runs one example (scripts/snow_cone.run)
# - Compares one output png to test data

set -e          #Exit if any command fails
set -o pipefail #Exit if any pipe fails

echo "Attempting to build Rescal-snow..."
cd ..
rm -rf build
mkdir build
cd build
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release ..
make -j 16
cd ..

# Test that Rescal-snow runs
echo "Attempting to run Rescal-snow using the snow_cone test case..."
echo "(This may take a few minutes)"
cd scripts
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
cd utilities
python3 param_space_exploration_example.py
cd ../..
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
cd scripts
./clean -f
rm -rf out
cd ..

echo "----Testing complete.----"
