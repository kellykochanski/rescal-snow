#!/bin/bash

# Kelly Kochanski Sept 27 2019
# Test script for Rescal-snow
# This script:
# - builds Rescal-snow
# - Runs one example (scripts/snow_cone.run)
# - Imports all the utilities scripts
# - Uses the heightmap utility to visualize some files
# - Uses the parameter space exploration utility
#    and rescal_utilities to set up some potential runs
# Success at each stage is confirmed by the user.


echo "Attempting to build Rescal-snow..."
cd ..
mkdir build
cd build
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release ..
make -j 16
cd ..
while true; do
	read -p "Did Rescal-snow build successfully? [y/n]" yn
	case $yn in
		[Yy]* ) echo "Good! Continuing."; break;;
		[Nn]* ) exit;;
		* ) echo "Please input [y/n] to continue.";;
	esac
done

# Test that Rescal-snow runs
echo "Attempting to run Rescal-snow using the snow_cone test case..."
echo "(This may take a few minutes)"
cd scripts
./snow_cone.run
eog *.png
while true; do
	read -p "Does the output look good? [y/n]" yn
	case $yn in
		[Yy]* ) echo "Good! Continuing."; break;;
		[Nn]* ) exit;;
		* ) echo "Please input [y/n] to continue.";;
	esac
done


# Test the visualization example
echo "Attempting visualization example..."
cp ../test/viz_example_test.py .
python3 viz_example_test.py
rm viz_example_test.py
eog out/*image.png
while true; do
	read -p "Does the output look good? [y/n]" yn
	case $yn in
		[Yy]* ) echo "Good! Continuing."; break;;
		[Nn]* ) exit;;
		* ) echo "Please input [y/n] to continue.";;
	esac
done

# Test rescal_utilities for setting up parallel runs
echo "Attempting rescal_utilities example..."
cd utilities
python3 param_space_exploration_example.py
cd ../..
echo "Printing parallel directories:"
ls test_parallel_runs
echo "Printing contents of one directory:"
ls test_parallel_runs/tauMin0_lambdaI0.001
rm -rf test_parallel_runs
while true; do
	read -p "Did Rescal-snow set up 10 parallel directories, each with rescal, genesis, real_data, run.run and run.par? [y/n]" yn
	case $yn in
		[Yy]* ) echo "Good! Continuing."; break;;
		[Nn]* ) exit;;
		* ) echo "Please input [y/n] to continue.";;
	esac
done

# Remove results of test
cd scripts
./clean -f
rm -rf out
cd ..

echo "----Testing complete.----"
