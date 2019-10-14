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


echo "Building Rescal-snow..."
cd ..
mkdir build
cd build
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release ..
make -j 16
cd ..

# Run snow-cone test to generate example data
echo "Attempting to run Rescal-snow using the snow_cone test case..."
echo "(This may take a few minutes)"
echo "(Output not displayed - use test-rescal.sh if necessary)"
cd scripts
./snow_cone.run

# Test the visualization example
echo "Attempting visualization example..."
cp ../test/viz_example_test.py .
python3 viz_example_test.py
rm viz_example_test.py
eog out/*image.png
echo "(Close image to continue)"
while true; do
	read -p "Does the output look good? [y/n]" yn
	case $yn in
		[Yy]* ) echo "Good! Continuing."; break;;
		[Nn]* ) exit;;
		* ) echo "Please input [y/n] to continue.";;
	esac
done


echo "----Testing complete.----"
