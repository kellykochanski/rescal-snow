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

# Catch common scripting errors with stricter checking
set -o errexit
set -o pipefail
set -o nounset

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

echo "Building Rescal-snow..."

# Use a subshell to automatically `cd` back to the right place
(cd "${test_dir%/}/.." || exit 1
mkdir build || true # Don't throw an error if the build exists
cd build
# ensure that we're not ignoring errors caused by recent changes in CMake build files
rm CMakeCache.txt
# Don't silence developer warnings during testing, fix any that may be triggered
cmake -Wdev -DCMAKE_BUILD_TYPE=Release ..
make clean # Make sure that we're actually rebuilding everything
make -j ${NCPU}
)

# Run snow-cone test to generate example data
echo "Attempting to run Rescal-snow using the snow_cone test case..."
echo "(This may take a few minutes)"
echo "(Output not displayed - use test-rescal.sh if necessary)"
cd "${test_dir%/}/../scripts" || exit 1
./snow_cone.run

# Test the visualization example
echo "Attempting visualization example..."
cp ../test/viz_example_test.py .
python3 viz_example_test.py
rm viz_example_test.py

# # Create a gif of the PNGs if we have imagemagick
# if command -v convert >/dev/null 2>&1 ; then
#     "${test_dir%/}/../scripts/visualization/png_to_gif.sh" out/ALTI_cone.gif out/*image.png
# fi

# Portably look for image viewing program, using system default application when possible
if command -v xdg-open > /dev/null 2>&1 ; then
    VIEWER=xdg-open
elif command -v open > /dev/null 2>&1 ; then
    VIEWER=open
elif command -v eog > /dev/null 2>&1 ; then
    VIEWER=eog
else
    no_viewer () {
        printf "No known image viewer found, please manually open and inspect %s.\n" "$*"
    }
    VIEWER=no_viewer
fi
$VIEWER out/*image.png
if [[ -f out/ALTI_cone.gif ]]; then
    $VIEWER out/ALTI_cone.gif
fi
echo "(Close image to continue)"
while true; do
    read -rp "Does the output look good? [y/n]" yn
    case $yn in
	[Yy]* ) echo "Good! Continuing."; break;;
	[Nn]* ) echo "Remember to manually run the clean script and remove the out directory"; exit 1 ;;
	* ) echo "Please input [y/n] to continue.";;
    esac
done

echo "----Testing complete.----"
