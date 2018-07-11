#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
echo This script will run analysis scripts on ReSCAL-snow Altitude files.
echo Option 1 - Run FFT analysis
echo Option 2 - Run XCorrelation
echo Option 3 - Run all analysis scripts
read -p 'Which option do you want? (1-3): ' option
read -p "Enter directory to analyze: " dir
if [ -d "input_data/$dir" ]; then
    if [ "$option" == "1" ]; then
	mkdir -p output_files/${dir}_OUT/fft
        read -p "Create png images and GIF animation? (y/n): " ans
        if [ "$ans" == "n" ]; then
      	    python fft2d_analysis.py "input_data/$dir" "output_files/${dir}_OUT/fft/" 0
        elif [ "$ans" == "y" ]; then
            fc=$(ls input_data/$dir/*.log | wc -l | xargs)
       	    read -p "Enter snapshot step size (1-$fc): " interval
            python fft2d_analysis.py "input_data/$dir" "output_files/${dir}_OUT/fft/" $interval
        fi
    elif [ "$option" == "2" ]; then
    	mkdir -p output_files/${dir}_OUT/xcor
	echo Performing x-correlation analysis... 
	python xcorr-analysis.py "input_data/$dir" "${dir}_OUT"
    elif [ "$option" == "3" ]; then
        mkdir -p output_files/${dir}_OUT/fft
	mkdir -p output_files/${dir}_OUT/xcor
	read -p "Create png images and GIF animation from fft analysis? (y/n): " ans
	if [ "$ans" == "n" ]; then
      	    python fft2d_analysis.py "input_data/$dir" "output_files/${dir}_OUT/fft/" 0 '.data'
        elif [ "$ans" == "y" ]; then
	    fc=$(ls input_data/$dir/*.log | wc -l | xargs)
       	    read -p "Enter snapshot step size (1-$fc): " interval
            python fft2d_analysis.py "input_data/$dir" "output_files/${dir}_OUT/fft/" $interval '.data'
        fi
	echo Performing x-correlation analysis...
	python xcorr-analysis.py "input_data/$dir" "${dir}_OUT"
    else
        echo You entered an invalid option.
    fi
else
    echo The directory: $parent_path/$dir was not found. Please try a different directory.
fi
