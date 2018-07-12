#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

in_dir="input_data"
out_dir="output_files"
ext=".data"
skip_files=5

echo This script will run analysis scripts on ReSCAL-snow Altitude files.
read -p "Enter directory to analyze: " dir
if [ -d "input_data/$dir" ]; then
    echo Option 1 - Run FFT analysis
    echo Option 2 - Run XCorrelation
    echo Option 3 - Run all analysis scripts
    read -p 'Which option do you want? (1-3): ' option
    if [ "$option" == "1" ]; then
	mkdir -p output_files/${dir}_OUT/fft
        read -p "Create png images and GIF animation? (y/n): " ans
        if [ "$ans" == "n" ]; then
      	    python fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/${dir}_OUT/fft/" 0 $ext $skip_files
        elif [ "$ans" == "y" ]; then
            fc=$(ls ${in_dir}/$dir/*.data | wc -l | xargs)
       	    read -p "Enter snapshot step size (1-$fc): " interval
            python fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/${dir}_OUT/fft/" $interval $ext $skip_files
        fi
    elif [ "$option" == "2" ]; then
    	mkdir -p ${out_dir}/${dir}_OUT/xcor
	echo Performing x-correlation analysis... 
	python xcorr-analysis.py "${in_dir}/$dir" "${dir}_OUT"
    elif [ "$option" == "3" ]; then
        mkdir -p ${out_dir}/${dir}_OUT/fft
	mkdir -p ${out_dir}/${dir}_OUT/xcor
	read -p "Create png images and GIF animation from fft analysis? (y/n): " ans
	if [ "$ans" == "n" ]; then
      	    python fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/${dir}_OUT/fft/" 0 $ext $skip_files
        elif [ "$ans" == "y" ]; then
	    fc=$(ls input_data/$dir/*.log | wc -l | xargs)
       	    read -p "Enter snapshot step size (1-$fc): " interval
            python fft2d_analysis.py "${}in_dir/$dir" "${out_dir}/${dir}_OUT/fft/" $interval $ext $skip_files
        fi
	echo Performing x-correlation analysis...
	python xcorr-analysis.py "{in_dir}/$dir" "${dir}_OUT"
    else
        echo You entered an invalid option.
    fi
else
    echo The directory: ${parent_path}/${in_dir}/$dir was not found. Please try a different directory.
fi
