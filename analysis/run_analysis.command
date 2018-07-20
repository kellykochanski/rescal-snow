#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

in_dir="input_data" #The main input directory relative to this script
out_dir="output_data" #The main output directory
fft_dir=""
xcor_dir=""
pref="ALTI"
skip_files=1

echo This script will run analysis scripts on ReSCAL-snow Altitude files.
read -p "Enter directory path to analyze: " dir
if [ -d "${in_dir}/$dir" ]; then
    echo Option 1 - Run FFT analysis
    echo Option 2 - Run XCorrelation
    echo Option 3 - Run all analysis scripts
    read -p 'Which option do you want? (1-3): ' option
    if [ "$option" == "1" ]; then
        read -p "Create png images and GIF animation? (y/n): " ans
        if [ "$ans" == "n" ]; then
      	    python ${fft_dir}fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/fft/" "$(basename $dir)_analysis" 0 $pref $skip_files
        elif [ "$ans" == "y" ]; then
            fc=$(ls ${in_dir}$dir/*.data | wc -l | xargs)
       	    read -p "Enter snapshot step size (1-$fc): " interval
            python ${fft_dir}/fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/fft/" "$(basename $dir)_analysis" $interval $pref $skip_files
        fi
    elif [ "$option" == "2" ]; then
	echo Performing x-correlation analysis... 
	python xcor-slices.py "${in_dir}/$dir" "${out_dir}/xcorr/" "$(basename $dir)"
    elif [ "$option" == "3" ]; then
	read -p "Create png images and GIF animation from fft analysis? (y/n): " ans
	if [ "$ans" == "n" ]; then
      	    python fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/fft/" "$(basename $dir)_analysis" 0 $pref $skip_files
        elif [ "$ans" == "y" ]; then
	    fc=$(ls input_data/$dir/*.log | wc -l | xargs)
       	    read -p "Enter snapshot step size (1-$fc): " interval
            python fft2d_analysis.py "${in_dir}/$dir" "${out_dir}/fft/" "$(basename $dir)_analysis" $interval $pref $skip_files
        fi
	echo Performing x-correlation analysis...
	python xcor-slices.py "{in_dir}/$dir" "${out_dir}/xcorr/" "$(basename $dir)"
    else
        echo You entered an invalid option.
    fi
else
    echo The directory: ${parent_path}/${in_dir}/$dir was not found. Please try a different directory.
fi
