"""
Rescal-snow: a cellular automaton model of self-organized snow
Copyright (C) 2019 Kelly Kochanski

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.
This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

in_dir="input_data" #The main input directory relative to this script
out_dir="output_data" #The main output directory
fft_dir="downie4_snowcal/rescal-snow/analysis/"
xcor_dir="robeson3_snowcal/rescal-snow/analysis/"
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
	read -p "Enter output name (leave blank if multiple directories): " outname
        if [[ -z "$outname" ]]; then
  	    python ${fft_dir}fft2d_analysis.py "${in_dir}/${dir}/" "${out_dir}/fft/${dir}/" "" 0
	else
	    read -p "Create png images and GIF animation? (y/n): " ans
      	    if [ "$ans" == "n" ]; then
	        python ${fft_dir}fft2d_analysis.py "${in_dir}/${dir}/" "${out_dir}/fft/" $outname 0
            elif [ "$ans" == "y" ]; then
                fc=$(ls ${in_dir}/$dir/*.data | wc -l | xargs)
       	        read -p "Enter snapshot step size (1-$fc): " interval
                python ${fft_dir}fft2d_analysis.py "${in_dir}/${dir}/" "${out_dir}/fft/" $outname $interval
            fi
	fi
    elif [ "$option" == "2" ]; then
	echo Performing x-correlation analysis... 
	python ${xcor_dir}xcor-slices.py "${in_dir}/$dir" "${out_dir}/xcorr/" "$(basename $dir)"
    elif [ "$option" == "3" ]; then
	read -p "Enter output name (leave blank if multiple directories): " outname
        if [ "$outname" != "" ]; then
  	    read -p "Create png images and GIF animation? (y/n): " ans
      	    if [ "$ans" == "n"]; then
	        python ${fft_dir}fft2d_analysis.py "${in_dir}/$dir/" "${out_dir}/fft/" $outname 0
            elif [ "$ans" == "y" ]; then
                fc=$(ls ${in_dir}$dir/*.data | wc -l | xargs)
       	        read -p "Enter snapshot step size (1-$fc): " interval
                python ${fft_dir}fft2d_analysis.py "${in_dir}/$dir/" "${out_dir}/fft/" $outname $interval
            fi
	else
	    python ${fft_dir}fft2d_analysis.py "${in_dir}/$dir/" "${out_dir}/fft/" "" 0
	fi
	echo Performing x-correlation analysis...
	python ${xcor_dir}xcor-slices.py "{in_dir}/$dir" "${out_dir}/xcorr/" "$(basename $dir)"
    else
        echo You entered an invalid option.
    fi
else
    echo The directory: ${parent_path}/${in_dir}/$dir was not found. Please try a different directory.
fi
