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

"""
# KK July 03 2018

This is a shell script designed to start organizing ReSCAL to perform a parameter space exploration;
it could be adapted to perform any set of runs with variable parameters

The "explore_parameter_space" function creates an array of subdirectories, each containing:
   - a run script (bash), generally called run.run
   - a parameter file (txt), generally caled run.par
 Note that the multi-directory structure is useful for several reasons:
   1. The output for each run stays in its own subdirectory, attached to the appropriate input files
   2. If multiple processors are available, they may simultaneously run separate instances of rescal-snow, one per directory,
       without any communication, race conditions, or other problems

"""

import rescal_utilities
import numpy as np
import datetime
import os

def explore_parameter_space(output_root, executable_location, parameter_ranges, fixed_parameters, n_runs):
	# Take in a dictionary of parameter ranges to explore,
	#   e.g. { 'Lambda_S' : [0.001, 0.2],
	#          'Lambda_F' : [3, 10] }

	# Setup n_runs runs to explore those, under directory output_ro
	# Rescal and genesis executables are located in directory executable_location

	if ( os.path.isdir(output_root) == False):
		os.mkdir(output_root)

	# Get some lists of parameter combos to explore
	parameter_combos = _uniform_random_search(parameter_ranges, n_runs)
 	
	# Write parameter files for each set of parameters
	now = datetime.datetime.now()
	for run in range(n_runs):
		run_output_dir 		= output_root + "/" + \
				str(now.year) + str(now.month) + str(now.day) + "_" + str(run)

		# Avoid overwriting previous results or using a messy directory
		while os.path.isdir(run_output_dir):
			run_output_dir += "A"

		this_run = rescal_utilities.Design_a_run()
		this_run.set_parameters(parameter_combos[run])
		this_run.set_parameters({'rescallocation' : executable_location,
			'parfile': 'run.par'})
		this_run.set_parameters(fixed_parameters)
		this_run.set_directory(run_output_dir)

		os.mkdir(run_output_dir)
		this_run.write()

def _uniform_random_search(parameter_ranges, n_runs):
	# For each run, each parameter takes a value drawn from a uniform distribution
	
	parameter_combinations = []

	for run in range(n_runs):
		these_parameters = {}
		for parameter in parameter_ranges.keys():
			valuerange	= parameter_ranges[parameter]
			these_parameters[parameter] = np.random.uniform(valuerange[0], valuerange[1])
	
		parameter_combinations.append(these_parameters)
	
	return parameter_combinations
