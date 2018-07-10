# KK July 03 2018
# Set up a bunch of runs to go in parallel

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
	parameter_combos = _parameter_search(parameter_ranges, n_runs, 'uniform')
 	
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


###---------------------------------------------------------------------------
###---------------------------------------------------------------------------

def _parameter_search(parameter_ranges, n_runs, search_type):
	# Return n_runs lists of parameter value
	# in lists with entries of forms { parameter_name : value }
	# (note that this is preferred format for scikit-learn searches)

	if search_type not in ['grid']:
		print "Search type " + search_type + " not yet implemented."
		print "Defaulting to a uniform random search."
		return _uniform_random_search(parameter_ranges, n_runs)
	
	elif search_type == 'grid':
		return _grid_search(parameter_ranges, n_runs)

	elif search_type == 'uniform':
		return _uniform_random_search(parameter_ranges, n_runs)


def _grid_search(parameter_ranges, n_runs):
	# Grid search of parameter spacee - TODO

	print("Grid search incomplete. Doing a uniform search instead.")
	return _uniform_random_search(parameter_ranges, n_runs)

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
