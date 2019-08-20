from rescal_utilities import *
import os
import shutil
import numpy as np

# Kelly Kochanski and Adam Rubin, 2018

# This is an example file showing how to run a parameter space exploration
# The exploration fixes most simluations parameters (the `parameters' dictionary')
#  then explores five snowfall rates, controlled by parameter Lambda_I,
#  six wind speeds, controlled by parameter Tau_min,
#  and all combinations thereof.

# The script creates new directories and manages the locations of input, output and executable files
#  so that all 30 simulations can be run in parallel.


executable_location = "../../build" # location of the compiled rescal executable
experiment_name = "test_gaussian_parallel" # header directory which all generated input goes into
experiment_directory = os.path.join("../../", experiment_name)
if not os.path.isdir(experiment_directory):
	os.mkdir(experiment_directory)

# Parameters held constant for all simulations
parameters = {
        'Model':  'SNO',
        'Output_directory': './out',
        'Csp_file': 'DUN.csp',
        'Csp_template': 'INPUT_ELEVATION(initial_elevation.data, 1)',
        'parfile': 'run.par',
	'Boundary':  'OPEN', #???
        'Time': 0.0, #???
        'H': 100, #100
        'L': 600, #600
        'D': 150, #400
        'Centering_delay': 0,
        'Phys_prop_file': 'real_data/sealevel_snow.prop',
        'Qsat_file': 'real_data/PDF.data',
        'Lambda_A': 1.0,
        'Lambda_E': 1.0,
        'Lambda_T': 1.5,
        'Lambda_C': 0.5, 
        'Lambda_G': 100000.0,
        'Lambda_D': 0.01,
        'Lambda_S': 0.0,
        'Lambda_F': 1.0,
        'Lambda_I': 0.0,
        'Tau_min': 100.0,
        'Tau_max': 1100.0,
        'Coef_A': 0.1,
        'Coef_B': 10.0,
        'Coef_C': 10.0, 
        'Prob_link_ET': 0.5,
        'Prob_link_TT': 1.0,
        'High_mobility': 1.0,
        'Ava_mode': 'TRANS',
        'Ava_angle': 35.0,
        'Ava_h_lim': 1.0,
        'Lgca_delay': 1.0,
        'Lgca_speedup': 1000.0,
	'rescallocation': '.' # this is where run script looks for rescal+genesis executables
}

num_trials = 360
#SEED = 94550
#np.random.seed(SEED)

curr_dirs = {}

while (True):
	if (len(curr_dirs) == num_trials):
		break

	Pile_height = np.random.random_integers(2, parameters['H'] // 3)
	Pile_width = np.random.random_integers(Pile_height, Pile_height * 5)
	
	curr_dirs[(Pile_height, Pile_width)] = curr_dirs.get((Pile_height, Pile_width), 0) + 1
	if (curr_dirs[(Pile_height, Pile_width)] > 1): 
		continue
	
	this_run = Design_a_run()
	this_run.set_parameters({'output interval' : '20t0', 'png interval' : '20t0'})
	this_run.set_header("Pile Height and Pile Width for Gaussian Terrain. ")
	this_run.set_name("run")
	# Where should the input for this single run go
	this_directory = os.path.join(experiment_directory, "Pile_height" + str(Pile_height) + "_Pile_width" + str(Pile_width))
	if not os.path.isdir(this_directory):
		os.mkdir(this_directory)
	this_run.set_directory(this_directory)
	#parameters['Lambda_I'] = Lambda_I
	#parameters['Tau_min'] = Tau_min
	#parameters['Tau_max'] = Tau_min + 1000
	this_run.set_parameters(parameters)
	this_run.write("gaussian", ["initial_elevation.data", parameters['L'], parameters['H'], parameters['D'], Pile_width, Pile_height])
	# copy executables for this run into this run directory
	shutil.copyfile(executable_location + "/rescal", this_directory+'/rescal')
	shutil.copyfile(executable_location + "/genesis", this_directory+'/genesis')
	shutil.copyfile(".." + "/gaussian.py", this_directory+'/gaussian.py')
	if not os.path.isdir(this_directory+'/real_data'):
		os.mkdir(this_directory+'/real_data')
	shutil.copyfile('../real_data/sealevel_snow.prop', this_directory+'/real_data/sealevel_snow.prop')
	shutil.copyfile('../real_data/PDF.data', this_directory+'/real_data/PDF.data')
	# while we're at it, make sure there's a place for output
	output_dir = parameters['Output_directory']
	if output_dir[0:1] == "./":
		output_dir = output_dir[2:]
		output_dir = os.path.join(experiment_directory, output_dir)
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)


