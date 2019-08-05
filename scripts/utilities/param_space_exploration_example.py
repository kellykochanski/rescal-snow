from rescal_utilities import *
import os
import shutil

# Kelly Kochanski and Adam Rubin, 2018

# This is an example file showing how to run a parameter space exploration
# The exploration fixes most simluations parameters (the `parameters' dictionary')
#  then explores five snowfall rates, controlled by parameter Lambda_I,
#  six wind speeds, controlled by parameter Tau_min,
#  and all combinations thereof.

# The script creates new directories and manages the locations of input, output and executable files
#  so that all 30 simulations can be run in parallel.


executable_location = ".." # location of the compiled rescal executable
experiment_name = "test_parallel_runs" # header directory which all generated input goes into
experiment_directory = os.path.join("../../", experiment_name)
if not os.path.isdir(experiment_directory):
	os.mkdir(experiment_directory)

# Parameters held constant for all simulations
parameters = {
        'Model':  'SNO',
        'Output_directory': './out',
        'Csp_file': 'DUN.csp',
        'Csp_template': 'SNOWFALL(4)',
        'parfile': 'run.par',
	'Boundary':  'OPEN',
        'Time': 0.0,
        'H': 50,
        'L': 400,
        'D': 100,
        'Centering_delay': 0,
        'Phys_prop_file': 'real_data/sealevel_snow.prop',
        'Qsat_file': 'real_data/PDF.data',
        'Lambda_A': 1,
        'Lambda_E': 1,
        'Lambda_T': 1.5,
        'Lambda_C': 0.5,
        'Lambda_G': 100000,
        'Lambda_D': 0.01,
        'Lambda_S': 0,
        'Lambda_F': 1,
        'Coef_A': 0.1,
        'Coef_B': 10,
        'Coef_C': 10,
        'Prob_link_ET': 0.5,
        'Prob_link_TT': 1.0,
        'High_mobility': 1,
        'Ava_mode': 'TRANS',
        'Ava_angle': 38,
        'Ava_h_lim': 1,
        'Lgca_delay': 1,
        'Lgca_speedup': 1000,
	'rescallocation': '.' # this is where run script looks for rescal+genesis executables
}

# These for loops vary the parameters for snowfall rate and wind speed
for Lambda_I in [0.001, 0.01] :
	for Tau_min in [0, 100, 200, 300, 1000] :
		this_run = Design_a_run()
		this_run.set_header("Baseline Lambda_I values at a specified tau_min")
		this_run.set_name("run")
		# Where should the input for this single run go?
		this_directory = os.path.join(experiment_directory, "tauMin" + str(Tau_min) + "_lambdaI" + str(Lambda_I))
		if not os.path.isdir(this_directory):
			os.mkdir(this_directory)
		this_run.set_directory(this_directory)
		parameters['Lambda_I'] = Lambda_I
		parameters['Tau_min'] = Tau_min
		parameters['Tau_max'] = Tau_min + 1000
		this_run.set_parameters(parameters)
		this_run.write()
		# copy executables for this run into this run directory
		shutil.copyfile(executable_location + "/rescal", this_directory+'/rescal')
		shutil.copyfile(executable_location + "/genesis", this_directory+'/genesis')
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

