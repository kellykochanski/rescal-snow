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

# Parameters held constant for all simulations
parameters = {
        'Model':  'SNO',
        'Output_directory': './out',
        'Csp_file': 'DUN.csp',
        'Csp_template': 'SNOWFALL(8)',
        'parfile': 'run.par',
	'Boundary':  'OPEN',
        'Time': 0.0,
        'H': 100,
        'L': 600,
        'D': 300,
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
        'Init_ncycl': 150,
	'rescallocation': '.' # this is where run script looks for rescal+genesis executables
}

# These for loops vary the parameters for snowfall rate and wind speed
for Lambda_I in [0.0, 0.0001, 0.001, 0.01, 0.1] :
	for Tau_min in [0, 10, 100, 200, 300, 1000] :
		this_run = Design_a_run()
		this_run.set_header("Baseline Lambda_I values at a specified tau_min")
		this_run.set_name("run")
		this_directory = "../runs/" + "tauMin" + str(Tau_min) + "_lambdaI" + str(Lambda_I)
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
		if not os.path.isdir(parameters['Output_directory']):
			os.mkdir(this_directory + '/out')

