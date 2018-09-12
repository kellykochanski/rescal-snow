import setup_runs

output_root = "../../new_runs3"

# Path from output_root/run_directory to the location of 'rescal' and 'genesis' executables
executable_location = "../../src"

parameter_ranges = {'Lambda_S' : [0,0.5],
		'Tau_min' : [0,40],
		'H' : [25, 200],
		'Lambda_F' : [0, 50],
		'Lambda_T' : [0, 30],
		'Lambda_D' : [0, 0.1],
		'Lambda_I' : [0, 0.5],
		'Ava_angle' : [10, 90],
		'Coef_A'   : [0, 3],
		'Csp_params' : [0, 25]}

fixed_parameters = {'Environment' : 'snow',
		'Csp_template' : 'LAYER(30)',
		'L' : 1500,
		'D' : 400,
		'Output_directory' : "."}

n_runs = 34

setup_runs.explore_parameter_space(output_root, executable_location, parameter_ranges, fixed_parameters, n_runs)
