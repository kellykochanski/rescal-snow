# KK Jun 25 20180
# Utilities to set up run scripts quickly and easily; these tools are designed to aid parameter space
# explorations, sensitivity analyses, and large batches of runs.
# Parameters : Class to hold, update, change, or write all the parameters that ReSCAL needs to run
# Run_Script : Class to hold, update, change, or write a ReSCAL run script with appropriate flags
# This includes almost all options and inputs, except those in the "real_data" file.

class Design_a_run():
	# Design a ReSCAL run
	# This interfaces with both the Parameters and Run_Script utilities
	# in order to produce run script and parameter files simultaneously
	# using a standardized input format

	def __init__(self):
		self.name 	= 'run'
		self.directory 	= '.' 
		self.parameters = Parameters()
		self.run_script = Run_Script()
		self.parameters.set_header("Test from write_tests.py")
         
	def set_parameters(self, param_dict):
        	# Sets all parameters, given as {name : value} pairs
        	# Allows user to mix things that go in parameter file
        	# with things that go in the run script, etc
        	# TODO make it possible to do compiler flags too?
                 
		for name in param_dict.keys():
			if self.__is_a_parameter(name): 
				self.parameters.set({name : param_dict[name]})
			elif self.__is_a_run_script_option(name):
				self.run_script.set({name : param_dict[name]})
			else:   
				print "Warning : skipped nonexistent parameter " + name
         
	def list_all(self):
		# Lists all available options/parameters
		return self.run_script.list_all(), self.parameters.list_all()
         
	def get(self, param):
		if param in self.parameters.list_all():
			return self.parameters.get(param)
		elif param in self.run_script.list_all():
			return self.run_script.get(param)
		else:
			print "Error: requested nonexistent parameter " + param
			return False

	def set_header(self, description):
        	self.parameters.set_header("Test from write_tests.py: " + description)
         
	def set_name(self, name):
		 self.name = name

	def set_directory(self, directory):
		self.directory = directory
         
	def __is_a_parameter(self, name):
        	# Checks whether the parameters dictionary contains 'name'
        	return ((name == 'Environment') or (name in self.parameters.list_all()))
         
	def __is_a_run_script_option(self, name):
        	# Checks whether 'name' is an option in the run script
        	return (name in self.run_script.list_all())
         
	def write(self, input_elevation_type=None, args=None):
        	self.run_script.write(self.directory + "/" + self.name + ".run", input_elevation_type, args)
        	self.parameters.write(self.directory + "/" + self.name + ".par")

##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

class Parameters():
	# Holds all the parameters that ReSCAL needs to run
	def __init__(self):
		self.parameters 		= {}
		self.parameter_descriptions 	= {}
		self.__special_parameter_descriptions = {}
		self.__set_special =			{}
		self.__default_header()
		self.__default_parameters()

	def new_parameter(self, name, description, default_value):
		# New parameter
		# Name should be a string (the parameter name expected by ReSCAL)
		# Description should be a string describing its use or purpose etc
		# Default value is the default value, generally a string or float
		self.parameters[name] 		  = default_value
		self.parameter_descriptions[name] = description

	def __special_parameter(self, name, description, setting_function):
		# Special parameters are not stored.
		# When set, they modify one or more other parameters in a way described by the setting_function()
		self.__special_parameter_descriptions[name]	= description
		self.__set_special[name]			= setting_function

	def __default_header(self):
		# Provides a basic description of this parameter set in form of a string
		self.header_start = "Parameter file for ReSCAL written from Parameters utility"
		self.header_body  = "Default parameters : sahara sand"

	def __default_parameters(self):
		# Default values for parameters
		# Initialized here for sahara sand
		# Basic model parameters and geometry of model run:
		self.new_parameter('Model', 'Defines a set of transitions for grains in the cellular automaton e.g. DUN or SNO', 'SNO')
		self.new_parameter('Output_directory', 'Directory where output files go', './out')

		self.new_parameter('Csp_file', 	'Genesis template?', 	'DUN.csp')
		self.new_parameter('Csp_template', 'Template containing initial shape for run', 'RCONE')
		self.new_parameter('Csp_params', 'Parameters for size, depth etc of template shape', [])
		self.new_parameter('Boundary', 	'Type of boundary conditions, e.g. OPEN, PERIODIC, REINJECTION', 'PERIODIC')
		self.new_parameter('Time', 	'Initial time', 	0.0)
		self.new_parameter('H',		'Height of model domain (cells)', 	80)
		self.new_parameter('L',		'Length, longitudinal, of model domain (cells)', 	300)
		self.new_parameter('D', 	'Depth, transverse to flow, of model domain (cells)', 	200)
		self.new_parameter('Centering_delay', 'Automatic re-centering of model domain on dune', 0)

		# Parameters describing the physical behavior of the grains in the model
		self.new_parameter('Phys_prop_file', 'File containing important physical properties e.g. gravity', 'real_data/desert_earth.prop')
		self.new_parameter('Qsat_file',	'File containing relationship between saturated and threshold flux', 'real_data/PDF.data')
		self.new_parameter('Lambda_E',  'Erosion rate (when threshold exceeded)', 1)
		self.new_parameter('Lambda_T', 	'Transition rate for mobile grains', 1.5)
		self.new_parameter('Lambda_C',	'Deposition rate for mobile grains on surface', 0.5)
		self.new_parameter('Lambda_G',	'Gravity - rate of fall of not-in-transport grains in air', 1000)
		self.new_parameter('Lambda_D', 	'Transition rate for diffusion', 0.01)
		self.new_parameter('Lambda_S',  'Transition rate for sintering/cohesion', 0)
		self.new_parameter('Lambda_F',  'Ratio of sintered:unsintered erosion thresholds', 3)
		self.new_parameter('Coef_A', 	'Ratio of vertical:horizontal transport of mobile grains', 0.1)
		self.new_parameter('Coef_B', 	'Ratio of deposition against an obstacle : deposition', 1)
		self.new_parameter('Coef_C',	'Ratio of deposition behind an obstacle : deposition', 3)
		self.new_parameter('Prob_link_ET', 'Probability of the transition links', 0.5)
		self.new_parameter('Prob_link_TT', 'Probability of the transition links', 1)
		self.new_parameter('High_mobility', 'Higher mobility of grains', 1)
		self.new_parameter('Lambda_I', 	'Injection rate', 0)
		self.new_parameter('Lambda_A', 	'??', 1)

		# Parameters describing avalanching
		self.new_parameter('Ava_mode', 	'Mode of avalanching', 'TRANS')
		self.new_parameter('Ava_angle', 'Angle of repose, steeper slopes avalanche (degrees)', 35)
		self.new_parameter('Ava_h_lim',	'Height limit in avalanches (cells)', 1)

		# Parameters for flow stabilization - how often the lgca runs and how strong is the flow
		self.new_parameter('Lgca_delay', 'Delay between flow cycles', 1)
		self.new_parameter('Lgca_speedup', 'Speedup of the flow stabilization', 1000)
		self.new_parameter('Tau_min', 	'Shear stress threshold for grain erosion, controls flow strength', 0)
		self.new_parameter('Tau_max', 	'Max shear stress, used in grain erosion/flow speed relation', 1000)

		# Special parameters - set by same interface as others, 
		# But are actually tools which modify one or more other parameters
		self.__special_parameter('Environment', 'Sets all defaults for an environment, e.g. sand or snow', self.__environment)
		self.__special_parameter('Csp_params',  'Sets parameters for a Csp_template', self.__update_template)

	def __environment(self, keyword):
		# Function to set special 'Environment' parameter
		#Update defaults fora  different enviroment, e.g. snow
		# (Could be added for subaqueous dunes, mars environment, etc)
		if keyword == "snow":
			self.header_body = "Default values - Niwot Ridge snow"
			self.set({'Model' 		: 'SNO',
				'Phys_prop_file' 	: 'real_data/niwot_snow.prop',
				'Ava_angle'		: 38,
				'Lambda_S'		: 0.001})

	def __update_template(self, template_parameters):
		# Function to set special 'Csp_params' parameter
		#  __update_template([20]) called on 'Csp_template = LAYER(10)' would change LAYER(10) to LAYER(20)
		#  __update_template([5,7]) called on 'Csp_template = FORSTEP(1,2)' would set Csp_template = FORSTEP(5,7)
		initial_template = self.get('Csp_template')
		base 		 = initial_template.split('(',1)[0]
		new_template     = base + "("
		if isinstance(template_parameters, list):
			for param in template_parameters:
				new_template = new_template + str(param) + ","
			new_template = new_template[:-1] + ")" # remove trailing comma, close parenthesis
		else:
			new_template = new_template + str(template_parameters) + ")"
		self.set({'Csp_template' : new_template})

	def set(self, name_value_dict):
		# Change the value of one or more parameters
		# input as a dictionary {"Parameter name" : value}

		for name in name_value_dict.keys():	
			# Paramters that can be overwritten by other parameters must go first
			if name in self.__special_parameter_descriptions.keys():
				setting_function = self.__set_special[name]
				setting_function(name_value_dict[name])
			# Non-special parameters are independent. Order doesn't matter.
			elif name in self.parameters.keys():
				self.parameters[name] 		= name_value_dict[name]
			# Error handling
			else:
				print("Parameter " + name + " has not been initialized and has no default value.")
				print("Parameter " + name + " skipped. Check spelling or add a default with new_parameter().")
	
	def get(self, name):
		return self.parameters[name]

	def list_all(self):
		return self.parameters.keys()

	def get_description(self, name):
		return self.parameter_descriptions[name]
	
	def set_header(self, new_header):
		# Give parameter file a useful descriptive header
		# Write function will append one comment (##) symbol; additional lines
		# must be commented appropriately.
		self.header_body = new_header

        def read(self, filename):
                #Read a parameter file into all parameters
                with open(filename, "r") as f:
                        prev_line = ""
                        for i, line in enumerate(f):
                                if line.strip() == "":
                                        continue
                                elif line.startswith("#"):
                                        prev_line = line
                                        continue
                                else:
                                        try:
                                                #Capture parameter description if it has one
                                                desc = ""
                                                if prev_line != "":
                                                        desc = prev_line.replace("#","").replace("\n",'')[1:]
                                                        prev_line = ""
                                                        
                                                par = line.split(' ')

                                                #Try reading string as int, then float before assuming string
                                                value = par[2].replace('\n','')
                                                try:
                                                        value = int(par[2])
                                                except:
                                                        try:
                                                                value = float(par[2])
                                                        except:
                                                                pass
                                                
                                                self.new_parameter(par[0], desc, value)
                                        except:
                                                print("Error occurred when trying to read parameter from file on line {}.\n".format(i+1))
                                                print("Line read: {}".format(line))
                
                
	def write(self, filename):
		# Write all parameters to a parameter file
		with open(filename, "w") as f:
			f.write("## " + self.header_start + "\n")
			f.write("## " + self.header_body  + "\n \n")
			for parameter in self.parameters.keys():
				f.write("# " + self.get_description(parameter) 	+ "\n")
				f.write(parameter + " = " + str(self.get(parameter)) + "\n \n")


##----------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------

class Run_Script():
	# Creates a ReSCAL run script
	# Yes, it's a python script that writes a bash script, and yes, that's a bit silly
	def __init__(self):
		self.options = {}
		self.__flag_options = {}
		self.__default_options()
		self.__default_header()

	def __default_header(self):
		self.header_start = "Script written with Run_Script utility (KK)"
		self.header_body  = "Default options"

	def __default_options(self):
		# set default options
		# ReSCAL lists a lot of command-line options in entry.c, show_general_options()
		self.options['clean'] 		= True
		self.options['backup'] 		= False
		self.options['parfile'] 	= "sno_cone.par"
		self.options['genesislog']	= "GENESIS.log"
		self.options['rescallog']	= "RESCAL.log"
		self.options['rescallocation']	= "../src"

		# automatically handle command line options here
		self.__set_flag_option('usage info',		'h',	True)
		self.__set_flag_option('show params', 		'hm', 	True)
		self.__set_flag_option('no video',		'nv',	True)
		self.__set_flag_option('info interval', 	'info',	True)
		self.__set_flag_option('output interval', 	'dcsp', '10t0')
		self.__set_flag_option('png interval',		'dpng', '10t0')
		self.__set_flag_option('stop after',		'stop', False)
		self.__set_flag_option('frame rate', 		'fr',	False)
		self.__set_flag_option('random seed',		's',	False)
		self.__set_flag_option('vel', 			'vel', 	True)
		self.__set_flag_option('vss', 			'vss', 	True)
		
		self.options['nice']		= False
		# locations for linking to physical properties files
		self.options['real_data_location'] 	= '../../scripts/real_data'

	def __set_flag_option(self, descriptive_name, flag, value):
		# Some options are command line flags. Some are not.
		# Try to hide this bit of extra complexity from the user.
		self.__flag_options[descriptive_name] 	= flag
		self.options[descriptive_name] 		= value

	def list_all(self):
		return self.options.keys()

	def get(self, option):
		return self.options[option]

	def set_header(self, header):
		self.header_body = header

	def set(self, option_dict):
		# Takes a dictionary of {"option name" : value} pairs
		for option in option_dict.keys():
			if option in self.options.keys():
				self.options[option] = option_dict[option]
			else:
				print "Skipping nonexistent option " + option
	
	def __write_run_rescal(self, f):
		# Sub-function of write() that writes the call to rescal
		# with appropriate flags

		f.write('# ----Rescal----\n')
		if self.options['nice']:
			f.write('nice ')
		f.write("./rescal $PAR_FILE")
		for option in self.__flag_options.keys():
			if (self.options[option] == True):
				# -flag
				f.write(" -" + self.__flag_options[option])
			elif self.options[option]: # == any value except False or True
				# -flag VALUE
				f.write(" -" + self.__flag_options[option] + " " + self.options[option])
		f.write("> $RESCAL_LOG_FILE\n \n")


	def write(self, filename, input_elevation_type=None, args=None):
		with open(filename, "w") as f:
			f.write("#!/bin/bash \n \n")
			f.write("################## \n ## ReSCAL run script ## \n################")
			f.write("\n \n ##" + self.header_start + "\n ##" + self.header_body + "\n \n")
			
			# Organizational options
			f.write("# ----Organizational tasks---- \n")
			if self.options['clean']:
				f.write("# Remove files from previous runs \n")
				f.write("./clean \n\n")
			if self.options['backup']:
				f.write("# Make an archive from the sources \n")
				f.write("./dobackup \n\n")

			#Linking - not optional
			f.write("if [ ! -e genesis ]; then \n  ln -s " + self.options['rescallocation'] +  "/genesis . \nfi \n")
			f.write("if [ ! -e rescal ]; then \n  ln -s "  + self.options['rescallocation'] + "/rescal . \nfi \n\n")
			#f.write("ln -s " + self.options['rescallocation'] + "/rescal-ui.xml . \n \n")

			#Parameter file
			f.write("# ----Parameter file----\n")
			f.write("PAR_FILE=\"" + self.options['parfile'] + "\"\n\n")
			f.write("echo PAR_FILE=$PAR_FILE\n\n")

			f.write('# -----Physical properties----\n')
			f.write("cp -r " + self.options['real_data_location'] + " . \n\n")

			# Run options
			f.write("# ----Run options----\n")
			f.write("export OMP_NUM_THREADS=1 \n")
			f.write("GENESIS_LOG_FILE=\"" + self.options['genesislog'] + "\"\n")
			f.write("RESCAL_LOG_FILE=\""  + self.options['rescallog']  + "\"\n\n")

			# Run optional input_elevation scripts
			if input_elevation_type == "gaussian": 
				f.write("# ----Gaussian----\n")
				f.write("python gaussian.py " + " ".join([str(i) for i in args]) + "\n\n")

			# Run genesis
			f.write("# ----Genesis----\n")
			f.write("./genesis -f $PAR_FILE -s 2000 > $GENESIS_LOG_FILE\n\n")

			# Run rescal
			self.__write_run_rescal(f)

			# Automatic analysis


