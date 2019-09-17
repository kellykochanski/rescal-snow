# KK Jun 25 20180
# Utilities to set up run scripts quickly and easily; these tools are designed to aid parameter space
# explorations, sensitivity analyses, and large batches of runs.
# Parameters : Class to hold, update, change, or write all the parameters that ReSCAL needs to run
# RunScript : Class to hold, update, change, or write a ReSCAL run script with appropriate flags
# This includes almost all options and inputs, except those in the "real_data" file.

import os
import re
import itertools
import subprocess
import glob

class DesignRun():
    """
    Design a ReSCAL run
    This interfaces with both the Parameters and RunScript utilities
    in order to produce run script and parameter files simultaneously
    using a standardized input format
    """

    def __init__(self):
        self.name       = 'run'
        self.directory  = '.'
        self.parameters = Parameters()
        self.run_script = RunScript()
        self.parameters.set_header("Test from write_tests.py")

    def set_parameters(self, param_dict):
        """
        Sets all parameters, given as {name : value} pairs
        Allows user to mix things that go in parameter file
        with things that go in the run script, etc
        """

        for name in param_dict.keys():
            if self._is_a_parameter(name):
                self.parameters.set({name : param_dict[name]})
            elif self._is_a_run_script_option(name):
                self.run_script.set({name : param_dict[name]})
            else:
                print("Warning : skipped nonexistent parameter " + name)

    def list_all(self):
        """Lists all available options/parameters"""
        return self.run_script.list_all(), self.parameters.list_all()

    def get(self, param):
        if param in self.parameters.list_all():
            return self.parameters.get(param)
        elif param in self.run_script.list_all():
            return self.run_script.get(param)
        else:
            print("Error: requested nonexistent parameter " + param)
            return False

    def set_header(self, description):
        self.parameters.set_header("Test from write_tests.py: " + description)

    def set_name(self, name):
        self.name = name

    def set_directory(self, directory):
        self.directory = directory

    def _is_a_parameter(self, name):
        # Checks whether the parameters dictionary contains 'name'
        return ((name == 'Environment') or (name in self.parameters.list_all()))

    def _is_a_run_script_option(self, name):
        # Checks whether 'name' is an option in the run script
        return (name in self.run_script.list_all())

    def write(self):
        self.run_script.write(self.directory + "/" + self.name + ".run")
        self.parameters.write(self.directory + "/" + self.name + ".par")

##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

class Parameters():
    """Holds all the parameters that ReSCAL needs to run"""
    def __init__(self):
        self.parameters                 = {}
        self.parameter_descriptions     = {}
        self._special_parameter_descriptions = {}
        self._set_special =                    {}
        self._default_header()
        self._default_parameters()

    def new_parameter(self, name, description, default_value):
        """
        New parameter
        Name should be a string (the parameter name expected by ReSCAL)
        Description should be a string describing its use or purpose etc
        Default value is the default value, generally a string or float
        """
        self.parameters[name]             = default_value
        self.parameter_descriptions[name] = description

    def _special_parameter(self, name, description, setting_function):
        """
        Special parameters are not stored.
        When set, they modify one or more other parameters in a way described by the setting_function()
        """
        self._special_parameter_descriptions[name]     = description
        self._set_special[name]                        = setting_function

    def _default_header(self):
        """Provides a basic description of this parameter set in form of a string"""
        self.header_start = "Parameter file for ReSCAL written from Parameters utility"
        self.header_body  = "Default parameters : sahara sand"

    def _default_parameters(self):
        """
        Default values for parameters
        Initialized here for sahara sand
        Basic model parameters and geometry of model run:
        """
        self.new_parameter('Model', 'Defines a set of transitions for grains in the cellular automaton e.g. DUN or SNO', 'SNO')
        self.new_parameter('Output_directory', 'Directory where output files go', './out')

        self.new_parameter('Csp_file',  'Genesis template?',    'DUN.csp')
        self.new_parameter('Csp_template', 'Template containing initial shape for run', 'RCONE')
        self.new_parameter('Csp_params', 'Parameters for size, depth etc of template shape', [])
        self.new_parameter('Boundary',  'Type of boundary conditions, e.g. OPEN, PERIODIC, REINJECTION', 'PERIODIC')
        self.new_parameter('Time',      'Initial time',         0.0)
        self.new_parameter('H',         'Height of model domain (cells)',       80)
        self.new_parameter('L',         'Length, longitudinal, of model domain (cells)',        300)
        self.new_parameter('D',         'Depth, transverse to flow, of model domain (cells)',   200)
        self.new_parameter('Centering_delay', 'Automatic re-centering of model domain on dune', 0)

        # Parameters describing the physical behavior of the grains in the model
        self.new_parameter('Phys_prop_file', 'File containing important physical properties e.g. gravity', 'real_data/desert_earth.prop')
        self.new_parameter('Qsat_file', 'File containing relationship between saturated and threshold flux', 'real_data/PDF.data')
        self.new_parameter('Lambda_E',  'Erosion rate (when threshold exceeded)', 1)
        self.new_parameter('Lambda_T',  'Transition rate for mobile grains', 1.5)
        self.new_parameter('Lambda_C',  'Deposition rate for mobile grains on surface', 0.5)
        self.new_parameter('Lambda_G',  'Gravity - rate of fall of not-in-transport grains in air', 1000)
        self.new_parameter('Lambda_D',  'Transition rate for diffusion', 0.01)
        self.new_parameter('Lambda_S',  'Transition rate for sintering/cohesion', 0)
        self.new_parameter('Lambda_F',  'Ratio of sintered:unsintered erosion thresholds', 3)
        self.new_parameter('Coef_A',    'Ratio of vertical:horizontal transport of mobile grains', 0.1)
        self.new_parameter('Coef_B',    'Ratio of deposition against an obstacle : deposition', 1)
        self.new_parameter('Coef_C',    'Ratio of deposition behind an obstacle : deposition', 3)
        self.new_parameter('Prob_link_ET', 'Probability of the transition links', 0.5)
        self.new_parameter('Prob_link_TT', 'Probability of the transition links', 1)
        self.new_parameter('High_mobility', 'Higher mobility of grains', 1)
        self.new_parameter('Lambda_I',  'Injection rate', 0)
        self.new_parameter('Lambda_A',  '??', 1)

        # Parameters describing avalanching
        self.new_parameter('Ava_mode',  'Mode of avalanching', 'TRANS')
        self.new_parameter('Ava_angle', 'Angle of repose, steeper slopes avalanche (degrees)', 35)
        self.new_parameter('Ava_h_lim', 'Height limit in avalanches (cells)', 1)

        # Parameters for flow stabilization - how often the lgca runs and how strong is the flow
        self.new_parameter('Lgca_delay', 'Delay between flow cycles', 1)
        self.new_parameter('Lgca_speedup', 'Speedup of the flow stabilization', 1000)
        self.new_parameter('Tau_min',   'Shear stress threshold for grain erosion, controls flow strength', 0)
        self.new_parameter('Tau_max',   'Max shear stress, used in grain erosion/flow speed relation', 1000)

        # Special parameters - set by same interface as others,
        # But are actually tools which modify one or more other parameters
        self._special_parameter('Environment', 'Sets all defaults for an environment, e.g. sand or snow', self._environment)
        self._special_parameter('Csp_params',  'Sets parameters for a Csp_template', self._update_template)

    def _environment(self, keyword):
        """
        Function to set special 'Environment' parameter
        Update defaults fora  different enviroment, e.g. snow
        (More environments could be added for subaqueous dunes, mars environment, etc)
        """
        if keyword == "snow":
            self.header_body = "Default values - Niwot Ridge snow"
            self.set({'Model'               : 'SNO',
                    'Phys_prop_file'        : 'real_data/niwot_snow.prop',
                    'Ava_angle'             : 38,
                    'Lambda_S'              : 0.001})

    def _update_template(self, template_parameters):
        """
        Function to set special 'Csp_params' parameter
         _update_template([20]) called on 'Csp_template = LAYER(10)' would change LAYER(10) to LAYER(20)
         _update_template([5,7]) called on 'Csp_template = FORSTEP(1,2)' would set Csp_template = FORSTEP(5,7)
         """
        initial_template = self.get('Csp_template')
        base             = initial_template.split('(',1)[0]
        new_template     = base + "("
        if isinstance(template_parameters, list):
            for param in template_parameters:
                new_template = new_template + str(param) + ","
            new_template = new_template[:-1] + ")" # remove trailing comma, close parenthesis
        else:
            new_template = new_template + str(template_parameters) + ")"
        self.set({'Csp_template' : new_template})

    def set(self, name_value_dict):
        """
        Change the value of one or more parameters
        input as a dictionary {"Parameter name" : value}
        """

        for name in name_value_dict.keys():
            # Paramters that can be overwritten by other parameters must go first
            if name in self._special_parameter_descriptions.keys():
                setting_function = self._set_special[name]
                setting_function(name_value_dict[name])
            # Non-special parameters are independent. Order doesn't matter.
            elif name in self.parameters.keys():
                self.parameters[name]           = name_value_dict[name]
            # Error handling
            else:
                print("Parameter " + name + " has not been initialized and has no default value.")
                print("Parameter " + name + " skipped. Check spelling or add a default with new_parameter().")

    # TODO, deal with all parameters and verify none are real_data/*
    def set_real_data_path(self, real_data_directory):
        """sets path for 'Phys_prop_file' and 'Q_sat_file"""
        keys = ['Phys_prop_file','Qsat_file']
        for key in keys:
            self.parameters[key] = os.path.join(real_data_directory, os.path.basename(self.parameters[key])) 
                
    def get(self, name):
        return self.parameters[name]

    def list_all(self):
        return self.parameters.keys()

    def get_description(self, name):
        return self.parameter_descriptions[name]

    def set_header(self, new_header):
        """
        Give parameter file a useful descriptive header
        Write function will append one comment (##) symbol; additional lines
        must be commented appropriately.
        """
        self.header_body = new_header

    def read(self, filename):
        """Read a parameter file into all parameters"""
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
        """Write all parameters to a parameter file"""
        with open(filename, "w") as f:
            f.write("## " + self.header_start + "\n")
            f.write("## " + self.header_body  + "\n \n")
            for parameter in self.parameters.keys():
                f.write("# " + self.get_description(parameter)  + "\n")
                f.write(parameter + " = " + str(self.get(parameter)) + "\n \n")


##----------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------

class RunScript():
    """
    Creates a ReSCAL run script
    Yes, it's a python script that writes a bash script.
    """
    def __init__(self):
        self.options = {}
        self._flag_options = {}
        self._default_options()
        self._default_header()

    def _default_header(self):
        self.header_start = "Script written with RunScript utility (KK)"
        self.header_body  = "Default options"

    def _default_options(self):
        """
        set default options
        ReSCAL lists a lot of command-line options in entry.c, show_general_options()
        """
        self.options['clean']           = True
        self.options['backup']          = False
        self.options['parfile']         = "sno_cone.par"
        self.options['genesislog']      = "GENESIS.log"
        self.options['rescallog']       = "RESCAL.log"
        self.options['rescallocation']  = "../src"

        # automatically handle command line options here
        self._set_flag_option('usage info',             'h',    True)
        self._set_flag_option('show params',            'hm',   True)
        self._set_flag_option('no video',               'nv',   True)
        self._set_flag_option('info interval',          'info', True)
        self._set_flag_option('output interval',        'dcsp', '10t0')
        self._set_flag_option('png interval',           'dpng', '10t0')
        self._set_flag_option('stop after',             'stop', False)
        self._set_flag_option('frame rate',             'fr',   False)
        self._set_flag_option('random seed',            's',    False)
        self._set_flag_option('vel',                    'vel',  True)
        self._set_flag_option('vss',                    'vss',  True)
        self._set_flag_option('quit',                   'q',   False)
        self._set_flag_option('alti only',              'altionly', False)
        self._set_flag_option('cellspace borders',      'csp_borders', False)
        self._set_flag_option('uncompressed cellspace', 'uncompressed_csp', False)
        self._set_flag_option('print performance',      'perf_print',  False)
        

        self.options['nice']            = False
        # locations for linking to physical properties files
        self.options['real_data_location']      = '../../scripts/real_data'

        # can use path to premade csp file instead of having
        # genesis make one, will still use 'Dun.csp' as the link
        self.options['premade_csp'] = False

    def _set_flag_option(self, descriptive_name, flag, value):
        """
        Some options are command line flags. Some are not.
        Try to hide this bit of extra complexity from the user.
        """
        self._flag_options[descriptive_name]   = flag
        self.options[descriptive_name]          = value

    def list_all(self):
        return self.options.keys()

    def get(self, option):
        return self.options[option]

    def set_header(self, header):
        self.header_body = header

    def set(self, option_dict):
        """Takes a dictionary of {"option name" : value} pairs"""
        for option in option_dict.keys():
            if option in self.options.keys():
                self.options[option] = option_dict[option]
            else:
                print("Skipping nonexistent option " + option)

    def _write_run_rescal(self, f):
        """
        Sub-function of write() that writes the call to rescal
        with appropriate flags
        """

        f.write('# ----Rescal----\n')
        if self.options['nice']:
            f.write('nice ')
        f.write("./rescal $PAR_FILE")
        for option in self._flag_options.keys():
            if (self.options[option] == True):
                # -flag
                f.write(" -" + self._flag_options[option])
            elif self.options[option]: # == any value except False or True
                # -flag VALUE
                f.write(" -" + str(self._flag_options[option]) + " " + str(self.options[option]))
        f.write("\n \n")



    def rescal_call_args(self):
        """Give a list of all the args for the actual call to rescal"""
        args = []
        if self.options['nice']:
            nice_list.append('nice')
        for option in self._flag_options.keys():
            if (self.options[option] == True):
                # -flag
                args.append("-" + self._flag_options[option])
            elif self.options[option]: # == any value except False or True
                # -flag VALUE
                args.append("-" + str(self._flag_options[option]))
                args.append(str(self.options[option]))
        return args
        
        

    def write(self, filename):
        """Write a rescal run script from a RunScript object"""
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
            f.write("ln -s " + self.options['rescallocation'] + "/rescal-ui.xml . \n \n")

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

            # Run genesis if no premade .csp file
            if  not self.options['premade_csp']:
                f.write("# ----Genesis----\n")
                f.write("./genesis -f $PAR_FILE -s 2000 > $GENESIS_LOG_FILE\n\n")
            else:
                # link to the premade .csp
                f.write('ln -s ' + self.options['premade_csp']  + ' DUN.csp\n\n')



            # Run rescal
            self._write_run_rescal(f)

        # make file executable
        os.chmod(filename, 0o764)


##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------



# TODO allow iteratables to be used
# not just lists
# TODO allow passing of arguments in dicts, but with no guaranteed ordering
def all_parameter_combos(named_parameter_lists):
    """
    Takes a parameter list of lists of lists such as:
       [['Coef_A', [0.1, 0.2, 0.3]], ['random seed', [1,2,3]]])
       so that each sublist has the variable name and set of parameters to vary over.
       Creates a generator that returns a dictionary with each possible parameter
       combination and a string that can be used to identify the varying parameters.
       If variable names have spaces they are replaced with underscores.
    """

    # separate parameter names from values lists, cast all values to strings
    names = [x[0] for x in named_parameter_lists]
    parameter_lists_original = [x[1] for x in named_parameter_lists]

    # cast all values to strings
    parameter_lists = [[str(i) for i in x] for x in parameter_lists_original]

    # create a generator that makes all combination of parameters lists (the cartesian product)
    # also creates a string for the file name
    for parameter_combination in itertools.product(*parameter_lists):
        # create a dict of the names and some value combination
        params_dict_to_add = dict(zip(names, parameter_combination))
        # creats a directory name based on the names and values
        # if any of the values are file paths, only the base name of the file is used
        directory_name_suffix = [name.replace(' ', '_') + '-' + \
                                 os.path.basename(str(params_dict_to_add[name])) for name in names]
        yield params_dict_to_add, '--'.join(directory_name_suffix)



#### NOTE: set up for power-lab machine usage, meaning on a single PC, not a cluster
def make_run_directories(fixed_params, variable_params, experiment_name, run_header, run_name):
    """
    Create a set of directories that contain all needed files for separate rescal runs.
       The directories exist in a top level directory of experiment_name.
       Each directory is names based on its varying parameters.
       The paths to the run scripts are returned for easy execution on the
       power-lab machines.
    """

    # make the top level directory, but don't overwrite a directory that
    # already exists
    experiment_directory_name = experiment_name
    if not os.path.isdir(experiment_name):
        os.mkdir(experiment_name)
    # TODO deal with case that directory already exists
    else:
        return
    
    # save the paths to the run scipts for running them later
    run_scripts = []

    # create full parameter set for each run
    # as well as the directory suffix for each run sub directory
    for current_variable_params in all_parameter_combos(variable_params):
        params_to_add, directory_suffix = current_variable_params
        parameters = {**fixed_params, **params_to_add}
        this_run = DesignRun()
        #this_run.set_header(run_header)
        this_run.set_name(run_name)
        # create a directory for each run inside the experiment directory
        run_directory = os.path.join(experiment_directory_name, directory_suffix)
        if not os.path.isdir(run_directory):
            os.mkdir(run_directory)

        this_run.set_directory(run_directory)
        this_run.set_parameters(parameters)
        this_run.write()

        # store all the run scripts for future use
        run_scripts.append(os.path.join(run_directory, run_name + '.run'))

    return run_scripts


#### NOTE: set up for power-lab machine usage, meaning on a single PC, not a cluster
def run_rescals(run_scripts):
    """
    Takes a list of paths to .run scripts that should already be in directories set up
    to run ReSCAL. An instance of ReSCAL is started using each run script and then
    the ReSCAL all run at the same time and asynchronously. This function waits for all the
    child processes to complete.
    """

    # if any spaces in path names, turns ' ' into '\ ' so the shell can understand them
    modded_scripts = []
    for run_script in run_scripts:
        modded_scripts.append(run_script.replace(' ', '\\ '))

    processes = []
    # move to the directory of each run_script and run it
    for modded_script in modded_scripts:
        processes.append(subprocess.Popen('cd ' + os.path.dirname(modded_script) + \
                                          ' && ' +  modded_script, shell=True))
    for p in processes:
        p.wait()


def par_to_dict(filename):
    """
    given a .par file, returns a dictionary of parameters
    TODO, could be more flexible, but works for the auto-generated ones
    """
    with open(filename, 'r') as f:
        data = f.read()
    
        # get all the assignment strings
        value_assignment = re.compile(r'([^#\n]*)[ ]*=[ ]*([^#\n]+)')
        return dict(re.findall(value_assignment, data))
        
def cmd_line_args(filename):
    """given a rescal .run file, get the command-line arguments for rescal"""
    with open(filename, 'r') as f:
        data = f.read()
        # get the rescal line
        
        rescal_line = re.search(r'(\./rescal.*)\n', data)
        # now get the args out
        flag_arg = re.compile(r'(-[^\s]*)+(?:[ \t]+|\n)([^-\s]*)')
        flag_value_pairs = flag_arg.findall(rescal_line.group(0))
        return dict(flag_value_pairs)

def get_files_to_process(top_dir, path_glob, exclude_globs):
    """
    given a top directory top_dir
    all of its subdirectories will be checked for files
    that match some the path_glob
    so any file of name top_dir/"any subdirectory of top_dir"/path_glob
    will match unless that file name matches some glob_exclude
    examples:
    for cellspace files path_glob='*.csp*' exclude_globs=['DUN.csp']
    for ALTI files      path_glob='ALTI*'  exclude_globs=[]
    """
    
    
    # get all the directories in top_directory (from stack overflow 973473)
    # make them absolute paths
    output_dirs = sorted([f.path for f in os.scandir(top_dir) if f.is_dir()])
    for i in range(len(output_dirs)):
        output_dirs[i] = os.path.abspath(output_dirs[i])

    # make 2D array or absolute paths to all the .csp* files
    paths = []
    for output_dir in output_dirs:
        # get all possible matches
        local_paths = glob.glob(output_dir + '/' + path_glob)
        # find the ones to exclude
        exclude_paths = []
        for exclude_glob in exclude_globs:
            exclude_paths += glob.glob(output_dir + '/' + exclude_glob)

        # now remove the bad paths from the good ones
        local_paths = sorted(list(set(local_paths) - set(exclude_paths)))
        
        # if this sub_dir has files to add, add them
        if local_paths:
            paths.append(local_paths)
    # if no files were added, paths will not have changed
    # however, always returs a 2D list
    if paths == []:
        return [[]]
    
    # truncate paths so that each row has the same number
    # it's possible that the output_dirs won't all have the
    # same number of .csp files if the simulation doesn't finish
    # so just limit row size to minimum of any output directory
    # there shouldn't be any empty lists in paths
    paths_truncated = []
    min_files_in_dir = min([len(p) for p in paths])
    for p in paths:
        paths_truncated.append(p[:min_files_in_dir])
    return paths_truncated

def random_initial_states(num_states, parameters, top_dir, run_header='run', run_name='run'):
    """create num_states random seeds"""
    seed_numbers = random.sample(range(1,1000000), num_states)
    seeds = [['random seed', seed_numbers]]
    
    run_files = make_run_directories(parameters, seeds, top_dir, run_header, run_name)

    run_rescals(run_files)

    paths = get_files_to_process(top_dir, )
