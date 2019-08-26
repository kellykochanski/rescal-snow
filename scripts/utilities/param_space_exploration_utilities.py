# KK July 03 2018

# This is a shell script designed to start organizing ReSCAL to perform a parameter space exploration;
#  it could be adapted to perform any set of runs with variable parameters

# The "explore_parameter_space" function creates an array of subdirectories, each containing:
#   - a run script (bash), generally called run.run
#   - a parameter file (txt), generally caled run.par
# Note that the multi-directory structure is useful for several reasons:
#   1. The output for each run stays in its own subdirectory, attached to the appropriate input files
#   2. If multiple processors are available, they may simultaneously run separate instances of rescal-snow, one per directory,
#       without any communication, race conditions, or other problems

# Set up a bunch of runs to go in parallel

import numpy as np
import datetime
import os
import itertools
import subprocess
import glob
import rescal_utilities


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
        run_output_file = '{year}{month}{day}_{run}'.format(year=str(now.year),
                                                            month=str(now.month),
                                                            day=str(now.day),
                                                            run=str(run))
        run_output_dir = os.path.join(output_root, run_ouput_file)

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
            valuerange      = parameter_ranges[parameter]
            these_parameters[parameter] = np.random.uniform(valuerange[0], valuerange[1])

        parameter_combinations.append(these_parameters)

    return parameter_combinations







# TODO allow iteratables to be used
# not just lists
# TODO allow passing of arguments in dicts, but with no guaranteed ordering
def all_parameter_combos(named_parameter_lists):
    '''Takes a parameter list of lists of lists such as:
       [['Coef_A', [0.1, 0.2, 0.3]], ['random seed', [1,2,3]]])
       so that each sublist has the variable name and set of parameters to vary over.
       Creates a generator that returns a dictionary with each possible parameter
       combination and a string that can be used to identify the varying parameters.
       If variable names have spaces they are replaced with underscores.'''

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
    '''Create a set of directories that contain all needed files for separate rescal runs.
       The directories exist in a top level directory of experiment_name.
       Each directory is names based on its varying parameters.
       The paths to the run scripts are returned for easy execution on the
       power-lab machines.'''

    # make the top level directory, but don't overwrite a directory that
    # already exists
    experiment_directory_name = experiment_name
    if not os.path.isdir(experiment_directory):
        os.mkdir(experiment_directory)
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
        this_run = rescal_utilities.Design_a_run()
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
    '''Takes path to a set of run_scripts that should be in the proper locations
       to run rescal and runs each of them. Runs should be asynchronous.'''

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


# given a top directory top_dir
# all of its subdirectories will be checked for files
# that match some the path_glob
# so any file of name top_dir/"any subdirectory of top_dir"/path_glob
# will match unless that file name matches some glob_exclude
# examples:
# for cellspace files path_glob='*.csp*' exclude_globs=['DUN.csp']
# for ALTI files      path_glob='ALTI*'  exclude_globs=[]
# TODO, 2 copies of this function, the other in param_space_exploration_utilities.py
def get_files_to_process(top_dir, path_glob, exclude_globs):
    
    # get all the directories in top_directory (from stack overflow 973473)
    # make them absolute paths
    output_dirs = sorted([f.path for f in os.scandir(top_dir) if f.is_dir()])
    for i in range(len(output_dirs)):
        output_dirs[i] = os.path.abspath(output_dirs[i])

    # make 2D array or absolute paths to all the .csp* files
    paths = []
    for output_dir in output_dirs:
        # get all possible matches
        local_paths = glob.glob(os.path.join(output_dir, path_glob))
        # find the ones to exclude
        exclude_paths = []
        for exclude_glob in exclude_globs:
            exclude_paths += glob.glob(os.path.join(output_dir, exclude_glob))

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


