
import numpy as np
import random
import os
import rescal_utilities
import cellspace
import pickle


# experiment directory setup
# expermient directory is a directory that contains all the files
# it is dealt with as an object in python
# RescalRun


# experiment_name
#   meta_data 
#
#



class RescalDataRunMetaData:

    def __init__(self):
        self.parameters = None
        self.intial_csps_created = False
        self.control_run_files_created = False        
        self.modded_csp_created = False
        self.moded_run_files_created = False
        self.main_experiment_started = False
        self.data_consolidated = False
        self.jobs_ids = []

    
# a class that represents a rescal data run
# uses directory in the top level directory of a ReSCAL install
class RescalDataRun:

    def __init__(self, experiment_directory, rescal_root=None, overwrite=False):
        # get the location of rescal
        if rescal_root is None:
            self.rescal_root = os.environ['RESCAL_ROOT']
        else:
            self.rescal_root = rescal_root


        # deal with experiment_directory
        if not os.path.isabs(experiment_directory):
            self.experiment_directory = os.path.join(self.rescal_root, experiment_directory)
        else:
            # TODO, this could break everything if user gives bad path
            self.experiment_directory = experiment_directory


        self.meta_data_path = os.path.join(self.experiment_directory, '.meta_data')

            
        # see if experiment_directory exists and is a directory
        # TODO deal with permissions
        if os.path.isdir(self.experiment_directory):
            # look for the meta_data file
            if os.path.isfile(self.meta_data_path):
                # if found, use it to initialize
                self.initialize_from_existing_experiment_directory()
                
        # failure, won't everwrite file
        elif os.path.isfile(self.experiment_name):
            print('A file exists at {ed}'.format(self.experiment_directory))
            print('Cannnot create experiment directory')

        # nothing exists at self.experiment directory, so make and setup new directory
        else:
            self.initialize_experiment_directory()
            


    # TODO maybe do more than just the meta_data
    def initialize_experiment_directory(self):
        self.create_meta_data()
        

    # create a meta_data object and file
    def create_meta_data(self):
        self.meta_data = RescalDataRunMetaData()
        pickle.dump(self.meta_data, self.meta_data_path)
    
            
    def initialize_from_existing_experiement(self):
        self.meta_data = pickle.load(self.meta_data_path)

        
        
    # setup directories, this includes running rescal
    # to randomly initialize
    def setup(self):
        pass

    













def flatten(x):
    if x == []:
        return []
    elif isinstance(x, list):
        return flatten(x[0]) + flatten(x[1:])
    else:
        return [x]
        



# common parameters for rescal
# this includes parameters that will be in the .par file
parameters_par_1 = {
    'Model':  'SNO',
    'Output_directory': './out',
    'Csp_file': 'DUN.csp',
    'Csp_template': 'SNOWFALL(4)',
    'parfile': 'run.par',
    'Boundary':  'OPEN',
    'Time': 0.0,
    'H': 80,
    'L': 400,
    'D': 200,
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
    'Lambda_I': 0.002,
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
    'Tau_min': 100,
    'Tau_max': 1100,
    'rescallocation': '../../scripts', # this is where run script looks for rescal+genesis executables
}

# command line arguments for rescal which will be in the .run file
parameters_run_1 = {
    'stop after' : '400_t0',
    'output interval' : '200_t0',
    'png interval' : False,
    'quit' : True,
    'random seed' : 6,
    'usage info' : False,
    'show params' : False,
    'info interval' : False,
}

# parameters that are not given to genesis or rescal
# but affect how these files are made
# e.g. 'premade_csp' will modify the run file so that genesis is not called and
# and premade .csp is used as the input to rescal
parameters_meta_1 = {
    'rescallocation': '../../scripts', # this is where run script looks for rescal+genesis executables
    #'premade_csp' : '/g/g13/defazio1/rescal-snow/scripts/SNO2500_t0_rng_4.csp'
}

# combine the 3 kinds of parameters into a single dict
# there should be any common keys
parameters_1 = {**parameters_par_1, **parameters_run_1, **parameters_meta_1}


# parameters that will differ between runs
# a list of lists
# each list has a key, which is some valid parameter name
# and a list of values for that key to take
# the total number of runs is the cardinality of the cartesian product
# of all the values lists
var_params_1 = [['random seed', [7,8,9,10,11,12,13,14,15,16]],
              ['Tau_min', [100, 110, 120]],
]


# parameters that wil change for the second set of runs
# these are the modified parameters
parameters_par_mod = {

}

parameters_run_mod = {

}

parameters_meta_mod = {

}

parameters_mod = {**parameters_par_mod, **parameters_run_mod, **parameters_meta_mod}



# makes randoms initial states
# sets up runs directories to make initial states
# runs the simulations to get the initial states
# returns the file paths to the .csp files created
def random_initial_states(num_states, parameters, top_dir, run_header='run', run_name='run'):

    # create num_states random seeds
    seed_numbers = random.sample(range(1,1000000), num_states)
    seeds = [['random seed', seed_numbers]]

    if run_files == None:
        run_files = rescal_utilities.make_run_directories(parameters, seeds, top_dir, run_header, run_name)

    rescal_utilities.run_rescals(run_files)

    paths = rescal_utilities.get_files_to_process(top_dir, cellspace.path_glob, cellspace.exclude_globs)
    paths_last = [x[-1] for x in paths]


    print(paths_last)
    
    return paths_last


# create modified versions of the outputs a bunch of .csp files
# paths to the .csp files
# makes a 3D list by [starting_csp][mod_type][mod_num]
def modify_outputs(top_dir, paths, mod_types, num_mods, exp_name='exp', ):


    # top dir to put the .csp files
    if not os.path.isdir(top_dir):
        os.mkdir(top_dir)

        
    # the 3D list of paths to all the mods
    mod_paths = []
    for i,path in enumerate(paths):
        c = cellspace.CellSpace(path)
        mod_paths_inner = []
        for mod_type in mod_types:
            # entire path name will be:
            # top_dir + '/' + exp_name + str(i) + '--' + mod_type '--' str(j) + '.csp'
            # where j indexes num_mods
            file_prefix = top_dir + '/' +  str(i) 
            mod_paths_inner.append(c.random_mods(file_prefix, mod_type,
                                                 num_mods=num_mods, temp_mod=True))
        mod_paths.append(mod_paths_inner)

    print(mod_paths)
        
    return mod_paths



# create runs scripts for the control run
def make_control_runs(num_runs, paths1D, parameters, top_dir, run_header='run', run_name='run', seeds=None):
    
    # top dir to put the .csp files
    if not os.path.isdir(top_dir):
        os.mkdir(top_dir)

    # create num_states random seeds
    if seeds == None:
        seed_numbers = random.sample(range(1,1000000), num_runs)
        seeds = [['random seed', seed_numbers]]

    for path in paths1D:    
        path_file_part = os.path.basename(path)[:-4]
        file_prefix = top_dir + '/' + path_file_part
        # set each run up to use a premade .csp
        csp_loc = {'premade_csp' : path}
        parameters_new = {**parameters, **csp_loc}
        all_paths.append(rescal_utilities.make_run_directories(parameters_new,
                                                               seeds, 
                                                               file_prefix,
                                                               run_header,
                                                               run_name))
    return all_paths
    
    




# takes the paths to the modded .csp files, creates a bunch of run folders for each with
# differnt random seeds
def make_modded_runs(num_runs, paths3D, parameters, top_dir, run_header='run', run_name='run', seeds=None):

    # top dir to put the .csp files
    if not os.path.isdir(top_dir):
        os.mkdir(top_dir)

    
    # create num_states random seeds
    if seeds == None:
        seed_numbers = random.sample(range(1,1000000), num_runs)
        seeds = [['random seed', seed_numbers]]

    # 4D array of paths once filled
    all_paths = []

    # paths comes from modify outputs, so it's 3D
    for paths2D in paths3D:
        for paths1D in paths2D:
            for path in paths1D:
                # set top_dir
                # remove .csp from paths names
                # TODO, only use filename of path for file_prefix
                path_file_part = os.path.basename(path)[:-4]
                file_prefix = top_dir + '/' + path_file_part
                # set each run up to use a premade .csp
                csp_loc = {'premade_csp' : path}
                parameters_new = {**parameters, **csp_loc}
                all_paths.append(rescal_utilities.make_run_directories(parameters_new,
                                                                       seeds, 
                                                                       file_prefix,
                                                                       run_header,
                                                                       run_name))
    return all_paths
    

# sets up a set of files to be run
# starts with some initial conditions
# creates num_rand_initials
# then |mod_types| x num_mods 
# then num_runs
# the result is lists of file paths, a 1D lists of controls
# and a 4D list of modded versions
# the 1D list is length rand_initials
# the 4D list is num_rand_initials x |mod_types| x num_mods x num_runs
def set_up_data_run(parameters_initial, parameters_after_mod,
                    num_rand_initials, mod_types, num_mods, num_runs,
                    top_dir_rand_initials, top_dir_mods, top_dir_experiment_runs,
                    rand_initial_files=None, chunk_size=128):

    # start with some initial configuration
    # run it to some time and get the csp files from the end of the run
    if rand_initial_files == None:
        rand_initial_files = random_initial_states(num_rand_initials,
                                                   parameters_initial,
                                                   top_dir_rand_initials)

    # create run files for the rand_initial_files
    # this is the control group
    make_control_runs()
        
    # create a bunch of modified versions of each .csp
    modded_csps = modify_outputs(top_dir_mods, rand_initial_files, mod_types, num_mods)

    # create run directories for all the modded .csp files, using differnt seeds
    run_paths = make_modded_runs(num_runs, modded_csps, parameters_after_mod,
                                 top_dir_experiment_runs)

    print(run_paths)
    make_job_files(run_paths, chunk_size)

    return rand_initial_files, modded_csps, run_paths
    



# does the data run
def do_data_run():




# makes pairs of .sbatch and submit files to split up jobs
def make_job_files(run_paths, chunk_size=36, sbatch_file='test', submit_file='submit'):
    # split up run_paths2d into a few chunks
    # first flatten and remove 'run.run'

    dirs = flatten(run_paths)
    
    # flatten run paths
    #for run_paths1D in run_paths2D:
    #for run_path in run_paths2D:
    #    dirs.append(os.path.dirname(run_path))
    # create the chunks
    chunks = [dirs[i:i+chunk_size] for i in range(0,len(dirs),chunk_size)]

    sbatch_files = []
    submit_files = []
    for i, chunk in enumerate(chunks):
        current_sbatch_file = sbatch_file + str(i) + '.sbatch'
        sbatch_files.append(current_sbatch_file)
        current_submit_file = submit_file + str(i) + '.sh'
        submit_files.append(current_submit_file)
        make_sbatch_file(chunk, current_sbatch_file, current_submit_file, run_number=i)
        make_submit_file(chunk, current_submit_file)
    return sbatch_files, submit_files
        

# the sbatch file that does the big set of parallel runs
def make_sbatch_file(run_paths, output_file='test', submit_file='submit', run_number=0): 
    
    with open(output_file, 'w') as f:
        # get the data needed
        email = 'defazio1@llnl.gov'
        ntasks = len(run_paths)
        time = '05:00:00'
        log_file = 'log' + str(run_number) + '.txt'
        
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --job-name r{rn}\n'.format(rn=run_number))
        f.write('#SBATCH --output={lf}\n'.format(lf=log_file))
        f.write('#SBATCH --qos normal\n')
        f.write('#SBATCH --mail-user={e}\n'.format(e=email))
        f.write('#SBATCH --mail-type=ALL\n')      
        f.write('#SBATCH --time={t}\n'.format(t=time))
        f.write('#SBATCH --ntasks {nt}\n'.format(nt=ntasks))
        f.write('srun --wait=0 --cpus-per-task=1 --ntasks={nt} {sf}'.format(nt=ntasks, sf=submit_file))
    # set permissions for output_file
    os.chmod(output_file, 0o764)
        

# the submit.h file that sbatch calls
def make_submit_file(dirs, output_file='submit.sh'):
    with open(output_file, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('run_dirs=(\n')
        print(dirs)
        for d in dirs:
            f.write(d+'\n')
        f.write(')\n')
        f.write('my_run_dir=\"{run_dirs[${PMI_RANK}]}\"\n')
        f.write('cd ${my_run_dir}\n')
        f.write('chmod u+rwx *\n')
        f.write('./run.run\n')
    # set permissions for output file
    os.chmod(output_file, 0o764)    



if __name__ == '__main__':

    p_mods = {
        'stop after' : '3600_t0',
        'output interval' : '100_t0',
        'alti only' : True,
    }


    
    parameters_2 = {**parameters_1, **p_mods}


    rif_s = ['/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-175495/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-349/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-447778/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-475580/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-495303/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-525175/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-582451/SNO00002_t0.csp.gz', '/g/g13/defazio1/summer_2019/rescal-snow/tdi/random_seed-909448/SNO00002_t0.csp.gz'] 
    
    # rif, mc, rp = set_up_data_run(parameters_1, parameters_2,
    #                               8, ['space_invader', 'sine'], 8, 8,
    #                               '../../tdi', '../../tdm', '../../tdr',
    #                               rand_initial_files=rif_s)
    
#    f = open('paths_saved')
#    li = f.read().splitlines()
    
    #print(li)

 #   sb, su = make_job_files(li, 128)
 #   print(sb)
 #   print(su)

