import numpy as np
import random
import os
import rescal_utilities
import cellspace
import pickle


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
    'stop after' : '50_t0',
    'output interval' : '50_t0',
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
}

# combine the 3 kinds of parameters into a single dict
# there should be any common keys
parameters_1 = {**parameters_par_1, **parameters_run_1, **parameters_meta_1}


# makes randoms initial states
# sets up runs directories to make initial states
# runs the simulations to get the initial states
# returns the file paths to the .csp files created
def random_initial_states(num_states, parameters, top_dir, run_header='run', run_name='run', run_files=None):

    # create num_states random seeds
    seed_numbers = random.sample(range(1,1000000), num_states)
    seeds = [['random seed', seed_numbers]]

    if run_files == None:
        run_files = rescal_utilities.make_run_directories(parameters, seeds, top_dir, run_header, run_name)

    rescal_utilities.run_rescals(run_files)

    paths = rescal_utilities.get_files_to_process(top_dir, cellspace.path_glob, cellspace.exclude_globs)
    paths_last = [x[-1] for x in paths]

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
    #make_control_runs()
    
    # create a bunch of modified versions of each .csp
    modded_csps = modify_outputs(top_dir_mods, rand_initial_files, mod_types, num_mods)

    # create run directories for all the modded .csp files, using differnt seeds
    run_paths = make_modded_runs(num_runs, modded_csps, parameters_after_mod,
                                 top_dir_experiment_runs)

    print(run_paths)
    make_job_files(run_paths, chunk_size)

    return rand_initial_files, modded_csps, run_paths
    



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
        time = '00:30:00'
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
            f.write(os.path.dirname(d)+'\n')
        f.write(')\n')
        f.write('my_run_dir=\"${run_dirs[${PMI_RANK}]}\"\n')
        f.write('cd ${my_run_dir}\n')
        f.write('chm od u+rwx *\n')
        f.write('./run.run\n')
    # set permissions for output file
    os.chmod(output_file, 0o764)    


# can iterate over run paths
def consolidate_runs(run_path):
    '''Given a directory for a single run,
    finds the ALTI files and consolidates them into a .npz file.'''

    # assumes that the output files are in run_path/out
    dir_path = os.path.dirname(run_path)
    output_dir = os.path.join(dir_path, 'out')

    altis = glob.glob(output_dir + '/ALTI*')
    
    ndarry_altis = []
    for alti in sorted(altis):
        ndarry_altis.append(np.loadtxt(alti))
    ndarray_altis = np.stack(ndarray_altis)
    save_file = os.path.join(output_dir, 'height_maps.npz')
    np.savez(save_file, height_maps=ndarray_altis)



  
    

if __name__ == '__main__':

    # modifications for the big runs
    p_mods = {
        'stop after' : '200_t0',
        'output interval' : '100_t0',
        'alti only' : True,
        'rescallocation': '../../../scripts',
        'real_data_location' : '../../../scripts/real_data'
    }

    
    parameters_2 = {**parameters_1, **p_mods}

    # will need to run the .sbatch file created
    rand_initial_cells, modded_cells, run_paths = set_up_data_run(parameters_1, parameters_2,
                                                                  2, ['space_invader', 'sine'], 1, 2,
                                                                  '../../tdi', '../../tdm', '../../tdr',)
        
    

    

    
