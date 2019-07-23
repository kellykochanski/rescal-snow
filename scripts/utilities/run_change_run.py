
import numpy as np
import random
import os
import rescal_utilities
import cellspace

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
    'stop after' : '200_t0',
    'output interval' : '200_t0',
    'png interval' : False,
    'quit' : True,
    'random seed' : 6,
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
            file_prefix = top_dir + '/' + exp_name + str(i) 
            mod_paths_inner.append(c.random_mods(file_prefix, mod_type,
                                                 num_mods=num_mods, temp_mod=True))
        mod_paths.append(mod_paths_inner)
    return mod_paths


# # TODO change function name
# # fixed_and_var_params  the fixed and variable parameters for the first set of runs
# # fixed_and_var_mods    the fixed and variable parameters for the second set of runs
# def mod_res(params_1, params_2, top_dir_name, var_params_1=None, var_params_2=None):

#     # will do a single run with seed set to 7 by default
#     if var_params_1 == None:
#         var_params_1 = [['random seed' [7]]]
#     if var_params_2 == None:
#         var_params_2 = [['random seed' [7]]]

#     # TODO maybe call it make run scripts
#     # create all the files and directories needed for the runs
#     run_files = param_space_exploration_utilities.make_run_directories(params_1,
#                                                                        var_params_1,
#                                                                        top_dir_name,
#                                                                        'run',
#                                                                        'run')
#     # run all the .run files, in parallel I think
#     param_space_exploration_utilities.run_rescals(run_files)

#     # get the paths to all the .csp file that were just creates
#     paths = param_space_exploration_utilities.get_files_to_process(top_dir_name,
#                                                                    cellspace.path_glob,
#                                                                    cellspace.exclude_globs)

#     return paths




# def modify_csp_and_restart():

    
#     files = pseu.make_run_directories(parameters, var_params, 'exp_csp', 'run', 'run')
#     pseu.run_rescals(files)
#     paths = pseu.get_files_to_process('../../exp_csp', '/*.csp.gz')
#     # paths should only contain one file
#     #print(paths)
#     c = cellspace.CellSpace(paths[0][1])
#  #   c.draw_height_map()

#     # edit the csp file


#     #barry_face = np.read('barry_b_w.npy').astype(np.uint8) * 6
    
#     pic = cellspace.make_gaussian(8, 30,30,11.0)
#     input_paths = c.multiple_random_pics(8, pic, 'hello')
#     #input_paths = c.multiple_random_pics(6, barry_face, 'hello')
    
#     #### make multiple edit and write them all out
#     #### make a list of the paths
#     #### then do a sun with premade_csp as a varparam


#     #c.draw_height_map()

#     # start a new run where a modified .csp is used as the start point

#     #print(input_paths)
#     parameters['stop after'] = '3000_t0'

    
#     new_var_params = var_params + [['premade_csp',input_paths]]

#     #print(new_var_params)
    
#     files = pseu.make_run_directories(parameters, new_var_params, 'exp_csp1', 'run', 'run')
#     pseu.run_rescals(files)
#     paths = pseu.get_files_to_process('../../exp_csp1', '/*.csp.gz')
#     # paths should only contain one file
#     #print(paths)
#     # for p in paths:
#     #     for q in p:
#     #         c = cellspace.CellSpace(q)
#     #         c.draw_height_map()
    

# def just_run_it():
#     files = pseu.make_run_directories(parameters, var_params, 'exp_xxx1', 'run', 'run')
#     pseu.run_rescals(files)
#     paths = pseu.get_files_to_process('../../exp_xxx1', '/*.csp.gz')
#     return paths
    

if __name__ == '__main__':
#    randos = random_initial_states(2, parameters_1, '../../rr')


    saved_randos = ['/g/g13/defazio1/summer_2019/rescal-snow/rr/random_seed-518542/SNO00001_t0.csp.gz',
                    '/g/g13/defazio1/summer_2019/rescal-snow/rr/random_seed-561509/SNO00001_t0.csp.gz']
    
 #   print(randos)

    
    modded_randos = modify_outputs('../../ss',
                                   saved_randos,
                                   ['space_invader'],
                                   2)
    print(modded_randos)

    saved_modded_randos = [[['/g/g13/defazio1/summer_2019/rescal-snow/ss/exp0--invader--0.csp',
                             '/g/g13/defazio1/summer_2019/rescal-snow/ss/exp0--invader--1.csp'],
                            ['/g/g13/defazio1/summer_2019/rescal-snow/ss/exp0--gaussian--0.csp',
                             '/g/g13/defazio1/summer_2019/rescal-snow/ss/exp0--gaussian--1.csp']],
                           
                           [['/g/g13/defazio1/summer_2019/rescal-snow/ss/exp1--invader--0.csp',
                             '/g/g13/defazio1/summer_2019/rescal-snow/ss/exp1--invader--1.csp'],
                            ['/g/g13/defazio1/summer_2019/rescal-snow/ss/exp1--gaussian--0.csp',
                             '/g/g13/defazio1/summer_2019/rescal-snow/ss/exp1--gaussian--1.csp']]]
    #modify_csp_and_restart()
    #just_run_it()

    
