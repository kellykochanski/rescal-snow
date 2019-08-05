import cellspace
# another exploration example varying the random seed
#from rescal_utilities import *
import time_lapse_compare as tlc
import time_lapse_visualize as tlv
import param_space_exploration_utilities as pseu
import numpy as np


# Parameters held constant for all simulations
parameters = {
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
    #'Lambda_I': 0.01 # original value was causing whole space to fill with snow
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
    'stop after' : '100_t0',
    'output interval' : '100_t0',
    'png interval' : False,
    'quit' : True,
    #'premade_csp' : '/g/g13/defazio1/rescal-snow/scripts/SNO2500_t0_rng_4.csp'
}

var_params = [['random seed', [7,8]]]

def modify_csp_and_restart():

    
    files = pseu.make_run_directories(parameters, var_params, 'exp_csp_s', 'run', 'run')
    pseu.run_rescals(files)
    paths = pseu.get_files_to_process('../../exp_csp_s', '/*.csp.gz')
    # paths should only contain one file
    #print(paths)
    c = cellspace.CellSpace(paths[0][1])
 #   c.draw_height_map()

    # edit the csp file


    #barry_face = np.read('barry_b_w.npy').astype(np.uint8) * 6
    
    pic = cellspace.make_gaussian(8, 30,30,11.0)
    input_paths = c.multiple_random_pics(2, pic, 'hello_s')
    #input_paths = c.multiple_random_pics(6, barry_face, 'hello')
    
    #### make multiple edit and write them all out
    #### make a list of the paths
    #### then do a sun with premade_csp as a varparam


    #c.draw_height_map()

    # start a new run where a modified .csp is used as the start point

    #print(input_paths)
    parameters['stop after'] = '200_t0' 

    
    new_var_params = var_params + [['premade_csp',input_paths]]

    #print(new_var_params)
    
    files = pseu.make_run_directories(parameters, new_var_params, 'exp_csp_s1', 'run', 'run')
    pseu.run_rescals(files)
    paths = pseu.get_files_to_process('../../exp_csp_s1', '/*.csp.gz')
    # paths should only contain one file
    #print(paths)
    # for p in paths:
    #     for q in p:
    #         c = cellspace.CellSpace(q)
    #         c.draw_height_map()
    

def just_run_it():
    files = pseu.make_run_directories(parameters, var_params, 'exp_xxx1', 'run', 'run')
    pseu.run_rescals(files)
    paths = pseu.get_files_to_process('../../exp_xxx1', '/*.csp.gz')
    return paths
    

if __name__ == '__main__':
    modify_csp_and_restart()
    #just_run_it()

    
