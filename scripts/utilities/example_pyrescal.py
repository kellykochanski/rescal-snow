#!/usr/bin/env python3

__author__ = """Gian-Carlo DeFazio"""
__doc__ = """This example is meant to test functionality of datarun.py
on slurm. This script can be called indirectly using sbatch.
More information is provided in the datarun tutorial in docs."""


import datarun
import os

parameters_par_1 = {
    'Model':  'SNO',
    'Output_directory': 'out',
    'Csp_file': 'DUN.csp',  # could just put the premade .csp in here
    'Csp_template': 'SNOWFALL(4)',
    'parfile': 'run.par',
    'Boundary':  'OPEN',
    'Time': 0.0,
    'H': 80,
    'L': 400,
    'D': 200,
    'Centering_delay': 0,
    'Phys_prop_file': 'real_data/sealevel_snow.prop', # need to use rescal home
    'Qsat_file': 'real_data/PDF.data', # need to use rescal home
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
}

# command line arguments for rescal which will be in the .run file
parameters_run_1 = {
    'stop after' : '200_t0',
    'output interval' : '50_t0',
    'png interval' : False,
    'quit' : True,
    'random seed' : 6,
    'usage info' : False,
    'show params' : False,
    'info interval' : False,
}

parameters_1 = {**parameters_par_1, **parameters_run_1}



if __name__ == '__main__':

    # will do a bunch of rescal runs with different random seeds

    # find the pmi_rank
    # (This is the unique identifier of each run on a parallel cluster.
    # It may have different names on different architectures.)
    if bool(os.environ.get('PMI_RANK', None)): 
    	pmi_rank = os.environ['PMI_RANK']
    else:
        print("Warning: process identifier PMI_RANK not found. If this run is parallel, collisions may occur.")
        pmi_rank = 0
        # safe to make sure data_runs directory exists
        rescal_root = os.environ['RESCAL_SNOW_ROOT']
        data_runs_dir = os.path.join(rescal_root, 'data_runs')
        if not os.path.isdir(data_runs_dir):
            os.mkdir(data_runs_dir)
        
    parameters_1['random seed'] = int(pmi_rank) + 1234
    output_dir = 'exp' + str(pmi_rank)
    
    dr = datarun.DataRun(parameters_1, output_dir)
    dr.setup()
    dr.run_simulation()
    #dr.run_without_piping()
