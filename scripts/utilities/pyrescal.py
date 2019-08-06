__author__ = '''Gian-Carlo DeFazio'''

__doc__ = r'''A python interface to rescal. A pyrescal process sets up the files needed,
then calls rescal. The data from rescal that would normally be written to .csp files is 
piped to the parent for processing.'''

import os
import subprocess
import struct

import cellspace
import rescal_utilities




# list of rescal parameters that are paths to files/directories that
# already exist
# constructor will try to verify that they exist 





class PyRescal:
    '''Runs a rescal process and gets its data.'''

    def __init__(self, parameters, experiment_directory, rescal_root=None,
                 keep_height_maps=True, keep_ffts=True, relative_parameter_paths=True):
        '''Set up all the file locations and lists to hold the data files.'''
        
        self.parameters = parameters

        # try to find RESCAL_SNOW_ROOT
        if rescal_root is None:
            self.rescal_root = os.path.expanduser(os.environ['RESCAL_SNOW_ROOT'])
        else:
            self.rescal_root = os.path.expanduser(rescal_root)
        self.rescal_root = os.path.abspath(self.rescal_root)

        
        
        # deal with experiment_directory
        # TODO deal with directory already existing
        if not os.path.isabs(experiment_directory):
            self.experiment_directory = os.path.join(self.rescal_root, experiment_directory)
        else:
            # TODO, this could break everything if user gives bad path
            self.experiment_directory = os.fspath(experiment_directory)
        self.experiment_directory = os.path.abspath(self.experiment_directory)
            
        # set up paths to other files/directories of interest
        self.scripts_directory = os.path.join(self.rescal_root, 'scripts')
        self.real_data_directory = os.path.join(self.scripts_directory, 'real_data')
        self.rescal_executable = os.path.join(self.scripts_directory, 'rescal')
        self.genesis_executable = os.path.join(self.scripts_directory, 'genesis')
        self.ui = os.path.join(self.scripts_directory, 'rescal_ui.xml')
        self.par = os.path.join(self.experiment_directory, 'run.par')

        # parameters needs to contain directory paths for .par file generation to work
        paths_to_add = {'rescallocation' : self.scripts_directory,
                        'real_data_location' : self.real_data_directory}
        
        
        
        # a list to hold all the data
        self.keep_height_maps = keep_height_maps
        self.keep_ffts = keep_ffts
        
        self.height_maps = []
        self.ffts = []

    def __str__(self):
        '''A very generic __str__ that just gets everythiing that not a built-in
        or a method.'''
        # get the attributes that matter
        # from stack-exchange
        # https://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
        attrs = [a for a in dir(self) if not a.startswith('__') and not callable(getattr(self,a))]
        s = []
        for attr in attrs:
            s.append(attr + ': ' + str(self.__dict__[attr]))
        return '\n'.join(s)
    
    
    def create_experiment_directory(self):
        ''' '''
        # make the experiment directory
        if not os.path.exists(self.experiment_directory):
            os.mkdir(self.experiment_directory)
            
        # create a run with the given parameters
        design = rescal_utilities.Design_a_run()
        design.set_parameters(self.parameters)

        # find file paths that are in the real data folder and fix them
        # using self.real_data_directory 
        design.parameters.set_real_data_path(self.real_data_directory)
        
        # write the par file
        design.parameters.write(self.par)

        # now get the cmd line args for rescal
        cmd_args = design.run_script.rescal_call_args()

        # put rescal and the par file in, after nice if nice exists
        if cmd_args[0].strip() == 'nice':
            cmd_args = cmd_args[0] + [self.rescal_executable, self.par] + cmd_args[1:]
        else:
            cmd_args = [self.rescal_executable, self.par] + cmd_args
        self.rescal_args = cmd_args
        
        # do genesis if needed or link a premade cellspace

        
            
    def receive_data(self):
        '''Get the data size, which will be in a 4 byte integer.
        Get the data itself and store it '''
        
        data_size_bytes = os.read(self.r, 4)
        data_size = struct.unpack('i', data_size_bytes)[0]

        self.data = os.read(self.r, data_size)
        


    def process_data(self):
        c = cellSpace(self.data)
        # get the actual heightmap and ffts arrays out
        if self.keep_height_maps:
            self.height_maps.append(c.height_map.heigt_map)
        if self.keep_ffts:
            self.ffts.append(c.height_map.fft_blur)
        


    
    

    def run_simulation(self):

        # setup pipe to rescal
        r,w = os.pipe()

        self.r = r
        
        # allow rescal to inherit the write side of the pipe
        os.set_inheritable(w, True)

        # add the pipe to self.rescal args
        self.rescal_args = self.rescal_args + ' -data_pipe ' + str(w)

        # close the write pipe on this side
        os.close(w)
        
        # start a rescal process that can send data back
        rescal = subprocess.Popen(self.rescal_args, pass_fds=([w]))

        # check if rescal is still running
        while rescal.poll() is None:

            # now get the data back and process it
            receive_data()
            process_data()

        
        
       


if __name__ == '__main__':

    parameters_par_1 = {
        'Model':  'SNO',
        'Output_directory': './out',
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

            
    p = PyRescal(parameters_1, 'pr_test')
