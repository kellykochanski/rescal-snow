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
        



    
pickled_experiment = 'pickled_experiment'



# hold the meta data for the data run
class MetaData:
    def __init__(self):
        self.intial_csps_created = False
        self.control_run_files_created = False        
        self.modded_csp_created = False
        self.moded_run_files_created = False
        self.main_experiment_started = False
        self.data_consolidated = False
        self.jobs = []

        # experiment data
        self.have_experiment_data = False
        self.parameters_initial = None
        self.parameters_after_mod = None
        self.num_rand_initial = None
        self.mod_types = None
        self.num_mods = None
        self.num_runs = None
        self.chunk_size = None


    def __str__(self):
        # get the attributes that matter
        # from stack-exchange
        # https://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
        attrs = [a for a in dir(self) if not a.startswith('__') and not callable(getattr(self,a))]
        s = []
        for attr in attrs:
            s.append(attr + ': ' + str(self.__dict__[attr]))
        return '\n'.join(s)
            


        
    
# holds the eperiment data for the data run
class ExperimentData:

    def __init__(self, parameters_initial=None, parameters_after_mod=None,
                 num_rand_initial=None, mod_types=None, num_mods=None,
                 num_runs=None, chunk_size=None):

        self.parameters_initial = parameters_initial
        self.parameters_after_mod = parameters_after_mod
        self.num_rand_initial = num_rand_initial 
        self.mod_types = mod_types
        self.num_mods = num_mods
        self.num_runs = num_runs
        self.chunk_size = chunk_size

                       
    # write pickled version of self to file
    def write(self, filename):
        with open(filename, 'wb+') as f:
            pickle.dump(self, f)

    
        
        
# a class that represents a rescal data run
# uses directory in the top level directory of a ReSCAL install
class DataRun:

    def __init__(self, experiment_directory, rescal_root=None, experiment_data=None, overwrite=False):
        # get the location of rescal
        if rescal_root is None:
            self.rescal_root = os.path.expanduser(os.environ['RESCAL_SNOW_ROOT'])
        else:
            self.rescal_root = os.path.expanduser(rescal_root)


        # deal with experiment_directory
        if not os.path.isabs(experiment_directory):
            self.experiment_directory = os.path.join(self.rescal_root, experiment_directory)
        else:
            # TODO, this could break everything if user gives bad path
            self.experiment_directory = os.fspath(experiment_directory)

        self.meta_data_path = os.path.join(self.experiment_directory, '.meta_data')
            
        # see if experiment_directory exists and is a directory
        # TODO deal with permissions
        if os.path.isdir(self.experiment_directory):
            # look for the meta_data file
            if os.path.isfile(self.meta_data_path):
                # if found, use it to initialize
                with open(self.meta_data_path, 'rb') as f: 
                    self.meta_data = pickle.load(f)
        # failure, won't everwrite file
        elif os.path.isfile(self.experiment_directory):
            print('A file exists at {ed}'.format(self.experiment_directory))
            print('Cannnot create experiment directory')

        # nothing exists at self.experiment_directory, so make and setup new directory
        else:
            # create the directory
            os.mkdir(self.experiment_directory)
            with open(self.meta_data_path, 'wb+') as f:
                self.meta_data = MetaData()
                pickle.dump(self.meta_data, f)

        self.setup_experiment(experiment_data)
            
        
    # write meta_data out to file
    def write_meta_data(self):
        with open(self.meta_data_path, 'wb') as f:
            pickle.dump(self.meta_data, f)
        

    # copies the experiment data into the meta data
    # from a ExperimentData object
    def copy_experiment_data(self, experiment_data):
        self.meta_data.parameters_initial = experiment_data.parameters_initial
        self.meta_data.parameters_after_mod = experiment_data.parameters_after_mod
        self.meta_data.num_rand_initial = experiment_data.num_rand_initial 
        self.meta_data.mod_types = experiment_data.mod_types
        self.meta_data.num_mods = experiment_data.num_mods
        self.meta_data.num_runs = experiment_data.num_runs
        self.meta_data.chunk_size = experiment_data.chunk_size
        self.meta_data.have_experiment_data = True
        # and save it
        self.write_meta_data()

            
    # get experiment parameters from a file path or
    # passed as an ExperimentData object
    def setup_experiment(self, experiment_data):

        if experiment_data is None and self.meta_data.have_experiment_data == False:
            print('WARNING: DataRun has no experiment data.')
        elif experiment_data is not None and self.meta_data.have_experiment_data == True:
            print('WARNING: DataRun experiment_data alreay exists.')
            print('cannot overwrite experiment data at instantiation time')
            print('use replace_experiment_data method')
        elif experiment_data is not None and self.meta_data.have_experiment_data == False:
            # determine if experiment_data is file name of ExperimentData
            if isinstance(experiment_data, ExperimentData):
                self.copy_experiment_data(experiment_data)
            elif isinstance(experiment_data, str):    
                # try to read in pickled ExperimentData
                with open(experiment_data, 'rb') as f: 
                    self.copy_experiment_data(pickle.load(f))
            else:
                print('WARNING: experiment_data must be of type str of ExperiementData')
                print('DataRun has no experiment_data.')
            #experiment_data is None and self.meta_data.have_experiment_data == True:
        else:
            # nothing to do, typical case when reloading from file
            pass

        
    # TODO maybe do more than just the meta_data
    def initialize_experiment_directory(self):
        self.create_meta_data()
        

    # create a meta_data object and file
    def create_meta_data(self):
        with open(self.meta_data_path, 'wb') as f:
            self.meta_data = MetaData()
            pickle.dump(self.meta_data, f)
    
            
    def initialize_from_existing_experiement(self):
        self.meta_data = pickle.load(self.meta_data_path)

        
        
    # setup directories, this includes running rescal
    # to randomly initialize
    def setup(self):
        pass


    # makes randoms initial states
    # sets up runs directories to make initial states
    # runs the simulations to get the initial states
    # saves the file paths to the .csp files created
    def random_initial_states(self):

        random_intials_paths = os.path.join(self.experiment_directory, 'random_initials_paths'))
        
        if self.meta_data.initial_csps_created = True:
            # read in paths from file
            with open(random_intials_paths, 'r') as f:
                self.random_initial_paths = f.read.splitlines()
        
        else:
            # will create directory called random initial to
            # do the intial random runs
            top_dir = os.path.join(self.experiment_directory, '.random_initials')

            # create num_states random seeds
            seed_numbers = random.sample(range(1,1000000), self.meta_data.num_states)
            seeds = [['random seed', seed_numbers]]

            run_files = rescal_utilities.make_run_directories(parameters, seeds, top_dir, run_header, run_name)

            rescal_utilities.run_rescals(run_files)

            paths = rescal_utilities.get_files_to_process(top_dir, cellspace.path_glob, cellspace.exclude_globs)
            paths_last = [x[-1] for x in paths]

            # save the paths
            with open(random_intials_paths, 'w+') as f:
                f.write('\n'.join(paths_last))

            self.random_initial_paths = paths_last

            self.meta_data.initial_csps_created = True





    
    












    
