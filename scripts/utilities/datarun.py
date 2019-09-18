__author__ = """Gian-Carlo DeFazio"""

__doc__ = r"""A python interface to Rescal-snow. A DataRun object takes in the 
parameters and meta-parameters required to run Rescal-snow. The DataRun can receive and process
the output of Rescal-snow while Rescal-snow is running.
To run Rescal-snow using a DataRun instance, the environment variable RESCAL_SNOW_ROOT
should be defined and be the path of a Rescal-snow installation. Also, a directory
for the output should be created. The default is RESCAL_SNOW_ROOT/data_runs.
"""

import os
import sys
import subprocess
import struct
import numpy as np

import cellspace
import rescal_utilities


def find_rescal_root(filename=None):
    """
    Tries to find RESCAL_SNOW_ROOT in the environment
    and make it into an absolute path.
    """
    if filename is None:
        if 'RESCAL_SNOW_ROOT' in os.environ:
            return os.path.expanduser(os.environ['RESCAL_SNOW_ROOT'])
        else:
            sys.stderr.write('Rescal root not found. Set environment variable RESCAL_SNOW_ROOT ' \
                             'to be the path to the top directory of a ReSCAL install or specify another path ' \
                             'for \'rescal_root\'.\n')
            sys.exit(-1)
    else:
        return os.path.expanduser(filename)




class DataRun:
    """
    Run ReSCAL with the option to receive and process the output 
    while ReSCAL is running. Save the processed output to a file.
    """

    def __init__(self, parameters, experiment_directory,
                 experiment_parent_directory='data_runs', rescal_root=None,
                 keep_height_maps=True, keep_ffts=True):
        """
        Collect and process the parameters for running ReSCAL and 
        processing its output. Set up the data structures to store ReSCAL outout.
        Determine where the executables and configuration files are or should be.
        For this initialization to work, the environment variable RESCAL_SNOW_ROOT 
        should exist, and its value should be the path to the top directory of an 
        installed ReSCAL. A 'data_runs' directory, or some other name specified by 
        experiment_parent_directory should already exist at RESCAL_SNOW_ROOT.
        
        arguments:

            parameters -- a dictionary containing the parameters for rescal 
            and meta-parameters for running ReSCAL.

            experiment_directory -- the directory that will hold 
            output data, and may have input files depending
            on the configuration.

        keyword arguments:

            experiment_parent_directory -- the parent directory of experiment_directory.
            the purpose of using this directory is to keep
            experiment directories from accumulating in the
            RESCAL_SNOW_ROOT directory. The directory must already exist.
            (default 'data_runs')

            rescal_root -- the path to a ReSCAL install that will be used 
            for the data run. If the value is None and the enviroment variable
            RESCAL_SNOW_ROOT exist, then the value of RESCAL_SNOW_ROOT will be 
            used. Otherwise, a user specified path will be used. If the value is 
            None and RESCAL_SNOW_ROOT is not defined, initialization will fail.
            (default None)

            keep_height_maps -- will keep and write out the height map data
            created by processing the cellspace outputs from ReSCAL.
            (default True)

            keep_ffts -- will keep and write out the fft blur data created by 
            processing cell space outputs from ReSCAL.
            (default True)
        """

        # set self.rescal_root using the environment or rescal_root argument
        # should be the parent directory of the ReSCAL src and scripts directories
        self.rescal_root = find_rescal_root(rescal_root)
        
        # find the experiment_parent directory, this needs to already exist
        # attempting to create will cause race conditions if doing multiple simultaneous runs
        if not os.path.isabs(experiment_parent_directory):
            self.experiment_parent_directory = os.path.join(self.rescal_root,
                                                            experiment_parent_directory)
        else:
            self.experiment_parent_directory = os.fspath(experiment_parent_directory)

        # deal with experiment_directory
        if not os.path.isabs(experiment_directory):
            self.experiment_directory = os.path.join(self.experiment_parent_directory,
                                                     experiment_directory)
        else:
            self.experiment_directory = os.fspath(experiment_directory)

        # deal with output directory which should be made by ReSCAL            
        ouput_directory = 'out'
        if 'Output_directory' in parameters.keys():
            output_directory = parameters['Output_directory']
        if not os.path.isabs(output_directory):
            self.output_directory = os.path.join(self.experiment_directory, output_directory)
        else:
            self.output_directory = os.fspath(output_directory)
            
        # set up paths to other files/directories of interest
        self.rescal_src_directory = os.path.join(self.rescal_root, 'src')
        self.scripts_directory = os.path.join(self.rescal_root, 'scripts')
        self.build_directory = os.path.join(self.rescal_root, 'build')
        self.real_data_directory = os.path.join(self.scripts_directory, 'real_data')
        self.rescal_executable = os.path.join(self.build_directory, 'rescal')
        self.genesis_executable = os.path.join(self.build_directory, 'genesis')
        self.ui_xml = os.path.join(self.rescal_src_directory, 'rescal-ui.xml')
        self.par = os.path.join(self.experiment_directory, 'run.par')
        self.meta_data = os.path.join(self.experiment_directory, 'meta_data')
        
        
        # parameters needs to contain directory paths for .par file generation to work
        paths_to_add = {'rescallocation' : self.build_directory,
                        'real_data_location' : self.real_data_directory}

        self.parameters = {**parameters, **paths_to_add}
        
        # for dealing with the data from ReSCAL
        self.keep_height_maps = keep_height_maps
        self.keep_ffts = keep_ffts
        self.height_maps = []
        self.ffts = []
        self.height_maps_path = os.path.join(self.output_directory, 'height_maps.npz')
        self.ffts_path = os.path.join(self.output_directory, 'ffts.npz')



    def check_paths(self):
        """
        Verify that the files and directories that should already exist 
        actually do exist.
        """

        dirs =  [self.rescal_root, self.scripts_directory,
                 self.real_data_directory, self.experiment_parent_directory,
                 self.rescal_src_directory]
        files = [self.rescal_executable, self.genesis_executable, self.ui_xml]
        for d in dirs:
            if not os.path.isdir(d):
                sys.stderr.write('ERROR: directory ' + d + ' not found.\n')
                sys.exit(-1)
        for f in files:
            if not os.path.isfile(f):
                sys.stderr.write('ERROR: file ' + f + ' not found.\n')
                sys.exit(-1)
                
        

    def __str__(self):
        """
        A very generic __str__ that just gets everythiing that's not a built-in
        or a method.
        """
        # get the attributes that matter
        # from stack-exchange
        # https://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
        attrs = [a for a in dir(self) if not a.startswith('__') and not callable(getattr(self,a))]
        s = []
        for attr in attrs:
            s.append(attr + ': ' + str(self.__dict__[attr]))
        return '\n'.join(s)
    
    
    def setup(self):
        """
        Create the experiment directory and the .par file. Create the .csp input 
        using genesis if neccesary. Deal with the file path to real_data.
        Get the args for the call to rescal.
        """

        # some files and directories should alreay exist
        self.check_paths()
        
        # make the experiment directory
        if not os.path.exists(self.experiment_directory):
            os.mkdir(self.experiment_directory)

        # deal with .csp files
        # if 'premade_csp' overwrite 'Csp_file'
        if 'premade_csp' in self.parameters.keys():
            # put the .csp name into the .par file as Csp_file
            design.set_parameters({'Csp_file' : self.parameters['premade_csp']})
        else:
            # if 'Csp_file' just a file name, put it into self.experiment_directory
            if 'Csp_file' in self.parameters.keys():
                # if it's an absolute path, leave it alone, otherwise
                # take the file name and put it into the experiment directory
                if not os.path.isabs(self.parameters['Csp_file']):
                    self.parameters['Csp_file'] = os.path.join(self.experiment_directory,
                                                               os.path.basename(self.parameters['Csp_file']))
            else:
                # give 'Csp_file' a generic name
                self.parameters['Csp_file'] = os.path.join(self.experiment_directory, 'Cells.csp')

        # set output directory
        # TODO allow users to specify their own output directory
        self.parameters['Output_directory'] = os.path.join(self.experiment_directory, 'out')
                
        # create a run with the given parameters
        design = rescal_utilities.DesignRun()
        design.set_parameters(self.parameters)

        # find file paths that are in the real data folder and fix them
        # using self.real_data_directory 
        design.parameters.set_real_data_path(self.real_data_directory)
            
        # write the par file
        design.parameters.write(self.par)

        # run genesis if neccesary
        if 'premade_csp' not in self.parameters.keys():
            # run genesis and put the output where specified by 'Csp_file' parameter
            p = subprocess.check_call([self.genesis_executable, '-f',  self.par, '-s', str(2000)],
                                 stdout=subprocess.DEVNULL)
        

        
        # now get the cmd line args for rescal
        cmd_args = design.run_script.rescal_call_args()

        # put rescal and the par file in, after nice if nice exists
        if cmd_args[0].strip() == 'nice':
            cmd_args = cmd_args[0] + [self.rescal_executable, self.par] + cmd_args[1:]
        else:
            cmd_args = [self.rescal_executable, self.par] + cmd_args
        self.rescal_args = cmd_args
        
            
    def receive_process_data(self):
        """
        Get the data size, which will be in a 4 byte integer.
        Get the data itself. Create a cellspace.CellSpace to
        process the data and keep the height_maps and ffts 
        based on flags.
        """

        # either block because there's nothing to read,
        # read the next .csp,
        # or see EOF and return
        data_size_bytes = os.read(self.r, 4)
        if not data_size_bytes:
            return
        data_size = struct.unpack('i', data_size_bytes)[0]

        self.data = os.read(self.r, data_size)
        
        c = cellspace.CellSpace(self.data)
        # get the actual heightmap and ffts arrays out
        if self.keep_height_maps:
            self.height_maps.append(c.height_map.height_map)
        if self.keep_ffts:
            self.ffts.append(c.height_map.fft_blur)

        
    def write_meta_data(self):
        """
        Write the rescal parameters and execuation call to a file
        in the experiment directory.
        """
        with open(self.meta_data, 'w+') as f:
            f.write('parameters = ' + str(self.parameters) + '\n')
            f.write('exec = ' + str(self.rescal_args) + '\n')


    def run_simulation(self):
        """
        Run ReSCAL and get the cellspace data
        that would otherwise be written to a file. Process and save the data
        to files and write out the meta-data. A new file descriptor is used,
        as opposed to just redirecting the stdot of ReSCAL, because the stdout
        of ReSCAL is full of text output that this function doesn't deal with.
        """

        # setup pipe to rescal
        r,w = os.pipe()

        # save the read side of the pipe
        self.r = r
        
        # allow rescal to inherit the write side of the pipe
        os.set_inheritable(w, True)

        # add the pipe to self.rescal args so rescal knows it's open
        self.rescal_args = self.rescal_args + ['-data_pipe',  str(w)]

        # start a rescal process that can send data back
        # data is send back one .csp at a time through the pipe
        rescal = subprocess.Popen(self.rescal_args,
                                  pass_fds=([w]),
                                  cwd=self.experiment_directory,
                                  stdout=subprocess.DEVNULL)

        # close the write pipe on this side
        os.close(w)

        # check if rescal is still running
        while rescal.poll() is None:

            # now get the data back and process it
            self.receive_process_data()

        # ensure that an output_directory exists, just in case ReSCAL didn't make it
        if not os.path.isdir(self.output_directory):
            os.mkdir(self.ouput_directory)
        # write the data
        if self.keep_height_maps:
            height_maps = np.array(self.height_maps)
            np.savez_compressed(self.height_maps_path, height_maps=height_maps)
        if self.keep_ffts:
            ffts = np.array(self.ffts) 
            np.savez_compressed(self.ffts_path, ffts=ffts)

        # write out the meta_data
        self.write_meta_data()

        

    def run_without_piping(self):
        """Do a rescal run that saves data to files like normal."""
        rescal = subprocess.Popen(self.rescal_args,
                                  cwd=self.experiment_directory,
                                  stdout=subprocess.DEVNULL)
        rescal.wait()

