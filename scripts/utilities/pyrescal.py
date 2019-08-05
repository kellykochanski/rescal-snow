__author__ = '''Gian-Carlo DeFazio'''

__doc__ = r'''A python interface to rescal. A pyrescal process sets up the files needed,
then calls rescal. The data from rescal that would normally be written to .csp files is 
piped to the parent for processing.'''

import os
import subprocess
import struct


class PyRescal:
    '''Runs a rescal process and gets its data.'''


    def __init__(self, parameters, processing_options, rescal_root=None):
        self.parameters = parameters
        self.processing_options = processing_options

        # try to find RESCAL_SNOW_ROOT
        if rescal_root is None:
            self.rescal_root = os.environ['RESCAL_SNOW_ROOT']
        else:
            self.rescal_root = rescal_root

        
            
    def receive_data(self):
        '''Get the data size, which will be in a 4 byte integer.
        Get the data itself and store it '''

        data_size_bytes = os.read(self.r, 4)
        data_size = struct.unpack('i', data_size_bytes)[0]

        self.data = os.read(self.r, data_size)
        

        


    
    

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
        rescal = subprocess.Popen([self.rescal_executable] + [self.rescal_args], pass_fds=([w]))

        # now get the data back and process it

    

    


    
