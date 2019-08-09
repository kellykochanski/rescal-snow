__doc__ = 'the cellspace module'
__author__ = 'Gian-Carlo DeFazio'
__date__ = 'July 31 2019'


import matplotlib.pyplot as plt
import matplotlib.colors as  colors

import pydoc

import struct
import numpy as np
import random
import gzip
import os
import sys
import argparse
import scipy.ndimage
import random

import heightmap



def make_gaussian(h, d_padding, l_padding, sigma):
    '''creates a heightmap for a guassian
    scales the max height to h
    set all values to integers'''
    
    # does it the cheap way, makes a dirac delta
    # meaning a 1 in the middle and 0 otherwise
    # then filters it
    dd = np.zeros([2 * l_padding + 1, 2 * d_padding + 1])
    l, d = dd.shape
    dd[l//2, d//2] = 1
    g = scipy.ndimage.gaussian_filter(dd, sigma)
    g = g * (h / g[l//2,d//2]) + 0.01
    return (np.round_(g)).astype(np.uint8)


_cell_colors_list = ['khaki',         # grain
                    'turquoise',     # mobile_grain
                    'red',       # vegetated_grain
                    'white',        # air
                    'forestgreen',  # vegetation
                    'lightcyan',   # boundary
                    'dimgray',      # neutral
                    'darkorange',   # input_of_sand
                    'orchid',       # output_of_sand
                    'red',          # tunnel
                    'lightcoral',   # EAUT
                    'lemonchiffon', # colored_grain
                    ]


cell_types = {'grain' : 0,
              'mobile_grain' : 1,
              'vegetated_grain' : 2,
              'air' : 3,
              'vegetation' : 4,
              'boundary' : 5,
              'neutral' : 6,
              'input_of_sand' : 7,
              'ouput_of_sand' : 8,
              'tunnel' : 9,
              'EAUT' : 10, # not sure, used only graphically
              'colored_grain' : 11 # used graphically
}


mods = ['sine', 'noise', 'guassian hill', 'space_invader', 'clip']

def _ints_to_colors(color_list):
    '''takes a list of colors and creates a colormap
    that maps the index of the color to the color when drawing'''
    length = len(color_list)
    c_rgba = []
    for c in color_list:
        c_rgba.append(colors.to_rgba_array(c))
    c_rgba = np.reshape(c_rgba, (c_rgba.shape[0], c_rgba.shape[2]))
    return colors.ListedColormap(c_rgba, N=12)

def _draw_cross_section(window, axes):
    n = colors.Normalize(vmin=0, vmax=11)
    axes.pcolormesh(window, cmap=rescal_color_map, norm=n)
    axes.invert_yaxis()
    return


# TODO: I think my depth, height, and length are oriented differently than ReSCAL :-(
class CellSpace:
    '''A ReSCAL cell space representation. This class deals with the .csp files
    made by ReSCAL and genesis. The .csp files contain a header followed by the 3D space
    of cells as a 1D array.

    The parts of the header that are used are organized as follows:
   
    the first 20 bytes:
        [0..3]   B3s  the magic number which is \\212 or 0o212
                      in the source and 'CSP' <<138 'CSP'>>
        [4..5]   2s   not sure what this is, but it's <<'@1'>>
        [6]      c    can be <<'_'>> for little endian and <<'b'>> for big endian
        [7]      c    a newline <<'\\n'>>
        [8..11]  i    a chunk size, this doesn't appear to be used,
                       but it's 8 + sizeof(int) = 12  <<12>>
        12..15]  4s   a string that indicates the header size is next <<'HDSZ'>>
        [16..19] i    the header size,
                      which includes these first 20 bytes <<some_number>>

    after the first 20 bytes, the format is in chunks
    each chunk starts with a number (int) that says its size in bytes
    this size includes the 4 bytes for the chunk size number
    byte values are relative to the start of the chunk
    the chunks that are used are the model, size, and cell size chunks

    Model chunk:  contains the model type, which I think is a 3 byte string
        [0-3]   i  the size the chunk  <<some_number>> (probably 12)
        [4-7]   4s the chunk type identifier <<'MODL'>>
        [8-10]  3s the model name, for example <<'SNO'>>
        [11]    B  I think this is pad <<'\\0'>>

    Size chunk:  contains the dimensions of the simulations space 
    height (H), length (L), and depth (D) as 3 ints
        [0-3]     i  the size the chunk  <<some_number>> (probably 20)
        [4-7]     4s the chunk type identifier <<'SIZE'>>
        [8-11]    i  the height H <<some_number>>
        [12-15]   i  the length L <<some_number>>
        [16-19]   i  the depth  D <<some_number>>

    Cell Size chunk: contains the number of bytes for each cell
    which seems to default to 8
        [0-3]     i  the size the chunk  <<some_number>> (probably 12)
        [4-7]     4s the chunk type identifier <<'CELL'>>
        [8-11]    i  the number of bytes per cell <<some_number>> (probably 8)

    Sizes with borders chunk: contains the dimension with border cells
    and a flag to indicate if there are the dimensions to be used
        [0-3]     i i  the size the chunk  <<some_number>> (probably 24)
        [4-7]     4s the chunk type identifier <<'BORD'>>
        [8-11]    i  the height H <<some_number>>
        [12-15]   i  the length L <<some_number>>
        [16-19]   i  the depth  D <<some_number>>
        [20-23]   i  flag to say if thses values should be used. 1 if yes, 0 if no.'''

    
    b0_19 = 'B3s2scci4si'
    model_chunk = 'i4s3sB'
    size_chunk = 'i4siii'
    cell_size_chunk = 'i4si'
    borders_chunk = 'i4siiii'

    
    def __init__(self, filename, keep_original=True):
        '''read in a .csp or .csp.gz file and create a CellSpace instance
        to represent the file contents.

        arguments:
            filename -- the path to the input file

        keyword arguments:
            keep_orignal -- if True, retain the original cell data from the
                            input file. If False discard the orginal cell data.
                            (default True)

        usage example:
            c = CellSpace('path_to_file.csp')
        '''
            
        # set within self._read
        self.header_size = None
        self.model = None
        self.height = None
        self.length = None
        self.depth = None
        self.cell_size = None
        self.cell_data_type = None

        # gets original data from .csp file if
        # keep_original == True
        self.cells_original = None
        self.header_binary = None
        self.keep_original = keep_original

        # determine if data is in a file or is alreay bytes
        if isinstance(filename, bytes):
            self._read(bs=filename)
        else:            
            self.input_file = filename
            
            # remove .gz by default
            if self.input_file.endswith('.gz'):
                self.output_file = filename
            else:
                self.output_file = filename[:-3]
            self._read()

        # potentially used later
        self.temp_mod_cells = None

        # made during analysis
        self.height_map = None
        self.surface_map = None

        # also makes surface_map
        self.make_height_map()
        self.make_ceiling_map()



    def _get_cell_data_type(self):
        '''get correct data type depending on cell size'''
        data_type = None
        if self.cell_size == 8:
            data_type = np.uint64
        elif self.cell_size == 4:
            data_type = np.uint32
        elif self.cell_size == 2:
            data_type = np.uint16
        elif self.cell_size == 1:
            data_type = np.uint8
        else:
            print('not sure of cell data type')
        return data_type


    def _read(self, bs=None):
        '''read in the file contents and store the data in a ndarray.
        store the header both as python variables, but also in the original binary
        form from the input file.'''

        if bs is None:
            f = None
            if self.input_file.endswith('.gz'):
                f = gzip.open(self.input_file)
            else:
                f = open(self.input_file, 'rb')
            bs = bytearray(f.read())
            f.close()
        

        self._read_header(bs)

        # get correct data type depending on cell size
        self.cell_data_type = self._get_cell_data_type()

        # read the cells in ndarray
        cells = np.frombuffer(bs, dtype=self.cell_data_type, offset=self.header_size)

        # detemine how to reshape the cells, depending on whether or not border
        # cells are present
        if self.use_borders:
            cells = np.reshape(cells, (self.depth_borders,self.height_borders,self.length_borders))
            # now do some math
            pad_depth = (self.depth_borders - self.depth) // 2
            pad_height = (self.height_borders - self.height) // 2
            pad_length = (self.length_borders - self.length) // 2
            # remove the borders
            cells = cells[pad_depth:-pad_depth,pad_height:-pad_height,pad_length:-pad_length]
        else:
            cells = np.reshape(cells, (self.depth,self.height,self.length))
        # may choose to discard original cells to save space
        if self.keep_original:
            self.cells_original = cells
        self.cells = cells.astype(np.uint8)
        # make copy of the header for easy writing to file
        self.header_binary = bs[:self.header_size]
        return


    # TODO deal with 'TIME'
    # TODO should all the chunk size checks be 12? I think not.
    # sets values of self.header_size
    # self.model, self.depth, self.length, self.height
    # self.cell_size
    def _read_header(self, bs):
        '''read the .csp header and store the pertinent values'''

        # unpack first 20 bytes
        header_20 = struct.unpack(CellSpace.b0_19, bs[:20])
        # get the size of the header
        header_size = header_20[-1]
        self.header_size = header_size

        # grab the rest of the header
        header_20_to_end = bs[20:header_size]

        # find the other pertinent data
        # iterate through all remaining chunks and extract the importaing data
        chunk_start = 20
        # should be at least 8 more bytes for the size and type
        while chunk_start + 8 <= header_size:

            # get the next chuck size and type, which should always be 8 bytes
            chunk_size, chunk_type = struct.unpack('i4s', bs[chunk_start:chunk_start+8])

            if chunk_start + chunk_size > header_size:
                print('bad chunk size , will exceed header size')

            bytes_slice = slice(chunk_start, chunk_start + chunk_size)

            if chunk_type == b'MODL':
                if chunk_size < 12:
                    print('bad chunk size for Model')
                else:
                    # extract and store model name
                    _, _, model, _ = struct.unpack(CellSpace.model_chunk,bs[bytes_slice])
                    self.model = model.decode('utf-8')

            elif chunk_type == b'SIZE':
                if chunk_size < 12:
                    print('bad chunk size for SIZE')
                else:
                    # extract and store height, length, and depth
                    _, _, h, l, d = struct.unpack(CellSpace.size_chunk, bs[bytes_slice])
                    self.height, self.length, self.depth = h, l, d

            elif chunk_type == b'CELL':
                if chunk_size < 12:
                    print('bad chunk size for CELL')
                else:
                    _, _, cell_size = struct.unpack(CellSpace.cell_size_chunk, bs[bytes_slice])
                    self.cell_size = cell_size

            elif chunk_type == b'TIME':
                pass

            elif chunk_type == b'BORD':
                if chunk_size < 12:
                    print('bad chunk size for BORD')
                else:
                    # extract and store height, length, and depth
                    _, _, h_b, l_b, d_b, u_b = struct.unpack(CellSpace.borders_chunk, bs[bytes_slice])
                    self.height_borders, self.length_borders, self.depth_borders, self.use_borders = h_b, l_b, d_b, u_b
                
            else:
                print('unrecognized chunk type {n}'.format(n=chunk_type.decode('utf-8')))

            chunk_start = chunk_start + chunk_size
        return


    def _overwrite_lsbs(self, temp_mod=False):
        '''creates cells that are identical to those read in (self.cells_original) except
        for the lsb which is rewritten with the values in self.cells
        can use a temporary modification to write out'''

        # expand bytes to type read in from the input file
        bytes_expanded = None
        if temp_mod:
            bytes_expanded = self.temp_mod_cells.astype(self.cells_original.dtype)
        else:
            bytes_expanded = self.cells.astype(self.cells_original.dtype)
        mask = np.array([0xffffffffffffff00]).astype(self.cells_original.dtype)
        return (self.cells_original & mask) | bytes_expanded


    # writes a .csp file
    # if no filename is given, the input file name is used
    def write(self, filename=None, temp_mod=False, compressed=False):
        ''' write the CellSpace to a file in the .csp format

        keyword arguments:
            filename --   the name of the file to write to.
                          by default, the name of the input file used to initialize
                          the CellSpace will be used. If a .csp.gz file was used, the .gz 
                          is removed. The default value is None, but this will cause
                          self.output_file to be used. Note that is is possible to
                          change self.output_file after initialization.
                          (default None)

            temp_mod --   if True, temp_mod_cells are written. if False cells are written. 
                          temp_mod_cells are a copy of the orginal cells that may be modified
                          while the original cell remain unchanged. 
                          (default False)
        
            compressed -- if True, a gzipped file is written and '.gz.' is added to 
                          the filename if not already present. If False, the output is
                          uncompressed and the filename is unchanged, even if it ends
                          with '.gz'
                          (default False)

        return value: None

        usage examples:
            c = CellSpace('input.csp')

            no arguments will overwrite the input file using cells and the file
            will no be compressed. 'input.csp' is overwritten.
                c.write()

            compressed=True will write to a file that has the same
            name as the input but with '.gz' added. 'input.csp' is not overwritten
            and 'input.csp.gz' is written or overwritten.
                c.write(compressed=True)

            temp_mod=True will use the temp_mod_cells, which will exist if some
            modification has been made which also set temp_mod=True.
                c.write(temp_mod=True)'''

        if filename is None:
            filename = self.output_file
        
        cells = None
        # get original bytes with lsbs overwritten
        if self.keep_original:
            cells = self._overwrite_lsbs(temp_mod)
        else:
            if temp_mod:
                cells = self.temp_mod_cells.astype(self.cell_data_type)
            else:
                cells = self.cells.astype(self.cell_data_type)
            
        if compressed:
            # if no .gz at end of file, add it
            if not filename.endswith('.gz'):
                filename = filname + '.gz'
            with gzip.open(filename) as f:
               f.write(self.header_binary)
               cells.tofile(f)
        else:
            with open(filename, 'wb') as f:
               f.write(self.header_binary)
               cells.tofile(f)




    # checks if a point is in the cell space
    def _in_cell_space(self, point):
        if len(point) != self.cells.ndim:
            return False
        for i, dim in zip(point, self.cells.shape):
            if i < 0 or i > dim - 1:
                return False
        return True



    # point can be either length 3 or 1 or a scalar
    # get a slice through the etire cell space at the given point and orientation
    def _cut_1(self, point, orientation):
        # deal with possible input types
        D,H,L = 0,0,0
        if type(point) == int:
            D,H,L = point, point, point
        elif len(point) == 1:
            D,H,L = point[0], point[0], point[0]
        elif len(point) == 3:
            D,H,L = point

        # verify slice is in range
        bounds_check_point = (0,0,0)
        if orientation == 'd':
            bounds_check_point = (D,0,0)
        elif orientation == 'h':
            bounds_check_point = (0,H,0)
        elif orientation == 'l':
            bounds_check_point = (0,0,L)

        # TODO handle the error or point out of bounds
        if not self._in_cell_space(bounds_check_point):
            return None

        if orientation == 'd':
            return self.cells[D,:,:].copy()
        elif orientation == 'h':
            return self.cells[:,H,:].copy()
        elif orientation == 'l':
            return self.cells[:, :,L].copy()
        # TODO handle error or invalid orientation
        else:
            return None


    # show cut from all 3 directions for a given point
    def _cut_3(self, point, reorient=False):
        d_cut = self._cut_1(point, 'd')
        h_cut = self._cut_1(point, 'h')
        l_cut = self._cut_1(point, 'l')

        if reorient:
            l_cut = np.transpose(l_cut)

        return d_cut, h_cut, l_cut


    # takes in csp column and makes a height map
    # find the first location of an air cell
    def _find_air(self, column):
        air_indices = np.where(column == 3)[0]
        if air_indices.size > 0:
            return air_indices[-1]
        else:
            return np.array([0])

    # takes in csp column and makes a height map
    def _find_air_or_mobile(self, column):
        air_indices = np.nonzero(column == 3)[0]
        mobile_sand_indices = np.nonzero(column == 1)[0]

        if air_indices.size > 0 and mobile_sand_indices.size > 0:
            return max(air_indices[-1], mobile_sand_indices[-1])
        elif air_indices.size > 0:
            return air_indices[-1]
        elif mobile_sand_indices.size > 0:
            mobile_sand_indices[-1]
        else:
            return 0

    # get the position of the bottom of the ceiling
    def _find_ceiling_one_column(self, column):

        # get the indices of each type
        boundary = cell_types['boundary'] 
        neutral = cell_types['neutral']
        input_of_sand =  cell_types['input_of_sand']
        boundary_indices = np.nonzero(column == boundary)[0]
        neutral_indices = np.nonzero(column == neutral)[0]
        input_of_sand_indices = np.nonzero(column == input_of_sand)[0]

        # iterate through the types of indices
        indices_lists = [boundary_indices, neutral_indices, input_of_sand_indices]
        min_position = 0
        for indices in indices_lists:
            if indices.size == 0:
                continue
            elif indices.size == 1:
                min_position = max(min_position, indices[0])
            # find the longes contiguous streak of ceiling cell types starting from 0
            # the highest index value is the bottom of the ceiling for this column
            else:
                for i in range(1,len(indices)):
                    if indices[i] == 1 + indices[i-1]:
                        min_position = max(min_position, indices[i])
                    else:
                        break
        return min_position
                
                # look for break in contiguous values
                
    # create a heightmap of the celing and get the lowest point                
    def make_ceiling_map(self):
        self.ceiling_map = np.apply_along_axis(self._find_ceiling_one_column, 1, self.cells)  
        
        

        
    # creates a surface map
    # like an inverted and offset heightmap, gives the positions
    # of the air cells at the surface of the sand
    # really it is just the air cell in each column with the highest
    # index (and idices start at the physical top)
    def make_surface_map(self):
        self.surface_map = np.apply_along_axis(self._find_air_or_mobile, 1, self.cells)


    # make a HeightMap object
    def make_height_map(self):
        self.make_surface_map()
        height_map = self.cells.shape[1] - 1 - self.surface_map
        # this is how good I am at OOP
        self.height_map = heightmap.HeightMap(height_map)


    ########################################
    #########      modification      #######
    ########################################


    # add sand grain at (depth, surface_map[depth, length], length)]
    def add_sand(self, depth, length):
        # get surface map or immobile sand
        sm = self.surface_map
        height = sm[depth, length]
        self.cells[depth, height, length] = cell_types['grain']

    # adds a randomly placed grain of sand to the surface
    def add_sand_random(self):
        d, _, l = self.cells.shape
        depth = random.randint(0,d)
        length = random.randint(0,l)
        self.add_sand(depth, length)


    # prototype to add grains to surface
    def add_square(self, depth, length, k, temp_mod=False):
        sm = self.surface_map
        for d in range(depth-5, depth+6):
            for l in range(length-5, length+6):
                h = sm[d, l]
                # was range(h-k, h+1):
                for i in range(h-k, h+1):
                    if temp_mod:
                        self.temp_mod_cells[d,i,l] = cell_types['grain']
                    else:
                        self.cells[d,i,l] = cell_types['grain']


    # perform multiple random edits and write out the edited versions
    # return the paths to the new files
    def multiple_random_edits(self, num_edits, file_prefix):
        # pick a random spot that won't be out of bounds
        d, _, l = self.cells.shape
        paths = []
        for i in range(num_edits):
            depth = random.randint(5,d-5)
            length = random.randint(5,l-5)
            self.temp_mod_cells = self.cells.copy()
            self.add_square(depth, length,  8, temp_mod=True)
            filename = file_prefix + str(i) + '.csp'
            self.write(filename=filename, temp_mod=True)
            paths.append(filename)

            # make path names absolute
            for i in range(len(paths)):
                paths[i] = os.path.abspath(paths[i])
        return paths


    def multiple_random_pics(self, num_edits, base_image, file_prefix):
        d, _, l = self.cells.shape
        paths = []
        #scale_factor = random.randint(4,8)
        #scaled_image = np.kron(base_image, make_scaler(scale_factor, scale_factor)) * 6
        scaled_image = base_image
        s_d, s_l = scaled_image.shape

        for i in range(num_edits):
            depth = random.randint(1, d - s_d - 1)
            length = random.randint(1, l - s_l - 1)
            self.temp_mod_cells = self.cells.copy()
            self.add_height_map((depth, length), scaled_image, temp_mod=True)
            filename = file_prefix + str(i) + '.csp'
            self.write(filename=filename, temp_mod=True)
            paths.append(filename)

            # make path names absolute
            for i in range(len(paths)):
                paths[i] = os.path.abspath(paths[i])
        return paths


    # TODO vectorize this
    # adds a height_map to the sand with heigh_map positioned by top_left
    def add_height_map(self, top_left, height_map, temp_mod=False):

        x_min, y_min = top_left
        x_len, y_len, = height_map.shape
        #image_window = self.cells[x_min:x_min + x, y_min:y_min + y]

        sm = self.surface_map
        for x in range(x_len):
            for y in range(y_len):
                h = sm[x_min + x, y_min + y]
                for i in range(h-height_map[x,y]+1, h+1):
                    if temp_mod:
                        self.temp_mod_cells[x_min+x,i,y_min+y] = cell_types['grain']
                    else:
                        self.cells[x_min+x,i,y_min+y] = cell_types['grain']

    # add a sin wave to the cellspace
    def add_sine(self, amp=5, f=10, temp_mod=False):
        # mgrid of (l,d) values
        a = np.mgrid[f*-1.0:f*1.0:f*2.0/self.length, f*-1.0:f*1.0:f*2.0/self.depth]
        s = np.sin(a[0]) + 1.0
        s = (s.T * amp)
        s = (np.round(s)).astype(np.uint8)
        #hm(s)
        self.add_height_map((0,0),s, temp_mod=temp_mod)


    # types of mods available
    mods = ['sine', 'noise', 'guassian hill', 'space_invader', 'clip']


    # applies a random modification of the given type
    def random_mods(self, file_prefix, mod_type='gaussian', num_mods=1, temp_mod=False):

        # get heightest sand to avoid going out of the space
        
        sand_top = np.min(self.surface_map)
        ceiling_bottom = np.max(self.ceiling_map)
        
        # max alteration height
        h = sand_top - ceiling_bottom - 1
        if h < 2:
            # TODO deal with not enough space to modify
            return None

        paths = []
        # do the modification
        for i in range(num_mods):

            if temp_mod:
                self.temp_mod_cells = self.cells.copy()

            if mod_type == 'sine':
                # superimpose a sine wave in the wind direction with random
                # amplitude and frequency
                # calculate frequency and amplitude values
                amp = random.randint(2,h)
                f = random.randint(1, 20)
                self.add_sine(amp, f, temp_mod)
            elif mod_type == 'gaussian':
                # create a gaussian
                
                pass
            elif mod_type == 'space_invader':
                # superimpose a space invader onto the cell space surface
                # random size and location
                # create the image
                image_depth, image_length = heightmap.invader.shape
                depth, _, length = self.cells.shape
                max_length_scaling = length // image_length
                max_depth_scaling = depth // image_depth
                length_scaling = random.randint(1, max_length_scaling)
                depth_scaling = random.randint(1, max_depth_scaling)
                
                height_scaling = random.randint(2,h)
                scaling_matrix = heightmap.make_scaler(depth_scaling, length_scaling)
                scaled_image = np.kron(heightmap.invader, scaling_matrix) * height_scaling

                # now place the image
                scaled_image_depth, scaled_image_length = scaled_image.shape
                #breakpoint()
                top_left_depth = random.randint(0, depth - scaled_image_depth)
                top_left_length = random.randint(0, length - scaled_image_length)
                top_left = (top_left_depth, top_left_length)
                self.add_height_map(top_left, scaled_image, temp_mod)
            else:
                pass

            

            # now write out the modification fo file
            file_out = file_prefix + '--' + mod_type + '--' + str(i) + '.csp'
            self.write(filename=file_out, temp_mod=temp_mod)
            paths.append(file_out)


        # make path names absolute
        for i in range(len(paths)):
            paths[i] = os.path.abspath(paths[i])
            
        return paths
            

    ###################################
    #########      drawing      #######
    ###################################



    def draw_height_map(self):
        '''updates and draws height_map'''
        self.make_height_map()
        self.height_map.draw()


    def draw_surface_map(self):
        '''updates and draws surface_map'''
        self.make_surface_map()
        plt.imshow(self.surface_map.astype(np.float32))
        plt.colorbar()
        plt.show()

    def draw_fft_blur(self):
        '''updates and draws the fft of the height_map'''
        self.make_height_map()
        self.height_map.draw_fft_blur()


    # put image on axes
    @staticmethod
    def _draw_cross_section(cut, axes):
        rescal_color_map = _ints_to_colors(_cell_colors_list)
        n = colors.Normalize(vmin=0, vmax=11)
        axes.pcolormesh(cut, cma=rescal_color_map, norm=n)
        axes.invert_yaxis()
        return

    @staticmethod
    def _draw_window(cut, center, axes, orientation='d',
                     pad_horizontal=5, pad_vertical=5, reorient=False):
        rescal_color_map = _ints_to_colors(_cell_colors_list)
        n = colors.Normalize(vmin=0, vmax=11)
        axes.pcolormesh(cut, cmap=rescal_color_map, norm=n)

        D,H,L = center
        pv = pad_vertical
        ph = pad_horizontal
        right_max, top_max = cut.shape
        if orientation == 'd':
            bottom, top = max(H-ph,0), min(H+ph+1,right_max)
            left, right = max(L-pv,0), min(L+pv+1,top_max)
            axes.set_xlim(left=left, right=right)
            axes.set_ylim(bottom=bottom, top=top)
            axes.invert_yaxis()
        elif orientation == 'h':
            bottom, top = max(D-ph,0), min(D+ph+1,right_max)
            left, right = max(L-pv,0), min(L+pv+1,top_max)
            axes.set_xlim(left=left, right=right)
            axes.set_ylim(bottom=bottom, top=top)
            axes.invert_yaxis()
        elif orientation == 'l':
            # D and H switched due to transpose
            if reorient:
                D, H = H, D
            bottom, top = max(D-ph,0), min(D+ph+1,right_max)
            left, right = max(H-pv,0), min(H+pv+1,top_max)
            axes.set_xlim(left=left, right=right)
            axes.set_ylim(bottom=bottom, top=top)
            axes.invert_yaxis()


    def edit_and_view(self, center, new_cell_type, pad_horizontal=5, pad_vertical=5):
        '''show before and after editing of cell.
        show entire view from by slicing at given depth, height, and length befre the edit.
        also show a zoomed in view before and after the edit.
        
        arguments:
            center -- a list-like object of 3 integers.
                      the location at which the zoomed windows are centered.
                      the depth, length, and height of the cuts are determined by center.

            new_cell_type -- the cell type that is written to the cell space at the position given
                             by center

        keyword arguments:
            pad_horizontal -- the padding on each side of center for the zoomed in windows
                              the total window width is (2 * pad_horizontal + 1)

            pad_vertical   -- the padding above and below center for the zoomed in windows
                              the total window height is (2 * pad_vertical + 1)
        
        returns: None'''

        fig, axs = plt.subplots(nrows=3, ncols=3)

        # draw the whole space at the 3 perspectives
        d_slice, h_slice, l_slice = self._cut_3(center, reorient=True)
        self._draw_cross_section(d_slice, axs[0][0])
        self._draw_cross_section(h_slice, axs[1][0])
        self._draw_cross_section(l_slice, axs[2][0])


        self._draw_window(d_slice, center, axs[0][1], orientation='d',
                          pad_horizontal=pad_horizontal, pad_vertical=pad_vertical)
        self._draw_window(h_slice, center, axs[1][1], orientation='h',
                          pad_horizontal=pad_horizontal, pad_vertical=pad_vertical)
        self._draw_window(l_slice, center, axs[2][1], orientation='l',
                          pad_horizontal=pad_horizontal, pad_vertical=pad_vertical)

        # do edit
        if new_cell_type is not None:
            self.cells[center] = new_cell_type

        d_slice, h_slice, l_slice = self._cut_3(center, reorient=True)

        self._draw_window(d_slice, center, axs[0][2], orientation='d',
                    pad_horizontal=pad_horizontal, pad_vertical=pad_vertical)
        self._draw_window(h_slice, center, axs[1][2], orientation='h',
                    pad_horizontal=pad_horizontal, pad_vertical=pad_vertical)
        self._draw_window(l_slice, center, axs[2][2], orientation='l',
                    pad_horizontal=pad_horizontal, pad_vertical=pad_vertical)
        plt.show()
        return

# make a main() that takes command line args
def _process_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--height_map', action='store_true',
                        help='draw height map of input file')
    parser.add_argument('-f', '--input_file', required=True,
                        help='input file path')
    parser.add_argument('--fft_blur', action='store_true',
                        help='draw colormap of part of low frequency values of 2Dfft magnitudes')
    return parser.parse_args()


def _main():    
    args = _process_args()
    c = CellSpace(args.input_file)
    if args.height_map:
        c.draw_height_map()
    if args.fft_blur:
        c.draw_fft_blur()
