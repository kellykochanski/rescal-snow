__doc__ = '''Utilities to read, write, modify, and visualize snow-rescal cell spaces.'''
__author__ = 'Gian-Carlo DeFazio'
__date__ = 'August 16 2019'

import matplotlib
import os
# Matplotlib will fail if no display is available (e.g. many high-performance computing environments)
if bool(os.environ.get('DISPLAY', None)) == False:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as  colors
import struct
import numpy as np
import random
import gzip
import sys
import argparse
import scipy.ndimage
import random

import heightmap


_cell_colors_list = [
    'khaki',        # grain
    'turquoise',     # mobile_grain
    'red',           # vegetated_grain
    'white',         # air
    'forestgreen',   # vegetation
    'lightcyan',     # boundary
    'dimgray',       # neutral
    'darkorange',    # input_of_sand
    'orchid',        # output_of_sand
    'red',           # tunnel
    'lightcoral',    # EAUT
    'lemonchiffon',  # colored_grain
]


cell_types = {
    'grain' : 0,
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


def in_bounds(base, overlay, position):
    '''Evaluates if overlay can be placed onto base at 
    position. position is the coordinates of base on which the top left 
    corner of overlay is placed.'''
    # verify that the dimensions match
    if base.ndim != len(position):
        return False
    if overlay.ndim > base.ndim:
        return False

    # verify that position is in bounds for base
    if not all([-b <= p and p < b for b,p in zip(base.shape, position)]):
        return False
    
    # change any negative values in position to
    # positive base on base.shape
    position_positive = [p if (p>=0) else b+p for b,p in zip(base.shape, position)]

    # now get the portions of dimensions that are applicable based
    # on overlay.ndim, get only the overlay.ndim last dimensions
    base_dims = base.shape[-overlay.ndim:]
    position_positive = position_positive[-overlay.ndim:]

    # verify that position_positive + overlay_dims
    # will not exceed base dims
    all_dims = zip(base_dims, overlay.shape, position_positive)
    return all([o + p <= b for b,o,p in all_dims])


def non_border(a):
    '''Get a slice of a that excludes border cells.'''
    # the stuff we care about has states 0,1,2,3
    valid = np.nonzero(a < 4)[0]
    return a[valid[0]: valid[-1]+1]


# replaced with find_air_or_mobile
def surface_position(a):
    '''find index of first solid sand starting from the top.'''
    sand_indices = np.nonzero(a == 0)[0]
    if sand_indices.size > 0:
        return sand_indices[0]
    else:
        return 0


def find_air_or_mobile(column):
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


def shift_fill(a, shift, fill=0):
    '''Shift an array to the right by shift and fill in with fill value.
    If shift value is negative, the shift will be to the left.'''
    # nothing to do
    if shift == 0:
        pass
    # no point in shifting, entire array will
    # be overwritten by fill value
    elif abs(shift) >= len(a):
        a[:] = fill
    elif shift > 0:
        a[shift:] = a[:-shift]
        a[:shift] = fill
    else:
        a[:shift] = a[-shift:]
        a[shift:] = fill
    

def change_surface_level(column, delta):
    '''Given a column in the 3D cellspace, modify the sand height by delta.
    Clip at the boundaries of the space. Try to keep mobile sand by moving it with
    the surface.'''

    # if nothing to do, do nothing
    if delta == 0:
        return

    # negate delta because the height and indices are reversed
    delta = -delta

    # ignore the boundary cell at the top and bottom
    column_in_bounds = non_border(column)
    
    # get surface height
    surface = find_air_or_mobile(column_in_bounds) + 1
    
    # now calculate how far to actually move the sand
    # the actual delta may be clipped by top or bottom
    # also set the fill values
    # then shift the columns
    if delta < 0:
        clipped_delta = max(delta, -surface)
        fill_value = cell_types['grain']
        shift_fill(column_in_bounds, clipped_delta, fill_value)
    else:
        clipped_delta = min(delta, len(column_in_bounds) - surface)
        fill_value = cell_types['air']
        shift_fill(column_in_bounds, clipped_delta, fill_value)


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


class CellSpace:
    '''A snow-rescal cell space representation. This class deals with the .csp files
    made by rescal and genesis. The .csp files contain a header followed by the 3D space
    of cells as a 1D array.'''

    
    def __init__(self, filename):
        '''read in a .csp or .csp.gz file and create a CellSpace instance
        to represent the file contents. filename can also be a bytearray object
        that contains the contents of a .csp file.

        arguments:
            filename -- the path to the input file

        the following attributes will be created:

        header_size: the number of bytes in the header of filename
        model: the snow-rescal model
        height: the height of the cell space 
        length the length of the cell space
        depth: the depth of the cell space
        cell_size: the number of bytes per cell

        cell_data_type: the data type used for each cell self.cells
        cells_original: the orginal cell values
        cell: a copy of cell_original that can be modified
        header_binary: the verbatim binary header of filename

        output_file: the default output file name to write to
        height_map: a heightmap.HeightMap instance with the height data
        dims_2d: the dimensions of the corresponding height map'''
            
        # determine if data is in a file or is already bytes
        if isinstance(filename, bytes):
            self._read(bs=filename)
            self.output_file = None
        else:            
            self.input_file = filename
            # remove .gz by default
            if self.input_file.endswith('.gz'):
                self.output_file = filename
            else:
                self.output_file = filename[:-3]
            self._read()
            
        # make_height_map also makes surface_map
        self.make_height_map()
        self.dims_2d = self.dims_2d()


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


    def dims_2d(self):
        '''Get the dimensions of a height map for cell.'''
        d, _, l = self.cells.shape
        return d,l
    

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
        self.cells_original = cells
        self.cells = cells.astype(np.uint8)
        # make copy of the header for easy writing to file
        self.header_binary = bs[:self.header_size]



    # TODO deal with 'TIME'
    # TODO should all the chunk size checks be 12? I think not.
    # sets values of self.header_size
    # self.model, self.depth, self.length, self.height
    # self.cell_size
    def _read_header(self, bs):
        '''Read the .csp header and store the pertinent values

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

        
        # unpack first 20 bytes
        header_20 = struct.unpack(b0_19, bs[:20])
        # get the size of the header
        header_size = header_20[-1]
        self.header_size = header_size

        # grab the rest of the header
        header_20_to_end = bs[20:header_size]

        # find the other pertinent data
        # iterate through all remaining chunks and extract the important data
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
                    _, _, model, _ = struct.unpack(model_chunk,bs[bytes_slice])
                    self.model = model.decode('utf-8')

            elif chunk_type == b'SIZE':
                if chunk_size < 12:
                    print('bad chunk size for SIZE')
                else:
                    # extract and store height, length, and depth
                    _, _, h, l, d = struct.unpack(size_chunk, bs[bytes_slice])
                    self.height, self.length, self.depth = h, l, d

            elif chunk_type == b'CELL':
                if chunk_size < 12:
                    print('bad chunk size for CELL')
                else:
                    _, _, cell_size = struct.unpack(cell_size_chunk, bs[bytes_slice])
                    self.cell_size = cell_size

            elif chunk_type == b'TIME':
                pass

            elif chunk_type == b'BORD':
                if chunk_size < 12:
                    print('bad chunk size for BORD')
                else:
                    # extract and store height, length, and depth
                    _, _, h_b, l_b, d_b, u_b = struct.unpack(borders_chunk, bs[bytes_slice])
                    self.height_borders, self.length_borders, self.depth_borders, self.use_borders = h_b, l_b, d_b, u_b
            else:
                print('unrecognized chunk type {n}'.format(n=chunk_type.decode('utf-8')))

            chunk_start = chunk_start + chunk_size
        return


    def _overwrite_lsbs(self):
        '''creates cells that are identical to those read in (self.cells_original) except
        for the lsb which is rewritten with the values in self.cells.'''
        # expand bytes to type read in from the input file
        bytes_expanded = self.cells.astype(self.cells_original.dtype)
        mask = np.array([0xffffffffffffff00]).astype(self.cells_original.dtype)
        return (self.cells_original & mask) | bytes_expanded



    def restore_original_cells(self):
        '''restore self.cells to its orginal value.'''
        self.cells = self.cells_original.astype(np.uint8)
    

    # writes a .csp file
    # if no filename is given, the input file name is used
    def write(self, filename=None, compressed=False):
        ''' write the CellSpace to a file in the .csp format

        keyword arguments:
            filename --   the name of the file to write to.
                          by default, the name of the input file used to initialize
                          the CellSpace will be used. If a .csp.gz file was used, the .gz 
                          is removed. The default value is None, but this will cause
                          self.output_file to be used. Note that is is possible to
                          change self.output_file after initialization.
                          (default None)

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
                c.write(compressed=True)'''

        if filename is None:
            filename = self.output_file
        
        if compressed:
            # if no .gz at end of file, add it
            if not filename.endswith('.gz'):
                filename = filename + '.gz'
            with gzip.open(filename, 'wb+') as f:
                # get the bytes of the header and cells
                all_bytes = self.header_binary + self._overwrite_lsbs().tobytes()
                f.write(all_bytes)
        else:
            with open(filename, 'wb') as f:
               f.write(self.header_binary)
               self._overwrite_lsbs().tofile(f)


    # checks if a point is in the cell space
    def _in_cell_space(self, point):
        if len(point) != self.cells.ndim:
            return False
        for i, dim in zip(point, self.cells.shape):
            if i < 0 or i > dim - 1:
                return False
        return True


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

                        

    def add_height_map(self, height_map, top_left_corner=(0,0)):
        # verify that the height_map can be placed at top_left
        if not in_bounds(self.cells[:,0,:], height_map, top_left_corner):
            # TODO raise an exception, bad input
            return
        # now grab the applicable portion of self.cells
        d_top_left, l_top_left = top_left_corner
        d_height_map, l_height_map = height_map.shape
        cells_to_mod = self.cells[d_top_left : d_top_left + d_height_map,
                                  :,
                                  l_top_left : l_top_left + l_height_map]

        # now iterate over the applicable cells, and add the heightmap
        for i in range(d_height_map):
            for j in range(l_height_map):
                change_surface_level(cells_to_mod[i,:,j], height_map[i,j])
                

    def add_height(self, height):
        '''Raise the entire surface level by height.
        If height is negative surface level will lower.'''
        # create heightmap of height
        h = np.full(self.dims_2d, height)
        self.add_height_map(h)
                
                        
    def add_sinusoid(self, amplitude=5, frequency=10, phase=0, wind_direction=True):
        '''Superimposes a sinusoidal wave onto the surface of the cell space.'''
        # create a 2D heightmap that will be superimposed onto the surface
        sinusoid_height_map = heightmap.make_sinusoid(amplitude, frequency, self.dims_2d, phase=phase)
        self.add_height_map(sinusoid_height_map)


    # applies a random modification of the given type
    def random_mods(self, file_prefix, mod_type='gaussian', num_mods=1):

        # get heightest sand to avoid going out of the space
        
        sand_top = np.min(self.surface_map)
        ceiling_bottom = np.max(self.ceiling_map)
        
        # max alteration height
        h = sand_top - ceiling_bottom - 1
        if h < 2:
            # TODO deal with not enough space to modify
            return None

        d, _, l = self.shape
        paths = []
        # do the modification
        for i in range(num_mods):

            if mod_type == 'sine':
                # calculate parameters to add_sinusoid
                amp = random.randint(2,h)
                f = random.randint(1, min(d,l)//10)
                wind_direction = random.choice([True, False])
                phase = random.uniform(0, 2*math.pi)
                self.add_sinusoid(amp, f, phase=phase, wind_direction=wind_direction)
            elif mod_type == 'gaussian':
                # create a gaussian
                
                pass
            elif mod_type == 'space_invader':
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
                self.add_height_map(top_left, scaled_image)
            else:
                pass

            

            # now write out the modification fo file
            file_out = file_prefix + '--' + mod_type + '--' + str(i) + '.csp'
            self.write(filename=file_out)
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



    # point can be either length 3 or 1 or a scalar
    # get a slice through the entire cell space at the given point and orientation
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


if __name__ == '__main__':
    c = CellSpace('DUN.csp')
