'''
Rescal-snow: a cellular automaton model of self-organized snow
Copyright (C) 2019 Kelly Kochanski
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''



import struct
import numpy as np
import random
import gzip
import os
import sys
import argparse
import scipy.ndimage
import heightmap

# import drawing modules if they exist
ocan_plot = True
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator
    from matplotlib import colors
    from mpl_toolkits.mplot3d import Axes3D
except:
    can_plot = False
    pass


# creates a heightmap for a guassian
# scales the max height to h
# set all values to integers
def make_gaussian(h, d_padding, l_padding, sigma):

    # does it the cheap way, makes a dirac delta
    # meaning a 1 in the middle and 0 otherwise
    # then filters it
    dd = np.zeros([2 * l_padding + 1, 2 * d_padding + 1])
    l, d = dd.shape
    dd[l//2, d//2] = 1
    g = scipy.ndimage.gaussian_filter(dd, sigma)
    g = g * (h / g[l//2,d//2]) + 0.01
    return (np.round_(g)).astype(np.uint8)

# map of cell numbers to colors
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


# takes a list of colors and creats a colormap
# that maps the index of the color to the color when drawing
def _ints_to_colors(color_list):
    length = len(color_list)
    # get the rgba for each color
    c_rgba = np.array([colors.to_rgba_array(c) for c in color_list])
    c_rgba = np.reshape(c_rgba,
                         (c_rgba.shape[0],
                          c_rgba.shape[2]))
    return ListedColormap(c_rgba, N=12)

#
def _draw_cross_section(window, axes):
    n = matplotlib.colors.Normalize(vmin=0, vmax=11)
    axes.pcolormesh(window, cmap=rescal_color_map, norm=n)
    axes.invert_yaxis()
    return

# for finding cellspace files
path_glob = '*.csp*'
exclude_globs = ['DUN.csp']


# A ReSCAL cell space, now in python form
# deals with the .csp used by ReSCAL
#### TODO: I think my depth, height, and length are oriented differently
#### than ReSCAL :-(
class CellSpace:

    # cell typew from rescal defs.h
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

    ## description of the header using the struct library formatting strings

    # the first 20 bytes of a .csp header
    # [0..3]   B3s  the magic number which is \212 or 0o212
    #               in the source and 'CSP' <<138 'CSP'>>
    # [4..5]   2s   not sure what this is, but it's <<'@1'>>
    # [6]      c    can be <<'_'>> for little endian and <<'b'>> for big endian
    # [7]      c    a newline <<'\n'>>
    # [8..11]  i    a chunk size, this doesn't appear to be used,
    #               but it's 8 + sizeof(int) = 12  <<12>>
    # [12..15] 4s   s string that indicates the header size is next <<'HDSZ'>>
    # [16..19] i    the header size,
    #               which includes these first 20 bytes <<some_number>>
    b0_19 = 'B3s2scci4si'

    # after the first 20 bytes, the format is in chunks
    # each chunk starts with a number (int) that says its size in bytes
    # this size includes the 4 bytes for the chunk size

    # byte values are relative to the start of the chunk

    # Model chunk:  contains the model type, which I think is a 3 byte string
    # [0-3]   i  the size the chunk  <<some_number>> (probably 12)
    # [4-7]   4s the chunk type identifier <<'MODL'>>
    # [8-10]  3s the model name, for example <<'SNO'>>
    # [11]    B  I think this is pad <<'\0'>>
    model_chunk = 'i4s3sB'

    # Size chunk:  contains the dimensions of the simulations space
    # height (H), length (L), and depth (D) as 3 ints
    # [0-3]     i  the size the chunk  <<some_number>> (probably 20)
    # [4-7]     4s the chunk type identifier <<'SIZE'>>
    # [8-11]    i  the height H <<some_number>>
    # [12-15]   i  the length L <<some_number>>
    # [16-19]   i  the depth  D <<some_number>>
    size_chunk = 'i4siii'

    # Cell Size chunk: contains the number of bytes for each cell
    # which seems to default to 8
    # [0-3]     i  the size the chunk  <<some_number>> (probably 12)
    # [4-7]     4s the chunk type identifier <<'CELL'>>
    # [8-11]    i  the number of bytes per cell <<some_number>> (probably 8)
    cell_size_chunk = 'i4si'


    # read from a .csp file output from ReSCAL
    def __init__(self, filename, keep_original=True):
        self.input_file = filename
        self.keep_original = keep_original
        # remove .gz by default
        if self.input_file.endswith('.gz'):
            self.output_file = filename
        else:
            self.output_file = filename[:-3]

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

        self._read()

        # potentially used later
        self.temp_mod_cells = None

        # made during analysis
        self.height_map = None
        self.surface_map = None
        self.make_height_map()



    def _get_cell_data_type(self):
        '''get correct data type depending on cell size.
           in header.'''
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


    # read in CSP file and store the data in a 3D ndarray
    # store the meta data in a dictionary
    def _read(self):

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
        cells = np.reshape(cells, (self.depth,self.height,self.length))
        # may choose to discard original cells to save space
        if self.keep_original:
            self.cells_original = cells
        self.cells = cells.astype(np.uint8)
        # make copy of the header for easy writing to file
        self.header_binary = bs[:self.header_size]
        return


    # TODO deal with 'TIME'
    # sets values of self.header_size
    # self.model, self.depth, self.length, self.height
    # self.cell_size
    def _read_header(self, bs):
        '''read the csp header.'''

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
            else:
                print('unrecognized chunk type {n}'.format(n=chunk_type.decode('utf-8')))

            chunk_start = chunk_start + chunk_size

        return




    # creates cells that are identical to those read in (self.cells_original) except
    # for the lsb which is rewritten with the values in self.cells
    # can use a temporary modification to write out
    def _overwrite_lsbs(self, temp_mod=False):

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
    def write(self, filename=None, temp_mod=False):

        if filename is None:
            filename = self.output_file

        f = open(filename, 'wb')

        # write the header
        f.write(self.header_binary)

        cells = None
        # get original bytes with lsbs overwritten
        if self.keep_original:
            cells = self._overwrite_lsbs(temp_mod)
        else:
            if temp_mod:
                cells = self.temp_mod_cells.astype(self.cell_data_type)
            else:
                cells = self.cells.astype(self.cell_data_type)
        # write the cells
        cells.tofile(f)
        f.close()



    # checks if a point is in the cell space
    def _in_cell_space(self, point):
        if len(point) != self.cells.ndim:
            return False
        for i, dim in zip(point, self.cells.shape):
            if i < 0 or i > dim - 1:
                return False
        return True



    # point can be either length 3 or 1 or a scalar
    # get a slice through the etire cell spcae at the given point and orientation
    def cut(self, point, orientation):
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
    def cut_3(self, point, reorient=False):
        d_cut = self.cut(point, 'd')
        h_cut = self.cut(point, 'h')
        l_cut = self.cut(point, 'l')

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
        sm = self.surface_map()
        height = sm[depth, length]
        self.cells[depth, height, length] = CellSpace.cell_types['grain']

    # adds a randomly placed grain of sand to the surface
    def add_sand_random(self):
        d, _, l = self.cells.shape
        depth = random.randint(0,d)
        length = random.randint(0,l)
        self.add_sand(depth, length)


    # prototype to add grains to surface
    def add_square(self, depth, length, k, temp_mod=False):
        sm = self.surface_map()
        for d in range(depth-5, depth+6):
            for l in range(length-5, length+6):
                h = sm[d, l]
                # was range(h-k, h+1):
                for i in range(h-k, h+1):
                    if temp_mod:
                        self.temp_mod_cells[d,i,l] = CellSpace.cell_types['grain']
                    else:
                        self.cells[d,i,l] = CellSpace.cell_types['grain']


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

        sm = self.surface_map()
        for x in range(x_len):
            for y in range(y_len):
                h = sm[x_min + x, y_min + y]
                for i in range(h-height_map[x,y]+1, h+1):
                    if temp_mod:
                        self.temp_mod_cells[x_min+x,i,y_min+y] = CellSpace.cell_types['grain']
                    else:
                        self.cells[x_min+x,i,y_min+y] = CellSpace.cell_types['grain']

    # add a sin wave to the cellspace
    def add_sine(self, amp=5, f=10):
        # mgrid of (l,d) values
        a = np.mgrid[f*-1.0:f*1.0:f*2.0/self.length, f*-1.0:f*1.0:f*2.0/self.depth]
        s = np.sin(a[0]) + 1.0
        s = (s.T * amp)
        s = (np.round(s)).astype(np.uint8)
        hm(s)
        self.add_height_map((0,0),s)


    ###################################
    #########      drawing      #######
    ###################################


    # sets standard colors for ReSCAL cell types
    _rescal_color_map = _ints_to_colors(_cell_colors_list)


    def draw_height_map(self):
        self.make_height_map()
        self.height_map.draw()


    def draw_surface_map(self):
        self.make_surface_map()
        plt.imshow(self.surface_map.astype(np.float32))
        plt.colorbar()
        plt.show()

    def draw_fft_blur(self):
        self.make_height_map()
        self.height_map.draw_fft_blur()


    # put image on axes
    @staticmethod
    def _draw_cross_section(cut, axes):
        n = matplotlib.colors.Normalize(vmin=0, vmax=11)
        axes.pcolormesh(cut, cmap=CellSpace._rescal_color_map, norm=n)
        axes.invert_yaxis()
        return

    @staticmethod
    def _draw_window(cut, center, axes, orientation='d',
                     pad_horizontal=5, pad_vertical=5, reorient=False):
        n = matplotlib.colors.Normalize(vmin=0, vmax=11)
        axes.pcolormesh(cut, cmap=CellSpace._rescal_color_map, norm=n)

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

    # show before and after editing of cell
    # show entire view
    # show zoomed in windows before and after edit
    def edit_and_view(self, center, new_cell_type, pad_horizontal=5, pad_vertical=5):
        fig, axs = plt.subplots(nrows=3, ncols=3)

        # draw the whole space at the 3 perspectives
        d_slice, h_slice, l_slice = self.cut_3(center, reorient=True)
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

        d_slice, h_slice, l_slice = self.cut_3(center, reorient=True)

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

# for running as shell command
def main():
    args = _process_args()
    c = CellSpace(args.input_file)
    if args.height_map:
        c.draw_height_map()
    if args.fft_blur:
        c.draw_fft_blur()
