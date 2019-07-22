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


# deals with the heightmaps from ReSCAL
import numpy as np
import argparse
import scipy.ndimage

can_plot = True
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import colors
    from mpl_toolkits.mplot3d import Axes3D
except:
    can_plot = False
    pass

# little example of making a height map and scaling it
# space invader height map
invader = np.array([[0,0,1,0,0,0,0,0,1,0,0],
                    [0,0,0,1,0,0,0,1,0,0,0],
                    [0,0,1,1,1,1,1,1,1,0,0],
                    [0,1,1,0,1,1,1,0,1,1,0],
                    [1,1,1,1,1,1,1,1,1,1,1],
                    [1,0,1,1,1,1,1,1,1,0,1],
                    [1,0,1,0,0,0,0,0,1,0,1],
                    [0,0,0,1,1,0,1,1,0,0,0]], dtype=np.uint8)

# make a lists of lists that can be used in np.kron
# to scale up images by a factor of 'h' in the horizontal direction
# and a factor of 'v' in the vertical direction
def make_scaler(h,v):
    return [[1 for j in range(v)] for i in range(h)]
invaders = np.kron(invader, make_scaler(6,6)) * 6


# for finding height_map files
# used with param_space_exploration_utilities.get_files_to_process
path_glob = 'out/ALTI*'
exclude_globs = []



#### This function should be imported from analysis ####
# performs fft2d analysis on the data taken from input file
# and returns single quadrant fft result.
# data -> data to perform fft on
def fft2d_analyze(data):
    #Data points for x and y axis
    dpx, dpy = data.shape

    #Create x, y axis for graphing 2d
    x = np.arange(0,dpx)
    y = np.arange(0,dpy)

    #Set DC frequency bin to 0
    fft_data = data - np.mean(data)

    #Get fft2d and resize to single quadrant
    # PY3 Note: need indeces to be integers
    return np.fft.fft2(fft_data)[0:int(dpx/2),0:int(dpy/2)]*2/(dpx*dpy)

def fft2d_analyze_map_pic(data):
    #Data points for x and y axis
    dpx, dpy = data.shape

    #Create x, y axis for graphing 2d
    x = np.arange(0,dpx)
    y = np.arange(0,dpy)

    #Set DC frequency bin to 0
    fft_data = data - np.mean(data)

    #Get fft2d and resize to single quadrant
    # PY3 Note: need indeces to be integers
    return np.fft.fftshift(np.fft.fft2(fft_data)/(dpx*dpy))


# gets fft2d (which drops all but top left quadrant)
# grabs a small piece on the top left
# apply gaussian blur
def fft2d_crop_blur(image):
    # get just the magnitudes
    fft2d = np.absolute(fft2d_analyze(image))
    x,y = fft2d.shape
    fft2d = fft2d[:x//6, :y//6]
    fft2d = scipy.ndimage.gaussian_filter(fft2d, sigma=1)
    fft2d[0][0] = 0.0
    return fft2d


# make pics with fft stuff in the middle
def fft2d_center_blur(image):
    # get just the magnitudes
    fft2d = np.absolute(fft2d_analyze_map_pic(image))
    x,y = fft2d.shape
    fft2d = fft2d[x//3:2*(x//3), y//3:2*(y//3)]
    return scipy.ndimage.gaussian_filter(fft2d, sigma=1)

def make_surface(height_map):
    # get indices
    x,y = height_map.shape
    xs = list(range(x))
    ys = list(range(y))
    xss, yss = np.meshgrid(xs, ys)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(45, 20)
    surf = ax.plot_surface(xss, yss, np.transpose(height_map), cmap=cm.plasma,
                           linewidth=0, antialiased=False)
    return surf

# from the matplotlib tutorials
def plot_3d(height_map):
    surf = make_surface(height_map)
    plt.show()


class HeightMap:
    # expect the heightmap files from rescal
    # to be 2D arrays of ints in text form
    # TODO should maybe do a try catch
    def __init__(self, height_input):
        # case that a filename is given
        if type(height_input) == str:
            self.input_file = height_input
            self.output_file = height_input
            self.height_map = np.loadtxt(height_input).astype(np.uint8)

        # otehrwise is should be a numpy array
        else:
            self.input_file = None
            self.output_file = None
            self.height_map = height_input

        self.length, self.depth = self.height_map.shape
        # more params for analysis
        self.make_summary_data()


    # draw height map
    # draw fft
    # say the summary data
    def display_summary(self):
        fig, axs = plt.subplots(nrows=2, ncols=1)
        axs[0].pcolormesh(self.height_map)
        axs[1].pcolormesh(self.fft_blur)
        axs[0].invert_yaxis()
        axs[1].invert_yaxis()
        print('depth={d} length={l}'.format(d=self.depth,
                                            l=self.length))
        print('min_height={min_h} max_height={max_h}'.format(min_h=self.min_height,
                                                             max_h=self.max_height))
        print('average_height={mh} var_height={vh}'.format(mh=self.average_height,
                                                           vh=self.var_height))
        plt.show()


    # write out the height_map
    # default should be in the same form as ReSCAL
    # can also write to .npy
    def write(self, filename=None, npy=False):
        if filename is None:
            filename = self.output_file

        # meaning self.output_file is None
        # and this instance wasn't made constructed using a file
        if filename is None:
            # TODO maybe do some default filename
            # or an error
            return
        else:
            if not npy:
                np.savetxt(filename, self.height_map, fmt='%s')
            else:
                np.save(filename, self.height_map)


    def save_as_pdf(self, filename, in_3d=False):
        ''' '''
        # make a figure
        if not in_3d:
            plt.imshow(self.height_map)
        else:
            make_surface(self.height_map)
        # save with transparent background and a small bounding box
        plt.savefig(filename, transparent=True, bbox_inches='tight')



    # draw a simple color_map of height_map
    def draw(self):
        plt.imshow(self.height_map)
        plt.colorbar()
        plt.show()

    # makes a blurred 2D fft of the height_map
    def make_fft_blur(self):
        return fft2d_crop_blur(self.height_map.astype(np.float32))


    # draw simple color map of fft_blur
    def draw_fft_blur(self, in_3d=False):
        if not in_3d:
            plt.imshow(self.fft_blur)
            plt.colorbar()
            plt.show()
        else:
            plot_3d(self.fft_blur)

    def save_fft_blur(self, filename, in_3d=False):
        if not in_3d:
            plt.imshow(self.fft_blur)
        else:
            surf = make_surface(self.fft_blur)
        plt.savefig(filename, transparent=True, bbox_inches='tight')


    # draw simple color map of fft_blur
    def draw_fft_center(self, in_3d=False):
        if not in_3d:
            plt.imshow(self.fft_center)
            plt.colorbar()
            plt.show()
        else:
            plot_3d(self.fft_center)



    # get some basic metrics
    def make_summary_data(self):
        hm = self.height_map.astype(np.float32)
        self.average_height = hm.mean()
        self.min_height = hm.min()
        self.max_height = hm.max()
        self.var_height = hm.var()
        self.fft_blur = fft2d_crop_blur(self.height_map.astype(np.float32))
        self.fft_center = fft2d_center_blur(self.height_map.astype(np.float32))

    # given a 2D list of file paths
    # returns 2D list of HeightMap objects
    def load_from_paths(paths):
        return [[HeightMap(path) for path in local_paths] for local_paths in paths]


def _process_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--height_map', action='store_true',
                        help='draw height map of input file')
    parser.add_argument('-f', '--input_file', required=True,
                        help='input file path')
    parser.add_argument('--fft_blur', action='store_true',
                        help='draw colormap of part of low frequency values of 2Dfft magnitudes')
    parser.add_argument('-s', '--summary', action='store_true',
                        help='draw colormap and fft blur and summary data')
    return parser.parse_args()

# for running as shell command
def main():
    args = _process_args()
    h = HeightMap(args.input_file)
    if args.height_map:
        h.draw()
    if args.fft_blur:
        h.draw_fft_blur()
    if args.summary:
        h.display_summary()
