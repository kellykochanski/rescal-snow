__author__ = 'Gian-Carlo DeFazio'
__date__ = '16 August 2019'
__doc__ = r"""
Utilities to create, read and visualize Rescal-snow height maps. These are managed through a HeightMap class.

Utilities for creating Rescal-snow heightmaps (these may be used as initial conditions for the Rescal-snow simulation; see the INPUT_ELEVATION CSP_TEMPLATE in `src/genesis.c`).
- invader_template
 - gaussian_hill
 - make_sinusoid
 - scale

Utilities for analyzing Rescal-snow heightmaps (most of these use fourier transforms, as snow/sand self-organization has strong emergent wavelengths):
 - fft2d_analyze
 - fft2d_analyze_map_pic
 - fft2d_crop_blur
 - fft2d_center_blur

 Utilities for visualizing height maps:
  - make_surface
"""

import sys
import os
import math
import numpy as np
import scipy.ndimage
import matplotlib
# Matplotlib will fail if no display is available (e.g. many high-performance computing environments)
if bool(os.environ.get('DISPLAY', None)) == False:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# space invader height map
invader_template = np.array([[0,0,1,0,0,0,0,0,1,0,0],
                             [0,0,0,1,0,0,0,1,0,0,0],
                             [0,0,1,1,1,1,1,1,1,0,0],
                             [0,1,1,0,1,1,1,0,1,1,0],
                             [1,1,1,1,1,1,1,1,1,1,1],
                             [1,0,1,1,1,1,1,1,1,0,1],
                             [1,0,1,0,0,0,0,0,1,0,1],
                             [0,0,0,1,1,0,1,1,0,0,0]], dtype=np.int32)


def gaussian_hill(amplitude, sigma, height_padding, width_padding):
    """
    Creates a guassian heightmap. The center of the guassian
    will have a value of amplitude. sigma can be either a scaler 
    or a 2-tuple. The padding determines the size of the array returned.
    The size is (2*height_padding+1, 2*width_padding+1)

    Inputs:
     - amplitude: height in cells of gaussian hills at center (highest point)
     - sigma: width (at standard deviation) of hills. 
          may be scalar (symmetric hills) or 2-tuple (asymmetric)
     - height_padding, width_padding: determine size of array returned
     """
    
    # does it the cheap way:
    # makes an impulse, meaning a 1 in the middle and 0 otherwise,
    # then filters it
    impulse = np.zeros([2 * height_padding + 1, 2 * width_padding + 1])
    height, width = impulse.shape
    impulse[height//2, width//2] = 1
    gaussian = scipy.ndimage.gaussian_filter(impulse, sigma)
    # scale so center has a height of amplitude
    gaussian = gaussian * (amplitude / gaussian[height//2,width//2])
    return np.round_(gaussian).astype(np.int32)


def make_sinusoid(height, frequency, dims, phase=0, wind_direction=True, no_negative=False):
    """
    Create a sinusoid on a 2D plane that has dimensions dims. 
    Does an inverse Fourier transform and discretizes. 
    The wave is resized so that its height ranges from 0 to height. 
    The frequency is relative to the space. If the frequency is x, wave will 
    complete x cycles across the space. If wind_direction is True, then the wave
    heights will vary in the horizontal direction, otherwise they will vary in the
    vertical direction.

    Inputs:
     - height: wave ranges from 0 to height
     - dims: dimensions of array
     - wind_direction: boolean, if True wave varies in horizontal, if False in vertical
     - no_negative: disallow negative heights, returning zeros as minimum

    """
    # create a 2D array of the same size as the dims
    grid = np.zeros(dims, dtype=np.complex)

    # convert phase to a complex number
    complex_val = math.cos(phase) + (math.sin(phase) * 1.0j)
    
    # add the frequency to the grid
    if wind_direction:
        grid[0,frequency] = complex_val
    else:
        grid[frequency,0] = complex_val
        
    # make the sinusoid
    height_map = np.fft.ifft2(grid).real

    # scale to amplitude
    if height_map.max() != 0.0:
        height_map = height_map *  height / (2 * height_map.max())
        
    # now get rid of the negatives and round
    if no_negative:
        height_map = height_map - np.min(height_map)
    return np.round_(height_map).astype(np.int32)


def scale(height_map, amplitude=1, vertical=1, horizontal=1):
    """
    Scales a height map. The heigh map dimensions are scaled in the 
    horizontal and vertical directions. The values are then scaled by a factor
    of amplitude.

    Inputs:
     - height_map
     - amplitude: scale factor for stretching height map
     - vertical, horizontal: relative scale factors in horizontal and vertical dimensions
    """
    return np.kron(height_map, np.full((vertical, horizontal), amplitude))


def fft2d_analyze(data):
    """
    This function should be imported from analysis
    performs fft2d analysis on the data taken from input file
    and returns single quadrant fft result.

    data -> data to perform fft on
    """

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
    """
    May be redundant with fft2d_analyze?
    """
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


def fft2d_crop_blur(image):
    """
    gets fft2d (which drops all but top left quadrant)
    grabs a small piece on the top left
    apply gaussian blur
    """

    # get just the magnitudes
    f_image = image.astype(np.float32)
    fft2d = np.absolute(fft2d_analyze(f_image))
    x,y = fft2d.shape
    # need to avoid setting array dimensions to zero here
    fft2d = fft2d[:max(x//6,1), :max(y//6,1)]
    fft2d = scipy.ndimage.gaussian_filter(fft2d, sigma=1)
    fft2d[0][0] = 0.0
    return fft2d


def fft2d_center_blur(image):
    """make pics with fft stuff in the middle"""
    
    # get just the magnitudes
    fft2d = np.absolute(fft2d_analyze_map_pic(image))
    x,y = fft2d.shape
    fft2d = fft2d[x//3:2*(x//3), y//3:2*(y//3)]
    return scipy.ndimage.gaussian_filter(fft2d, sigma=1)


def make_surface(height_map):
    """Plot and return a height map as a 3D surface using matplotlib.pyplot.plot_surface, plasma colors"""
    # get indices
    x,y = height_map.shape
    xs = list(range(x))
    ys = list(range(y))
    xss, yss = np.meshgrid(xs, ys)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.view_init(45, 20)
    return ax.plot_surface(xss, yss, np.transpose(height_map), cmap=cm.plasma,
                           linewidth=0, antialiased=False)


def draw(a):
    """draw a simple color_map of an array"""
    plt.imshow(a)
    plt.colorbar()
    plt.show()


def plot_3d(a):
    """draw and show a 3D surface of from a 2D array of heights"""
    surf = make_surface(a)
    plt.show()


class HeightMap:

    # TODO should maybe do a try catch
    def __init__(self, height_input):
        """
        Input
         - height_input, a rescal elevation array (a 2D array of integers)
             or a filename pointing to such an array in a text file
        """ 
        # case that a filename is given
        if type(height_input) == str:
            self.input_file = height_input
            self.output_file = height_input
            try:
                self.height_map = np.loadtxt(height_input).astype(np.uint8)
            except:
                print("HeightMap couldn't read file " + height_input + ".")
                print("This probably happened because the input file is an irregular shape.")

        # otehrwise is should be a numpy array
        else:
            self.input_file = None
            self.output_file = None
            self.height_map = height_input

        self.length, self.depth = self.height_map.shape
        # more params for analysis
        self.make_summary_data()


    def display_summary(self):
        """
        Create plot of height map,
        its fourier transform,
        and summary data
        """
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


    def write(self, filename=None, npy=False):
        """
        write out the height_map
        default should be in the same form as ReSCAL
        can also write to .npy
        """
        
        if filename is None:
            if self.output_file is None:
                sys.stderr.write('ERROR: No filename given and no default file name. Cannot write HeightMap to file.\n')
                return
            else:
                filename = self.output_file
        else:
            if not npy:
                np.savetxt(filename, self.height_map, fmt='%s')
            else:
                np.save(filename, self.height_map)


    def save_as_pdf(self, filename, in_3d=False):
        """Save a HeightMap figure"""
        plt.xticks([])
        plt.yticks([])
        
        if not in_3d:
            plt.imshow(self.height_map)
        else:
            make_surface(self.height_map)
        # save with transparent background and a small bounding box
        plt.savefig(filename, transparent=True, bbox_inches='tight')

    def save_color_map(self, filename):
        """Draw a color map and save as png"""
        plt.imshow(self.height_map)
        plt.colorbar()
        plt.savefig(filename, bbox_inches='tight', transparent=True)
        plt.close()

    def draw(self):
        """draw a simple color_map of height_map""" 
        plt.imshow(self.height_map)
        plt.colorbar()
        plt.show()

    
    def make_fft_blur(self):
        """makes a blurred 2D fft of the height_map"""        
        return fft2d_crop_blur(self.height_map.astype(np.float32))

    
    def draw_fft_blur(self, in_3d=False):
        """draw simple color map of fft_blur"""
        if not in_3d:
            plt.imshow(self.fft_blur)
            plt.colorbar()
            plt.show()
        else:
            plot_3d(self.fft_blur)

            
    def save_fft_blur(self, filename, in_3d=False):
        """save fftblur to file"""
        if not in_3d:
            plt.imshow(self.fft_blur)
        else:
            surf = make_surface(self.fft_blur)
        plt.savefig(filename, transparent=True, bbox_inches='tight')


    def draw_fft_center(self, in_3d=False):
        """draw simple color map of fft_blur"""
        if not in_3d:
            plt.imshow(self.fft_center)
            plt.colorbar()
            plt.show()
        else:
            plot_3d(self.fft_center)

    
    def make_summary_data(self):
        """get some basic analysis metrics"""
        hm = self.height_map.astype(np.float32)
        self.average_height = hm.mean()
        self.min_height = hm.min()
        self.max_height = hm.max()
        self.var_height = hm.var()
        self.fft_blur = fft2d_crop_blur(hm)
        self.fft_center = fft2d_center_blur(hm)


