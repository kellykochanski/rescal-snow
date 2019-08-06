import numpy as np
import os
import sys
import glob
from scipy.ndimage import gaussian_filter
#import fft2d_analysis



#### Should be gotten by importing  fft2d_analysis
#### but no pandas on shattuck, so just getting these functions by pasting

#### Start fft2d_analysis ####

#Performs fft2d analysis on the data taken from input file and returns single quadrant fft result.
#data -> data to perform fft on
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

# vectorized version of fft2d_analyze
# signature means 2D to 2D
fft2d_analyze_vec = np.vectorize(fft2d_analyze, signature='(a,b)->(c,d)')


#### End   fft2d_analysis ####




# crop a 2D ndarray, cuts 1/4 off all edges
def crop(image):
    x,y = image.shape
    return image[x//4:(3*(x//4)), y//4:(3*(y//4))]
crop_vec = np.vectorize(crop, signature='(a,b)->(c,d)')

# set mean to zero for ndarry x
def zero_mean(x):
    return x - x.mean()

zero_mean_vec = np.vectorize(zero_mean, signature='(a,b)->(a,b)')



# easy metrics
def L2(x, y):
    return np.mean(np.square(y - x), axis=0)

def L1(x, y):
    return np.mean(np.absolute(y - x))

#def LN


def fft_comp(x, y, n=10):
    ''' finds largest n values in each, then compares them
    and averages over the numeber used.
    x and y are 2D ndarrays assumed to be the same size and much
    larger than n'''

    # get amplitudes of x and y, put into linear array
    fx = np.absolute(fft2d_analyze(x)).ravel()
    fy = np.absolute(fft2d_analyze(y)).ravel()

    # get locations of n largest values in x and y
    lx = fx.argsort()[:n]
    ly = fy.argsort()[:n]
    lxy = np.union1d(lx,ly)

    # do a L2 means suared differnce on the largest values
    return L2(fx[lxy], fy[lxy])


def fft_comp_all(x,y):

    # get amplitudes of x and y, put into linear array
    fx = np.absolute(fft2d_analyze(x))
    fy = np.absolute(fft2d_analyze(y))

    return L2(fx,fy)

standard_metrics = [L1, L2]

# gets fft2d (which drops all but top left quadrant)
# grabs top left quadrant (so now top left 1/16 of original)
# applie gaussian blur
def fft2d_crop_blur(image):
    # get just the magnitudes
    fft2d = np.absolute(fft2d_analyze(image))
    x,y = fft2d.shape
    fft2d = fft2d[:x//2, :y//2]
    return gaussian_filter(fft2d, sigma=1)

fft2d_crop_blur_vec = np.vectorize(fft2d_crop_blur, signature='(a,b)->(c,d)')


# generates a num_dirs x num_dirs ndarray of comparisons
def all_ways_compare(metric, height_maps):
    comparisons = np.empty([height_maps.shape[0], height_maps.shape[0]])
    for i in range(height_maps.shape[0]):
        for j in range(height_maps.shape[0]):
            comparisons[i,j] = metric(height_maps[i], height_maps[j])
    return comparisons

# generate a 1 x num_times array of comparisons
def initial_compare(metric, height_maps):
    comparisons = np.empty([height_maps.shape[0]])
    for j in range(height_maps.shape[0]):
        comparisons[j] = metric(height_maps[0], height_maps[j])
    return comparisons

# does an all_way_compare for all metrics and all time_steps
def do_all_comparisons(height_maps, metrics):
    num_dirs = height_maps.shape[0]
    num_times = height_maps.shape[1]
    all_comparisons = np.empty([len(metrics), num_times, num_dirs, num_dirs])
    for i in range(len(metrics)):
        for j in range(num_times):
            all_comparisons[i,j] = all_ways_compare(metrics[i], height_maps[:, j])
    return all_comparisons


def do_initial_comparisons(height_maps, metrics):
    num_dirs = height_maps.shape[0]
    num_times = height_maps.shape[1]
    initial_comparisons = np.empty([num_dirs, len(metrics), num_times])
    for i in range(num_dirs):
        for m in range(len(metrics)):
            initial_comparisons[i,m] = initial_compare(metrics[m], height_maps[i,:])
    return initial_comparisons

# given a top directory, which contains sub directories of height maps
# puts all height maps and comparisons into a big picture
# TODO deal with no files found
def make_heights_and_comparisons(top_directory, metrics, zero_mean=False, crop=False):

    #metrics = standard_metrics + metrics

    # 2D array of absolute paths to all the height maps
    # gives dimensions of height maps
    directories, x, y = get_alti_files(top_directory)
    num_dirs = len(directories)
    num_times = len(directories[0])

    # read in all files to giant ndarray
    # height_maps = np.empty([num_dirs, num_times, x, y], dtype=np.uint8)
    height_maps = np.empty([num_dirs, num_times, x, y], dtype=np.float32)
    for i in range(num_dirs):
        for j in range(num_times):
            height_maps[i,j] = np.loadtxt(directories[i][j])

    # set average height to zero for each height map
    if zero_mean:
        height_maps = zero_mean_vec(height_maps)

    if crop:
        height_maps = crop_vec(height_maps)

    initial_comparisons = do_initial_comparisons(height_maps, metrics)
    all_comparisons = do_all_comparisons(height_maps, metrics)

    return height_maps, all_comparisons, initial_comparisons

def make_compressed_heights_and_comparisons(filename, top_directory, metrics):
    '''Returns a .npz archive of the heights and comparison data.'''

    h, a, i = make_heights_and_comparisons(top_directory, metrics)
    compress_heights_and_comparisons(filename, h, a, i)

def compress_heights_and_comparisons(filename, height_maps, all_comparisons, initial_comparisons):
    '''Creates a file <filename>.npz which is a compressed archive of
   height_maps, all_comparisons, initial_comparisons.'''

    np.savez_compressed(filename,
                        height_maps=height_maps,
                        all_comparisons=all_comparisons,
                        initial_comparions=initial_comparisons)










def extract_heights_and_comparisons(filename):
    '''Extracts height_maps, all_comparisons, initial_comparisons from an .npz archive.'''
    data = np.load(filename)
    return data['height_maps'], data['all_comparisons'], data['initial_comparisons']

# TODO ensure alti_paths are absolute paths
def get_alti_files(top_directory):
    '''Finds all the ALTI* files which are the height
       map outputs of rescal.
       Assumes that the top directory is full of directories,
       each of which contains a directory called 'out',
       and each 'out' contains the same set of ALTI* file names
       and each ALTI file has the same dimensions.'''

    # get all the directories in top_directory (from stack overflow 973473)
    output_dirs = [f.path for f in os.scandir(top_directory) if f.is_dir()]
    for i in range(len(output_dirs)):
        output_dirs[i] = os.path.abspath(output_dirs[i])

    # make 2D array or absolute paths to all the ALTI* files
    alti_paths = []
    #print(output_dirs)
    for output_dir in sorted(output_dirs):
        local_alti_files = glob.glob(output_dir + '/out/ALTI*')
        # if directory has no alti files, don't append empty list
        if local_alti_files:
            alti_paths.append(sorted(local_alti_files))

    # find the dimensions of one of the files
        x = 0
        y = 0
        if alti_paths and alti_paths[0]:
            # alti_paths[0][0] must exist because alti_paths should not contain empty lists
            test_load = np.loadtxt(alti_paths[0][0])
            x, y = test_load.shape
    return alti_paths, x, y
