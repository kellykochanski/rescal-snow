import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


# like decorate axes but takes in height_map ndarry
# was 'PiYG'
def decorate_axes(height_map, axes, fig, color_map='PiYG', levels=None, title=''):
    '''Creates an image of height_map and attaches it to axes.
       The default color sheme is pink, white, green.
       The color gradation levels are calculated from height_map unless
       other levels are passed in.
       No image title unless specified.
       A corresponding color bar is added to the figure.'''

    # get min and max value based on height_map if not passed in
    if levels is None:
        levels = MaxNLocator(nbins=15).tick_values(height_map.min(), height_map.max())

    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels
    # (this or pretty much from the matplotlib tutorial)
    cmap = plt.get_cmap(color_map)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # make the image and attach it to the axes, also add the colorbar
    im = axes.pcolormesh(height_map, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=axes)

    # put zero at top of y-axis
    axes.invert_yaxis()

    # title blank be default
    axes.set_title(title)


# for some ndarray a, get num_bins color binning values
def get_color_bins(a, num_bins=30):
    return MaxNLocator(nbins=num_bins).tick_values(a.min(),
                                                   a.max())






def show_interspersed(image_arrays):
    '''intersperse images from list of image arrays.
       x1, x2 etc are 2D arrays, diferent letters
       may have differnt sizes
       x = np.array([[x0,x1],
                     [x2,x3]])
       y = np.array([[y0,y1],
                     [y2,y3]])
       z = np.array([[z0,z1],
                     [z2,z3]])

      output:  x0 x1
               y0 y1
               z0 z1
               x2 x3
               y2 y3
               z2 z3
       '''


    color_bins = [get_color_bins(x) for x in image_arrays]

    fig, axes = plt.subplots(nrows=(image_arrays[0].shape[0] * len(image_arrays)),
                             ncols=(image_arrays[0].shape[1]))

    # get dims to iterate over
    rows, cols = image_arrays[0].shape[0], image_arrays[0].shape[1]
    for a in range(len(image_arrays)):
        for i in range(rows):
            for j in range(cols):
                decorate_axes(image_arrays[a][i,j], axes[i*len(image_arrays)+a, j],
                              fig, levels=color_bins[a])

    plt.show()
    return


def draw_comparison_chart(height_maps, all_comparisons, initial_comparisons):
    '''Creates a figure and associated axes's and processes the height_maps
       and comparison to create images and arranges those images on the figure.'''

    num_dirs = height_maps.shape[0]
    num_times = height_maps.shape[1]

    # make the figure and axes's
    fig, axes = plt.subplots(nrows=(num_dirs+len(all_comparisons)), ncols=(num_times+1))

    # get min and max values for height maps
    # and for the comparison maps to do consistent color grading
    height_maps_levels         = MaxNLocator(nbins=15).tick_values(height_maps.min(),
                                                                   height_maps.max())
    all_comparisons_levels     = MaxNLocator(nbins=15).tick_values(all_comparisons.min(),
                                                                   all_comparisons.max())
    initial_comparisons_levels = MaxNLocator(nbins=15).tick_values(initial_comparisons.min(),
                                                                   initial_comparisons.max())

    # makes the height_map images and place on figure
    for i in range(num_dirs):
        for j in range(num_times):
            decorate_axes(height_maps[i,j], axes[i][j],
                          fig, levels=height_maps_levels)

    # makes the all-way comparison images and place on figure
    for i in range(len(all_comparisons)):
        for j in range(num_times):
            decorate_axes(all_comparisons[i,j], axes[i+num_dirs][j],
                          fig, levels=all_comparisons_levels)

    # makes the initial comparison images and place on figure
    for i in range(num_dirs):
        decorate_axes(initial_comparisons[i], axes[i][-1],
                      fig, levels=initial_comparisons_levels)

    # make the unused axes's invisible
    for i in range(num_dirs, num_dirs + len(all_comparisons)):
        axes[i][-1].axis('off')

    fig.suptitle('all the datums', fontsize=14)

    plt.show()

    return fig
