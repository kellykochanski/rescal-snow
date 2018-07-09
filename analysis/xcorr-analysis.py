import math
import time
import glob
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit

if len(sys.argv) < 3:
    print('Please provide a path to altitude files, and a output directory name')
    exit()
else:
    alti_path = str(sys.argv[1])
    output_dir = str(sys.argv[2])

# Run options:
show_graph = False
alti_file_offset = 100
time_correlation = True
output_directory = 'output_files/' + output_dir + '/xcor/'

# Output files
cross_correlations_file = open( output_directory + '/cross-correlations.txt','w')
#log_file = open(output_directory + '/xcorr.log', 'w')
meta_analysis_file = open(output_directory + '/meta_analysis.txt', 'w')

current_file_id = 00000
next_file_id = current_file_id+alti_file_offset
alti_file_1 = alti_path + '/ALTI%05.0f_t0.log' % current_file_id
alti_file_2 = alti_path + '/ALTI%05.0f_t0.log' % next_file_id

# Vector to store the calculated velocities
y_values = []
x_values = []

correlation_times = []

print(alti_file_2)

# Compute the velocity of dunes between altitude files while they exist
while (os.path.exists(alti_file_1) and os.path.exists(alti_file_2)):
    
    # Load the alittudes
    altitudes_1_matrix = np.loadtxt(alti_file_1).T
    altitudes_2_matrix = np.loadtxt(alti_file_2).T
    
    print('Analyising t0: ' + str(current_file_id))
    
    # Get the width and length of cell space
    signal_length = altitudes_1_matrix.shape[0]
    signal_width = altitudes_1_matrix.shape[1]
    
    # Repeat the second altitudes in a tile pattern for accurate cross correlation
    altitudes_2_matrix_tiled = np.tile(altitudes_2_matrix, [3,3])
    
    # Compute Cross Correlation
    start = time.time()
    xcorr = signal.correlate2d(altitudes_1_matrix, altitudes_2_matrix_tiled, mode='same')
    correlation_times.append(time.time() - start)
    
    # Save the cross correlation for future use
    cross_correlations_file.write(str(current_file_id) + ':\n')
    np.savetxt(cross_correlations_file, xcorr, delimiter=',')
    cross_correlations_file.write('\n\n')

    # Get the indices maximum values along width and length
    max_width = np.argmax(np.max(xcorr, axis=0))+1
    max_length = np.argmax(np.max(xcorr, axis=1))+1
    
    # Adjust the indices so that they are zero-centered from the middle of the cross correlation matrix
    offset_width = max_width - signal_width/2
    offset_length = max_length - signal_length/2
    
    # Save the velocity and invert its sign
    y_values.append(np.float64(-offset_length))
    x_values.append(np.float64(current_file_id))
    
    # Increment file ids
    current_file_id += alti_file_offset
    next_file_id += alti_file_offset
    alti_file_1 = alti_path + '/ALTI%05.0f_t0.log' % current_file_id
    alti_file_2 = alti_path + '/ALTI%05.0f_t0.log' % next_file_id

meta_analysis_file.write('Average correlation time: ' + str(np.mean(correlation_times)))

# Remove the first velocity as it is garbage data
y_values = y_values[1:]
x_values = x_values[1:]

# Save the velocities for future use
meta_analysis_file.write('\n\nVelocities omiting first, step size of ' + str(alti_file_offset) + ': ' + str(y_values))

############# Try to fit the velocity points #############
def function_to_fit(t, A, b, m, c):
    # Cast all elements of t to integers
    return ((A) * math.e**((t)/(-b))) + ((m) * (t) + (c))
try:
    popt, pcov = curve_fit(function_to_fit, x_values, y_values, p0 = [1,50,1,1])
    A, b, m, c = popt
    calc_y_values = function_to_fit(np.asarray(x_values), float(A), float(b), float(m), float(c))

    # Calculate R^2
    y_bar = np.mean(y_values)
    SS_tot = np.sum((y_values - y_bar)**2)
    SS_res = np.sum((y_values - calc_y_values)**2)
    r_squared = 1 - (SS_res/SS_tot)

    # Save and print out calculated values
    meta_analysis_file.write('\n')
    meta_analysis_file.write('Initial state coefficient:    ' + str(A))
    meta_analysis_file.write('Rate of approach:        ' + str(b))
    meta_analysis_file.write('Slope of quasi-stable state:    ' + str(m))
    meta_analysis_file.write('State intercept:        ' + str(c))
    meta_analysis_file.write('Initial velocity:        ' + str(y_values[0]))
    meta_analysis_file.write('Final velocity:            ' + str(y_values[-1]))
    meta_analysis_file.write('R Squared value:        ' + str(r_squared))
    meta_analysis_file.write('\n')

#    print('\n')
#    print('Initial state coefficient:    ' + str(A))
#    print('Rate of approach:        ' + str(b))
#    print('Slope of quasi-stable state:    ' + str(m))
#    print('State intercept:        ' + str(c))
#    print('Initial velocity:        ' + str(y_values[0]))
#    print('Final velocity:            ' + str(y_values[-1]))
#    print('R Squared value:        ' + str(r_squared))
#    print('\n')

    # Graphing
    plt.figure('plot')
    plt.plot(x_values,y_values)
    plt.plot(x_values, calc_y_values)
    plt.savefig(output_directory + '/velocity_with_fit.png')
    if show_graph:
        plt.show()

except:
    print('Error fitting function to velocity data, probably too few examples.')
    meta_analysis_file.write('\n')
    meta_analysis_file.write('Initial velocity:        ' + str(y_values[0]))
    meta_analysis_file.write('Final velocity:            ' + str(y_values[-1]))
    meta_analysis_file.write('\n')
