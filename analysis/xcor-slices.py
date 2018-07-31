import math
import time
import glob
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit

if len(sys.argv) < 4:
    print('Please provide a path to altitude files, an output directory name, and a name to prepend output files with.')
    exit()
else:
    alti_path = str(sys.argv[1])
    output_dir = str(sys.argv[2])
    output_filename = str(sys.argv[3])

# Run options:
show_graph = False
alti_file_offset = 100
time_correlation = True

# Output files
cross_correlations_file = open( output_dir + '/' + output_filename + '_cross_correlations.txt','w')
meta_analysis_file = open(output_dir + '/' + output_filename + '_meta_analysis.txt', 'w')

# Add file offset to skip the first file as cross correlation with flat surface is ill defined
current_file_id = 00000 + alti_file_offset
next_file_id = current_file_id+alti_file_offset
alti_file_1 = alti_path + '/ALTI%05.0f_t0.log' % current_file_id
alti_file_2 = alti_path + '/ALTI%05.0f_t0.log' % next_file_id

# Vector to store the calculated velocities
y_values = []

# Vector to store the t0 values for each velocity
x_values = []

correlation_times = []

# Function to compute the velocity between
def velocity_between(signal_1, signal_2):
    velocities = []
    if (signal_1.shape != signal_2.shape):
        print('\x1b[6;37;41m' + 'XCOR ERROR: trying to correlate differnt shaped signals' + '\x1b[0m')
        exit()
    
    for slice_index in range(0,altitudes_1_matrix.shape[1]):

        slice_1 = signal_1[:, slice_index]
        slice_2_tiled = np.tile(signal_2[:, slice_index], [3])

        one_d_correlation = signal.correlate(slice_1, slice_2_tiled, mode='same')
        
        plt.figure('1D xcor')
        plt.plot(one_d_correlation)
        plt.show()
        
        location_of_max = np.argmax(one_d_correlation)
        offset = location_of_max - (altitudes_1_matrix.shape[0] / 2)
        velocities.append(offset)

    return abs(np.mean(velocities))

# Check that it is able to open at least the first
if not(os.path.exists(alti_file_1)):
    print('\x1b[6;37;41m' + 'XCOR ERROR: unable to load file at path: ' + alti_file_1  + '\x1b[0m')
    exit()


# Compute the velocity of dunes between altitude files while they exist
while (os.path.exists(alti_file_1) and os.path.exists(alti_file_2)):
    
    # Load the alittudes
    altitudes_1_matrix = np.loadtxt(alti_file_1).T
    altitudes_2_matrix = np.loadtxt(alti_file_2).T
    
    print(altitudes_1_matrix.shape)

    # Cross correlate slices
    start = time.time()
    velocity = velocity_between(altitudes_1_matrix, altitudes_2_matrix)
    exit()
    # Scale the velocity to the time steps between each measurement
    velocity = float(velocity) / alti_file_offset
    correlation_times.append(time.time() - start)
    
    # Print output after each correlation
    # print(str(current_file_id) + '    Velocity: ' + str(velocity) + '    Time: ' + str(correlation_times[-1]))

    y_values.append(velocity)
    x_values.append(np.float64(current_file_id))
    
    # Increment file ids
    current_file_id += alti_file_offset
    next_file_id += alti_file_offset
    alti_file_1 = alti_path + '/ALTI%05.0f_t0.log' % current_file_id
    alti_file_2 = alti_path + '/ALTI%05.0f_t0.log' % next_file_id

# Check that there are at least 2 velocities captured
if len(y_values) < 2:
    print('\x1b[6;37;41m' + 'XCOR ERROR: too few files to correlate' + '\x1b[0m')
    exit()

meta_analysis_file.write('Average correlation time: ' + str(np.mean(correlation_times)))
meta_analysis_file.write('\nTotal correlation time: ' + str(np.sum(correlation_times)))

# Save the velocities for future use
meta_analysis_file.write('\n\nVelocities omiting first, step size of ' + str(alti_file_offset) + ': ' + str(y_values))

############# Try to fit the velocity points #############
def function_to_fit(t, A, b, m, c):
    # Cast all elements of t to integers
    return ((A) * math.e**((t)/(-b))) + ((m) * (t) + (c))

try:
    popt, pcov = curve_fit(function_to_fit, x_values, y_values, p0 = [1,50,-1,1], sigma=np.full(len(y_values), 0.005))
    perr = np.sqrt(np.diag(pcov))
    print(perr)
    A, b, m, c = popt
    calc_y_values = function_to_fit(np.asarray(x_values), float(A), float(b), float(m), float(c))

    # Calculate R^2
    y_bar = np.mean(y_values)
    SS_tot = np.sum((y_values - y_bar)**2)
    SS_res = np.sum((y_values - calc_y_values)**2)
    r_squared = 1 - (SS_res/SS_tot)

    # Save and print out calculated values
    meta_analysis_file.write('\n\n')
    meta_analysis_file.write('Initial state coefficient:    ' + str(A) + '\n')
    meta_analysis_file.write('Rate of approach:        ' + str(b) + '\n')
    meta_analysis_file.write('Slope of quasi-stable state:    ' + str(m) + '\n')
    meta_analysis_file.write('State intercept:        ' + str(c) + '\n')
    meta_analysis_file.write('Initial velocity:        ' + str(y_values[0]) + '\n')
    meta_analysis_file.write('Final velocity:            ' + str(y_values[-1]) + '\n')
    meta_analysis_file.write('R Squared value:        ' + str(r_squared) + '\n')
    meta_analysis_file.write('Error bars:           ' + str(perr) + '\n')
    meta_analysis_file.write('\n\n')

    # Graphing
    plt.figure('plot')
    plt.plot(x_values,y_values)
    plt.plot(x_values, calc_y_values)
    plt.savefig(output_dir + '/' + output_filename + '_velocity_with_fit.png')

except:
    plt.figure('Velocities')
    plt.plot(x_values,y_values)
    plt.savefig(output_dir + '/' + output_filename + '_velocity.png')
    meta_analysis_file.write('\n')
    meta_analysis_file.write('Initial velocity:        ' + str(y_values[0]))
    meta_analysis_file.write('Final velocity:            ' + str(y_values[-1]))
    meta_analysis_file.write('\n')

    print('Error, unable to fit velocities with function.')
