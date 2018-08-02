import math
import time
import glob
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit
import rescal_utilities as ru

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

# Get the timestamp for filenames
timestamp_string = str(time.time())

# Output files
log_file = open( output_dir + '/' + output_filename + '_' + timestamp_string + '_log.txt','w')
analysis_file = open(output_dir + '/' + output_filename + '_' + timestamp_string +  '_analysis.txt', 'w')

# Write the parameters from the par file
# get par file
print(alti_path + '/*.par')
par_file_name = glob.glob(alti_path + '/*.par')[0]

params = ru.Parameters()
try:
    params.read(par_file_name)
except:
    print("Parameter file not found at: {}".format(par_file_path))
    params = None

try:
    #Write summary data
    if params:
        lam_s = params.get('Lambda_S')
        tau = params.get('Tau_min')
        Ava = params.get('Ava_angle')
        h = params.get('H')
        d = params.get('D')
        l = params.get('L')
        meta_analysis_file.write('Parameter info:\nLambda S: {}\nTau_min: {}\nAva_angle: {}\nHeight: {} Depth: {} Length: {}\nFrequencies logged: {}\n\n'.format(lam_s,tau,Ava,h,d,l,freqs))
except:
    log_file.write('\nFailed to write parameters to output file\n')
    print('\x1b[6;37;43m' + 'XCOR WARNING: Falide to write parameters to output file' + '\x1b[0m')

# Add file offset to skip the first file as cross correlation with flat surface is ill defined
current_file_id = 00000 + alti_file_offset
next_file_id = current_file_id+alti_file_offset
alti_file_1 = alti_path + '/ALTI%05.0f_t0.data' % current_file_id
alti_file_2 = alti_path + '/ALTI%05.0f_t0.data' % next_file_id

# Vector to store the calculated velocities
y_values = []

# Vector to store the t0 values for each velocity
x_values = []

correlation_times = []

# Function to compute the velocity between
def velocity_between(signal_1, signal_2):
    velocities = []
    
    for slice_index in range(0,altitudes_1_matrix.shape[1]):

        slice_1 = signal_1[:, slice_index]
        slice_2_tiled = np.tile(signal_2[:, slice_index], [3])

        one_d_correlation = signal.correlate(slice_1, slice_2_tiled, mode='same')
        
        #plt.figure('1D xcor')
        #plt.plot(one_d_correlation)
        #plt.show()
        
        location_of_max = np.argmax(one_d_correlation)
        offset = location_of_max - (altitudes_1_matrix.shape[0] / 2)
        velocities.append(offset)

    return abs(np.mean(velocities))

# Check that it is able to open at least the first
if not(os.path.exists(alti_file_1)):
    log_file.write('\nERROR: unable to load file at path: ' + alti_file_1 + '\n')
    print('\x1b[6;37;41m' + 'XCOR ERROR: unable to load file at path: ' + alti_file_1  + '\x1b[0m')
    exit()


# Compute the velocity of dunes between altitude files while they exist
while (os.path.exists(alti_file_1) and os.path.exists(alti_file_2)):
    
    # Load the alittudes
    altitudes_1_matrix = np.loadtxt(alti_file_1).T
    altitudes_2_matrix = np.loadtxt(alti_file_2).T
    
    if (altitudes_1_matrix.shape != altitudes_2_matrix.shape):
        log_file.write('\nWARNING: Trying to correlate differently shaped signals, terminating early\n')
        print('\x1b[6;37;43m' + 'XCOR WARNING: trying to correlate differnt shaped signals early termination' + '\x1b[0m')
        break
    
    print('Analyzing: ' + str(current_file_id))

    # Cross correlate slices
    start = time.time()
    velocity = velocity_between(altitudes_1_matrix, altitudes_2_matrix)
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
    alti_file_1 = alti_path + '/ALTI%05.0f_t0.data' % current_file_id
    alti_file_2 = alti_path + '/ALTI%05.0f_t0.data' % next_file_id

# Check that there are at least 2 velocities captured
if len(y_values) < 2:
    log_file.write('\nERROR: Too few output files to correlate\n')
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
    
    log_file.write('\nERROR: Unable to fit velocities with function\n')
    print('\x1b[6;37;41m' + 'XCOR ERROR: unable to fit velocities with function' + '\x1b[0m')
