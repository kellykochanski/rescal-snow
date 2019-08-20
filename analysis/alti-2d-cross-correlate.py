"""
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
"""


################################################################
# Use cross-correlation to find the offset between two signals #
################################################################
import sys
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib import cm
from math import pi
from mpl_toolkits.mplot3d import Axes3D
import time

## Set up variables for signal generation
#signal_length = 1000.0
#signal_width = 20.0
#TAU = 2*pi
#resolution = 1.0
#signal_1_values = np.arange(0,signal_length,resolution)
#signal_2_values = np.arange(0,signal_length,resolution)
#
## Generate signals and turn them into surfaces
#
#signal1_harmonic = np.cos((signal_1_values/100)*20)
#signal1 = signal1_harmonic + np.cos(TAU*(signal_1_values/100))
##plt.plot(signal1_harmonic)
##plt.show
#signal2 = np.cos(TAU*(signal_2_values/100))
#signal1_surface = np.reshape(np.repeat(signal1, signal_width), newshape=(len(signal1), -1))
#signal2_surface = np.reshape(np.repeat(signal2, signal_width), newshape=(len(signal2), -1))
#
## Add cross dimensional sine wave to surfaces
##signal1_surface += signal1.T
##signal2_surface += signal2.T

# Load altitudes from files

if len(sys.argv) < 3:
    print("Please provide the numbers of the altitude files to analize e.g. 0004 0008")
    exit()
else:
    altitude_number = sys.argv[1]
    alti_file_1 = 'ALTI' + str(sys.argv[1]) + '_t0.log'
    alti_file_2 = 'ALTI' + str(sys.argv[2]) + '_t0.log'

print("loading: " + str(alti_file_1) + " " + str(alti_file_2))

#alti_file_1 = 'ALTI0285.data'
#alti_file_2 = 'ALTI0300.data'
altitudes_1_matrix = np.loadtxt(alti_file_1).T
altitudes_2_matrix = np.loadtxt(alti_file_2).T

# Set width and length of cells
signal_length = altitudes_1_matrix.shape[0]
signal_width = altitudes_1_matrix.shape[1]
print('signal width:' + str(signal_width))
print('signal length:' + str(signal_length))


# Compute Cross Correlation
start = time.time()
xcorr = signal.correlate2d(altitudes_1_matrix, altitudes_2_matrix, mode='same')
print('Time to cross correlate: ' + str(time.time() - start) + ' seconds')

# Find the indices of the max
max_width = np.argmax(np.max(xcorr, axis=0))+1
max_length = np.argmax(np.max(xcorr, axis=1))+1

print('Max location X:' + str(max_width) + ' Y:' + str(max_length))

# Scale the found location to it's distance from the center
offset_width = max_width - signal_width/2
offset_length = max_length - signal_length/2

print('The offset x (width) between the two signals is ' + str(offset_width) + ' cells.')
print('The offset y (length) between the two signals is ' + str(offset_length) + ' cells.')

#################################################################
############################ GRAPHING ###########################
#################################################################
#
## Create mesh grid for graphing
#X, Y = np.meshgrid(np.arange(0, signal_width), np.arange(0, signal_length))
#
### Test plot signal 1
##plt.figure("Test")
##plt.imshow(signal1_surface)
##plt.colorbar()
##
### Test plot alti
##plt.figure("alti1")
##plt.imshow(altitudes_1_matrix)
##plt.colorbar()
##plt.figure("alti2")
##plt.imshow(altitudes_2_matrix)
##plt.colorbar()
#
### Plot first signal
##fig = plt.figure('Signal 1')
##ax = fig.gca(projection='3d')
##surf = ax.plot_surface(X, Y, signal1_surface, cmap=cm.coolwarm)
##fig.colorbar(surf, shrink=0.5, aspect=5)
#
#
## Test plot alti
#X2, Y2 = np.meshgrid(np.arange(0, altitudes_1_matrix.shape[1]), np.arange(0, altitudes_1_matrix.shape[0]))
#fig = plt.figure('Alti 1')
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X2, Y2, altitudes_1_matrix, cmap=cm.coolwarm)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#fig = plt.figure('Alti 2')
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X2, Y2, altitudes_2_matrix, cmap=cm.coolwarm)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
## Test plot correlation
#fig = plt.figure('Cross Corrleation')
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X2, Y2, xcorr, cmap=cm.coolwarm)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
### Plot second signal
##fig = plt.figure('Signal 2')
##ax = fig.gca(projection='3d')
##surf = ax.plot_surface(X, Y, signal2_surface, cmap=cm.coolwarm)
##fig.colorbar(surf, shrink=0.5, aspect=5)
##
### Plot the correlation
##fig = plt.figure('Cross Correlation')
##ax = fig.gca(projection='3d')
##surf = ax.plot_surface(X, Y, xcorr, linewidth=0, cmap=cm.coolwarm)
##fig.colorbar(surf, shrink=0.5, aspect=5)
#
#plt.show()

