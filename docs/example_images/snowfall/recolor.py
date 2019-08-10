#!/usr/bin/env python3

## This script recolors 'ALTI00xxx.log' files with matplotlib color schemes
# It allows users to change rescal-snow visualizations when the default colors are not useful.
# This particular script is configured to work with the example output in docs/example_images/snowfall

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import os

dirname = '..'
for filename in os.listdir(dirname):
  if filename.startswith('ALTI') and filename.endswith('0_t0.log'):
    data = np.loadtxt(os.path.join(dirname, filename))
    vmx = 20
    if vmx < data.max():
      print("Warning: some values saturated")
    plt.imshow(data, cmap='magma', vmin=2, vmax=vmx)
    plt.colorbar()
    plt.savefig(filename[:-4]+'_recolored.png')
    print('Recolored ' + filename)
    plt.close()
