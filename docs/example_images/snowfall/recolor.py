## This script recolors 'ALTI00xxx.log' files with matplotlib color schemes
# It allows users to change rescal-snow visualizations when the default colors are not useful.
# This particular script is configured to work with the example output in docs/example_images/snowfall

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import os

dirname = '.'
for filename in os.listdir(dirname):
	if filename[:4] == 'ALTI':
		if filename[-8:] == '0_t0.log':
			raw = pd.read_csv(os.path.join(dirname, filename))
			data = [[int(num) for num in filter(None, value[0].split(' '))] for value in raw.values]
			vmx = 20
			if vmx < max(data):
				print("Warning: some values saturated")
			plt.imshow(data, cmap='magma', vmin=2, vmax=vmx)
			plt.colorbar()
			plt.savefig(filename[:-4]+'_recolored.png')
			print('Recolored ' + filename)
			plt.close()
