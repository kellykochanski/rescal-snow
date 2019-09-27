import sys
import os

sys.path.append('utilities')
import heightmap
import cellspace # import these to check quickly for syntax errors
import datarun

filenames = ['out/'+f for f in os.listdir('out') if f[:4]=='ALTI']
filenames = [f for f in filenames if f[-4:] == '.log']
filenames.sort()

for name in filenames:
    hm = heightmap.HeightMap(name)
    hm.save_color_map(name[:-4] + '_image.png')

