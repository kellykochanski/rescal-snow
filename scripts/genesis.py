#!/usr/bin/python

## Generation of a 3d cellular space (raw binary format)

## From Rescal v1.6  (2016) -- not tested with rescal-snow (2019)

import sys

#model DUN
from cellspace import *

#model AVA
#from cellspace_ava import *

def usage():
  print "genesis.py - " + MOD_NAME + " model"
  print "Usage: genesis.py <H> <L> <D>"
  print "   or: genesis.py -f <PARAMETERS_FILE>"
  print "  H\theight"
  print "  L\tlength"
  print "  D\tdepth"
  quit()

def read_params(file):
  """Read parameters file"""
  print "genesis.py - " + MOD_NAME + " model"
  print "reading parameters file: " + file
  fpar = open(file,'r')
  while True:
    line = fpar.readline()
    if line=="":
      break
    if len(line)==1 or line[0]=='#':
      continue
    l=line.replace(' ','').split('=')
    par=l[0]
    val=l[1].replace('\n','')
    #print par + " = " + val
    if par=="H":
      H=int(val)
    elif par=="L":
      L=int(val)
    elif par=="D":
      D=int(val)
    elif par=="Bin_file":
      bin_filename=val;
    elif par=="Csp_file":
      print "Failed to create CSP file: CSP format not supported yet ..."
      sys.exit()
    #else:
    #  print "unkown parameter: " + par
  return [H,L,D,bin_filename]

## parse arguments
if len(sys.argv)==3 and sys.argv[1]=="-f":
  H,L,D,bin_filename = read_params(sys.argv[2])
elif len(sys.argv)==4:
  ## read dimensions
  H=int(sys.argv[1])
  L=int(sys.argv[2])
  D=int(sys.argv[3])
  bin_filename = "%s_%d_%d_%d.bin" % (MOD_NAME,H,L,D)
else:
  usage()

print "H =",H
print "L =",L
print "D =",D

## generating cellular space
cs=cellspace(H,L,D)

## writing binary file
print "writing binary file: " + bin_filename
fbin=open(bin_filename,"wb")
print fbin
fbin.write(cs)
fbin.close()
print "done"

