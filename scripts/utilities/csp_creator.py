#!/usr/bin/python

from numpy import *

## DUN model
MOD_NAME="DUN"

## predefined cell types
GR=0
EAUC=3
DUM=6
IN=7
OUT=8

###############################################################
def create_simple_space(H,L,D):
  """creates a simple cellular space without sand"""

  ## create cellular space
  cs=zeros((D,H,L),int8)

  ## air volume
  cs[:,:,:]=EAUC

  ## ceiling
  cs[:,0,:]=DUM

  ## floor
  cs[:,H-1,:]=DUM
  
  return cs


###############################################################
def sand_layer(H,L,D):
  """generates a cellular space with uniform sand layer"""
  cs=create_simple_space(H,L,D)
  
  print "generating uniform sand layer"
  
  ## sand layer
  cs[:,H-2,:]=GR
  
  return cs

###############################################################
def sand_block(H,L,D):
  """generates a cellular space with a block of sand"""
  cs=create_simple_space(H,L,D)  
  
  print "generating a block of sand"
  
  ## sand block 20x20x20
  l=10
  L2=L/2
  D2=D/2
  cs[L2-l:L2+l,H-2-2*l:H-2,D2-l:D2+l]=GR
  
  return cs

###############################################################
def cellspace(H,L,D):
  """generates cellular space"""
  cs=sand_layer(H,L,D)
  #cs=sand_block(H,L,D)
  return cs

