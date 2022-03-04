#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division

import numpy as np
import matplotlib.pyplot as mp
import pyramses.amr as new
import pyradmc3d.amr as old

# build amr grid for testing purposes

def build_grid():
  
#   ngridx, ngridy, ngridz = 3, 3, 3
#   
#   edgex, edgey, edgez = np.linspace(0,1,(ngridx+1)), np.linspace(0,1,(ngridy+1)), np.linspace(0,1,(ngridz+1))
#   
#   xx, yy, zz = 0.5*(edgex[1:] + edgex[:-1]), 0.5*(edgey[1:] + edgey[:-1]), 0.5*(edgez[1:] + edgez[:-1])
#   dxx, dyy, dzz = np.diff(edgex), np.diff(edgey), np.diff(edgez)
#   
#   z, y, x = np.meshgrid(zz,yy,xx,indexing='ij')
#   dz, dy, dx = np.meshgrid(dzz,dyy,dxx,indexing='ij')
#   
#   dz, dy, dx = dz.ravel(), dy.ravel(), dx.ravel()
#   
#   pos = np.vstack((x.ravel(),y.ravel(),z.ravel())).T
#   
#   amrtree = np.zeros(dx.size,dtype=np.int8)
#   
#   pos,dx,amrtree = refine(pos,dx,amrtree,13)
  
  ngridx, ngridy, ngridz = 5, 5, 5
  
  edgex, edgey, edgez = np.linspace(0,1,(ngridx+1)), np.linspace(0,1,(ngridy+1)), np.linspace(0,1,(ngridz+1))
  
  xx, yy, zz = 0.5*(edgex[1:] + edgex[:-1]), 0.5*(edgey[1:] + edgey[:-1]), 0.5*(edgez[1:] + edgez[:-1])
  dxx, dyy, dzz = np.diff(edgex), np.diff(edgey), np.diff(edgez)
  
  z, y, x = np.meshgrid(zz,yy,xx,indexing='ij')
  dz, dy, dx = np.meshgrid(dzz,dyy,dxx,indexing='ij')
  
  dz, dy, dx = dz.ravel(), dy.ravel(), dx.ravel()
  
  pos = np.vstack((x.ravel(),y.ravel(),z.ravel())).T
  
  amrtree = np.zeros(dx.size,dtype=np.int8)
  
  for i in range(300):
    index = np.random.randint(0,dx.size)
    pos,dx,amrtree = refine(pos,dx,amrtree,index)
  
  return pos,dx,amrtree

def refine(pos,dx,amrtree,index):
  
  x_refiner = np.array([-0.25,0.25,-0.25,0.25,-0.25,0.25,-0.25,0.25])
  y_refiner = np.array([-0.25,-0.25,0.25,0.25,-0.25,-0.25,0.25,0.25])
  z_refiner = np.array([-0.25,-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25])
  
  dx_refiner = np.ones(8)/2.
  
  refiner = np.vstack((x_refiner,y_refiner,z_refiner)).T
  amr_refiner = np.array([1,0,0,0,0,0,0,0,0],dtype=np.int8)
  
  # amrtree indexing
  indarr = np.arange(amrtree.size)
  leaf_index = amrtree == 0
  amrtree_index = indarr[leaf_index][index]
  
  if index == 0:
    pos = np.vstack((refiner*dx[0]+pos[0,:],pos[1:,:]))
    dx  = np.hstack((dx_refiner*dx[0],dx[1:]))
  elif index == (len(dx)-1):
    pos = np.vstack((pos[:-1,:],refiner*dx[-1]+pos[-1,:]))
    dx  = np.hstack((dx[:-1],dx_refiner*dx[-1]))
  else:
    pos = np.vstack((pos[:index,:],refiner*dx[index]+pos[index,:],pos[(index+1):,:]))
    dx  = np.hstack((dx[:index],dx_refiner*dx[index],dx[(index+1):]))
  
  if amrtree_index == 0:
    amrtree = np.hstack((amr_refiner,amrtree[1:]))
  elif amrtree_index == (amrtree.size -1):
    amrtree = np.hstack((amrtree[:-1],amr_refiner))
  else:
    amrtree = np.hstack((amrtree[:amrtree_index],amr_refiner,amrtree[(amrtree_index+1):]))
  
  return pos,dx,amrtree

def shuffle(pos,dx):
  
  ind = np.arange(dx.size)
  np.random.shuffle(ind)
  return pos[ind], dx[ind]

if __name__ == '__main__':
  
  pos, dx, amrtree = build_grid()
  
  spos, sdx = shuffle(pos,dx)
  
  nind, namrtree = new.amrsort(spos,sdx, strip_edge=True)
  oind, oamrtree = old.amrsort(spos,sdx)
  
  print((nind == oind).all())
  print((namrtree == oamrtree).all())
  print(namrtree.size,namrtree.sum(),amrtree.sum())
  