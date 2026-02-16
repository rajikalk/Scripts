#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division

import numpy as np

def extrapolate1d(x,xp,fp,left=None,right=None):
  """
  Extrapolation routine. Piecewise linear interpolation with extrapolation as specified
  by user
  Input:
    x     : Interpolated x values. Must be monotonically increasing
    xp    : x-coordinates of data points
    fp    : y-values of data points
  Keywords:
    left  : Directions for left extrapolation. None is no extrapolation (an error is thrown
            if values exceed range), 'constant' is constant extrapolation. 'linear' gives
            linear extrapolation
    right : Like left but for right values
  Returns:
        f : Interpolated y values
  """
  
  x = np.asarray(x)
  
  assert np.all(np.diff(xp)), "xp must be monotonically increasing"
  assert np.all(np.isfinite(fp)), "All elements in fp must be finite"
  
  f = np.interp(x,xp,fp)
  
  rindex = x > xp[-1]
  if rindex.sum() > 0: # Extrapolation needed
    if right is None:
      raise ValueError("Right interpolation not specified but {} points are beyond right edge".format(rindex.sum()))
    elif right == 'linear':
      # calculate linear extrapolation from last two grid points
      a = (fp[-1] - fp[-2])/(xp[-1] - xp[-2])
      b = fp[-1] - a*xp[-1]
      f[rindex] = a*x[rindex] + b
    elif right == 'constant':
      pass
    else:
      raise ValueError("{} type extrapolation not recognized".format(right))
  
  lindex = x < xp[0]
  if lindex.sum() > 0: # Extrapolation needed
    if left is None:
      raise ValueError("Left interpolation not specified but {} points are beyond right edge".format(lindex.sum()))
    elif left == 'linear':
      # calculate linear extrapolation from last two grid points
      a = (fp[1] - fp[0])/(xp[1] - xp[0])
      b = fp[0] - a*xp[0]
      f[lindex] = a*x[lindex] + b
    elif left == 'constant':
      pass
    else:
      raise ValueError("{} type extrapolation not recognized".format(left))
  
  return f

def calculate_dm(star,starsdat,dt=1):
  """
  Function for calculating sink accretion rates
  star is single snapshot sink data
  starsdat is sink data from multiple snapshot or stars.dat file, gathered into a list
  dt is the accretion window in indices. Default is to calculate accretion over one index
  returns the mass difference over the accretion window and the accretion window itself
  in code units
  """
  
  snapshot_times = np.array([s['snapshot_time'] for s in starsdat])
  
  nstar = len(star['m'])
  
  idx        = extrapolate1d(star['snapshot_time'],snapshot_times,np.arange(snapshot_times.size),left='linear',right='linear')
  delta_time = extrapolate1d([idx-dt/2,idx+dt/2],np.arange(snapshot_times.size),snapshot_times,left='linear',right='linear')
  delta_time = delta_time[1]-delta_time[0]
  indices = np.arange(np.floor(idx-dt),np.ceil(idx+dt)+0.1,dtype=np.int)
  
  masses  = np.zeros((nstar,indices.size))
  
  # Fill masses array
  for i,index in enumerate(indices):
    if index < 0:
      masses[:,i] = 0.
    elif index >= len(starsdat):
      masses[:,i] = np.abs(starsdat[-1]['m'][:nstar])
    else:
      s = starsdat[index]['m']
      if len(starsdat[index]['m']) > nstar:
        masses[:,i] = np.abs(starsdat[index]['m'][:nstar])
      elif len(starsdat[index]['m']) < nstar:
        s = np.zeros(nstar,dtype=np.float)
        s[:starsdat[index]['m'].size] = np.abs(starsdat[index]['m'])
        masses[:,i] = s
      else:
        masses[:,i] = np.abs(starsdat[index]['m'])
    
  dm = np.zeros(nstar,dtype=np.float)
  for i in range(nstar):
    imass = np.interp([idx-dt/2,idx+dt/2],indices,masses[i,:])
    dm[i] = imass[1]-imass[0]
  
  return dm, delta_time