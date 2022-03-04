#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import numpy as np

from scipy.spatial.distance import pdist

# personal modules
from cyclic import dist
from binned_statistic import binned_statistic

def sphereBoxVolume(r,L=1):
  """
  Volume of a sphere with radius r intersected by a cube with side length L.
  The centre of the cube and sphere coindice and only the spherical volume inside the cube
  is acounted for
  """
  
  if r < 0:
    return np.nan
  elif r == 0:
    return 0
  elif r <= L/2: # sphere does not intersect with cube
    return 4/3. * np.pi * r**3
  elif r < np.sqrt(2)/2. * L: # six spherical caps escape the cube faces
    return 4/3. * np.pi * r**3 - np.pi/4*(2*r - L)**2*(4*r + L)
  elif r < np.sqrt(3)/2 * L: # the spherical caps overflow
    # solved by calculating the region integral in mathematica
    # syntax:
    # Integrate[If[x^2 + y^2 + z^2 < r, 1, 0], {x, -1, 1}, {y, -1, 1}, {z, -1, 1}, Assumptions -> {r > 2, r < 3}]
    # The box size was fixed to L=2 and r is assumed to be squared so the units have to be
    # scaled before applying the result
    r2 = (2*r/L)**2
    
    result = 1/3*(-6*np.pi+24*np.sqrt(-2+r2)+6*np.pi*r2+3*np.arctan(1/np.sqrt(-2+r2))+9*r2*np.arctan(1/np.sqrt(-2+r2))+12*r2**(3/2)*np.arcsin(1/(1-r2))-3*np.arcsin(1/np.sqrt(-1+r2))+9*r2*np.arcsin(1/np.sqrt(-1+r2))-np.arcsin(1/np.sqrt(-1+r2))+3*r2*np.arcsin(1/np.sqrt(-1+r2))+np.arctan(1/np.sqrt(-2+r2))+3*r2*np.arctan(1/np.sqrt(-2+r2))+24*np.arctan(np.sqrt(-2+r2))-48*r2*np.arctan(np.sqrt(-2+r2))+16*r2**(3/2)*np.arctan(np.sqrt((-2+r2)/r2))-4*r2**(3/2)*np.arctan(1/np.sqrt((-2+r2)*r2)))
    
    return (L/2)**3*result
  else: # the inside of the sphere fills the cube
    return L**3

def sphereBoxShellVolume(r,L=1):
  """
  Volume of shells
  """
  result = np.empty(r.size-1,dtype=np.float64)
  for i in range(1,r.size):
    result[i-1] = sphereBoxVolume(r[i],L=L) - sphereBoxVolume(r[i-1],L=L)
  
  return result

def circleBoxArea(r,L=1):
  """
  Two dimensional version of sphereBoxVolume
  """
  if r < 0:
    return np.nan
  elif r == 0:
    return 0
  elif r <= L/2:
    return np.pi*r**2
  elif r <= np.sqrt(2)/2 * L:
    return np.pi*r**2 - 4*r**2*np.arccos(L/(2*r)) + 4*(L/2)**2*np.sqrt((2*r/L)**2 - 1)
  else:
    return L**2

def circleBoxShellArea(r,L=1):
  """
  Area of rings
  """
  result = np.empty(r.size-1,dtype=np.float64)
  for i in range(1,r.size):
    result[i-1] = circleBoxArea(r[i],L=L) - circleBoxArea(r[i-1],L=L)
  return result

def twoPointCorrelationFunction(positions,**kw):
  """
  Calculate the two point correlation function for a number of point particles
  """
  # ------- Handle keywords
  Lbox   = kw.get('Lbox',1)
  cyclic = kw.get('cyclic',True)
  rmin   = kw.get('rmin',None)
  rmax   = kw.get('rmax',None)
  nbin   = kw.get('nbin',10)
  xlog   = kw.get('xlog',True)
  ylog   = kw.get('ylog',True)
  
  if cyclic is False:
    raise NotImplementedError("Not yet implemented this feature")
  
  Nstar, ndim = positions.shape
  
  separations = pdist(positions,lambda x1,x2: dist(x1,x2,period=Lbox,cyclic=cyclic))
  
  # ------- set edges
  rmin = separations.min() if rmin is None else rmin
  rmax = separations.max() if rmax is None else rmax
  if xlog is True:
    redges = np.logspace(np.log10(rmin),np.log10(rmax),nbin+1)
  else:
    redges = np.linspace(rmin,rmax,nbin+1)
  rcentre = 0.5*(redges[1:] + redges[:-1])
  
  count,_ = np.histogram(separations,redges)
  scount = np.sqrt(count)
  
  if ndim == 2:
    nstar = Nstar/Lbox**2
    dV    = circleBoxShellArea(redges,L=Lbox)
  elif ndim == 3:
    nstar = Nstar/Lbox**3
    dV    = sphereBoxShellVolume(redges,L=Lbox)
  
  # pair correlation function
  g  = count/(Nstar*nstar*dV/2)
  sg = scount/(Nstar*nstar*dV/2)
  
  return rcentre, g, sg