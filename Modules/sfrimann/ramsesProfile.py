#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import numpy as np
import astropy.units as unit
import astropy.constants as cnst

import weights as w

def densityProfile(cell,**kw):
  """
  make a radial density profile from pyramses cell data
  """
  # ------- Handle keywords
  rmin = kw.get('rmin',None)
  rmax = kw.get('rmax',None)
  nr   = kw.get('nr',100)
  
  # ------- set rmin and rmax
  if rmin is None:
    rmin = np.sqrt(3*(cell.dx.min()/2)**2).to('AU').value
  
  if rmax is None:
    xmi, xma = cell.xlim
    ymi, yma = cell.ylim
    zmi, zma = cell.zlim
    rmax = min(abs(xmi),xma,abs(ymi),yma,abs(zmi),zma).to('AU').value
  
  # ------- set up radial grid
  r   = np.logspace(np.log10(rmin),np.log10(rmax),nr)
  re, mass, nbin = [], [], []
  
  # ------- read cell quantities
  cellRadius = cell.radius.to('AU').value
  cellDensity = cell.density.cgs.value
  cellVolume = cell.volume.cgs.value
  cellMass   = cell.mass.cgs.value
  
  # ------- loop over radial cells taking all values inside a given radius
  for i in range(r.size):
    index, = np.where(cellRadius <= r[i])
    if index.size == 0:
      continue
    ivolume = cellVolume[index].sum()
    equivalent_radius = (3/4./np.pi*ivolume)**(1/3) * unit.cm.to('AU')
    if len(re) > 0:
      if equivalent_radius == re[-1]:
        continue
    
    re.append(equivalent_radius)
    imass  = cellMass[index].sum()
    mass.append(imass)
    nbin.append(index.size)
  
  # ------- convert lists to numpy arrays
  re = np.asarray(re)
  nbin = np.asarray(nbin)
  mass = np.asarray(mass)
  
  nbin = np.hstack((nbin[0],np.diff(nbin))) # cumulative nbin to shell nbin
  re   = np.hstack((0.,re)) # add zero to edges
  
  rc = 0.5*(re[1:] + re[:-1]) # central coordinates
  
  density = np.empty(rc.size,dtype=np.float64)
  sdensity = np.empty((2,rc.size),dtype=np.float64)
  nb = 0
  isort = np.argsort(cellRadius)
  for i,n in enumerate(nbin):
    dens = cellDensity[isort][nb:(nb+n)]
    weights = cellVolume[isort][nb:(nb+n)]
    density[i] = w.wmedian(dens,weights=weights)
    if dens.size > 1:
      sdensity[:,i] = w.wpercentile(dens,[0.159,0.841],weights=weights)
    else:
      sdensity[:,i] = np.nan
    nb += n
  sdensity[0] = density - sdensity[0]
  sdensity[1] = sdensity[1] - density
  
  dic = dict(radius = rc,
             mass   = mass,
             nbin   = nbin,
             density = density,
             sdensity = sdensity)
  
  return dic