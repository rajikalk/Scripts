#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def two_power(theta,xdata):
  """
  Simple power law model
  """
  power1, power2, xloc = theta
  return np.where(xdata <= xloc,(xdata/xloc)**(-power1),(xdata/xloc)**(-power2))

def offset(theta,r,z,scales):
  """
  constant offset
  """
  c0 = 10**theta if scales == 'log' else theta
  array = np.empty(r.shape) ; array.fill(float(c0))
  return array

def disk_model_classic(theta,r,z,scales):
  """
  Classical disk model
  """
  sigma0, h0, rpower, hpower = np.where(scales == 'log',10**theta,theta)
  r0 = 100. # sample distance in AU
  hr =  h0 * (r/r0)**hpower # pressure scale height
  return sigma0 * (r/r0)**(-rpower)/np.sqrt(2*np.pi)/hr * np.exp(-.5*(z/hr)**2)

def disk_model_plane(theta,r,z,scales):
  """
  disk model with central density being determined from the plane
  """
  sigma0, h0, rloc, rpower1, rpower2, hpower = np.where(scales == 'log',10**theta,theta)
  r0 = 100. # sample distance in AU
  hr =  h0 * (r/r0)**hpower
#   hr = np.where(hr < 20/1.15, 20/1.15, hr)
  return sigma0 * two_power([rpower1,rpower2,rloc],r) * np.exp(-.5*(z/hr)**2)

def disk_model_with_cutoff(theta,r,z,scales):
  """
  disk model with cutoff
  """
  sigma0, h0, rloc, rpower1, rpower2, hpower, rout, fac = np.where(scales == 'log',10**theta,theta)
  distance = np.sqrt(r**2 + z**2)
  
  r0 = 100. # AU
  hr =  h0 * (r/r0)**hpower
  return sigma0 * two_power([rpower1,rpower2,rloc],r) * np.exp(-.5*(z/hr)**2) * np.where(distance < rout,1.,np.exp(-fac*(distance-rout)))

def disk_model_only_cutoff(theta,r,z,scales,sigmar,hr):
  """
  disk model with central density being determined from the plane
  """
  rcut = 10**theta if scales == 'log' else theta
  diskmodel = sigmar * np.exp(-.5*(z/hr)**2)
  return np.where(r <= rcut,diskmodel,diskmodel*np.exp(-.5*r/10*2))

# ------- envelope models
def envelope_model_sphere(theta,r,z,scales):
  """
  spherical envelope model
  """
  sigma0 = 10**theta[0] if scales[0] == 'log' else theta[0]
  power  = 10**theta[1] if scales[1] == 'log' else theta[1]

  power = theta[0] ; sigma0 = 10**theta[1]
  distance = np.sqrt(r**2 + z**2)
  r0 = 1000.
  return sigma0*(distance/r0)**(-power)