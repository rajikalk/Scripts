#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
#from numpy.polynomial.hermite import hermval #numpy 1.6.1

# ---- Error definitions ----
class InputError(Exception):
  pass

class Gauss:
  """Class to create one-dimensional Gaussian functions"""
  def __init__(self,*args):
    if len(args) == 4:
      self.c = np.array(args,dtype=np.float)
    elif len(args) == 3:
      self.c = np.array(args,dtype=np.float)
      self.c = np.append(self.c,0.)
    elif len(args) == 2:
      self.c = np.array(args,dtype=np.float)
      self.c = np.append(self.c,[0.,0.])
    elif len(args) == 1:
      self.c = np.array(args,dtype=np.float)
      self.c = np.append(self.c,[1.,0.,0.])
    else:
      self.c = np.array([1.,1.,0.,0.],dtype=np.float)
    self.sigma, self.amp, self.mean, self.offset = self.c

  def __call__(self,x,deriv=0):
    argument = (x-self.mean)/(np.sqrt(2)*self.sigma)
    if deriv == 0:
      return self.amp*np.exp(-argument**2) + self.offset
    elif deriv == 1:
      pre = -2*argument
    elif deriv == 2:
      pre = 4*argument**2 - 2
    else:
      raise IOError('deriv must be a 0, 1, or 2')
    
    return self.amp*pre*np.exp(-argument**2)

    #The below will only work with numpy 1.6.1 and above
    #  n = zeros(deriv+1)
    #  n[-1] = 1 if deriv%2 == 0 else -1

    #  return hermval(argument,n)*exp(-argument**2)