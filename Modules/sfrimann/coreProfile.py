#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import numpy as np

from scipy.optimize import leastsq

from binned_statistic import binned_statistic


def powerLaw1DFloor(xdata,*params):
  """
  power law envelope model of the form
    density = a * xdata**(-b) + c
  """
  if len(params) == 1:
    a, b, c = params[0]
  else:
    a, b, c = params
  return a * xdata**(-b) + c

def fitCoreProfile(radius,density,**kw):
  """
  fit core profile
  """
  
  # ------- Handle keywords
  volumes  = kw.get('volumes',None) # array of volumes. Useful for cell data
  masses   = kw.get('masses',None)  # array of particle masses. Useful for particle data
  rmin     = kw.get('rmin',0)       # minimum r value used in fit
#   relative = kw.get('relative',True)
  nbin     = kw.get('nbin',51)       # 0 means no binning
  xlog     = kw.get('xlog',True)
  ylog     = kw.get('ylog',True)
  
  r = np.log10(radius) if xlog is True else radius
  d = np.log10(density) if ylog is True else density
  
  if nbin > 0:
    range = (np.log10(max(rmin,radius.min())),np.log10(radius.max())) if xlog is True else (max(rmin,radius.min()),radius.max())
    xc, y, yl, yu, ys, count = binned_statistic(r,d,range=range,nbin=nbin)
  else:
    raise NotImplementedError('Not yet implemented non-binned data')
  
  index = np.isfinite(1/ys) # only use bins containing more than one data point
  
  modelfun = powerLaw1DFloor
  
  # ------- set error function
  if xlog is True:
    if ylog is True:
      errfun = lambda x: ((y - np.log10(modelfun(10**xc,x)))/ys)[index]
    else:
      errfun = lambda x: ((y - modelfun(10**xc,x))/ys)[index]
  else:
    if ylog is True:
      errfun = lambda x: ((y - np.log10(modelfun(xc,x)))/ys)[index]
    else:
      errfun = lambda x: ((y - modelfun(xc,x))/ys)[index]
  
  # starting guesses
  p0 = [10**(-12),2,10**(-19)]
  
  # do fit
  pfit, pcov, infodict, errmsg, success = leastsq(errfun,p0,full_output=1)
  
  if not success in [1,2,3,4]:
    raise ValueError("leastsq did not converge to a solution")
  
  if (len(y) > len(p0)) and pcov is not None:
    s_sq = (errfun(pfit)**2).sum()/(len(y)-len(p0))
    pcov = pcov * s_sq
  else:
    pcov = np.inf
  
  return pfit, pcov