#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import numpy as np

def mad(x,normed=True,axis=None):
  """
  Median absolute deviation
  if normed is True it is put on std scale
  """
  factor = 1.4826 if normed is True else 1.
  return np.median(np.abs(x - np.median(x,axis=axis)),axis=axis) * factor

def wmad(x,normed=True,axis=None,weights=None):
  """
  Weighted median aboslute deviation
  """
  if weights is None:
    return mad(x,normed=normed,axis=axis)
  
  factor = 1.4826 if normed is True else 1.
  return wmedian(np.abs(x - wmedian(x,axis=axis,weights=weights)),axis=axis,weights=weights)*factor
  
def wmedian(x,weights=None,axis=None):
  """
  Calculate weighted median
  """
  if weights is None:
    return np.median(x,axis=axis)
  else:
    return wpercentile(x,0.5,weights=weights,axis=axis)
#   ix = np.argsort(x)
#   xx, ww = x[ix], w[ix]
#   nx = len(xx)
#   Sn = np.cumsum(ww)
#   pn = 1./Sn[-1] * (Sn - ww/2) # weighted percentiles
# #   print pn
#   index = np.searchsorted(pn,0.5)
# #   print index
#   if index == 0:
#     return xx[index]
#   elif index == nx:
#     return xx[index-1]
#   elif np.abs(pn[index]-0.5) < np.abs(pn[index-1]-0.5):
#     return xx[index]
#   elif np.abs(pn[index]-0.5) > np.abs(pn[index-1]-0.5):
#     return xx[index-1]
#   elif np.abs(pn[index]-0.5) == np.abs(pn[index-1]-0.5):
#     return np.mean(xx[index-1:index+1])

def wpercentile(x,q,weights=None,axis=None):
  """
  Calculate weighted percentile
  """
  if weights is None:
    return np.percentile(x,q,axis=axis)
  
  if axis is not None:
    raise NotImplementedError("I thought this would be easy to implement. Turns out it is a pain. Will look at it if it is ever needed. The issue comes with axis in np.argsort")
  
  ix = np.argsort(x,axis=axis)
  if axis is None:
    xx, ww = x.ravel()[ix], weights.ravel()[ix]
  else:
    pass
  nx = len(xx)
  Sn = np.cumsum(ww)
  pn = 1./Sn[-1] * (Sn - ww/2) # weighted percentiles
#   print pn
  index = np.searchsorted(pn,q)
#   print index
  if not hasattr(index,'__iter__'):
    if index == 0:
      return xx[index]
    elif index == nx:
      return xx[index-1]
    elif np.isclose(np.abs(pn[index]-q),np.abs(pn[index-1]-q)):
      return np.mean(xx[index-1:index+1])
    elif np.abs(pn[index]-q) < np.abs(pn[index-1]-q):
      return xx[index]
    elif np.abs(pn[index]-q) > np.abs(pn[index-1]-q):
      return xx[index-1]
  else:
    res = np.empty(len(index))
    for i,ii in enumerate(index):
      if ii == 0:
        res[i] = xx[ii]
      elif ii == nx:
        res[i] = xx[ii-1]
      elif np.isclose(np.abs(pn[ii]-q[i]),np.abs(pn[ii-1]-q[i])):
        res[i] = np.mean(xx[ii-1:ii+1])
      elif np.abs(pn[ii]-q[i]) < np.abs(pn[ii-1]-q[i]):
        res[i] = xx[ii]
      elif np.abs(pn[ii]-q[i]) > np.abs(pn[ii-1]-q[i]):
        res[i] = xx[ii-1]
    return res