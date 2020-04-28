#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import numpy as np

def circMap(x,min=0,max=1):
  """
  Map points into a cyclic domain bounded by min and max
  """
  assert np.all(max > min), "maximum must be larger than minimum"
  size = max - min
  return (x - min) % size + min

def circmean(x,weights=None,axis=None,min=0,max=1,arithmetic=False):
  """
  Calculate circular mean
  Calling Sequence:
    m = circmean(x,min=0,max=1,weights=None,axis=None,arithmetic=False)
  Input:
    x  : Array containing data to be averaged. If x is not an array, a conversion is attempted
  Keywords:
    min        : minimum value in the range where x is cyclic
    max        : maximum value in the range where x is ciclic
    weights    : An array of weights associated with the values in x. Each value in a contributes 
                 to the average according to its associated weight. The weights array can either 
                 be 1-D (in which case its length must be the size of x along the given axis) or 
                 of the same shape as x. If weights=None, then all data in a are assumed to have 
                 a weight equal to one.
    axis       : Axis along which to average x. If None, averaging is done over the flattened array
    arithmetic : Do an arithmetic correction
  Notes:
    The function first maps the values of x onto the cyclic domain defined by min and max
    using circMap
    
    The circular mean is calculated by first mapping x onto a unit circle and finding the 
    average angle before mapping back to physical space
    .. math::
       
       x_\mathrm{map} = \frac{x}{\mathrm{max}-\mathrm{min}} * 2 * \pi - \pi
       
       \bar{x_\mathrm{map}} = \mathrm{arctan2} \left(\frac(\sum w_i \sin x_i}{\sum w_i},\frac(\sum w_i \cos x_i}{\sum w_i}\right)
       
       \bar{x} = \frac{\bar{x_\mathrm{map}} + \pi}{2\pi} \times (\mathrm{max} - \mathrm{min}) + \mathrm{min}
    
    The circular mean is usually close to the arithmetic mean, but not exactly the same.
    If the arithmetic keyword is True the values will be centered around the circular mean
    and the relative arithemtic mean will be calculated, assuming no circularity. Since the
    referene is close to the true mean, the assumed non-circularity should not introduce any
    errors and the arithmetic value is returned
  """
  
  # map points to cyclic space in case some fall outside
  x = circMap(x,min=min,max=max)
  
  #map to +/- pi
  x_map = (x-min)/(max-min) * 2*np.pi - np.pi
  # circular mean
  xmean_map = np.arctan2(np.average(np.sin(x_map),weights=weights,axis=axis),np.average(np.cos(x_map),weights=weights,axis=axis))
  #map back
  xmean = (xmean_map + np.pi)/(2*np.pi) * (max-min) + min
  
  if arithmetic:
    x = circMap(x,min=xmean-0.5*(max-min),max=xmean+0.5*(max-min))
    xmean = np.average(x,weights=weights,axis=axis)
  
  return xmean

def relativePosition(x1,x2,cyclic=True,period=1.):
  """
  Relative positions between two points
  
  Input:
    x1, x2 : The two points for which to calculate the relative position. Should be of
             identical sizes. The position is calculated relative to x1
  Keywords:
    cyclic : Use cyclic position cordinates (Default: True)
    period : period of cyclic coordinates. Only relevant if cyclic is True (Default: 1)
  """

  # ------- Sanity checks
  if np.ndim(x1) != np.ndim(x2):
    raise ValueError("Dimensionality of x1 and x2 should be the same")

  if np.ndim(x1) > 0:
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    if x1.shape != x2.shape:
      raise ValueError("x1 and x2 must have the same shape")

  # ------- Relative Euclidian distances
  if cyclic:
    return (((x1 - x2 + period/2.) % period) - period/2.)
  else:
    return x1 - x2

def dist(x1,x2,**kw):
  """
  Euclidian distance between two points
  
  Input:
    x1, x2 : The two points for which to calculate the distance. Should be zero- or one-
             or two-dimensional and of identical shapes. If the array is two-dimensional 
             different points are assumed to be stacked along the zeroth dimension
  Keywords:
    cyclic : Use cyclic position cordinates (Default: True)
    period : period of cyclic coordinates. Only relevant if cyclic is True (Default: 1)
  """
  relPos = relativePosition(x1,x2,**kw)
  if np.ndim(relPos) <= 1:
    return np.sqrt(np.sum(relPos**2))
  else:
    return np.sqrt(np.sum(relPos**2,axis=-1))
