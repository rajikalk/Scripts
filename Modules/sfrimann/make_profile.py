#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from weights import wmedian

def make_profile(X,Y,**kw):
  """
  Make profile of data
  """
  
  # ------- Handle Keywords
  nbin  = kw.get('nbin',100)
  nbinx = kw.get('nbinx',nbin)
  nbiny = kw.get('nbiny',nbin)
  xlog  = kw.get('xlog',True)
  ylog  = kw.get('ylog',True)
  weights = kw.get('weights',np.ones(len(X)))
  sigma = kw.get('sigma',1.)
  
  # ------- Sort by distance from centre
  index   = np.argsort(X)
  xx      = X[index]
  yy      = Y[index]
  ww      = weights[index]
  
  if ylog is True:
    index = np.where(yy > 0)
    xx    = xx[index]
    yy    = yy[index]
    ww    = ww[index]
  
  # ------- log or not
  xx     = np.log10(xx) if xlog is True else xx
  yy     = np.log10(yy) if ylog is True else yy

  # ------- Bin the data
  bin_edges = np.linspace(xx.min(),xx.max(),num=nbinx+1)
  bin_index = np.searchsorted(bin_edges,xx,side='right')
  bin_spacing = np.median(np.diff(bin_edges))
  
  # ------- Initialise Arrays
  x  = np.interp(np.arange(nbinx)+0.5,np.arange(nbinx+1),bin_edges)
  
  # ------- Calculate median in each bin
  y, s0 = np.empty(len(x)), np.empty(len(x))
  for i,bin_center in enumerate(x):
    cnt, const = 0, 1.
    while cnt <= 5:
      index, = np.where((xx >= (bin_center-const*bin_spacing/2.)) & (xx <= (bin_center+const*bin_spacing/2.)))
      cnt = len(index)
      const += 1
    y[i]  = wmedian(yy[index],ww[index])
    s0[i] = sigma*1.4826*wmedian(np.abs(yy[index] - y[i]), ww[index])
  
  # ------- If logarithmic transform back to linear coordinates
  x = 10**x if xlog is True else x
  if ylog is True:
    upper_error = 10**(y + s0)
    lower_error = 10**(y - s0)
    y = 10**y
    s0 = 10**s0
  else:
    upper_error = y + s0
    lower_error = y - s0
  
  # ------- Calculate 2d Histogram
  histxy, yedges, xedges = np.histogram2d(yy,xx,bins=[nbiny,nbinx],weights=weights)
  histx = np.interp(np.arange(nbinx)+0.5,np.arange(nbinx+1),xedges)
  histy = np.interp(np.arange(nbiny)+0.5,np.arange(nbiny+1),yedges)
  histx = 10**histx if xlog is True else histx
  histy = 10**histy if ylog is True else histy

  # ------- normalize ysum to 1
  sumy = histxy.sum(axis=0)
  for i,s in enumerate(sumy):
    if s > 0.:
      histxy[:,i] /= s

  return {'x':x,'y':y,'s0':s0,'lower_error':lower_error,'upper_error':upper_error,\
          'nbin':nbin,'xlog':xlog,'ylog':ylog,'histxy':histxy,'histx':histx,\
          'histy':histy}
