#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

from scipy.stats import binned_statistic as bin_stat
import numpy as np

def binned_statistic(xdata,ydata,range=None,xlog=False,ylog=False,nbin=10,sigma=1.):
  if range is None:
    if xlog:
      range = [np.log10(xdata.min()),np.log10(xdata.max())]
    else:
      range = [xdata.min(),xdata.max()]
  else:
    if xlog:
      range = [np.log10(range[0]),np.log10(range[1])]

  xval = np.log10(xdata) if xlog else xdata
  yval = np.log10(ydata) if ylog else ydata
  
  def mad(x):
    return 1.4826*np.median(np.abs(np.median(x) - x))
  
  count, xe, binno = bin_stat(xval,yval,statistic='count',bins=nbin,range=range)
  y , xe, binno = bin_stat(xval,yval,statistic='median',bins=nbin,range=range)
  ys, xe, binno = bin_stat(xval,yval,statistic=mad,bins=nbin,range=range)
  xc = 0.5*(xe[1:] + xe[:-1])
  
  yl, yu = y-sigma*ys, y+sigma*ys
  
  y  = 10**y if ylog else y
  yl = 10**yl if ylog else yl
  yu = 10**yu if ylog else yu
  xc = 10**xc if xlog else xc
  
  return xc, y, yl, yu, ys, count
