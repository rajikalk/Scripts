#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division

from numpy import poly1d, polyfit, argsort, median, diff, where, logical_and, abs, empty
import matplotlib.pyplot as mp
import numpy as np
from time import time

# ---- Error definitions ----
class InputError(Exception):
  pass

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
def polydiv(x,y,order,plot=True,xlabel='x',ylabel='y',filename='polydata.pdf'):
  """Calling Sequence:
    yc = polydiv(x,y,order,plot=True,xlabel='x',ylabel='y',filename='polydata.pdf')

  Input:
    x, y     : Data arrays
    order    : Order of the fitted polynomial

  Output:
    yc       : Corrected data

  Keywords:
    plot     : Make a plot (Default: True)
    xlabel   : xlabel (Default: 'x')
    ylabel   : ylabel (Default: 'y')
    filename : filename (Default: 'polydata.pdf')

  Description:
    Takes the inputdata, fits a polynomial and divides.
  """
  p    = poly1d(polyfit(x,y,order)) #Fit polynomial
  poly = p(x)                       #Tabulated polynomial
  yc   = y/poly                     #Corrected data
  
  if plot is True:
    mp.subplot(211)
    mp.plot(x,y)
    mp.plot(x,poly,'r')
    mp.title('Original Data')
    mp.ylabel(ylabel)
    
    mp.subplot(212)
    mp.plot(x,yc)
    mp.title('Corrected Data')
    mp.xlabel(xlabel)
    mp.ylabel(ylabel)
    
    mp.savefig(filename)
    mp.close()
  
  return yc

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def datafold(x,y,period,x0=0.,plot=True,xlabel='Phase',ylabel='y',filename='folddata.pdf',phase=True):
  """Calling Sequence:
    xf, yf = datafold(x,y,period,t0=0.,plot=True,xlabel='x',ylabel='y',filename='polydata.pdf')

  Input:
    x, y     : Data arrays
    period   : Period to fold with

  Output:
    xf,yf   : Folded data

  Keywords:
    x0       : x-value for zero phase
    plot     : Make a plot (Default: True)
    xlabel   : xlabel (Default: 'x')
    ylabel   : ylabel (Default: 'y')
    filename : filename (Default: 'polydata.pdf')

  Description:
    Take the input data and fold with period.
  """
  xf = (x - x0) % period #Folded x-data
  si = argsort(xf)       #Sorting index
  xf = xf[si]            #Sorting x-data
  yf = y[si]             #Folded y-data
  
  if phase is True: xf /= period
  
  if plot is True:
    mp.plot(xf,yf,',')
    mp.title('Folded Data')
    mp.xlabel(xlabel)
    mp.ylabel(ylabel)
    
    mp.savefig(filename)
    mp.close()
  
  return xf, yf

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
def hampel(x,y,threshold=3.,dx=None,full=False,weights=None,timeit=False,percentile=0.):
  """Calling Sequence:
    yy, s0, gp, n, no = hampel(x,y,threshold=3.,dx=None.,plot=True,xlabel='x',ylabel='y',filename='filterdata.pdf')

  Input:
    x, y     : Data arrays (Row or column vectors -- y should be normally distributed)

  Output:
    yy       : Filtered data
    s0       : Array of estimated uncertainies
    gp       : Array with indices of good points (non-filtered values)
    n        : Arrray of points included for every filter iteration
    no       : Number of outliners

  Keywords:
    threshold : threshold value (Default: 3 -- 0 corresponds to a median filter)
    dx        : Window half-size of the window (Default: 3*median(diff(x)) )

  Description:
    hampel returns the Hampel filtered values of the elements in Y. It was developed
    to detect outliers in a time series, but it can also be used as an alternative to
    the standard median filter.

    References:
    Chapters 1.4.2, 3.2.2 and 4.3.4 in Mining Imperfect Data: Dealing with
    Contamination and Incomplete Records by Ronald K. Pearson.

    http://exploringdatablog.blogspot.com/2012/01/moving-window-filters-and-pracma.html

    Adapted from MATLAB code originally written by Michael L. Nielsen
    http://www.mathworks.com/matlabcentral/fileexchange/34795
  """
  
  t0 = time()
  
  #Default value for dx
  if dx is None:
    dx = 3*median(diff(x))

  #Initialize filtered data as copy of original data
  yy = y.copy()

  #Initializing
  y0 = empty(len(y))
  s0 = y0.copy()
  n  = y0.copy()
  
  xs    = np.argsort(x)
  ilow  = np.searchsorted(x[xs],x-dx)
  ihigh = np.searchsorted(x[xs],x+dx,side='right')
  
  #Main loop
  for i in range(len(y)):
    #index = where(logical_and(x[i] - dx <= x, x <= x[i] + dx)) #Window for this iteration
    #n[i]  = index[0].size #Number of points inside window
    index = xs[ilow[i]:ihigh[i]]
    n[i] = index.size
    if percentile > 0:
      y0[i] = np.percentile(y[index][np.isfinite(y[index])],percentile) #Determine median
      s0[i] = 1.4826*median(abs(y[index][np.isfinite(y[index])] - y0[i])) #Median absolute deviation
    elif weights == None:
      y0[i] = median(y[index][np.isfinite(y[index])]) #Determine median
      s0[i] = 1.4826*median(abs(y[index][np.isfinite(y[index])] - y0[i])) #Median absolute deviation
    else:
      y0[i] = wmedian(y[index][np.isfinite(y[index])],weights[index][np.isfinite(y[index])]) #Determine median
      s0[i] = 1.4826*wmedian(abs(y[index][np.isfinite(y[index])] - y0[i]),weights[index][np.isfinite(y[index])]) #Median absolute deviation
  
  gp = where(abs(y - y0) <= threshold*s0) #Index of good points
  ol = where(abs(y - y0) > threshold*s0)  #Index of outliers
  
  yy[ol] = y0[ol]     #Replace outliners
  no     = ol[0].size #Number of outliners
  
  if timeit == True:
    
    hours   = (time()-t0)/3600.
    minutes = (hours % 1)*60.
    seconds = (minutes % 1)*60.
    
    print("hampel finished in %.0f h %.0f m %.0f s" % (hours, minutes, seconds))
  
  if full is True:
    return yy, s0, gp, n, no
  else:
    return yy, s0

# ------- Helper Function
def wmedian(x,w):
  """
  Calculate weighted median
  """
  
  if len(x) == 0:
    return np.nan
  
  if w == None:
    return median(x)
  
  assert len(w) == len(x), 'wmedian: data and weights must be same shape'
  
  ix = np.argsort(x)
  xx, ww = x[ix], w[ix]
  nx = len(x)
  Sn = np.cumsum(ww)
  pn = 1./Sn[-1] * (Sn - ww/2) # weighted percentiles
  index = np.searchsorted(pn,0.5)
  if index == 0:
    return xx[index]
  elif index == nx:
    return xx[index-1]
  elif np.abs(pn[index]-0.5) < np.abs(pn[index-1]-0.5):
    return xx[index]
  elif np.abs(pn[index]-0.5) > np.abs(pn[index-1]-0.5):
    return xx[index-1]
  elif np.abs(pn[index]-0.5) == np.abs(pn[index-1]-0.5):
    return np.mean(xx[index-1:index+1])

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=x#np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y
