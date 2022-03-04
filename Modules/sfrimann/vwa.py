#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os.path import abspath, exists
import numpy as np
import matplotlib.pyplot as mp
from timeseries import smooth, hampel

# ---- Error definitions ----
class NumberError(Exception):
  pass

class FileError(Exception):
  pass

class ProgramError(Exception):
  pass


def import_observations(filename):
  """Calling Sequence:
    wl, fl = import_observations(filename)

  Input:
    filename : path to the obs-file containing the spectrum

  Output:
    wl       : array containing wavelength
    fl       : array containing flux measurements
    
  Keywords:
    None

  Description:
    Reads spectra in the obs format used in vwa (developed for idl)
  """

  #Path to file
  filepath = abspath(filename)

  #Check if file exists
  if not exists(filepath):
    raise FileError("Can't find obs-file")

  #Read file
  lines = open(filepath).read().split('\n')

  #Number of points
  n0    = int(lines[1])

  #Create wavelength and flux lists
  wl = lines[2:n0+2]
  fl = lines[n0+2:-1]

  #Convert to numpy arrays
  wl = np.array(wl,dtype=np.float64)
  fl = np.array(fl,dtype=np.float64)

  return wl, fl

def ptp_scatter(data):
  """Calling Sequence:
    noise = ptp_scater(data)

  Input:
    data  : Data array for the point to point scatter

  Output:
    noise : Estimate noise

  Keywords:
    None

  Description:
    Robust estimater of point to point scatter in a spectrum
  """
  nmin = 20 #Minimum number of points in the spectrum to get estimation

  #Median absolute deviation of data (Normalized to standard deviation)
  mad = 1.4826*np.median(np.abs(data - np.median(data)))

  #Index of data more than thee sigma away
  index = np.where(np.abs(data - np.median(data)) < 3.*mad)

  #Number of points available for estimation
  n     = np.where(np.abs(data - np.median(data)) < 3.*mad,1,0).sum()

  #Check number of legimate points
  if n < nmin:
    raise NumberError('n < nmin')

  #Noise estimate
  noise = np.std(data[index]) #np.sqrt((np.diff(data[index])**2).sum()/((n-1)))

  return noise

def estimate_sn(wl,fl,dx=None,threshold=0.95):

  index = np.where(fl > threshold)
  nindex = np.where(fl <= threshold)

  if dx is None:
    noise = ptp_scatter(fl[index])
    return 1./noise
  else:
    wli, fli = wl[index], fl[index]
    #mfl = smooth(fli,window='flat')
    mfl, s0, gp, n, no = hampel(wli,fli,dx=2.,plot=False,threshold=0)
    wli, fli = wl[index], fl[index]

#  index = ~np.isnan(noise)
#  wli, fli, noise, = wli[index], fli[index], noise[index]

#  sn = 1./noise

#  p    = np.poly1d(np.polyfit(wli,sn,2)) #Fit polynomial
#  poly = p(wl)                       #Tabulated polynomial

#  mp.plot(wli,sn,',')
#  mp.plot(wl,poly,'r')
#  mp.show()
  return wli, fli, mfl, s0

  
def estimate_sn2(wl,fl,dx=None,threshold=0.95):
  
  index = np.where(fl > threshold)
#  index = np.argsort(fl)
#  index = index[-int((1-threshold)*len(index)):]
  print len(fl), len(index[0])
  
  if dx is None:
    noise = ptp_scatter(fl[index])
    return 1./noise
  else:
    wli, fli = wl[index], fl[index]
    noise = np.empty(len(fli))
    for i in xrange(len(fli)):
      noise[i] = window(wli,fli,dx,i)

  index = ~np.isnan(noise)
  wli, fli, noise, = wli[index], fli[index], noise[index]
  
  sn = 1./noise
  
  p    = np.poly1d(np.polyfit(wli,sn,2)) #Fit polynomial
  poly = p(wl)                       #Tabulated polynomial

  mp.plot(wli,sn,',')
  mp.plot(wl,poly,'r')
  mp.show()
  return wli, fli, sn, poly
  
  
def window(x,y,dx,i,threshold=5.):

  dn    = np.median(np.diff(x))
  ddx   = np.diff(x)
  index = np.where(ddx > threshold*dn)
  xborder = x[index]
  ddx   = ddx[index]
  if len(xborder) == 0:
    lx = x[i] - dx
    ux = x[i] + dx
#    print 'case0'
  elif len(xborder[np.where(xborder < x[i])]) == 0:
    lx = x[i] - dx
    ux = xborder[0]
#    print 'case1'
  elif len(xborder[np.where(xborder >= x[i])]) == 0:
    lx = xborder[-1]
    ux = x[i] + dx
#    print 'case2'
  elif len(xborder[np.where(xborder < x[i])]) and len(xborder[np.where(xborder >= x[i])]) > 0:
    lx = xborder[np.where(xborder < x[i])][-1] + ddx[np.where(xborder < x[i])][-1]
    ux = xborder[np.where(xborder >= x[i])][0]
    #print 'case3'
  else:
    raise ProgramError('Something went wrong')

  if lx < (x[i] - dx): lx = x[i] - dx
  if ux > (x[i] + dx): ux = x[i] + dx

  index = np.where(np.logical_and(lx <= x, x <= ux)) #Window for this iteration

  try:
    noise = ptp_scatter(y[index])
  except NumberError:
    noise = np.nan
  return noise