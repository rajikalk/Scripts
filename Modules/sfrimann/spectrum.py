#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import array, float64, mean, nan_to_num, empty, median, diff, max as nmax, min as nmin, arange, where, logical_and, sin, cos
from time import time as systemtime
import matplotlib.pyplot as mp
from math import pi
from numexpr import evaluate as neeval

# ---- Error definitions ----
class InputError(Exception):
  pass

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def spectrum_core(time, data, freq, weights=None):
  """Calling Sequence:
    alpha, beta = spectrum_core(time, data, freq, Weights=None)

  Input:
    time    : Time array
    data    : Data array
    freq    : Frequency array with frequencies to be evaluated

  Output:
    alpha   : alpha coefficient (see below)
    beta    : beta coefficient (see below)

  Keywords:
    Weights : Array with weights to be used. If this array is not present no weighting will be used

  Description:
    Computes the Lomb-Periodogram as defined in Lomb 1975, with added weights as described in Zechmeiter and Kürster 2008.
    Returns the linear coefficients of the fitted expression:

      y_fit = alpha*sin(2*pi*freq*time) + beta*cos(2*pi*freq*time)

    spectrum_core is an internal function and is not supposed to be called directly
  """

  #Number of frequencies
  n = freq.size

  #Create matrix of frequency multiplied by time (time varying along 0-axis i.e. vertically)
  arg  = array([time],dtype=float64).repeat(n, axis=0).T * array(freq,dtype=float64)
  #Create matrix of data
  y    = array([data],dtype=float64).repeat(n, axis=0).T

  #Calculate sines and cosines
  s, c = neeval('sin(2*pi*arg)'), neeval('cos(2*pi*arg)')

  #Calculate weighed weighed elements depending on whether weights are defined or not
  if weights is None:
    wy     = y
    ws, wc = s, c
  else:
    ww     = (array([weights]).repeat(n, axis=0)).T #weights
    wy     = neeval('ww * y')
    ws, wc = neeval('ww * s'), neeval('ww * c')

  #Calculating  weighted s, c, ss, cc, sc for all frequencies
  ys = neeval('wy * s').sum(axis = 0)
  yc = neeval('wy * c').sum(axis = 0)
  ss = neeval('ws * s').sum(axis = 0)
  cc = neeval('wc * c').sum(axis = 0)
  sc = neeval('ws * c').sum(axis = 0)

  #alpha and beta
  D     = neeval('ss*cc-sc*sc')
  alpha = nan_to_num(neeval('(ys*cc - yc*sc)/D'))
  beta  = nan_to_num(neeval('(yc*ss - ys*sc)/D'))
  
  return alpha, beta

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def chunkeval(time,data,freq,maxentries=10000000,weights=None):
  """Calling Sequence:
    alpha, beta, power = chunkeval(time, data, freq, maxentires=10000000, weights=None)

  Input:
    time       : Time array
    data       : Data array
    freq       : Frequency array with frequencies to be evaluated

  Output:
    alpha      : alpha coefficient (see spectrum_core)
    beta       : beta coefficient (see spectrum_core)
    power      : Power Spectrum

  Keywords:
    maxentries : Maximum number of entries in the combined time/frequency matrix (for memory management)
    weights    : Array with weights to be used. If this array is not present no weighting will be used

  Description:
    Super-routine for spectrum_core to help with evaluation of the Lomb-Scargle periodogram.
    The routine chops up the freqency array to keep the memory usage down.
    The routine also calculates the power spectrum using the alpha and beta coefficients.
  """

  if freq.size*time.size < maxentries: #time and frequency arrays are small enough, no need to chop up.
    result = array(spectrum_core(time,data,freq,weights=weights))
  else: #Chop up using maxentries
    result = empty([2,len(freq)],dtype=float64)  #initialize result array
    nfreq  = maxentries/time.size  #number of frequencies per iteration
    for i in range(0, freq.size, nfreq):  #iterate over frequencies
      #Call spectrum_core and fill in result array
      result[0, i:i+nfreq], result[1, i:i+nfreq] = spectrum_core(time,data,freq[i:i+nfreq],weights=weights)
  return result[0], result[1], (result * result).sum(axis = 0)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def powerspectrum(*args,**kw):
  """Calling Sequence:
    alpha, beta, power, [freq] = powerspectrum(time, data, [freq,] weights=None, timeit=True, ofac=4)

  Input:
    time       : Time array
    data       : Data array
    freq       : Frequency array with frequencies to be evaluated. If not present one will be generated

  Output:
    alpha      : alpha coefficient (see spectrum_core)
    beta       : beta coefficient (see spectrum_core)
    power      : Power Spectrum
    freq       : Generated frequency array (if not present in input)

  Keywords:
    weights    : Array with weights to be used. If this array is not present no weighting will be used
    timeit     : If true prints timing information
    ofac       : Oversampling parameter

  Description:
    Main Routine for evaluation of power-spectra
  """

  # ------ Starting time
  t0 = systemtime()
  
  # ------- Handle keywords
  weights = kw.get('weights',None)
  timeit  = kw.get('timeit',True)
  ofac    = kw.get('ofac',4.)

  # ------ Handle arguments
  if len(args) == 2:
    time = args[0]
    data = args[1]
  elif len(args) == 3:
    time = args[0]
    data = args[1]
    freq = args[2]
  else:
    raise InputError('Wrong number of inputs')

  time = array(time,dtype=float64).squeeze()
  data = array(data,dtype=float64).squeeze()

  # ------ Subtract zero frequency
  ddata = data - mean(data)

  # ------ Handle frequency array
  if len(args) == 2:
    dt    = median(diff(time))
    nyq   = 1./(2*dt)
    tdiff = nmax(time) - nmin(time)
    freq  = arange(0,nyq,1./(ofac*tdiff),dtype=float64)
  else:
    freq = array(freq,dtype=float64).squeeze()

  # ------ Handle weights
  if weights is not None:
    weights = array(weights,dtype=float64).squeeze()

  # ------ Calculate Power Spectrum
  alpha, beta, power = chunkeval(time,ddata,freq,weights=weights)

  # ------ Ending time
  if timeit:
    print 'spectrum.py: powerspectrum finished in %3.1f seconds' % (systemtime() - t0)
  
  # ------ Return
  if len(args) == 2:
    return alpha, beta, power, freq
  else:
    return alpha, beta, power

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def windowfunction(*args,**kw):
  """Calling Sequence:
    power, [freq] = windowfunction(time, [freq,] weights=None, timeit=True, ofac=4, wcf=None)

  Input:
    time       : Time array
    freq       : Frequency array with frequencies to be evaluated. If not present one will be generated

  Output:
    power      : Power Spectrum
    freq       : Generated frequency array (if not present in input)

  Keywords:
    weights    : Array with weights to be used. If this array is not present no weighting will be used
    timeit     : If true prints timing information
    ofac       : Oversampling parameter
    wcf        : If given this will the be frequency for the pseudo window. If not given a central frequency is used

  Description:
    The window function for a Lomb periodogram is not strictly defined. This functions calculates a pseudo-window.
  """
  # ------ Starting time
  t0 = systemtime()

  # ------- Handle keywords
  weights = kw.get('weights',None)
  timeit  = kw.get('timeit',True)
  wcf     = kw.get('wcf',None) #Window center frequency
  ofac    = kw.get('ofac',4.)
  
  # ------ Handle arguments
  if len(args) == 1:
    time = args[0]
  elif len(args) == 2:
    time = args[0]
    freq = args[1]
  else:
    raise InputError('Wrong number of inputs')

  time = array(time,dtype=float64).squeeze()

  # ------ Handle frequency array
  if len(args) == 1:
    dt    = median(diff(time))
    nyq   = 1./(2*dt)
    tdiff = nmax(time) - nmin(time)
    freq  = arange(0.01*nyq,0.99*nyq,1./(ofac*tdiff),dtype=float64)
  else:
    freq = array(freq,dtype=float64).squeeze()
  
  # ------ Handle weights
  if weights is not None:
    weights = array(weights,dtype=float64).squeeze()

  # ------ Estimate Window Center Frequency
  if wcf is None:
    wcf = median(freq)

  # ------ Make data
  arg  = 2*pi*wcf*time
  #data = 0.5 * (np.sin(arg) + np.cos(arg))
  
  # ------ Calculate Power Spectrum (Only power output)
  power = 0.5*(chunkeval(time,sin(arg),freq,weights=weights)[-1] + chunkeval(time,cos(arg),freq,weights=weights)[-1])

  # ------ Ending time
  if timeit:
    print 'spectrum.py: windowfunction finished in %3.1f seconds' % (systemtime() - t0)

  # ------ Return
  if len(args) == 1:
    return power, freq
  else:
    return power

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  """Calling Sequence:
    bpdata = bandpass(time, data, flow, fhigh, [freq,] weights=None, invert=False, ofac=4, timeit=True)

  Input:
    time       : Time array
    data       : Data array
    flow       : Lower freqyency boundary
    fhigh      : Upper frequency boundary
    freq       : frequency array (optional)

  Output:
    bpdata     : bandpassed data

  Keywords:
    weights    : Array with weights to be used. If this array is not present no weighting will be used
    timeit     : If true prints timing information
    ofac       : Oversampling parameter
    invert     : Invert range

  Description:
    The window function for a Lomb periodogram is not strictly defined. This functions calculates a pseudo-window.
  """
def bandpass(*args,**kw):
  # ------ Starting time
  t0 = systemtime()

  # ------- Handle keywords
  weights = kw.get('weights',None)
  invert  = kw.get('invert',False)
  timeit  = kw.get('timeit',True)
  ofac    = kw.get('ofac',4.)

  # ------ Handle arguments
  if len(args) == 4:
    time  = args[0]
    data  = args[1]
    flow  = args[2]
    fhigh = args[3]
  elif len(args) == 5:
    time  = args[0]
    data  = args[1]
    flow  = args[2]
    fhigh = args[3]
    freq  = args[4]
  else:
    raise InputError('Wrong number of inputs')

  time = array(time,dtype=float64).squeeze()
  data = array(data,dtype=float64).squeeze()
  
  # ------ Handle weights
  if weights is not None:
    weights = array(weights,dtype=float64).squeeze()
  
  # ------ Calculate Power Spectrum
  if len(args) == 4:
    alpha, beta, power, freq = powerspectrum(time,data,weights=weights,timeit=False,ofac=ofac)
  if len(args) == 5:
    freq = array(freq,dtype=float64).squeeze()
    alpha, beta, power = powerspectrum(time,data,freq,weights=weights,timeit=False,ofac=ofac)

  # ------ Calculate Window function
  dt     = median(diff(time))
  nyq    = 1./(2*dt)
  dfreq  = median(diff(freq))
  wfreq  = arange(0.05*nyq,0.95*nyq,dfreq,dtype=float64)
  winsum = windowfunction(time,wfreq,weights=weights,timeit=True).sum()
  
  # ------ Bandpass range
  bins_inside  = where(logical_and(freq >= flow, freq <= fhigh),1,0).sum()
  bins_outside = freq.size - bins_inside
  if (bins_inside <= bins_outside and not invert) or (bins_inside >= bins_outside and invert):
    normal = True
  else:
    normal = False
  print 'normal:', normal
  print 'invert:', invert
  if (invert and normal) or (not invert and not normal):
    print 'forkert'
    index = where(~logical_and(freq >= flow, freq <= fhigh))
  elif (~invert and normal) or (invert and ~normal):
    print 'rigtig'
    index = where(logical_and(freq >= flow, freq <= fhigh))
  
  n        = time.size
  alphamat = array([alpha[index]],dtype=float64).repeat(n,axis=0).T
  betamat  = array([beta[index]],dtype=float64).repeat(n,axis=0).T
  #Create matrix of frequency multiplied by time (frequency varying along 0-axis i.e. vertically)
  arg      = array([freq[index]],dtype=float64).repeat(n, axis=0).T * array(time,dtype=float64)
  #Calculate sines and cosines
  s, c = neeval('sin(2*pi*arg)'), neeval('cos(2*pi*arg)')
  
  if normal:
    bpdata = neeval('alphamat * s + betamat * c').sum(axis=0)/winsum + mean(data)
  else:
    bpdata = data - neeval('alphamat * s + betamat * c').sum(axis=0)/winsum
  
  # ------ Ending time
  if timeit:
    print 'spectrum.py: bandpass finished in %3.1f seconds' % (systemtime() - t0)

  # ------ Return
  return bpdata

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

if __name__ == '__main__':

  time = np.arange(0,30,0.1)

  amp = [1.3, 5.1, 2.1, 1.0]
  f   = [1.0, 2.2, 0.7, 3.3]
  ph  = [0.5, 3.1, 2.2, 1.1]
  m   = 0

  data = np.zeros(time.size) + m
  for i in range(len(amp)):
    data += amp[i]*np.sin(2*pi*time*f[i]+ph[i])
  
#  data = 1.2*np.sin(2*pi*time*1.3)
  freq = np.arange(0.05,4.5,0.001)
  weights = None# np.random.random(len(time))
  #print weights

  alpha, beta, power= powerspectrum(time,data,freq,ofac=10)
  bpdata = bandpass(time,data,1.9,2.5,freq,ofac=10)
  
  alpha, beta, power2 = powerspectrum(time,bpdata,freq,ofac=10)

  print np.max(np.sqrt(power)), np.max(np.sqrt(power2))
  
  mp.figure(1)
  mp.plot(time,data)

  mp.figure(2)
  mp.plot(freq,np.sqrt(power))
  mp.plot(freq,np.sqrt(power2),'r')
  mp.ylim((0,6))
  
  mp.figure(3)
  mp.plot(time,bpdata)

  mp.figure(4)
  mp.plot(freq,np.sqrt(power2))
  mp.ylim((0,6))
  
  
  mp.show()