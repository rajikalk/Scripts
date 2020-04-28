#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import interp1d
import matplotlib.pyplot as mp
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
def zucker(wlf,f,wlg,g,res=None,ran=100.,fittype='gauss',nchp=1,ofac=2.,returnML=False,normalize=True,plot=False,title='CrossCorr.pdf'):
  """s_hat, s_hat_sigma, dv, C, [Ne, a_hat, sigma_hat] = ZUCKER(wlf,f,wlg,g,res=None,ran=100,fittype='gauss',nchp=1,returnML=False,plot=False,title='CrossCorr.pdf')

ZUCKER returns the cross-correlation function (CCF) between a spectrum, f, and a template spectrum, g.

REFERENCES
Zucker, Mon. Not. R. Astron. Soc. 342, 1291-1298 (2003)

INPUT
     wlf, f : wavelength and values of the original spectrum
     wlg, g : wavelength and values of the template spectrum

KEYWORDS
        res : resolution of the CCF in km/s (if not given a value is calculated)
        ran : range of the Cross Correlatoin in km/s. The CCF go from -ran to ran in steps of res
       nchp : number of pieces in which to chop the spectrum. Default: 1 (the spectrum is not chopped)
       ofac : over sampling parameter. Default: 2
    fittype : fit a function to the CCF to estimate midpoint. 'gauss' or 'polynomial'.
   returnML : if spectrum is chopped into pieces return a Maximum Likelihooed estimate (True/False)
  normalize : normalize the error estimate (Errors will be reduced to to linear interpolation of spectrum). Default: True
       plot : make a plot (True/False)
      title : title of the plot (Default: 'CrossCorr.pdf')

OUTPUT
         dv : vector of velocity shifts
          C : CCF
      s_hat : midpoint of the CCF
s_hat_sigma : uncertainty of the midpoint
      a_hat : estimated amplitude factor of the template
  sigma_hat : estiamted error in the original spectrum

DESCRIPTION
ZUCKER calculates the CCF between a spectrum and a template using the Maximum Likelihood approach described by Zucker (2003). The function returns the CCF itself, as well as a an estimate of the velocity shift and error.

The spectrum is resampled on a new wavelength grid corresponding to constant velocity steps, with a resolution either determined by input or calculated by the program.

The spectrum can be chopped up into a number of equally sized pieces set by nchp. This produces nchp CCF's which can either be examined individually or combined using the Maximum Likelihood approach from Zucker (2003). The advantage of chopping the spectrum up is that the S/N ratio along the spectrum may be variable, and that regions with many narrow lines give a better determination of the shifts that few wide lines. By chopping the spectrum and combining using the Maximum likelihood approach the best parts of the spectrum is weighed higher. The Maximum Likelihood approach is used if nchp > 1 and returnML is True.

To determine the velocity shift a gauss curve or parabola. This is set by the keyword fittype, which needs to be either 'gauss' or 'polynomial'. Default is 'gauss'.

The spectrum is put onto a new wavelength grid using linear interpolation. This will artificially reduce the scatter of the spectrum, which will under estimate the errors. If Normalize keyword is set a normalization factor is calculated and applied to the error esimates

VERSION 0.1
LAST UPDATED 6/10 2012

(C) Søren Frimann
"""

  #Convert inputs to numpy arrays
  wlf = np.asarray(wlf)
  f   = np.asarray(f)
  wlg = np.asarray(wlg)
  g   = np.asarray(g)
  
  #check arguments
  if wlf.ndim != 1:
    raise TypeError, "expected 1D vector for wlf"
  if wlg.ndim != 1:
    raise TypeError, "expected 1D vector for wlg"
  if wlf.size == 0:
    raise TypeError, "expected non-empty vector for wlf"
  if wlg.size == 0:
    raise TypeError, "expected non-empty vector for wlg"
  if wlf.shape != f.shape:
    raise TypeError, "wlf and f must have same dimension"
  if wlg.shape != g.shape:
    raise TypeError, "wlg and g must have same dimension"
  if ran < 0:
    raise ValueError, "ran must be non-negative"
  if nchp != int(nchp) or nchp < 1:
    raise TypeError, "nchp must be a non-negative integer"
  if ofac < 0:
    raise ValueError, "ofac must be non-negative"
  if not fittype in ['gauss', 'polynomial']:
    raise ValueError, "fittype must be either 'gauss' or 'polynomial'"
  if returnML is True and nchp == 1:
    raise TypeError, "returning Maximum Likelihood estimate not not sensible when not combining spectra"
  
  #Create interpolation objects
  f_int = interp1d(wlf,f,kind='linear')
  g_int = interp1d(wlg,g,kind='linear')

  #white noise on same wavelength grid as spectrum. Used to estimate normalization factor
  noise = np.random.normal(size=len(wlf),scale=1.)
  noise_int = interp1d(wlf,noise)
  
  #Determine best resolution
  dlambda = np.median(np.diff(wlf))
  c       = 299792.458  #speed of light in km/s
  lambda0 = wlf[0]
  lambdan = wlf[-1]
  nres    = dlambda/lambda0 * c  #natural resolution in km/s
  if res is None: res = nres / ofac

  #put spectra on new grid uniform in velocity
  n     = np.ceil((np.log(lambdan) - np.log(lambda0))/np.log(1+res/c))  #number of points in new grid
  bound = int(np.floor(ran/res))  #Number of shifts (in bins) to make in the cross correlation function
  dn    = int(np.floor(n/nchp)) #number of bins in each piece
  n_beg = np.arange(nchp)*dn #vector of starting bins for each piece

  #initialize arrays
  C           = np.empty((nchp,2*bound+1))
  Ne          = np.empty(nchp)
  s_hat       = Ne.copy()
  s_hat_sigma = Ne.copy()
  a_hat       = Ne.copy()
  sigma_hat   = Ne.copy()
  norm_factor = Ne.copy()

  #-#-#-#-#-# iterate over pieces #-#-#-#-#-#
  if plot is True and not returnML:
    pdf = PdfPages(title)
  
  for i,nb in enumerate(n_beg):
    wl_low, wl_high = (1 + res/c)**(nb-bound)*lambda0, (1 + res/c)**(nb+dn+bound)*lambda0  #high and low wavelengths

    #if first wavelength in template > wl_low and last wavelength is < wl_high
    #the template will be 'bound' bins smaller than the original spectrum in both ends
    if wlg[0] > wl_low or wlg[-1] < wl_high:
      nn = np.arange(nb,nb+dn)
      wwlf = (1 + res/c)**nn*lambda0
      wwlg = wwlf[bound:-bound]
    else:  #else the template will be 'bound' bins larger than original spectrum in both ends
      nn = np.arange(nb-bound,nb+dn+bound)
      wwlg = (1 + res/c)**nn*lambda0
      wwlf = wwlg[bound:-bound]

    dv, C[i,:], Ne[i], s_hat[i], s_hat_sigma[i], a_hat[i], sigma_hat[i], C_fun = zucker_core(f_int(wwlf),g_int(wwlg),res,nres,fittype=fittype)

    #Normalization factor for error estimate
    norm_factor[i] = np.std(noise)/np.std(noise_int(wwlf))
    
    #plotting
    if plot is True and not returnML:
      mp.subplot(311)
      mp.plot(wwlf,f_int(wwlf))
      mp.title('Original Spectrum')
      mp.xlabel('Wavelength (AA)')
      mp.ylabel('Relative Flux')

      mp.subplot(312)
      mp.plot(wwlg,g_int(wwlg))
      mp.title('Template Spectrum')
      mp.xlabel('Wavelength (AA)')
      mp.ylabel('Relative Flux')
      
      mp.subplot(313)
      mp.plot(dv,C[i,:],',')
      mp.plot(dv, C_fun(dv),'r')
      mp.figtext(0.65,0.3,'dv = %.3f +/- %.3f' % (s_hat[i], s_hat_sigma[i]))
      mp.title('Cross Correlation Function')
      mp.xlabel('Velocity (km/s)')
      mp.ylabel('C(dv)')
      pdf.savefig()
      mp.close()

  if plot is True and not returnML:
    pdf.close()

  #Normalize error estimate
  if normalize is True:
    s_hat_sigma *= norm_factor
    sigma_hat   *= norm_factor

  if not returnML:
    return dv, C, s_hat, s_hat_sigma, a_hat, sigma_hat
  
  if nchp > 1:
    ML = np.sqrt(1 - (1 - C**2).prod(axis=0)**(1/nchp))

    #Roughly determine velocity shift
    mid    = dv[np.argmax(ML)]
    binmin = 10 #minimum number of velocity bins used for fitting
    sigma  = binmin*res #calculate sigma from binmin
    
    #fit parabola or gauss to better determine velocity shift
    if fittype is 'polynomial':
      weights, ML_fun, s_hat, MLs, d2MLs = fitpol(dv, ML, mid, sigma)
    elif fittype is 'gauss':
      weights, ML_fun, s_hat, MLs, d2MLs = fitgauss(dv, ML, mid, sigma)

    s_hat_sigma = np.sqrt(-1/Ne.sum()*MLs/d2MLs*(1-MLs**2)/MLs**2)

    #Normalize error estimate
    if normalize is True:
      norm_factor  = (norm_factor/sigma_hat**2).sum()/(1/sigma_hat**2).sum() #Normalize norm factor
      s_hat_sigma *= norm_factor
    
    if plot is True:
      ddv = np.linspace(dv[0], dv[-1], num=500)
      mp.plot(dv,ML,',')
      mp.plot(ddv,ML_fun(ddv),'r')
      mp.title('Maximum Likelihood Cross Correlation Function')
      mp.xlabel('Velocity (km/s)')
      mp.ylabel('ML(dv)')
      mp.figtext(0.65,0.85,'dv = %.3f +/- %.3f' % (s_hat, s_hat_sigma))
      mp.savefig(title)
      mp.close()

  return dv, ML, s_hat, s_hat_sigma, a_hat, sigma_hat

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
def zucker_core(f,g,res,nres,fittype='gauss'):
  #f is the spectrum
  #g is the template
  
  sf, sg = np.std(f), np.std(g) #standard deviation of spectra
  mf, mg = np.mean(f), np.mean(g) #mean of spectra
  N  = np.min((len(f),len(g))) #number of overlapping points
  Ne = N*res/nres #effective number of overlapping points
  
  #Cross correlation function
  C = np.correlate(f-mf,g-mg,mode='valid')/(N*sf*sg)
  
  #Create velocity axis
  bound = (len(C)-1)/2
  dv    = np.arange(-bound,bound+1)*res
  
  #Roughly determine velocity shift
  mid    = dv[np.argmax(C)]
  binmin = 10 #minimum number of velocity bins used for fitting
  sigma  = binmin*res #calculate sigma from binmin

  #fit parabola or gauss to better determine velocity shift
  if fittype is 'polynomial':
    weights, C_fun, s_hat, Cs, d2Cs = fitpol(dv, C, mid, sigma)
  elif fittype is 'gauss':
    weights, C_fun, s_hat, Cs, d2Cs = fitgauss(dv, C, mid, sigma)
  
  #determine main parameters
  sigma_hat = np.sqrt(sf**2*(1-Cs**2))
  a_hat = Cs*sf/sg
  s_hat_sigma = np.sqrt(-sigma_hat**2/(a_hat*Ne*d2Cs*sf*sg))
  sigma_hat_sigma = np.sqrt(sigma_hat**2/(2*N))
  a_hat_sigma = np.sqrt(sigma_hat**2/(N*sg**2))

  return dv, C, Ne, s_hat, s_hat_sigma, a_hat, sigma_hat, C_fun


#-#-#-#-#-# helper functions #-#-#-#-#-#
def gauss(x,amp,mid,sigma,offset):  #Gaussian helper function
  return amp*np.exp(-(x-mid)**2/(2*sigma**2)) + offset

def fitpol(dv, C, mid, sigma): #fit parabola
  weights = gauss(dv,1.,mid,sigma,0) #weigths for the fit
  C_fun   = np.poly1d(polyfit(dv,C,2,w=weights)[::-1]) #fitted polynomial

  s_hat = C_fun.deriv(1).r #shift (s_hat) that maximize Cross correlation function (root of the parabola)
  Cs = C_fun(s_hat)        #C(s_hat)
  d2Cs = C_fun.deriv(2)(s_hat)  #d2C(s_hat)/ds_hat2

  return weights, C_fun, s_hat, Cs, d2Cs

def fitgauss(dv, C, mid, sigma): #fit gauss
  weights = gauss(dv,0.99,mid,sigma,0.01)  #weights of the fit
  ww = 1./np.sqrt(weights)  #convert weights to sigma
  p  = curve_fit(gauss,dv,C,sigma=ww,p0=(1.,mid,10.,0.2))[0] #coefficients of the gaussfit
  C_fun = lambda x: gauss(x,p[0],p[1],p[2],p[3]) #fitted gaussian

  s_hat = p[1]  #shift (s_hat) that maximize Cross correlation function (midpoint of gauss)
  Cs    = p[0] + p[3]   #C(s_hat)
  d2Cs  = -p[0]/p[2]**2 #d2C(s_hat)/ds_hat2

  return weights, C_fun, s_hat, Cs, d2Cs