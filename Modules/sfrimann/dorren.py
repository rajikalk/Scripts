#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
try:
  from numexpr import evaluate as neeval
  inumexpr = True
except ImportError:
  inumexpr = False
from math import pi
from time import time as systemtime
import matplotlib.pyplot as mp

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

class dorren(object):
  """spotobject = dorren(xdata=np.arange(0,1,0.001),ydata=None,sigma=1,nspot=1,
                         numexpr=False,longref=0,brightspot=False,
                         returnflux=False,constantval=[],constantpar=[])

dorren gives a spot class that can be used for calculating spot models using
Dorren (1987), as well as calculating likelihood given measured data

REFERENCES
Dorren, Astrophysical Journal, 320, 756-767 (1987)
Generalization to Quadratic Limb Darkening calculated by H.E. Frölich

KEYWORDS
  xdata       : Array of x-values. These x-values are also used when calculating models
  ydata       : Array of observations.
  sigma       : one-sigma errors. Default = 1. (No error)
  nspot       : Number of spots to be modeled (Default: 1)
  numexpr     : Use numexpr? (Default: False)
  longref     : Reference longitude. (Default: 0)
  brightspots : Model bright spots. (Default: False)
  returnflux  : Return flux values when calculating models. (Default: magnitude)
  contantval  : Array of parameter values to be kept constant
  constantpar : Array of paramter names to be kept constant (Couples to constantval)

OUTPUT
  spot object

DESCRIPTION
Creates spot object that can be used for calculating models, and used in fit-
ting routines. Calling½  the object returns ln(likelihood)

VERSION 0.2
LAST UPDATED 1/4  2012
             28/8 2012

(C) Søren Frimann
"""
  # ------ Initializing -------------------------------------------------------
  def __init__(self,**kw):
    self.xdata       = kw.get('xdata',np.arange(0,1,0.001))  # xdata (create if not present)
    self.ydata       = kw.get('ydata',None)        # ydata
    self.sigma       = kw.get('sigma',1.)          # measurement errors
    self.nspot       = kw.get('nspot',1)           # number of spots to model
    self.numexpr     = kw.get('numexpr',False)     # use numexpr?
    self.longref     = kw.get('longref',0.)        # reference longitude
    self.brightspot  = kw.get('brightspot',False)  # model bright spots?
    self.returnflux  = kw.get('returnflux',False)  # If true return flux. (Default: magnitude)
    self.vsini       = kw.get('vsini',1.)          # stellar vsini
    self.constantval = kw.get('constantval',[])
    constantpar      = kw.get('constantpar',[])
    nspot            = self.nspot

    # ------ Handle observational data ----------------------------------------
    self.xdata       = np.asanyarray(self.xdata)
    if self.ydata != None:
      self.ydata     = np.asanyarray(self.ydata)
    if np.size(self.sigma) == 1:
      self.sigma     = np.repeat(self.sigma,len(self.xdata))
    else:
      self.sigma     = np.asanyarray(self.sigma)

    # ------ check arguments --------------------------------------------------
    if self.ydata is not None:
      if self.ydata.shape != self.xdata.shape:
        raise TypeError, "xdata and ydata must have same dimension"
    if nspot != int(nspot) or nspot < 1:
      raise ValueError, "nspot must be a non-negative integer"
    if len(constantpar) != len(np.unique(constantpar)):
      raise ValueError, "don't repeat elements in constantpar"
    if len(self.constantval) != len(constantpar):
      raise TypeError, "constantval and constantpar must be same length"

    # ------ std of y data ----------------------------------------------------
    if self.ydata != None:
      self.ystd = np.std(self.ydata)
    else:
      self.ystd = 0.
    
    # ------ define paramter ranges -------------------------------------------
    slatmin , slatmax , slatname    = -1.  , 1.    , 'sin_lat'
    longmin , longmax , longname    = -180., 180.  , 'long'
    sizemin , sizemax , sizename    = -25. , 0.    , 'size'
    permin  , permax  , pername     = 0.   , 1.    , 'period'
    cinclmin, cinclmax, cinclname   = 0.   , 1.    , 'cos_incl'
    uamin   , uamax   , uaname      = 0.   , 1.    , 'ua'
    ubmin   , ubmax   , ubname      = 0.   , 1.    , 'ub'
    smin    , smax    , sname       = 0.   , 4.    , 's'
    if self.brightspot == True:
      Fspotmin, Fspotmax, Fspotname = 1.   , 2.    , 'Fspot'
    else:
      Fspotmin, Fspotmax, Fspotname = 0.   , 1.    , 'Fspot'

    # ------ parnames ---------------------------------------------------------
    num   = np.arange(1,nspot+1).astype('S2')
    pname = np.hstack((np.core.defchararray.add(np.repeat(slatname,nspot),num),
                       np.core.defchararray.add(np.repeat(longname,nspot),num),
                       np.core.defchararray.add(np.repeat(sizename,nspot),num),pername,cinclname,uaname,ubname,Fspotname,sname))
    pmin  = np.hstack((np.repeat(slatmin,nspot),np.repeat(longmin,nspot),np.repeat(sizemin,nspot),permin,cinclmin,uamin,ubmin,Fspotmin,smin))
    pmax  = np.hstack((np.repeat(slatmax,nspot),np.repeat(longmax,nspot),np.repeat(sizemax,nspot),permax,cinclmax,uamax,ubmax,Fspotmax,smax))

    # ------ Error ------------------------------------------------------------
    if np.in1d(constantpar,pname).all() == False:
      raise ValueError, "one or more elements in constantpar not in pname"

    # ------ Set constants ---------------------------------------------------
    self.wh    = np.squeeze(np.array([np.where(item == pname) for item in constantpar]))
    self.mask  = np.in1d(pname,constantpar)
    self.pname = np.ma.array(pname,mask=self.mask)
    self.pmin  = np.ma.array(pmin,mask=self.mask)
    self.pmax  = np.ma.array(pmax,mask=self.mask)

    self.npar  = self.pmin.count() #number of paramters
    self.ndim  = self.pmin.size
    self.N     = self.xdata.size
    
  # ------ priors -------------------------------------------------------------
  def flatprior(self,par,pmin,pmax,invcum=False):
    """Calculates flat prior probability given par, pmin, pmax
       invcum keyword is used in draw_from_prior"""
    if invcum is False:
      if (par > pmax) or (par < pmin):
        return 0
      else:
        return 1/(pmax-pmin)
    else:
      return par*(pmax-pmin)+pmin
        
  # ------ Transform parameters -----------------------------------------------
  def fullpar(self,par_in):
    """Combine constant and variable paramters to get all paramters to be
    passed onto the model function"""
    if len(self.constantval) == 0:
      par_out = np.asanyarray(par_in)
    elif len(par_in) == 0:
      par_out = np.asanyarray(self.constantval)
    else:
      par_out = np.empty(self.ndim)
      par_out[self.wh], par_out[~self.mask] = self.constantval, par_in

    if len(par_out) != 3*self.nspot+6:
      raise ValueError, 'len(par_out) != number of paramters for spotmodel'

    self.slat   = par_out[0:self.nspot]
    self.long0  = par_out[self.nspot:2*self.nspot]
    self.size   = par_out[2*self.nspot:3*self.nspot]
    self.period = par_out[3*self.nspot]
    self.cincl  = par_out[3*self.nspot+1]
    self.ua     = par_out[3*self.nspot+2]
    self.ub     = par_out[3*self.nspot+3]
    self.Fspot  = par_out[3*self.nspot+4]
    self.s      = par_out[3*self.nspot+5]
    
    self.long0 = (self.long0 + 180.)%360. - 180.
    
    return par_out

  # ------ Draw samples from prior distribution -------------------------------
  def draw_from_prior(self):
    """Draw samples from the prior distribution """
    par = np.random.rand(self.npar)
    for i,(pmi,pma) in enumerate(zip(self.pmin.compressed(),self.pmax.compressed())):
      par[i] = self.flatprior(par[i],pmi,pma,invcum=True)
    return par

  # ------ Model function -----------------------------------------------------
  def model(self,*arg,**kw):
    """Wrapper to actual model function
    Input is paramter array
    If called with fulloutput keyword it also returns light curves for individual
    spots"""
    
    fulloutput = kw.get('fulloutput',False) # also output result from individual spots
    returnrv   = kw.get('returnrv',False)
    if len(arg) == 0:
      par_in = []
    elif len(arg) == 1:
      par_in = arg[0]
    else:
      raise ValueError, "model takes maximum one input"
    par    = self.fullpar(par_in)

    # ------ Calculate models -------------------------------------------------
    if returnrv:
      temp = np.empty((self.nspot,self.N))
      for i in range(self.nspot):
        temp[i,:] = spotmodel(self.xdata,self.slat[i],self.long0[i],self.size[i],self.period,self.cincl,ua=self.ua,ub=self.ub,Fspot=self.Fspot,numexpr=self.numexpr,longref=self.longref,returnflux=True)
      return rv(self.xdata,temp,self.slat,self.long0,self.period,self.vsini,longref=self.longref)
    if not fulloutput:
      return spotmodel(self.xdata,self.slat,self.long0,self.size,self.period,self.cincl,ua=self.ua,ub=self.ub,Fspot=self.Fspot,numexpr=self.numexpr,longref=self.longref,returnflux=self.returnflux)
    else:
      result = np.empty((self.nspot+1,self.N))
      result[0,:] = spotmodel(self.xdata,self.slat,self.long0,self.size,self.period,self.cincl,ua=self.ua,ub=self.ub,Fspot=self.Fspot,numexpr=self.numexpr,longref=self.longref,returnflux=self.returnflux)
      for i in range(self.nspot):
        result[i+1,:] = spotmodel(self.xdata,self.slat[i],self.long0[i],self.size[i],self.period,self.cincl,ua=self.ua,ub=self.ub,Fspot=self.Fspot,numexpr=self.numexpr,longref=self.longref,returnflux=self.returnflux)
      return result

  # ------ Logprior -----------------------------------------------------------
  def logprior(self,par):
    """Calculate ln(prior prob) """
    if len(par) == 0: return 0
    lp = 0
    for p,(pn,(pmi,pma)) in zip(par,zip(self.pname.compressed(),zip(self.pmin.compressed(),self.pmax.compressed()))):
      if 'long' in pn:
        p = (p+pma)%(pma-pmi) + pmi
      lp += np.log(self.flatprior(p,pmi,pma))
    return lp

  # ------ Log likelihood -----------------------------------------------------
  def loglike(self,par_in):
    """Calculates ln(likelihood) assuming IID Gaussian uncertainties
    Ignores any offset by integrating it out analytically"""
    ymodel = self.model(par_in)
    a      = (1/(2*(self.sigma**2 + (self.ystd*self.s)**2))).sum() #constants for integration
    b      = ((ymodel - self.ydata)/(self.sigma**2 + (self.ystd*self.s)**2)).sum()
    return (- self.N/2*np.log(2*pi) - 0.5*np.log(self.sigma**2 + (self.ystd*self.s)**2).sum() - ((self.ydata - ymodel)**2./(2*(self.sigma**2 + (self.ystd*self.s)**2))).sum() #classical likelihood
           + 0.5*np.log(pi/a) + b**2/(4*a) - 1) #analytical integration of constant offset

  # ------ convert longitudes -------------------------------------------------
  def convert_longitude(self,longitudes):
    """ Converts a list of longitudes from radian to degree.
    Also transforms cycially according to the formula long = (long+180.)%360. - 180."""

    return (np.asanyarray(longitudes)*360/2/pi+180.)%360. - 180.
    

  # ------ Return posterior ---------------------------------------------------
  def __call__(self,*arg):
    """Calculates posterior probability"""
    if self.logprior(*arg) == -np.inf:
      return -np.inf
    else:
      return self.logprior(*arg) + self.loglike(arg)

  # ------ Speed tester for different models ----------------------------------
  def test(self):
    """test speed of using numexpr """
    par = [self.draw_from_prior() for i in xrange(100)]
    numexpr_temp = self.numexpr
    for self.numexpr in [True, False]:
      t0 = systemtime()
      for i in par:
        dump = self.__call__(i)
      print "numexpr: %s. %f seconds" % (self.numexpr, (systemtime() - t0)/100)
    self.numexpr, self.simple = numexpr_temp, simple_temp

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
def rv(x,result,sin_lat,long0_d,period,vsini,longref=0.):

  # ------ Initializing -------------------------------------------------------
  ns    = np.asarray(sin_lat).size  # Number of spots
  nd    = np.asarray(x).size        # Number of data points

  rv = np.zeros(nd)
  for i in range(ns):
    F   = result[i,:] - 1
    phi = 2*pi*(x - longref)/period + long0_d[i]*2*pi/360
    rv += -F*vsini*np.cos(np.arcsin(sin_lat[i]))*np.sin(phi)

  return rv

def spotmodel(x,sin_lat,long0_d,size,period,cos_incl,ua=0.,ub=0.,Fspot=0.,
              numexpr=True,longref=0,returnflux=False):
  """ymodel = spotmodel(x,sin_lat,long0_d,size,period,cos_incl,ua=0.,ub=0.,
                        Fspot=0.,numexpr=True,longref=0,returnflux=False)

spotmodel calculates the synthetic light curve for a number of circular spots
with given parameters.

REFERENCES
Dorren, Astrophysical Journal, 320, 756-767 (1987)
Generalization to Quadratic Limb Darkening calculated by H.E. Frölich

INPUT
  x          : xdata
  sin_lat    : sin(spot latitude). Input either scalar or array of length nspot
  long0_d    : Spot latitude in degree (reference longitude given by longref)
  size       : Spot size = ln(sin^2(alpha))
  period     : Period of the spot(s) in same units as the xdata. If a scalar all
               spots have the same period. If an 1D array with same number of en-
               tries as spots the spots will have different periods
  cos_incl   : cos(stellar inclination)
KEYWORDS
  ua         : Linear limb-darkening coefficient of the star, and spot
               (Default: 0, no limb-darkening)
  ub         : Quadratic limb-darkening coefficient of the star, and spot
               (Default: 0, no limb-darkening)
  Fspot      : Flux from the spots relative to the star (Default: 0, dark spot)
  numexpr    : Use the numexpr module for calculations which is expected to be
               faster for large arrays, but might not always be available
               (Default: False)
  longref    : Reference longitude (Default = 0)
  returnflux : If True return flux (Default: magnitude)

OUTPUT
     ymodel : 1D array with the model, same length as x. In magnitude.

DESCRIPTION
Calculates a spot model as described by Dorren (1987). For the time being the
model can only handle constant spot sizes. Generalized to quadratic limb dark-
ening. Limb darkening of spots is set equal to limb darkening of star.

VERSION 0.2
LAST UPDATED 19/4 2012
             24/8 2012 (Added Quadratic Limb Darkening)

(C) Søren Frimann
"""
  # ------ Check if numexpr is available --------------------------------------
  if numexpr is True and inumexpr is False:
    numexpr = False
  # ------ Initializing -------------------------------------------------------
  ns    = np.asarray(sin_lat).size  # Number of spots
  nd    = np.asarray(x).size        # Number of data points
  Fstar = 1.  # Stellar flux

  # ------ Converting from trigonometric functions to angles ------------------
  lat   = np.arcsin(np.asanyarray([sin_lat]))
  incl  = np.arccos(cos_incl)
  # ------ Converting from degree to radians-----------------------------------
  long0 = np.asanyarray([long0_d])*2*pi/360

  # ------ Creating arrays with correct dimensions ----------------------------
  xx    = np.asanyarray([(x - longref)]).repeat(ns,axis=0)
  lat   = lat.repeat(nd,axis=0).T
  long0 = long0.repeat(nd,axis=0).T
  size  = np.asanyarray([size],dtype=np.float64).repeat(nd,axis=0).T
  if np.asarray(period).size == ns:
    period = np.asanyarray([period],dtype=np.float64).repeat(nd,axis=0).T

  # ------ Stellar limb darkening = Spot limb darkening -----------------------
  ua_star, ub_star, ua_spot, ub_spot = ua, ub, ua, ub

  # ------ Calculate models ---------------------------------------------------
  # ------ Determining a, b, (and c) ------------------------------------------
  a = (1. - ua_star -   ub_star) - (1 - ua_spot -   ub_spot)*Fspot/Fstar
  b = (     ua_star + 2*ub_star) - (    ua_spot + 2*ub_spot)*Fspot/Fstar
  c = (             -   ub_star) - (            -   ub_spot)*Fspot/Fstar
  
  if numexpr:
    az = neeval('xx*2*pi/period + long0')

    # ------ Spot angular size ------------------------------------------------
    alpha = neeval('arcsin(sqrt(exp(size)))')  # Spot radius

    # ------ Determining angles -----------------------------------------------
    beta  = neeval('arccos(cos(incl)*sin(lat) + sin(incl)*cos(lat)*cos(az))')
    delta = np.nan_to_num(neeval('where(beta - alpha >= pi/2,pi,arccos(1./(tan(alpha)*tan(beta))))'))
    zeta  = neeval('arcsin(sin(delta)*sin(alpha))')

    # ------ Determining T ----------------------------------------------------
    T = neeval('where(beta <= pi/2,arctan(sin(zeta)*tan(beta)),pi - arctan(-sin(zeta)*tan(beta)))')

    # ------ Trigonometric arrays ---------------------------------------------
    sbeta , cbeta  = neeval('sin(beta)') , neeval('cos(beta)')
    salpha, calpha = neeval('sin(alpha)'), neeval('cos(alpha)')
    
    # ------ Determining A, B, (and C) ----------------------------------------
    if ub_star != 0:  # Only calculate C if ub != 0
      A01  = neeval('sbeta*sin(delta)')
      B011 = neeval('(2*((5*sbeta**2-2)*salpha**2+4*cbeta**2)*(delta-pi)*salpha**2+(sbeta+2)*(sbeta-2)*delta)*cbeta+2*(A01**2-3)*A01*alpha')
      B012 = neeval('-(3*(2*(A01**2-3)*salpha**2-3*(2*salpha**4-1)*cbeta**2)+4*(A01**2-12*sbeta**2+9)*salpha**4)')
      B01  = neeval('3*salpha*B011+A01*calpha*B012')
      B021 = neeval('-((A01**2-3)*A01*alpha-(3*cbeta*delta-2*zeta))')
      B022 = neeval('(sbeta**2+2)*delta*salpha-A01*calpha*cbeta')
      B02  = neeval('6*salpha*cbeta**2*B021-3*cbeta**3*B022')

    Aint = neeval('zeta + (pi - delta)*cbeta*salpha**2 - sin(zeta)*sbeta*calpha')
    Bint = neeval('1./3.*(pi - delta)*(-2*calpha**3 - 3*sbeta**2*calpha*salpha**2) + 2./3.*(pi - T) + 1./6.*sin(zeta)*sin(2*beta)*(2-3*calpha**2)')
    
    if ub_star != 0:
      Cint = neeval('-(cbeta**2*B01+B02)/(24*cbeta**2*salpha)')
    else:
      Cint = 0
  else:  # numexpr not available or not wanted
    np.seterr(invalid='ignore')
    az = xx*2*pi/period + long0

    # ------ Spot angular size ------------------------------------------------
    alpha = np.arcsin(np.sqrt(np.exp(size))) #Spot radius
    
    # ------ Determining angles -----------------------------------------------
    beta  = np.arccos(np.cos(incl)*np.sin(lat) + np.sin(incl)*np.cos(lat)*np.cos(az))
    delta = np.nan_to_num(np.where(beta - alpha >= pi/2,pi,np.arccos(1./(np.tan(alpha)*np.tan(beta)))))
    zeta  = np.arcsin(np.sin(delta)*np.sin(alpha))

    # ------ Determining T ----------------------------------------------------
    T = np.where(beta <= pi/2,np.arctan(np.sin(zeta)*np.tan(beta)),pi - np.arctan(-np.sin(zeta)*np.tan(beta)))

    # ------ Trigonometric arrays ---------------------------------------------
    sbeta , cbeta  = np.sin(beta) , np.cos(beta)
    salpha, calpha = np.sin(alpha), np.cos(alpha)

    # ------ Determining A, B, (and C) ----------------------------------------
    if ub_star != 0:  # Only calculate C if ub != 0
      A01  = sbeta*np.sin(delta)
      B011 = (2*((5*sbeta**2-2)*salpha**2+4*cbeta**2)*(delta-pi)*salpha**2+(sbeta+2)*(sbeta-2)*delta)*cbeta+2*(A01**2-3)*A01*alpha
      B012 = -(3*(2*(A01**2-3)*salpha**2-3*(2*salpha**4-1)*cbeta**2)+4*(A01**2-12*sbeta**2+9)*salpha**4)
      B01  = 3*salpha*B011+A01*calpha*B012
      B021 = -((A01**2-3)*A01*alpha-(3*cbeta*delta-2*zeta))
      B022 = (sbeta**2+2)*delta*salpha-A01*calpha*cbeta
      B02  = 6*salpha*cbeta**2*B021-3*cbeta**3*B022

    Aint = zeta + (pi - delta)*cbeta*salpha**2 - np.sin(zeta)*sbeta*calpha
    Bint = 1./3.*(pi - delta)*(-2*calpha**3 - 3*sbeta**2*calpha*salpha**2) + 2./3.*(pi - T) + 1./6.*np.sin(zeta)*np.sin(2*beta)*(2-3*calpha**2)
    
    if ub_star != 0:
      Cint = -(cbeta**2*B01+B02)/(24*cbeta**2*salpha)
    else:
      Cint = np.array([0,0])

      np.seterr(invalid='print')

  # ------ Writing and returning lightcurve -----------------------------------
  if returnflux is True:
    return 1 - (a*Aint.sum(axis=0) + b*Bint.sum(axis=0) + c*Cint.sum(axis=0) )/(pi*(1 - ua_star/3 - ub_star/6))
  else:
    return -2.5 * np.log10( 1 - (a*Aint.sum(axis=0) + b*Bint.sum(axis=0) + c*Cint.sum(axis=0) )/(pi*(1 - ua_star/3 - ub_star/6)) )
