#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np

import os

arcsec2rad = np.pi/648000.

def uvamp(u,v,visibilities,uvrange=None,offset=[0.0,0.0]):
  """
  uv amplitude function
  u and v are in wavelengths
  offset is in arcsec
  """
  
  # setup uv bins
  if uvrange is None:
    uvrange = np.arange(0,max(u.max()/1000.,v.max())/1000.,5)
  else:
    uvrange = np.asarray(uvrange)
  
  nbin = uvrange.size - 1
  
  uvamp, sigmean, sn, expect, npoint = np.empty(nbin), np.empty(nbin), np.empty(nbin), np.empty(nbin), np.empty(nbin)
  
  uu, vv = u/1000., v/1000.
  
  vis = visibilities * np.exp(-2*np.pi*1j*(offset[0]*uu+offset[1]*vv)*arcsec2rad*1000)
  
  for i in range(nbin):
    index     = (uu**2 + vv**2 >= uvrange[i]**2) & (uu**2 + vv**2 < uvrange[i+1]**2)
    
    npoint[i] = index.sum()
    if npoint[i] == 0:
      uvamp[i], sigmean[i], sn[i], expect[i] = 0, 0, 0, 0
      continue
    
    uvamp[i]   = np.absolute(np.mean(vis[index]))
    varr2      = ((vis[index].real**2).sum() - npoint[i]*np.mean(vis[index].real)**2)/(npoint[i]-1)
    vari2      = ((vis[index].imag**2).sum() - npoint[i]*np.mean(vis[index].imag)**2)/(npoint[i]-1)
    sigtot     = vis[index].real.mean()**2/uvamp[i]**2*varr2 + vis[index].imag.mean()**2/uvamp[i]**2*vari2
    sigmean[i] = np.sqrt(sigtot/(npoint[i]-2))
    sn[i]      = uvamp[i]/sigmean[i]
    expect[i]  = np.sqrt(np.pi/2)*sigmean[i]
  
  uvpoints = 0.5*(uvrange[1:]+uvrange[:-1])
  
  return {'uvrange':uvpoints,'uvamp':uvamp,'sigma':sigmean,'sn':sn,'expect':expect,'npoint':npoint}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def uvpoint(u,v=None,amp=1.,x0=0.,y0=0.,grid=False):
  """
  point source in uv-plane
  u and v are in wavelenghts
  """
  
  v = u if v is None else v
  
  x0    *= arcsec2rad
  y0    *= arcsec2rad
  
  if grid:
    uu, vv = np.meshgrid(u,v,sparse=True)
  else:
    uu, vv = u, v
  
  return amp * np.exp(2*np.pi*1j*(x0*uu+y0*vv))

def uvgauss(u,v=None,amp=1.,fwhmx=1.,fwhmy=None,x0=0.,y0=0.,theta=0.,offset=0.,fluxdens=True,grid=False):
  """
  two dimensional gaussian in uv-plane
  theta is in degree
  u and v are in wavelenghts
  """
  
  v = u if v is None else v
  
  theta *= np.pi/180
  fwhmx *= arcsec2rad
  x0    *= arcsec2rad
  y0    *= arcsec2rad
  
  if grid:
    uu, vv = np.meshgrid(u,v,sparse=True)
  else:
    uu, vv = u, v
  
  fwhmy = fwhmx if fwhmy is None else fwhmy*arcsec2rad
  
  sigx, sigy = fwhmx/(2*np.sqrt(2*np.log(2))), fwhmy/(2*np.sqrt(2*np.log(2)))
  Sigx, Sigy = 1./(2*np.pi*sigx), 1./(2*np.pi*sigy)
  
  a =   np.cos(theta)**2/(2*Sigx**2) + np.sin(theta)**2/(2*Sigy**2)
  b = - np.sin(2*theta)/(4*Sigx**2)  + np.sin(2*theta)/(4*Sigy**2)
  c =   np.sin(theta)**2/(2*Sigx**2) + np.cos(theta)**2/(2*Sigy**2)
  
  amp = amp / (2*np.pi*sigx*sigy) if fluxdens else amp
  
  return 2 * np.pi * amp * sigx * sigy * np.exp(-(a*uu**2 + 2*b*uu*vv + c*vv**2)) * np.exp(2*np.pi*1j*(x0*uu+y0*vv)) + offset
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class Uvfit(object):

  def __init__(self,u,v,visibility,modelfun,thetamin,thetamax,**kw):  
    uvmin         = kw.get('uvmin',15)
    uvmax         = kw.get('uvmax',200)
    u2            = kw.get('u2',None)
    v2            = kw.get('v2',None)
    vis2          = kw.get('vis2',None)
#     self.niter    = kw.get('niter',1000)
#     self.nwalkers = kw.get('nwalker',100)
#     self.threads  = kw.get('threads',1)
  
    self.modelfun = modelfun
    self.thetamin = thetamin
    self.thetamax = thetamax
  
    self.ndim       = len(thetamax)
    self.prior_type = np.empty(self.ndim,dtype='U4')
    self.prior_type.fill('flat')
  
    index = (np.sqrt((u/1000.)**2 + (v/1000.)**2) >= uvmin) & (np.sqrt((u/1000)**2 + (v/1000)**2) <= uvmax)
    
    self.uu  = u[index]
    self.vv  = v[index]
    self.vis = visibility[index]
    
    if vis2 is None:
      self.more_data = False
    else:
      self.more_data = True
      
      index = (np.sqrt((u2/1000.)**2 + (v2/1000.)**2) >= uvmin) & (np.sqrt((u2/1000)**2 + (v2/1000)**2) <= uvmax)
      
      self.uu2  = u2[index]
      self.vv2  = v2[index]
      self.vis2 = vis2[index]
  
  def lnprior(self,theta):
    """
    log prior.
    Calculates the prior probability given an array of parameters.
    Uses the type of priors and value ranges from __init__
    """
  
    result = 0.
    for t,tmin,tmax,ptype in zip(theta,self.thetamin,self.thetamax,self.prior_type):
      if (t > tmax) or (t < tmin):
        return -np.inf
      else:
        if ptype == 'flat':
          result += 1./(tmax-tmin)
        else:
          raise ValueError('Do not recotnize type of prior %s' % ptype)
  
    return result

#   def lnlike(self,theta):
#     """
#     calculate log likelihood
#     Give array of parameters as input
#     """
#     
#     lns = theta[-1]
#     
#     model1    = self.modelfun(self.uu,self.vv,theta[:-1])
#     residual1 = np.hstack((self.vis.real,self.vis.imag)) - np.hstack((model1.real,model1.imag))
#     sigma1    = 10**(lns)
#     lnlike1   = -0.5*np.sum(residual1**2/sigma1**2 + np.log(2*np.pi*sigma1**2))
#     
#     if self.more_data:
#       model2    = self.modelfun(self.uu2,self.vv2,theta[:-1])
#       residual2 = np.hstack((self.vis2.real,self.vis2.imag)) - np.hstack((model2.real,model2.imag))
#       sigma2    = 10**(lns)
#       lnlike2   = -0.5*np.sum(residual2**2/sigma2**2 + np.log(2*np.pi*sigma2**2))
#     else:
#       lnlike2   = 0
# 
#     return lnlike1 + lnlike2

  def lnlike(self,theta):
    """
    calculate log likelihood
    Give array of parameters as input
    """
    
    model1    = self.modelfun(self.uu,self.vv,theta)
    residual1 = np.hstack((self.vis.real,self.vis.imag)) - np.hstack((model1.real,model1.imag))
    lnlike1   = -residual1.size/2.*np.log(np.sum(residual1**2))
    
    if self.more_data:
      model2    = self.modelfun(self.uu2,self.vv2,theta)
      residual2 = np.hstack((self.vis2.real,self.vis2.imag)) - np.hstack((model2.real,model2.imag))
      lnlike2   = -residual2.size/2.*np.log(np.sum(residual2**2))
    else:
      lnlike2   = 0

    return lnlike1 + lnlike2

  def lnprob(self,theta):
    """
    posterior probability
    """
    lp = self.lnprior(theta)
    if not np.isfinite(lp):
      return -np.inf
    return lp + self.lnlike(theta)
#     return self.lnlike(theta)