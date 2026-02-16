#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division

import numpy as np
import itertools
import os

import astropy.units as unit

from scipy.integrate import odeint
from scipy.integrate import ode

class FreezeoutSublimationIterator(object):
  
  def __init__(self,mu,max_abund,min_abund=None,Tevap=None,Eb=None,kcrd=0.):
    self.mu    = mu
    self.Tevap = Tevap
    self.kcrd  = kcrd
    self.Eb    = Eb
    self.max_abund = max_abund
    
    self.min_abund = self.max_abund/100 if min_abund is None else min_abund
    
    if (self.Eb is None) & (self.Tevap is None):
      raise ValueError('Tevap or Eb has to be given')

  def kdep(self, Tgas, nH, mu=None):
    """
    Follows Rodgers & Charnley (2003, ApJ, 585, 355) with atomic weights
    given by NIST (http://www.nist.gov/pml/data/comp.cfm).
    See the discussion surrouding eq. (3) of Rodgers & Charnley for
    more information
    
    Input:
      Tgas : Gas temperature in Kelvin
      nH   : Hydrogen number density
    Keywords:
      mu   : molecular weight. If not given __init__ mu is used (Default: None)
    Returns:
      ksub : Sublimation rate in s^-1
    """
    if mu is None:
      return 4.55e-18 * np.sqrt(Tgas / self.mu) * nH
    else:
      return 4.55e-18 * np.sqrt(Tgas / mu) * nH

  def ksub(self, Tdust, nu_vibr=2.0e12,Eb=None,Tevap=None):
    """
    Calculate the thermal evaporation time scale following Rodgers &
    Charnley (2003, ApJ, 585, 355); the time scale is simply the
    reciprocal of the evaporation rate (eq. 4 of Rodgers & Charnley).
    
    Input:
      Tdust   : Dust temperature in Kelvin
    Keywords:
      nu_vibr : The default value for the vibrational frequency of the input molecule 
                (`nu_vibr`) also comes from Rodgers & Charnley.
      Eb      : Binding energy of the molecule in Kelvin (Eb/kB). If None value from
                __init__ function is used (Default: None)
    Returns:
      kdep    : Depletion (freeze-out) rate in s^-1
    """
    
    if Eb is not None:
      return nu_vibr * np.exp(-Eb/Tdust)
    
    if self.Eb is None:
      self.Eb = -self.Tevap*np.log(unit.s.to('yr')/nu_vibr)
    if self.Tevap is None:
      self.Tevap = -self.Eb/np.log(unit.s.to('yr')/nu_vibr)
    
    return nu_vibr * np.exp(-self.Eb/Tdust)
  
  def initialise(self,*args,equilibrium=True):
    """
    Initialise abundance values.
    Input:
      nH    : Hydrogen number density
      Tgas  : Gas temperature in Kelvin.
      Tdust : Dust temperature in Kelvin. If not given assumed equal to Tgas
    Keywords:
      equilibrium : If True initialise abundances according to the ordinary differential
                    equation. If False set abundance equal to max abundance when evapora-
                    tion time scale < 1 yr (Default: True)
    Returns:
      abundance   : abundance values
    """
    
    if len(args) == 2:
      nH    = args[0]
      Tgas  = args[1]
      Tdust = args[1]
    elif len(args) == 3:
      nH    = args[0]
      Tgas  = args[1]
      Tdust = args[2]
    else:
      raise ValueError("initialise takes two or three inputs. {} given".format(len(args)))
    
    if equilibrium:
      abundance = (self.kcrd+self.ksub(Tdust))/(self.kdep(Tgas,nH) + self.kcrd + self.ksub(Tdust))*self.max_abund
    else:
      abundance = np.where(self.ksub(Tdust) > unit.s.to('yr'), self.max_abund, self.min_abund)
#     abundance = np.where(self.ksub(Tdust) > unit.s.to('yr'), self.max_abund, self.kcrd/(self.kdep(Tgas,nH2) + self.kcrd)*self.max_abund)
    try:
      abundance[abundance < self.min_abund] = self.min_abund
    except:
      if abundance < self.min_abund:
        abundance = self.min_abund
    
    return abundance
  
  def analyticalStep(self,*args,dt=10.,equilibrium=True,nstep=1):
    """
    Do one step of integration 
    Input:
      abund       : Starting abundances
      nH          : Hydrogen number density
      Tgas        : Gas temperature in Kelvin
      Tdust       : Dust temperature in Kelvin. Assumed equal to Tgas if not given
    Keywords:
      dt          : Time step in yr (is broken down into smaller steps of 1 yr in the 
                    solver)
      equilibrium : 
    """
    
    if len(args) == 3:
      abund = args[0]
      nH    = args[1]
      Tgas  = args[2]
      Tdust = args[2]
    elif len(args) == 4:
      abund = args[0]
      nH    = args[1]
      Tgas  = args[2]
      Tdust = args[3]
    else:
      raise ValueError('step takes three or four inputs. {} given'.format(len(args)))
    
    if isinstance(nH,list):
      nH0 = nH[0]
      nH1 = nH[1]
    else:
      nH0 = nH1 = nH

    if isinstance(Tgas,list):
      Tgas0 = Tgas[0]
      Tgas1 = Tgas[1]
    else:
      Tgas0 = Tgas1 = Tgas

    if isinstance(Tdust,list):
      Tdust0 = Tdust[0]
      Tdust1 = Tdust[1]
    else:
      Tdust0 = Tdust1 = Tdust
    
    t  = np.linspace(0,1,num=nstep+1)
    
    # ------- depletion and sublimation rates
    ksub0 = self.ksub(Tdust0) * dt * unit.yr.to('s')
    ksub1 = self.ksub(Tdust1) * dt * unit.yr.to('s')
    kdep0 = self.kdep(Tgas0,nH0) * dt * unit.yr.to('s')
    kdep1 = self.kdep(Tgas1,nH1) * dt * unit.yr.to('s')
    
    p0 = ksub0 + kdep0 + self.kcrd
    p1 = ksub1 + kdep1 + self.kcrd
    q0 = ksub0 + self.kcrd
    q1 = ksub1 + self.kcrd

    def solution(x,init):
      from scipy.special import erf
      from numpy.lib.scimath import sqrt
      part1 = np.exp(-((p1-p0)*x**2/2 + x*p0))
      part2 = (q0 - q1)/(p0-p1)*np.exp((p1-p0)*x**2/2 + x*p0)
      part3 = np.sqrt(np.pi/2)*np.exp(p0**2/(2*p0-2*p1))*erf((p0*(x-1)-p1*x)/(np.sqrt(2)*sqrt(p0-p1)))/((p0-p1)+0j)**1.5
      
      print(q1)
      
      res  = part1*(self.max_abund*(part2 + part3*(p0*q1 - p1*q0)) + init )
        
      return res
    
#     def solution(time,init):
#       from scipy.special import erfi
#       from numpy.lib.scimath import sqrt
#       part1 = np.exp(-((p1-p0)*(time-t0)**2/(2*(t1-t0)) + (time-t0)*p0))
#       erfi1 = erfi((p0*(time-t1)+p1*(t0-time))/(np.sqrt(2)*sqrt(p0-p1)*sqrt(t0-t1)))
#       erfi2 = erfi((p0*sqrt(t0-t1))/(np.sqrt(2)*sqrt(p0-p1)))
#       part2 = (np.sqrt(2*np.pi) * np.complex(t0-t1)**1.5 * np.exp(p0**2*(t1-t0)/(2*p0-p1)) * (erfi1 - erfi2) +
#                self.max_abund*sqrt(p0-p1)*(time-t0)*(q0*(time+t0-2*t1)+q1*(t0-time)))/(2*sqrt(p0-p1)*(t0-t1))
#       
#       return part1*(part2 + init)

#     def solution(time,init):
#       def func(x,i):
# #         return (np.exp((p1-p0)[i]*(x-t0)**2/(2*(t1-t0)) + (x-t0)*p0[i]) + ((x-t0)*(q1-q0)[i]/(t1-t0)*self.max_abund))
#         return (np.exp((p1-p0)[i]*(x-t0)**2/(2*(t1-t0)) + (x-t0)*p0[i]) + ((x-t0)*(q1-q0)[i]/(t1-t0)*self.max_abund))
#       from scipy.special import erfi
#       from numpy.lib.scimath import sqrt
#       from scipy.integrate import quad
#       part1 = np.exp(-((p1-p0)*(time-t0)**2/(2*(t1-t0)) + (time-t0)*p0))
#       
#       part2 = np.array([quad(func,t0,time,args=(ii))[0] for ii in range(p0.size)])
#       print(part2)
#       
#       
#       return part1*(part2 + init)
    
    sol = [solution(tt,abund) for tt in t]
    
    sol = np.vstack(sol)
    
    return sol, t*dt
  
  def step(self,*args,dt=10.,equilibrium=True,nstep=1):
    """
    Do one step of integration 
    Input:
      abund       : Starting abundances
      nH          : Hydrogen number density
      Tgas        : Gas temperature in Kelvin
      Tdust       : Dust temperature in Kelvin. Assumed equal to Tgas if not given
    Keywords:
      dt          : Time step in yr (is broken down into smaller steps of 1 yr in the 
                    solver)
      equilibrium : 
    """
    
    if len(args) == 3:
      abund = args[0]
      nH    = args[1]
      Tgas  = args[2]
      Tdust = args[2]
    elif len(args) == 4:
      abund = args[0]
      nH    = args[1]
      Tgas  = args[2]
      Tdust = args[3]
    else:
      raise ValueError('step takes three or four inputs. {} given'.format(len(args)))
    
    if isinstance(nH,list):
      nH0 = nH[0]
      nH1 = nH[1]
    else:
      nH0 = nH1 = nH

    if isinstance(Tgas,list):
      Tgas0 = Tgas[0]
      Tgas1 = Tgas[1]
    else:
      Tgas0 = Tgas1 = Tgas

    if isinstance(Tdust,list):
      Tdust0 = Tdust[0]
      Tdust1 = Tdust[1]
    else:
      Tdust0 = Tdust1 = Tdust
    
    t  = np.linspace(0,dt,num=nstep+1)*unit.yr.to('s')
    t0 = t[0]
    t1 = t[-1]
    
    # ------- depletion and sublimation rates
    ksub0 = self.ksub(Tdust0)
    ksub1 = self.ksub(Tdust1)
    kdep0 = self.kdep(Tgas0,nH0)
    kdep1 = self.kdep(Tgas1,nH1)
    
#     solution = [abund]
#     for i in range(1,len(t)):
#       ksub = (t[i] - t0)*(ksub1 - ksub0)/(t1 - t0) + ksub0
#       kdep = (t[i] - t0)*(kdep1 - kdep0)/(t1 - t0) + kdep0
#       dt   = t[i] - t[i-1]
#       
#       nabund = (1 - kdep*dt)*solution[-1]
#       if equilibrium:
#         nabund = (self.max_abund - nabund)*(self.kcrd + ksub)*dt
#       else:
#         nabund = np.where(ksub > unit.s.to('yr'), self.max_abund, nabund)
# 
#       try:
#         nabund[nabund > self.max_abund] = self.max_abund
#         nabund[nabund < self.min_abund] = self.min_abund
#       except:
#         if nabund > self.max_abund:
#           nabund = self.max_abund
#         elif nabund < self.min_abund:
#           nabund = self.min_abund
#       
#       solution.append(nabund)
#     solution = np.vstack(solution)

    def fun(t,y,ksub0,ksub1,kdep0,kdep1):
      ksub = (t - t0)*(ksub1 - ksub0)/(t1 - t0) + ksub0
      kdep = (t - t0)*(kdep1 - kdep0)/(t1 - t0) + kdep0
      
      return (ksub + self.kcrd)*(self.max_abund - y) - kdep*y

#     def jac(t,y,ksub0,ksub1,kdep0,kdep1):
#       ksub = (t - t0)*(ksub1 - ksub0)/(t1 - t0) + ksub0
#       kdep = (t - t0)*(kdep1 - kdep0)/(t1 - t0) + kdep0
#       
#       return np.vstack((np.zeros(ksub.size),-(ksub + self.kcrd + kdep),np.zeros(ksub.size)))

    
#     solution = odeint(fun,abund,t,full_output=False,Dfun=Dfun,args=(ksub0,ksub1,kdep0,kdep1),ml=0,mu=0)
    
    r = ode(fun).set_integrator('vode',method='bdf',uband=0,lband=0,nsteps=3000)
#     r = ode(fun).set_integrator('lsoda',lband=1,uband=1,nsteps=3000)
    r.set_initial_value(abund,t0).set_f_params(ksub0,ksub1,kdep0,kdep1)#.set_jac_params(ksub0,ksub1,kdep0,kdep1)
    solution = [abund]
    success  = True
    for i in range(1,len(t)):
#       print(r.t*unit.s.to('yr'),t[i]*unit.s.to('yr'),(t[i]-t[i-1])*unit.s.to('yr'))
      r.integrate(t[i])
      solution.append(r.y)
      if not r.successful():
        success = False
    solution = np.vstack(solution)
    
#     ksub = (t[:,np.newaxis] - t0)*(ksub1[np.newaxis,:] - ksub0[np.newaxis,:])/(t1 - t0) + ksub0[np.newaxis,:]
#     kdep = (t[:,np.newaxis] - t0)*(kdep1[np.newaxis,:] - kdep0[np.newaxis,:])/(t1 - t0) + kdep0[np.newaxis,:]
    
    return solution, t*unit.s.to('yr'), success
      
#   def step(self,abund,nH,Tdust,Tgas,dt):
#     
#     abundance_out = abundance_in - abundance_in*self.kdep(Tgas,nH2)*dt
#     
#     if self.equilibrium:
#       abundance_out = (self.max_abund - abundance_out)*(self.kcrd+self.ksub(Tdust))*dt
#     else:
#       abundance_out = np.where(self.ksub(Tdust) > co.s2yr, self.max_abund, abundance_out)
#     
#     try:
#       abundance_out[abundance_out > self.max_abund] = self.max_abund
#       abundance_out[abundance_out < self.min_abund] = self.min_abund
#     except:
#       if abundance_out > self.max_abund:
#         abundance_out = self.max_abund
#       elif abundance_out < self.min_abund:
#         abundance_out = self.min_abund
# 
#     return abundance_out
