#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from cross_correlate import zucker
from scipy.interpolate import interp1d
import matplotlib.pyplot as mp
import sys

c = 299792.458  #speed of light in km/s
v = -30

wl, fl = np.loadtxt('../t7400_g360_p000.dat',unpack=True,dtype=np.float64)

bound = 130
wlf, flf = wl[bound:-bound].copy(), fl[bound:-bound].copy()

#dv, C, s_hat, s_hat_sigma, a_hat, sigma_hat, norm_factor = zucker(wlf,f,wlg,g,nchp=1,returnML=False)
#dv, ML, s_hat, s_hat_sigma, a_hat, sigma_hat, norm_factor = zucker(wlf,f,wlg,g,nchp=1,returnML=True)

n     = 500
ofac  = 2.
noise = 0.1
ofac  = 2.
nchp  = 20.

s_hat = np.empty(n)
s_hat_sigma = np.empty(n)
norm_factor = np.empty(n)
sigma_hat = np.empty(n)
a_hat = np.empty(n)
for i in range(n):
  noise = np.random.normal(size=len(wlf),scale=1.) * np.linspace(0.05,0.2,num=len(flf))
  fflf = flf + noise
  dv, C, s_hat[i], s_hat_sigma[i], a_hat, sigma_hat, norm_factor = zucker(wlf,fflf,wl,fl,nchp=nchp,ofac=ofac,returnML=True)
  #s_hat[i], s_hat_sigma[i], a_hat[i], sigma_hat[i] = np.mean(ss_hat), np.mean(ss_hat_sigma), np.mean(aa_hat), np.mean(ssigma_hat)

#print (s_hat/s_hat_sigma**2).sum()/(1/s_hat_sigma**2).sum()
#print (1/s_hat_sigma).sum()/(1/s_hat_sigma**2).sum()
print norm_factor
print 'shift     : ', np.mean(s_hat), '+/-', np.std(s_hat)
print 'sigma     : ', np.mean(s_hat_sigma), '+/-', np.std(s_hat_sigma)
#print 'sigma_hat : ', np.mean(sigma_hat), '+/-', np.std(sigma_hat)
#print 'a_hat     : ', np.mean(a_hat), '+/-', np.std(a_hat)
print 'norm_fac  : ', np.mean(norm_factor), '+/-', np.std(norm_factor)
print '------------------------------------'
print 'sigma     : ', np.mean(norm_factor)*np.mean(s_hat_sigma), '+/-', np.sqrt(np.var(s_hat_sigma)+np.var(norm_factor))
print 'sigma_hat : ', np.mean(norm_factor)*np.mean(sigma_hat), '+/-', np.sqrt(np.var(sigma_hat)+np.var(norm_factor))

mp.hist(s_hat,bins=30)
mp.show()
  