#!/usr/bin/python
# -*- coding: utf-8 -*-

#Library of covariance functions. 
#Also containes helper-function cvm that calculates covariance matrices

import numpy as np
from math import sqrt

#Calculates Covariance matrix of matrix with dimensions len(y) x len(x).
#x and y are vectors with input values. **kw should contain hyperparameters
def cvm(hyp,covar,y,x):
  if np.size(y) > 1 and np.size(x) > 1:
    xx, yy = np.meshgrid(x,y)
    return covar(hyp,xx,yy)
  else:
    return covar(hyp,x,y)

#Squared Exponential covariance function: a^2*exp(-(x-xp)^2/(2*ell^2))
#Takes two hyperparameters: ln(a) and ln(ell)
#Logarithms are used when we deal with scale parameters (which are natually positive)
def SE(hyp,x,xp):
  a2 = np.exp(2*hyp[0]); ell2 = np.exp(2*hyp[1])
  return a2*np.exp(-(x - xp)**2/(2.*ell2))

def expon(hyp,x,xp):
  a2 = np.exp(2*hyp[0]); ell = np.exp(hyp[1])
  return a2*np.exp(-np.abs(x - xp)/ell)

def matern32(hyp,x,xp):
  ell = np.exp(hyp[0])
  return (1 + sqrt(3)*np.abs(x - xp)/ell)*np.exp(-sqrt(3.)*np.abs(x - xp)/ell)

#Periodic covariance function a^2*exp(-2sin(pi(x-xp)/period)^2/ell^2)
#Takes three hyperparameters: ln(period), ln(a), and ln(ell)
#Logarithms are used when we deal with scale parameters (which are natually positive)
def period(hyp,x,xp):
  per = np.exp(hyp[0]); a2 = np.exp(2*hyp[1]); ell2 = np.exp(2*hyp[2])
  return a2*np.exp(-2*np.sin(np.pi*(x - xp)/per)**2/ell2)

def RQ(hyp,x,xp):
  a2 = np.exp(2*hyp[0]); ell2 = np.exp(2*hyp[1]); shape = np.exp(2*hyp[2])
  return a2*(1 + (x-xp)**2/(2*shape*ell2))**(-shape)