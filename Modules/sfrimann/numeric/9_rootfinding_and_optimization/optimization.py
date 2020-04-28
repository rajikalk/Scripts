#!/usr/bin/python
# -*- coding: utf-8 -*-

# Implementation of downhill simplex (amoeba)

from __future__ import division
from numpy import dot, array, float64, sum, argmin, argmax
from math import sqrt
import sys
sys.path.append('../6_QR')
from QR import qrdec, qrback

# ---- Error definitions ----
class DimensionError(Exception):
  pass

#Downhill simplex routine.
#fun is the function to be minimized
#p is the starting simplex. It must be an array with dimensions n x (n+1), where n is the number of parameters taken by fun
#acc is the desired accuracy of the minimization (default at 1e-5)
def amoeba(fun,p,acc=1e-6):
  
  p = array(p, float64)          #make starting simplex into numpy array and use double precision
  
  if p.shape[0] != p.shape[1]-1: #check error condition
    raise DimensionError('Error: p must have dimensions n x (n+1)')
  
  fval = array([fun(p[:,i]) for i in range(p.shape[1])]) #calculate initial function values
  iter = 0                       #initialize iteration
  while size(p) > acc:           #main loop with stopping criterium
    iter += 1                    #step forward in iteration
    
    #Find indexes for min and max simplex as well as their values
    mn = argmin(fval); mx = argmax(fval); plo = p[:,mn]; phi = p[:,mx]
    
    pce = centroid(p,mx,p.shape[1]) #calculate centroid
    pre = pce + (pce - phi)         #calculate reflection
    if fun(pre) < fun(phi):         #try reflection
      fval[mx] = fun(pre)
      p[:,mx] = pre
      if fun(pre) < fun(plo):       #f(reflection) < f(lowest) try expansion
        pex = pce + 2*(pce - phi)   #calculate expansion
        if fun(pex) < fun(pre):     #try expansion
          fval[mx] = fun(pex)
          p[:,mx] = pex
    else:
      pco = pce + 0.5*(phi - pce)   #if reflection didn't work try contraction
      if fun(pco) < fun(phi):
        fval[mx] = fun(pco)
        p[:,mx] = pco
      else:
        for i in range(p.shape[1]): #if nothing worked do reduction
          if i != mn:
            p[:,i] = 0.5*(p[:,i] + p[:,mn])
            fval[i] = fun(p[:,i])
  
  #return minimum simplex
  return p[:,mn], iter

#Definition of auxiliary functions
def centroid(p,mx,n):
  if mx == 0:
    pce = sum(p[:,1:],1)/(n-1)
  elif mx == (n-1):
    pce = sum(p[:,:n-1],1)/(n-1)
  else:
    pce = (sum(p[:,:mx],1) + sum(p[:,mx+1:],1))/(n-1)
  return pce

def norm(v):
  return sqrt(dot(v,v))

def dist(v1,v2):
  return norm(v1 - v2)

def size(p):
  return norm([dist(p[:,i],p[:,0]) for i in range(1,p.shape[1])])