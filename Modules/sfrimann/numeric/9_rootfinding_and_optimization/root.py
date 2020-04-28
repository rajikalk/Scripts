#!/usr/bin/python
# -*- coding: utf-8 -*-

# Implementation of Newton's method for root finding

from __future__ import division
from numpy import empty, dot, array, float64
from math import sqrt
import sys
sys.path.append('../6_QR')
from QR import qrdec, qrback

# ---- Error definitions ----
class DimensionError(Exception):
  pass

#fun is array of functions
#x is array of starting parameters (len(x) must be equal to len(fun))
#acc is tolerance for solution
#dx is finite difference used when approximating Jacobi matrix
def newton(fun,x,acc=1e-6,dx=1e-3):
  
  x = array(x, float64) #Make copy of starting guess and use double precision
  
  if len(fun) != len(x):   #Error condition
    raise DimensionError('Error: Number of functions must equal number of parameters!')
  
  J = empty((len(fun),len(x)))   #Initialize Jacobi matrix
  minusfx = empty(len(x))        #Initialize minusfx
  
  for i,f in enumerate(fun):     #Fill in minusfx
    minusfx[i] = -f(x)
  
  iter = 0                       #Iteration initialization
  while sqrt(dot(minusfx,minusfx)) > acc: #Outer while loop. Checking stopping condition
    
    iter += 1                    #Iterate up
    
    for i,f in enumerate(fun):            #Evaluate Jacobi matrix using finite difference
      for j in range(len(x)):
        x[j] += dx
        J[i,j] = (f(x) + minusfx[i])/dx
        x[j] -= dx

    Q,R = qrdec(J); Dx = qrback(Q,R,minusfx) #Solve system using QR-Decomposition to find Newton step
    
    s = 2   #Initialize s
    minusfz = minusfx.copy()  #Initialize minusfz
    #Inner while loop. Adjust Newton step
    while sqrt(dot(minusfz,minusfz)) > (1-s/2)*sqrt(dot(minusfx,minusfx)) and s > 1/128:
      s /= 2
      z = x + s*Dx
      for i,f in enumerate(fun):
        minusfz[i] = -f(z)
    
    minusfx = minusfz; x = z  #Set step
  
  return x, iter  #Return root