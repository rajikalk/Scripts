#!/usr/bin/python
# -*- coding: utf-8 -*-

from numpy import array, ones, zeros, dot, float64
from math import sqrt
import sys
sys.setrecursionlimit(1200)


def adaptint(func,a,b,acc,eps,oldfs=None):
  
  x  = array([1./6, 2./6, 4./6, 5./6],dtype='float64') #Abscissas
  wh = array([2./6, 1./6, 1./6, 2./6],dtype='float64') #Weights of higher order quadrature
  wl = 1./4*ones(4,dtype='float64')                    #Weights of lower order quadrature
  p  = [True, False, False, True]                      #New points to be calculated for each recursion
  n  = len(x)                                          #Number of abscissas
  h  = float64(b - a)                                  #Integration interval
  a  = float64(a)
  b  = float64(b)
  
  if oldfs == None:                    #First call
    fs = array([ func(a+i*h) for i in x])
  else:
    fs = zeros(n)
    k  = 0
    for i in range(n):
      if p[i] == True: fs[i] = func(a+x[i]*h)   #Make new points
      else:            fs[i] = oldfs[k]; k += 1 #Reuse old points

  q4  = h*dot(wh,fs)       #Four point estimate
  q2  = h*dot(wl,fs)       #Two point estimate
  tol = acc + eps*abs(q4)  #Tolerance
  err = abs(q4 - q2)/(n-1) #Estimated error
  
  if err < tol:            #Error estimate fulfilled, return
    return [q4, err, 1]
  else:                    #Error estimate not fulfilled, recursion
    acc   = acc/sqrt(2)    #Rescale accuracy goal
    mid   = float64((a+b)/2) #Mid-point
    left  = fs[:2]         #Store left points
    right = fs[2:]         #Store right points
    
    #Recursive calls
    [ql, erl, ml] = adaptint(func,a,mid,acc,eps,oldfs=left)
    [qr, err, mr] = adaptint(func,mid,b,acc,eps,oldfs=right)
    
    return [ql+qr, sqrt(erl**2 + err**2), ml+mr+1]