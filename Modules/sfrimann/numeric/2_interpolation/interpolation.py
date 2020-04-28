#!/usr/bin/python
# -*- coding: utf-8 -*-

#Routines for polynomial, linear and quadratic spline interpolation, and binary search.

from numpy import array, arange, exp, zeros

class RangeError(Exception):
  pass

#-----------------------------------------------------------------------------#
def binsearch(xx,x):
  #Binary search for value x in array xx. if x lies between xx[i] and xx[i+1]
  #the program returns the index i.
  
  #Determine nmin and nmac from number of elements in "xx"
  nmin = 0
  nmax = len(xx)-1
  
  #Check if searchvalue is in array
  if x < xx[nmin] or xx[nmax] < x:
    raise RangeError("Error: Search value is outside array range")
  
  #Locate position of searchvalue in array
  while (nmax - nmin) > 1:
    nmid = int((nmax + nmin)/2)
    if x >= xx[nmid]:
      nmin = nmid
    else:
      nmax = nmid

  return nmin

#-----------------------------------------------------------------------------#
def polyinterpol(xin,yin,xtab):
   #Polynomial interpolation of a tabulated dataset

   #Output array
   yout = zeros(len(xtab))

   for i in range(len(xin)):
      prod = 1.
      for j in range(len(xin)):
         if j != i:
            prod = prod * ((xtab-xin[j])/(xin[i]-xin[j]))
      yout = yout + yin[i]*prod
   
   return yout

#-----------------------------------------------------------------------------#
def lininterpol(xin,yin,xtab):
   #Linear spline interpolation
   yout = zeros(len(xtab)) #Initialize

   for i in range(len(xtab)):
     #Find out between which indicies of xin xtab[i] is located
     ii      = binsearch(xin,xtab[i])
     yout[i] = yin[ii] + ((yin[ii+1]-yin[ii])/(xin[ii+1] - xin[ii]))*(xtab[i] - xin[ii])
   
   return yout

#-----------------------------------------------------------------------------#
def quadinterpol(xin,yin,xtab):
   #Quadratic spline interpolation
   
   #Initialization
   yout = zeros(len(xtab))
   ai   = zeros(len(xin)-1)
   aii  = zeros(len(xin)-1)
   
   #Recursion up
   for i in range(len(aii)-1):
     #Defining deltas
     delx  = xin[i+1] - xin[i]
     delxx = xin[i+2] - xin[i+1]
     dely  = yin[i+1]-yin[i]
     delyy = yin[i+2]-yin[i+1]
     
     aii[i+1] = (1./delxx) * (delyy/delxx - dely/delx - aii[i]*delx)   

   #Recursion down 
   for i in range(len(ai)-2, -1,-1):
     #Defining deltas
     delx  = xin[i+1] - xin[i]
     delxx = xin[i+2] - xin[i+1]
     dely  = yin[i+1]-yin[i]
     delyy = yin[i+2]-yin[i+1]
     
     ai[i] = (1./delx) * (delyy/delxx - dely/delx - ai[i+1]*delxx)

   acomb = (ai + aii)/2.

   for i in range(len(xtab)):
     #Interpolation
     ii      = binsearch(xin,xtab[i])
     yout[i] = yin[ii] + ((yin[ii+1]-yin[ii])/(xin[ii+1] - xin[ii]))*(xtab[i] - xin[ii]) + acomb[ii]*(xtab[i]- xin[ii])*(xtab[i]-xin[ii+1])
   
   return yout