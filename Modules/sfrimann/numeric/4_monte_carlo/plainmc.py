#!/usr/bin/python
# -*- coding: utf-8 -*-

#Implementation of plain multidimensional Monte Carlo Integration

from numpy import array
from math import sqrt
from random import uniform

#Uniform random number in range between a and b
def randompar(a,b):
  return array([ uniform(a[i],b[i]) for i in a ])

#Multidimensional Monte Carlo Integration
#func is the function to be integrated
#a and b are arrays with the integration limits (length of the array is equal to the number of dimensions)
#N is the number of sample points
def plainmc(func,a,b,N):
  
  vol = 1.
  for i in a:
    vol *= b[i] - a[i]              #Calculate volume
  
  sum1 = 0.
  sum2 = 0.
  for i in range(N):                #Calculate sums
    f = func(randompar(a,b))
    sum1 += f
    sum2 += f**2
  
  average  = sum1/N
  variance = sum2/N - average**2
  integral = vol*average            #Integral estimate
  error    = vol*sqrt(variance/N)   #Error estimate
  
  return [integral, error]