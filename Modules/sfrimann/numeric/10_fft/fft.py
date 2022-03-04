#!/usr/bin/python
# -*- coding: utf-8 -*-

# Implementation of fast fourier transform

from __future__ import division
from cmath import pi, exp
from math import sqrt
from numpy import concatenate, arange, array, float64, arange


# ---- Error definitions ----
class InputError(Exception):
  pass

#fft_core function
#Not meant for direct calling. Use the wrapper fft.
#x is the input array which needs to be numpy, and *must* have lenght equal to a power of 2.
#sign gives either forward or inverse fft (-1 for forward, +1 for inverse)
#note that the transform is unitary so you don't need to normalize when doing the inverse
def fft_core(x,sign):
  
  N = len(x)
  if N != 1:
    #Recursive calls. Divide and conquer
    ceven = fft_core(x[::2] ,sign) #Even fft
    codd  = fft_core(x[1::2],sign) #Odd fft
    
    M = int(N/2)
    
    w = exp(sign*2*pi*1j/N)

    #Concatenate and normalize
    return concatenate((ceven + w**arange(M)*codd, ceven - w**arange(M)*codd))/sqrt(2)
  else: #when we're down to the one-element level
    return x

def fft(x,sign=-1,dx=None):
  
  if abs(sign) != 1: #if sign is something else than 1 scale it to 1
    sign /= abs(sign)
  
  if sign == -1: #When we do reverse fft x has an extra element. 
    N = len(x)
  else:
    N = len(x)-1
  
  if (N & (N - 1)) != 0: #check if power of two
    raise InputError('Error: Input vector *must* have a length that is a power of two. Zero pad your data!')
  
  if sign == 1: #If we do the inverse recreate correct format for x
    x = concatenate((x[N/2:],x[1:N/2]))
  
  ft = fft_core(x,sign)
  
  if sign == -1: #make format that makes sense if you want to plot frequency
    ft = concatenate((ft[N/2:],ft[:N/2+1]))
  
  if dx != None: #if dx is given make frequency vector
    freq = arange(-N/2,N/2+1)/(N*dx)
    return ft, freq
  else:
    return ft
  