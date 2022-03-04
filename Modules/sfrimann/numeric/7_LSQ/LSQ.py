#!/usr/bin/python
# -*- coding: utf-8 -*-

# Linear least square fitter

from __future__ import division
from numpy import array, float64, empty, dot, sqrt, zeros, diag
import sys
sys.path.append('../6_QR')
from QR import qrdec, qrback, inverse

# ---- Error definitions ----
class DimensionError(Exception):
  pass



# ---- Wrapper for lsfit_raw ----
# x is abscissa values. y is data values. dy is data errors, and fun is array of functions.
# keyword w is statistical weights - if present dy input is ignored. If not preset w is set to 1/dy.
# w and dy has to be either scalars (uniform weights/errors) or len(y) element vectors.
# returns array, c, of coefficients, standard deviation of coefficients, covariance matrix, and vector of y values for the best fit
def lsfit(x,y,dy,fun,w=None):
  
  x  = array(x , dtype=float64)                 #Convert inputs to double
  y  = array(y , dtype=float64)
  dy = array(dy, dtype=float64)
  
  if w == None:                                 #Check if w is set
    if len([dy]) != 1 and len([dy]) != len(y):  #Error condition
      raise DimensionError('Error: dy must either have one element (uniform error), or same number of elements as y.')
    w = 1/dy                                    #Setting Weight
  else:
    w = array(w, dtype=float64)
    if len([w]) != 1 and len([w]) != len(y):    #Error condition
      raise DimensionError('Error: w must either have one element (uniform weight), or same number of elements as y.')
  
  c, covar = lsfit_raw(x,y,w,fun)               #Run lsfit_raw
  
  error = sqrt(diag(covar))                  #Get errors from diagonal of covariance matrix
  
  yfit = zeros(len(x))
  for i,f in enumerate(fun): yfit += f(x)*c[i]  #y-values of best fit
  
  return c, error, covar, yfit

# ---- Raw least square routine ----
# x is abscissa values. y is data values. w is weights, and fun is array of functions
# returns array, c, of coefficients to be coupled with the functions, and covariance matrix
def lsfit_raw(x,y,w,fun):
  A = empty([len(x),len(fun)])    #Initialize A
  
  for i,f in enumerate(fun):      #Fill in A (note that the functions must be able to take vectors)
    A[:,i] = f(x)*w
    
  b = y*w                         #Make b

  Q,R = qrdec(A)                 #QR decomposition
  c = qrback(Q,R,b)              #Back-substitute for c
  
  covar = inverse(dot(R.T,R))    #Covariance matrix
  
  return c, covar                #Return