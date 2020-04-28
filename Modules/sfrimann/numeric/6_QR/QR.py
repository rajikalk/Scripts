#!/usr/bin/python
# -*- coding: utf-8 -*-

#Solution of linear equations by QR decomposition

from __future__ import division
from math import sqrt
import numpy as np



# ---- Error Definitions
class Singular(Exception):
  pass

class DimensionError(Exception):
  pass



# ---- QR Decomposition Routine
# A must be a multidimensional array with m >= n where m is number of rows and n is number of columns
# tol is a user specified tolerance level, defining the limit for when a matrix is considered singular
def qrdec(A,tol=0):              #QR Decomposition
  A = A.astype(np.float64)       #Convert A to Double
  
  m,n = A.shape                  #Check Shape
  
  if m < n:                      #Check Error Condition
    raise DimensionError('Error: Number of rows must be >= number of columns.')
  
  R = np.zeros((n,n),np.float64) #Initialize R and Q
  Q = np.empty((m,n),np.float64)
  
  for i in range(n):             #Outer loop
    R[i,i] = sqrt(np.dot(A[:,i],A[:,i])) #Calculate R diagonals
    if R[i,i] < tol:             #Check if singular
      raise Singular('Error: Input matrix is probably singular.')
    
    Q.T[i] = A.T[i] / R[i,i]     #Orthogonalization
    for j in range(i+1,n):       #Inner loop
      R[i,j]  = np.dot(Q[:,i],A[:,j]) #Off-diagonal R
      A.T[j] -= Q.T[i]*R[i,j]    #Remove everything that points in direction already orthogonalized
  
  return Q, R                    #Return Q and R



# ---- QR Back Substitution Routine
# Q and R are output from qrdec
# b is the vector on the right hand side
def qrback(Q,R,b):
  b = b.astype(np.float64)        #Convert b to Double
  
  if len(b.shape) > 1:            #Check dimension
    raise DimensionError('Error: b must be a vector. Use qrsystem if you wish to solve several systems at once')
  
  c = np.dot(Q.T,b)               #Calculate c
  x = np.zeros(len(c),np.float64) #Initialize x
  
  for i in range(len(x)-1,-1,-1): #Back substitution loop
    x[i] = (c[i] - np.dot(R[i],x))/R[i,i]
  
  return x



# ---- QR Back Substitution Wrapper for several solutions
# Q and R are output from qrdec
# B is multidimensional array on the right hand side
def qrsystem(Q,R,B):
  B = np.array(B)                #Convert to numpy array
  
  size1 = B.shape                #Get shape

  if len(size1) == 1:            #if vector just run qrback
    return qrback(Q,R,B)
  
  size2 = R.shape                #Shape of R
  
  if size1[0] != size2[0]:       #Check if B is compatiable with R
    raise DimensionError('Error: B must have same height as R')
  
  X = np.empty(size1,np.float64) #Initialize
  for i in range(size1[1]):      #Loop over Columns and call qrback
    X.T[i] = qrback(Q,R,B.T[i])
  
  return X



# ---- Inverse matrix calculation using qrsystem
def inverse(A):
  m,n = A.shape
  
  if m != n:
    raise DimensionError('Error: A is not square')
  
  Q,R = qrdec(A)
  I   = np.eye(n)
  
  return qrsystem(Q,R,I)



# ---- Calculation of Absolute Determinant using R
def absdet(R):
  return np.prod(np.diag(R))