#!/usr/bin/python
# -*- coding: utf-8 -*-

# Finding Eigen Values and Vectors for real symmetric matrices

from __future__ import division
from numpy import array, float64, eye, diag
from math import atan2, cos, sin

# ---- Error definitions ----
class MatrixError(Exception):
  pass

class DimensionError(Exception):
  pass

# jacobi does cyclic Jacobi rotations to matrix A
# eps is limit for zeroing (1e-12 as standard)
# returns list of eigenvalues, E, and matrix of eigen-
# vectors, V.
def jacobi(A, eps=1e-12):

  A = array(A, dtype=float64) #Copy A
  
  #Error conditions
  if A.shape[0] != A.shape[1]:
    raise DimensionError('Error: Matrix must be square')

  if (A.T != A).any():
    print A.T == A
    raise MatrixError('Error: Matrix is not symmetric')
  
  n = A.shape[0]        #length
  V = eye(n)            #Initialize V as Idendity matrix
  nsum = 0              #Initialize nsum variable
  ssum = eps + 1        #Initialize ssum variable
  while ssum > eps:     #Condition for continuing rotations
    for p in range(n):  #Iterate above the diagonal
      for q in range(p+1,n):
        if abs(A[p,q]) > eps/n**2:  #Only rotate if value is above eps
          A,V = rotate(p,q,A,V)     #Rotation subroutine
          nsum += 1                 #Update nsum
    ssum = 0
    for i,row in enumerate(A):      #Check if summed elements are below eps
      ssum += sum(abs(row[i+1:]))
  
  E = diag(A)  #Get eigenvalues from diagonal of A
  
  return E, V, nsum

#Subroutine for jacobi
#Given p and q it rotates A, and updates V
def rotate(p,q,A,V):
  #Initial values
  n = A.shape[0]; App = A[p,p]; Aqq = A[q,q]; Apq = A[p,q]
  phi = 0.5*atan2(2*Apq, Aqq-App)   #Find phi
  c = cos(phi); s = sin(phi)        #Calculate sin and cos
  
  A[p,p] = c*c*App + s*s*Aqq - 2*s*c*Apq  #Update diagonal
  A[q,q] = s*s*App + c*c*Aqq + 2*s*c*Apq
  A[p,q] = 0  #This is zero by construction
  
  #Iterate over and update remaining off-diagonal elements
  for i in range(p):
    Aip = A[i,p]; Aiq = A[i,q]
    A[i,p] = c*Aip - s*Aiq
    A[i,q] = c*Aiq + s*Aip
    
  for i in range(p+1,q):
    Api = A[p,i]; Aiq = A[i,q]
    A[p,i] = c*Api - s*Aiq
    A[i,q] = c*Aiq + s*Api
  
  for i in range(q+1,n):
    Api = A[p,i]; Aqi= A[q,i]
    A[p,i] = c*Api - s*Aqi
    A[q,i] = c*Aqi + s*Api
  
  #Update V
  for i in range(n):
    Vip = V[i,p]; Viq = V[i,q]
    V[i,p] = c*Vip - s*Viq
    V[i,q] = s*Vip + c*Viq
  
  return A, V
