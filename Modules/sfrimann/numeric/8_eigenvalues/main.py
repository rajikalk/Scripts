#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from jacobi import jacobi

A1 = np.array([[1,-2,4],[-2,3,5],[4,5,1]])
E1,V1,nsum1 = jacobi(A1)

A2=np.array([[3,4,-6,2],[4,1,3,-8],[-6,3,-9,1],[2,-8,1,7]],dtype=np.float64)
E2,V2,nsum2 = jacobi(A2)

print "Example 1. 3x3 matrix"
print "Original matrix:"
print A1
print " "
print "Number of summations: %3i" % nsum1
print " "
print "Eigenvalues:"
print E1
print " "
print "Eigenvectors (columns):"
print V1
print " "
print "Test that this solves the eigensystem (Evaluating A*V*diag(E)^-1)"
print np.dot(np.dot(A1,V1),np.diag(1/E1))

print ""
print "Example 2. 4x4 matrix"
print "Original matrix:"
print A2
print " "
print "Number of summations: %3i" % nsum2
print " "
print "Eigenvalues:"
print E2
print " "
print "Eigenvectors (columns):"
print V2
print " "
print "Test that this solves the eigensystem (Evaluating A*V*diag(E)^-1)"
print np.dot(np.dot(A2,V2),np.diag(1/E2))
