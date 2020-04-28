#!/usr/bin/python
# -*- coding: utf-8 -*-

from numpy import array, dot
from QR import qrdec, qrback, qrsystem, inverse, absdet

# ---- Example -----

tolerance = 0.000001		#set Tolerance

A = array([[1,4,7,3,4],		#Input matrix
	   [3,5,2,4,1],
	   [8,4,1,7,2],
	   [9,4,8,2,8],
	   [8,1,4,4,6]])

b = array([4,2,7,1,0])		#Solution
print A
Q,R = qrdec(A,tol=tolerance)	#QR Decomposition
res = qrback(Q,R,b)		#QR Back Substitution
INV = inverse(A)	   	#Inverse
print A
print dot(Q,R)

# ----- Tests -----

print 'Solution of system: Ax = b'
print ''
print 'Printing matrices'
print ''
print 'This is the R matrix:'
print R
print ''
print 'This is the Q matrix:'
print Q
print ''
print 'Q^TQ is the idendity matrix:'
print dot(Q.T,Q).round()
print ''
print 'Printing results'
print ''
print 'Solution of the system:'
print res
print ''
print 'Test if this is a solution to the problem:'
print dot(A,res)
print 'Here is b for reference:'
print b
print ''
print 'Inverse of A:'
print INV
print ''
print 'Test that this is the inverse of A:'
print dot(INV,A).round()
print ''
print 'Absolute determinant of A:'
print absdet(R)
