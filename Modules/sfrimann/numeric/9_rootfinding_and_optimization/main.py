#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from random import normalvariate
from root import newton
from optimization import amoeba

# --- Example 1: Root Finding using Newton's method ---

#Coupled functions of which to find the root:
fun = [lambda x: 1-x[0], lambda x: 10*(x[1]-x[0]*x[0])]

#Initial guess for the parameters
x = np.random.uniform(-10,10,len(fun))

#Find roots
root, iter = newton(fun,x)

#Print Results
print ''
print '           Example 1: Root Finding of the coupled equations (1-x)=0 and 10(y-x^2)=0'
print '   Analytical result: x =  1.00        , y =  1.00'
print '       Initial guess: x = % .10f, y = % .10f' % (x[0], x[1])
print '              Result: x = % .10f, y = % .10f' % (root[0], root[1])
print 'Number of iterations: %2i' % iter


# --- Example 2: Minimization of Rosenbrock function using downhill simplex ---

#Function to be minimized
fun = lambda x: (1-x[0])**2+100*(x[1]-x[0]*x[0])**2

#Initial simplex
init = np.random.uniform(-100,100,(2,3))

#Find minimum
res, iter = amoeba(fun,init)

#Print Results
print ''
print '           Example 2: Minimization of Rosenbrock function = (1-x)^2 + 100(y-x^2)^2'
print '   Analytical result: x =  1.00        , y =  1.00'
print '              Result: x = % .10f, y = % .10f' % (res[0], res[1])
print 'Number of iterations: %2i' % iter


# --- Example 3: Minimization of Himmelblau's function using downhill simplex ---

#Function to be minimized
fun = lambda x: ((x[0]**2+x[1]-11)**2+(x[0]+x[1]**2-7)**2)

#Initial simplex
init = np.random.uniform(-100,100,(2,3))

#Find minimum
res, iter = amoeba(fun,init)

#Print Results
print ''
print '           Example 3: Minimization of Himmelblau\'s function = (x^2 + y - 11)^2 + (x + y^2 -7)^2'
print '  Analytical results: x =  3.000000 , y =  2.000000 (Source: Wikipedia)'
print '                    : x = -2.805118 , y =  3.131312'
print '                    : x = -3.779310 , y = -3.283186'
print '                    : x =  3.584428 , y = -1.848126'
print ''
print '              Result: x = % .6f , y = % .6f' % (res[0], res[1])
print 'Number of iterations: %2i' % iter
