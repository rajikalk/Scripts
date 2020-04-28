#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import e
from numpy import empty
import matplotlib.pyplot as mp
from ODE import rkdrive

# ---- Example 1: Exponential

def func1(x,y):
  return y #Defining ODE for e

a   = 0
b   = 1
y0  = 1
acc = 0.00001
eps = 0.00001
h   = 0.00001

xlist, ylist = rkdrive(func1,a,b,y0,acc,eps,h)

#print results
print ''
print 'Equation: dy/dx = y in range 0 to 1'
print 'Analytical result = e = %.5f' % e
print '          rk45 result = %.5f' % ylist[-1]
print ' Number of iterations = %5i' % len(ylist)

# ---- Example 2: Damped oscillator

def func2(x,y,c):
  yp0 = y[1]
  yp1 = -c[0]*y[1] - y[0]
  yp2 = y[3]
  yp3 = -c[1]*y[3] - y[2]
  yp4 = y[5]
  yp5 = -c[2]*y[5] - y[4]
  yp6 = y[7]
  yp7 = -c[3]*y[7] - y[6]
  yp8 = y[9]
  yp9 = -c[4]*y[9] - y[8]
  return [yp0, yp1, yp2, yp3, yp4, yp5, yp6, yp7, yp8, yp9]

a   = 0
b   = 10
y0  = 5*[0, 1]
acc = 0.00001
eps = 0.00001
h   = 0.00001
c   = [0.3, 0.9, 1.5, 2.1, 2.9]

xlist, ylist = rkdrive(func2,a,b,y0,acc,eps,h,Extra=c)

y = empty((5, len(ylist)))
for i, arr in enumerate(ylist):
  y[:,i] = arr[1::2]

#Plot results
mp.figure(1)
for i in range(5):
  mp.plot(xlist,y[i].T)
mp.legend( ('c = 0.3','c = 0.9','c = 1.5','c = 2.1','c = 2.9'), loc='upper right')
mp.suptitle('Damped Oscillator')
mp.title('$d^2y/dt^2 = -c dy/dt - y$')
mp.xlabel('t')
mp.ylabel('y')
mp.savefig('damped_oscillator.png', format='png')