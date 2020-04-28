#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as mp
from random import normalvariate
from LSQ import lsfit

# --- Fit 1 ---

x  = np.arange(-10,10.01,1)   #Create abscissa value
y  = 3*x - 0.5*x**2           #Create data values
dy = np.empty(len(y))         #Initialize errors

std = 10                      #Standard Deviation for errors

for i in range(len(y)):       #Fill errors
  dy[i] = normalvariate(0,std)

y += dy                       #Data with errors

fun = [lambda x: x, lambda x: x**2] #Array of functions

c, error, covar, yfit = lsfit(x,y,dy,fun) #Call to lsfit

yfit1 = np.zeros(len(x))
yfit2 = np.zeros(len(x))
for i,f in enumerate(fun):
  yfit1 += f(x)*(c[i] - error[i])
  yfit2 += f(x)*(c[i] + error[i])

#Print results
print ''
print 'Model: c[0]*x + c[1]*x^2 (actual data follows y = 3x - 0.5x^2)'
for i in range(len(c)):
  print 'c[%1i] = %.5f +/- %.5f' % (i,c[i],error[i])

mp.figure(1)
mp.errorbar(x, y, yerr=dy, xerr=None,fmt='o',mfc='blue',ms=3)
mp.plot(x,yfit,'r')
#mp.plot(x,yfit1,'g')
#mp.plot(x,yfit2,'k')
#mp.legend( ('Estimated Error','Actual Error','$1/\sqrt{N}$'), loc='upper right')
#mp.suptitle('Error convergence as function of N')
mp.title('Fit of $f(x) = ax + bx^2$')
mp.xlabel('x')
mp.ylabel('f(x)')
mp.savefig('fit1.png', format='png')

# --- Fit 2 ---

x  = np.arange(0,10.01,0.5)   #Create abscissa value
y  = 1.2*x + 2.3              #Create data values
dy = np.empty(len(y))         #Initialize errors

std = 30                      #Standard Deviation for errors

for i in range(len(y)):       #Fill errors
  dy[i] = normalvariate(0,std)

y += dy                       #Data with errors

fun = [lambda x: x, lambda x: 1] #Array of functions

c, error, covar, yfit = lsfit(x,y,dy,fun) #Call to lsfit

yfit1 = np.zeros(len(x))
yfit2 = np.zeros(len(x))
for i,f in enumerate(fun):
  yfit1 += f(x)*(c[i] - error[i])
  yfit2 += f(x)*(c[i] + error[i])

#Print results
print ''
print 'Model: c[0]*x + c[1] (actual data follows y = 1.2x + 2.3)'
for i in range(len(c)):
  print 'c[%1i] = %.5f +/- %.5f' % (i,c[i],error[i])

mp.figure(2)
mp.plot(x,yfit,'r')
mp.plot(x,yfit1,'g')
mp.plot(x,yfit2,'k')
mp.errorbar(x, y, yerr=dy, xerr=None,fmt='o',mfc='blue',ms=3)
mp.title('Fit of $f(x) = ax + b$')
mp.legend( ('Best Fit', 'Lower Bound (1$\sigma$)', 'Upper Bound (1$\sigma$)'))
mp.xlabel('x')
mp.ylabel('f(x)')
mp.savefig('fit2.png', format='png')
