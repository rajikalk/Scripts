#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import log, sqrt, pi
from adaptint import adaptint

def func1(x):
  return log(x)/sqrt(x)

ran = [0, 1]
acc = 0.001
eps = 0.001
[q, err, m] = adaptint(func1,ran[0],ran[1],acc,eps)

print ''
print '         Integral: ln(x)/sqrt(x) in range %01d to %01d' %(ran[0], ran[1])
print '           Result: %f' % q
print '            Error:  %.16f' % err
print '   Required Error:  %.16f' % eps
print '  Number of calls:  %d'  % m

def func2(x):
  return 4*sqrt(1-(1-x)**2)

ran = [0, 1]
acc = 0.0000000000001
eps = 0.0000000000001
[q, err, m] = adaptint(func2,ran[0],ran[1],acc,eps)

print ''
print '         Integral: 4*sqrt(1-(1-x)^2) in range %01d to %01d' %(ran[0], ran[1])
print '           Result: %.16f' % q
print '               Pi: %.16f' % pi
print '            Error: %.16f' % err
print '   Required Error: %.16f' % eps
print '  Number of calls: %d'  % m

def func3(x):
  return x**(-x)

ran = [0, 1]
acc = 0.0000000000001
eps = 0.0000000000001
[q, err, m] = adaptint(func3,ran[0],ran[1],acc,eps)

print ''
print ' Sophomores Dream 1'
print '         Integral: x^-x in range %01d to %01d' %(ran[0], ran[1])
print '           Result: %.16f' % q
print ' Wikipedia Result: 1.291285997'
print '            Error: %.16f' % err
print '   Required Error: %.16f' % eps
print '  Number of calls: %d'  % m

def func4(x):
  return x**x

ran = [0, 1]
acc = 0.0000000000001
eps = 0.0000000000001
[q, err, m] = adaptint(func4,ran[0],ran[1],acc,eps)

print ''
print ' Sophomores Dream 2'
print '         Integral: x^-x in range %01d to %01d' %(ran[0], ran[1])
print '           Result: %.16f' % q
print ' Wikipedia Result: 0.783430510712'
print '            Error: %.16f' % err
print '   Required Error: %.16f' % eps
print '  Number of calls: %d'  % m