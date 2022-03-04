#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import log, sqrt, pi, cos
from numpy import zeros, ones, arange, median, abs, round, dot
import matplotlib.pyplot as mp
from plainmc import plainmc

class OddError(Exception):
  pass

#Find the indicy of an array corresponding to the median.
#If array has even number of elements return the indicies of the two middle elements.
#Array must be numpy array.
def argmedian(array):
  length = len(array)
  srt    = array.argsort()
  if length % 2 == 0:
    return srt[length/2-1:length/2+1]
  else:
    return srt[length/2]

# ---- Integral of (1-cos(x)cos(y)cos(z))^-1*pi^-3 ------------------------------------- #
def func1(x):                                               #Define function
  return 1/((1 - cos(x[0])*cos(x[1])*cos(x[2]))*pi**3)

a   = zeros(3, dtype='float')                               #Limits
b   = pi*ones(3, dtype='float')
N   = round(10**arange(1,5.1,0.1))                          #Number of points (equally spaced on logaritmic axis)
av  = 21                                                    #Number of samples from which to draw out median (must be odd)
res = 1.3932039296856768591842462603255                     #Analytical Result

if av % 2 == 0:
  raise OddError("Error: av must be odd number to in order to get true median")

integral        = zeros(av)                                 #Initialize
esterror        = zeros(av)
median_integral = zeros(len(N))
median_esterror = zeros(len(N))


for i in range(len(N)):                                     #Integrate with different N
  for j in range(av):                                       #Integrate av times for same N to avoid outliners
    [integral[j], esterror[j]] = plainmc(func1,a,b,N[i])    #Integrate
  median_integral[i] = integral[argmedian(integral)]        #Pull out medians
  median_esterror[i] = esterror[argmedian(integral)]

#Plot results
mp.figure(1)
mp.loglog(N,median_esterror,'ob')
mp.loglog(N,abs(median_integral-res),'or')
mp.loglog(N,N**(-0.5),'k')
mp.legend( ('Estimated Error','Actual Error','$1/\sqrt{N}$'), loc='upper right')
mp.suptitle('Error convergence as function of N')
mp.title('$(1-\cos(x)\cos(y)\cos(z))^{-1} \cdot \pi^{-3}$')
mp.xlabel('N')
mp.ylabel('Error')
mp.savefig('MCerror_cos.png', format='png')

#Print results
print ''
print 'Integral: (1-cos(x)cos(y)cos(z))^-1*pi^-3 | Analytical result = 1.3932039296856768591842462603255'
for i in range(len(N)):
  print 'N = %6i | Result = %.5f | Estimated error = %.5f | Actual error = %.5f' % (N[i], median_integral[i], median_esterror[i], abs(median_integral[i]-res))
  



# ---- Volume of x^2 + y^2 + z^2 <= 1 ------------------------------------- #
def func2(x):                                               #Define function
  if dot(x,x) <= 1:
    return 1.
  else:
    return 0.

a   = -1*ones(3, dtype='float')                               #Limits
b   = ones(3, dtype='float')
N   = round(10**arange(1,5.1,0.1))                          #Number of points (equally spaced on logaritmic axis)
av  = 21                                                    #Number of samples from which to draw out median (must be odd)
res = pi*4./3,                                              #Analytical Result

if av % 2 == 0:
  raise OddError("Error: av must be odd number to in order to get true median")

integral        = zeros(av)                                 #Initialize
esterror        = zeros(av)
median_integral = zeros(len(N))
median_esterror = zeros(len(N))


for i in range(len(N)):                                     #Integrate with different N
  for j in range(av):                                       #Integrate av times for same N to avoid outliners
    [integral[j], esterror[j]] = plainmc(func2,a,b,N[i])    #Integrate
  median_integral[i] = integral[argmedian(integral)]        #Pull out medians
  median_esterror[i] = esterror[argmedian(integral)]

#Plot results
mp.figure(2)
mp.loglog(N,median_esterror,'ob')
mp.loglog(N,abs(median_integral-res),'or')
mp.loglog(N,N**(-0.5),'k')
mp.legend( ('Estimated Error','Actual Error','$1/\sqrt{N}$'), loc='upper right')
mp.suptitle('Error convergence as function of N')
mp.title('$x^2 + y^2 + z^2 \leq 1$')
mp.xlabel('N')
mp.ylabel('Error')
mp.savefig('MCerror_sphere.png', format='png')

#Print results
print ''
print 'Volume of x^2 + y^2 + z^2 <= 1 | Analytical result = %.5f' % res
for i in range(len(N)):
  print 'N = %6i | Result = %.5f | Estimated error = %.5f | Actual error = %.5f' % (N[i], median_integral[i], median_esterror[i], abs(median_integral[i]-res))