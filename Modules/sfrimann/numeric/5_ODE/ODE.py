#!/usr/bin/python
# -*- coding: utf-8 -*-

#Solution of ordinary differential equations

from __future__ import division
from math import sqrt
from numpy import array

#The Extra keyword is present in case the function takes extra (constant) parameters
def rk45step(fun,x,y,h,Extra=None): #Embedded fifth order Runge-Kutta Method. Stepper function
  # ---- Dormand-Prince Parameters for embedded Runge-Kutta method (Press et. al 2007, p930)
  c = [0, 1/5, 3/10, 4/5, 8/9, 1]
  
  a1 = [1/5                                                       ]
  a2 = [3/40      ,  9/40                                         ]
  a3 = [44/45     , -56/15     , 32/9                             ]
  a4 = [19372/6561, -25360/2187, 64448/6561, -212/729             ]
  a5 = [9017/3168 , -355/33    , 46732/5247,  49/176 , -5103/18656]
  
  b    = [35/384    , 0, 500/1113  , 125/192, -2187/6784   , 11/84   ]
  bemb = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100]

  # ---- Implementation (Press et. al 2007, eq. 17.2.4)
  _fun = fun
  def fun(*x):
    if x[-1] is None:
      return array(_fun(*x[:-1]))
    else:
      return array(_fun(*x))
    
  k0 = h*fun(x, y, Extra)
  k1 = h*fun(x + c[1]*h, y + a1[0]*k0, Extra)
  k2 = h*fun(x + c[2]*h, y + a2[0]*k0 + a2[1]*k1, Extra)
  k3 = h*fun(x + c[3]*h, y + a3[0]*k0 + a3[1]*k1 + a3[2]*k2, Extra)
  k4 = h*fun(x + c[4]*h, y + a4[0]*k0 + a4[1]*k1 + a4[2]*k2 + a4[3]*k3, Extra)
  k5 = h*fun(x + c[5]*h, y + a5[0]*k0 + a5[1]*k1 + a5[2]*k2 + a5[3]*k3 + a5[4]*k4, Extra)
    
  # ---- Step forward (Press et. al. 2007, eq. 17.2.5)
  ynew = y + b[0]*k0 + b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4 + b[5]*k5
  
  # ---- Estimate error from embedded fourth order (Press et. al. 2007, eq. 17.2.6)
  dy   = (b[0] - bemb[0])*k0 + (b[1] - bemb[1])*k1 + (b[2] - bemb[2])*k2 + (b[3] - bemb[3])*k3 + (b[4] - bemb[4])*k4 + (b[5] - bemb[5])*k5
  
  return ynew, dy


def rkdrive(fun,a,b,y0,acc,eps,h,Extra=None):
  xlist = [a]; ylist = [y0]                        #Write output into lists
  
  while xlist[-1] < b:                             #While x is still less than b
    if xlist[-1] + h > b: h = b - xlist[-1]        #Last step has to land on b
    ynew, dy = rk45step(fun,xlist[-1],ylist[-1],h,Extra=Extra) #Stepper routine
    
    if isinstance(ynew,float):                     #Special case if ynew is float (non-coupled equations)
      err = abs(dy/(acc + abs(ynew)*eps))          #Calculate error for float case
    else:
      err = 0                                      #Calculate error for list case
      for i in range(len(ynew)): err += (dy[i]/(acc + abs(ynew[i])*eps))**2
      err = sqrt(err/len(ynew))
    
    if err <= 1: xlist.append(xlist[-1]+h); ylist.append(ynew) #Check if step is accepted
    
    hnew = 0.95 * h * err**(-0.2)                  #New step size
    if hnew/h > 10: h*=10                          #Limit step jumps
    elif h/hnew > 5:  h/=5
    else: h = hnew
    
  return xlist, ylist