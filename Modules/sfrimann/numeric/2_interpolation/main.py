#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as mp
from interpolation import *

#Gaussian interpolation
xi   = np.arange(-12,14,2)
yi   = exp(-(xi**2.)/9.)
xtab = arange(-12,12.01,0.01)

ypol  = polyinterpol(xi,yi,xtab)
ylin  = lininterpol(xi,yi,xtab)
yquad = quadinterpol(xi,yi,xtab)

mp.figure(1)
mp.plot(xi,yi,'o')
mp.plot(xtab,ypol,'r')
mp.plot(xtab,ylin,'g')
mp.plot(xtab,yquad,'k')
mp.axis([-15,15,-0.2,1.2])
mp.legend( ('Data','Poly interp','Lin interp','Quad interp'), loc='upper right')
mp.title('Interpolation of Gaussian')
mp.savefig('Gaussian_interp.png', format='png')

#Step function
yi = np.concatenate((np.zeros(7),np.ones(6)))

ypol  = polyinterpol(xi,yi,xtab)
ylin  = lininterpol(xi,yi,xtab)
yquad = quadinterpol(xi,yi,xtab)

mp.figure(2)
mp.plot(xi,yi,'o')
mp.plot(xtab,ypol,'r')
mp.plot(xtab,ylin,'g')
mp.plot(xtab,yquad,'k')
mp.axis([-15,15,-0.2,1.2])
mp.legend( ('Data','Poly interp','Lin interp','Quad interp'), loc='upper left')
mp.title('Interpolation of Step function')
mp.savefig('Step_interp.png', format='png')