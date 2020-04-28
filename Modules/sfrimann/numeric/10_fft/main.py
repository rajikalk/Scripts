#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from random import gauss
from fft import fft
import pylab as mp
import math as m


#Example 1: Simple harmonic oscillator

bit = 10
f   = 0.1
dx  = 0.01
A   = 5

x = np.arange(2**bit,dtype=np.float64)*dx
y = A*np.sin(2*m.pi*f*x)

ft, freq = fft(y,dx=dx)

mp.figure(1)
mp.plot(freq+10,ft.real,'r')
mp.plot(freq-10,ft.imag,'b')
mp.plot(freq,abs(ft),'k')
mp.legend( ('real part','Imaginary part','Power Spectrum'), loc='upper right')
mp.suptitle('FFT of Simple Harmonic Oscillator')
mp.title('(red and blue curve shifted along freq for clarity)')
mp.xlabel('frequency')
mp.ylabel('Fourier Transform')
mp.savefig('FourierTransform.png')


#Example 2: Noise filter

bit = 15
f   = [0.2,1,5,10]
dx  = 0.0001
A   = [3,0.2,1,1]
ph  = [0.2,1, -0.5, 2]


x = np.arange(2**bit,dtype=np.float64)*dx
yn = np.array([gauss(0,1) for i in range(len(x))])
yo = np.zeros(len(x))
for i in range(len(ph)):
  yo += A[i]*np.sin(2*m.pi*f[i]*x + ph[i])

y = yn + yo

ft, freq = fft(y,dx=dx)

#Filtering
ftfil = ft.copy()
ftfil[np.where((abs(ftfil) < 4*np.median(abs(ftfil))))] = 0

yfil = fft(ftfil,1)

#plotting
mp.figure(2)
mp.plot(freq,abs(ft))
mp.plot(freq,abs(ftfil),'r')
mp.xlim(0,15)
mp.title('Power Spectrum of Composite oscillator')
mp.legend( ('Original PS', 'Filtered PS'), loc='upper right')
mp.xlabel('Fequency')
mp.ylabel('Power')
mp.savefig('CompositePower.png')

mp.figure(3)
mp.plot(x,y)
mp.plot(x,yfil.real+5,'r')
mp.plot(x,yo+8,'g')
mp.legend( ('Noisy Signal','Original Signal (offset)', 'Filtered Signal (offset)'), loc='upper right')
mp.title('Signal')
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('Signal.png')

