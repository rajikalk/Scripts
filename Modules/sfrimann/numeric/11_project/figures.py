#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mp
import numpy as np
import scipy as sp
import random as rnd
from math import pi
from matplotlib.backends.backend_pdf import PdfPages
import covar as cv
import gp

#figure 1

def func(x,sigma,seed=None):
  np.random.seed(seed=seed)
  return (-(1-x)**2 + x + 1.1)/2 + sigma*np.random.randn(len(x))

covar      = cv.SE    #Covariance matrix to be used
x_training = np.array([0.1, 0.31, 0.35, 0.5, 0.9, 1.1])
x_grid     = np.arange(0,1.55,0.05)
x_test     = 1.4
sigma      = 0.1
y          = func(x_training,sigma,seed=23)
f          = func(x_grid,0)
hypdict    = {'hyp':[0,0],'mean':0,'sigma':sigma}
f_test,f_test_sigma,fcovar,alpha,nlml = gp.gp(hypdict,cv.SE,x_training,y,x_test)
f_grid,f_grid_sigma,fcovar,alpha,nlml = gp.gp(hypdict,cv.SE,x_training,y,x_grid)


mp.errorbar(x_training,y,fmt='o',yerr=sigma)
mp.annotate('?',(1.4,1.0))
mp.axis([0, 1.5, 0, 1.5])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig1.eps')
mp.clf()

#figure 2

mp.fill_between(x_grid,f_grid+2*f_grid_sigma,y2=f_grid-2*f_grid_sigma,color='0.75')
mp.errorbar(x_training,y,fmt='.',yerr=sigma)
mp.errorbar(x_test,f_test,fmt='.',yerr=f_test_sigma)
#mp.plot(x_grid,f_grid,'--k')
mp.plot(x_grid,func(x_grid,0),'--k')
mp.axis([0, 1.5, 0, 1.5])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig2.eps')
mp.clf()

#figure 3

x_grid     = np.arange(0,5.0,0.05)
x_training = np.array([1.0, 4.5])
y          = np.array([1.2, 0.4])
hypdict['sigma'] = 0
fm, fe, fc, alpha, nlml = gp.gp(hypdict,covar,x_training,y,x_grid)

prior_sample       = gp.draw_sample(hypdict,covar,x_grid,ns=4,seed=567)
conditioned_sample = gp.draw_sample(hypdict,covar,x_training,y,x_grid,ns=4,seed=567)

mp.fill_between(x_grid,2,y2=-2,color='0.75')
mp.plot(x_grid,prior_sample.T)
mp.axis([0, 5.0, -2.1, 2.1])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig3.eps')
mp.clf()

#figure 4

mp.fill_between(x_grid,fm+2*fe,y2=fm-2*fe,color='0.75')
mp.plot(x_grid,fm,'--k')
mp.plot(x_grid,conditioned_sample.T)
mp.plot(x_training,y,'ok')
mp.axis([0, 5.0, -1.6, 2.4])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig4.eps')
mp.clf()

#figure 5

covar      = cv.SE    #Covariance matrix to be used
x_training = np.array([0.1, 0.31, 0.35, 0.5, 0.9, 1.1])
x_grid     = np.arange(-10,10.05,0.05)
x_test     = 1.4
sigma      = 0.1
y          = func(x_training,sigma,seed=23)
f          = func(x_grid,0)
hypdict    = {'hyp':[0,0],'mean':0,'sigma':sigma}
f_test,f_test_sigma,fcovar,alpha,nlml = gp.gp(hypdict,cv.SE,x_training,y,x_test)
f_grid,f_grid_sigma,fcovar,alpha,nlml = gp.gp(hypdict,cv.SE,x_training,y,x_grid)

w, h = mp.figaspect(0.5)
fig = mp.figure(figsize=(w,h),facecolor='w') 
mp.fill_between(x_grid,f_grid+2*f_grid_sigma,y2=f_grid-2*f_grid_sigma,color='0.75')
mp.errorbar(x_training,y,fmt='.',yerr=sigma)
mp.errorbar(x_test,f_test,fmt='.',yerr=f_test_sigma)
#mp.plot(x_grid,f_grid,'--k')
mp.plot(x_grid,func(x_grid,0),'--k')
mp.axis([-5, 5, -2.1, 2.1])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig5.eps')
mp.clf()

#figure 6

covar      = cv.SE    #Covariance matrix to be used
xt         = np.random.rand(28)*14-7
xgrid      = np.arange(-7,7.01,0.01)
sigma      = 0.1
hypdict    = {'hyp':[0,0],'mean':0,'sigma':sigma}
y          = gp.draw_sample(hypdict,covar,xt,seed=89756)
fm,fe,fc,alpha,nlml = gp.gp(hypdict,cv.SE,xt,y,xgrid)

print 'Negative log marginal likelihood'
print nlml

w, h = mp.figaspect(1)
fig = mp.figure(figsize=(w,h),facecolor='w') 
mp.fill_between(xgrid,fm+2*fe,y2=fm-2*fe,color='0.75')
mp.plot(xt,y,'+')
mp.plot(xgrid,fm,'--k')
mp.axis([-7, 7, -4, 4])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig6.eps')
mp.clf()

#figure 7
hypdict    = {'hyp':list(np.log([0.5,0.3])),'mean':0,'sigma':0.0005}
fm,fe,fc,alpha,nlml = gp.gp(hypdict,cv.SE,xt,y,xgrid)

print 'Negative log marginal likelihood'
print nlml

w, h = mp.figaspect(1)
fig = mp.figure(figsize=(w,h),facecolor='w') 
mp.fill_between(xgrid,fm+2*fe,y2=fm-2*fe,color='0.75')
mp.plot(xt,y,'+')
mp.plot(xgrid,fm,'--k')
mp.axis([-7, 7, -4, 4])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig7.eps')
mp.clf()

#figure 8
hypdict    = {'hyp':list(np.log([1,3.0])),'mean':0,'sigma':0.9}
fm,fe,fc,alpha,nlml = gp.gp(hypdict,cv.SE,xt,y,xgrid)

print 'Negative log marginal likelihood'
print nlml

w, h = mp.figaspect(1)
fig = mp.figure(figsize=(w,h),facecolor='w') 
mp.fill_between(xgrid,fm+2*fe,y2=fm-2*fe,color='0.75')
mp.plot(xt,y,'+')
mp.plot(xgrid,fm,'--k')
mp.axis([-7, 7, -4, 4])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig8.eps')
mp.clf()

#figure 9

def func(x,sigma,seed=None):
  np.random.seed(seed=seed)
  return np.exp(x/10.) + np.sin(2*pi*x) + 2 + sigma*np.random.randn(len(x))

def covar(hyp,x,xp): #covar is this time a sum of two other covariance functions
  return cv.SE(hyp[:2],x,xp) + cv.period(hyp[2:],x,xp)

x_training = np.arange(0,20,0.1)
sigma      = 0.2
y          = func(x_training,sigma,seed=44)
x_grid     = np.arange(0,30,0.05)
f          = func(x_grid,0)
hypdict    = {'hyp':[0,0,0,0,0],'mean':sum(y)/float(len(y)),'sigma':sigma}
mask       = {'hyp':True,       'mean':False,               'sigma':True}

hypdictop, nlml = gp.hyper_optimize(hypdict,covar,x_training,y,mask)

print 'List of optimized hyperparamters:'
print np.exp(hypdictop['hyp'])
print 'mean:'
print np.exp(hypdictop['mean'])
print 'Optimized sigma:'
print np.exp(hypdictop['sigma'])

fm,fe,fc,al,nlml = gp.gp(hypdictop,covar,x_training,y,x_grid)

gd = np.where(x_grid >= 20)
fm = fm[gd]
fe = fe[gd]

mp.fill_between(x_grid[gd],fm+2*fe,y2=fm-2*fe,color='0.75')
mp.plot(x_training,y)
mp.axis([0, 30, 0, 15])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig9.eps')
mp.clf()

#figure 10

SEpred  = np.dot(cv.cvm(hypdictop['hyp'][:2],cv.SE,x_grid,x_training),al) + hypdictop['mean']
perpred = np.dot(cv.cvm(hypdictop['hyp'][2:],cv.period,x_grid,x_training),al)

mp.plot(x_grid,SEpred,'b')
mp.plot(x_grid,perpred,'r')
mp.axis([0, 30, -2, 15])
mp.xlabel('x')
mp.ylabel('y')
mp.savefig('fig10.eps')
mp.clf()
