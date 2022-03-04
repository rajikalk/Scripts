#!/usr/bin/python
# -*- coding: utf-8 -*-

from StringIO import StringIO
from datetime import datetime
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

#generator for sunspot group number
def generator(file):
    lines = file.readlines()
    total = len(lines) / 39
    for n in range(total):
      lns = lines[n*39:(n+1)*39]
      io = StringIO(''.join(lns[6:-2]))
      year = int(lns[1].split()[-1])
      array = np.loadtxt(io)
      yield year, array

#read sunspot group number
def read_sgn():
  
  day   = range(1,32)
  month = range(1,13)
  
  obs = list([])
  time = list([])
  
  for year, data in generator(open('dailyrg.dat')):
    #number_of_days = float(datetime(year,12,31).timetuple().tm_yday)
    for m in month:
      temp = list([])
      for d in day:
        if not data[d-1][m] == -99:
          temp.append(data[d-1][m])
          #obs.append(data[d-1][m])
          #time.append(year + datetime(year,m,d).timetuple().tm_yday/number_of_days)
      if len(temp) == 1:
        obs.append(temp[0])
        time.append(year + m/12.)
      elif len(temp) > 1:
        obs.append(sum(temp)/float(len(temp)))
        time.append(year + m/12.)
  
  return time,obs

#read relative sunspot number
def read_rsn():
  file = open('MONTHLY.PLT')
  lines = file.readlines()
  io = StringIO(''.join(lines))
  array = np.loadtxt(io)
  time = list([])
  obs = list([])
  for i in array:
    time.append(i[0]+i[1]/12.)
    obs.append(i[2])
  return time, obs

if __name__ == '__main__':
  time,obs = read_rsn()
  
  w, h = mp.figaspect(0.5)
  fig = mp.figure(figsize=(w,h),facecolor='w') 
  #mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  mp.plot(time,obs,'.')
  #mp.plot(timegrid,fm,'k')
  mp.xlabel('year')
  mp.ylabel('Sunspot Number')
  #mp.axis([1950,2050,0,300])
  mp.savefig('sunspot.eps')
  mp.clf()
  
  #Model 1
  covar      = cv.period
  time       = np.array(time)
  obs        = np.array(obs)
  mask       = {'hyp':True ,'mean':False ,'sigma':True}
  timegrid   = np.arange(1950,2050,0.5)
  points     = 600
  ntrial     = 20
  
  hyphigh   = np.log([13, 400, 5   ])
  hyplow    = np.log([9 , 10 , 0.01])
  sigmahigh = 100
  sigmalow  = 1
  
  #hyp0 = np.random.uniform(hyplow[0],hyphigh[0],(1,ntrial))
  #for i in range(1,len(hyphigh)):
  #  hyp0       = np.concatenate((hyp0,np.random.uniform(hyplow[i],hyphigh[i].T,(1,ntrial))),axis=0)
  #sigma0 = np.random.uniform(sigmalow,sigmahigh,ntrial)
  #hyp0 = hyp0.T
  #
  #nlmlarr    = list([])
  #hypdictarr = list([])
  #for i in range(ntrial):
  #  hypdict0   = {'hyp':list(hyp0[i]),'mean':np.mean(obs),'sigma':sigma0[i]}
  #  hypdict, nlml = gp.hyper_optimize(hypdict0,covar,time[-points::2],obs[-points::2],mask)
  #  nlmlarr.append(nlml)
  #  hypdictarr.append(hypdict.copy())
  #nlmlarr    = np.array(nlmlarr)
  #hypdictarr = np.array(hypdictarr)
  #np.savez('sunspot1', nlmlarr=nlmlarr, hypdictarr=hypdictarr)
  
  npzfile = np.load('sunspot1.npz')
  nlmlarr = npzfile['nlmlarr']
  hypdictarr = npzfile['hypdictarr']
  
  u, indices = np.unique(np.around(nlmlarr,decimals=3), return_index=True)
  mi = np.argmin(u)
  nlml = u[mi]
  hypdict = hypdictarr[indices[mi]]
  
  print '############################ Model 1 - Periodic ############################'
  print '                          Number of trials = %3i' % ntrial
  print '             Number of unique convergences = %3i' % len(u)
  print 'Best log marginal likelihood for the model = %.5f' % -nlml
  print '                                 Best mean = %.5f' % hypdict['mean']
  print '                                Best sigma = %.5f' % hypdict['sigma']
  print '                      Best hyperparameters =', np.exp(hypdict['hyp'])

  
  #fm, fe, fc, al, nlml = gp.gp(hypdict,covar,time,obs,timegrid)
  
  #fe = np.sqrt(fe**2 + hypdict['sigma']**2)
  
  #w, h = mp.figaspect(0.5)
  #fig = mp.figure(figsize=(w,h),facecolor='w') 
  #mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  #mp.plot(time,obs,'.')
  #mp.plot(timegrid,fm,'k')
  #mp.xlabel('year')
  #mp.ylabel('Sunspot Number')
  #mp.axis([1950,2050,0,300])
  #mp.savefig('sunspot1.eps')
  #mp.clf()
  
  #---------------------------------------------------------------- Model 2 ------------------------------------
  
  def covar(hyp,x,xp):
    return cv.period(hyp[0:3],x,xp)*cv.SE([1,hyp[3:]],x,xp)
  
  #covar      = cv.perioddecay
  time       = np.array(time)
  obs        = np.array(obs)
  mask       = {'hyp':True ,'mean':False ,'sigma':True}
  timegrid   = np.arange(1950,2050,0.5)
  points     = 600
  ntrial     = 20
  
  hyphigh   = np.log([13, 400, 5   ,100])
  hyplow    = np.log([9 , 10 , 0.01,10 ])
  sigmahigh = 100
  sigmalow  = 1
  
  #hyp0 = np.random.uniform(hyplow[0],hyphigh[0],(1,ntrial))
  #for i in range(1,len(hyphigh)):
  #  hyp0       = np.concatenate((hyp0,np.random.uniform(hyplow[i],hyphigh[i].T,(1,ntrial))),axis=0)
  #sigma0 = np.random.uniform(sigmalow,sigmahigh,ntrial)
  #hyp0 = hyp0.T
  #
  #nlmlarr    = list([])
  #hypdictarr = list([])
  #for i in range(ntrial):
  #  hypdict0   = {'hyp':list(hyp0[i]),'mean':np.mean(obs),'sigma':sigma0[i]}
  #  hypdict, nlml = gp.hyper_optimize(hypdict0,covar,time[-points::2],obs[-points::2],mask)
  #  nlmlarr.append(nlml)
  #  hypdictarr.append(hypdict.copy())
  #nlmlarr    = np.array(nlmlarr)
  #hypdictarr = np.array(hypdictarr)
  #np.savez('sunspot2', nlmlarr=nlmlarr, hypdictarr=hypdictarr)
  
  npzfile = np.load('sunspot2.npz')
  nlmlarr = npzfile['nlmlarr']
  hypdictarr = npzfile['hypdictarr']
  
  u, indices = np.unique(np.around(nlmlarr,decimals=3), return_index=True)
  mi = np.argmin(u)
  nlml = u[mi]
  hypdict = hypdictarr[indices[mi]]
  
  print '############################ Model 2 - Periodic with decay ############################'
  print '                          Number of trials = %3i' % ntrial
  print '             Number of unique convergences = %3i' % len(u)
  print 'Best log marginal likelihood for the model = %.5f' % -nlml
  print '                                 Best mean = %.5f' % hypdict['mean']
  print '                                Best sigma = %.5f' % hypdict['sigma']
  print '                      Best hyperparameters =', np.exp(hypdict['hyp'])

  
  #fm, fe, fc, al, nlml = gp.gp(hypdict,covar,time,obs,timegrid)
  
  #fe = np.sqrt(fe**2 + hypdict['sigma']**2)
  
  #w, h = mp.figaspect(0.5)
  #fig = mp.figure(figsize=(w,h),facecolor='w') 
  #mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  #mp.plot(time,obs,'.')
  #mp.plot(timegrid,fm,'k')
  #mp.xlabel('year')
  #mp.ylabel('Sunspot Number')
  #mp.axis([1950,2050,0,300])
  #mp.savefig('sunspot2.eps')
  #mp.clf()
  
    #---------------------------------------------------------------- Model 3 ------------------------------------
  
  def covar(hyp,x,xp):
    return cv.period(hyp[0:3],x,xp)*cv.period([hyp[3],1,hyp[4]],x,xp)*cv.SE([1,hyp[5:]],x,xp)
  
  time       = np.array(time)
  obs        = np.array(obs)
  mask       = {'hyp':True ,'mean':False ,'sigma':True}
  timegrid   = np.arange(time[0],2050,1)
  points     = 600
  ntrial     = 20
  
  hyphigh   = np.log([13, 400, 5   ,500,10  ,100])
  hyplow    = np.log([9 , 10 , 0.01,50 ,0.01,10 ])
  sigmahigh = 100
  sigmalow  = 1
  
  hyp0 = np.random.uniform(hyplow[0],hyphigh[0],(1,ntrial))
  for i in range(1,len(hyphigh)):
    hyp0       = np.concatenate((hyp0,np.random.uniform(hyplow[i],hyphigh[i].T,(1,ntrial))),axis=0)
  sigma0 = np.random.uniform(sigmalow,sigmahigh,ntrial)
  hyp0 = hyp0.T
  
  nlmlarr    = list([])
  hypdictarr = list([])
  for i in range(ntrial):
    hypdict0   = {'hyp':list(hyp0[i]),'mean':np.mean(obs),'sigma':sigma0[i]}
    hypdict, nlml = gp.hyper_optimize(hypdict0,covar,time[0::8],obs[0::8],mask)
    nlmlarr.append(nlml)
    hypdictarr.append(hypdict.copy())
  nlmlarr    = np.array(nlmlarr)
  hypdictarr = np.array(hypdictarr)
  np.savez('sunspot3', nlmlarr=nlmlarr, hypdictarr=hypdictarr)
  
  #npzfile = np.load('sunspot3.npz')
  #nlmlarr = npzfile['nlmlarr']
  #hypdictarr = npzfile['hypdictarr']
  
  u, indices = np.unique(np.around(nlmlarr,decimals=3), return_index=True)
  mi = np.argmin(u)
  nlml = u[mi]
  hypdict = hypdictarr[indices[mi]]
  
  print '############################ Model 3 - Double Periodic with decay ############################'
  print '                          Number of trials = %3i' % ntrial
  print '             Number of unique convergences = %3i' % len(u)
  print 'Best log marginal likelihood for the model = %.5f' % -nlml
  print '                                 Best mean = %.5f' % hypdict['mean']
  print '                                Best sigma = %.5f' % hypdict['sigma']
  print '                      Best hyperparameters =', np.exp(hypdict['hyp'])

  
  fm, fe, fc, al, nlml = gp.gp(hypdict,covar,time,obs,timegrid)
  
  #fe = np.sqrt(fe**2 + hypdict['sigma']**2)
  
  w, h = mp.figaspect(0.5)
  fig = mp.figure(figsize=(w,h),facecolor='w') 
  mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  mp.plot(time,obs,'.')
  mp.plot(timegrid,fm,'k')
  mp.xlabel('year')
  mp.ylabel('Sunspot Number')
  mp.axis([time[0],2050,0,300])
  mp.savefig('sunspot3.eps')
  mp.clf()