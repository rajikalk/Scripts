#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mp
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import covar as cv
import gp

#generator for sunspot group number
def generator(file):
  lines = file.readlines()
  for n in lines:
    time = list([])
    csv = open(n[:-1])
    day = int(n[-14:-12])
    tide = np.loadtxt(csv,delimiter=',',skiprows=1,usecols=(7,8))
    csv = open(n[:-1])
    lns = csv.readlines()
    for i in lns[1:]:
      h = int(i[11:13])
      m = int(i[14:16])
      time.append(day + (h*60 + m)/1440.)
    yield time, list(tide[:,1])
      

def read_july():
  tide = list([])
  time = list([])
  
  for ttime, ttide in generator(open('/home/frimann/numeric/11_project/tides/july/list')):
    time.extend(ttime)
    tide.extend(ttide)
  
  time = np.array(time)
  tide = np.array(tide)
  index = np.argsort(time)
  time = time[index]
  tide = tide[index]
  
  index = np.nonzero(np.diff(tide))
  index = list(index[0])
  i = 0
  while len(index) > i+1:
    if index[i] != index[i+1]+1:
      del index[i+1]
    i += 1
  time = time[index]
  tide = tide[index]

  return time, tide

if __name__ == '__main__':
  time, obs = read_july()
  
  w, h = mp.figaspect(0.5)
  fig = mp.figure(figsize=(w,h),facecolor='w') 
  #mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  mp.plot(time[-1000::1],obs[-1000::1],'.')
  #mp.plot(timegrid,fm,'k')
  mp.xlabel('Day in july 2011')
  mp.ylabel('Tide Height (m)')
  #mp.axis([1950,2050,0,300])
  mp.savefig('tideheight.eps')
  mp.clf()
  
  #Model 1
  covar      = cv.period
  mask       = {'hyp':True ,'mean':False ,'sigma':True}
  timegrid   = np.arange(20,30,0.01)
  points     = 500
  delta      = 2
  ntrial     = 30
  
  hyphigh   = np.log([2   ,10   , 10   ])
  hyplow    = np.log([0.5 ,0.1 , 0.01])
  sigmahigh = 1
  sigmalow  = 0.001
  
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
  #  hypdict, nlml = gp.hyper_optimize(hypdict0,covar,time[-delta*points::delta],obs[-delta*points::delta],mask)
  #  nlmlarr.append(nlml)
  #  hypdictarr.append(hypdict.copy())
  #nlmlarr    = np.array(nlmlarr)
  #hypdictarr = np.array(hypdictarr)
  #np.savez('tideheight1', nlmlarr=nlmlarr, hypdictarr=hypdictarr)
  
  npzfile = np.load('tideheight1.npz')
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

  
  fm, fe, fc, al, nlml = gp.gp(hypdict,covar,time[-delta*points:],obs[-delta*points:],timegrid)
  
  fe = np.sqrt(fe**2 + hypdict['sigma']**2)
  
  w, h = mp.figaspect(0.5)
  fig = mp.figure(figsize=(w,h),facecolor='w') 
  mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  mp.plot(time,obs,'.')
  mp.plot(timegrid,fm,'k')
  mp.xlabel('Day in july 2011')
  mp.ylabel('Tide Height (m)')
  mp.axis([20,30,0.5,5])
  mp.savefig('tideheight1.eps')
  mp.clf()
  
  #Model 2
  def covar(hyp,x,xp):
    return cv.period(hyp[0:3],x,xp)*cv.period([hyp[3],1,hyp[4]],x,xp)
  
  covar      = cv.period
  mask       = {'hyp':True ,'mean':False ,'sigma':True}
  timegrid   = np.arange(20,30,0.01)
  points     = 500
  delta      = 2
  ntrial     = 5
  
  hyphigh   = np.log([10   ,10   , 10,  10, 10  ])
  hyplow    = np.log([0.5 ,0.1 , 0.01, 0.5, 0.01])
  sigmahigh = 1
  sigmalow  = 0.001
  
  hyp0 = np.random.uniform(hyplow[0],hyphigh[0],(1,ntrial))
  for i in range(1,len(hyphigh)):
    hyp0       = np.concatenate((hyp0,np.random.uniform(hyplow[i],hyphigh[i].T,(1,ntrial))),axis=0)
  sigma0 = np.random.uniform(sigmalow,sigmahigh,ntrial)
  hyp0 = hyp0.T
  
  nlmlarr    = list([])
  hypdictarr = list([])
  for i in range(ntrial):
    hypdict0   = {'hyp':list(hyp0[i]),'mean':np.mean(obs),'sigma':sigma0[i]}
    hypdict, nlml = gp.hyper_optimize(hypdict0,covar,time[-delta*points::delta],obs[-delta*points::delta],mask)
    nlmlarr.append(nlml)
    hypdictarr.append(hypdict.copy())
  nlmlarr    = np.array(nlmlarr)
  hypdictarr = np.array(hypdictarr)
  np.savez('tideheight2', nlmlarr=nlmlarr, hypdictarr=hypdictarr)
  
  npzfile = np.load('tideheight2.npz')
  nlmlarr = npzfile['nlmlarr']
  hypdictarr = npzfile['hypdictarr']
  
  u, indices = np.unique(np.around(nlmlarr,decimals=3), return_index=True)
  mi = np.argmin(u)
  nlml = u[mi]
  hypdict = hypdictarr[indices[mi]]
  
  print '############################ Model 2 - Double Periodic ############################'
  print '                          Number of trials = %3i' % ntrial
  print '             Number of unique convergences = %3i' % len(u)
  print 'Best log marginal likelihood for the model = %.5f' % -nlml
  print '                                 Best mean = %.5f' % hypdict['mean']
  print '                                Best sigma = %.5f' % hypdict['sigma']
  print '                      Best hyperparameters =', np.exp(hypdict['hyp'])

  
  fm, fe, fc, al, nlml = gp.gp(hypdict,covar,time[-delta*points:],obs[-delta*points:],timegrid)
  print nlml
  fe = np.sqrt(fe**2 + hypdict['sigma']**2)
  
  w, h = mp.figaspect(0.5)
  fig = mp.figure(figsize=(w,h),facecolor='w') 
  mp.fill_between(timegrid,fm+2*fe,y2=fm-2*fe,color='0.75')
  mp.plot(time,obs,'.')
  mp.plot(timegrid,fm,'k')
  mp.xlabel('Day in july 2011')
  mp.ylabel('Tide Height (m)')
  mp.axis([20,30,0.5,5])
  mp.savefig('tideheight2.eps')
  mp.clf()