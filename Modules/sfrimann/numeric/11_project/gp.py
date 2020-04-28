#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import scipy.optimize as sp
from math import pi
from matplotlib.backends.backend_pdf import PdfPages
from covar import cvm

# ---- Error definitions ----
class InputError(Exception):
  pass


#---------------------------------------------------------------------------------------------------#
#main function. Calculates the Gaussian Process
def gp(*args,**kw):
  hypdict = args[0]  #hyperparameters
  
  tol = kw.get('tol',1e-4)
  #Sanity checks
  if not isinstance(hypdict,dict):
    raise InputError('hyperparameters must be given in dictionary')
  elif not (hypdict.has_key('hyp') or hypdict.has_key('mean') or hypdict.has_key('sigma')):
    raise InputError('hyperparameter dictionary must have keys hyp, mean, and sigma')
  else:
    hyp   = hypdict['hyp']
    mean  = hypdict['mean']
    sigma = hypdict['sigma']
    print hyp
  if len(args) == 3: #Three inputs: prior
    covar = args[1]  #callable covariance function
    x     = args[2]  #vector of x-values (must be one-dimensional)
    
    f_prior_mean  = np.zeros(len(x)) + mean
    if sigma != 0:
      f_prior_covar = cvm(hyp,covar,x,x) + sigma*sigma*np.eye(len(x))
    else:
      f_prior_covar = cvm(hyp,covar,x,x)
    return f_prior_mean, f_prior_covar
  elif len(args) == 4: #Four inputs: training
    covar = args[1]    #callable covariance function
    x     = args[2]    #vector of training x-values (must be one-dimensional)
    y     = args[3]    #vector of training y-values
    
    Kxx = cvm(hyp,covar,x,x)  #covariance matrix
    
    #If sigma is given add the variances to the diagonal of Kxx
    if sigma != 0 and sigma > tol:
      if isinstance(sigma,(tuple,list,np.ndarray)):
        if len(sigma) == len(x):
          if isinstance(sigma,(list,tuple)):
            for i in range(len(sigma)):
              Kxx[i,i] += sigma[i]*sigma[i]
          else:
            Kxx += np.diag(sigma*sigma)
        else:
          raise InputError('sigma must be same length as x')
      elif isinstance(sigma,(int,long,float)):
        Kxx += sigma*sigma*np.eye(len(x))
    else:
      Kxx += tol*np.eye(len(x)) #if sigma is zero add small number to diagonal to avoid numerical problems
    
    L     = np.linalg.cholesky(Kxx)  #Cholesky decomposition
    alpha = np.linalg.solve(np.transpose(L),np.linalg.solve(L,y-mean))  #alpha value
    nlml  = np.dot(y-mean,alpha)/2. + sum(np.log(np.diag(L))) + len(y)/2.*np.log(2*pi)  #negative log marginal likelihood
    return nlml
  elif len(args) == 5: #Five inputs: prediction
    covar = args[1]  #callable covariance function
    x     = args[2]  #vector of training x-values (must be one-dimensional)
    y     = args[3]  #vector of training y-values
    xs    = args[4]  #vector of prediction x-values
    
    Kxx    = cvm(hyp,covar,x,x)   #Covariance of training input`
    Kxsx   = cvm(hyp,covar,xs,x)  #Cross covariance
    Kxsxs  = cvm(hyp,covar,xs,xs) #Covariance of prediction input
    
    #If sigma is given add the variances to the diagonal of Kxx
    if sigma != 0 and sigma > tol:
      if isinstance(sigma,(tuple,list,np.ndarray)):
        if len(sigma) == len(x):
          if isinstance(sigma,(list,tuple)):
            for i in range(len(sigma)):
              Kxx[i,i] += sigma[i]*sigma[i]
          else:
            Kxx += np.diag(sigma*sigma)
        else:
          raise InputError('sigma must be same length as x')
      elif isinstance(sigma,(int,long,float)):
        Kxx += sigma*sigma*np.eye(len(x))
    else:
      Kxx += tol*np.eye(len(x)) #if sigma is zero add small number to diagonal to avoid numerical problems
    
    L     = np.linalg.cholesky(Kxx)                                    #Cholesky decomposition of Kxx
    alpha = np.linalg.solve(np.transpose(L),np.linalg.solve(L,y-mean)) #alpha vector
    fsm   = np.dot(Kxsx,alpha) + mean                                  #predictive mean of latent function
    v     = np.linalg.solve(L,np.transpose(Kxsx))                      #v vector
    fsc   = Kxsxs - np.dot(np.transpose(v),v)                          #predictive covariance of latent function
    nlml  = np.dot(y-mean,alpha)/2. + sum(np.log(np.diag(L))) + len(y)/2.*np.log(2*pi) #negative log marginal likelihood
    
    if np.size(fsc) > 1:
      fse = np.sqrt(np.diag(fsc))
    else:
      fse = np.sqrt(fsc)
    
    return fsm, fse, fsc, alpha, nlml
  else: #Else wrong number of inputs: Raise error
    raise InputError('Wrong number of inputs')
  
#---------------------------------------------------------------------------------------------------#
# Draws a sample from either a prior or conditioned distribution given appropiate input
def draw_sample(*args,**kw):
  
  #Get keywords
  ns   = kw.get('ns',1)
  seed = kw.get('seed',None)
  
  if len(args) == 3:
    fm, fc = gp(*args) #Calculate mean vector and covariance matrix of the GP
  else:
    fm, fe, fc, al, nlml = gp(*args)
  
  np.random.seed(seed=seed) #seed (default is None)
  #if ns == 1 only draw one sample. Otherwise draw ns samples
  if ns == 1:
    return np.random.multivariate_normal(fm,fc)
  else:
    return np.random.multivariate_normal(fm,fc,(ns))



#---------------------------------------------------------------------------------------------------#
# Optimizes hyperparameters using downhill simplex
def hyper_optimize(hypdict0,covar,x,y,mask):

  #Wrapper function
  def wrapper_gp(hyp0,covar,x,y,hypdict0,mask):
    if (mask['mean'] and mask['sigma']) is True:
      hypdict0['hyp']   = hyp0[2:]
      hypdict0['mean']  = hyp0[1]
      hypdict0['sigma'] = np.exp(hyp0[0])
    elif mask['mean'] is True:
      hypdict0['hyp']   = hyp0[1:]
      hypdict0['mean']  = hyp0[0]
    elif mask['sigma'] is True:
      hypdict0['hyp']   = hyp0[1:]
      hypdict0['sigma'] = np.exp(hyp0[0])
    else:
      hypdict0['hyp']   = hyp0
    return gp(hypdict0,covar,x,y)
  
  #Sanity Checks
  if not isinstance(hypdict0,dict):
    raise InputError('hyperparameters must be given in dictionary')
  elif not (hypdict0.has_key('hyp') or hypdict0.has_key('mean') or hypdict0.has_key('sigma')):
    raise InputError('hyperparameter dictionary must have keys hyp, mean, and sigma')
  
  hypdict = hypdict0.copy() #Copy input dictionary
  
  if not isinstance(mask,dict):
    raise InputError('masking must be given in dictionary')
  elif not (mask.has_key('hyp') or mask.has_key('mean') or mask.has_key('sigma')):
    raise InputError('mask dictionary must have keys hyp, mean, and sigma')
  elif (mask['mean'] and mask['sigma']) is True: #Optimize hyp, mean, and sigma
    hyp0      = [np.log(hypdict['sigma']),hypdict['mean']] + hypdict['hyp'][:] #Create hyp0 vector for scipy.fmin
    hyp,nlml,a1,a2,a3 = sp.fmin(wrapper_gp,hyp0,args=(covar,x,y,hypdict,mask),full_output=True)
  elif mask['mean'] is True: #Optimize hyp and mean
    hyp0      = hypdict['hyp'][:].insert(0,hypdict['mean'])  #Create hyp0 vector for scipy.fmin
    hyp,nlml,a1,a2,a3 = sp.fmin(wrapper_gp,hyp0,args=(covar,x,y,hypdict,mask),full_output=True)
  elif mask['sigma'] is True: #Optimize hyp and sigma
    hyp0      = hypdict['hyp'][:]  #Create hyp0 vector for scipy.fmin
    hyp0.insert(0,np.log(hypdict['sigma']))
    hyp,nlml,a1,a2,a3 = sp.fmin(wrapper_gp,hyp0,args=(covar,x,y,hypdict,mask),full_output=True)
  else: #Optimize hyp
    hyp0 = hypdict['hyp'][:]  #Create hyp0 vector for scipy.fmin
    hyp,nlml,a1,a2,a3 = sp.fmin(wrapper_gp,hyp0,args=(covar,x,y,hypdict,mask),full_output=True)
  
  #Return optimized parameters and negative log marginal likelihood for the best solution
  return hypdict, nlml
