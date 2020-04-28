#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import numpy as np
import os

from scipy.interpolate import griddata
from glob import glob

def _read_BM96(filename='./pms_grids/prems/canonical/tabloid020c'):
  """
  Function for reading Bernasconi & Maeder (1996) pms tracks; part of the Geneva database
  of stellar models.
  
  Reference:
    Bernasconi & Maeder (1996) A&A 307, 829
  
  Files downloaded from:
    
  
  The Default model read is canonical (no accretion luminosity) for solar abundance
  
  The columns of the returnd grid are
  
  [0]  Mass    : Stellar mass in Msun
  [1]  NB      : Number of selected point
  [2]  age     : Age of star in yr
  [3]  logL    : log luminosity of star
  [4]  logTe   : log effective temperature in log(K)
  [5]  R       : Toral radius in Rsun
  [6]  Lgrav   : Gravitational luminosity in Lsun
  [7]  Lnucl   : Nuclear luminosity in Lsun
  [8]  Mcore   : Mass of the convective core in fraction of total mass
  [9]  logTb   : log temperature at the base of the convective envelope in log(K)
  [10] Menv    : Mass of the convective envelope in fraction of total mass
  [11] H2      : 2H surface abundance
  [12] Li6     : 6Li surface abundance
  [13] Li7     : 7Li surface abundance
  """

  with open(os.path.join(os.path.dirname(__file__),filename),'r') as f:
    mass_list, table = [], []
    # step through the lines of the file
    for line in f:
      if 'STELLAR MODEL' in line:
        _, _, _, mass, _, _, _, _, Z = line.split()
        mass = float(mass)
      if len(line.split()) > 9:
        if 'NB' in line:
          continue
        table.append(line.rstrip('\n'))
        mass_list.append(mass)
      
  grid = np.loadtxt(StringIO('\n'.join(table)))
  grid = np.vstack((np.array(mass_list),grid.T)).T
  
  return grid


def _read_SDF2000(filename='./pms_grids/SDF2000/*.hrd'):
  """
  Function for reading Siess 2000 pms tracks by
  
    Siess, Dufour & Forestini (2000) A&A 358 593
  
  Files downloaded from:
  
    http://www.astro.ulb.ac.be/~siess/pmwiki/pmwiki.php/StellarModels/PMS
    
  The models are for solar abundance without overshoot
  
  The columns of the returned grid are
  [0]  model   : the model number
  [1]  phase   : the phase number (1 for PMS, 2 for MS, 3 for Post-MS, ...)
  [2]  Lum     : the total surface luminosity (in solar units)
  [3]  Mbol    : the corresponding bolometric absolute magnitude
  [4]  Reff    : the radius at photosphere (in solar units)
  [5]  Rstar   : the total surface radius [at optical depth = 0.005] (R, in solar units)
  [6]  Teff    : the effective temperature (in K)
  [7]  rho_eff : the volumic mass at photosphere (in cgs)
  [8]  log g   : the gravity at photosphere (in log10 of cgs)
  [9]  M       : the total mass (M, in solar units)
  [10] age     : the age (in yr)
  """
  
  filenames = sorted(glob(os.path.join(os.path.dirname(__file__),filename)))
  
  grid = []
  for filename in filenames:
    grid.append(np.loadtxt(filename))
  
  grid = np.vstack(grid)
  
  return grid

def _read_DM97(filenames=['./pms_grids/DAntona1998.txt','./pms_grids/DAntona1997.txt']):
  """
  Function for reading DAntona 1997+1998 pms tracks by 
    
    D'Antona & Mazzitelli (1997) Mem. S.A.It., 68, 807
  
  Files downloaded from: http://www.mporzio.astro.it/dantona/prems.html
  
  Due to column mismatches only four columns are returned in the grid although more values
  are followed

  The columns of the returned grid are
  [0]  M/Ms     : mass of the star in units of solar mass 
  [1]  log t    : log age of the star (in yr)
  [2]  L/Ls     : log luminosity in units of solar luminosity (value used Ls=3.839d+33)
  [3]  log Teff : log effective temperature (in K)
  """
  mass_list, table = [], []
  for filename in filenames:
    with open(os.path.join(os.path.dirname(__file__),filename),'r') as f:
      # step through the lines of the file
      for line in f:
        if line[0] == '#':
          if 'M/Msun=' in line:
            _, mass, _ = line.split()
            mass = float(mass)
        else:
          table.append(line.rstrip('\n'))
          mass_list.append(mass)
  
  grid = np.loadtxt(StringIO('\n'.join(table)),usecols=[1,2,3])
  grid = np.vstack((np.array(mass_list),grid.T)).T
  
  return grid

def _read_BHAC15(filename='./pms_grids/BHAC15_tracks+structure.txt'):
  """
  Function for reading BHAC15 pms grid by Baraffe et al. (2013) A&A 577 A42
  
  The file that is read here can be downloaded from
  
  http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks+structure

  The columns of the returned grid are
  [0]  M/Ms    : mass of the star in units of solar mass 
  [1]  log t   : log age of the star (in yr)
  [2]  Teff    : effective temperature (in K)
  [3]  L/Ls    : log luminosity in units of solar luminosity (value used Ls=3.839d+33)
  [4]  g       : log g  (surface gravity)
  [5]  R/Rs    : radius of the star in units of solar radius    (value used Rs=6.96d10)
  [6]  Li/Li0  : ratio of surface lithium abundance to initial abundance
  [7]  log Tc  : log of central temprature
  [8]  log ROc : log of central density (in gr/cc)
  [9]  Mrad    : mass of radiative core (in solar mass)
  [10] Rrad    : radius of radiatif core (in solar radius)
  [11] k2conv  : convective gyration radius 
  [12] k2rad   : radiative gyration radius
  """
  
  # open file
  table = []
  with open(os.path.join(os.path.dirname(__file__),filename),'r') as f:
    i = 0 # initialize counter
    # step through the lines of the file
    for line in f:
      # there's always a table after three lines beginning with !, so count them
      if line[0] == '!':
        i += 1
        continue
      # blank lines come right after tables so reset counter
      if line == '\n':
        i = 0
        continue
      # if i == 3 it means we are inside a table so save the line
      if i == 3:
        table.append(line)
  
  data = np.loadtxt(StringIO(''.join(table)))
  
  return data

def pms_interpolate(xval,yval,**kw):
  """
  Description
    Function for interpolating the grid of pms models by Baraffe et al. (2013) A&A 577 A42
  Syntax:
    interp value = pms_model(xi,yi,indexin=[0,1],indexout=2,method='linear',expandgrid=False)
  Input
    xi, yi     : points at which to interpolate the data
  Keywords
    indexin    : indices of the interpolation points (legend below; default: [0,1])
    indexout   : index of the data to be interpolated (legend below; default: 2)
    method     : method of interpolation ('linear' or 'nearest'; default: 'linear')
  Returns
    zi         : interpolated data
  Comments
    For indices 0, 2, and 5 (M/Ms, Teff, and R/Rs) interpolation is done in log-space
  index legend
  """
  
  # ------- Handle keywords
  index         = kw.get('index',[0,1,2])
  method        = kw.get('method','linear')
  interpolation = kw.get('interpolation',['linear','linear','linear'])
  grid          = kw.get('grid','DM97')
  returnlimits  = kw.get('returnlimits',False)
  
  assert method in ['linear','nearest'], "only nearest and linear interpolation allowed"
  for interpol in interpolation:
    assert interpol in ['linear','log'], "Interpolation method has to be either 'linear' or 'log'"
  
  if grid == 'BM96':
    data = _read_BM96()
    for i in [3,4,9]:
      data[:,i] = 10**data[:,i] # Go from logspace to linspace    
  elif grid == 'DM97':
    data = _read_DM97()
    data[:,1:] = 10**data[:,1:] # Go from logspace to linspace
  elif grid == 'BHAC15':
    data = _read_BHAC15()
    for i in [1,3,7,8]:
      data[:,i] = 10**data[:,i] # Go from logspace to linspace
  elif grid == 'SDF2000':
    data = _read_SDF2000()
  else:
    raise ValueError("grid must be BM96, DM97, SDF2000, or BHAC15. It is {}".format(grid))
  
  xgrid, ygrid, zgrid = data[:,index[0]], data[:,index[1]], data[:,index[2]]
  
  if returnlimits:
    return xgrid.min(), xgrid.max(), ygrid.min(), ygrid.max()
  
  # make some values log, depending on the index
  xxval = np.asarray(np.log10(xval)) if interpolation[0] == 'log' else np.asarray(xval)
  yyval = np.asarray(np.log10(yval)) if interpolation[1] == 'log' else np.asarray(yval)
  
  xxgrid = np.log10(xgrid) if interpolation[0] == 'log' else np.asarray(xgrid)
  yygrid = np.log10(ygrid) if interpolation[1] == 'log' else np.asarray(ygrid)
  zzgrid = np.log10(zgrid) if interpolation[2] == 'log' else np.asarray(zgrid)
    
  # interpolate
  zval = griddata(np.vstack((xxgrid,yygrid)).T,zzgrid,(xxval,yyval),method=method)
  
  zval = 10**zval if interpolation[2] == 'log' else zval
  
  return zval

def ageMassInterpolate(mass,age,addtime=0,grid='DM97'):
  
  # indices for age, mass, and Lum/Teff
  if grid == 'BM96':
    indices = [[0,2,3],[0,2,4]]
  elif grid == 'DM97':
    indices = [[0,1,2],[0,1,3]]
  elif grid == 'BHAC15':
    indices = [[0,1,3],[0,1,2]]
  elif grid == 'SDF2000':
    indices = [[9,10,2],[9,10,6]]
  else:
    raise ValueError("grid must be DM97, BHAC15 or SDF2000. It is {}".format(grid))
  
  luminosity = pms_interpolate(mass,age,method='linear',grid=grid,index=indices[0],interpolation = ['log','log','log'])
  Teff       = pms_interpolate(mass,age,method='linear',grid=grid,index=indices[1],interpolation = ['log','log','log'])

  if np.isfinite(luminosity).all(): # no points fall outside grid
    return luminosity, Teff, np.ones(luminosity.size,dtype=np.bool)
  
  inside = ~np.isnan(luminosity)
  
  luminosity[np.isnan(luminosity)] = pms_interpolate(mass,age,method='nearest',grid=grid,index=indices[0],interpolation = ['log','log','log'])[np.isnan(luminosity)]
  Teff[np.isnan(Teff)]             = pms_interpolate(mass,age,method='nearest',grid=grid,index=indices[1],interpolation = ['log','log','log'])[np.isnan(Teff)]
  
  massmin, massmax, agemin, agemax = pms_interpolate(mass,age,grid=grid,returnlimits=True)
  
  index = np.where(age < agemin)
  luminosity[index] = luminosity[index]*(age[index]/agemin)**5
  
  return luminosity, Teff, inside