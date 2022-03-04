#!/usr/bin/env python
# -*- coding: utf-8 -*-



import numpy as np
import os
import matplotlib.pyplot as mp
import re

from glob import glob

from mirpy import miriad

def uvfit(datadir='.',uvmin=15,gauss1d=False):
  
  curdir = os.getcwd()
  
  try:
    os.chdir(datadir)
    visdata = glob('./*.uvdata')
    if len(visdata) == 0:
      raise IOError("could not find %s" % visdata)
    elif len(visdata) > 1:
      raise IOError("Multiple visibility files found %i" % (len(visdata)))

    visdata = visdata[0]
    source, void = os.path.splitext(visdata)
    
    imout, rms, beamMajor, beamMinor, beampa = read_imfit(''.join(open('imfit.txt','r').readlines()))
    
    os.system('rm -r *.blct.uv')
    subfiles = sorted(glob('*sub*.aver.uv'))
    uvname = source+'.aver.uv' if len(subfiles) == 0 else subfiles[-1]
#     uvname = source+'.aver.uv'
    miriad.uvcat(vis=uvname, select='uvrange(%i,100)' % uvmin, out=source+'.blct.uv')
    
    if gauss1d:
      uvfit = miriad.uvfit(vis=source+'.blct.uv', fix='c', object='gaussian',spar='%f,%f,%f,%f' % (imout['peak'],imout['x'],imout['y'],np.sqrt(imout['dmaj']*imout['dmin'])))
    else:
      uvfit = miriad.uvfit(vis=source+'.blct.uv', object='gaussian',spar='%f,%f,%f,%f,%f,%f' % (imout['peak'],imout['x'],imout['y'],imout['dmaj'],imout['dmin'],imout['dpa']))
    with open('uvfit.txt','w') as f:
      f.write(uvfit)
    
    print(uvfit)
    
    uvfit = re.findall("([-+]?\d+[\.]?\d*[-+Ee]*\d*)",uvfit)
    
    if gauss1d:
      uvpar = {'flux':float(uvfit[6]),'sflux':float(uvfit[7]),
               'x':float(uvfit[8]),'y':float(uvfit[9]),
               'sx':float(uvfit[10]),'sy':float(uvfit[11]),
               'fwhmx':float(uvfit[12]),'fwhmy':float(uvfit[13]),
               'sfwhmx':float(uvfit[14]),'sfwhmy':float(uvfit[15])}
    else:
      uvpar = {'flux':float(uvfit[6]),'sflux':float(uvfit[7]),
               'x':float(uvfit[8]),'y':float(uvfit[9]),
               'sx':float(uvfit[10]),'sy':float(uvfit[11]),
               'fwhmx':float(uvfit[12]),'fwhmy':float(uvfit[13]),
               'sfwhmx':float(uvfit[14]),'sfwhmy':float(uvfit[15]),
               'pa':float(uvfit[16]),'spa':float(uvfit[17])}
    
  finally:
    os.chdir(curdir)
  
  return uvpar

def clean(v0,datadir='.'):

  curdir = os.getcwd()

  try:
    os.chdir(datadir)
    visdata = glob('./*.uvdata')
    if len(visdata) == 0:
      raise IOError("could not find %s" % visdata)
    elif len(visdata) > 1:
      raise IOError("Multiple visibility files found %i" % (len(visdata)))

    visdata = visdata[0]
    source, void = os.path.splitext(visdata)

    #set C18O rest velocity
    void = miriad.puthd(_in=visdata+'/restfreq',value=219.560368)

    dv = 0.2
    Dv = 4.
    nchan = int(Dv/dv)
#     v0 = 9.1 #km/s
    # cut out all but +/- 2 km/s from rest
    # velocity syntax is number of channels, first channel, channel spacing, 
    #         channel width
    try:
      miriad.uvaver(vis=visdata,line='velocity,%i,%.1f,%.1f,%.1f' % (nchan,v0-Dv/2,dv,dv), out=source+'.vcut.uv')
    except:
      pass

    # integrate (average) down to one channel
    # will make integrated intensity, but in wrong units
    # to get to integrated intensity, multiply by (n*dv)
    try:
      miriad.uvaver(vis=source+'.vcut.uv',line='velocity,%i,%.1f,%.1f,%.1f' % (1,v0-Dv/2,Dv,Dv), out=source+'.aver.uv')
    except:
      pass

    os.system('rm -rf *.dirtymap* *.beam*')

    output = miriad.invert(vis=source+'.aver.uv', map=source+'.dirtymap', beam=source+'.beam', cell=0.8, imsize=100, robust=1.0, options='systemp,double')

    rms = float(re.findall("Theoretical rms noise: ([-+]?\d+[\.]?\d*[-+Ee]*\d*)",output)[0])

    #clean the dirty map, just do a very basic clean for now
    #restore image
    os.system('rm -rf *.clean* *.cleanmap*')
    cleanoutput = miriad.clean(map=source+'.dirtymap', beam=source+'.beam', out=source+'.clean', cutoff=rms*2, niters=10000, gain=0.05, options='positive')

    print(cleanoutput)

    miriad.restor(map=source+'.dirtymap', beam=source+'.beam', model=source+'.clean',out=source+'.cleanmap')

    miriad.fits(_in=source+'.dirtymap', op='xyout', out=source+'.dirtymap.fits')
    miriad.fits(_in=source+'.beam', op='xyout', out=source+'.beam.fits')
    miriad.fits(_in=source+'.cleanmap', op='xyout', out=source+'.cleanmap.fits')
  finally:
    os.chdir(curdir)

def uvamp(datadir='.',uvdelta=5,factor=1.,offset=[],visdata=None):
  
  curdir = os.getcwd()
  
  try:
    
    assert os.path.exists(visdata), "Could not find file %s" % visdata
  
    source, void = os.path.splitext(visdata)

    #select the region for imfit, rename the region file 
    #fit a gaussian in the image plane, use the resulting source position 
    #     to calculate amplitude vs uvdistance
    #uvamp bin syntax is number of bins, bin spacing, and units

    if len(offset) == 2:
      offx, offy = offset
    else:
      if not os.path.exists('imfit.txt'):
        os.system('rm imfit.region')
        miriad.cgcurs(_in=source+'.cleanmap', device='/xs', options='wedge,cursor,region')
        os.system('mv cgcurs.region imfit.region')
        output = miriad.imfit(_in=source+'.cleanmap', object='gaussian', region='@imfit.region')
        with open('imfit.txt','w') as f:
          f.write(output)
        print(output)
      else:
        output = open('imfit.txt','r').readlines()
        output = ''.join(output)

      offx, offy = list(map(float,re.findall(r"Offset Position \(arcsec\):\s*([-+]?\d+[\.]?\d*[-+Ee]*\d*)\s*([-+]?\d+[\.]?\d*[-+Ee]*\d*)",output)[0]))      
    
    uvamp = miriad.uvamp(vis=visdata,offset='%f,%f' % (offx,offy), bin='200,%i,klam' % uvdelta).decode('ascii')
    with open('uvamp.txt','w') as f:
      f.write(uvamp)
    
    def load_uvamp(filename):
      uvmin, uvmax, amplitude, sigma, sn, expect, pnts = np.genfromtxt(filename,skip_header=7,skip_footer=3,unpack=True)
  
      uvrad = 0.5*(uvmin+uvmax)
      index = pnts > 0
      uvrad, amplitude, sigma, sn, expect = uvrad[index], amplitude[index], sigma[index], sn[index], expect[index]
      amplitude = amplitude*factor
      sigma = sigma*factor
      expect = expect*factor
  
      return uvrad, amplitude, sigma, sn, expect
    
    uvrad, amplitude, sigma, sn, expect = load_uvamp('uvamp.txt')
    
  finally:
    os.chdir(curdir)
  
  return uvrad, amplitude, sigma, sn, expect

def read_uvfit(input):
  """
  read output string from miriad uvfit routine and extract results
  """
  
  strings = input.split('Source')
  rms     = re.findall(r'RMS residual is (.*)',strings[0])
  assert len(rms) == 1, "String longer than 1. %i" % len(rms)
  rms     = float(rms[0])
  
  strings = strings[1:]
  result = []
  for string in strings:
    object = re.findall(r'Object type: ([a-z]+)',string)[0]
    if object == 'gaussian':
      # get flux
      try:
        flux, sflux = re.findall(r'Flux:\s*(\S*)\s*\+\/- (\S*)',string)[0]
        flux, sflux = float(flux), float(sflux)
      except:
        flux = re.findall(r'Flux:\s*(\S*)',string)[0]
        flux, sflux = float(flux), np.nan
      # get offset positions
      x, y = re.findall(r'Offset Position \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
      x, y = float(x), float(y)
      # get offset position errors
      try:
        sx, sy = re.findall(r'Positional errors \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
        sx, sy = float(sx), float(sy)
      except:
        sx, sy = np.nan, np.nan
      # get major and minor axes
      maj, min = re.findall(r'Major,minor axes \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
      maj, min = float(maj), float(min)
      # get axes errors
      try:
        smaj, smin = re.findall(r'Axes errors \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
        smaj, smin = float(smaj), float(smin)
      except:
        smaj, smin = np.nan, np.nan
      # get position angle
      pos = re.findall(r'Position angle \(degrees\):\s*(\S*)',string)[0]
      pos = float(pos)
      # get position angle error
      try:
        spos = re.findall(r'Pos\. angle error \(degrees\):\s*(\S*)',string)[0]
        spos = float(spos)
      except:
        spos = np.nan    
      
      result.append(dict(flux=flux,sflux=sflux,x=x,y=y,sx=sx,sy=sy,maj=maj,
                         min=min,smaj=smaj,smin=smin,pos=pos,spos=spos,object=object))
    elif object == 'point':
      # get flux
      flux, sflux = re.findall(r'Flux:\s*(\S*)\s*\+\/- (\S*)',string)[0]
      flux, sflux = float(flux), float(sflux)
      # get offset positions
      x, y = re.findall(r'Offset Position \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
      x, y = float(x), float(y)
      # get offset position errors
      try:
        sx, sy = re.findall(r'Positional errors \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
        sx, sy = float(sx), float(sy)
      except:
        sx, sy = np.nan, np.nan
      
      result.append(dict(flux=flux,sflux=sflux,x=x,y=y,sx=sx,sy=sy,object=object))
      
    else:
      raise ValueError('Object name %s does not match any known' % object)
  
  return result, rms

def read_imfit(input,filename=False):
  """
  read output string from miriad imfit routine and extract results
  """
  
  if filename is True:
    input = open(input,'r').readlines()
    input = ''.join(input)
  
  strings = input.split('Source')
  rms     = re.findall(r'RMS residual is (\S*)',strings[0])
  assert len(rms) == 1, "String longer than 1. %i" % len(rms)
  rms     = float(rms[0])
  
  # beam major minor axes
  try:
    beamMaj, beamMin = re.findall(r'Beam Major, minor axes \(arcsec\):\s*(\S*)\s*(\S*)',strings[0])[0]
    beamMaj, beamMin = float(beamMaj), float(beamMin)
  
    # beam position angle
    beampa = re.findall(r'Beam Position angle \(degrees\):\s*(\S*)',strings[0])[0]
    beampa = float(beampa)
  except:
    beamMaj, beamMin, beampa = np.nan, np.nan, np.nan
  
  strings = strings[1:]
  result = []
  for string in strings:
    object = re.findall(r'Object type: ([a-z]+)',string)[0]
    if object == 'beam':
      try:
        # get major axis
        maj, smaj = re.findall(r'Major axis \(arcsec\):\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]
        maj, smaj = float(maj), float(smaj)
      except:
        maj       = re.findall(r'Major axis \(arcsec\):\s*(\S*)',string)[0]
        maj, smaj = float(maj), np.nan
      
      try:
        # get minor axis
        min, smin = re.findall(r'Minor axis \(arcsec\):\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]
        min, smin = float(min), float(smin)
      except:
        # get minor axis
        min       = re.findall(r'Minor axis \(arcsec\):\s*(\S*)',string)[0]
        min, smin = float(min), np.nan
      
      # get position angle
      pa , spa  = re.findall(r'Position angle \(degrees\):\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]  
      pa , spa  = float(pa), float(spa)
      
      result.append(dict(maj=maj,smaj=smaj,min=min,smin=smin,pa=pa,spa=spa))
    elif object == 'gaussian':
      # get peak
      peak, speak = re.findall(r'Peak value:\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]
      peak, speak = float(peak), float(speak)
      
      # get integrated flux
      flux = re.findall(r'Total integrated flux:\s*(\S*)',string)[0]
      flux = float(flux)
      
      # get offset positions
      x, y = re.findall(r'Offset Position \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
      x, y = float(x), float(y)
      # get offset position errors
      try:
        sx, sy = re.findall(r'Positional errors \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
        sx, sy = float(sx), float(sy)
      except:
        sx, sy = np.nan, np.nan
      
      # get major axis
      maj, smaj = re.findall(r'Major axis \(arcsec\):\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]
      maj, smaj = float(maj), float(smaj)
      
      # get minor axis
      min, smin = re.findall(r'Minor axis \(arcsec\):\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]
      min, smin = float(min), float(smin)
      
      # get position angle
      pa , spa  = re.findall(r'Position angle \(degrees\):\s*(\S*)\s*\+\/-\s*(\S*)',string)[0]  
      pa , spa  = float(pa), float(spa)
      
      # deconvolved major, minor axes
      dmaj, dmin = re.findall(r'Deconvolved Major, minor axes \(arcsec\):\s*(\S*)\s*(\S*)',string)[0]
      dmaj, dmin = float(dmaj), float(dmin)
      
      # deconvolved position angle
      dpa = re.findall(r'Deconvolved Position angle \(degrees\):\s*(\S*)',string)[0]
      dpa = float(dpa)
      
      result.append(dict(peak=peak,speak=speak,flux=flux,x=x,y=y,sx=sx,sy=sy,maj=maj,
                         smaj=smaj,min=min,smin=smin,pa=pa,spa=spa,dmaj=dmaj,dmin=dmin,
                         dpa=dpa))
    else:
      raise ValueError('Object name %s does not match any known' % object)
  
  return result, rms, beamMaj, beamMin, beampa

def subtract(datadir='.',time=0):
  
  curdir = os.getcwd()
  
  try:
    os.chdir(datadir)
    
    visdata = glob('*.uvdata')
    if len(visdata) == 0:
      raise IOError("could not find %s" % visdata)
    elif len(visdata) > 1:
      raise IOError("Multiple visibility files found %i" % (len(visdata)))

    visdata = visdata[0]
    source, void = os.path.splitext(visdata)
    
    cleanext = '.cleanmap' if time == 0 else '.sub%i.cleanmap' % (time-1)
    
    if not os.path.exists('imfit_sub%i.txt' % time):
      os.system('rm sub%i.region' % time)
      miriad.cgcurs(_in=source+cleanext, device='/xs', options='wedge,cursor,region')
      os.system('mv cgcurs.region sub%i.region' % time)
    
      output = miriad.imfit(_in=source+cleanext, object='gaussian', region='@sub%i.region' % time)
    
      with open('imfit_sub%i.txt' % time,'w') as f:
        f.write(output)
    else:
      output = open('imfit_sub%i.txt' % time,'r').readlines()
      output = ''.join(output)
    
    print(output)
    
    imout, rms, beamMajor, beamMinor, beampa = read_imfit(output)
    
    # create zero image
    if not os.path.exists(source+'.zero'):
      miriad.maths(exp=source+cleanext+'*0', out=source+'.zero')
    
    os.system('rm -r *sub%i.model* *sub%i.dirtymap* *sub%i.beam *sub%i.aver.uv' % (time,time,time,time))
    # make model gaussian
    miriad.imgen(_in=source+'.zero',out=source+'.sub%i.model' % time, object='gaussian',spar='%f,%f,%f,%f,%f,%f' % (imout['peak'],imout['x'],imout['y'],imout['dmaj'],imout['dmin'],imout['dpa']))
    
    #convert to Jy/pixel
    factor = 0.8**2/(2*np.pi*beamMajor*beamMinor/2.35482**2)
    print(1/factor)
    miriad.maths(exp=source+'.sub%i.model*%f' % (time,factor), out=source+'.sub%i.model.perpixel' % time)
    
    # convert model gaussian to uv space
    uvext = '.aver.uv' if time == 0 else '.sub%i.aver.uv' % (time-1)
    miriad.uvmodel(vis=source+uvext,model=source+'.sub%i.model.perpixel' % time,options='subtract',out=source+'.sub%i.aver.uv' % time)
    
    # create new dirty image with subtracted uv data
    miriad.invert(vis=source+'.sub%i.aver.uv' % time,map=source+'.sub%i.dirtymap' % time, beam=source+'.sub%i.beam' % time, cell=0.8, imsize=100, robust=1.0, options='systemp,double')

    miriad.fits(_in=source+'.sub%i.dirtymap' % time, op='xyout', out=source+'.sub%i.dirtymap.fits' % time)
    
  finally:
    os.chdir(curdir)