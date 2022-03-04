#!/usr/bin/env python
# encoding: utf-8
"""
clfind2d.py

Created by José Ramón Sánchez-Gallego on 2010-06-07.

This script is a translation of the clfind2d.pro code by J. P. Williams
http://www.ifa.hawaii.edu/users/jpw/clumpfind.shtml.

Needs: pyfits, numpy and the routine search2d provided in a separated file
Usage: clfind2d.py fileIn levels [-l|--nolog] [-o|--out]
        fileIn: The original FITS file
        levels: The levels of the contours [l1,l2,l3,...]
        -l|--nolog: Disable the logging system
        -o|--out: The root of the output files. If not present, the [fileIn] root
                 will be used
        -n|-npixmin: Minimum number of pixels for a region to be accepted
        -i|--interval: Interval of levels [start,end,step]. Disables the levels argument

"""


import numpy as np
import sys
import os
from time import strftime, gmtime, clock
from scipy.ndimage import generate_binary_structure, label

try:
  from astropy.io import fits as pf
except:
  import pyfits as pf

def printLog(text, log, silent):
  if silent is False:
    sys.stdout.write(text)
    sys.stdout.flush()
  if log is not None:
    unit = open(log, 'a')
    print >>unit, text,
    unit.close()
  return

# This routine identifies all the clumps between two consecutive contour levels.
def defReg(data, levMin, levMax, diagonal=False):

# Pixels with values between the two levels
  levelMask = (data >= levMin) & (data < levMax)
  nPixInLevel = np.sum(levelMask)

# Is no pixels between [levMin,levMax) are found, returns
  if nPixInLevel == 0:
    return nPixInLevel, None, 0

  if diagonal:
    struct = generate_binary_structure(2, 2)
  else:
    struct = generate_binary_structure(2, 1)
  reg, nReg = label(levelMask, structure=struct)
  return nPixInLevel, reg, nReg

# Decides if a set of pixels is a new clump or belong to a previous clump.
# In that case, merges the two clumps or mask the pixels to the closest
# clump
def defClump(data, levMin, reg, nReg, mask, clumpPeak, nCl):

  nNew = 0
  if nReg == 0: return 0
  
  nPix2, pix2, nReg2 = defReg(data, levMin, data.max()+1000, diagonal=True)

  for nR in range(1, nReg + 1):
    pix1 = np.where(reg == nR)
    nPix1 = pix1[0].size

    peak = np.argmax(data[pix1])
    iPeak = pix1[0][peak]
    jPeak = pix1[1][peak]

    maskPix2 = pix2 == pix2[iPeak, jPeak]
    aPix2 = mask[maskPix2]

    if np.amax(aPix2) == 0:
      # Found new clump!
      nNew += 1
      mask[maskPix2] = nCl + nNew
      peak = np.argmax(data[pix1])
      iPeak = pix1[0][peak]
      jPeak = pix1[1][peak]
      clumpPeak[nCl+nNew, 0] = iPeak
      clumpPeak[nCl+nNew, 1] = jPeak
    else:
      nC = np.unique(np.sort(aPix2, axis=None))
      if np.amin(nC) == 0:
        if nC.size == 2:
          # Extending clump
          mask[maskPix2] = nC[1]
        else:
          clumpMerge = nC[1:]
          # Merging clumps
          iC = clumpPeak[clumpMerge, 0]
          jC = clumpPeak[clumpMerge, 1]
          for nR1 in range(0, nPix1):
            i = pix1[0][nR1]
            j = pix1[1][nR1]
            d = (i-iC)**2 + (j-jC)**2
            m = np.argmin(d)
            mask[i, j] = clumpMerge[m]

  return nNew

# Clumps with fewer pixels than nPixMin are removed
# The labels of the clumps are resorted accordingly
def testBad(data, nPixMin, nCl, clumpPeak, mask):

  dMax = data[clumpPeak[1:nCl+1, 0], clumpPeak[1:nCl+1, 1]]
  newOrder = np.argsort(dMax)[::-1]

  nClNew = 0
  nBad = 0
  mask0 = mask.copy()

  for n0 in newOrder:
    iClp = mask0 == n0+1
    count = np.sum(iClp)
    if count <= nPixMin:
      nBad += 1
      mask[iClp] = 0
    else:
      nClNew += 1
      mask[iClp] = nClNew
  return nBad

# Main routine
def clfind2d(data, levels, log=False, nPixMin=20, silent=False):

  # If log == True, sets the name of the log file
  if log is True:
    logFile = 'clfind2d.log'
    if os.path.exists(logFile): os.remove(logFile)
  else:
    logFile = None

  tStart = clock()  # The clockwatch is on
  currentTime = strftime('%a, %d %b %Y %H:%M:%S GMT', gmtime())
  
  printLog('----------------------------------------------------------------\n', logFile, silent)
  printLog('CLFIND2d: %s\n' % currentTime, logFile, silent)
  if log is True:
    printLog('Log file = %s\n' % logFile, logFile, silent)
  printLog('Minimum number of pixels = %d\n' % nPixMin, logFile, silent)
  levStr = np.round(levels, 4)
  levStr = ','.join(map(str, levStr))
  printLog('Contour levels at = %s\n' % levStr, logFile, silent)
  printLog('----------------------------------------------------------------\n', logFile, silent)

  iSize, jSize = data.shape
  
  # Bad pixels are set to -999.9
  data[np.isnan(data)] = -999.9

  # Adds a top level to the contours list
  levels = list(levels)
  levels = np.append(levels, data.max()+1000)
  nLev = levels.size

  max_clumps = 1e6

  # The array which will contain the mask
  mask = np.zeros([iSize, jSize], np.int)
  clumpPeak = np.zeros([max_clumps, 2], np.int)

  # For each contour level, identifies the clumps and sorts them
  nCl = 0
  for nWork in range(nLev-2, -1, -1):
    levMin = levels[nWork]
    levMax = levels[nWork + 1]

    nPix, reg, nReg = defReg(data, levMin, levMax, diagonal=True)
    
    printLog('Contour level %f: %i pixels, %i regions, ' % (levels[nWork], nPix, nReg), logFile, silent)
    nNew = defClump(data, levels[nWork], reg, nReg, mask, clumpPeak, nCl)
    printLog('%i new clumps\n' % (nNew), logFile, silent)
    nCl += nNew
  # Tests for bad pixels
  nBad = testBad(data, nPixMin, nCl, clumpPeak, mask)
  nCl -= nBad

  printLog('%d clumps found (%d rejected)\n' % (nCl, nBad), logFile, silent)
  printLog('================================================================\n', logFile, silent)

  tEnd = clock()
  interval = (tEnd - tStart)/60.0
  printLog('%.1f minutes elapsed\n' % interval, logFile, silent)
  printLog('\n', logFile, silent)

  return mask

def clstat2d(mask,data):
  
  unique = np.unique(mask)
  # number of clumps
  nclump = len(unique)-1
  if nclump == 0:
    return
  unique = unique[1:]
  
  ny, nx = data.shape
  # ------- loop over clumps
  ppeakx, ppeaky, ppeak, eedge, ssumflux, ssx, ssy, rradius = [], [], [], [], [], [], [], []
  for i in unique:
    index = np.where(mask == i)
    flux  = np.where(mask == i,data,0)
    npix  = np.count_nonzero(flux)
    
    iy, ix = np.mgrid[:ny,:nx]
    
    # ------- check if clump is at edge
    edge = True if any(map(any,[ix[index] == 0,ix[index] == nx-1,iy[index] == 0,iy[index] == ny-1])) else False
    
    peak = flux.max()
    peaky, peakx = np.unravel_index(flux.argmax(),flux.shape)
    
    sumflux = flux.sum()
    sx      = 2.355*np.sqrt((ix*ix*flux).sum()/sumflux - ((ix*flux).sum()/sumflux)**2)
    sy      = 2.355*np.sqrt((iy*iy*flux).sum()/sumflux - ((iy*flux).sum()/sumflux)**2)
    radius  = np.sqrt(npix/np.pi)
    
    ppeakx.append(peakx)
    ppeaky.append(peaky)
    ppeak.append(peak)
    eedge.append(edge)
    ssumflux.append(sumflux)
    ssx.append(sx)
    ssy.append(sy)
    rradius.append(radius)

  return ppeakx, ppeaky, ppeak, eedge, ssumflux, ssx, ssy, rradius