#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as mp
from math import pi

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

class spotplot(object):
  # ------ Initializing -------------------------------------------------------------------------------------------------------
  def __init__(self,**kw):
    self.cos_incl    = kw.get('cos_incl', 0)                    #spot paramters
    self.sin_lat     = kw.get('sin_lat' , 0)
    self.long0       = kw.get('long0'   , 0)
    self.size        = kw.get('size'    ,-3)
    self.ua          = kw.get('ua'      , 0)
    self.ub          = kw.get('ub'      , 0)
    self.Fspot       = kw.get('Fspot'   , 0)
    self.nspot       = kw.get('nspot'   , 1)                    #number of spots to model
    ngrid            = kw.get('ngrid'   , 501)
    nspot            = self.nspot

    self.cos_incl    = np.asanyarray(self.cos_incl)          #convert to numpy arrays
    self.sin_lat     = np.asanyarray(self.sin_lat )
    self.long0       = np.asanyarray(self.long0   )
    self.size        = np.asanyarray(self.size    )

    # ------ min/max-values for paramters ---------------------------------------------------------------------------------------
    slatmin , slatmax  = -1.  , 1.    #sin_lat
    longmin , longmax  = -180., 180.  #long
    sizemin , sizemax  = -25. , 0.    #size
    cinclmin, cinclmax = 0.   , 1.    #cos_incl
    uamin   , uamax    = 0.   , 1.    #ua
    ubmin   , ubmax    = 0.   , 1.    #ub
    Fspotmin, Fspotmax = 0.   , 2.    #Fspot
    
    # ------ check arguments ----------------------------------------------------------------------------------------------------
    if (np.array((self.sin_lat.size,self.long0.size,self.size.size)) != self.nspot).any():
      raise TypeError, "Number of elements in sin_lat, long0, size arrays not equal to number of spots"
    if (self.sin_lat > slatmax).any() or (self.sin_lat < slatmin).any():
      raise ValueError, "sin_lat outside allowed range."
    if (self.size > sizemax).any() or (self.size < sizemin).any():
      raise ValueError, "size outside allowed range."
    if self.cos_incl > cinclmax or self.cos_incl < cinclmin:
      raise ValueError, "cos_incl outside allowed range."
    if self.ua > uamax or self.ua < uamin:
      raise ValueError, "ua outside allowed range."
    if self.ub > ubmax or self.ub < ubmin:
      raise ValueError, "ub outside allowed range."
    if self.Fspot > Fspotmax or self.Fspot < Fspotmin:
      raise ValueError, "Fspot outside allowed range."

    # ------ Initialize disk ----------------------------------------------------------------------------------------------------
    xy      = np.linspace(-1,1,num=ngrid) #xy coordinate
    x, y    = np.meshgrid(xy,xy)
    self.xy = xy

    radius = np.where( x**2 + y**2 > 1)

    self.x_disk, self.y_disk = x.copy(), y.copy() 
    self.x_disk[radius], self.y_disk[radius] = np.nan, np.nan
    self.z_disk = np.sqrt( 1 - self.x_disk**2 - self.y_disk**2)
    self.disk   = (1 - self.ua*(1 - self.z_disk) - self.ub*(1 - self.z_disk)**2)/(pi*(1 - self.ua/3 - self.ub/6))

  def spotmodel(self,phase):
    # ------ Initializing ----------------------------------------------------------------------
    Fstar            = 1.

    # ------ Converting angles from degree to radians ------------------------------------------
    lat   = np.arcsin(self.sin_lat)
    incl  = np.arccos(self.cos_incl)
    long0 = self.long0*2*pi/360

    # ------ Calculate size --------------------------------------------------------------------
    alpha = np.arcsin( np.sqrt( np.exp(self.size) ) ) #Angular size of spot

    # ------ Coefficients for plane a*x + b*y + c*z = d ----------------------------------------
    a = np.cos(lat) * np.sin(long0 + 2*pi*phase)
    b = np.sin(lat)
    c = np.cos(lat) * np.cos(long0 + 2*pi*phase)
    d = np.sqrt(1 - np.sin(alpha)**2)

    # ------ Create disk coordinates after rotation of inclination -----------------------------
    xx =  self.x_disk.copy()
    yy =  self.y_disk*np.sin(incl) + self.z_disk*np.cos(incl)
    zz = -self.y_disk*np.cos(incl) + self.z_disk*np.sin(incl)

    disk = self.disk.copy()
    
    # ------ Initialize data arrays ------------------------------------------------------------
    if self.nspot == 1:
      plane       = a*xx + b*yy + c*zz
      index       = np.where(np.logical_and(plane >= d, plane <= 1.))
      zzz         = np.sqrt(1 - xx[index]**2 - yy[index]**2)
      disk[index] = self.Fspot * (1 - self.ua*(1 - zzz) - self.ub*(1 - zzz)**2)/(pi*(1 - self.ua/3 - self.ub/6))
    else:
      for i in range(self.nspot): #loop over spots
        plane       = a[i]*xx + b[i]*yy + c[i]*zz
        index       = np.where(np.logical_and(plane >= d[i], plane <= 1.))
        zzz         = np.sqrt(1 - xx[index]**2 - yy[index]**2)
        disk[index] = self.Fspot * (1 - self.ua*(1 - zzz) - self.ub*(1 - zzz)**2)/(pi*(1 - self.ua/3 - self.ub/6))

    return disk

  def plot(self,phases,**kw):
    # ------ Handle Input ----------------------------------------------------------------------
    nx       = kw.get('nx',4)
    ny       = kw.get('ny',2)
    filename = kw.get('filename','spot_map')
    lmin     = kw.get('lmin',0)
    lmax     = kw.get('lmax',1)
    lnum     = kw.get('lnum',100)
    gray     = kw.get('gray',False)

    phases = np.asanyarray(phases)
    level  = np.linspace(lmin,lmax,num=lnum)

    # ------ helper function -------------------------------------------------------------------
    def plotcore(phase):
      disk = self.spotmodel(phase)
      print "phase: %1.2f, min: %1.2f, max: %1.2f" % (phase, np.nanmin(disk), np.nanmax(disk))
      ax.contourf(self.xy,self.xy,disk,cmap=cmap,levels=level)
      ax.plot(0,np.sqrt(1 - self.cos_incl**2),'xk')
      ax.text(0.02,0.02,"$\phi = %1.2f$" % (phase),color='white',size='x-small',transform=ax.transAxes)
      ax.axis('equal')
      ax.set_xticks([])
      ax.set_yticks([])
      return

    # ------ check arguments -------------------------------------------------------------------
    if phases.size > nx*ny:
      raise ValueError, "Number of phase points greater than nx*ny"
    
    if gray is True:
      cmap = mp.cm.gist_gray
    else:
      cmap = mp.cm.YlOrBr_r

    # ------ Plot initialization ---------------------------------------------------------------
    fig_width_cm  = 10
    inches_per_cm = 1.5/2.54                    # Convert cm to inch
    fig_width     = fig_width_cm*inches_per_cm  # width in inches
    fig_height    = fig_width*1.1*ny/nx         # height in inches
    fig_size      =  [fig_width,fig_height]
    params        = {'savefig.dpi'   : 300      ,
                     'font.family'   : 'serif'  ,
                     'text.usetex'   : True     ,
                     'figure.figsize': fig_size ,
                     'figure.subplot.wspace' : 0,
                     'figure.subplot.hspace' : 0}
    mp.rcParams.update(params)

    # ------ Start working ---------------------------------------------------------------------
    fig = mp.figure()
    if phases.size == 1:
      ax = fig.add_subplot(1,1,1,axisbg='black')
      plotcore(phases)
    else:
      for i,phase in enumerate(phases):
        ax = fig.add_subplot(ny,nx,i+1,axisbg='black')
        plotcore(phase)
    mp.savefig(filename+'.png',bbox_inches='tight')
    mp.close()
    mp.clf() #Clear figure
    
  # ------ Call --------------------------------------------------------------------------------
  def __call__(self,phases,**kw):
    self.plot(phases,**kw)
    return