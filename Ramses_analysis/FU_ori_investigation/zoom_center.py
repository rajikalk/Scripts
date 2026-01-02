#!/usr/bin/env python

import os
import time
import numpy as np

fname = '/lustre/astro/troels/IMF_512/stars_red_512.sav' 
fname = '/lustre/astro/troels/IMF_512_cores/stars_red_512.sav' 
if len(s) == 0 : 
    restore, fname 
 
# Sink to zoom in on 
is=17 
 
# Extract first sink record where sink exists 
ii , = np.where(s[is,:].m > 0) & ss=reform(s[is,ii[0]]) 
print( 'Coordinate at time of formation: ', (ss.t - ss.time)*[ss.ux,ss.uy,ss.uz]+[ss.x,ss.y,ss.z] 
print( 'Time difference between first snapshot and time of formation: ', (ss.t - ss.time) 
 
# Time of snapshot just before formation -- put by hand 
# Found by looking at data/output_XXXXX/info_XXXXX.txt 
tsnap=0.104721058347503E+01 
 
tmid = 0.5*(tsnap + ss.t) 
dt = tmid - ss.time 
print( "Zoom center :", dt*[ss.ux,ss.uy,ss.uz]+[ss.x,ss.y,ss.z] 
print( tsnap, tmid, ss.t, ss.time, dt, ss.t - ss.time 
 
 
