#!/usr/bin/env python

import os
import time
import numpy as np
import pickle

fname = '/home/100/rlk100/gdata/RAMSES/Global/stars_red_512.pkl'
glob_data = open(fname, 'rb')
sink_data = pickle.load(glob_data)
glob_data.close()
 
# Sink to zoom in on 
i_sink=17
 
# Extract first sink record where sink exists 
form_it = np.where(sink_data['m'].T[i_sink]>0)[0][0]
tflush = sink_data['tflush'].T[i_sink][form_it]
tform = sink_data['time'].T[i_sink][form_it]
print( 'Coordinate at time of formation: ', (ss.t - ss.time)*[ss.ux,ss.uy,ss.uz]+[ss.x,ss.y,ss.z]
print( 'Time difference between first snapshot and time of formation: ', (ss.t - ss.time) 
 
# Time of snapshot just before formation -- put by hand 
# Found by looking at data/output_XXXXX/info_XXXXX.txt 
tsnap=0.104721058347503E+01 
 
tmid = 0.5*(tsnap + ss.t) 
dt = tmid - ss.time 
print( "Zoom center :", dt*[ss.ux,ss.uy,ss.uz]+[ss.x,ss.y,ss.z] 
print( tsnap, tmid, ss.t, ss.time, dt, ss.t - ss.time 
 
 
