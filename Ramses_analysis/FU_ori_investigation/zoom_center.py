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
tform = sink_data['tflush'].T[i_sink][form_it]
tfile = sink_data['time'].T[i_sink][form_it]
vel = np.array([sink_data['ux'].T[i_sink][form_it], sink_data['uy'].T[i_sink][form_it], sink_data['uz'].T[i_sink][form_it]])
pos = np.array([sink_data['x'].T[i_sink][form_it], sink_data['y'].T[i_sink][form_it], sink_data['z'].T[i_sink][form_it]])
form_pos = (tfile-tform)*vel + pos
print( 'Coordinate at time of formation: ', form_pos)
print( 'Time difference between first snapshot and time of formation: ', (tfile-tform))
 
# Time of snapshot just before formation -- put by hand 
# Found by looking at data/output_XXXXX/info_XXXXX.txt 
tsnap=0.104721058347503E+01 
 
tmid = 0.5*(tsnap + tfile)
dt = tmid - tform
print( "Zoom center :", dt*vel+pos)
print( tsnap, tmid, tfile, tform, dt, tfile - tform)
 
 
