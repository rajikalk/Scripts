#!/usr/bin/env python
import yt
yt.enable_parallelism()
import numpy as np
import sys
import os
import pickle

vel_bins = np.linspace(0, 20, 21).tolist()
bin_centers = ((np.array(vel_bins[1:]) + np.array(vel_bins[:-1]))/2).tolist()
bin_centers.append(bin_centers[-1]+((np.array(vel_bins[1:]) - np.array(vel_bins[:-1]))[0]))
vel_bins.append(100)
        
file = open('gathered_burst_analysys_sink_'+str(sink_id)+'.pkl', 'rb')
sink_all = pickle.load(file)
file.close()
    
for ind in range(len(sink_all['time'])):
    time_val = sink_all['time'][ind].in_units('yr')
    
    hist = sink_all['outflow_distribution'][ind]
    
    plt.clf()
    plt.bar(bin_centers, hist, width=0.5, color='b')
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('#')
    plt.xlim([0, 21])
    plt.ylim([bottom=0])
    save_name = "dist_" + ("%06d" % frames[file_int-1]) + ".jpg"
    plt.savefig(save_name, bbox_inches='tight', pad_inches=0.02)
