import csv
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
import yt

pickle_files = ['/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L08.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L09.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L10.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L11.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L12.pkl']
label = ["Lvl=8 (10.07AU)", "Lvl=9 (5.04AU)", "Lvl=10 (2.52AU)", "Lvl=11 (1.26AU)", "Lvl=12 (0.63AU)"]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(two_col_width, single_col_width*2), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

cit = -1
for pick_file in pickle_files:
    print("read pickle", pick_file)
    cit = cit + 1
    file_open = open(pick_file, 'rb')
    sink_data, line_counter = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    
    if len(list(sink_data.keys())) != 2:
        print("Number fo sink != 2")
    else:
        sys_form_time = -1*np.inf
        for key in list(sink_data.keys()):
            if sink_data[key]['time'][0] > sys_form_time:
                sys_form_time = sink_data[key]['time'][0]
                secondary_key = key
        primary_key = np.nan
        for key in list(sink_data.keys()):
            if key != secondary_key:
                primary_key = key
        import pdb
        pdb.set_trace()
        
        #KEEP EDITING HERE. PLOT THE BINARY SEPARATION.
    
