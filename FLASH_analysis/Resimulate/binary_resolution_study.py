import csv
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
import yt

pickle_files = ['/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L08.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L09.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L10.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L11.pkl', '/scratch/ek9/rlk100/Analyisis/Sink_particle_pickles/L12.pkl']
label = ["Lvl=8 (10.07AU)", "Lvl=9 (5.04AU)", "Lvl=10 (2.52AU)", "Lvl=11 (1.26AU)", "Lvl=12 (0.63AU)"]
r_acc = 5*np.array([10.07, 5.04, 2.52, 1.26, 0.63])

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(two_col_width, single_col_width*1.5), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

cit = -1
for pick_file in pickle_files:
    print("read pickle", pick_file)
    cit = cit + 1
    file_open = open(pick_file, 'rb')
    file_open = open(pick_file, 'rb')
    sink_data, line_counter = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    
    if len(list(sink_data.keys())) != 2:
        print("Number of sink != 2")
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
        
        form_ind = np.where(sink_data[primary_key]['time']==sys_form_time)[0][0]
        dx = sink_data[primary_key]['posx'][form_ind:] - sink_data[secondary_key]['posx']
        dy = sink_data[primary_key]['posy'][form_ind:] - sink_data[secondary_key]['posy']
        dz = sink_data[primary_key]['posz'][form_ind:] - sink_data[secondary_key]['posz']
        sep = np.sqrt(dx**2 + dy**2 + dz**2)
        sep = yt.YTArray(sep, 'cm').in_units('au')
        time = yt.YTArray((sink_data[secondary_key]['time'] - sys_form_time), 's').in_units('yr')
        axs.axhline(y=r_acc[cit], color=proj_colours[cit], linestyle='--')
        axs.semilogy(time, sep, color=proj_colours[cit], label=label[cit])
    print("read pickle", pick_file)
    
axs.set_xlabel('Time (yr)')
axs.set_ylabel('Separation (au)')
axs.legend()
axs.set_xlim([0, 50000])
plt.savefig('separation_resolution.png', bbox_inches='tight', dpi=300)
