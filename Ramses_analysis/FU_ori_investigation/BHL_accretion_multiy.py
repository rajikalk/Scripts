#!/usr/bin/env python
import glob
import sys
import argparse
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os

cmap=plt.cm.gist_heat
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "mathtext.fontset": "stixsans"  # Force math to use a sans-serif look
    })

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#-----------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------
#Ploting parameters
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 9

plt.clf()
fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))

time_bounds = [[3800, 4900],[5575, 5700], [6580, 6720], [7295, 7365], [7850, 7900]]
colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

#Start by loading pickel data and then deleting what we don't need

event_it = args.event_identifier
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_burst = 7340
    end_time = 7340
    
try:
    sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    sys.stdout.flush()
except:
    sink_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
    print("read pickle", sink_pickle)
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    sys.stdout.flush()
    
plt.title("Burst event "+str(event_it), y=0.8)
start_ind = np.argmin(abs(particle_data['time']-start_time))
end_ind = np.argmin(abs(particle_data['time']-end_time))
#axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
lns_all = []
lns1 = plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='k', ls='-', label="Accretion rate")
lns_all.append(lns1)
Radii = [3, 4, 5, 6]

for rad in Radii:
    directory = "/home/100/rlk100/rlk/RAMSES/Analysis/BHL_analytical_calc/Event_"+str(args.event_identifier)+"/Radius_"+str(rad) +"/"
    plot = True
    if os.path.exists(directory+'BHL_accretion.pkl'):
        file_open = open(directory+'BHL_accretion.pkl', 'rb')
        save_dict = pickle.load(file_open)
        file_open.close()
    elif os.path.exists(directory+'BHL_accretion_0.pkl'):
        pickle_files = sorted(glob.glob(directory+"BHL_accretion_*.pkl"))
        save_dict = {}
        save_dict.update({"Time": np.array([])})
        save_dict.update({"BHL_Acc_acc_low": np.array([])})
        save_dict.update({"BHL_Acc_acc_high": np.array([])})
        save_dict.update({"Density": np.array([[]])})
        save_dict.update({"Rel_kep": np.array([[]])})
        for pickle_file in pickle_files:
            file_open = open(pickle_file, 'rb')
            save_dict_r = pickle.load(file_open)
            file_open.close()
            for key in save_dict_r.keys():
                save_dict[key] = np.append(save_dict[key],save_dict_r[key])
        sorted_inds = np.argsort(save_dict["Time"])
        for key in save_dict.keys():
            save_dict[key] = save_dict[key][sorted_inds]
    else:
        print("No data right now for Radius", rad)
        plot = False
    
    if plot == True:
        BHL_mean = (save_dict["BHL_Acc_acc_low"]+save_dict["BHL_Acc_acc_high"])/2
        lns3 = plt.semilogy(save_dict["Time"], BHL_mean, color=colours[pickle_files.index(pickle_file)], ls=':', label="r = "+str(rad)+"AU")
        lns_all.append(lns3)
        plt.ylim([np.min(particle_data['mdot'].T[1][start_ind:end_ind]), np.max(particle_data['mdot'].T[1][start_ind:end_ind])])
        plt.fill_between(save_dict["Time"], save_dict["BHL_Acc_acc_low"], save_dict["BHL_Acc_acc_high"], color=colours[pickle_files.index(pickle_file)], alpha=0.5)
axes_1_twin = plt.twinx()
lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")
lns_all.append(lns2)
#Plot accretion and separation. This should be loaded from a pickle

plt.xlabel('Time (yr)', labelpad=-0.2, fontsize=font_size) #($yr$)
plt.ylabel('Accretion rate (M$_\odot$/yr)', labelpad=-0.2, fontsize=font_size)# (M$_\odot/yr$)
axes_1_twin.set_ylabel('Separation (au)', fontsize=font_size)
plt.tick_params(axis='x', which='major', direction='in', color='k', top=True)
plt.tick_params(axis='y', which='major', direction='in', color='k', right=True)
#plt.xaxis.label.set_color('black')
#plt.yaxis.label.set_color('black')
plt.tick_params(axis='both', labelsize=font_size)
plt.xlim([start_time, end_time])
plt.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
labs = [l.get_label() for l in lns_all]
plt.legend(lns_all, labs, loc='upper left')
plt.savefig("BHL_Event_"+str(event_it)+"_multi_rad.png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
print('Saved figure with BHL Accretion')
sys.stdout.flush()
