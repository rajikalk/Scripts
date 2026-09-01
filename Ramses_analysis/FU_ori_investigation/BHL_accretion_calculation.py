#!/usr/bin/env python
import glob
import sys
import argparse
import numpy as np
import pickle
import my_ramses_module as mym
import matplotlib.pyplot as plt
import os
import yt
import my_ramses_fields_short as myf
import gc
#-----------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument("-ax", "--axis", default='xy', type=str)
parser.add_argument("-sph_rad", "--measuring_sphere_radius", )
parser.add_argument('files', nargs='*')
args = parser.parse_args()


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "mathtext.fontset": "stixsans"  # Force math to use a sans-serif look
})

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#Ploting parameters
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 9
mym.set_global_font_size(font_size)


#------------------------------------------------------
time_bounds = [[3800, 4900],[5575, 5700], [6580, 6720], [7295, 7365], [7850, 7900]]
cmap=plt.cm.gist_heat

#Start by loading pickel data and then deleting what we don't need
try:
    sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
except:
    sink_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
    print("read pickle", sink_pickle)
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    

event_it = args.event_identifier
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_burst = 7340
    end_time = 7340
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)


plt.clf()
fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
            
plt.title("Burst event "+str(event_it), y=0.8)
start_ind = np.argmin(abs(particle_data['time']-start_time))
end_ind = np.argmin(abs(particle_data['time']-end_time))
#axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
lns1 = plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='b', ls='-', label="Accretion rate")
axes_1_twin = plt.twinx()
lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")
lns = lns1+lns2
labs = [l.get_label() for l in lns]
plt.legend(lns, labs, loc='upper left')
            
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

plt.savefig("BHL_Event_"+str(event_it)+"_"+args.axis+".pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)

sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/'
files = sorted(glob.glob(sim_data_dir+"*/info*.txt"))

time_arr = np.array([])
BHL_Acc_acc = np.array([])
sink_id = 45

for file in files:
    ds = yt.load(file, units_override=units_override)
    
    sink_mass = ds.r["gas", "sink_particle_mass"][sink_id]
    
    #Get sink position
    sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
    sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
    sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
    sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
    
    #Cet sink velocity
    sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
    sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
    sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
    sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
    
    #Define measuring sphere:
    
    
    
    

    import pdb
    pdb.set_trace()
    
