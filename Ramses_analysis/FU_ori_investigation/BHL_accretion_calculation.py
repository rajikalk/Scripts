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
yt.enable_parallelism()
import my_ramses_fields_short as myf
import gc
from mpi4py.MPI import COMM_WORLD as CW

#-----------------------------------------------------
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument("-ax", "--axis", default='xy', type=str)
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=5, type=float)
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
    particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
except:
    sink_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
    print("read pickle", sink_pickle)
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    
if rank != 0:
    del particle_data, counter
    gc.collect()

event_it = args.event_identifier
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_burst = 7340
    end_time = 7340
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

if rank == 0:
    plt.clf()
    fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
                
    plt.title("Burst event "+str(event_it), y=0.8)
    start_ind = np.argmin(abs(particle_data['time']-start_time))
    end_ind = np.argmin(abs(particle_data['time']-end_time))
    #axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
    lns1 = plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='b', ls='-', label="Accretion rate")
    axes_1_twin = plt.twinx()
    lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")
                
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

    del particle_data, counter
    gc.collect()
    
sys.stdout.flush()
CW.Barrier()

sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/'
files = sorted(glob.glob(sim_data_dir+"*/info*.txt"))

if os.path.exists('BHL_accretion.pkl'):
    file_open = open('BHL_accretion.pkl', 'rb')
    time_arr, BHL_Acc_acc_low, BHL_Acc_acc_high = pickle.load(file_open)
    file_open.close()
    if len(time_arr) != len(files):
        files = files[len(time_arr):]
else:
    time_arr = np.array([])
    BHL_Acc_acc_low = np.array([])
    BHL_Acc_acc_high = np.array([])
sink_id = 45
#sink_form_time = np.nan

sys.stdout.flush()
CW.Barrier()

if len(files)>0:
    rit = -1
    for file in files:
        rit = rit+1
        if rit == size:
            rit = 0
        if rank == rit:
            ds = yt.load(file, units_override=units_override)
            #if np.isnan(sink_form_time):
            #    sink_form_time = ds.r["sink_particle_form_time"][sink_id]
            time_val = ds.current_time.in_units('yr').value - sink_form_time.in_units('yr').value
            time_arr = np.append(time_arr, time_val)
            
            sink_mass = ds.r["gas", "sink_particle_mass"][sink_id]
            gc.collect()
            
            #Get sink position
            sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
            sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
            sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
            sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
            del sink_particle_posx, sink_particle_posy, sink_particle_posz
            gc.collect()
            
            #get sink velocity
            sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
            sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
            sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
            sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
            del sink_particle_velx, sink_particle_vely, sink_particle_velz
            gc.collect()
            
            #Define measuring sphere:
            radius = yt.YTQuantity(args.measuring_sphere_radius, 'au')
            
            #Get inds in measuring sphere
            dd = ds.all_data()
            dx = dd['x'].in_units('au') - sink_pos[0].in_units('au')
            dy = dd['y'].in_units('au') - sink_pos[1].in_units('au')
            dz = dd['z'].in_units('au') - sink_pos[2].in_units('au')
            sep = np.sqrt(dx**2 + dy**2 + dz**2)
            del dx, dy, dz
            gc.collect()
            
            sphere_inds = np.where(sep<radius)[0]
            del sep
            gc.collect()
            
            mean_density = np.mean(ds.r["gas", "Density"][sphere_inds])
            
            #Calculate bulk velocity of the sphere
            sph_velx = np.mean(ds.r["ramses", "x-velocity"][sphere_inds].in_units('km/s'))
            sph_vely = np.mean(ds.r["ramses", "y-velocity"][sphere_inds].in_units('km/s'))
            sph_velz = np.mean(ds.r["ramses", "z-velocity"][sphere_inds].in_units('km/s'))
            bulk_velocity = yt.YTArray([sph_velx, sph_vely, sph_velz])
            del sph_velx, sph_vely, sph_velz
            gc.collect
            rel_vel = bulk_velocity - sink_vel
            rel_speed = np.sqrt(np.sum(rel_vel**2))
            del rel_vel
            gc.collect
            
            sound_speed = np.mean(np.sqrt((ds.r["gas", "Gamma"][sphere_inds]*ds.r["gas", "Pressure"][sphere_inds])/ds.r["gas", "Density"][sphere_inds]).in_units('km/s'))
            
            alpha = yt.YTArray([1, 2], '')
            BHL_top = 2*np.pi* (sink_mass.in_cgs()*yt.units.gravitational_constant_cgs)**2 * mean_density.in_cgs()
            del sink_mass, mean_density
            gc.collect
            BHL_bot = (rel_speed.in_cgs()**2 + sound_speed.in_cgs()**2)**(3./2.)
            del rel_speed, sound_speed
            gc.collect
            BHL = (alpha * (BHL_top/BHL_bot)).in_units('msun/yr')
            del alpha, BHL_top, BHL_bot
            gc.collect
            BHL_Acc_acc_low = np.append(BHL_Acc_acc_low, BHL[0])
            BHL_Acc_acc_high = np.append(BHL_Acc_acc_high, BHL[1])

            #Save BHL Calculation
            file_open = open('BHL_accretion_'+str(rank)+'.pkl', 'wb')
            pickle.dump((time_arr, BHL_Acc_acc_low, BHL_Acc_acc_high), file_open)
            file_open.close()
            print("RANK "+str(rank)+": Calculated BHL for file", file)
        
sys.stdout.flush()
CW.Barrier()

if rank == 0:
    time_arr = np.array([])
    BHL_Acc_acc_low = np.array([])
    BHL_Acc_acc_high = np.array([])
    for rit in range(size):
        file_open = open('BHL_accretion_'+str(rit)+'.pkl', 'rb')
        time_arr_r, BHL_Acc_acc_low_r, BHL_Acc_acc_high_r = pickle.loadfile_open)
        file_open.close()
        time_arr = np.append(time_arr,time_arr_r)
        BHL_Acc_acc_low = np.append(BHL_Acc_acc_low, BHL_Acc_acc_low_r)
        BHL_Acc_acc_high = np.append(BHL_Acc_acc_high, BHL_Acc_acc_high_r)
    sorted_inds = np.argsort(time_arr)
    time_arr = time_arr[sorted_inds]
    BHL_Acc_acc_low = BHL_Acc_acc_low[sorted_inds]
    BHL_Acc_acc_high = BHL_Acc_acc_high[sorted_inds]

    lns3 = plt.fill_between(time_arr, BHL_Acc_acc_low, BHL_Acc_acc_high, color='g', alpha=0.5, label="BHL prediction")
    lns = lns1+lns2+ln3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='upper left')
    plt.savefig("BHL_Event_"+str(event_it)+"_"+args.axis+".pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)
