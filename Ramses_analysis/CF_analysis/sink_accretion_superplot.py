import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
from mpi4py.MPI import COMM_WORLD as CW
import sys
import collections
import os
import matplotlib.collections as mcoll
import matplotlib.path as mpath

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot_type", "--plot_accretion_type", help="which plot do you want to make?", default='Shut_off', type=str)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

Shut_off = [0, 1, 17, 20, 21, 26, 49, 59, 64, 65, 80, 101, 105, 115, 117, 132, 151, 159, 167, 173, 191]
Continuous = [2, 8, 11, 84, 89, 90, 96, 118, 126, 127, 129, 139, 140, 157, 162, 164, 170, 172, 186, 208]
Late_accretion = [13, 24, 27, 30,  82, 94, 112, 156, 202]
Slow_drop_off = [28, 35, 38, 42, 46, 51, 57, 62, 66, 68, 69, 72, 77, 79, 95, 97, 119, 130, 136, 144, 158, 160, 161, 163, 168, 169, 171, 179, 189]
Sporadic = [31, 32, 37, 40, 41, 44, 47, 48, 53, 54, 55, 67, 75, 75, 83, 92, 109, 116, 123, 154, 174, 182, 187, 188, 205]
Weird = [33, 71, 86, 91, 93, 98, 111, 114, 124, 133, 135, 143, 149, 177, 178, 182, 192, 193, 198, 201, 223, 235]

rank = CW.Get_rank()
size = CW.Get_size()

global_data_pickle_file = '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl'

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
    
file_open = open(global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
    
dm = global_data['dm']*units['mass_unit'].in_units('Msun')
dt = (global_data['time'] - global_data['tflush'])*units['time_unit'].in_units('yr')
Accretion_array = (dm/dt).T

Time_arr = (global_data['time']*units['time_unit'].in_units('yr')).T

window = yt.YTQuantity(500, 'yr')
rit = -1
make_pickle = False
if make_pickle:
    for sink_it in range(len(Time_arr)):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            print("smoothing sink", sink_it)
            start_ind = np.where(Time_arr[sink_it]>0)[0][0]
            Time = Time_arr[sink_it][start_ind:] - Time_arr[sink_it][start_ind]
            Accretion = Accretion_array[sink_it][start_ind:]
            mean_time = []
            mean_acc = []
            for time in Time:
                start_time = time - window/2.
                end_time = time + window/2.
                start_ind = np.argmin(abs(Time - start_time))
                end_ind = np.argmin(abs(Time - end_time))
                if end_ind != start_ind:
                    non_nan_inds = np.argwhere(np.isnan(Accretion[start_ind:end_ind])==False)
                    mean_t = np.mean(Time[start_ind:end_ind][non_nan_inds])
                    mean_a = np.mean(Accretion[start_ind:end_ind][non_nan_inds])
                    mean_time.append(mean_t)
                    mean_acc.append(mean_a)
                    
            plt.clf()
            plt.figure(figsize=(6,2))
            plt.semilogy(mean_time, mean_acc, linewidth=0.25)
            plt.xlim([0.0, 1500000])
            plt.ylim(bottom=1e-9)
            plt.xlabel('Time since formation')
            plt.ylabel('Accretion rate')
            plt.axvspan(0, 10000, alpha=0.5, color='grey')
            plt.title('Sink No.' + str(sink_it))
            plt.savefig('sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+".pdf")
            pickle_name = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
            file = open(pickle_name, 'wb')
            pickle.dump((mean_time, mean_acc),file)
            file.close()
            print("created "+pickle_name)
        
#gather pickles
if args.plot_accretion_type == 'Continuous':
    Continuous = [2, 8, 11, 89, 90, 127, 129, 139, 140, 157, 162, 164, 172, 186, 208]# Good
    plt.clf()
    #plt.figure(figsize=(6,2))
    for sink_it in Continuous:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        plt.semilogy(mean_time, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
        plt.xlim([0.0, 1500000])
        plt.ylim(bottom=1e-9)
        plt.ylabel('Smoothed Accretion rate (M$_\odot$)')
        plt.xlabel('Time since star formation')
        plt.legend(loc='upper right')
        plt.savefig("Continuous.png")
        print("saved Continuous.png")
elif args.plot_accretion_type == 'Shut_off':
    Shut_off = [0, 1, 17, 21, 49, 59, 101, 105, 117, 167, 191] #Good!
    plt.clf()
    #plt.figure(figsize=(6,2))
    for sink_it in Shut_off:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        plt.semilogy(mean_time, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
        plt.xlim([0.0, 1500000])
        plt.ylim(bottom=1e-9)
        plt.ylabel('Smoothed Accretion rate (M$_\odot$)')
        plt.xlabel('Time since star formation')
        plt.legend(loc='upper right')
        plt.savefig("Shut_off.png")
        print("saved Shut_off.png")
elif args.plot_accretion_type == 'Slow_drop_off':
    Slow_drop_off = [28, 35, 38, 42, 46, 51, 57, 66, 68, 72, 77, 95, 119, 144, 158, 160, 163, 168, 169, 171, 179, 189]
    plt.clf()
    #plt.figure(figsize=(6,2))
    for sink_it in Slow_drop_off:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        plt.semilogy(mean_time, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
        plt.xlim([0.0, 1500000])
        plt.ylim(bottom=1e-9)
        plt.ylabel('Smoothed Accretion rate (M$_\odot$)')
        plt.xlabel('Time since star formation')
        plt.legend(loc='upper right')
        plt.savefig("Slow_drop_off.png")
        print("saved Slow_drop_off.png")
elif args.plot_accretion_type == 'Late_accretion':
    Late_accretion = [13, 24, 27, 30,  82, 112, 156]# Good!
    plt.clf()
    #plt.figure(figsize=(6,2))
    for sink_it in Late_accretion:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        plt.semilogy(mean_time, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
        plt.xlim([0.0, 1500000])
        plt.ylim(bottom=1e-9)
        plt.ylabel('Smoothed Accretion rate (M$_\odot$)')
        plt.xlabel('Time since star formation')
        plt.legend(loc='upper right')
        plt.savefig("Late_accretion.png")
        print("saved Late_accretion.png")
elif args.plot_accretion_type == 'Weird':
    Weird = [33, 55, 86, 98, 114, 135, 143, 178, 182, 192, 193, 201, 223]# Good!
    plt.clf()
    #plt.figure(figsize=(6,2))
    for sink_it in Weird:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        plt.semilogy(mean_time, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
        plt.xlim([0.0, 1500000])
        plt.ylim(bottom=1e-9)
        plt.ylabel('Smoothed Accretion rate (M$_\odot$)')
        plt.xlabel('Time since star formation')
        plt.legend(loc='upper right')
        plt.savefig("Weird.png")
        print("saved Weird.png")
else:
    panel_tag = ['(a)', '(b)', '(c)', '(d)', '(e)']
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(7.5, 6), sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0.0)
    Continuous = [2, 8, 11, 89, 90, 127, 129, 139, 140, 157, 162, 164, 172, 186, 208]
    for sink_it in Continuous:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        axs[0].semilogy(np.array(mean_time)/1.e6, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
    axs[0].set_xlim([0.0, 1.5])
    axs[0].set_ylim(bottom=1e-9)
    axs[0].text(1.44,1.e-7,panel_tag[0])
    plt.savefig("Accretion_evol.pdf", format='pdf', bbox_inches='tight')
    print("saved Accretion_evol.pdf")
    xticklabels = axs[0].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    Shut_off = [0, 1, 17, 21, 49, 59, 101, 105, 117, 167, 191] #Good!
    for sink_it in Shut_off:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        axs[1].semilogy(np.array(mean_time)/1.e6, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
    axs[1].text(1.44,1.e-7,panel_tag[1])
    plt.savefig("Accretion_evol.pdf", format='pdf', bbox_inches='tight')
    print("saved Accretion_evol.pdf")
    xticklabels = axs[1].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    Slow_drop_off = [28, 35, 38, 42, 46, 51, 57, 66, 68, 72, 77, 95, 119, 144, 158, 160, 163, 168, 169, 171, 179, 189]
    for sink_it in Slow_drop_off:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        axs[2].semilogy(np.array(mean_time)/1.e6, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
    axs[2].set_ylabel('Smoothed accretion rate (M$_\odot$/yr)')
    axs[2].text(1.44,1.e-7,panel_tag[2])
    plt.savefig("Accretion_evol.pdf", format='pdf', bbox_inches='tight')
    print("saved Accretion_evol.pdf")
    xticklabels = axs[2].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    Late_accretion = [13, 24, 27, 30,  82, 112, 156]# Good!
    for sink_it in Late_accretion:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        axs[3].semilogy(np.array(mean_time)/1.e6, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
    axs[3].text(1.44,1.e-7,panel_tag[3])
    plt.savefig("Accretion_evol.pdf", format='pdf', bbox_inches='tight')
    print("saved Accretion_evol.pdf")
    xticklabels = axs[3].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    Weird = [33, 55, 86, 98, 114, 135, 143, 178, 182, 192, 193, 201, 223]# Good!
    for sink_it in Weird:
        pickle_file = 'sink_no_'+str(sink_it)+'_smoothed_'+str(window.value)+'.pkl'
        file = open(pickle_file, 'rb')
        mean_time, mean_acc = pickle.load(file)
        file.close()
        axs[4].semilogy(np.array(mean_time)/1.e6, mean_acc, linewidth=0.25, label=str(sink_it))
        print('plotted from', pickle_file)
    axs[4].set_xlabel('Time since star formation (Myr)')
    axs[4].text(1.44,1.e-7,panel_tag[4])
    plt.savefig("Accretion_evol.pdf", format='pdf', bbox_inches='tight')
    print("saved Accretion_evol.pdf")
    

'''
pickles = glob.glob('sink_no_*_smoothed_'+str(window.value)+'.pkl')
plt.clf()
plt.figure(figsize=(6,2))
max_acc_dot = []
for pickle_file in pickles:
    file = open(pickle_file, 'rb')
    mean_time, mean_acc = pickle.load(file)
    da = np.array(mean_acc)[2:] - np.array(mean_acc)[:-2]
    dt = np.array(mean_time)[2:] - np.array(mean_time)[:-2]
    max_acc_dot.append(np.nanmax(da/dt))
    file.close()
    plt.semilogy(mean_time, mean_acc, linewidth=0.25, alpha=0.5)
    print('plotted from', pickle_file)
    plt.xlim([0.0, 1500000])
    plt.ylim(bottom=1e-9)
    plt.ylabel('Smoothed Accretion rate (M$_\odot$)')
    plt.xlabel('Time since star formation')
    plt.savefig(pickle_file.split('.pkl')[0]+".pdf")
    print("saved", pickle_file.split('.pkl')[0]+".pdf")
'''
"""
plt.clf()
plt.figure(figsize=(12,4))
window = yt.YTQuantity(10000, 'yr')
for sink_it in range(len(Time_arr)):
    start_ind = np.where(Time_arr[sink_it]>0)[0][0]
    Time = Time_arr[sink_it][start_ind:] - Time_arr[sink_it][start_ind]
    Accretion = Accretion_array[sink_it][start_ind:]
    da = Accretion[2:] - Accretion[:-2]
    dt = Time[2:] - Time[:-2]
    if np.nanmax(da/dt) > 1.e-7:
        '''
        mean_time = []
        mean_acc = []
        for time in Time:
            start_time = time - window/2.
            end_time = time + window/2.
            start_ind = np.argmin(abs(Time - start_time))
            end_ind = np.argmin(abs(Time - end_time))
            if end_ind != start_ind:
                non_nan_inds = np.argwhere(np.isnan(Accretion[start_ind:end_ind])==False)
                mean_t = np.mean(Time[start_ind:end_ind][non_nan_inds])
                mean_a = np.mean(Accretion[start_ind:end_ind][non_nan_inds])
                mean_time.append(mean_t)
                mean_acc.append(mean_a)
        
        plt.semilogy(mean_time, mean_acc, linewidth=0.25, alpha=0.25)
        '''
        plt.semilogy(Time, Accretion, linewidth=0.25, alpha=0.25)
        plt.xlim([0.0, 1700000])
        plt.ylim(bottom=1e-9)
        print('plotted sink', sink_it)
        plt.savefig("test_acc.png")
"""
