#!/usr/bin/env python
import numpy as np
import h5py
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import commands
import yt
import re
import os
import argparse

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="What is the input directory?", required=True)
    parser.add_argument("-o", "--output_dir", help="what is the output directory? if this is not defined it will be the same as the input directory")
    parser.add_argument("-of", "--output_file", help="Do you have a file ot write out to?", default="out.csv")
    parser.add_argument("-ff", "--file_frequency", help="every how many files do you want to use?", default=100, type=int)
    #parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()
output_file = args.output_file

#define files to read
dir = args.input_dir
if dir[-1] != '/':
    dir = dir + '/'
if args.output_dir == None:
    split_cur = dir.split('/O')
    save_dir = split_cur[0] + '/YT_O' + split_cur[1]
else:
    save_dir = args.output_dir
if save_dir[-1] != '/':
    save_dir = save_dir + '/'

'''
den_pert = sys.argv[3]
if '.0' in den_pert:
    den_pert = 'SingleStar'
'''
file_name = 'WIND_hdf5_plt_cnt_*'

#define measured quantities:
time = []
mass = []
maximum_speed = []
momentum = []
ang_momentum = []
'''
#create list of directories
dirs = commands.getoutput('ls ' + dir).split('\n')
lref = 0
del_dirs = []
for d in dirs:
    if den_pert not in d or 'lref' not in d[-7:]:
        del_dirs.append(d)
    else:
        time.append([])
        mass.append([])
        maximum_speed.append([])
        momentum.append([])
        ang_momentum.append([])

print "Found directories"

for d in del_dirs:
    dirs.remove(d)
dirs.sort(key=natural_keys)
print "loaded dirs"

#define time selection:
dt = 20.0

#define analysis cylinder:
radius = 500.0
#region:
upper_cylinder_bound = 250.0
lower_cylinder_bound = -250.0
height = 500.0
'''
#define analysis cylinder:
radius = 500.0
#region:
upper_cylinder_bound = 250.0
lower_cylinder_bound = -250.0
height = 500.0

#define constants
ts = yt.DatasetSeries(dir +file_name)
ts = ts[1::10]
domain = ts[0].domain_width[0]
AU = ts[0].unit_registry['au'][0]
Year = ts[0].unit_registry['yr'][0]

#define prevs
prev_mass = [0.0, 0.0]
prev_maximum_speed = [0.0, 0.0]
prev_momentum = [0.0, 0.0]
prev_ang_momentum = 0.0
prev_time = np.inf
'''
dit = 0
for d in dirs:

ts = DatasetSeries(dir + file_name)
'''
#read in sink creation time
files = commands.getoutput('ls ' + dir + 'WIND_proj*').split('\n')
sink_form = 0.0
for file in files:
    f = h5py.File(file, 'r')
    if len(f.keys()) > 12:
        sink_form = f['time'][0]/Year
        print "found sink formation as time", sink_form
        break
'''
#find times:
start_year = 0.0
max_time = 3100.0
t_prev = start_year
ts_usable = []
files = commands.getoutput('ls ' + dir + d + '/' + file_name).split('\n')
for it in range(len(files)):
    f = h5py.File(files[it], 'r')
    time_val = f['real scalars'][0][1]/Year - sink_form
    if time_val -t_prev > dt and t_prev < max_time:
        t_prev = t_prev + dt
        ts_usable.append(ts[it])
    f.close()
print "Found usable files"
'''
#open files to read data out to:
if os.path.isfile(save_dir + output_file) == False:
    f = open(save_dir + output_file, 'a+')
    f.close()
f = open(save_dir + output_file, 'w')
f.write('Lref ' + dir.split('_')[-1] + ': Time, Mass, Momentum, Angular Momentum, Max speed \n')
f.close()

file = 0
for pf in ts:
    time_val = pf['time']/Year - sink_form
    
    if time_val > 0.0:
        #define cylinders:
        tube_1 = pf.disk([0.0, 0.0, (upper_cylinder_bound+(height/2.))*AU], [0.0, 0.0, 1.0], (radius, 'au'), (height/2., 'au'))
        tube_2 = pf.disk([0.0, 0.0, (lower_cylinder_bound-(height/2.))*AU], [0.0, 0.0, -1.0], (radius, 'au'), (height/2., 'au'))

        #calculate mass in cylinders
        # for region 1
        pos_pos = np.where(tube_1['velocity_z'] > 0.0)[0]
        if len(pos_pos) == 0:
            mass_1 = 0.0
            speed_1 = 0.0
            max_speed_1 = 0.0
            mom_1 = 0.0
        else:
            mass_1 = tube_1['cell_mass'][pos_pos].in_units('Msun').value
            speed_1 = tube_1['velocity_magnitude'][pos_pos].in_units("km/s").value
            max_speed_1 = np.max(speed_1)
            mom_1 = speed_1*mass_1
        
        #for region 2
        neg_pos = np.where(tube_2['velocity_z'] < 0.0)[0]
        if len(neg_pos) == 0:
            mass_2 = 0.0
            speed_2 = 0.0
            max_speed_2 = 0.0
            mom_2 = 0.0
        else:
            mass_2 = tube_2['cell_mass'][neg_pos].in_units('Msun').value
            speed_2 = tube_2['velocity_magnitude'][neg_pos].in_units("km/s").value
            max_speed_2 = np.max(speed_2)
            mom_2 = speed_2*mass_2
        
        #Calculate outflow mass:
        outflow_mass = np.sum(mass_1) + np.sum(mass_2)
        
        #Calculate max outflow speed:
        if max_speed_1 > max_speed_2:
            max_speed = abs(max_speed_1)
        else:
            max_speed = abs(max_speed_2)

        #Calculate outflow momentum:
        mom = np.sum(mom_1) + np.sum(mom_2)

        #Calculate outflow angular momentum:
        L_x = np.sum(tube_1['angular_momentum_x'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['angular_momentum_x'][neg_pos].in_units("Msun * km**2 / s").value)
        L_y = np.sum(tube_1['angular_momentum_y'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['angular_momentum_y'][neg_pos].in_units("Msun * km**2 / s").value)
        L_z = np.sum(tube_1['angular_momentum_z'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['angular_momentum_z'][neg_pos].in_units("Msun * km**2 / s").value)
        L = np.sqrt((L_x)**2. + (L_y)**2. + (L_z)**2.)
        
        #Append quantities
        if outflow_mass == 0.0:
            outflow_mass = np.nan
            max_speed = np.nan
            mom = np.nan
            L = np.nan
            
        time.append(time_val)
        mass.append(outflow_mass)
        maximum_speed.append(max_speed)
        momentum.append(mom)
        ang_momentum.append(L)
        print "OUTPUT=", time_val, outflow_mass, mom, L, max_speed
        f = open(save_dir + output_file, 'a')
        f.write(str(time_val) + ',' + str(outflow_mass) + ',' + str(mom) + ',' + str(L) + ',' + str(max_speed) + '\n')
        f.close()
        prev_time = time_val

        file = file + 1
        print "Progress:", file, "/", len(ts)
#dit = dit + 1
'''
time = np.array(time)
mass = np.array(mass).clip(min=0.000001)
momentum = np.array(momentum).clip(min=0.000001)
ang_momentum = np.array(ang_momentum).clip(min=0.000001)
maximum_speed = np.array(maximum_speed).clip(min=0.000001)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)
ls = ['-', '--', '-.', ':', '+', 'x']


#for m in range(len(mass)):
ax1.semilogy(time, mass)
#ax1.plot(time, mass)
ax1.set_xlabel('time since protostar formation [yr]')
ax1.set_ylabel('outflow mass [M$_\odot$]')
ax1.set_xlim([0, 3000])
ax1.set_ylim([1.e-5, 1.0])

#for p in range(len(momentum)):
ax2.semilogy(time, momentum)
#ax2.plot(time, momentum)
ax2.set_xlabel('time since protostar formation [yr]')
ax2.set_ylabel('outflow momentum[M$_\odot$kms$^{-1}$]')
ax2.set_xlim([0, 3000])
ax2.set_ylim([1.e-5, 1.0])

#for l in range(len(ang_momentum)):
ax3.semilogy(time, ang_momentum)
#ax3.plot(time, ang_momentum)
ax3.set_xlabel('time since protostar formation [yr]')
ax3.set_ylabel('outflow ang. mom. [M$_\odot$km$^2$s$^{-1}$]')
ax3.set_xlim([0, 3000])
ax3.set_ylim([1.e5, 1.e10])

#for s in range(len(maximum_speed)):
ax4.semilogy(time, maximum_speed, ls)
#ax4.plot(time, maximum_speed)
ax4.set_xlabel('time since protostar formation [yr]')
ax4.set_ylabel('maximum outflow speed [kms$^{-1}$]')
ax4.set_xlim([0, 3000])
ax4.set_ylim([0.1, 1.e4])

plt.tight_layout()

filename = 'Mass_outflow'
plt.savefig(save_dir + filename + ".eps", format='eps', bbox_inches='tight')
'''
