#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import csv
import gc
import collections


def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

#=======MAIN=======
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(3000,"Msun")}
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

sys.stdout.flush()
CW.Barrier()

#-------------------------------------------------------
Sim_path = '/lustre/astro/troels/IMF_256_fixed_dt/data/'
txt_files = sorted(glob.glob(Sim_path+"*/stars_output.snktxt"))
sim_file_times = []

for output_txt in txt_files:
    with open(output_txt, 'r') as txt_file:
        reader = csv.reader(txt_file)
        for row in reader:
            time_val = float(row[0].split('   ')[-2])
            sim_file_times.append(time_val)
            break

gc.collect()
dt_min = np.min((np.array(sim_file_times[1:]) - np.array(sim_file_times[:-1])))*units['time_unit'].in_units('yr')

sys.stdout.flush()
CW.Barrier()


#-------------------------------------
#Find system candidates:

birth_con_pickle = "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl"

file = open(birth_con_pickle, 'rb')
Sink_birth_all = pickle.load(file)
file.close()

Bound_core_frag_candidates = []
Unbound_core_frag_candidates = []
Dynamical_capture_candidates = []

sink_id = 0
while sink_id < len(Sink_birth_all.keys())-1:
    sink_id = sink_id + 1
    if Sink_birth_all[str(sink_id)][2] != 'nan':
        if Sink_birth_all[str(sink_id)][0] == True:
            if '[' not in Sink_birth_all[str(sink_id)][2]:
                Bound_core_frag_candidates.append((sink_id, Sink_birth_all[str(sink_id)][1]))
        else:
            if Sink_birth_all[str(sink_id)][1] == Sink_birth_all[str(sink_id)][2] and Sink_birth_all[str(sink_id)][-2] > dt_min:
                Unbound_core_frag_candidates.append((sink_id, Sink_birth_all[str(sink_id)][1]))
            elif Sink_birth_all[str(sink_id)][1] not in flatten(eval(Sink_birth_all[str(sink_id)][2])) and Sink_birth_all[str(sink_id)][-2] > dt_min:
                Dynamical_capture_candidates.append((sink_id, (Sink_birth_all[str(sink_id)][1], Sink_birth_all[str(sink_id)][2])))

global_pickle = '/groups/astro/rlk/rlk/Global_sink_pickles/G100_full.pkl'
file = open(global_pickle, 'rb')
global_data = pickle.load(file)
file.close()

rm_pair = []
for pair in Bound_core_frag_candidates:
    center_sink = pair[0]
    form_ind = np.where(global_data['m'].T[center_sink]>0)[0][0]
    secondary_form_time = global_data['time'].T[center_sink][form_ind]
    unbound_sink = int(pair[1])
    form_ind = np.where(global_data['m'].T[unbound_sink]>0)[0][0]
    primary_form_time = global_data['time'].T[unbound_sink][form_ind]
    dt = secondary_form_time - primary_form_time
    if dt < dt_min:
        rm_pair.append(pair)

import pdb
pdb.set_trace()
Bound_core_frag_candidates = list(set(Bound_core_frag_candidates).symmetric_difference(set(rm_pair)))

rm_pair = []
for pair in Unbound_core_frag_candidates:
    center_sink = pair[0]
    form_ind = np.where(global_data['m'].T[center_sink]>0)[0][0]
    if '[' in pair[1]:
        unbound_sink = flatten(eval(pair[1]))[np.argmax(global_data['m'][form_ind][flatten(eval(pair[1]))])]
    else:
        unbound_sink = int(pair[1])
    form_pos = np.array([global_data['x'].T[center_sink][form_ind], global_data['y'].T[center_sink][form_ind], global_data['z'].T[center_sink][form_ind]])*units['length_unit'].in_units('au')
    unbound_sink_pos = np.array([global_data['x'].T[unbound_sink][form_ind], global_data['y'].T[unbound_sink][form_ind], global_data['z'].T[unbound_sink][form_ind]])*units['length_unit'].in_units('au')
    d_pos = abs(form_pos-unbound_sink_pos)
    if True in (d_pos>10000):
        rm_pair.append(pair)

Unbound_core_frag_candidates = list(set(Unbound_core_frag_candidates).symmetric_difference(set(rm_pair)))

rm_pair = []
for pair in Dynamical_capture_candidates:
    center_sink = pair[0]
    unbound_sink = pair[1][0]
    form_ind = np.where(global_data['m'].T[center_sink]>0)[0][0]
    form_pos = np.array([global_data['x'].T[center_sink][form_ind], global_data['y'].T[center_sink][form_ind], global_data['z'].T[center_sink][form_ind]])*units['length_unit'].in_units('au')
    unbound_sink_pos = np.array([global_data['x'].T[unbound_sink][form_ind], global_data['y'].T[unbound_sink][form_ind], global_data['z'].T[unbound_sink][form_ind]])*units['length_unit'].in_units('au')
    d_pos = abs(form_pos-unbound_sink_pos)
    if True in (d_pos>10000):
        rm_pair.append(pair)
        
Dynamical_capture_candidates = list(set(Dynamical_capture_candidates).symmetric_difference(set(rm_pair)))

#Set systems to use:
Bound_core_frag_system = Bound_core_frag_candidates[2]
Unbound_core_frag_system = Unbound_core_frag_candidates[0]
Dynamical_capture_system = Dynamical_capture_candidates[0]

#Bound core fragmentation
Bound_primary_form_time = global_data['time'].T[Bound_core_frag_system[1]][np.where(global_data['m'].T[Bound_core_frag_system[1]]>0)[0][0]]
Bound_secondary_form_time = global_data['time'].T[Bound_core_frag_system[0]][np.where(global_data['m'].T[Bound_core_frag_system[0]]>0)[0][0]]
Bound_m_times = [Bound_secondary_form_time, Bound_primary_form_time]

Unbound_secondary_form_time = global_data['time'].T[Unbound_core_frag_system[0]][np.where(global_data['m'].T[Unbound_core_frag_system[0]]>0)[0][0]]
Unbound_bound_time = (Unbound_secondary_form_time*units['time_unit'].in_units('yr').value + Sink_birth_all[str(Unbound_core_frag_system[0])][-2])/units['time_unit'].in_units('yr').value
Unbound_m_times = [Unbound_bound_time, Unbound_secondary_form_time]

Dynamical_secondary_form_time = global_data['time'].T[Dynamical_capture_system[0]][np.where(global_data['m'].T[Dynamical_capture_system[0]]>0)[0][0]]
Dynamical_bound_time = (Dynamical_secondary_form_time*units['time_unit'].in_units('yr').value + Sink_birth_all[str(Dynamical_capture_system[0])][-2])/units['time_unit'].in_units('yr').value
Dynamical_m_times = [Dynamical_bound_time, Dynamical_secondary_form_time]

del Sink_birth_all
del global_data
gc.collect()
#----------------------------------------------------------------------
#Bound core fragmentation pathway
usable_files = []

for m_time in Bound_m_times:
    match_time_ind = np.argmin(abs(np.array(sim_file_times) - m_time))
    if sim_file_times[match_time_ind] < m_time:
        match_time_ind = match_time_ind + 1
    #use string manipulation to get the relative info file
    star_file = txt_files[match_time_ind]
    info_file = star_file.split('stars_output.snktxt')[0] + 'info*.txt'
    usable_files.append(glob.glob(info_file)[0])
    
#Add presink frame
match_time_ind = match_time_ind - 1
star_file = txt_files[match_time_ind]
info_file = star_file.split('stars_output.snktxt')[0] + 'info*.txt'
usable_files.append(glob.glob(info_file)[0])
    
center_sink = Bound_core_frag_system[1]
gc.collect()

sys.stdout.flush()
CW.Barrier()

from pyramses import rsink
center_positions = []
pickle_file_preffix = 'bound_core_frag_'
pit = 4
Core_frag_sinks = sorted(flatten(Bound_core_frag_system))
max_seps = []
for fn_it in range(len(usable_files)):
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
    fn = usable_files[fn_it]
    file_no = int(fn.split('output_')[-1].split('/')[0])
    datadir = fn.split('output_')[0]
    loaded_sink_data = rsink(file_no, datadir=datadir)
    try:
        center_pos = yt.YTArray([loaded_sink_data['x'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['y'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['z'][center_sink]*units['length_unit'].in_units('au')])
        sink_creation_time = loaded_sink_data['tcreate'][center_sink]*units['time_unit'].in_units('yr')
        center_positions.append(center_pos)
    except:
        center_pos = center_positions[-1]
        center_positions.append(center_pos)
    #x_lim = [center_pos[0] - thickness/2, center_pos[0] + thickness/2]
    #y_lim = [center_pos[1] - thickness/2, center_pos[1] + thickness/2]
    #z_lim = [center_pos[2] - thickness/2, center_pos[2] + thickness/2]
    #sinks_in_box = np.where((loaded_sink_data['x']*units['length_unit'].in_units('au')>x_lim[0])&(loaded_sink_data['x']*units['length_unit'].in_units('au')<x_lim[1])&(loaded_sink_data['y']*units['length_unit'].in_units('au')>y_lim[0])&(loaded_sink_data['y']*units['length_unit'].in_units('au')<y_lim[1])&(loaded_sink_data['z']*units['length_unit'].in_units('au')>z_lim[0])&(loaded_sink_data['z']*units['length_unit'].in_units('au')<z_lim[1]))[0]
    if len(loaded_sink_data['m'])>Core_frag_sinks[-1]:
        particle_masses = loaded_sink_data['m'][Core_frag_sinks]*units['mass_unit'].in_units('Msun')
        particle_x_pos = loaded_sink_data['x'][Core_frag_sinks]*units['length_unit'].in_units('au')
        particle_y_pos = loaded_sink_data['y'][Core_frag_sinks]*units['length_unit'].in_units('au')
    elif len(loaded_sink_data['m'])>Core_frag_sinks[0]:
        particle_masses = loaded_sink_data['m'][Core_frag_sinks[0]]*units['mass_unit'].in_units('Msun')
        particle_x_pos = loaded_sink_data['x'][Core_frag_sinks[0]]*units['length_unit'].in_units('au')
        particle_y_pos = loaded_sink_data['y'][Core_frag_sinks[0]]*units['length_unit'].in_units('au')
    else:
        particle_masses = yt.YTArray([], 'Msun')
        particle_x_pos = yt.YTArray([], 'au')
        particle_y_pos = yt.YTArray([], 'au')
    try:
        dx = np.max(abs(particle_x_pos-particle_x_pos[0]))
        dy = np.max(abs(particle_y_pos-particle_y_pos[0]))
        if dx > dy:
            max_seps.append(dx)
        else:
            max_seps.append(dy)
    except:
        pass
    gc.collect()
    #particle_masses = dd['sink_particle_mass']

    #if np.remainder(rank,48) == 0:
    file = open(pickle_file, 'wb')
    #pickle.dump((image, time_val, particle_positions, particle_masses), file)
    pickle.dump((particle_x_pos, particle_y_pos, particle_masses), file)
    file.close()
    print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
    #del x_lim
    #del y_lim
    #del z_lim
    gc.collect()

max_sep = np.max(max_seps)
thickness = yt.YTQuantity(np.ceil(max_sep/100)*100+500, 'au')

#del units
gc.collect()
pit = 4

sys.stdout.flush()
CW.Barrier()
cit = -1
import os
for usuable_file in usable_files:
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    if os.path.exists(pickle_file) == False:
        cit = cit + 1
        ds = yt.load(usuable_file, units_override=units_override)
        #dd = ds.all_data()

        center_pos = center_positions[cit]
        time_val = ds.current_time.in_units('yr') - sink_creation_time
        
        left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
        gc.collect()
        
        axis_ind = 2
        proj = yt.ProjectionPlot(ds, axis_ind, ("ramses", "Density"), width=thickness, data_source=region, method='integrate', center=(center_pos, 'AU'))
        proj_array = np.array(proj.frb.data[("ramses", "Density")])/thickness.in_units('cm')
        image = proj_array*units['density_unit'].in_units('g/cm**3')
        del proj
        del proj_array
        gc.collect()
        
        file = open(pickle_file, 'wb')
        pickle.dump((image, time_val), file)
        file.close()
        print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)

sys.stdout.flush()
CW.Barrier()
   
#Make frames
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
pickle_files = sorted(glob.glob('bound_core_frag_*_part.pkl'))
pit = -1
cit = 0
for pickle_file in pickle_files:
    pit = pit + 1
    cit = cit + 1
    center_pos = center_positions[-1*cit]
    
    file = open(pickle_file, 'rb')
    particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()
    
    file = open("".join(pickle_file.split('_part')), 'rb')
    image, time_val = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()

    file_name = pickle_file_preffix + '_' + ("%06d" % pit)
    
    xlim = [-1*thickness, thickness]
    ylim = [-1*thickness, thickness]
    x = np.linspace(xlim[0], xlim[1], 800)
    y = np.linspace(ylim[0], ylim[1], 800)
    X, Y = np.meshgrid(x, y)
    
    X = X + center_pos[0]
    Y = Y + center_pos[1]
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    
    annotate_space = (xlim[1] - xlim[0])/31
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 31:
        val = annotate_space*counter + annotate_space/2. + xlim[0]
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
          
    has_particles = True
    xabel = "X (AU)"
    yabel = "Y (AU)"
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xlabel(xabel, labelpad=-1, fontsize=10)
    ax.set_ylabel(yabel, fontsize=10) #, labelpad=-20
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=1.e-19, vmax=1.e-16), rasterized=True)
    plt.gca().set_aspect('equal')
    #plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
    cbar = plt.colorbar(plot, pad=0.0)
    #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)
    #ax.scatter(particle_x_pos, particle_y_pos, color='c', s=1)
    try:
        mym.annotate_particles(ax, np.array([particle_x_pos, particle_y_pos]), 200, limits=[xlim, ylim], annotate_field=particle_masses, particle_tags=Core_frag_sinks)
    except:
        try:
            mym.annotate_particles(ax, np.array([[particle_x_pos], [particle_y_pos]]), 200, limits=[xlim, ylim], annotate_field=[particle_masses], particle_tags=Core_frag_sinks)
        except:
            pass
    
    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)

    plt.tick_params(axis='both', which='major')# labelsize=16)
    for line in ax.xaxis.get_ticklines():
        line.set_color('white')
    for line in ax.yaxis.get_ticklines():
        line.set_color('white')

    try:
        plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
        try:
            plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        except:
            print("Couldn't outline time string")
    except:
        print("Couldn't plot time string")

    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
    print('Created frame ' + file_name + '.png')

#------------------------------------------------------------------------------------------------------------------------------------------------
#Unbound core frag
usable_files = []

for m_time in Unbound_m_times:
    match_time_ind = np.argmin(abs(np.array(sim_file_times) - m_time))
    if sim_file_times[match_time_ind] < m_time:
        match_time_ind = match_time_ind + 1
    #use string manipulation to get the relative info file
    star_file = txt_files[match_time_ind]
    info_file = star_file.split('stars_output.snktxt')[0] + 'info*.txt'
    usable_files.append(glob.glob(info_file)[0])
    
#Add presink frame
match_time_ind = match_time_ind - 1
star_file = txt_files[match_time_ind]
info_file = star_file.split('stars_output.snktxt')[0] + 'info*.txt'
usable_files.append(glob.glob(info_file)[0])
    
gc.collect()

sys.stdout.flush()
CW.Barrier()

from pyramses import rsink
max_seps = []
center_positions = []
pickle_file_preffix = 'unbound_core_frag_'
pit = 4
Unbound_core_frag_system = [Unbound_core_frag_system[0], eval(Unbound_core_frag_system[1])]
center_sink = Unbound_core_frag_system[1]
Core_frag_sinks = sorted(flatten(Unbound_core_frag_system))
for fn_it in range(len(usable_files)):
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
    fn = usable_files[fn_it]
    file_no = int(fn.split('output_')[-1].split('/')[0])
    datadir = fn.split('output_')[0]
    loaded_sink_data = rsink(file_no, datadir=datadir)
    try:
        center_pos = yt.YTArray([loaded_sink_data['x'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['y'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['z'][center_sink]*units['length_unit'].in_units('au')])
        sink_creation_time = loaded_sink_data['tcreate'][center_sink]*units['time_unit'].in_units('yr')
        center_positions.append(center_pos)
    except:
        center_pos = center_positions[-1]
        center_positions.append(center_pos)
    #x_lim = [center_pos[0] - thickness/2, center_pos[0] + thickness/2]
    #y_lim = [center_pos[1] - thickness/2, center_pos[1] + thickness/2]
    #z_lim = [center_pos[2] - thickness/2, center_pos[2] + thickness/2]
    #sinks_in_box = np.where((loaded_sink_data['x']*units['length_unit'].in_units('au')>x_lim[0])&(loaded_sink_data['x']*units['length_unit'].in_units('au')<x_lim[1])&(loaded_sink_data['y']*units['length_unit'].in_units('au')>y_lim[0])&(loaded_sink_data['y']*units['length_unit'].in_units('au')<y_lim[1])&(loaded_sink_data['z']*units['length_unit'].in_units('au')>z_lim[0])&(loaded_sink_data['z']*units['length_unit'].in_units('au')<z_lim[1]))[0]
    if len(loaded_sink_data['m'])>Core_frag_sinks[-1]:
        particle_masses = loaded_sink_data['m'][Core_frag_sinks]*units['mass_unit'].in_units('Msun')
        particle_x_pos = loaded_sink_data['x'][Core_frag_sinks]*units['length_unit'].in_units('au')
        particle_y_pos = loaded_sink_data['y'][Core_frag_sinks]*units['length_unit'].in_units('au')
    elif len(loaded_sink_data['m'])>Core_frag_sinks[0]:
        particle_masses = loaded_sink_data['m'][Core_frag_sinks[0]]*units['mass_unit'].in_units('Msun')
        particle_x_pos = loaded_sink_data['x'][Core_frag_sinks[0]]*units['length_unit'].in_units('au')
        particle_y_pos = loaded_sink_data['y'][Core_frag_sinks[0]]*units['length_unit'].in_units('au')
    else:
        particle_masses = yt.YTArray([], 'Msun')
        particle_x_pos = yt.YTArray([], 'au')
        particle_y_pos = yt.YTArray([], 'au')
    try:
        dx = np.max(abs(particle_x_pos-particle_x_pos[0]))
        dy = np.max(abs(particle_y_pos-particle_y_pos[0]))
        if dx > dy:
            max_seps.append(dx)
        else:
            max_seps.append(dy)
    except:
        pass
    gc.collect()
    #particle_masses = dd['sink_particle_mass']

    #if np.remainder(rank,48) == 0:
    file = open(pickle_file, 'wb')
    #pickle.dump((image, time_val, particle_positions, particle_masses), file)
    pickle.dump((particle_x_pos, particle_y_pos, particle_masses), file)
    file.close()
    print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
    #del x_lim
    #del y_lim
    #del z_lim
    del particle_masses
    del particle_x_pos
    del particle_y_pos
    gc.collect()

max_sep = np.max(max_seps)
thickness = yt.YTQuantity(np.ceil(max_sep/100)*100+500, 'au')

#del units
gc.collect()
pit = 4

sys.stdout.flush()
CW.Barrier()
cit = -1
import os
for usuable_file in usable_files:
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    if os.path.exists(pickle_file) == False:
        cit = cit + 1
        ds = yt.load(usuable_file, units_override=units_override)
        #dd = ds.all_data()

        center_pos = center_positions[cit]
        time_val = ds.current_time.in_units('yr') - sink_creation_time
        
        left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
        gc.collect()
        
        axis_ind = 2
        proj = yt.ProjectionPlot(ds, axis_ind, ("ramses", "Density"), width=thickness, data_source=region, method='integrate', center=(center_pos, 'AU'))
        proj_array = np.array(proj.frb.data[("ramses", "Density")])/thickness.in_units('cm')
        image = proj_array*units['density_unit'].in_units('g/cm**3')
        del proj
        del proj_array
        gc.collect()
        
        file = open(pickle_file, 'wb')
        pickle.dump((image, time_val), file)
        file.close()
        print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)

sys.stdout.flush()
CW.Barrier()
   
#Make frames
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
pickle_files = sorted(glob.glob('unbound_core_frag_*_part.pkl'))
pit = -1
cit = 0
for pickle_file in pickle_files:
    pit = pit + 1
    cit = cit + 1
    center_pos = center_positions[-1*cit]
    
    file = open(pickle_file, 'rb')
    particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()
    
    file = open("".join(pickle_file.split('_part')), 'rb')
    image, time_val = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()

    file_name = pickle_file_preffix + '_' + ("%06d" % pit)
    
    xlim = [-1*thickness, thickness]
    ylim = [-1*thickness, thickness]
    x = np.linspace(xlim[0], xlim[1], 800)
    y = np.linspace(ylim[0], ylim[1], 800)
    X, Y = np.meshgrid(x, y)
    
    X = X + center_pos[0]
    Y = Y + center_pos[1]
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    
    annotate_space = (xlim[1] - xlim[0])/31
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 31:
        val = annotate_space*counter + annotate_space/2. + xlim[0]
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
          
    has_particles = True
    xabel = "X (AU)"
    yabel = "Y (AU)"
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xlabel(xabel, labelpad=-1, fontsize=10)
    ax.set_ylabel(yabel, fontsize=10) #, labelpad=-20
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=np.min(image), vmax=np.max(image)), rasterized=True)
    plt.gca().set_aspect('equal')
    #plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
    cbar = plt.colorbar(plot, pad=0.0)
    #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)
    #ax.scatter(particle_x_pos, particle_y_pos, color='c', s=1)
    try:
        mym.annotate_particles(ax, np.array([particle_x_pos, particle_y_pos]), 200, limits=[xlim, ylim], annotate_field=particle_masses, particle_tags=Core_frag_sinks)
    except:
        try:
            mym.annotate_particles(ax, np.array([[particle_x_pos], [particle_y_pos]]), 200, limits=[xlim, ylim], annotate_field=[particle_masses], particle_tags=Core_frag_sinks)
        except:
            pass
    
    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)

    plt.tick_params(axis='both', which='major')# labelsize=16)
    for line in ax.xaxis.get_ticklines():
        line.set_color('white')
    for line in ax.yaxis.get_ticklines():
        line.set_color('white')

    try:
        plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
        try:
            plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        except:
            print("Couldn't outline time string")
    except:
        print("Couldn't plot time string")
                
    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
    print('Created frame ' + file_name + '.png')

#-------------------------------------------------------------------------------------------------------------------------------------------

#Dynamical capture
usable_files = []


for m_time in Dynamical_m_times:
    match_time_ind = np.argmin(abs(np.array(sim_file_times) - m_time))
    if sim_file_times[match_time_ind] < m_time:
        match_time_ind = match_time_ind + 1
    #use string manipulation to get the relative info file
    star_file = txt_files[match_time_ind]
    info_file = star_file.split('stars_output.snktxt')[0] + 'info*.txt'
    usable_files.append(glob.glob(info_file)[0])
    
#Add presink frame
match_time_ind = match_time_ind - 1
star_file = txt_files[match_time_ind]
info_file = star_file.split('stars_output.snktxt')[0] + 'info*.txt'
usable_files.append(glob.glob(info_file)[0])
    
center_sink = Dynamical_capture_system[0]
gc.collect()

sys.stdout.flush()
CW.Barrier()

from pyramses import rsink
max_seps = []
center_positions = []
pickle_file_preffix = 'dynamical_capture_'
pit = 4
Core_frag_sinks = sorted(flatten(Dynamical_capture_system))
for fn_it in range(len(usable_files)):
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
    fn = usable_files[fn_it]
    file_no = int(fn.split('output_')[-1].split('/')[0])
    datadir = fn.split('output_')[0]
    loaded_sink_data = rsink(file_no, datadir=datadir)
    try:
        center_pos = yt.YTArray([loaded_sink_data['x'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['y'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['z'][center_sink]*units['length_unit'].in_units('au')])
        sink_creation_time = loaded_sink_data['tcreate'][center_sink]*units['time_unit'].in_units('yr')
        center_positions.append(center_pos)
    except:
        center_pos = center_positions[-1]
        center_positions.append(center_pos)
    #x_lim = [center_pos[0] - thickness/2, center_pos[0] + thickness/2]
    #y_lim = [center_pos[1] - thickness/2, center_pos[1] + thickness/2]
    #z_lim = [center_pos[2] - thickness/2, center_pos[2] + thickness/2]
    #sinks_in_box = np.where((loaded_sink_data['x']*units['length_unit'].in_units('au')>x_lim[0])&(loaded_sink_data['x']*units['length_unit'].in_units('au')<x_lim[1])&(loaded_sink_data['y']*units['length_unit'].in_units('au')>y_lim[0])&(loaded_sink_data['y']*units['length_unit'].in_units('au')<y_lim[1])&(loaded_sink_data['z']*units['length_unit'].in_units('au')>z_lim[0])&(loaded_sink_data['z']*units['length_unit'].in_units('au')<z_lim[1]))[0]
    if len(loaded_sink_data['m'])>Core_frag_sinks[-1]:
        particle_masses = loaded_sink_data['m'][Core_frag_sinks]*units['mass_unit'].in_units('Msun')
        particle_x_pos = loaded_sink_data['x'][Core_frag_sinks]*units['length_unit'].in_units('au')
        particle_y_pos = loaded_sink_data['y'][Core_frag_sinks]*units['length_unit'].in_units('au')
    elif len(loaded_sink_data['m'])>Core_frag_sinks[0]:
        particle_masses = loaded_sink_data['m'][Core_frag_sinks[0]]*units['mass_unit'].in_units('Msun')
        particle_x_pos = loaded_sink_data['x'][Core_frag_sinks[0]]*units['length_unit'].in_units('au')
        particle_y_pos = loaded_sink_data['y'][Core_frag_sinks[0]]*units['length_unit'].in_units('au')
    else:
        particle_masses = yt.YTArray([], 'Msun')
        particle_x_pos = yt.YTArray([], 'au')
        particle_y_pos = yt.YTArray([], 'au')
    try:
        dx = np.max(abs(particle_x_pos-particle_x_pos[0]))
        dy = np.max(abs(particle_y_pos-particle_y_pos[0]))
        if dx > dy:
            max_seps.append(dx)
        else:
            max_seps.append(dy)
    except:
        pass
    gc.collect()
    #particle_masses = dd['sink_particle_mass']

    #if np.remainder(rank,48) == 0:
    file = open(pickle_file, 'wb')
    #pickle.dump((image, time_val, particle_positions, particle_masses), file)
    pickle.dump((particle_x_pos, particle_y_pos, particle_masses), file)
    file.close()
    print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
    #del x_lim
    #del y_lim
    #del z_lim
    del particle_masses
    del particle_x_pos
    del particle_y_pos
    gc.collect()

max_sep = np.max(max_seps)
thickness = yt.YTQuantity(np.ceil(max_sep/100)*100+500, 'au')

#del units
gc.collect()
pit = 4

sys.stdout.flush()
CW.Barrier()
cit = -1
import os
for usuable_file in usable_files:
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    if os.path.exists(pickle_file) == False:
        cit = cit + 1
        ds = yt.load(usuable_file, units_override=units_override)
        #dd = ds.all_data()

        center_pos = center_positions[cit]
        time_val = ds.current_time.in_units('yr') - sink_creation_time
        
        left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
        gc.collect()
        
        axis_ind = 2
        proj = yt.ProjectionPlot(ds, axis_ind, ("ramses", "Density"), width=thickness, data_source=region, method='integrate', center=(center_pos, 'AU'))
        proj_array = np.array(proj.frb.data[("ramses", "Density")])/thickness.in_units('cm')
        image = proj_array*units['density_unit'].in_units('g/cm**3')
        del proj
        del proj_array
        gc.collect()
        
        file = open(pickle_file, 'wb')
        pickle.dump((image, time_val), file)
        file.close()
        print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)

sys.stdout.flush()
CW.Barrier()
   
#Make frames
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
pickle_files = sorted(glob.glob('dynamical_capture_*_part.pkl'))
pit = -1
cit = 0
for pickle_file in pickle_files:
    pit = pit + 1
    cit = cit + 1
    center_pos = center_positions[-1*cit]
    
    file = open(pickle_file, 'rb')
    particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()
    
    file = open("".join(pickle_file.split('_part')), 'rb')
    image, time_val = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()

    file_name = pickle_file_preffix + '_' + ("%06d" % pit)
    
    xlim = [-1*thickness, thickness]
    ylim = [-1*thickness, thickness]
    x = np.linspace(xlim[0], xlim[1], 800)
    y = np.linspace(ylim[0], ylim[1], 800)
    X, Y = np.meshgrid(x, y)
    
    X = X + center_pos[0]
    Y = Y + center_pos[1]
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    
    annotate_space = (xlim[1] - xlim[0])/31
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 31:
        val = annotate_space*counter + annotate_space/2. + xlim[0]
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
          
    has_particles = True
    xabel = "X (AU)"
    yabel = "Y (AU)"
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xlabel(xabel, labelpad=-1, fontsize=10)
    ax.set_ylabel(yabel, fontsize=10) #, labelpad=-20
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=np.min(image), vmax=np.max(image)), rasterized=True)
    plt.gca().set_aspect('equal')
    #plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
    cbar = plt.colorbar(plot, pad=0.0)
    #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)
    #ax.scatter(particle_x_pos, particle_y_pos, color='c', s=1)
    try:
        mym.annotate_particles(ax, np.array([particle_x_pos, particle_y_pos]), 200, limits=[xlim, ylim], annotate_field=particle_masses, particle_tags=Core_frag_sinks)
    except:
        try:
            mym.annotate_particles(ax, np.array([[particle_x_pos], [particle_y_pos]]), 200, limits=[xlim, ylim], annotate_field=[particle_masses], particle_tags=Core_frag_sinks)
        except:
            pass
    
    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)

    plt.tick_params(axis='both', which='major')# labelsize=16)
    for line in ax.xaxis.get_ticklines():
        line.set_color('white')
    for line in ax.yaxis.get_ticklines():
        line.set_color('white')

    try:
        plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
        try:
            plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        except:
            print("Couldn't outline time string")
    except:
        print("Couldn't plot time string")
                
    

    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
    print('Created frame ' + file_name + '.png')
