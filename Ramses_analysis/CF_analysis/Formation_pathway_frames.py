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
import os

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

if rank == 0 and os.path.exists('candidates.pkl') == False:
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

    print("sorted all systems into formation pathways")

    global_pickle = '/groups/astro/rlk/rlk/High_cadence/G100.pkl'
    file = open(global_pickle, 'rb')
    global_data = pickle.load(file)
    file.close()
    
    print("read in global data")

    #rm_pair = []
    Bound_core_frag_candidates_reduced = []
    for pair in Bound_core_frag_candidates:
        center_sink = pair[0]
        form_ind = np.where(global_data['m'].T[center_sink]>0)[0][0]
        secondary_form_time = global_data['time'].T[center_sink][form_ind]
        unbound_sink = int(pair[1])
        form_ind = np.where(global_data['m'].T[unbound_sink]>0)[0][0]
        primary_form_time = global_data['time'].T[unbound_sink][form_ind]
        dt = (secondary_form_time - primary_form_time)*units['time_unit'].in_units('yr')
        if dt > dt_min:
            #if '[' in pair[1]:
            #    other_ind = np.max(flatten(eval(pair[1])))
            #else:
            #    other_ind = int(pair[1])
            #Save times formation
            other_ind = int(pair[1])
            Bound_primary_form_time = global_data['time'].T[other_ind][np.where(global_data['m'].T[other_ind]>0)[0][0]]
            Bound_secondary_form_time = global_data['time'].T[pair[0]][np.where(global_data['m'].T[pair[0]]>0)[0][0]]
            Bound_m_times = [Bound_secondary_form_time, Bound_primary_form_time]
            Bound_core_frag_candidates_reduced.append([pair, Bound_m_times])

    Bound_core_frag_candidates = Bound_core_frag_candidates_reduced
    del Bound_core_frag_candidates_reduced

    print('removed all bound core fragmentation candidates with formation times separated by < dt_min')

    #rm_pair = []
    Unbound_core_frag_candidates_reduced = []
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
        if False in (d_pos<10000):
            if '[' in pair[1]:
                other_ind = np.max(flatten(eval(pair[1])))
            else:
                other_ind = int(pair[1])
            Unbound_primary_form_time = global_data['time'].T[other_ind][np.where(global_data['m'].T[other_ind]>0)[0][0]]
            Unbound_secondary_form_time = global_data['time'].T[pair[0]][np.where(global_data['m'].T[pair[0]]>0)[0][0]]
            Unbound_m_times = [Unbound_secondary_form_time, Unbound_primary_form_time]
            Unbound_core_frag_candidates_reduced.append([pair, Unbound_m_times])
            

    Unbound_core_frag_candidates = Unbound_core_frag_candidates_reduced #list(set(Unbound_core_frag_candidates).symmetric_difference(set(rm_pair)))
    del Unbound_core_frag_candidates_reduced

    print('removed all unbound core fragmentation candidates with initial separation > 10000 au')

    #rm_pair = []
    Dynamical_capture_candidates_reduced = []
    for pair in Dynamical_capture_candidates:
        center_sink = pair[0]
        unbound_sink = pair[1][0]
        form_ind = np.where(global_data['m'].T[center_sink]>0)[0][0]
        form_pos = np.array([global_data['x'].T[center_sink][form_ind], global_data['y'].T[center_sink][form_ind], global_data['z'].T[center_sink][form_ind]])*units['length_unit'].in_units('au')
        unbound_sink_pos = np.array([global_data['x'].T[unbound_sink][form_ind], global_data['y'].T[unbound_sink][form_ind], global_data['z'].T[unbound_sink][form_ind]])*units['length_unit'].in_units('au')
        d_pos = np.sqrt(np.sum((form_pos-unbound_sink_pos)**2))
        if d_pos<20000:
            if '[' in pair[1][1]:
                other_ind = np.max(flatten(eval(pair[1][1])))
            else:
                other_ind = int(pair[1][1])
            #print('removing', pair, 'because d_min=', d_pos)
            #rm_pair.append(pair)
            Dynamical_secondary_form_time = global_data['time'].T[pair[0]][np.where(global_data['m'].T[pair[0]]>0)[0][0]]
            Dynamical_bound_time = (Dynamical_secondary_form_time*units['time_unit'].in_units('yr').value + Sink_birth_all[str(pair[0])][-2])/units['time_unit'].in_units('yr').value
            Dynamical_m_times = [Dynamical_bound_time, Dynamical_secondary_form_time]
            Dynamical_capture_candidates_reduced.append([pair, Dynamical_m_times])
            
    Dynamical_capture_candidates = Dynamical_capture_candidates_reduced#list(set(Dynamical_capture_candidates).symmetric_difference(set(rm_pair)))
    del Dynamical_capture_candidates_reduced
    
    print('removed dynamical capture candidates that had birth separations >20000 au')

    del Sink_birth_all
    del global_data
    gc.collect()
    
    #Save candidates in pickle
    candidate_pickles = 'candidates.pkl'
    file = open(candidate_pickles, 'wb')
    pickle.dump((Bound_core_frag_candidates, Unbound_core_frag_candidates, Dynamical_capture_candidates), file)
    file.close()

sys.stdout.flush()
CW.Barrier()


file = open('candidates.pkl', 'rb')
Bound_core_frag_candidates, Unbound_core_frag_candidates, Dynamical_capture_candidates = pickle.load(file)
file.close()

sys.stdout.flush()
CW.Barrier()

#investigate system with sink 88
for system in yt.parallel_objects(Bound_core_frag_candidates, njobs=int(size/(3))):# Bound_core_frag_candidates: #3 projections
    print('Processing system', system, 'on rank', rank)
    
    Bound_m_times = system[1]
    
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

    print("usable files for Bound core fragmentation are", usable_files)

    
    if type(system[0][1]) == str:
        if '[' in system[0][1]:
            import pdb
            pdb.set_trace()
        else:
            center_sink = int(system[0][1])
    else:
        center_sink = system[0][1]
    gc.collect()

    sys.stdout.flush()
    CW.Barrier()

    #Find Sink positions
    pickle_file_preffix = 'bound_core_frag_'+str(system[0]) + '_'
    
    from pyramses import rsink
    center_positions = []
    pit = 4
    try:
        try:
            system[0][1] = int(system[0][1])
        except:
            system[0][1] = flatten(eval(system[0][1]))
        Core_frag_sinks = [system[0][0]] + [system[0][1]]
    except:
        Core_frag_sinks = list(system[0])
    all_max_seps = []
    for fn in usable_files:#yt.parallel_objects(usable_files, njobs=int(3)): #range(len(usable_files)):
        pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
        max_seps = []
        if os.path.exists(pickle_file) == False:
            print('Getting sink positions from', fn, 'on rank', rank)
            pit = pit - 1
            #fn = usable_files[fn_it]
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
            existing_sinks = list(set(Core_frag_sinks).intersection(np.arange(len(loaded_sink_data['m']))))
            if len(existing_sinks)>0:
                particle_masses = loaded_sink_data['m'][existing_sinks]*units['mass_unit'].in_units('Msun')
                particle_x_pos = loaded_sink_data['x'][existing_sinks]*units['length_unit'].in_units('au')
                particle_y_pos = loaded_sink_data['y'][existing_sinks]*units['length_unit'].in_units('au')
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

            if np.remainder(rank, 3) == 0:
                #if np.remainder(rank,48) == 0:
                all_max_seps = all_max_seps + max_seps
                file = open(pickle_file, 'wb')
                #pickle.dump((image, time_val, particle_positions, particle_masses), file)
                pickle.dump((particle_x_pos, particle_y_pos, particle_masses, max_seps), file)
                file.close()
                print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
            #del x_lim
            #del y_lim
            #del z_lim
            gc.collect()
        else:
            file = open(pickle_file, 'rb')
            particle_x_pos, particle_y_pos, particle_masses, max_seps = pickle.load(file)
            file.close()
            all_max_seps = all_max_seps + max_seps
            
    max_sep = np.max(flatten(all_max_seps))
    thickness = yt.YTQuantity(np.ceil(max_sep/100)*100+500, 'au')

    #del units
    gc.collect()
    pit = 4

    sys.stdout.flush()
    CW.Barrier()
    cit = -1
    for usable in yt.parallel_objects(usable_files, njobs=int(3)):
        print('making projection of', usable, 'on rank', rank)
        #for usable in usable_files:
        pit = pit - 1
        pickle_file = pickle_file_preffix + str(pit) + '.pkl'
        if os.path.exists(pickle_file) == False:
            cit = cit + 1
            ds = yt.load(usable, units_override=units_override)
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
    pickle_files = sorted(glob.glob(pickle_file_preffix + '*_part.pkl'))
    #cit = 0
    #for pickle_file in pickle_files:
    for pickle_file in yt.parallel_objects(pickle_files, njobs=int(3)):
        pit = pickle_files.index(pickle_file)
        #cit = cit + 1
        center_pos = center_positions[::-1][pit]
        
        file = open(pickle_file, 'rb')
        particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()
        
        file = open("".join(pickle_file.split('_part')), 'rb')
        image, time_val = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()

        file_name = pickle_file_preffix + ("%06d" % pit)
        
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
        print('Created frame ' + file_name + '.png, on rank', rank)

sys.stdout.flush()
CW.Barrier()

for system in yt.parallel_objects(Unbound_core_frag_candidates, njobs=int(size/(3))):# Bound_core_frag_candidates: #3 projections
    #Bound core fragmentation
    print('Processing system', system)
    
    Unbound_m_times = system[1]
    
    #Bound core fragmentation pathway
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

    print("usable files for Unbound core fragmentation are", usable_files)

    if type(system[0][1]) == str:
        if '[' in system[0][1]:
            import pdb
            pdb.set_trace()
        else:
            center_sink = int(system[0][1])
    else:
        center_sink = system[0][1]
    gc.collect()
    gc.collect()

    sys.stdout.flush()
    CW.Barrier()

    #Find Sink positions
    from pyramses import rsink
    center_positions = []
    pickle_file_preffix = 'unbound_core_frag_'+str(system[0]) + '_'
    pit = 4
    try:
        try:
            system[0][1] = int(system[0][1])
        except:
            system[0][1] = flatten(eval(system[0][1]))
        Core_frag_sinks = [system[0][0]] + [system[0][1]]
    except:
        Core_frag_sinks = list(system[0])
    max_seps = []
    for fn in usable_files:#yt.parallel_objects(usable_files, njobs=int(3)): #range(len(usable_files)):
        print('Getting sink positions from', fn, 'on rank', rank)
        pit = pit - 1
        pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
        #fn = usable_files[fn_it]
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
        existing_sinks = list(set(Core_frag_sinks).intersection(np.arange(len(loaded_sink_data['m']))))
        if len(existing_sinks)>0:
            particle_masses = loaded_sink_data['m'][existing_sinks]*units['mass_unit'].in_units('Msun')
            particle_x_pos = loaded_sink_data['x'][existing_sinks]*units['length_unit'].in_units('au')
            particle_y_pos = loaded_sink_data['y'][existing_sinks]*units['length_unit'].in_units('au')
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

        if np.remainder(rank, 3) == 0:
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
    for usable in yt.parallel_objects(usable_files, njobs=int(3)):
        #for usable in usable_files:
        pit = pit - 1
        pickle_file = pickle_file_preffix + str(pit) + '.pkl'
        if os.path.exists(pickle_file) == False:
            cit = cit + 1
            ds = yt.load(usable, units_override=units_override)
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
    pickle_files = sorted(glob.glob(pickle_file_preffix + '*_part.pkl'))
    #cit = 0
    #for pickle_file in pickle_files:
    for pickle_file in yt.parallel_objects(pickle_files, njobs=int(3)):
        pit = pickle_files.index(pickle_file)
        #cit = cit + 1
        center_pos = center_positions[::-1][pit]
        
        file = open(pickle_file, 'rb')
        particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()
        
        file = open("".join(pickle_file.split('_part')), 'rb')
        image, time_val = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()

        file_name = pickle_file_preffix + ("%06d" % pit)
        
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


sys.stdout.flush()
CW.Barrier()

for system in yt.parallel_objects(Dynamical_capture_candidates, njobs=int(size/(3))):# Bound_core_frag_candidates: #3 projections
    #Bound core fragmentation
    print('Processing system', system)

    Dynamical_m_times = system[1]
    
    #Bound core fragmentation pathway
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

    print("usable files for Bound core fragmentation are", usable_files)

    if type(system[0][1]) == str:
        if '[' in system[0][1]:
            import pdb
            pdb.set_trace()
        else:
            center_sink = int(system[0][1])
    else:
        center_sink = system[0][1]
    gc.collect()

    sys.stdout.flush()
    CW.Barrier()

    #Find Sink positions
    from pyramses import rsink
    center_positions = []
    pickle_file_preffix = 'dynamical_capt_'+str(system[0]) + '_'
    pit = 4
    try:
        try:
            system[0][1] = int(system[0][1])
        except:
            system[0][1] = flatten(eval(system[0][1]))
        Core_frag_sinks = [system[0][0]] + [system[0][1]]
    except:
        Core_frag_sinks = list(system[0])
    max_seps = []
    for fn in usable_files:#yt.parallel_objects(usable_files, njobs=int(3)): #range(len(usable_files)):
        print('Getting sink positions from', fn, 'on rank', rank)
        pit = pit - 1
        pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
        #fn = usable_files[fn_it]
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
        existing_sinks = list(set(Core_frag_sinks).intersection(np.arange(len(loaded_sink_data['m']))))
        if len(existing_sinks)>0:
            particle_masses = loaded_sink_data['m'][existing_sinks]*units['mass_unit'].in_units('Msun')
            particle_x_pos = loaded_sink_data['x'][existing_sinks]*units['length_unit'].in_units('au')
            particle_y_pos = loaded_sink_data['y'][existing_sinks]*units['length_unit'].in_units('au')
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

        if np.remainder(rank, 3) == 0:
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
    for usable in yt.parallel_objects(usable_files, njobs=int(3)):
        #for usable in usable_files:
        pit = pit - 1
        pickle_file = pickle_file_preffix + str(pit) + '.pkl'
        if os.path.exists(pickle_file) == False:
            cit = cit + 1
            ds = yt.load(usable, units_override=units_override)
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
    pickle_files = sorted(glob.glob(pickle_file_preffix + '*_part.pkl'))
    #cit = 0
    #for pickle_file in pickle_files:
    for pickle_file in yt.parallel_objects(pickle_files, njobs=int(3)):
        pit = pickle_files.index(pickle_file)
        #cit = cit + 1
        center_pos = center_positions[::-1][pit]
        
        file = open(pickle_file, 'rb')
        particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()
        
        file = open("".join(pickle_file.split('_part')), 'rb')
        image, time_val = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()

        file_name = pickle_file_preffix + ("%06d" % pit)
        
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

#Unbound_secondary_form_time = global_data['time'].T[Unbound_core_frag_system[0]][np.where(global_data['m'].T[Unbound_core_frag_system[0]]>0)[0][0]]
#Unbound_bound_time = (Unbound_secondary_form_time*units['time_unit'].in_units('yr').value + Sink_birth_all[str(Unbound_core_frag_system[0])][-2])/units['time_unit'].in_units('yr').value
#Unbound_m_times = [Unbound_bound_time, Unbound_secondary_form_time]

#Dynamical_secondary_form_time = global_data['time'].T[Dynamical_capture_system[0]][np.where(global_data['m'].T[Dynamical_capture_system[0]]>0)[0][0]]
#Dynamical_bound_time = (Dynamical_secondary_form_time*units['time_unit'].in_units('yr').value + Sink_birth_all[str(Dynamical_capture_system[0])][-2])/units['time_unit'].in_units('yr').value
#Dynamical_m_times = [Dynamical_bound_time, Dynamical_secondary_form_time]

#print('found times for the plotting candidates')
'''
del Sink_birth_all
del global_data
gc.collect()
#----------------------------------------------------------------------

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
for usable in usable_files:
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    if os.path.exists(pickle_file) == False:
        cit = cit + 1
        ds = yt.load(usable, units_override=units_override)
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
for usable in usable_files:
    pit = pit - 1
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    if os.path.exists(pickle_file) == False:
        cit = cit + 1
        ds = yt.load(usable, units_override=units_override)
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
'''
