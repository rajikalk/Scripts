#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import my_ramses_module as mym
import my_ramses_fields as myf
import csv
import gc

#=======MAIN=======
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = '100'

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

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
mym.set_units(units_override)

sys.stdout.flush()
CW.Barrier()

#------------------------------
Sim_path = '/lustre/astro/troels/IMF_256_fixed_dt/data/'
files = sorted(glob.glob(Sim_path+"*/info*.txt"))
txt_files = sorted(glob.glob(Sim_path+"*/stars_output.snktxt"))
sim_file_times = []

for output_txt in txt_files:
    with open(output_txt, 'rU') as txt_file:
        reader = csv.reader(txt_file)
        for row in reader:
            time_val = float(row[0].split('   ')[-2])
            sim_file_times.append(time_val)
            break

del txt_files
gc.collect()

sys.stdout.flush()
CW.Barrier()

Interested_sinks = [36, 14, 2]
Other_sink = [4, [10, [5, 9]], [1, 3]]

#----------------------------------------------------------------------
#Bound core fragmentation pathway
Primary_form_time = 1.0365265956563827
Secondary_form_time = 1.0460617956407776

m_times = [Secondary_form_time, Primary_form_time]
usuable_file_inds = []

for m_time in m_times:
    match_time_ind = np.argmin(abs(np.array(sim_file_times) - m_time))
    if sim_file_times[match_time_ind] < m_time:
        match_time_ind = match_time_ind + 1
    usuable_file_inds.append(match_time_ind)

usuable_file_inds = [16, 5, 4]
#usuable_file_inds.append(usuable_file_inds[-1]-1)
usable_files = np.array(files)[usuable_file_inds]
center_sink = Other_sink[0]
del usuable_file_inds
gc.collect()

sys.stdout.flush()
CW.Barrier()

from pyramses import rsink
thickness = yt.YTQuantity(5000, 'au')
center_positions = []
pit = 4
for fn_it in range(len(usable_files)):
    pit = pit - 1
    pickle_file = "movie_frame_" + ("%06d" % pit) + "_part.pkl"
    fn = usable_files[fn_it]
    file_no = int(fn.split('output_')[-1].split('/')[0])
    datadir = fn.split('output_')[0]
    loaded_sink_data = rsink(file_no, datadir=datadir)
    try:
        center_pos = yt.YTArray([loaded_sink_data['x'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['y'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['z'][center_sink]*units['length_unit'].in_units('au')])
        import pdb
        pdb.set_trace()
        sink_creation_time = loaded_sink_data['tcreate'][center_sink]*units['time_unit'].in_units('yr')
        center_positions.append(center_pos)
    except:
        center_pos = center_positions[-1]
        center_positions.append(center_pos)
    x_lim = [center_pos[0] - thickness/2, center_pos[0] + thickness/2]
    y_lim = [center_pos[1] - thickness/2, center_pos[1] + thickness/2]
    z_lim = [center_pos[2] - thickness/2, center_pos[2] + thickness/2]
    sinks_in_box = np.where((loaded_sink_data['x']*units['length_unit'].in_units('au')>x_lim[0])&(loaded_sink_data['x']*units['length_unit'].in_units('au')<x_lim[1])&(loaded_sink_data['y']*units['length_unit'].in_units('au')>y_lim[0])&(loaded_sink_data['y']*units['length_unit'].in_units('au')<y_lim[1])&(loaded_sink_data['z']*units['length_unit'].in_units('au')>z_lim[0])&(loaded_sink_data['z']*units['length_unit'].in_units('au')<z_lim[1]))[0]
    particle_masses = loaded_sink_data['m'][sinks_in_box]*units['mass_unit']
    particle_x_pos = loaded_sink_data['x'][sinks_in_box]*units['length_unit']
    particle_y_pos = loaded_sink_data['y'][sinks_in_box]*units['length_unit']
    gc.collect()
    #particle_masses = dd['sink_particle_mass']

    #if np.remainder(rank,48) == 0:
    file = open(pickle_file, 'wb')
    #pickle.dump((image, time_val, particle_positions, particle_masses), file)
    pickle.dump((particle_x_pos, particle_y_pos, particle_masses), file)
    file.close()
    print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
    del x_lim
    del y_lim
    del z_lim
    del particle_masses
    del particle_x_pos
    del particle_y_pos
    gc.collect()

prev_center = np.nan
sink_creation_time = np.nan
pickle_file_preffix = 'bound_core_frag_'
pit = 4

sys.stdout.flush()
CW.Barrier()
cit = -1
for usuable_file in usable_files:
    pit = pit - 1
    cit = cit + 1
    ds = yt.load(usuable_file, units_override=units_override)
    #dd = ds.all_data()

    center_pos = center_positions[cit]
    time_val = ds.current_time.in_units('yr') - sink_creation_time
    
    axis_ind = 2
    left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
    right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
    region = ds.box(left_corner, right_corner)
    del left_corner
    del right_corner
    gc.collect()
    
    proj = yt.ProjectionPlot(ds, axis_ind, ("ramses", "Density"), width=thickness, data_source=region, method='integrate', center=(center_pos, 'AU'))
    proj_array = np.array(proj.frb.data[("ramses", "Density")]/units['length_unit'].in_units('cm'))
    image = proj_array*units['density_unit'].in_units('g/cm**3')
    del proj
    del proj_array
    gc.collect()
    
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    file = open(pickle_file, 'wb')
    pickle.dump((image, time_val), file)
    file.close()
    print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
            
    print('FINISHED MAKING YT PROJECTIONS ON RANK', rank)

sys.stdout.flush()
CW.Barrier()
import pdb
pdb.set_trace()
   
#Make frames
pickle_files = sorted(glob.glob('bound_core_frag_*.pkl'))
pit = -1
for pickle_file in pickle_files:
    pit = pit + 1
    file = open(pickle_file, 'rb')
    image, part_info, time_val = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()

    file_name = save_dir + "movie_frame_" + ("%06d" % pit)
    
    xlim = [-2500, 2500]
    ylim = [-2500, 2500]
    x = np.linspace(xlim[0], xlim[1], 800)
    y = np.linspace(ylim[0], ylim[1], 800)
    X, Y = np.meshgrid(x, y)
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

    
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
    plt.gca().set_aspect('equal')
    #plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
    cbar = plt.colorbar(plot, pad=0.0)
    #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)

    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
    
    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)

    plt.tick_params(axis='both', which='major')# labelsize=16)
    for line in ax.xaxis.get_ticklines():
        line.set_color('white')
    for line in ax.yaxis.get_ticklines():
        line.set_color('white')

    try:
        plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
        try:
            plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        except:
            print("Couldn't outline time string")
    except:
        print("Couldn't plot time string")
                
    

    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
    print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')

#Delay core frag pathway
Primary_form_time = 1.0387929956526736
Secondary_form_time = 1.040190745650386
System_bound_time = 1.0405948456497247

m_times = [Primary_form_time, Secondary_form_time, System_bound_time]

#Dynamical capture
Star_form_time = 1.0358556956574807
Capture_time = 1.036077995657117

m_times = [Star_form_time, Capture_time]
