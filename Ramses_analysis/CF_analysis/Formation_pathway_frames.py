#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
import my_ramses_fields as myf
import csv

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-22')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-19)
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xy")
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-G_mass", "--Gas_mass", type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

args = parse_inputs()

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
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
mym.set_units(units_override)

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
usuable_file_inds.append(usuable_file_inds[-1]-1)
usuable_files = np.array(files)[usuable_file_inds]
center_sink = Other_sink[0]

thickness = yt.YTQuantity(5000, 'au')
x = np.linspace(-1*thickness/2, thickness/2, 800)
y = np.linspace(-1*thickness/2, thickness/2, 800)
X, Y = np.meshgrid(x, y)
annotate_space = (-1*thickness/2 - thickness/2)/31
x_ind = []
y_ind = []
counter = 0
while counter < 31:
    val = annotate_space*counter + annotate_space/2. - thickness/2
    x_ind.append(int(val))
    y_ind.append(int(val))
    counter = counter + 1
X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
xlim = [-2500, 2500]
ylim = [-2500, 2500]
x_width = (xlim[1] -xlim[0])
y_width = (ylim[1] -ylim[0])

cbar_max = args.colourbar_max
try:
    cbar_min = float(args.colourbar_min)
except:
    cbar_min = float(args.colourbar_min[1:])

prev_center = np.nan
sink_creation_time = np.nan
pickle_file_preffix = 'bound_core_frag_'
pit = 4
for usuable_file in usuable_files:
    pit = pit -1
    pickle_file = pickle_file_preffix + str(pit) + '.pkl'
    ds = yt.load(usuable_file, units_override=units_override)
    dd = ds.all_data()

    try:
        center_pos = yt.YTArray([dd['sink_particle_posx'][center_sink], dd['sink_particle_posy'][center_sink], dd['sink_particle_posz'][center_sink]]).in_units('au')
    except:
        center_pos = prev_center
    if np.isnan(prev_center):
        prev_center = center_pos
        
    if np.isnan(sink_creation_time):
        sink_creation_time = dd['sink_particle_form_time'][Interested_sinks[0]]
        time_val = ds.current_time.value*scale_t.in_units('yr') - sink_creation_time
        time_val = np.round(time_val)
        
    if args.axis == 'xy':
        axis_ind = 2
        left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
    elif args.axis == 'xz':
        axis_ind = 1
        left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.5*thickness), center_pos[2]-(0.75*thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.5*thickness), center_pos[2]+(0.55*thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
    elif args.axis == 'yz':
        axis_ind = 0
        left_corner = yt.YTArray([center_pos[0]-(0.5*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.75*thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.5*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.75*thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
        
    TM = np.sum(region['cell_mass'].in_units('g'))
    x_top = np.sum(region['cell_mass'].in_units('g')*region['x-velocity'].in_units('cm/s'))
    y_top = np.sum(region['cell_mass'].in_units('g')*region['y-velocity'].in_units('cm/s'))
    z_top = np.sum(region['cell_mass'].in_units('g')*region['z-velocity'].in_units('cm/s'))
    com_vel = [(x_top/TM), (y_top/TM), (z_top/TM)]
    center_vel = yt.YTArray(com_vel, 'cm')
        
    X_image = X
    Y_image = Y
    X_image_vel = X_vel
    Y_image_vel = Y_vel
    
    part_info = mym.get_particle_data(ds, axis=args.axis, sink_id=center_sink, region=region)
    
    if args.axis == 'xy':
        part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0].value
        part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[1].value
    elif args.axis == 'xz':
        part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0].value
        part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2].value
    elif args.axis == 'yz':
        part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[1].value
        part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2].value
            
    vel1_field = args.axis[0] + '-velocity'
    vel2_field = args.axis[1] + '-velocity'
    mag1_field = 'mag' + args.axis[0]
    mag2_field = 'mag' + args.axis[1]
    proj_dict = {'density':[], vel1_field:[], vel2_field:[], mag1_field:[], mag2_field:[]}
    proj_dict_keys = str(proj_dict.keys()).split("['")[1].split("']")[0].split("', '")
    proj_field_list =[('gas', 'density'), ('ramses', vel1_field), ('ramses', vel2_field), ('gas', mag1_field), ('gas', mag2_field)]
    proj_root_rank = int(rank/len(proj_field_list))*len(proj_field_list)
    
    for field in yt.parallel_objects(proj_field_list):
        proj = yt.ProjectionPlot(ds, axis_ind, field, width=(x_width,'au'), data_source=region, method='integrate', center=(center_pos, 'AU'))
        if 'mag' in str(field):
            if args.axis == 'xz':
                proj_array = np.array(proj.frb.data[field].T.in_units('cm*gauss')/thickness.in_units('cm'))
            else:
                proj_array = np.array(proj.frb.data[field].in_units('cm*gauss')/thickness.in_units('cm'))
        elif field == ('gas', 'density'):
            if args.axis == 'xz':
                proj_array = np.array((proj.frb.data[field].T/thickness.in_units('cm')).in_units("g/cm**3"))
            else:
                proj_array = np.array((proj.frb.data[field]/thickness.in_units('cm')).in_units("g/cm**3"))
        else:
            if args.axis == 'xz':
                proj_array = np.array(proj.frb.data[field].T.in_cgs()/thickness.in_units('cm'))
            else:
                proj_array = np.array(proj.frb.data[field].in_cgs()/thickness.in_units('cm'))
        if rank == proj_root_rank:
            proj_dict[field[1]] = proj_array
        else:
            file = open(pickle_file.split('.pkl')[0] + '_proj_data_' + str(proj_root_rank)+ str(proj_dict_keys.index(field[1])) + '.pkl', 'wb')
            pickle.dump((field[1], proj_array), file)
            file.close()
        
    if rank == proj_root_rank and size > 1:
        for kit in range(1,len(proj_dict_keys)):
            file = open(pickle_file.split('.pkl')[0] + '_proj_data_' +str(proj_root_rank) +str(kit)+'.pkl', 'rb')
            key, proj_array = pickle.load(file)
            file.close()
            proj_dict[key] = proj_array
            os.remove(pickle_file.split('.pkl')[0] + '_proj_data_' +str(proj_root_rank) +str(kit)+'.pkl')
            
    if rank == proj_root_rank:
        image = proj_dict[proj_dict_keys[0]]
        velx_full = proj_dict[proj_dict_keys[1]]
        vely_full = proj_dict[proj_dict_keys[2]]
        magx = proj_dict[proj_dict_keys[3]]
        magy = proj_dict[proj_dict_keys[4]]
    
    if rank == proj_root_rank:
        velx, vely, velz = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=center_vel, axis=args.axis)
        del velx_full
        del vely_full

        args_dict = {}
        if args.annotate_time == "True":
            args_dict.update({'annotate_time': r"$t$="+str(int(time_val))+"yr"})
        args_dict.update({'field':'density'})
        args_dict.update({'annotate_velocity': True})
        args_dict.update({'time_val': time_val})
        args_dict.update({'cbar_min': cbar_min})
        args_dict.update({'cbar_max': cbar_max})
        args_dict.update({'xabel': "X (AU)"})
        args_dict.update({'yabel': "Y (AU)"})
        args_dict.update({'axlim':args.ax_lim})
        args_dict.update({'xlim':xlim})
        args_dict.update({'ylim':ylim})
        args_dict.update({'has_particles':has_particles})

        if args.absolute_image != "False":
            image = abs(image)
        file = open(pickle_file, 'wb')
        pickle.dump((X_image, Y_image, image, magx, magy, X_image_vel, Y_image_vel, velx, vely, part_info, args_dict, []), file)
        file.close()
        print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
        del image
        del magx
        del magy
        del velx
        del vely
        del args_dict
    del has_particles
    del time_val
    del center_vel
    del part_info
    del X_image
    del Y_image
    del X_image_vel
    del Y_image_vel
        
    print('FINISHED MAKING YT PROJECTIONS ON RANK', rank)

    sys.stdout.flush()
    CW.Barrier()
    import pdb
    pdb.set_trace()
   



#Delay core frag pathway
Primary_form_time = 1.0387929956526736
Secondary_form_time = 1.040190745650386
System_bound_time = 1.0405948456497247

m_times = [Primary_form_time, Secondary_form_time, System_bound_time]

#Dynamical capture
Star_form_time = 1.0358556956574807
Capture_time = 1.036077995657117

m_times = [Star_form_time, Capture_time]
