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
    parser.add_argument("-G_mass", "--Gas_mass", type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
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

region_thickness = 5000
for usuable_file in usuable_files:
    ds = yt.load(usuable_file, units_override=units_override)
    dd = ds.all_data()
    
    center_pos = yt.YTArray([dd['sink_particle_posx'][center_sink], dd['sink_particle_posy'][center_sink], dd['sink_particle_posz'][center_sink]]).in_units('au')
        
    if args.axis == 'xy':
        axis_ind = 2
        left_corner = yt.YTArray([center_pos[0]-(0.75*region_thickness), center_pos[1]-(0.75*region_thickness), center_pos[2]-(0.5*region_thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*region_thickness), center_pos[1]+(0.75*region_thickness), center_pos[2]+(0.5*region_thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
    elif args.axis == 'xz':
        axis_ind = 1
        left_corner = yt.YTArray([center_pos[0]-(0.75*region_thickness), center_pos[1]-(0.5*region_thickness), center_pos[2]-(0.75*region_thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.75*region_thickness), center_pos[1]+(0.5*region_thickness), center_pos[2]+(0.55*region_thickness)], 'AU')
        region = ds.box(left_corner, right_corner)

        del left_corner
        del right_corner
    elif args.axis == 'yz':
        axis_ind = 0
        left_corner = yt.YTArray([center_pos[0]-(0.5*region_thickness), center_pos[1]-(0.75*region_thickness), center_pos[2]-(0.75*region_thickness)], 'AU')
        right_corner = yt.YTArray([center_pos[0]+(0.5*region_thickness), center_pos[1]+(0.75*region_thickness), center_pos[2]+(0.75*region_thickness)], 'AU')
        region = ds.box(left_corner, right_corner)
        del left_corner
        del right_corner
        
    if args.update_ax_lim == 'True':
        if args.axis == 'xy':
            X_image = X + center_pos[0]
            Y_image = Y + center_pos[1]
            X_image_vel = X_vel + center_pos[0]
            Y_image_vel = Y_vel + center_pos[1]
        elif args.axis == 'xz':
            X_image = X + center_pos[0]
            Y_image = Y + center_pos[2]
            X_image_vel = X_vel + center_pos[0]
            Y_image_vel = Y_vel + center_pos[2]
        elif args.axis == 'yz':
            X_image = X + center_pos[1]
            Y_image = Y + center_pos[2]
            X_image_vel = X_vel + center_pos[1]
            Y_image_vel = Y_vel + center_pos[2]
    else:
        X_image = X
        Y_image = Y
        X_image_vel = X_vel
        Y_image_vel = Y_vel
    
    if args.update_ax_lim == 'False':
        if args.axis == 'xy':
            part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0]
            part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[1]
        elif args.axis == 'xz':
            part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0]
            part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2]
        elif args.axis == 'yz':
            part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[1]
            part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2]
            
    vel1_field = args.axis[0] + '-velocity'
    vel2_field = args.axis[1] + '-velocity'
    vel3_field = perp_vel + '-velocity'
    mag1_field = 'mag' + args.axis[0]
    mag2_field = 'mag' + args.axis[1]
    proj_dict = {simfo['field'][1]:[], vel1_field:[], vel2_field:[], vel3_field:[], mag1_field:[], mag2_field:[]}
    proj_dict_keys = str(proj_dict.keys()).split("['")[1].split("']")[0].split("', '")
    proj_field_list =[simfo['field'], ('ramses', vel1_field), ('ramses', vel2_field), ('ramses', vel3_field), ('gas', mag1_field), ('gas', mag2_field)]
    proj_root_rank = int(rank/len(proj_field_list))*len(proj_field_list)
    
    for field in yt.parallel_objects(proj_field_list):
        proj = yt.ProjectionPlot(ds, axis_ind, field, width=(x_width,'au'), weight_field=weight_field, data_source=region, method='integrate', center=(center_pos, 'AU'))
        proj.set_buff_size([args.resolution, args.resolution])
        if 'mag' in str(field):
            if weight_field == None:
                if args.axis == 'xz':
                    proj_array = np.array(proj.frb.data[field].T.in_units('cm*gauss')/thickness.in_units('cm'))
                else:
                    proj_array = np.array(proj.frb.data[field].in_units('cm*gauss')/thickness.in_units('cm'))
            else:
                if args.axis == 'xz':
                    proj_array = np.array(proj.frb.data[field].T.in_units('gauss'))
                else:
                    proj_array = np.array(proj.frb.data[field].in_units('gauss'))
        elif args.field in str(field):
            if weight_field == None:
                if args.axis == 'xz':
                    if args.divide_by_proj_thickness == "True":
                        proj_array = np.array((proj.frb.data[field].T/thickness.in_units('cm')).in_units(args.field_unit))
                    else:
                        proj_array = np.array(proj.frb.data[field].T.in_units(args.field_unit+"*cm"))
                else:
                    if args.divide_by_proj_thickness == "True":
                        proj_array = np.array((proj.frb.data[field]/thickness.in_units('cm')).in_units(args.field_unit))
                    else:
                        proj_array = np.array(proj.frb.data[field].in_units(args.field_unit+"*cm"))
            else:
                if args.axis == 'xz':
                    proj_array = np.array(proj.frb.data[field].T.in_units(args.field_unit))
                else:
                    proj_array = np.array(proj.frb.data[field].in_units(args.field_unit))
        else:
            if weight_field == None:
                if args.axis == 'xz':
                    proj_array = np.array(proj.frb.data[field].T.in_cgs()/thickness.in_units('cm'))
                else:
                    proj_array = np.array(proj.frb.data[field].in_cgs()/thickness.in_units('cm'))
            else:
                if args.axis == 'xz':
                    proj_array = np.array(proj.frb.data[field].T.in_cgs())
                else:
                    proj_array = np.array(proj.frb.data[field].in_cgs())
        if str(args.field) in field and 'velocity' in str(args.field):
            proj_array = proj_array + com_vel[-1].in_units(args.field_unit).value
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
        velz_full = proj_dict[proj_dict_keys[3]]
        magx = proj_dict[proj_dict_keys[4]]
        magy = proj_dict[proj_dict_keys[5]]
    
    part_info = mym.get_particle_data(ds, axis=args.axis, sink_id=sink_id, region=region)
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
