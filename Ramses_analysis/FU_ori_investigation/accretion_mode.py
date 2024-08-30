#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import my_ramses_fields as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def has_sinks(ds):
    '''
    Checks particle file to see if particles exists, or tries the plot file.
    '''
    dd = ds.all_data()
    if len(dd['sink_particle_tag'][myf.get_centred_sink_id():].astype(int)) != 0:
        del dd
        return True
    else:
        del dd
        return False

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

#File files
usable_files = sorted(glob.glob(input_dir+"*/info*.txt"))

sys.stdout.flush()
CW.Barrier()

#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
mym.set_units(units_override)


#find sink particle to center on and formation time
del units_override['density_unit']
ds = yt.load(usable_files[-1], units_override=units_override)
#try:
dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
del dd
    
no_data_points = len(usable_files)
    
sys.stdout.flush()
CW.Barrier()

if args.make_frames_only == 'False':
    #Trying yt parallelism
    file_int = -1
    for fn in yt.parallel_objects(usable_files, njobs=int(size/6)):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
            if usable_files[file_int] == usable_files[file_int-1]:
                os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        make_pickle = False
        if args.plot_time is None:
            pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl"
        else:
            pickle_file = save_dir + "time_" + str(args.plot_time) +".pkl"
        if os.path.isfile(pickle_file) == False:
            make_pickle = True
        elif os.path.isfile(pickle_file) == True:
            if os.stat(pickle_file).st_size == 0:
                make_pickle = True
        if usable_files[file_int] == usable_files[file_int-1]:
            os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
            print("copied", save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl", "to",  save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl")
        if make_pickle == True:
            
            ds = yt.load(fn, units_override=units_override)
            dd = ds.all_data()
            has_particles = has_sinks(ds)
            
            #Define box:
            center_pos = dd['Center_Position'].in_units('au').value
            if args.axis == 'xy':
                axis_ind = 2
                left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[1]-(0.75*y_width), center_pos[2]-(0.5*args.slice_thickness)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[1]+(0.75*y_width), center_pos[2]+(0.5*args.slice_thickness)], 'AU')
                region = ds.box(left_corner, right_corner)
                del left_corner
                del right_corner
            elif args.axis == 'xz':
                axis_ind = 1
                left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[1]-(0.5*args.slice_thickness), center_pos[2]-(0.75*y_width)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[1]+(0.5*args.slice_thickness), center_pos[2]+(0.55*y_width)], 'AU')
                region = ds.box(left_corner, right_corner)

                del left_corner
                del right_corner
            elif args.axis == 'yz':
                axis_ind = 0
                left_corner = yt.YTArray([center_pos[0]-(0.5*args.slice_thickness), center_pos[1]-(0.75*x_width), center_pos[2]-(0.75*y_width)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+(0.5*args.slice_thickness), center_pos[1]+(0.75*x_width), center_pos[2]+(0.75*y_width)], 'AU')
                region = ds.box(left_corner, right_corner)
                del left_corner
                del right_corner
            
            if has_particles:
                part_info = mym.get_particle_data(ds, axis=args.axis, sink_id=sink_id, region=region)
            else:
                part_info = {}
            
            try:
                time_val = m_times[file_int]
            except:
                sink_creation_time = np.min(dd['particle_creation_time'].value)
                time_real = yt.YTQuantity(ds.current_time.value - sink_creation_time, 's')
                time_val = np.round(time_real.in_units('yr'))
                del sink_creation_time
                del time_real
                
            if args.use_angular_momentum != 'False':
                if len(part_info['particle_mass']) == 1:
                    left_corner_test = yt.YTArray([dd['sink_particle_posx'][sink_id].in_units('AU').value - 100, dd['sink_particle_posy'][sink_id].in_units('AU').value - 100, dd['sink_particle_posz'][sink_id].in_units('AU').value - 100], 'AU')
                    right_corner_test = yt.YTArray([dd['sink_particle_posx'][sink_id].in_units('AU').value + 100, dd['sink_particle_posy'][sink_id].in_units('AU').value + 100, dd['sink_particle_posz'][sink_id].in_units('AU').value + 100], 'AU')
                    region = ds.box(left_corner_test, right_corner_test)
                    L_x = np.sum(region['Angular_Momentum_x'].value)
                    L_y = np.sum(region['Angular_Momentum_y'].value)
                    L_z = np.sum(region['Angular_Momentum_z'].value)
                    L = np.array([L_x, L_y, L_z])/np.sum(region['Angular_Momentum'].value)
                    myf.set_normal(L)
                    print("L =", L)
                    del left_corner_test
                    del right_corner_test
                    del region
                    del L_x
                    del L_y
                    del L_z
                else:
                    L_x = np.sum(dd['Orbital_Angular_Momentum_x'].value)
                    L_y = np.sum(dd['Orbital_Angular_Momentum_y'].value)
                    L_z = np.sum(dd['Orbital_Angular_Momentum_z'].value)
                    L = np.array([L_x, L_y, L_z])/np.sum(dd['Orbital_Angular_Momentum'].value)
                    myf.set_normal(L)
                    print("L =", L)
            
            #mass_array = dd['sink_particle_mass']

            #print('Center Pos=' + str(center_pos))
            
            #Update X and Y to be centered on center position
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
            
            #print("initialised fields")
            '''
            if args.axis == 'xy':
                axis_ind = 2
                left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[1]-(0.75*y_width), center_pos[2]-(0.5*args.slice_thickness)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[1]+(0.75*y_width), center_pos[2]+(0.5*args.slice_thickness)], 'AU')
                region = ds.box(left_corner, right_corner)
                del left_corner
                del right_corner
            elif args.axis == 'xz':
                axis_ind = 1
                left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[1]-(0.5*args.slice_thickness), center_pos[2]-(0.75*y_width)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[1]+(0.5*args.slice_thickness), center_pos[2]+(0.55*y_width)], 'AU')
                region = ds.box(left_corner, right_corner)
                del left_corner
                del right_corner
            elif args.axis == 'yz':
                axis_ind = 0
                left_corner = yt.YTArray([center_pos[0]-(0.5*args.slice_thickness), center_pos[1]-(0.75*x_width), center_pos[2]-(0.75*y_width)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+(0.5*args.slice_thickness), center_pos[1]+(0.75*x_width), center_pos[2]+(0.75*y_width)], 'AU')
                region = ds.box(left_corner, right_corner)
                del left_corner
                del right_corner
            '''
            if args.use_angular_momentum != 'False':
                region = yt.disk(center_pos, L, (np.sqrt((0.5*x_width)**2 + (0.5*y_width)), 'AU'), (args.slice_thickness/2, 'AU'))
            weight_field = args.weight_field
            
            myf.set_center_pos_ind(args.calculation_center)

            if args.image_center == 0 and args.use_gas_center_calc=='False':
                #If not using the gas the calculate the Center posiiton and velocity in the fields
                TM = np.sum(region['cell_mass'].in_units('g'))
                x_top = np.sum(region['cell_mass'].in_units('g')*region['x-velocity'].in_units('cm/s'))
                y_top = np.sum(region['cell_mass'].in_units('g')*region['y-velocity'].in_units('cm/s'))
                z_top = np.sum(region['cell_mass'].in_units('g')*region['z-velocity'].in_units('cm/s'))
                com_vel = [(x_top/TM), (y_top/TM), (z_top/TM)]
                center_vel = yt.YTArray(com_vel, 'cm')
                del TM
                del x_top
                del y_top
                del z_top
                #del com_vel
            else:
                center_vel = region['Center_Velocity'].in_units('cm/s').value
            #myf.set_center_pos_ind(args.image_center)
            #print("center_vel =", center_vel, "on rank", rank, "for", ds)
            
            if args.axis == 'xy':
                center_vel_plane = np.array([center_vel[0], center_vel[1]])
                perp_vel = 'z'
            elif args.axis == 'xz':
                center_vel_plane = np.array([center_vel[0], center_vel[2]])
                perp_vel = 'y'
            elif args.axis == 'yz':
                center_vel_plane = np.array([center_vel[1], center_vel[2]])
                perp_vel = 'x'
            
            if args.use_angular_momentum == 'False':
                vel1_field = args.axis[0] + '-velocity'
                vel2_field = args.axis[1] + '-velocity'
                vel3_field = perp_vel + '-velocity'
                mag1_field = 'mag' + args.axis[0]
                mag2_field = 'mag' + args.axis[1]
                proj_dict = {simfo['field'][1]:[], vel1_field:[], vel2_field:[], vel3_field:[], mag1_field:[], mag2_field:[]}
                proj_dict_keys = str(proj_dict.keys()).split("['")[1].split("']")[0].split("', '")
                proj_field_list =[simfo['field'], ('ramses', vel1_field), ('ramses', vel2_field), ('ramses', vel3_field), ('gas', mag1_field), ('gas', mag2_field)]
                proj_root_rank = int(rank/len(proj_field_list))*len(proj_field_list)
                
                proj_dict = {}
                for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
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
                    sto.result_id = field[1]
                    sto.result = proj_array
                    '''
                    if rank == proj_root_rank:
                        proj_dict[field[1]] = proj_array
                    else:
                        file = open(pickle_file.split('.pkl')[0] + '_proj_data_' + str(proj_root_rank)+ str(proj_dict_keys.index(field[1])) + '.pkl', 'wb')
                        pickle.dump((field[1], proj_array), file)
                        file.close()
                    '''
                '''
                if rank == proj_root_rank and size > 1:
                    for kit in range(1,len(proj_dict_keys)):
                        file = open(pickle_file.split('.pkl')[0] + '_proj_data_' +str(proj_root_rank) +str(kit)+'.pkl', 'rb')
                        key, proj_array = pickle.load(file)
                        file.close()
                        proj_dict[key] = proj_array
                        os.remove(pickle_file.split('.pkl')[0] + '_proj_data_' +str(proj_root_rank) +str(kit)+'.pkl')
                '''
                if rank == proj_root_rank:
                    image = proj_dict[proj_dict_keys[0]]
                    velx_full = proj_dict[proj_dict_keys[1]]
                    vely_full = proj_dict[proj_dict_keys[2]]
                    velz_full = proj_dict[proj_dict_keys[3]]
                    magx = proj_dict[proj_dict_keys[4]]
                    magy = proj_dict[proj_dict_keys[5]]
                        
            elif args.use_angular_momentum != 'False':
                proj_root_rank = int(rank/7)*7
                #proj_dict = {simfo['field'][1]:[]}
                proj_dict = {simfo['field'][1]:[], 'Projected_Velocity_x':[], 'Projected_Velocity_y':[], 'Projected_Velocity_z':[], 'Projected_Magnetic_Field_x':[], 'Projected_Magnetic_Field_y':[], 'Projected_Magnetic_Field_z':[]}
                proj_dict_keys = str(proj_dict.keys()).split("['")[1].split("']")[0].split("', '")
                #proj_field_list =[simfo['field']]
                proj_field_list =[simfo['field'], ('gas', 'Projected_Velocity_x'), ('gas', 'Projected_Velocity_y'), ('gas', 'Projected_Velocity_z'), ('gas', 'Projected_Magnetic_Field_x'), ('gas', 'Projected_Magnetic_Field_y'), ('gas', 'Projected_Magnetic_Field_z')]
                
                for field in yt.parallel_objects(proj_field_list):
                    proj = yt.OffAxisProjectionPlot(ds, L, field, width=(x_width/2, 'AU'), weight_field=weight_field, method='integrate', center=(center_pos, 'AU'), depth=(args.slice_thickness, 'AU'))
                    if 'mag' in str(field):
                        if weight_field == None:
                            proj_array = np.array(proj.frb.data[field].in_units('cm*gauss')/thickness.in_units('cm'))
                        else:
                            proj_array = np.array(proj.frb.data[field].in_units('gauss'))
                    else:
                       if weight_field == None:
                            proj_array = np.array(proj.frb.data[field].in_cgs()/thickness.in_units('cm'))
                       else:
                           proj_array = np.array(proj.frb.data[field].in_cgs())
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
                
                #Figure out vectors projections onto axes
                y_axis_vector = proj.data_source.orienter.north_vector
                projected_velocity = yt.YTArray([proj_dict['Projected_Velocity_x'], proj_dict['Projected_Velocity_y'], proj_dict['Projected_Velocity_z']], 'cm/s')
                #projected_velocity = yt.YTArray([proj_dict['x-velocity'], proj_dict['y-velocity'], proj_dict['z-velocity']], 'cm/s')
                proj_x_vel_1 = (np.dot(projected_velocity.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[0]
                proj_x_vel_2 = (np.dot(projected_velocity.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[1]
                proj_x_vel_3 = (np.dot(projected_velocity.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[2]
                proj_x_vel = yt.YTArray([proj_x_vel_1.T, proj_x_vel_2.T, proj_x_vel_3.T], 'cm/s')
                proj_x_vel_mag = np.sqrt(proj_x_vel_1**2 + proj_x_vel_2**2 + proj_x_vel_3**2)
                proj_y_vel = (projected_velocity - proj_x_vel)
                proj_y_vel_mag = np.sqrt(proj_y_vel[0]**2 + proj_y_vel[1]**2 + proj_y_vel[2]**2)
                #center_vel = [np.mean(proj_x_vel_mag), np.mean(proj_y_vel_mag)]
                center_vel = [0.0, 0.0]
                
                projected_B = yt.YTArray([proj_dict['Projected_Magnetic_Field_x'], proj_dict['Projected_Magnetic_Field_y'], proj_dict['Projected_Magnetic_Field_z']], 'G')
                #projected_B = yt.YTArray([proj_dict['magx'], proj_dict['magy'], proj_dict['magz']], 'G')
                proj_x_B_1 = (np.dot(projected_B.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[0]
                proj_x_B_2 = (np.dot(projected_B.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[1]
                proj_x_B_3 = (np.dot(projected_B.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[2]
                proj_x_B = yt.YTArray([proj_x_B_1.T, proj_x_B_2.T, proj_x_B_3.T], 'G')
                proj_x_B_mag = np.sqrt(proj_x_B_1**2 + proj_x_B_2**2 + proj_x_B_3**2)
                proj_y_B = (projected_B - proj_x_B)
                proj_y_B_mag = np.sqrt(proj_y_B[0]**2 + proj_y_B[1]**2 + proj_y_B[2]**2)

                #Update particle data
                y_axis_vector = proj.data_source.orienter.north_vector
                dd = ds.all_data()
                projected_position = yt.YTArray([dd['Projected_Particle_Posx'].in_units('AU').value, dd['Projected_Particle_Posy'].in_units('AU').value, dd['Projected_Particle_Posz'].in_units('AU').value], 'AU')
                proj_x_pos_1 = (np.dot(projected_position.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[0]
                proj_x_pos_2 = (np.dot(projected_position.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[1]
                proj_x_pos_3 = (np.dot(projected_position.T, y_axis_vector)/np.dot(y_axis_vector,y_axis_vector))*y_axis_vector[2]
                proj_x_pos = yt.YTArray([proj_x_pos_1.T, proj_x_pos_2.T, proj_x_pos_3.T], 'G')
                proj_x_pos_mag = np.sqrt(proj_x_pos_1**2 + proj_x_pos_2**2 + proj_x_pos_3**2)
                proj_y_pos = (projected_position - proj_x_pos)
                proj_y_pos_mag = np.sqrt(proj_y_pos[0]**2 + proj_y_pos[1]**2 + proj_y_pos[2]**2)
                positions = np.array([proj_x_pos_mag.value,proj_y_pos_mag.value])
                part_info['particle_position'] = positions
                
                if rank == proj_root_rank:
                    image = proj_dict[proj_dict_keys[0]]
                    velx_full = proj_x_vel_mag
                    vely_full = proj_y_vel_mag
                    magx = proj_x_B_mag.value
                    magy = proj_y_B_mag.value
                    part_info['particle_position'] = positions
            
            if rank == proj_root_rank:
                velx, vely, velz = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=center_vel, velz_full=velz_full, axis=args.axis)
                del velx_full
                del vely_full
                del velz_full

                args_dict = {}
                if args.annotate_time == "True":
                    args_dict.update({'annotate_time': r"$t$="+str(int(time_val))+"yr"})
                args_dict.update({'field':simfo['field']})
                args_dict.update({'annotate_velocity': args.plot_velocity_legend})
                args_dict.update({'time_val': time_val})
                args_dict.update({'cbar_min': cbar_min})
                args_dict.update({'cbar_max': cbar_max})
                args_dict.update({'title': title})
                args_dict.update({'xabel': xabel})
                args_dict.update({'yabel': yabel})
                args_dict.update({'axlim':args.ax_lim})
                args_dict.update({'xlim':xlim})
                args_dict.update({'ylim':ylim})
                args_dict.update({'has_particles':has_particles})

                if args.absolute_image != "False":
                    image = abs(image)
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, image, magx, magy, X_image_vel, Y_image_vel, velx, vely, velz, part_info, args_dict, simfo), file)
                file.close()
                print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
                del image
                del magx
                del magy
                del velx
                del vely
                del velz
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

