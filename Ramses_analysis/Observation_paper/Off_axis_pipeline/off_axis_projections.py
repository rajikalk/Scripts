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
    #Projection field information
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="Number_Density")
    parser.add_argument("-f_unit", "--field_unit", help="What units would you like to plot the field?", default="cm**-3")
    parser.add_argument("-wf", "--weight_field", help="Do you want to have a weighted projection plot?", type=str, default=None)
    
    #Projection properties
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 500AU", type=float, default=300.)
    parser.add_argument("-proj_sep", "--projected_separation", help="if you want to make a projection such that the separation is a particular ammount, what is that? (in AU)", type=float, default=200.0)
    parser.add_argument("-threshold", "--density_threshold", help="What number density threshold would you like to use?", type=float, default=0.0)
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-ic", "--image_center", help="where would you like to center the image?", type=int, default=0) #0 = center of mass, 1=primary companion, 2=secondary companion
    parser.add_argument("-use_gas", "--use_gas_center_calc", help="Do you want to use gas when calculating the center position and veloity?", type=str, default='True')
    parser.add_argument("-use_part_for_vel", "--use_particle_for_center_vel_calc", help="Do you want to use the particles to calculate center velocity?", type=str, default='True')
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=150)
    
    #projection processing
    parser.add_argument("-div_by_thickness", "--divide_by_proj_thickness", help="Would you like to divide the field by the thickness of the projection?", default="True", type=str)
    parser.add_argument("-res", "--resolution", help="define image resolution", default=800, type=int)
    
    #Time inputs
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 100., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default=0, type=int)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-spec_file", "--specific_file", help="Do you want to use a specific file", type=str, default=None)
    
    #simulation information
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="what is the density so that the units can be calculated correctly?", type=str, default='100')
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
        
def xy_rotation_matrix(theta):
    """
    creates rotation matrix for given angle theta along the xy plane
    [cos(theta), -sin(theta), 0]
    [sin(theta), cos(theta) , 0]
    [0         , 0          , 1]
    """
    rot = np.array([[np.cos(theta), -1*np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0,0,1]])
    return rot
    
def z_rotation_matrix(phi):
    """
    creates rotation matrix for given angle phi along the yz plane
    [1, 0       , 0        ]
    [0, cos(phi), -sin(phi)]
    [0, sin(phi),  cos(phi)]
    """
    rot = np.array([[1,0,0], [0, np.cos(phi), -1*np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    return rot
    
def projected_vector(vector, proj_vector):
    """
    Calculates the position of vector projected onto proj_vector
    """
    vector_units = vector.units
    proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
    proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
    proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
    proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
    return proj_v

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
files = sorted(glob.glob(input_dir+"*/info*.txt"))

sys.stdout.flush()
CW.Barrier()

#Define units to override:
#BEWARE THIS IS HARD CODED
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.simulation_density_id

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
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')
mym.set_units(units_override)

#find sink particle to center on and formation time
ds = yt.load(files[-1], units_override=units_override)
field_it = [i for i, v in enumerate(ds.derived_field_list) if v[1] == args.field][0]
field = ds.derived_field_list[field_it]

#Sets center for calculating center position and velocity
myf.set_center_pos_ind(args.image_center)
myf.set_center_vel_ind(0)

if args.use_particle_for_center_vel_calc == 'True':
    myf.set_com_vel_use_part(True)
    myf.set_com_vel_use_gas(False)
else:
    myf.set_com_vel_use_part(False)

if args.use_gas_center_calc == 'True':
    myf.set_com_pos_use_gas(True)
else:
    myf.set_com_pos_use_gas(False)
    
#Make sure to only use gas when calculating the center velocity

sys.stdout.flush()
CW.Barrier()

#Sets some projection parameters
thickness = yt.YTQuantity(args.slice_thickness, 'AU')
if args.weight_field == 'None':
    weight_field = None
else:
    weight_field = args.weight_field

       
sys.stdout.flush()
CW.Barrier()

#Find sink particle that your are interested in, and get the formation
dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
sink_form_companion = dd['sink_particle_form_time'][sink_id+1]#Assumes the companion has the sink id of the primary sink id +1
if args.start_frame == 0 and args.plot_time == None:
    args.start_frame = int((sink_form_companion - sink_form_time)/(args.time_step))+1
del dd

sys.stdout.flush()
CW.Barrier()


verbatim = False
if rank == 0:
    verbatim = True
if args.specific_file == None:
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=args.end_time, form_time=sink_form_time)
    m_times = m_times[args.start_frame:]
    usable_files = mym.find_files(m_times, files, sink_form_time,sink_id, verbatim=False)
else:
    usable_files = [args.specific_file]
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=args.end_time, form_time=sink_form_time)
if rank == 0:
    print("Usable_files=", usable_files)
    
no_frames = len(m_times)
frames = list(range(args.start_frame, no_frames))
del sink_form_time
del files
    
sys.stdout.flush()
CW.Barrier()

#Trying yt parallelism
file_int = -1
for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/(8*4))): #8 projection and 4 fields.
    fn = usable_files[fn_it]
    print("File", fn, "is going to rank", rank)
    if size > 1:
        file_int = usable_files.index(fn)
    else:
        file_int = file_int + 1 #This was implemented for the specific cases where dt < simulation dump frequency, and the same file might be used for multiple times.
    
    if args.plot_time != None:
        if os.path.exists(save_dir + "time_" + (str(int(args.plot_time)))) == False:
            try:
                os.makedirs(save_dir + "time_" + (str(int(args.plot_time))))
            except:
                print(save_dir + "time_" + (str(int(args.plot_time))), "Already exists")
    else:
        if os.path.exists(save_dir + "movie_frame_" + ("%06d" % frames[file_int])) == False:
            try:
                os.makedirs(save_dir + "movie_frame_" + ("%06d" % frames[file_int]))
            except:
                print(save_dir + "movie_frame_" + ("%06d" % frames[file_int]), "Already exists")
    
    if len(glob.glob(save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/*.pkl")) == 8:
        make_pickle = False
        print("All projections for this time have been made")
    else:
        make_pickle = True
        
    if args.plot_time is None:
        pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/"
    else:
        pickle_file = save_dir + "time_" + (str(int(args.plot_time))) + "/"
        
    if make_pickle == True:
        ds = yt.load(fn, units_override=units_override)
        dd = ds.all_data()
        
        part_info = mym.get_particle_data(ds, sink_id=sink_id)
        
        time_val = m_times[file_int]
        
        center_pos = dd['Center_Position'].in_units('AU')
        center_vel = dd['Center_Velocity'].in_units('cm/s')
        
        part_posx = dd['sink_particle_posx'][sink_id:].in_units('AU') - center_pos[0]
        part_posy = dd['sink_particle_posy'][sink_id:].in_units('AU') - center_pos[1]
        part_posz = dd['sink_particle_posz'][sink_id:].in_units('AU') - center_pos[2]
        
        #Calculatae separation vector
        pos_array = yt.YTArray([part_posx, part_posy, part_posz]).T
        separation = pos_array[0] - pos_array[1]
        separation_magnitude = np.sqrt(separation[0]**2 + separation[1]**2 + separation[2]**2)
        sepration_unit = separation/separation_magnitude
        
        #Angle around xy plane:
        sep_xy_mag = np.sqrt(separation[0]**2 + separation[1]**2)
        sep_xy_unit = yt.YTArray([separation[0]/sep_xy_mag, separation[1]/sep_xy_mag, 0], 'AU')
        theta = np.arccos(np.dot(sep_xy_unit, np.array([0,1,0])))
        #Angle from z axis:
        phi = np.arccos(np.dot(sepration_unit, [0,0,1]))
        
        #Calculate angle alpha between projection and separation vector such that the projected separation is what you selected (default 200AU)
        projected_separation = yt.YTQuantity(args.projected_separation, 'AU')
        alpha = np.arcsin(projected_separation/separation_magnitude)
        proj_length = np.sqrt(separation_magnitude**2 - projected_separation**2)
        radius = proj_length*(projected_separation/separation_magnitude)
        sep_z = np.sqrt(proj_length**2 - radius**2)
        
        vectors_along_cone = np.array([[radius, 0, sep_z],\
                                       [radius/np.sqrt(2), radius/np.sqrt(2), sep_z],\
                                       [0, radius, sep_z],\
                                       [-1*radius/np.sqrt(2), radius/np.sqrt(2), sep_z],\
                                       [-radius, 0, sep_z],\
                                       [-1*radius/np.sqrt(2), -1*radius/np.sqrt(2), sep_z],\
                                       [0, -radius, sep_z],\
                                       [radius/np.sqrt(2), -1*radius/np.sqrt(2), sep_z]])
        
        
        #Figure out how to rotate the projection vectors to be the same orientation as the separation vector.
        #Rotate around XY plane
        z_rot = z_rotation_matrix(-1*phi)
        #z_rot_rev = z_rotation_matrix(phi)
        #xy_rot_rev = xy_rotation_matrix(-1*theta)
        theta_pos = yt.YTArray(np.dot(xy_rotation_matrix(theta), np.dot(z_rot, [0,0,1])), 'AU').value - sepration_unit
        theta_neg = yt.YTArray(np.dot(xy_rotation_matrix(-1*theta), np.dot(z_rot, [0,0,1])), 'AU').value - sepration_unit
        theta_pos_len = np.sqrt(np.sum(theta_pos**2))
        theta_neg_len = np.sqrt(np.sum(theta_neg**2))
        if theta_pos_len<theta_neg_len:
            xy_rot = xy_rotation_matrix(theta)
        elif theta_pos_len>theta_neg_len:
            xy_rot = xy_rotation_matrix(-1*theta)
        else:
            print("PROBLEM WITH FINDING CORRECT ROTATION")
            import pdb
            pdb.set_trace()
        
        #This section transforms the vectors_along_cone to be orientation for the separation vector
        projection_vectors = []
        north_vectors = []
        for vector in vectors_along_cone:
            if separation_magnitude > projected_separation:
                proj_vector = yt.YTArray(np.dot(xy_rot, np.dot(z_rot, vector)), 'AU')
            
                Proj_sep_proj = projected_vector(separation,proj_vector)
                sep_on_proj_length = np.sqrt(Proj_sep_proj[0]**2 + Proj_sep_proj[1]**2 + Proj_sep_proj[2]**2)
                calculated_proj_separation = np.sqrt(separation_magnitude**2 - sep_on_proj_length**2)
                if np.round(calculated_proj_separation) != projected_separation:
                    print("CALCULATED PROJECTED SEPARATION IS", np.round(calculated_proj_separation), "FOR FILE", fn, "AT TIME", time_val)
                    import pdb
                    pdb.set_trace()
                projection_vectors.append(proj_vector)
                north_vector = separation - Proj_sep_proj
                north_vectors.append(north_vector)
            else:
                projection_vectors.append(np.array([np.nan, np.nan, np.nan]))
                north_vectors.append(np.array([np.nan, np.nan, np.nan]))
        
        #Now that projection and north vectors have been generated, lets create the projection
        for proj_it in yt.parallel_objects(range(len(projection_vectors)), njobs=int(8)):# range(len(projection_vectors)):
            if True in np.isnan(projection_vectors[proj_it]):
                print("Skipping projection because  projection vector is Nan")
            elif args.plot_time == None and os.path.exists(save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it) + ".pkl"):
                print("Skipping because", save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it) + ".pkl", "exists")
            elif args.plot_time != None and os.path.exists(save_dir + "time_" + (str(int(args.plot_time))) + "/projection_" + str(proj_it) + ".pkl"):
                print("Skipping because", save_dir + "time_" + (str(int(args.plot_time))) + "/projection_" + str(proj_it) + ".pkl", "exists")
            else:
                #Calculate projected particle positions
                projected_particle_posy = projected_vector(pos_array, north_vectors[proj_it])
                proj_part_y_mag = np.sqrt(np.sum((projected_particle_posy**2), axis=1))
                proj_part_y_unit = (projected_particle_posy.T/proj_part_y_mag).T
                north_mag = np.sqrt(north_vectors[proj_it][0]**2 + north_vectors[proj_it][1]**2 + north_vectors[proj_it][2]**2)
                north_unit = north_vectors[proj_it]/north_mag
                north_sign = np.dot(north_unit, proj_part_y_unit.T)
                proj_part_y = proj_part_y_mag*north_sign
                proj_part_y = np.nan_to_num(proj_part_y)
                
                proj_vector_mag = np.sqrt(np.sum(projection_vectors[proj_it]**2))
                proj_vector_unit = projection_vectors[proj_it]/proj_vector_mag
                east_unit_vector = np.cross(north_unit, proj_vector_unit)
                
                projected_particle_posx = projected_vector(pos_array, east_unit_vector)
                proj_part_x_mag = np.sqrt(np.sum((projected_particle_posx**2), axis=1))
                proj_part_x_unit = (projected_particle_posx.T/proj_part_x_mag).T
                east_sign = np.dot(east_unit_vector, proj_part_x_unit.T)
                proj_part_x = proj_part_x_mag*east_sign
                proj_part_x = np.nan_to_num(proj_part_x)
                
                projected_particle_posz = projected_vector(pos_array, proj_vector_unit)
                proj_part_z_mag = np.sqrt(np.sum((projected_particle_posz**2), axis=1))
                proj_part_z_unit = (projected_particle_posz.T/proj_part_z_mag).T
                proj_sign = np.dot(proj_vector_unit, proj_part_z_unit.T)
                proj_part_z = proj_part_z_mag*proj_sign
                proj_part_z = np.nan_to_num(proj_part_z)
                part_info.update({'particle_position_z':proj_part_z.value})
                
                part_info['particle_position'] = np.array([[proj_part_x[0].value, proj_part_x[1].value],[proj_part_y[0].value, proj_part_y[1].value]])
                
                #Calculate center velocity
                center_vel_proj_y = projected_vector(center_vel, north_vectors[proj_it])
                center_vel_y = np.sqrt(center_vel_proj_y.T[0]**2 + center_vel_proj_y.T[1]**2 + center_vel_proj_y.T[2]**2).in_units('cm/s')
                
                center_vel_proj_x = projected_vector(center_vel, east_unit_vector)
                center_vel_x = np.sqrt(center_vel_proj_x.T[0]**2 + center_vel_proj_x.T[1]**2 + center_vel_proj_x.T[2]**2).in_units('cm/s')
                
                center_vel_proj_rv = projected_vector(center_vel, proj_vector_unit)
                center_vel_rv_mag = np.sqrt(np.sum(center_vel_proj_rv**2))
                center_vel_rv_unit = center_vel_proj_rv/center_vel_rv_mag
                rv_sign = np.dot(proj_vector_unit, center_vel_rv_unit)
                center_vel_rv = center_vel_rv_mag*rv_sign
                
                center_vel_image = np.array([center_vel_x, center_vel_y])
            
                #set vectors in my fields:
                myf.set_normal(proj_vector_unit)
                myf.set_east_vector(east_unit_vector)
                myf.set_north_vector(north_unit)
                
                #These numbers are just creative pickle file management
                div_32 = int(rank/32)
                rem_32 = np.remainder(rank,32)
                div_4 = int(rem_32/4)
                proj_root_rank = div_32*32 + div_4*4

                #Fields to project
                field_list = [field, ('gas', 'Radial_Velocity'), ('gas', 'Proj_x_velocity'), ('gas', 'Proj_y_velocity')]
                proj_dict = {field[1]:[], 'Radial_Velocity':[], 'Proj_x_velocity':[], 'Proj_y_velocity':[]}
                
                proj_dict_keys = str(proj_dict.keys()).split("['")[1].split("']")[0].split("', '")
                for field in yt.parallel_objects(field_list):
                    print("Calculating projection with normal", proj_vector_unit, "for field", field, "on rank", rank)
                    proj = yt.OffAxisProjectionPlot(ds, proj_vector_unit, field, width=(2*args.ax_lim, 'AU'), weight_field=weight_field, method='integrate', center=(center_pos.value, 'AU'), depth=(args.slice_thickness, 'AU'), north_vector=north_unit)
                    
                    if args.resolution != 800:
                        proj.set_buff_size([args.resolution, args.resolution])
                
                    if args.field in str(field):
                        if weight_field == None:
                            if args.divide_by_proj_thickness == "True":
                                proj_array = np.array((proj.frb.data[field]/thickness.in_units('cm')).in_units(args.field_unit))
                            else:
                                proj_array = np.array(proj.frb.data[field].in_units(args.field_unit+"*cm"))
                        else:
                            if args.divide_by_proj_thickness == "True":
                                proj_array = np.array(proj.frb.data[field].in_units(args.field_unit))
                            else:
                                proj_array = np.array(proj.frb.data[field].in_units(args.field_unit)*thickness.in_units('cm'))
                    else:
                        if weight_field == None:
                            proj_array = np.array(proj.frb.data[field].in_cgs()/thickness.in_units('cm'))
                        else:
                            proj_array = np.array(proj.frb.data[field].in_cgs())
                            
                    if rank == proj_root_rank:
                        proj_dict[field[1]] = proj_array
                    else:
                        print("Dumping projection data into", pickle_file + 'proj_data_' + str(proj_root_rank)+ str(proj_dict_keys.index(field[1])) + '.pkl')
                        file = open(pickle_file + 'proj_data_' + str(proj_root_rank)+ str(proj_dict_keys.index(field[1])) + '.pkl', 'wb')
                        pickle.dump((field[1], proj_array), file)
                        file.close()
                
                if rank == proj_root_rank and size > 1:
                    #check to see if all proj files exist yet
                    for kit in range(1,len(proj_dict_keys)):
                        file = open(pickle_file + 'proj_data_' +str(proj_root_rank) +str(kit)+'.pkl', 'rb')
                        key, proj_array = pickle.load(file)
                        file.close()
                        proj_dict[key] = proj_array
                        os.remove(pickle_file + 'proj_data_' +str(proj_root_rank) +str(kit)+'.pkl')
                        
                if rank == proj_root_rank:
                    image = yt.YTArray(proj_dict[proj_dict_keys[0]], args.field_unit)
                    vel_rad = yt.YTArray(proj_dict[proj_dict_keys[1]], 'cm/s')
                    velx_full = yt.YTArray(proj_dict[proj_dict_keys[2]], 'cm/s')
                    vely_full = yt.YTArray(proj_dict[proj_dict_keys[3]], 'cm/s')
                    
                    if args.image_center == 2:
                        image = np.flip(np.flip(image, axis=0), axis=1)
                        vel_rad = np.flip(np.flip(vel_rad, axis=0), axis=1)
                        velx_full = np.flip(np.flip(velx_full, axis=0), axis=1)
                        vely_full = np.flip(np.flip(vely_full, axis=0), axis=1)
                        part_info['particle_position'][1][0] = part_info['particle_position'][1][0]*-1

                    image_dict = {'time': time_val, 'center_vel_image':center_vel_image, 'proj_vector': proj_vector_unit, 'field': field}
                    
                    pickle_file = pickle_file + 'projection_' + str(proj_it) + '.pkl'
                    file = open(pickle_file, 'wb')
                    pickle.dump((image, vel_rad, velx_full, vely_full, part_info, image_dict), file)
                    file.close()
                    print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
                    del time_val
                    del dd
                    del proj
                    del image
                    del vel_rad
                    del velx_full
                    del vely_full
sys.stdout.flush()
CW.Barrier()

#Section to plot figures:
print("Finished generating projection pickles")
