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
import matplotlib.pyplot as plt
import matplotlib as cm
from matplotlib.colors import LogNorm
from matplotlib import transforms
import matplotlib.patheffects as path_effects

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="Number_Density")
    parser.add_argument("-f_unit", "--field_unit", help="What units would you like to plot the field?", default="cm**-3")
    parser.add_argument("-div_by_thickness", "--divide_by_proj_thickness", help="Woudl you like to divide the field by the thickness of the projection?", default="True", type=str)
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 100., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default=0, type=int)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-vaf", "--velocity_annotation_frequency", help="how many velocity vectors do you want annotated across one side?", type=float, default=31.)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-mt", "--movie_times", help="What movies times would you like plotted?", type=list, default=[])
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default=None)#'1.e-16')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=None)#1.e-14)
    parser.add_argument("-ic", "--image_center", help="where would you like to center the image?", type=int, default=0)
    parser.add_argument("-cc", "--calculation_center", help="where would you like to calculate center positionn and velocity?", type=int, default=0)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=250)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-proj_or", "--projection_orientation", help="Do you want to set the projection orientation? give as angle (in degrees) from positive y-axis", default=None, type=float)
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 300AU", type=float, default=500.)
    parser.add_argument("-wf", "--weight_field", help="Do you want to have a weighted projection plot?", type=str, default=None)
    parser.add_argument("-use_gas", "--use_gas_center_calc", help="Do you want to use gas when calculating the center position adn veloity?", type=str, default='True')
    parser.add_argument("-all_files", "--use_all_files", help="Do you want to make frames using all available files instead of at particular time steps?", type=str, default='False')
    parser.add_argument("-update_alim", "--update_ax_lim", help="Do you want to update the axes limits by taking away the center position values or not?", type=str, default='False')
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-frames_only", "--make_frames_only", help="do you only want to make frames?", default='False', type=str)
    parser.add_argument("-debug", "--debug_plotting", help="Do you want to debug why plotting is messing up", default='False', type=str)
    parser.add_argument("-res", "--resolution", help="define image resolution", default=800, type=int)
    parser.add_argument("-active_rad", "--active_radius", help="within what radius of the centered sink do you want to consider when using sink and gas for calculations", type=float, default=10000.0)
    parser.add_argument("-proj_sep", "--projected_separation", help="if you want to make a projection such that the separation is a particular ammount, what is that?", type=float, default=200.0)
    parser.add_argument("-use_part_for_vel", "--use_particle_for_center_vel_calc", help="Do you want to use the particles to calculate center velocity?", type=str, default='True')
    parser.add_argument("-image_dir", "--image_directory", help="The directory where the projection pickles are stored", type=str, default="./")
    parser.add_argument("-threshold", "--density_threshold", help="What number density threshold would you like to use?", type=float, default=0.0)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def sim_info(ds,args):
    """
    Finds particle info, relevant to frame size and such. NOTE ACCRETION RADIUS IS GIVEN FROM PARTICLE INFO FUNCTION
    """
    dd = ds.all_data()
    field_it = [i for i, v in enumerate(ds.derived_field_list) if v[1] == args.field][0]
    field = ds.derived_field_list[field_it]
    dim = args.resolution
    if args.ax_lim == None:
        xmin = -1000
        xmax = 1000
        ymin = -1000
        ymax = 1000
    else:
        xmin = -1*args.ax_lim
        xmax = args.ax_lim
        ymin = -1*args.ax_lim
        ymax = args.ax_lim
    cl = (xmax-xmin)/dim
    annotate_freq = dim/args.velocity_annotation_frequency
    smoothing = annotate_freq/2
    unit_string = str(dd[field[1]].in_cgs().units)
    split_string = unit_string.split('**')
    unit_string = "^".join(split_string)
    split_string = unit_string.split('*')
    unit_string = " ".join(split_string)
    sim_info = {'field': field,
                'dimension': dim,
                'xmin': xmin,
                'xmax': xmax,
                'ymin': ymin,
                'ymax': ymax,
                'cell_length': cl,
                'annotate_freq': annotate_freq,
                'smoothing': smoothing,
                'unit_string': unit_string
                }
    del field_it
    del field
    del dim
    del xmin
    del xmax
    del ymin
    del ymax
    del cl
    del annotate_freq
    del smoothing
    del unit_string
    del dd
    return sim_info

def image_properties(X, Y, args, sim_info):
    xlabel = 'Distance from center (AU)'
    ylabel = 'Distance from center (AU)'
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    return xlabel, ylabel, xlim, ylim

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
    Calculates the position of vecter projected onto vector
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

#Set some plot variables independant on data files
cbar_max = args.colourbar_max
if args.colourbar_min == None:
    cbar_min = args.colourbar_min
else:
    try:
        cbar_min = float(args.colourbar_min)
    except:
        cbar_min = float(args.colourbar_min[1:])
    
title_parts = args.title.split('_')
title = ''
for part in title_parts:
    if part != title_parts[-1]:
        title = title + part + ' '
    else:
        title = title + part
mym.set_global_font_size(args.text_font)

sys.stdout.flush()
CW.Barrier()

#File files
files = sorted(glob.glob(input_dir+"*/info*.txt"))

sys.stdout.flush()
CW.Barrier()

#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif args.simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif args.simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
mym.set_units(units_override)

#find sink particle to center on and formation time
ds = yt.load(files[-1], units_override=units_override)

#Get simulation information
if rank == 0:
    print("loading fields")
simfo = sim_info(ds, args)
x = np.linspace(simfo['xmin'], simfo['xmax'], simfo['dimension'])
y = np.linspace(simfo['ymin'], simfo['ymax'], simfo['dimension'])
X, Y = np.meshgrid(x, y)
annotate_space = (simfo['xmax'] - simfo['xmin'])/args.velocity_annotation_frequency
x_ind = []
y_ind = []
counter = 0
while counter < args.velocity_annotation_frequency:
    val = annotate_space*counter + annotate_space/2. + simfo['xmin']
    x_ind.append(int(val))
    y_ind.append(int(val))
    counter = counter + 1
xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
if args.ax_lim != None:
    xlim = [-1*args.ax_lim, args.ax_lim]
    ylim = [-1*args.ax_lim, args.ax_lim]
x_width = (xlim[1] -xlim[0])
y_width = (ylim[1] -ylim[0])
thickness = yt.YTQuantity(args.slice_thickness, 'AU')
rv_channels = np.arange(2,40,2)
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

if args.weight_field == 'None':
    weight_field = None
else:
    weight_field = args.weight_field

if args.plot_time != None:
    if args.weight_field == 'None':
        weight_field = None
        pickle_file = save_dir + "time_" + (str(int(args.plot_time))) + "/"
    else:
        weight_field = args.weight_field
        pickle_file = save_dir + "time_" + (str(int(args.plot_time))) + "_weight_" + args.field + "/"
sys.stdout.flush()
CW.Barrier()

dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
sink_form_companion = dd['sink_particle_form_time'][sink_id+1]
if args.start_frame == 0 and args.plot_time == None:
    args.start_frame = int((sink_form_companion - sink_form_time)/(args.time_step))+1
del dd

if rank == 0:
    print("Finding times")
if args.plot_time != None:
    m_times = [args.plot_time]
else:
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=args.end_time, form_time=sink_form_time)
if rank == 0:
    print("Found times")
    
no_frames = len(m_times)
m_times = m_times[args.start_frame:]
frames = list(range(args.start_frame, no_frames))
    
sys.stdout.flush()
CW.Barrier()

if args.make_frames_only == 'False':
    verbatim = False
    if rank == 0:
        verbatim = True
    usable_files = mym.find_files(m_times, files, sink_form_time,sink_id, verbatim=False)
    del sink_form_time
    del files
    
sys.stdout.flush()
CW.Barrier()

if args.make_frames_only == 'False':
    #Trying yt parallelism
    file_int = -1
    #rit = 32
    #for fn in usable_files:
    for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/(18*8))):
        fn = usable_files[fn_it]
        print("File", fn, "is going to rank", rank)
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
            '''
            #if usable_files[file_int] == usable_files[file_int-1]:
                #os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
            '''
        if args.plot_time != None:
            if os.path.exists(pickle_file) == False:
                try:
                    os.makedirs(pickle_file)
                except:
                    print(pickle_file, "Already exists")
        
        else:
            if os.path.exists(save_dir + "movie_frame_" + ("%06d" % frames[file_int])) == False:
                try:
                    os.makedirs(save_dir + "movie_frame_" + ("%06d" % frames[file_int]))
                except:
                    print(save_dir + "movie_frame_" + ("%06d" % frames[file_int]), "Already exists")
        if args.plot_time == None:
            if len(glob.glob(save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/*/*.pkl")) == (18*8):
                make_pickle = False
                print("All channel maps for this file have been created")
            else:
                make_pickle = True
        else:
            if len(glob.glob(pickle_file + "*/*.pkl")) == (18*8):
                make_pickle = False
                print("All channel maps for this file have been created")
            else:
                make_pickle = True
        if make_pickle == True:
            ds = yt.load(fn, units_override=units_override)
            dd = ds.all_data()
            
            has_particles = has_sinks(ds)
            part_info = mym.get_particle_data(ds, sink_id=sink_id)
            
            time_val = m_times[file_int]
            
            center_pos = dd['Center_Position'].in_units('AU')
            center_vel = dd['Center_Velocity'].in_units('cm/s')
            
            part_posx = dd['sink_particle_posx'][sink_id:].in_units('AU') - center_pos[0]
            part_posy = dd['sink_particle_posy'][sink_id:].in_units('AU') - center_pos[1]
            part_posz = dd['sink_particle_posz'][sink_id:].in_units('AU') - center_pos[2]
            
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
            
            #Calculate angle alpha between projection and separation vector such that the projected sepation is what you selected (default 200AU)
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
                                           
            #cone_vector_length_no_tan = np.sqrt(vectors_along_cone_no_tan[:,0]**2 + vectors_along_cone_no_tan[:,1]**2 + vectors_along_cone_no_tan[:,2]**2)
            #vectors_along_cone_no_tan = vectors_along_cone_no_tan/cone_vector_length_no_tan[0]
            '''
            vectors_along_cone = np.array([[separation_magnitude*np.tan(alpha), 0, separation_magnitude],\
                                           [-1*separation_magnitude*np.tan(alpha), 0, separation_magnitude],\
                                           [0, separation_magnitude*np.tan(alpha), separation_magnitude],\
                                           [0, -1*separation_magnitude*np.tan(alpha), separation_magnitude],\
                                           [(separation_magnitude*np.tan(alpha))/np.sqrt(2), (separation_magnitude*np.tan(alpha))/np.sqrt(2), separation_magnitude],\
                                           [-1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), ((separation_magnitude*np.tan(alpha))/np.sqrt(2)), separation_magnitude],\
                                           [(separation_magnitude*np.tan(alpha))/np.sqrt(2), -1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), separation_magnitude],\
                                           [-1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), -1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), separation_magnitude]])
            '''
            #cone_vector_length = np.sqrt(vectors_along_cone[:,0]**2 + vectors_along_cone[:,1]**2 + vectors_along_cone[:,2]**2)
            #vectors_along_cone = vectors_along_cone/cone_vector_length[0]
            
            
            #Figure out how to rotate the projection vectors to be the same reference as the separation vector.
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
            for proj_it in yt.parallel_objects(range(len(projection_vectors)), njobs=8):# range(len(projection_vectors)):
                if args.plot_time == None:
                    try:
                        if os.path.exists(save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it)) == False:
                            os.makedirs(save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it))
                            print("Created directory:", save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it))
                    except:
                        print("save directory:", save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it), "exists")
                    pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it)
                    if len(glob.glob(save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + "/projection_" + str(proj_it) + "/*.pkl")) == 18:
                        make_proj = False
                        print("All channels for this projection have been made")
                    else:
                        make_proj = True
                elif True in np.isnan(projection_vectors[proj_it]):
                    print("Skipping projection because vector is Nan")
                    make_proj = False
                else:
                    try:
                        if os.path.exists(pickle_file + "/projection_" + str(proj_it)) == False:
                            os.makedirs(pickle_file + "/projection_" + str(proj_it))
                            print("Created directory:", pickle_file + "/projection_" + str(proj_it))
                    except:
                        print("save directory:", pickle_file + "/projection_" + str(proj_it), "exists")
                    if len(glob.glob(pickle_file + "/projection_" + str(proj_it) + "/*.pkl")) == 18:
                        make_proj = False
                        print("All channels for this projection have been made")
                    else:
                        make_proj = True
                if make_proj:
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
                    east_unit_vector = np.cross(proj_vector_unit, north_unit)
                    #east_unit_vector = np.cross(north_unit, proj_vector_unit)
                    
                    projected_particle_posx = projected_vector(pos_array, east_unit_vector)
                    proj_part_x_mag = np.sqrt(np.sum((projected_particle_posx**2), axis=1))
                    proj_part_x_unit = (projected_particle_posx.T/proj_part_x_mag).T
                    east_sign = np.dot(east_unit_vector, proj_part_x_unit.T)
                    proj_part_x = proj_part_x_mag*east_sign
                    proj_part_x = np.nan_to_num(proj_part_x)
                    
                    part_info['particle_position'] = np.array([[proj_part_x[0].value, proj_part_x[0].value],[proj_part_y[0].value, proj_part_y[1].value]])
                    
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
                
                    #set vectors:
                    myf.set_normal(proj_vector_unit)
                    myf.set_east_vector(east_unit_vector)
                    myf.set_north_vector(north_unit)
                    
                    #Caculate pojections!
                    '''
                    div_32 = int(rank/32)
                    rem_32 = np.remainder(rank,32)
                    div_4 = int(rem_32/4)
                    proj_root_rank = div_32*32 + div_4*4
                    '''
                    for rv_channel_it in yt.parallel_objects(range(len(rv_channels)-1)):
                        if os.path.exists(save_dir +'/projection_' +str(proj_it) + '/channel_' + str(int(rv_channels[rv_channel_it])) + '_' + str(int(rv_channels[rv_channel_it+1])) + '.pkl') == False:
                            make_channel_proj = True
                        else:
                            make_channel_proj = False
                            print("Channel map of", str(int(rv_channels[rv_channel_it])) + '_' + str(int(rv_channels[rv_channel_it+1])), "for projeciton", proj_it, "Has been made, so skipping this")
                        if make_channel_proj:
                            #rv_cut_region = dd.cut_region(["(abs(obj['Radial_Velocity'].in_units('km/s')) < " + str(rv_channels[rv_channel_it+1]) + ") & (abs(obj['Radial_Velocity'].in_units('km/s')) > " + str(rv_channels[rv_channel_it]) + ") & (obj['"+ str(simfo['field'][1])+"'].in_units('"+str(args.field_unit)+"') > " + str(args.density_threshold) +")"])
                            rv_cut_region = dd.cut_region(["(abs(obj['Radial_Velocity'].in_units('km/s')) < " + str(rv_channels[rv_channel_it+1]) + ") & (abs(obj['Radial_Velocity'].in_units('km/s')) > " + str(rv_channels[rv_channel_it]) + ")"])
                        
                            print("Calculating projection with normal", projection_vectors[proj_it], "for rv channel", [rv_channels[rv_channel_it], rv_channels[rv_channel_it+1]], "on rank", rank)
                            field_list = [simfo['field'], ('gas', 'Radial_Velocity')]
                            try:
                                proj = yt.OffAxisProjectionPlot(ds, proj_vector_unit, field_list, width=(x_width, 'AU'), weight_field=weight_field, method='integrate', center=(center_pos.value, 'AU'), depth=(args.slice_thickness, 'AU'), north_vector=north_unit, data_source=rv_cut_region)
                                proj.set_buff_size([args.resolution, args.resolution])
                                
                                for field in field_list:
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
                                    if field == simfo['field']:
                                        image = yt.YTArray(proj_array, args.field_unit)
                                    else:
                                        sign_array = proj_array
                                
                                del proj
                            except:
                                print("Couldn't make projection so returning arrays of NaN on rank", rank)
                                image = yt.YTArray(np.ones((800,800))*np.nan,'g/cm**2')
                                sign_array = yt.YTArray(np.ones((800,800))*np.nan,'g/cm**2')
                            
                            args_dict = {}
                            if args.annotate_time == "True":
                                args_dict.update({'annotate_time': '$t$='+str(int(time_val))+'yr'})
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
                            
                            pickle_file = pickle_file + '/channel_' + str(int(rv_channels[rv_channel_it])) + '_' + str(int(rv_channels[rv_channel_it+1])) + '.pkl'
                            
                            file = open(pickle_file, 'wb')
                            pickle.dump((X, Y, image, sign_array, part_info, args_dict, simfo), file)
                            file.close()
                            print("Created Pickle:", pickle_file, "for  file:", str(ds))
                            del has_particles
                            del time_val
                            del dd
                            del center_vel
                            del image
                            del part_info
        sys.stdout.flush()
        CW.Barrier()
            
sys.stdout.flush()
CW.Barrier()
            
#Section to plot figures:
print("Finished generating channel pickles")
if args.plot_time is None:
    pickle_files = sorted(glob.glob(save_dir+"*/*/channel*.pkl"))
else:
    pickle_files = sorted(glob.glob(save_dir+"time_" + str(int(args.plot_time)) + "*/*/channel*.pkl"))
rit = -1
for pickle_file in pickle_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        proj_number = pickle_file.split('_')[-3]
        print("on rank", rank, "using pickle_file", pickle_file)
        file = open(pickle_file, 'rb')
        X, Y, image, sign_array, part_info, args_dict, simfo = pickle.load(file)
        image[image==0] = np.nan
        file.close()
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        
        time_val = args_dict['time_val']
        
        if np.round(np.mean(args_dict['xlim'])) == np.round(np.mean(X)):
            xlim = args_dict['xlim']
            ylim = args_dict['ylim']
        else:
            xlim = args_dict['xlim'] + np.mean(X)
            ylim = args_dict['ylim'] + np.mean(Y)
        has_particles = args_dict['has_particles']
        time_val = args_dict['time_val']
        xabel = args_dict['xabel']
        yabel = args_dict['yabel']
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
        ax.set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        rv_channel = [pickle_file.split('_')[-2], pickle_file.split('_')[-1].split('.')[0]]
        
        if len(m_times) > 1:
            if args.output_filename == None:
                file_name = pickle_file.split('.pkl')[0]
            else:
                file_name = args.output_filename + "_" + str(int(time_val))
        else:
            if args.output_filename != None:
                file_name = args.output_filename
            else:
                file_name = pickle_file.split('.pkl')[0]
                
        try:
            if args.plot_time != None:
                image_file=os.getcwd().split('Channel_maps')[0] + 'Time_Series' + os.getcwd().split('Channel_maps')[1] + '/time_' + str(int(args.plot_time)) + '/' + pickle_file.split('/')[2]+ '.pkl'
            else:
                image_file = os.getcwd().split('Channel_maps')[0] + 'Time_Series' + os.getcwd().split('Channel_maps')[1] + '/' + pickle_file.split('/')[-3] + '/' + pickle_file.split('/')[-2] + '.pkl'
            file = open(image_file, 'rb')
            stuff = pickle.load(file)
            file.close()
            background = stuff[2]
            if None in (cbar_min, cbar_max):
                plot = ax.pcolormesh(X, Y, background, cmap=plt.cm.Greys, rasterized=True, alpha=0.5)
            else:
                plot = ax.pcolormesh(X, Y, background, cmap=plt.cm.Greys, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, alpha=0.1)
            print("plotted image")
        except:
            print("Couldn't get image data")
        '''
        if None in (cbar_min, cbar_max):
            plot = ax.pcolormesh(X, Y, sign_array, cmap=plt.cm.seismic_r, rasterized=True, alpha=0.1)
        elif 0.0 in (cbar_min, cbar_max):
            plot = ax.pcolormesh(X, Y, sign_array, cmap=plt.cm.seismic_r, rasterized=True, vmin=-1*cbar_max, vmax=cbar_max, alpha=0.1)
        else:
            plot = ax.pcolormesh(X, Y, sign_array, cmap=plt.cm.seismic_r, norm=LogNorm(vmin=-1*cbar_max, vmax=cbar_max), rasterized=True, alpha=0.1)
        '''
        #plot_rv = ax.pcolormesh(X, Y, sign_array, cmap=plt.cm.seismic_r, rasterized=True, alpha=0.05)#, vmin=-1*cbar_max, vmax=cbar_max)
        #plt.gca().set_aspect('equal')
        #cbar = plt.colorbar(plot_rv, pad=0.0)
        non_nan_inds = np.where(np.isnan(image) == False)
        if len(non_nan_inds[0]) > 0:
            std = np.std(image[non_nan_inds])
            max = np.max(image[non_nan_inds])
            mean = np.mean(image[non_nan_inds])
            #mean = np.mean(image[non_nan_inds])
            try:
                level_b = np.arange(1, max.value, std.value*3)
                level_r = -1*np.arange(1, max.value, std.value*3)[::-1]
                CS_b = ax.contour(X,Y,np.nan_to_num(image*np.sign(sign_array)), levels=level_b, linewidths=0.5, colors='b', alpha=1.0)
                #CS_b = ax.contour(X,Y,np.nan_to_num(image), levels=level_b, linewidths=0.5, colors='b', alpha=1.0)
                CS_r = ax.contour(X,Y,np.nan_to_num(image*np.sign(sign_array)), levels=level_r, linewidths=0.5, colors='r', alpha=1.0)
                #ax.clabel(CS,inline=1)
                print("Contours plotted")
            except:
                print("Couldn't plot contours")
        plt.gca().set_aspect('equal')
        #cbar = plt.colorbar(plot, pad=0.0)
        
        if has_particles:
            if args.annotate_particles_mass == True:
                mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
            else:
                mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)
        
        '''
        if 'Density' in simfo['field']:
            if args.divide_by_proj_thickness == "True":
                cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
            else:
                cbar.set_label(r"Column Density (g$\,$cm$^{-2}$)", rotation=270, labelpad=14, size=args.text_font)
        elif 'Number_Density' in simfo['field']:
            if args.divide_by_proj_thickness == 'True':
                cbar.set_label(r"Number Density (cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
            else:
                cbar.set_label(r"Number Density (cm$^{-2}$)", rotation=270, labelpad=14, size=args.text_font)
        else:
            label_string = simfo['field'][1] + ' ($' + args.field_unit + '$)'
            cbar.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=args.text_font)
        '''
        if len(title) > 0:
            title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
            title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

        plt.tick_params(axis='both', which='major')# labelsize=16)
        for line in ax.xaxis.get_ticklines():
            line.set_color('white')
        for line in ax.yaxis.get_ticklines():
            line.set_color('white')

        if args.annotate_time == "True":
            try:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                time_string = "$t$="+str(int(time_val))+"yr"
                time_string_raw = r"{}".format(time_string)
                time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                except:
                    print("Couldn't outline time string")
            except:
                print("Couldn't plot time string")
        
        if size > 1:
            try:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame of projection', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
            except:
                print("couldn't save for the dviread.py problem. Make frame " + str(proj_number) + " on ipython")
        else:
            try:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame of projection', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
            except:
                print("Couldn't save figure " + file_name + ".jpg")

print("Completed making frames on rank", rank)
