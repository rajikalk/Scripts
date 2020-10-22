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
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="Density")
    parser.add_argument("-f_unit", "--field_unit", help="What units would you like to plot the field?", default="g/cm**3")
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
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-16')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
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
    parser.add_argument("-proj_sep", "--projected_separation", help="if you want to make a projection such that the separation is a particular ammount, what is that?", type=float, default=300.0)
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
    proj_v = yt.YTArray(np.array([proj_v_x, proj_v_y,proj_v_z]).T, vector_units)
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
units_override = {"length_unit":(4.0,"pc"), "mass_unit":(2998,"Msun"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "density_unit":(46.84375, "Msun/pc**3")}
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
X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
if args.ax_lim != None:
    xlim = [-1*args.ax_lim, args.ax_lim]
    ylim = [-1*args.ax_lim, args.ax_lim]
x_width = (xlim[1] -xlim[0])
y_width = (ylim[1] -ylim[0])
thickness = yt.YTQuantity(args.slice_thickness, 'AU')
#Sets center for calculating center position and velocity
myf.set_center(args.image_center)

#Set to make sure that particles aren't used to calculate the center velocity
myf.set_com_vel_use_part(False)

if args.use_gas_center_calc == 'True':
    myf.set_com_pos_use_gas(True)
else:
    myf.set_com_pos_use_gas(False)
    
#Make sure to only use gas when calculating the center velocity

sys.stdout.flush()
CW.Barrier()

if args.plot_time != None:
    if args.weight_field == 'None':
        weight_field = None
        pickle_file = save_dir + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
    else:
        weight_field = args.weight_field
        pickle_file = save_dir + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + ".pkl"
       
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
del dd

if args.plot_time != None:
    m_times = [args.plot_time]
else:
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=args.presink_frames, end_time=args.end_time, form_time=sink_form_time)
    
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
    for fn in yt.parallel_objects(usable_files, njobs=int(size/5)):
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
        if make_pickle == True:
            ds = yt.load(fn, units_override=units_override)
            dd = ds.all_data()
            
            time_val = m_times[file_int]
            has_particles = has_sinks(ds)
            if has_particles:
                part_info = mym.get_particle_data(ds, sink_id=sink_id)
            else:
                part_info = {}
            
            center_pos = dd['Center_Position'].in_units('AU')
            
            part_posx = dd['sink_particle_posx'][sink_id:].in_units('AU') - center_pos[0]
            part_posy = dd['sink_particle_posy'][sink_id:].in_units('AU') - center_pos[1]
            part_posz = dd['sink_particle_posz'][sink_id:].in_units('AU') - center_pos[2]
            
            pos_array = yt.YTArray([part_posx, part_posy, part_posz]).T
            separation = pos_array[0] - pos_array[1]
            separation_magnitude = np.sqrt(separation[0]**2 + separation[1]**2 + separation[2]**2)
            
            #Angle around xy plane:
            sep_xy_mag = np.sqrt(separation[0]**2 + separation[1]**2)
            sep_xy_unit = yt.YTArray([separation[0]/sep_xy_mag, separation[1]/sep_xy_mag, 0], 'AU')
            theta = np.arccos(np.dot(sep_xy_unit, np.array([0,1,0])))
            #Angle from z axis:
            phi = np.arccos(np.dot(separation/separation_magnitude, [0,0,1]))
            
            #Calculate angle alpha between projection and separation vector such that the projected sepation is what you selected (default 200AU)
            projected_separation = yt.YTQuantity(args.projected_separation, 'AU')
            #DOUBLE CHECK ALPHA CALCULATION
            alpha = np.arcsin(projected_separation/separation_magnitude)
            
            vectors_along_cone = np.array([[separation_magnitude*np.tan(alpha), 0, separation_magnitude],\
                                           [-1*separation_magnitude*np.tan(alpha), 0, separation_magnitude],\
                                           [0, separation_magnitude*np.tan(alpha), separation_magnitude],\
                                           [0, -1*separation_magnitude*np.tan(alpha), separation_magnitude],\
                                           [(separation_magnitude*np.tan(alpha))/np.sqrt(2), (separation_magnitude*np.tan(alpha))/np.sqrt(2), separation_magnitude],\
                                           [-1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), ((separation_magnitude*np.tan(alpha))/np.sqrt(2)), separation_magnitude],\
                                           [(separation_magnitude*np.tan(alpha))/np.sqrt(2), -1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), separation_magnitude],\
                                           [-1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), -1*((separation_magnitude*np.tan(alpha))/np.sqrt(2)), separation_magnitude]])
            
            #Figure out how to rotate the projection vectors to be the same reference as the separation vector.
            #Rotate around XY plane
            z_rot = z_rotation_matrix(-1*phi)
            xy_rot = xy_rotation_matrix(theta)
            projection_vectors = []
            north_vectors = []
            for vector in vectors_along_cone:
                z_rot_vector = np.dot(z_rot, vector)
                proj_vector = yt.YTArray(np.dot(xy_rot, z_rot_vector), 'AU')
                proj_mag = np.sqrt(proj_vector[0]**2 + proj_vector[1]**2 + proj_vector[2]**2)
                proj_unit = proj_vector/proj_mag
                projection_vectors.append(proj_vector)
                
                Proj_sep_proj = (np.dot(separation,proj_vector)/np.dot(proj_vector, proj_vector))* proj_vector
                north_vector = separation - Proj_sep_proj
                north_vectors.append(north_vector)
                
                sep_unit = separation/separation_magnitude
                
                angle_sep_proj = np.degrees(np.arccos(np.dot(proj_unit, sep_unit)))
                print("Alpha =", np.degrees(alpha), "but angle between sep and proj is", angle_sep_proj)
            
            #Now that projection and north vectors have been generated, lets create the projection
            for proj_it in range(len(projection_vectors)):
                #Calculate projected particle positions
                projected_particle_positions = projected_vector(pos_array, north_vectors[proj_it])
                proj_part_mag = np.sqrt(np.sum((projected_particle_positions**2), axis=1))
                proj_part_unit = (projected_particle_positions.T/proj_part_mag).T
                
                north_mag = np.sqrt(north_vectors[proj_it][0]**2 + north_vectors[proj_it][1]**2 + north_vectors[proj_it][2]**2)
                north_unit = north_vectors[proj_it]/north_mag
                
                y_sign = np.dot(proj_part_unit,north_unit)
                particle_y_plot = y_sign*proj_part_mag
                part_info['particle_position'] = np.array([[0, 0],[particle_y_plot[0].value, particle_y_plot[1].value]])
                
                #Caculate pojections!
                proj_vector_mag = np.sqrt(np.sum(projection_vectors[proj_it]**2))
                proj_vector_unit = projection_vectors[proj_it]/proj_vector_mag
                field_ist = [simfo['field'], ('ramses', 'x-velocity'), ('ramses', 'y-velocity'), ('ramses', 'z-velocity')] # ('gas', mag1_field), ('gas', mag2_field)
                proj = yt.OffAxisProjectionPlot(ds, proj_vector_unit, field_ist, width=(x_width/2, 'AU'), weight_field=weight_field, method='integrate', center=(center_pos.value, 'AU'), depth=(args.slice_thickness, 'AU'), north_vector=north_unit)
                
                image = proj.frb.data[simfo['field']]#/thickness.in_units('cm').value
                vel_x = proj.frb.data[('ramses', 'x-velocity')].in_cgs()/thickness.in_units('cm')
                vel_y = proj.frb.data[('ramses', 'y-velocity')].in_cgs()/thickness.in_units('cm')
                vel_z = proj.frb.data[('ramses', 'z-velocity')].in_cgs()/thickness.in_units('cm')
                
                vel_array = yt.YTArray([vel_x, vel_y, vel_z]).T
                RV = projected_vector(vel_array, projection_vectors[proj_it])
                import pdb
                pdb.set_trace()
                '''
                velx_full = proj.frb.data[('gas', 'Projected_Velocity')].in_units('cm**2/s').value/thickness.in_units('cm').value
                vely_full = proj.frb.data[('flash', 'velz')].in_units('cm**2/s').value/thickness.in_units('cm').value
                magx = proj.frb.data[('gas', 'Projected_Magnetic_Field')].in_units('cm*gauss').value/thickness.in_units('cm').value
                magy = proj.frb.data[('flash', 'magz')].in_units('cm*gauss').value/thickness.in_units('cm').value
                
                velx, vely = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=center_vel)
                '''
                
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
                
                file = open(pickle_file, 'wb')
                if args.absolute_image != "False":
                    image = abs(image)
                pickle.dump((X_image, Y_image, image, magx, magy, X_image_vel, Y_image_vel, velx, vely, part_info, args_dict, simfo), file)
                file.close()
                print("Created Pickle:", pickle_file, "for  file:", str(ds))
                del has_particles
                del time_val
                del x_width
                del y_width
                del thickness
                del dd
                del center_vel
                del proj
                del image
                del magx
                del magy
                del velx
                del vely
                del part_info
                
                #Create figure section. ONCE THIS IS FINALISED MOVE TO BE OUTSIDE OF THE FOR LOOP SO THEY CAN BE CREATE FROM PICKLE
                print("on rank,", rank, "using pickle_file", pickle_file)
                file = open(pickle_file, 'rb')
                X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo = pickle.load(file)
                #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                file.close()
                
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
                
                if len(usable_files) > 1:
                    if args.output_filename == None:
                        file_name = save_dir + "movie_frame_" + ("%06d" % frames[frame_val])
                    else:
                        file_name = args.output_filename + "_" + str(int(time_val))
                else:
                    if args.output_filename != None:
                        file_name = args.output_filename
                    else:
                        file_name = save_dir + "time_" + str(args.plot_time)
                
                if 0.0 in (cbar_min, cbar_max):
                    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.brg, rasterized=True, vmin=cbar_min, vmax=cbar_max)
                else:
                    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
                plt.gca().set_aspect('equal')
                
                import pdb
                pdb.set_trace()
