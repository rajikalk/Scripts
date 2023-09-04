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
    parser.add_argument("-abs", "--absolute_image", help="do you want to get the absolute value of the image field?", default="False")
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="Density")
    parser.add_argument("-f_unit", "--field_unit", help="What units would you like to plot the field?", default="g/cm**3")
    parser.add_argument("-div_by_thickness", "--divide_by_proj_thickness", help="Woudl you like to divide the field by the thickness of the projection?", default="True", type=str)
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 50., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default=0, type=int)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", type=int, default = 25)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-pzl", "--plot_z_velocities", help="do you want to plot the z velocity?", type=str, default='False')
    parser.add_argument("-vaf", "--velocity_annotation_frequency", help="how many velocity vectors do you want annotated across one side?", type=float, default=31.)
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-mt", "--movie_times", help="What movies times would you like plotted?", type=list, default=[])
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-16)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("-ic", "--image_center", help="where would you like to center the image?", type=int, default=0)
    parser.add_argument("-cc", "--calculation_center", help="where would you like to calculate center positionn and velocity?", type=int, default=0)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=250)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default='True')
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-proj_or", "--projection_orientation", help="Do you want to set the projection orientation? give as angle (in degrees) from positive y-axis", default=None, type=float)
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 300AU", type=float, default=500.)
    parser.add_argument("-use_disk", "--use_disk_angular_momentum", help="Do you want to use the disk angular momentum to define the normal vector for a projection?", default="False")
    parser.add_argument("-wf", "--weight_field", help="Do you want to have a weighted projection plot?", type=str, default=None)
    parser.add_argument("-use_gas", "--use_gas_center_calc", help="Do you want to use gas when calculating the center position adn veloity?", type=str, default='True')
    parser.add_argument("-all_files", "--use_all_files", help="Do you want to make frames using all available files instead of at particular time steps?", type=str, default='False')
    parser.add_argument("-update_alim", "--update_ax_lim", help="Do you want to update the axes limits by taking away the center position values or not?", type=str, default='False')
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-frames_only", "--make_frames_only", help="do you only want to make frames?", default='False', type=str)
    parser.add_argument("-use_L", "--use_angular_momentum", help="use disc angular momentum to find normal", default='False', type=str)
    parser.add_argument("-debug", "--debug_plotting", help="Do you want to debug why plotting is messing up", default='False', type=str)
    parser.add_argument("-res", "--resolution", help="define image resolution", default=4096, type=int)
    parser.add_argument("-active_rad", "--active_radius", help="within what radius of the centered sink do you want to consider when using sink and gas for calculations", type=float, default=10000.0)
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
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
    if args.axis == "xy":
        xlabel = r"$x$ (AU)"
        ylabel = r"$y$ (AU)"
    elif args.axis == "xz":
        xlabel = r"$x$ (AU)"
        ylabel = r"$z$ (AU)"
    elif args.axis == "yz":
        xlabel = r"$y$ (AU)"
        ylabel = r"$z$ (AU)"
    else:
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
if args.plot_time != None:
    m_times = [args.plot_time]
else:
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=args.presink_frames, end_time=args.end_time, form_time=sink_form_time)
    
no_frames = len(m_times)
m_times = m_times[args.start_frame:]

if args.make_frames_only == 'False':
    """
    except:
        files = files[:-1]
        ds = yt.load(files[-1], units_override=units_override)
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
    """
        
    sys.stdout.flush()
    CW.Barrier()

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
    if args.projection_orientation != None:
        y_val = 1./np.tan(np.deg2rad(args.projection_orientation))
        L = [1.0, y_val, 0.0]
    else:
        if args.axis == 'xy':
            L = [0.0, 0.0, 1.0]
        elif args.axis == 'xz':
            L = [0.0, 1.0, 0.0]
        elif args.axis == 'yz':
            L = [1.0, 0.0, 0.0]
    myf.set_normal(L)
    xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
    if args.ax_lim != None:
        xlim = [-1*args.ax_lim, args.ax_lim]
        ylim = [-1*args.ax_lim, args.ax_lim]
    x_width = (xlim[1] -xlim[0])
    y_width = (ylim[1] -ylim[0])
    thickness = yt.YTQuantity(args.slice_thickness, 'AU')
    #Sets center for calculating center position and velocity
    myf.set_center_pos_ind(args.image_center)

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
            pickle_file = save_dir + args.axis + '_' + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
        else:
            weight_field = args.weight_field
            pickle_file = save_dir + args.axis + '_' + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + ".pkl"
           
    sys.stdout.flush()
    CW.Barrier()

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

    sys.stdout.flush()
    CW.Barrier()

    if rank == 0:
        for sim_file in usable_files:
            file_int = file_int + 1
            try:
                if usable_files[file_int] == usable_files[file_int-1]:
                    os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
            except:
                continue
        print('Finished copying pickles that use the same file for the same frame')
    
    del usable_files
    del frames
sys.stdout.flush()
CW.Barrier()

#Copy pickle files:
if rank == 0 and args.plot_time == None:
    pickle_files = sorted(glob.glob(save_dir+"movie_frame_*.pkl"))
    file_counter = 0
    end_file_number = int(pickle_files[-1].split('_')[-1].split('.')[0]) + 1
    while file_counter < end_file_number:
        pickle_file = save_dir+"movie_frame_"+("%06d" % file_counter)+".pkl"
        if pickle_file not in pickle_files:
            os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % int(file_counter-1)) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % file_counter) + ".pkl ")
            print("copied file to "+ save_dir + "movie_frame_" + ("%06d" % file_counter) + ".pkl")
        file_counter = file_counter + 1

import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects

sys.stdout.flush()
CW.Barrier()

if args.plot_time != None:
    pickle_files = [save_dir + 'time_' + str(float(args.plot_time)) +'.pkl']
else:
    pickle_files = sorted(glob.glob(save_dir+"movie_frame_*.pkl"))
    start_file_int = pickle_files.index(save_dir+"movie_frame_"+("%06d" % args.start_frame)+".pkl")
    pickle_files = pickle_files[start_file_int:]

sys.stdout.flush()
CW.Barrier()

rit = args.working_rank - 1
for pickle_file in pickle_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        print('making frame from', pickle_file, 'on rank', rank)
        frame_no = int(pickle_file.split('_')[-1].split('.')[0])
        file = open(pickle_file, 'rb')
        X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()
        
        time_val = args_dict['time_val']
        
        if len(m_times) > 1:
            if args.output_filename == None:
                file_name = save_dir + "movie_frame_" + ("%06d" % frame_no)
            else:
                file_name = args.output_filename + "_" + str(int(time_val))
        else:
            if args.output_filename != None:
                file_name = args.output_filename
            else:
                file_name = save_dir + "time_" + str(args.plot_time)
        
        make_plot = True
        if os.path.isfile(file_name + ".jpg"):
            make_plot = False
            
        if make_plot == True:
            print(file_name + ".jpg" + " doesn't exist, so plotting image")
            #import pdb
            #pdb.set_trace()
            xlim = args_dict['xlim']
            ylim = args_dict['ylim']
            if np.round(np.mean(args_dict['xlim'])) != np.round(np.mean(X)):
                shift_x = np.mean(X)
                shift_y = np.mean(Y)
                X = X - shift_x
                Y = Y - shift_y
                X_vel = X_vel - shift_x
                Y_vel = Y_vel - shift_y
                part_info['particle_position'][0] = part_info['particle_position'][0] - shift_x
                part_info['particle_position'][1] = part_info['particle_position'][1] - shift_y
                  
            has_particles = args_dict['has_particles']
            xabel = args_dict['xabel']
            yabel = args_dict['yabel']
            plt.clf()
            fig, ax = plt.subplots()
            ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
            ax.set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if args.debug_plotting != 'False':
                plt.savefig("Test_776.jpg", format='jpg', bbox_inches='tight')
            
            if 0.0 in (cbar_min, cbar_max) or len(np.where(np.array([cbar_min, cbar_max]) < 0)[0]) > 0:
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.bwr, rasterized=True, vmin=cbar_min, vmax=cbar_max, zorder=1)
            else:
                cmap=plt.cm.gist_heat
                plot = ax.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            plt.gca().set_aspect('equal')
            if args.debug_plotting != 'False':
                plt.savefig("Test_784.jpg", format='jpg', bbox_inches='tight')
            if frame_no > 0 or time_val > -1.0:
                # plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            if args.debug_plotting != 'False':
                plt.savefig("Test_790.jpg", format='jpg', bbox_inches='tight')
            cbar = plt.colorbar(plot, pad=0.0)
            if args.debug_plotting != 'False':
                plt.savefig("Test_793.jpg", format='jpg', bbox_inches='tight')
            if args.plot_z_velocities == 'False':
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)#velz)
            else:
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)
            if args.debug_plotting != 'False':
                plt.savefig("Test_796.jpg", format='jpg', bbox_inches='tight')
                
            if has_particles:
                if args.ax_lim > 5000:
                    ax.scatter((part_info['particle_position'][0]), (part_info['particle_position'][1]), color='c', s=0.5)
                else:
                    if args.annotate_particles_mass == 'True':
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7)
                    else:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None, zorder=7)
            if args.debug_plotting != 'False':
                plt.savefig("Test_804.jpg", format='jpg', bbox_inches='tight')
            
            if 'Density' in simfo['field']:
                if args.divide_by_proj_thickness == "True":
                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
                else:
                    cbar.set_label(r"Column Density (g$\,$cm$^{-2}$)", rotation=270, labelpad=14, size=args.text_font)
            else:
                label_string = simfo['field'][1] + ' ($' + args.field_unit + '$)'
                cbar.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=args.text_font)
            if args.debug_plotting != 'False':
                plt.savefig("Test_812.jpg", format='jpg', bbox_inches='tight')

            if len(title) > 0:
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            if args.debug_plotting != 'False':
                plt.savefig("Test_818.jpg", format='jpg', bbox_inches='tight')

            plt.tick_params(axis='both', which='major')# labelsize=16)
            for line in ax.xaxis.get_ticklines():
                line.set_color('white')
            for line in ax.yaxis.get_ticklines():
                line.set_color('white')
            if args.debug_plotting != 'False':
                plt.savefig("Test_826.jpg", format='jpg', bbox_inches='tight')
            
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
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
                    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                except:
                    print("couldn't save for the dviread.py problem. Make frame " + str(frame_no) + " on ipython")
            else:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
                #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
        
sys.stdout.flush()
CW.Barrier()

print("completed making movie frames on rank", rank)

