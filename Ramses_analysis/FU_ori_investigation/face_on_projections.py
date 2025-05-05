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
from matplotlib import ticker
from scipy.ndimage import gaussian_filter
plt.rcParams["text.usetex"] = False

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="Density")
    parser.add_argument("-f_unit", "--field_unit", help="What units would you like to plot the field?", default="g/cm**-3")
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
    parser.add_argument("-use_part_for_vel", "--use_particle_for_center_vel_calc", help="Do you want to use the particles to calculate center velocity?", type=str, default='True')
    parser.add_argument("-all_files", "--use_all_files", help="Do you want to make frames using all available files instead of at particular time steps?", type=str, default='False')
    parser.add_argument("-update_alim", "--update_ax_lim", help="Do you want to update the axes limits by taking away the center position values or not?", type=str, default='False')
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-frames_only", "--make_frames_only", help="do you only want to make frames?", default='False', type=str)
    parser.add_argument("-debug", "--debug_plotting", help="Do you want to debug why plotting is messing up", default='False', type=str)
    parser.add_argument("-res", "--resolution", help="define image resolution", default=800, type=int)
    parser.add_argument("-active_rad", "--active_radius", help="within what radius of the centered sink do you want to consider when using sink and gas for calculations", type=float, default=10)

    parser.add_argument("-threshold", "--density_threshold", help="What number density threshold would you like to use?", type=float, default=0.0)
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("-sm", "--skip_made", help="do you want to skip frames already made?", type=str, default='True')
    parser.add_argument("-spec_file", "--specific_file", help="Do you want to use a specific file", type=str, default=None)
    parser.add_argument("-conv", "--convolve", help="Do you want to convolve the image with a gaussian beam?", type=str, default='True')
    parser.add_argument("-transpose", "--transpose_image", help="Do you want to transport image?", type=str, default='False')
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
    unit_string = str(dd[field].in_cgs().units)
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

#units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
#scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
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
        pickle_file = save_dir + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
    else:
        weight_field = args.weight_field
        pickle_file = save_dir + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + "/"
       
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

if rank == 0:
    print("Finding times")
if args.plot_time != None:
    m_times = [args.plot_time]
elif len(args.movie_times) > 0:
    m_times = eval(args.movie_times)
    no_frames = len(m_times)
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
    if args.specific_file == None:
        usable_files = mym.find_files(m_times, files, sink_form_time,sink_id, verbatim=False)
    else:
        usable_files = [args.specific_file]
    if rank == 0:
        if size == 0:
            import pdb
            pdb.set_trace()
        print("Usable_files=", usable_files)
    del sink_form_time
    del files
    
sys.stdout.flush()
CW.Barrier()

if args.make_frames_only == 'False':
    #Trying yt parallelism
    file_int = -1
    for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size)):
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
        if args.plot_time is None:
            pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[file_int])
        else:
            pickle_file = save_dir + "time_" + (str(int(args.plot_time)))
        if os.path.isfile(pickle_file) == False:
            ds = yt.load(fn, units_override=units_override)
            dd = ds.all_data()
            
            has_particles = has_sinks(ds)
            part_info = mym.get_particle_data(ds, sink_id=sink_id)
            
            time_val = m_times[file_int]
            
            center_pos = dd['Center_Position'].in_units('AU')
            center_vel = dd['Center_Velocity'].in_units('cm/s')
            
            measuring_sphere_primary = ds.sphere(center_pos, yt.YTQuantity(args.active_radius, 'au'))
            cell_mass = measuring_sphere_primary['mass'].in_units('g')
            
            sph_dx = measuring_sphere_primary['x'].in_units('cm') - center_pos[0].in_units('cm')
            sph_dy = measuring_sphere_primary['y'].in_units('cm') - center_pos[1].in_units('cm')
            sph_dz = measuring_sphere_primary['z'].in_units('cm') - center_pos[2].in_units('cm')
            
            sph_radial_vector = yt.YTArray([sph_dx, sph_dy, sph_dz]).T
            
            sph_dvx = measuring_sphere_primary['velocity_x'].in_units('cm/s') - center_vel[0].in_units('cm/s')
            sph_dvy = measuring_sphere_primary['velocity_y'].in_units('cm/s') - center_vel[1].in_units('cm/s')
            sph_dvz = measuring_sphere_primary['velocity_z'].in_units('cm/s') - center_vel[2].in_units('cm/s')
            
            sph_radial_vector = yt.YTArray([sph_dvx, sph_dvy, sph_dvz]).T
            
            L_x = cell_mass*(sph_dy*sph_dvz - sph_dz*sph_dvy)
            L_y = cell_mass*(sph_dx*sph_dvz - sph_dz*sph_dvx)
            L_z = cell_mass*(sph_dx*sph_dvy - sph_dy*sph_dvx)
            
            L_tot = yt.YTArray([np.sum(L_x), np.sum(L_y), np.sum(L_z)])
            L_mag = np.sqrt(np.sum(L_tot**2))
            L_unit = L_tot/L_mag
            
            field = simfo['field']
            proj = yt.OffAxisProjectionPlot(ds, L_unit, simfo['field'], width=(x_width, 'AU'), weight_field=weight_field, method='integrate', center=(center_pos.value, 'AU'), depth=(args.slice_thickness, 'AU'))#data_source=rv_cut_region
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
                    
            image = proj_array
            
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
            args_dict.update({'proj_vector': L_unit})
            
            pickle_file = pickle_file + '.pkl'
            file = open(pickle_file, 'wb')
            pickle.dump((X, Y, image, None, X_vel, Y_vel, None, None, part_info, args_dict, simfo, None), file)
            file.close()
            print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
sys.stdout.flush()
CW.Barrier()

#Section to plot figures:
print("Finished generating projection pickles")
if args.plot_time is None:
    pickle_files = sorted(glob.glob(save_dir + "movie_frame_*.pkl"))
else:
    pickle_files = [save_dir + "time_" + (str(int(args.plot_time)))+".pkl"]
#if rank==0:
#    print("pickle_files =", pickle_files)
rit = -1
for pickle_file in pickle_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        if len(m_times) > 1:
            if args.output_filename == None:
                file_name = pickle_file.split('.pkl')[0]
            else:
                file_name = args.output_filename + "_" + str(int(time_val))
        else:
            if args.output_filename != None:
                file_name = args.output_filename
            else:
                file_name = save_dir + "time_" + str(int(args.plot_time))
        if os.path.exists(file_name + ".jpg") and args.skip_made == 'True':
            if os.stat(file_name + ".jpg").st_size == 0:
                make_frame = True
            else:
                make_frame = False
                print(file_name + ".jpg already exists, so skipping")
        else:
            make_frame = True
        if make_frame:
            proj_number = pickle_file.split('.pkl')[0][-1]
            print("on rank", rank, "using pickle_file", pickle_file)
            file = open(pickle_file, 'rb')
            X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
            
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
            #title_str = "proj vec:["+str(np.round(args_dict['proj_vector'][0].value*100)/100)+ "," + str(np.round(args_dict['proj_vector'][1].value*100)/100) +"," +str(np.round(args_dict['proj_vector'][2].value*100)/100)+"], companion LOS pos:"+str(int(part_info['particle_position_z'][1]))+"AU"
            ax.set_title("Time:"+str(time_val)+"yr")
            
            if None in (cbar_min, cbar_max):
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.magma, rasterized=True)
            elif 0.0 in (cbar_min, cbar_max):
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.magma, rasterized=True, vmin=cbar_min, vmax=cbar_max)
            else:
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.magma, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
            plt.gca().set_aspect('equal')
            cbar = plt.colorbar(plot, pad=0.0)
            
            #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)

            if has_particles:
                if size > 1:
                    args.annotate_particles_mass = False
                if args.annotate_particles_mass == True:
                    mym.annotate_particles(ax, np.array([[0], [0]]), part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
                else:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)

            if 'Density' in simfo['field']:
                if args.divide_by_proj_thickness == "True":
                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
                else:
                    cbar.set_label(r"Column Density (g$\,$cm$^{-2}$)", rotation=270, labelpad=14, size=args.text_font)

            if len(title) > 0:
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

            plt.tick_params(axis='both', which='major')# labelsize=16)
            for line in ax.xaxis.get_ticklines():
                line.set_color('white')
            for line in ax.yaxis.get_ticklines():
                line.set_color('white')

            if args.annotate_time == "True":
                #time_string = "t="+str(int(time_val))+"yr"
                #ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string, va="center", ha="left", color='w', fontsize=args.text_font)
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    #time_string = "t="+str(int(time_val))+"yr"
                    time_string = "$t$="+str(int(time_val))+"yr"
                    time_string_raw = r"{}".format(time_string)
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string, va="center", ha="left", color='w', fontsize=args.text_font)
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
                    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    print('Created frame of projection', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                except:
                    print("couldn't save for the dviread.py problem. Make frame " + str(proj_number) + " on ipython")
            else:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame of projection', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')

print("Completed making frames on rank", rank)
