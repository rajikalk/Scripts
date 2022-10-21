#!/usr/bin/env python
import h5py
import numpy as np
#import numpy.ma as ma
#from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import csv
import sys
import os
from matplotlib import transforms
import glob
import my_module as mym
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.patheffects as path_effects
import yt
import my_fields as myf
from PIL import Image

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-abs", "--absolute_image", help="do you want to get the absolute value of the image field?", default="False")
    parser.add_argument("-zt", "--zoom_times", help="0 is default zoom", default=0, type=float)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    #parser.add_argument("-st", "--start_time", help="What time would you like to start calculating times from?", type=float, default=0.0)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", type=int, default = 25)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-vaf", "--velocity_annotation_frequency", help="how many velocity vectors do you want annotated across one side?", type=float, default=31.)
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-mt", "--movie_times", help="What movies times would you like plotted?", type=list, default=[])
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-16')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("-ic", "--image_center", help="where would you like to center the image?", type=int, default=0)
    parser.add_argument("-cc", "--calculation_center", help="where would you like to center the image?", type=int, default=0)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=None)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-yt", "--yt_proj", help="Do you want to use yt to create projections as opposed to the movie files?", default=False)
    parser.add_argument("-proj_or", "--projection_orientation", help="Do you want to set the projection orientation? give as angle (in degrees) from positive y-axis", default=None, type=float)
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 300AU", type=float, default=300.)
    parser.add_argument("-use_disk", "--use_disk_angular_momentum", help="Do you want to use the disk angular momentum to define the normal vector for a projection?", default="False")
    parser.add_argument("-wf", "--weight_field", help="Do you want to have a weighted projection plot?", type=str, default=None)
    parser.add_argument("-db", "--debug", help="Wanting to use the debugger where you inserted it?", type=str, default='False')
    parser.add_argument("-use_gas", "--use_gas_center_calc", help="Do you want to use gas when calculating the center position adn veloity?", type=str, default='True')
    parser.add_argument("-all_files", "--use_all_files", help="Do you want to make frames using all available files instead of at particular time steps?", type=str, default='False')
    parser.add_argument("-update_alim", "--update_ax_lim", help="Do you want to update the axes limits bu taking away the center position values or not?", type=str, default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def get_files(path, args):
    '''
    Gets movie files (proj or slice) files if possible, if not it find the plt files
    '''
    if args.yt_proj == False:
        if args.axis == "xz":
            type = "proj"
        else:
            type = "slice"
        source_directory = sorted(glob.glob(path + 'WIND_' + type + '*'))[:-1]
    else:
        source_directory = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))[:-1]
    return source_directory

def has_sinks(file):
    '''
    Checks particle file to see if particles exists, or tries the plot file.
    '''
    try:
        part_file = file[:-12] + 'part' + file[-5:]
        f = h5py.File(part_file, 'r')
        if f[list(f.keys())[1]][-1][-1] > 0:
            f.close()
            return True
        else:
            f.close()
            return False
    except:
        f = h5py.File(file, 'r')
        if "particlepositions" in list(f.keys()):
            f.close()
            return True
        else:
            f.close()
            return False

def sim_info(path, file, args):
    """
    Finds particle info, relevant to frame size and such. NOTE ACCRETION RADIUS IS GIVEN FROM PARTICLE INFO FUNCTION
    """
    try:
        f = h5py.File(file, 'r')
        for key in list(f.keys()):
            if args.field in key:
                field = key
        dim = np.shape(f[field])[0]
        if args.zoom_times > 0:
            zoom_cell = int((dim - dim/float(args.zoom_times))/2.)
        else:
            zoom_cell = 0
        xmin = f['minmax_xyz'][0][0]/yt.units.au.in_units('cm').value
        xmax = f['minmax_xyz'][0][1]/yt.units.au.in_units('cm').value
        xmin_full = xmin
        cl = (xmax-xmin)/dim
        cell_positions = np.arange(xmin, xmax-1, cl)
        xmin = f['minmax_xyz'][0][0]/yt.units.au.in_units('cm').value + zoom_cell*cl
        xmax = f['minmax_xyz'][0][1]/yt.units.au.in_units('cm').value - zoom_cell*cl
        if args.axis == "xy":
            ymin = f['minmax_xyz'][1][0]/yt.units.au.in_units('cm').value + zoom_cell*cl
            ymax = f['minmax_xyz'][1][1]/yt.units.au.in_units('cm').value - zoom_cell*cl
        else:
            ymin = f['minmax_xyz'][2][0]/yt.units.au.in_units('cm').value + zoom_cell*cl
            ymax = f['minmax_xyz'][2][1]/yt.units.au.in_units('cm').value - zoom_cell*cl
        f.close()
        annotate_freq = ((xmax/cl) - (xmin/cl))/args.velocity_annotation_frequency
    except:
        f = h5py.File(file, 'r')
        f.close()
        if args.field == 'dens':
            field = ('flash', 'dens')
        else:
            part_file = file[:-12] + 'part' + file[-5:]
            f = yt.load(file, particle_filename=part_file)
            field = f.derived_field_list[[x[1] for x in f.derived_field_list].index(args.field)]
            f.close()
        dim = 800
        zoom_cell = 0.0
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
        xmin_full = xmin
        annotate_freq = dim/args.velocity_annotation_frequency
    smoothing = annotate_freq/2
    if args.axis == "xz":
        type = "proj"
    else:
        type = "slice"
    sim_info = {'field': field,
                'dimension': dim,
                'zoom_cell': zoom_cell,
                'movie_file_type': type,
                'xmin': xmin,
                'xmax': xmax,
                'ymin': ymin,
                'ymax': ymax,
                'cell_length': cl,
                'annotate_freq': annotate_freq,
                'smoothing': smoothing,
                'xmin_full':  xmin_full
                }
    f.close()
    return sim_info

def get_image_arrays(f, field, simfo, args, X, Y):
    '''
    Gets image array from MOVIE files, not plot files
    '''
    dim = int(simfo['dimension'])
    image = []
    xpos = int(np.round(np.min(X) - simfo['xmin_full'])/simfo['cell_length'])
    ypos = int(np.round(np.min(Y) - simfo['xmin_full'])/simfo['cell_length'])
    print("XPOS =", xpos)
    for x in range(ypos, ypos+len(X[0])):
        image_val = f[field][x]
        if np.shape(image_val)[0] == 1:
            image_val = image_val.transpose()
        image_val = image_val[xpos:xpos+len(X[0])]
        if simfo['movie_file_type'] == "proj":
            image_val = image_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
        image.append(image_val)
    image = np.array(image)
    if np.shape(image)[-1] == 1:
        image = image[:, :, 0]
    return image

def image_properties(X, Y, args, sim_info):
    if args.axis == "xy":
        ylabel = '$y$ (AU)'
    else:
        ylabel = '$z$ (AU)'
    if args.yt_proj == False:
        xlabel = '$x$ (AU)'
    else:
        if args.axis == "xy":
            xlabel = '$x$ (AU)'
        else:
            xlabel = 'Distance from center (AU)'
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    return xlabel, ylabel, xlim, ylim

def rainbow_text(x,y,ls,lc,**kw):
    """
    Take a list of strings ``ls`` and colors ``lc`` and place them next to each
    other, with text ls[i] being shown in color lc[i]
    
    This example shows how to do both vertical and horizontal text, and will
    pass all keyword arguments to plt.text, so you can set the font size,
    family, etc.
    """
    t = plt.gca().transData
    fig = plt.gcf()
    plt.show()
    
    #horizontal version
    for s,c in zip(ls,lc):
        text = plt.text(x,y," "+s+" ",color=c, transform=t, **kw)
        text.draw(fig.canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(text._transform, x=ex.width, units='dots')

#=======MAIN=======
def main():
    
    rank = CW.Get_rank()
    size = CW.Get_size()
    #comm = MPI.COMM_WORLD
    
    # Read in directories:
    path = sys.argv[1]
    save_dir = sys.argv[2]
    if rank == 0 and os.path.exists(save_dir) == False:
        os.makedirs(save_dir)
    
    args = parse_inputs()
    
    mym.set_global_font_size(args.text_font)
    files = get_files(path, args)
    simfo = sim_info(path, files[0], args)
    if args.yt_proj == False:
        X, Y, X_vel, Y_vel, cl = mym.initialise_grid(files[-1], zoom_times=args.zoom_times, num_of_vectors=args.velocity_annotation_frequency)
        L=None
    else:
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
            L = [0.0, 0.0, 1.0]
        myf.set_normal(L)
        print("SET PROJECTION ORIENTATION L=", myf.get_normal())
    if args.yt_proj == False and args.image_center != 0:
        sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
    elif args.yt_proj != False and args.image_center != 0:
        sim_files = files
    xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
    if args.ax_lim != None:
        xlim = [-1*args.ax_lim, args.ax_lim]
        ylim = [-1*args.ax_lim, args.ax_lim]
    #Sets center for calculating center position and velocity
    myf.set_center(args.image_center)
    
    #Sets whether to use gas for center posiiton and velocity calculation
    if args.use_gas_center_calc == 'True':
        myf.set_use_gas(True)
    else:
        myf.set_use_gas(False)
   
    title_parts = args.title.split('_')
    title = ''
    for part in title_parts:
        if part != title_parts[-1]:
            title = title + part + ' '
        else:
            title = title + part
    
    # Initialise Grid and build lists
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=args.presink_frames, end_time=args.end_time)
    sys.stdout.flush()
    CW.Barrier()

    if args.use_all_files == 'False':
        no_frames = len(m_times)
        m_times = m_times[args.start_frame:]
        usable_files = mym.find_files(m_times, files)
        frames = list(range(args.start_frame, no_frames))
    elif args.use_all_files != 'False' and args.plot_time != None:
        usable_files = mym.find_files([args.plot_time], files)
        start_index = files.index(usable_files[0])
        args.plot_time = None
        end_file = mym.find_files([args.end_time], files)
        end_index = files.index(end_file[0])
        usable_files = files[start_index:end_index]
        frames = list(range(len(usable_files)))
        no_frames = len(usable_files)
    else:
        start_file = mym.find_files([0], files)
        usable_files = files[files.index(start_file[0]):]
        frames = list(range(len(usable_files)))
        no_frames = len(usable_files)
    if args.image_center != 0 and args.yt_proj == False:
        usable_sim_files = mym.find_files(m_times, sim_files)
        del sim_files
    sys.stdout.flush()
    CW.Barrier()

    # Define colourbar bounds
    cbar_max = args.colourbar_max
    try:
        cbar_min = float(args.colourbar_min)
    except:
        cbar_min = float(args.colourbar_min[1:])
    
    sys.stdout.flush()
    CW.Barrier()

    if args.yt_proj:
        thickness = yt.YTArray(args.slice_thickness, 'AU')
        if args.plot_time != None:
            if args.weight_field == 'None':
                weight_field = None
                pickle_file = save_dir + args.axis + '_' + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
            else:
                weight_field = args.weight_field
                pickle_file = save_dir + args.axis + '_' + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + ".pkl"
        rit = args.working_rank - 1
        ts = yt.load(usable_files)#, parallel=1)
        file_int = -1
        #yt.enable_parallelism()
        #ds = yt.load(usable_files[0])
        #dd = ds.all_data()
        #print("loaded initial file")
        #center_pos = dd['Center_Position'].in_units('au').value
        #for file_int in range(len(usable_files)):
        for dataset in ts.piter():
            ds = dataset
            if size > 1:
                file_int = usable_files.index(path + str(ds))
            else:
                file_int = file_int + 1
                if usable_files[file_int] == usable_files[file_int-1]:
                    os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
            #rit = rit + 1
            #print('rit=', rit, ', file_int=', file_int)
            #if rit == size:
            #    rit = 0
            #rit = rank
            #if rank == rit:
            #ds = yt.load(usable_files[file_int])
            if args.plot_time is None:
                pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl"
            make_pickle = False
            if os.path.isfile(pickle_file) == False:
                make_pickle = True
            elif os.path.isfile(pickle_file) == True:
                if os.stat(pickle_file).st_size == 0:
                    make_pickle = True
            if make_pickle == True:
                sys.stdout.flush()
                CW.Barrier()
                print("PICKLE:", pickle_file,"DOESN'T EXIST. MAKING PROJECTION  FOR FRAME", frames[file_int], "ON RANK", rank, "USING FILE", ds)
                sys.stdout.flush()
                CW.Barrier()
                    
                has_particles = has_sinks(path + str(ds)) #(usable_file)#(path + str(ds))
                if has_particles:
                    part_info = mym.get_particle_data(path + str(ds), args.axis, proj_or=L)
                else:
                    part_info = {}
                    #part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2]
                    
                '''
                if args.ax_lim != None:
                    if has_particles and args.image_center != 0:
                        xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center - 1], args.ax_lim + part_info['particle_position'][0][args.image_center - 1]]
                        ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center - 1], args.ax_lim + part_info['particle_position'][1][args.image_center - 1]]
                '''
                
                x_width = (xlim[1] -xlim[0])
                y_width = (ylim[1] -ylim[0])
                thickness = yt.YTQuantity(args.slice_thickness, 'AU')
                '''
                if args.image_center != 0 and has_particles:
                    original_positions = [X, Y, X_vel, Y_vel]
        
                    #x_pos = np.round(part_info['particle_position'][0][args.image_center - 1]/simfo['cell_length'])*simfo['cell_length']
                    #y_pos = np.round(part_info['particle_position'][1][args.image_center - 1]/simfo['cell_length'])*simfo['cell_length']
                    x_pos = 0
                    y_pos = 0
                    X = X + x_pos
                    Y = Y + y_pos
                    X_vel = X_vel + x_pos
                    Y_vel = Y_vel + y_pos
                '''

                print("in time series loop")
                #ds = yt.load(usable_file)
                dd = ds.all_data()
                print("loaded all data")
                
                try:
                    time_val = m_times[file_int]
                except:
                    sink_creation_time = np.min(dd['particle_creation_time'].value)
                    time_real = yt.YTQuantity(ds.current_time.value - sink_creation_time, 's')
                    time_val = np.round(time_real.in_units('yr'))
                
                center_pos = dd['Center_Position'].in_units('au').value
                print('Center Pos=' + str(center_pos))
                
                #Update X and Y to be centered on center position
                if args.update_ax_lim == 'True':
                    if args.axis == 'xy':
                        X_image = X + center_pos[0]
                        Y_image = Y + center_pos[1]
                        X_image_vel = X_vel + center_pos[0]
                        Y_image_vel = Y_vel + center_pos[1]
                    else:
                        X_image = X + center_pos[0]
                        Y_image = Y + center_pos[2]
                        X_image_vel = X_vel + center_pos[0]
                        Y_image_vel = Y_vel + center_pos[2]
                else:
                    if args.axis == 'xy':
                        X_image = X
                        Y_image = Y
                        X_image_vel = X_vel
                        Y_image_vel = Y_vel
                    else:
                        X_image = X
                        Y_image = Y
                        X_image_vel = X_vel
                        Y_image_vel = Y_vel
                
                if has_particles:
                    if args.update_ax_lim == 'False':
                        part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0]
                        if args.axis == 'xy':
                            part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[1]
                        else:
                            part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2]
                
                print("initialised fields")
                myf.set_center(args.calculation_center)
                myf.set_use_gas(True)
                center_vel = dd['Center_Velocity'].value
                if args.use_gas_center_calc == 'True':
                    myf.set_use_gas(True)
                else:
                    myf.set_use_gas(False)
                myf.set_center(args.image_center)
                print("center_vel =", center_vel, "on rank", rank, "for", ds)
                
                if args.axis == 'xy':
                    center_vel=center_vel[:2]
                    if args.use_disk_angular_momentum != "False":
                        disk = ds.disk(center_pos, L, (args.ax_lim*2, 'au'), (args.slice_thickness/2., 'au'))
                        tot_vec = [np.sum(disk['Angular_Momentum_x']).value, np.sum(disk['Angular_Momentum_y']).value, np.sum(disk['Angular_Momentum_z']).value]
                        tot_mag = np.sqrt(tot_vec[0]**2. + tot_vec[1]**2. + tot_vec[2]**2.)
                        L = tot_vec/tot_mag
                        print("SET PROJECTION ORIENTATION L=", L)
                else:
                    center_vel=center_vel[::2]
                
                if args.axis == 'xy':
                    axis_ind = 2
                    left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[1]-(0.75*y_width), center_pos[2]-(0.5*args.slice_thickness)], 'AU')
                    right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[1]+(0.75*y_width), center_pos[2]+(0.5*args.slice_thickness)], 'AU')
                    region = ds.box(left_corner, right_corner)
                else:
                    axis_ind = 1
                    #left_corner = yt.YTArray([center_pos[2]-(0.5*args.slice_thickness), center_pos[0]-(0.75*x_width), center_pos[1]-(0.75*y_width)], 'AU')
                    #right_corner = yt.YTArray([center_pos[2]+(0.5*args.slice_thickness), center_pos[0]+(0.75*x_width), center_pos[1]+(0.75*y_width)], 'AU')
                    left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[2]-(0.5*args.slice_thickness),  center_pos[1]-(0.75*y_width)], 'AU')
                    right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[2]+(0.5*args.slice_thickness), center_pos[1]+(0.75*y_width)], 'AU')
                    region = ds.box(left_corner, right_corner)
                 #(usable_file)#(path + str(ds))
        
                #proj.set_buff_size(1024)
                
                #Set weight field
                weight_field = args.weight_field
                
                if weight_field == None and args.projection_orientation == None:
                    if args.axis == 'xy':
                        proj = yt.ProjectionPlot(ds, axis_ind, [simfo['field'], 'Corrected_velx', 'Corrected_vely', 'magx', 'magy'], width=(x_width,'au'), weight_field=weight_field, data_source=region, method='integrate', center=dd['Center_Position'].in_units('cm'))
                        proj_depth = yt.ProjectionPlot(ds, axis_ind, ['z', 'Neg_z', 'dz', 'Neg_dz'], width=(x_width,'au'), weight_field=None, data_source=region, method='mip', center=dd['Center_Position'].in_units('cm'))
                        image = proj.frb.data[simfo['field']].value/(proj_depth.frb.data[('gas', 'Neg_z')].in_units('cm') + proj_depth.frb.data[('gas', 'z')].in_units('cm')).value
                        print("IMAGE VALUES ARE BETWEEN:", np.min(image), np.max(image))
                        velx_full = proj.frb.data[('gas', 'Corrected_velx')].in_units('cm**2/s')/((proj_depth.frb.data[('gas', 'Neg_z')].in_units('cm') + proj_depth.frb.data[('gas', 'Neg_dz')].in_units('cm')/2.) + (proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'dz')].in_units('cm')/2.))
                        print("IMAGE VELX ARE BETWEEN:", np.min(velx_full), np.max(velx_full))
                        vely_full = proj.frb.data[('gas', 'Corrected_vely')].in_units('cm**2/s')/((proj_depth.frb.data[('gas', 'Neg_z')].in_units('cm') + proj_depth.frb.data[('gas', 'Neg_dz')].in_units('cm')/2.) + (proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'dz')].in_units('cm')/2.))
                        print("IMAGE VELY ARE BETWEEN:", np.min(vely_full), np.max(vely_full))
                        magx = proj.frb.data[('flash', 'magx')].in_units('cm*gauss')/((proj_depth.frb.data[('gas', 'Neg_z')].in_units('cm') + proj_depth.frb.data[('gas', 'Neg_dz')].in_units('cm')/2.) + (proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'dz')].in_units('cm')/2.))
                        print("IMAGE MAGX ARE BETWEEN:", np.min(magx), np.max(magx))
                        magy = proj.frb.data[('flash', 'magy')].in_units('cm*gauss')/((proj_depth.frb.data[('gas', 'Neg_z')].in_units('cm') + proj_depth.frb.data[('gas', 'Neg_dz')].in_units('cm')/2.) + (proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'dz')].in_units('cm')/2.))
                        print("IMAGE MAGY ARE BETWEEN:", np.min(magy), np.max(magy))
                    else:
                        proj = yt.ProjectionPlot(ds, axis_ind, [simfo['field'], 'velx', 'velz', 'magx', 'magz'], width=(x_width,'au'), weight_field=weight_field, data_source=region, method='integrate', center=dd['Center_Position'].in_units('cm'))
                        image = proj.frb.data[simfo['field']].T/thickness.in_units('cm')
                        print("IMAGE VALUES ARE BETWEEN:", np.min(image), np.max(image))
                        velx_full = proj.frb.data[('flash', 'velx')].T.in_units('cm**2/s')/thickness.in_units('cm')
                        print("IMAGE VELX ARE BETWEEN:", np.min(velx_full), np.max(velx_full))
                        vely_full = proj.frb.data[('flash', 'velz')].T.in_units('cm**2/s')/thickness.in_units('cm')
                        print("IMAGE VELY ARE BETWEEN:", np.min(vely_full), np.max(vely_full))
                        magx = proj.frb.data[('flash', 'magx')].T.in_units('cm*gauss')/thickness.in_units('cm')
                        print("IMAGE MAGX ARE BETWEEN:", np.min(magx), np.max(magx))
                        magy = proj.frb.data[('flash', 'magz')].T.in_units('cm*gauss')/thickness.in_units('cm')
                        print("IMAGE MAGY ARE BETWEEN:", np.min(magy), np.max(magy))
                    
                    velx_full = np.array(velx_full.value)
                    vely_full = np.array(vely_full.value)
                    magx = np.array(magx.value)
                    magy = np.array(magy.value)
                elif args.projection_orientation == None:
                    image = proj.frb.data[simfo['field']].value
                    velx_full = proj.frb.data[('gas', 'Corrected_velx')].in_units('cm/s').value
                    vely_full = proj.frb.data[('gas', 'Corrected_vely')].in_units('cm/s').value
                    magx = proj.frb.data[('flash', 'magx')].in_units('gauss').value
                    magy = proj.frb.data[('flash', 'magy')].in_units('gauss').value
                    
                elif args.projection_orientation != None:
                    
                    #part_info = mym.get_particle_data(path + str(ds), args.axis, proj_or=L)
                    #part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0]
                    #part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[2]
                    #print("particle_pos=", part_info['particle_position'])
                    
                    proj = yt.OffAxisProjectionPlot(ds, L, [simfo['field'], 'Projected_Velocity', 'velz', 'Projected_Magnetic_Field', 'magz'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'), weight_field=weight_field)
                    if weight_field == None:
                        image = proj.frb.data[simfo['field']].value/thickness.in_units('cm').value
                        velx_full = proj.frb.data[('gas', 'Projected_Velocity')].in_units('cm**2/s').value/thickness.in_units('cm').value
                        vely_full = proj.frb.data[('flash', 'velz')].in_units('cm**2/s').value/thickness.in_units('cm').value
                        magx = proj.frb.data[('gas', 'Projected_Magnetic_Field')].in_units('cm*gauss').value/thickness.in_units('cm').value
                        magy = proj.frb.data[('flash', 'magz')].in_units('cm*gauss').value/thickness.in_units('cm').value
                    else:
                        image = proj.frb.data[simfo['field']].value
                        velx_full = proj.frb.data[('gas', 'Projected_Velocity')].in_units('cm/s').value
                        vely_full = proj.frb.data[('flash', 'velz')].in_units('cm/s').value
                        magx = proj.frb.data[('gas', 'Projected_Magnetic_Field')].in_units('gauss').value
                        magy = proj.frb.data[('flash', 'magz')].in_units('gauss').value
                else:
                    if args.use_disk_angular_momentum == "False":
                        #if has_particles:
                        #    part_info = mym.get_particle_data(path + str(ds), args.axis, proj_or=L)
                        #else:
                        #    part_info = {}
                        
                        #part_info = mym.get_particle_data(path + str(ds), args.axis, proj_or=L)
                        #if args.image_center == 0:
                        #part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0]
                        #part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[1]
                    
                        print("particle_pos=", part_info['particle_position'])
                        
                        left_corner = yt.YTArray([center_pos[0]-(0.75*x_width), center_pos[1]-(0.75*y_width), center_pos[2]-(0.5*args.slice_thickness)], 'AU')
                        right_corner = yt.YTArray([center_pos[0]+(0.75*x_width), center_pos[1]+(0.75*y_width), center_pos[2]+(0.5*args.slice_thickness)], 'AU')
                        region = ds.box(left_corner, right_corner)
                        proj = yt.ProjectionPlot(ds, 2, [simfo['field'], 'velx', 'vely', 'magx', 'magy'], width=(x_width,'au'), weight_field=weight_field, data_source=region, method='integrate')
                        if weight_field == None:
                            proj_depth = yt.ProjectionPlot(ds, 2, ['z', 'Neg_z'], width=(x_width,'au'), weight_field=None, data_source=region, method='mip')
                            image = proj.frb.data[simfo['field']].value/(proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'z')].in_units('cm')).value
                            print("IMAGE VALUES ARE BETWEEN:", np.min(image), np.max(image))
                            velx_full = proj.frb.data[('flash', 'velx')].in_units('cm**2/s').value/(proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'z')].in_units('cm')).value
                            print("IMAGE VELX ARE BETWEEN:", np.min(velx_full), np.max(velx_full))
                            vely_full = proj.frb.data[('flash', 'vely')].in_units('cm**2/s').value/(proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'z')].in_units('cm')).value
                            print("IMAGE VELY ARE BETWEEN:", np.min(vely_full), np.max(vely_full))
                            magx = proj.frb.data[('flash', 'magx')].in_units('cm*gauss')/(proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'z')].in_units('cm'))
                            print("IMAGE MAGX ARE BETWEEN:", np.min(magx), np.max(magx))
                            magy = proj.frb.data[('flash', 'magy')].in_units('cm*gauss')/(proj_depth.frb.data[('gas', 'z')].in_units('cm') + proj_depth.frb.data[('gas', 'z')].in_units('cm'))
                            print("IMAGE MAGY ARE BETWEEN:", np.min(magy), np.max(magy))
                        else:
                            image = proj.frb.data[simfo['field']].value
                            velx_full = proj.frb.data[('flash', 'velx')].in_units('cm/s').value
                            vely_full = proj.frb.data[('flash', 'vely')].in_units('cm/s').value
                            magx = proj.frb.data[('flash', 'magx')].in_units('gauss').value
                            magy = proj.frb.data[('flash', 'magy')].in_units('gauss').value
                    else:

                        #part_info = mym.get_particle_data(path + str(ds), args.axis, proj_or=L)
                        #if args.image_center == 0:
                        #    part_info['particle_position'][0] = part_info['particle_position'][0] - center_pos[0]
                        #    part_info['particle_position'][1] = part_info['particle_position'][1] - center_pos[1]
                        #else:
                        #    part_info['particle_position'][0] = part_info['particle_position'][0] - part_info['particle_position'][0][args.image_center-1]
                        #    part_info['particle_position'][1] = part_info['particle_position'][1] - part_info['particle_position'][1][args.image_center-1]
                        print("particle_pos=", part_info['particle_position'])
                        
                        proj = yt.OffAxisProjectionPlot(ds, L, [simfo['field'], 'velx', 'vely', 'magx', 'magy'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'), weight_field=weight_field)
                        if weight_field == None:
                            image = proj.frb.data[simfo['field']].value/thickness.in_units('cm').value
                            velx_full = proj.frb.data[('flash', 'velx')].in_units('cm**2/s').value/thickness.in_units('cm').value
                            vely_full = proj.frb.data[('flash', 'vely')].in_units('cm**2/s').value/thickness.in_units('cm').value
                            magx = proj.frb.data[('flash', 'magx')].in_units('cm*gauss').value/thickness.in_units('cm').value
                            magy = proj.frb.data[('flash', 'magy')].in_units('cm*gauss').value/thickness.in_units('cm').value
                        else:
                            image = proj.frb.data[simfo['field']].value
                            velx_full = proj.frb.data[('flash', 'velx')].in_units('cm/s').value
                            vely_full = proj.frb.data[('flash', 'vely')].in_units('cm/s').value
                            magx = proj.frb.data[('flash', 'magx')].in_units('gauss').value
                            magy = proj.frb.data[('flash', 'magy')].in_units('gauss').value
                
                velx, vely = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=center_vel)
                del velx_full
                del vely_full

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
            
        print('FINISHED MAKING YT PROJECTIONS')
                
    sys.stdout.flush()
    CW.Barrier()

    rit = args.working_rank
    for frame_val in range(len(frames)):
        if rank == rit:
            #print "creating frame", frames[frame_val], "on rank", rank
            if args.plot_time != None:
                if args.weight_field == 'None':
                    weight_field = None
                    pickle_file = save_dir + args.field + "_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
                else:
                    weight_field = args.weight_field
                    #pickle_file = path + args.field + "_movie_time_" + (str(args.plot_time)) + ".pkl"
                    pickle_file = save_dir + args.axis + '_' + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + ".pkl"
            else:
                pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
            if os.path.isfile(pickle_file) == False:
                print("NOT USING PICKLE ON RANK", rank)
                time_val = m_times[frame_val]
                print("FILE =", usable_files[frame_val])
                has_particles = has_sinks(usable_files[frame_val])
                if has_particles:
                    part_info = mym.get_particle_data(usable_files[frame_val], args.axis, proj_or=L)
                else:
                    part_info = {}
                if args.ax_lim != None:
                    if has_particles and args.image_center != 0:
                        xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center - 1], args.ax_lim + part_info['particle_position'][0][args.image_center - 1]]
                        ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center - 1], args.ax_lim + part_info['particle_position'][1][args.image_center - 1]]
                center_vel = [0.0, 0.0, 0.0]
                if args.image_center != 0 and has_particles:
                    original_positions = [X, Y, X_vel, Y_vel]
                    
                    x_pos = 0#np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                    y_pos = 0#np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                    X = X + x_pos
                    Y = Y + y_pos
                    X_vel = X_vel + x_pos
                    Y_vel = Y_vel + y_pos
                    
                    part_file = usable_files[frame_val][:-13] + 'part' + usable_files[frame_val][-6:]
                    f = h5py.File(part_file, 'r')
                    ordered_inds = np.argsort(f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['tag                     '])[0][0]])
                    center_vel_x = f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['velx                    '])[0][0]][ordered_inds][args.image_center - 1]
                    center_vel_y = f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['vely                    '])[0][0]][ordered_inds][args.image_center - 1]
                    center_vel_z = f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['velz                    '])[0][0]][ordered_inds][args.image_center - 1]
                    f.close()
                    center_vel = [center_vel_x, center_vel_y, center_vel_z]
                    print("CENTER_VEL=", center_vel)
                    f.close()
                xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
                #xlim = xlim - center_pos[0]
                #ylim = ylim - center_pos[1]
                if args.axis == 'xy':
                    center_vel=center_vel[:2]
                else:
                    center_vel=center_vel[::2]

                if args.yt_proj == False:
                    f = h5py.File(usable_files[frame_val], 'r')
                    image = get_image_arrays(f, simfo['field'], simfo, args, X, Y)
                    print("image shape=", np.shape(image))
                    print("grid shape=", np.shape(X))
                    magx = get_image_arrays(f, 'mag'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args, X, Y)
                    magy = get_image_arrays(f, 'mag'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args, X, Y)
                    x_pos_min = int(np.round(np.min(X) - simfo['xmin_full'])/simfo['cell_length'])
                    y_pos_min = int(np.round(np.min(Y) - simfo['xmin_full'])/simfo['cell_length'])
                    if np.shape(f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis]) == (2048, 2048):
                        velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis]]
                    elif args.axis == 'xy':
                        try:
                            velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0]]
                        except:
                            velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis]]
                    else:
                        velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:]]
                    velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X, velocity_data[0], velocity_data[1], center_vel=center_vel)
                    f.close()
                    #makine movie frame pickle
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
                    print("Opened pickle file")
                    #pickle.dump((usable_files[frame_val], X, Y, X_vel, Y_vel, image, velx, vely, part_info, args_dict, simfo, args, magx, magy), file)
                    pickle.dump((X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo), file)
                    print("Dumped data into pickle")
                    file.close()
            
            if args.pickle_dump == False:
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
                if frame_val > 0 or time_val > -1.0:
                    plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                else:
                    plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5)
                cbar = plt.colorbar(plot, pad=0.0)
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)

                if has_particles:
                    if args.annotate_particles_mass == True:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'],depth_array=part_info['depth_position'])
                    else:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None,depth_array=part_info['depth_position'])

                if args.annotate_time == "True":
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(time_val))+'yr', va="center", ha="left", color='w', fontsize=args.text_font)
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                    
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                if 'dens' in simfo['field']:
                    cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                else:
                    cbar.set_label(simfo['field'], rotation=270, labelpad=14, size=args.text_font)

                plt.tick_params(axis='both', which='major')# labelsize=16)
                for line in ax.xaxis.get_ticklines():
                    line.set_color('white')
                for line in ax.yaxis.get_ticklines():
                    line.set_color('white')

                plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight')
                #Convert to jpeg
                eps_image = Image.open(file_name + ".eps")
                eps_image.load(scale=4)
                eps_image.save(file_name + ".jpg")
            
                #os.system('convert -antialias -quality 100 -density 200 -resize 100% -flatten ' + file_name + '.eps ' + file_name + '.jpg')
                #os.remove(file_name + '.eps')
                print('Created frame', (frames[frame_val]+1), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.eps')
                #del image
                #del magx
                #del magy
                #del velx
                #del vely
                '''
                if args.image_center != 0 and has_particles:
                    X, Y, X_vel, Y_vel = original_positions
                '''
            
                sys.stdout.flush()
                CW.Barrier()

            else:
                print("Creating pickle")
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
                print("Built dictionary")
                try:
                    pickle_file = save_dir + str(args.axis) +'_' + args_dict['field'].split('_')[0] + "_thickness_" + str(int(args.slice_thickness)) +'_AU_movie_time_' + str(float(args_dict['time_val'])) + '.pkl'
                except:
                    pickle_file = save_dir + str(args.axis) +'_' + args_dict['field'][1] + "_thickness_" + str(int(args.slice_thickness)) +'_AU_movie_time_' + str(float(args_dict['time_val'])) + '.pkl'
                print("Got pickle file name:", pickle_file)
                file = open(pickle_file, 'wb')
                print("Opened pickle file")
                #pickle.dump((usable_files[frame_val], X, Y, X_vel, Y_vel, image, velx, vely, part_info, args_dict, simfo, args, magx, magy), file)
                pickle.dump((X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo), file)
                print("Dumped data into pickle")
                file.close()
                print("Created pickle:", pickle_file)
        
        rit = rit + 1
        if rit == size:
            rit = 0

    print("completed making movie frames on rank", rank)

if __name__ == '__main__': main()


