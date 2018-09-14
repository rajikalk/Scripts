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
from subprocess import call
import csv
import sys
import os
from matplotlib import transforms
import glob
import my_module as mym
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.patheffects as path_effects
import yt
import my_fields as myf

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-zt", "--zoom_times", help="0 is default zoom", default=0)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    #parser.add_argument("-st", "--start_time", help="What time would you like to start calculating times from?", type=float, default=0.0)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", type=int, default = 25)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-plr", "--plot_lref", help="would you like to annotate the refinement level?", default=False)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-sc", "--smooth_cells", help="how many cells would you like the smooth the velocities over? If not defined it is set to half the annotatation frequency")
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-mt", "--movie_times", help="What movies times would you like plotted?", type=list, default=[])
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-16)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("-ic", "--image_center", help="where would you like to center the image?", type=int, default=0)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=None)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=5000, type=int)
    parser.add_argument("-yt", "--yt_proj", help="Do you want to use yt to create projections as opposed to the movie files?", default=False)
    parser.add_argument("-proj_or", "--projection_orientation", help="Do you want to set the projection orientation? give as angle (in degrees) from positive y-axis", default=None, type=float)
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 100AU", type=float, default=300.)
    parser.add_argument("-use_disk", "--use_disk_angular_momentum", help="Do you want to use the disk angular momentum to define the normal vector for a projection?", default="False")
    parser.add_argument("-wf", "--weight_field", help="Do you want to have a weighted projection plot?", type=str, default='dens')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def get_files(path, args):
    if args.yt_proj == False:
        if args.axis == "xz":
            type = "proj"
        else:
            type = "slice"
        source_directory = sorted(glob.glob(path + 'WIND_' + type + '*'))
    else:
        source_directory = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
    return source_directory

def has_sinks(file):
    try:
        part_file = file[:-12] + 'part' + file[-5:]
        f = h5py.File(part_file, 'r')
        if f[f.keys()[1]][-1][-1] > 0:
            f.close()
            return True
        else:
            f.close()
            return False
    except:
        f = h5py.File(file, 'r')
        if "particlepositions" in f.keys():
            f.close()
            return True
        else:
            f.close()
            return False

def sim_info(path, file, args):
    """
    Finds particle info, relevant to frame size and such. NOTE ACCRETION RADIUS IS GIVEN FROM PARTICLE INFO FUNCTION
    """
    path_split = path.split('/')
    for p_s in path_split:
        if 'omega' in p_s:
            ang_val = p_s.split('_')[-1]
        else:
            ang_val = 0.2
        if 'lref' in p_s:
            temp_split = p_s.split('_')
            for t_s_it in range(len(temp_split)):
                if 'lref' == temp_split[t_s_it]:
                    lref = temp_split[t_s_it+1]
                if '0.' in temp_split[t_s_it]:
                    den_pert = temp_split[t_s_it]
                else:
                    den_pert = '0.0'
    lref = path.split('lref')[-1].split('/')[0].split('_')[-1]
    if '0.' in  path.split('lref')[0].split('/')[-1]:
        den_pert = path.split('lref')[0].split('/')[-1].split('_')[-2]
    else:
        den_pert = '0.0'
    movie_type = path.split('/')[-2].split('_')[0]

    try:
        f = h5py.File(file, 'r')
        for key in f.keys():
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
        if args.axis == "xz":
            type = "proj"
        else:
            type = "slice"
        annotate_freq = ((xmax/cl) - (xmin/cl))/31.
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
        type = "hdf5_plt_cnt"
        annotate_freq = dim/31.
    if args.smooth_cells == None:
        smoothing = annotate_freq/2
    else:
        smoothing = int(args.smooth_cells)
    sim_info = {'angular_momentum':ang_val,
                'movie_type':movie_type,
                'field': field,
                'dimension': dim,
                'zoom_cell': zoom_cell,
                'movie_file_type': type,
                'xmin': xmin,
                'xmax': xmax,
                'ymin': ymin,
                'ymax': ymax,
                'cell_length': cl,
                'annotate_freq': annotate_freq,
                'type': type,
                'smoothing': smoothing,
                'refinement_level': lref,
                'den_pert': den_pert,
                'xmin_full':  xmin_full
                }
    f.close()
    return sim_info

def get_image_arrays(f, field, simfo, args, X, Y):
    dim = int(simfo['dimension'])
    image = []
    xpos = int(np.round(np.min(X) - simfo['xmin_full'])/simfo['cell_length'])
    ypos = int(np.round(np.min(Y) - simfo['xmin_full'])/simfo['cell_length'])
    print "XPOS =", xpos
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
    other, with text ls[i] being shown in color lc[i].
    
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
    if os.path.exists(save_dir) == False:
        os.makedirs(save_dir)
    
    args = parse_inputs()
    mym.set_global_font_size(args.text_font)
    files = get_files(path, args)
    simfo = sim_info(path, files[0], args)
    if args.yt_proj == False:
        X, Y, X_vel, Y_vel, cl = mym.initialise_grid(files[-1], zoom_times=args.zoom_times)
        L=None
    else:
        x = np.linspace(simfo['xmin'], simfo['xmax'], simfo['dimension'])
        y = np.linspace(simfo['ymin'], simfo['ymax'], simfo['dimension'])
        X, Y = np.meshgrid(x, y)
                
        annotate_space = (simfo['xmax'] - simfo['xmin'])/31.
        x_ind = []
        y_ind = []
        counter = 0
        while counter < 31:
            val = annotate_space*counter + annotate_space/2. + simfo['xmin']
            x_ind.append(int(val))
            y_ind.append(int(val))
            counter = counter + 1
        X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
        if args.projection_orientation != None:
            y_val = 1./np.tan(np.deg2rad(args.projection_orientation))
            L = [1.0, y_val, 0.0]
        else:
            if args.axis == 'xz':
                L = [1.0, 0.0, 0.0]
            else:
                L = [0.0, 0.0, 1.0]
        myf.set_normal(L)
        print "SET PROJECTION ORIENTATION L=", myf.get_normal()
    if args.yt_proj == False and args.image_center != 0:
        sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
    elif args.yt_proj != False and args.image_center != 0:
        sim_files = files
    xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
    if args.ax_lim != None:
        xlim = [-1*args.ax_lim, args.ax_lim]
        ylim = [-1*args.ax_lim, args.ax_lim]
    myf.set_center(args.image_center)

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
    no_frames = len(m_times)
    m_times = m_times[args.start_frame:]
    sys.stdout.flush()
    CW.Barrier()

    usable_files = mym.find_files(m_times, files)
    del files
    if args.image_center != 0 and args.yt_proj == False:
        usable_sim_files = mym.find_files(m_times, sim_files)
        del sim_files
    sys.stdout.flush()
    CW.Barrier()
    frames = range(args.start_frame, no_frames)

    # Define colourbar bounds
    cbar_max = args.colourbar_max
    cbar_min = args.colourbar_min
    
    sys.stdout.flush()
    CW.Barrier()

    if args.yt_proj:
        yt.enable_parallelism()
        ts = yt.DatasetSeries(usable_files, parallel=size/16.)
        center_pos = np.array([0.0, 0.0, 0.0])
        thickness = yt.YTArray(args.slice_thickness, 'AU')
        if args.plot_time != None:
            if args.weight_field == 'None':
                weight_field = None
                pickle_file = path + args.field + "_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
            else:
                weight_field = args.weight_field
                pickle_file = path + args.field + "_movie_time_" + (str(args.plot_time)) + ".pkl"
        for ds in ts.piter():
            #for usable_file in usable_files:
            file_int = usable_files.index(path + str(ds))# (usable_file) #(path + str(ds))
            if args.plot_time is None:
                pickle_file = path + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl"
            if os.path.isfile(pickle_file) == False:
                print "PICKLE:", pickle_file,"DOESN'T EXIST. MAKING PROJECTION  FOR FRAME", frames[file_int], "ON RANK", rank
                time_val = m_times[file_int]
                has_particles = has_sinks(path + str(ds)) #(usable_file)#(path + str(ds))
                if has_particles:
                    part_info = mym.get_particle_data(path + str(ds), args.axis, proj_or=L)
                else:
                    part_info = {}
                if args.ax_lim != None:
                    if has_particles and args.image_center != 0:
                        xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center - 1], args.ax_lim + part_info['particle_position'][0][args.image_center - 1]]
                        ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center - 1], args.ax_lim + part_info['particle_position'][1][args.image_center - 1]]
                x_width = (xlim[1] -xlim[0])
                y_width = (ylim[1] -ylim[0])
                thickness = yt.YTArray(args.slice_thickness, 'AU')

                if args.image_center != 0 and has_particles:
                    original_positions = [X, Y, X_vel, Y_vel]
        
                    x_pos = np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                    y_pos = np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                    X = X + x_pos
                    Y = Y + y_pos
                    X_vel = X_vel + x_pos
                    Y_vel = Y_vel + y_pos

                print "in time series loop"
                #ds = yt.load(usable_file)
                dd = ds.all_data()
                print "loaded all data"
                center_pos = dd['Center_Position'].value
                if args.image_center != 0 and has_particles == True:
                    center_pos = dd['Center_Position'].value
                else:
                    center_pos = np.array([0.0, 0.0, 0.0])
                center_vel = dd['Center_Velocity'].value
                part_pos = dd['All_Particle_Positions']
                part_mass = dd['All_Particle_Masses']
                print "np.mean(dd['Particle_Potential']):", np.mean(dd['Particle_Potential'])
                print "center_vel =", myf.get_center_vel(), "on rank", rank, "for", ds
                if args.axis == 'xy':
                    center_vel=center_vel[:2]
                    if args.use_disk_angular_momentum != "False":
                        disk = ds.disk(center_pos, L, (args.ax_lim*2, 'au'), (args.slice_thickness*2, 'au'))
                        tot_vec = [np.sum(disk['Angular_Momentum_x']).value, np.sum(disk['Angular_Momentum_y']).value, np.sum(disk['Angular_Momentum_z']).value]
                        tot_mag = np.sqrt(tot_vec[0]**2. + tot_vec[1]**2. + tot_vec[2]**2.)
                        L = tot_vec/tot_mag
                        print "SET PROJECTION ORIENTATION L=", L
                else:
                    center_vel=center_vel[::2]
                if args.axis == "xz":
                    proj = yt.OffAxisProjectionPlot(ds, L, [simfo['field'], 'Projected_Velocity', 'velz', 'Projected_Magnetic_Field', 'magz'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'), weight_field=weight_field)
                    if weight_field == None:
                        image = proj.frb.data[simfo['field']]
                        velx_full = proj.frb.data[('gas', 'Projected_Velocity')].in_units('cm**2/s').value
                        vely_full = proj.frb.data[('flash', 'velz')].in_units('cm**2/s').value
                        magx = proj.frb.data[('gas', 'Projected_Magnetic_Field')].in_units('cm*gauss')
                        magy = proj.frb.data[('flash', 'magz')].in_units('cm*gauss')
                    else:
                        image = proj.frb.data[simfo['field']].value
                        velx_full = proj.frb.data[('gas', 'Projected_Velocity')].in_units('cm/s').value
                        vely_full = proj.frb.data[('flash', 'velz')].in_units('cm/s').value
                        magx = proj.frb.data[('gas', 'Projected_Magnetic_Field')].in_units('gauss').value
                        magy = proj.frb.data[('flash', 'magz')].in_units('gauss').value
                else:
                    proj = yt.OffAxisProjectionPlot(ds, L, [simfo['field'], 'velx', 'vely', 'magx', 'magy'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'), weight_field=weight_field)
                    if weight_field == None:
                        image = proj.frb.data[simfo['field']]
                        velx_full = proj.frb.data[('flash', 'velx')].in_units('cm**2/s').value
                        vely_full = proj.frb.data[('flash', 'vely')].in_units('cm**2/s').value
                        magx = proj.frb.data[('flash', 'magx')].in_units('cm*gauss')
                        magy = proj.frb.data[('flash', 'magy')].in_units('cm*gauss')
                    else:
                        image = proj.frb.data[simfo['field']].value
                        velx_full = proj.frb.data[('flash', 'velx')].in_units('cm/s').value
                        vely_full = proj.frb.data[('flash', 'vely')].in_units('cm/s').value
                        magx = proj.frb.data[('flash', 'magx')].in_units('gauss').value
                        magy = proj.frb.data[('flash', 'magy')].in_units('gauss').value

                velx, vely = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=center_vel)
                del velx_full
                del vely_full

                print "Creating pickle"
                args_dict = {}
                if args.annotate_time == "True":
                    args_dict.update({'annotate_time': '$t$='+str(int(time_val))+'yr'})
                if args.plot_lref == True:
                    args_dict.update({'annotate_lref': '$r_{acc}$='+str(r_acc)+'AU'})
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

                if rank == 0:
                    file = open(pickle_file, 'w+')
                    pickle.dump((X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo), file)
                    file.close()
                print "Created Pickle:", pickle_file, "for  file:", str(ds)
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
                
    sys.stdout.flush()
    CW.Barrier()


    rit = args.working_rank
    for frame_val in range(len(frames)):
        if rank == rit:
            print "creating frame", frames[frame_val], "on rank", rank
            if args.plot_time != None:
                if args.weight_field == 'None':
                    weight_field = None
                    pickle_file = path + args.field + "_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
                else:
                    weight_field = args.weight_field
                    pickle_file = path + args.field + "_movie_time_" + (str(args.plot_time)) + ".pkl"
            else:
                pickle_file = path + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
            if args.yt_proj and os.path.isfile(pickle_file):
                print "USING PICKLED FILE:", pickle_file
                file = open(pickle_file, 'r')
                X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo = pickle.load(file)
                #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                file.close()
                xlim = args_dict['xlim']
                ylim = args_dict['ylim']
                has_particles = args_dict['has_particles']
                time_val = args_dict['time_val']
                xabel = args_dict['xabel']
                yabel = args_dict['yabel']
            else:
                print "NOT USING PICKLE ON RANK", rank
                time_val = m_times[frame_val]
                print "FILE =", usable_files[frame_val]
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
                    
                    x_pos = np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                    y_pos = np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                    X = X + x_pos
                    Y = Y + y_pos
                    X_vel = X_vel + x_pos
                    Y_vel = Y_vel + y_pos
                    
                    part_file = usable_files[frame_val][:-12] + 'part' + usable_files[frame_val][-5:]
                    f = h5py.File(part_file, 'r')
                    ordered_inds = np.argsort(f[f.keys()[11]][:,np.where(f[f.keys()[5]][:] == ['tag                     '])[0][0]])
                    center_vel_x = f[f.keys()[11]][:,np.where(f[f.keys()[5]][:] == ['velx                    '])[0][0]][ordered_inds][args.image_center - 1]
                    center_vel_y = f[f.keys()[11]][:,np.where(f[f.keys()[5]][:] == ['vely                    '])[0][0]][ordered_inds][args.image_center - 1]
                    center_vel_z = f[f.keys()[11]][:,np.where(f[f.keys()[5]][:] == ['velz                    '])[0][0]][ordered_inds][args.image_center - 1]
                    f.close()
                    center_vel = [center_vel_x, center_vel_y, center_vel_z]
                    print "CENTER_VEL=", center_vel
                    f.close()
                xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
                if args.axis == 'xy':
                    center_vel=center_vel[:2]
                else:
                    center_vel=center_vel[::2]

                if args.yt_proj == False:
                    f = h5py.File(usable_files[frame_val], 'r')
                    image = get_image_arrays(f, simfo['field'], simfo, args, X, Y)
                    print "image shape=", np.shape(image)
                    print "grid shape=", np.shape(X)
                    magx = get_image_arrays(f, 'mag'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args, X, Y)
                    magy = get_image_arrays(f, 'mag'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args, X, Y)
                    x_pos_min = int(np.round(np.min(X) - simfo['xmin_full'])/simfo['cell_length'])
                    y_pos_min = int(np.round(np.min(Y) - simfo['xmin_full'])/simfo['cell_length'])
                    if np.shape(f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis]) == (2048, 2048):
                        velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis]]
                    elif args.axis == 'xy':
                        velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0]]
                    else:
                        velocity_data = [f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,0,:], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,0,:]]
                    velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X, velocity_data[0], velocity_data[1], center_vel=center_vel)
                    f.close()
            
            if args.pickle_dump == False:
                plt.clf()
                fig, ax = plt.subplots()
                ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
                ax.set_ylabel(yabel, labelpad=-20, fontsize=args.text_font)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                if 0.0 in (cbar_min, cbar_max):
                    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.brg, rasterized=True, vmin=cbar_min, vmax=cbar_max)
                else:
                    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
                plt.gca().set_aspect('equal')
                if frame_val > 0 or time_val > -1.0:
                    plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                else:
                    plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5)
                
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)

                if has_particles:
                    if args.annotate_particles_mass == True:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'],depth_array=part_info['depth_position'])
                    else:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None,depth_array=part_info['depth_position'])

                if args.plot_lref == True:
                    r_acc = np.round(part_info['accretion_rad'])
                    ax.annotate('$r_{acc}$='+str(r_acc)+'AU', xy=(0.98*simfo['xmax'], 0.93*simfo['ymax']), va="center", ha="right", color='w', fontsize=args.text_font)

                if args.annotate_time == "True":
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(time_val))+'yr', va="center", ha="left", color='w', fontsize=args.text_font)
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                    
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                cbar = plt.colorbar(plot, pad=0.0)
                if simfo['field'] == 'dens':
                    cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                else:
                    cbar.set_label(simfo['field'], rotation=270, labelpad=14, size=args.text_font)

                plt.tick_params(axis='both', which='major')# labelsize=16)
                for line in ax.xaxis.get_ticklines():
                    line.set_color('white')
                for line in ax.yaxis.get_ticklines():
                    line.set_color('white')
        
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

                plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight')
            
                call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
                os.remove(file_name + '.eps')
                print 'Created frame', (frames[frame_val]+1), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.eps'
                del image
                del magx
                del magy
                del velx
                del vely
                
                if args.image_center != 0 and has_particles:
                    X, Y, X_vel, Y_vel = original_positions
            
                sys.stdout.flush()
                CW.Barrier()

            else:
                print "Creating pickle"
                args_dict = {}
                if args.annotate_time == "True":
                    args_dict.update({'annotate_time': '$t$='+str(int(time_val))+'yr'})
                if args.plot_lref == True:
                    args_dict.update({'annotate_lref': '$r_{acc}$='+str(r_acc)+'AU'})
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
                print "Built dictionary"
                pickle_file = save_dir + 'movie_pickle.pkl'
                print "Got pickle file name"
                file = open(pickle_file, 'w+')
                print "Opened pickle file"
                #pickle.dump((usable_files[frame_val], X, Y, X_vel, Y_vel, image, velx, vely, part_info, args_dict, simfo, args, magx, magy), file)
                pickle.dump((X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo), file)
                print "Dumped data into pickle"
                file.close()
                print "Created pickle"
    
        rit = rit +1
        if rit == size:
            rit = 0

    print "completed making movie frames on rank", rank

if __name__ == '__main__': main()


