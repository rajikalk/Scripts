#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import h5py
import numpy as np
import numpy.ma as ma
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
from subprocess import call
import csv
import sys
import os
from matplotlib import transforms
import glob
import my_module as mym
import pickle
import matplotlib.patheffects as path_effects
import yt
import my_fields as myf
import ast
import matplotlib.gridspec as gridspec
import argparse

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

def define_constants():
    constants = {'year':31557600.0, 'au':1.496e13, 'Msun':1.98841586e+33}
    return constants

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-zt", "--zoom_times", help="0 is default zoom", default=0)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10, type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    parser.add_argument("-st", "--start_time", help="What time woudl you like to start calculating times from?", type=int)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", default = 25)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-plr", "--plot_lref", help="would you like to annotate the refinement level?", default=False)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-sc", "--smooth_cells", help="how many cells would you like the smooth the velocities over? If not defined it is set to half the annotatation frequency")
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
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
    parser.add_argument("-in_file", "--input_file", help="file in input data about the mosaic plot")
    parser.add_argument("-save_dir", "--save_directory", help="where do you want to save the frames?", default='./')
    parser.add_argument("-sx", "--share_x", help="do you want to share the x axis?", default=False)
    parser.add_argument("-sy", "--share_y", help="do you want to share the y axis?", default=True)
    parser.add_argument("-sa", "--share_ax", help="do you want to share axes?", default=True)
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 100AU", type=float, default=300.)
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
    c = define_constants()
    type = get_files(path, args)[1]
    path_split = path.split('/')
    for p_s in path_split:
        if 'omega' in p_s:
            ang_val = p_s.split('_')[-1]
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
        if 'r_accretion' in list(f.keys()):
            racc = f['r_accretion'][0]/yt.units.AU.in_units('cm').value
        else:
            racc = 0.0
        for key in list(f.keys()):
            if args.field in key:
                field = key
        dim = np.shape(f[field])[0]
        if args.zoom_times > 0:
            zoom_cell = int((dim - dim/float(args.zoom_times))/2.)
        else:
            zoom_cell = 0
        xmin = f['minmax_xyz'][0][0]/c['au']
        xmax = f['minmax_xyz'][0][1]/c['au']
        xmin_full = xmin
        cl = (xmax-xmin)/dim
        cell_positions = np.arange(xmin, xmax-1, cl)
        xmin = f['minmax_xyz'][0][0]/c['au'] + zoom_cell*cl
        xmax = f['minmax_xyz'][0][1]/c['au'] - zoom_cell*cl
        if args.axis == "xy":
            ymin = f['minmax_xyz'][1][0]/c['au'] + zoom_cell*cl
            ymax = f['minmax_xyz'][1][1]/c['au'] - zoom_cell*cl
        else:
            ymin = f['minmax_xyz'][2][0]/c['au'] + zoom_cell*cl
            ymax = f['minmax_xyz'][2][1]/c['au'] - zoom_cell*cl
        if args.axis == "xz":
            type = "proj"
        else:
            type = "slice"
        annotate_freq = ((xmax/cl) - (xmin/cl))/31.
    except:
        f = h5py.File(file, 'r')
        if has_sinks(file):
            racc = f[list(f.keys())[18]][109][-1]/yt.units.AU.in_units('cm').value
        else:
            racc = 0.0
        f.close()
        if args.field == 'dens':
            field = ('flash', 'dens')
        else:
            part_file = file[:-12] + 'part' + file[-5:]
            f = yt.load(file, particle_filename=part_file)
            import pdb
            pdb.set_trace()
            field = f.derived_field_list[[x[1] for x in f.derived_field_list].index(args.field)]
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
    print("PLOT FIELD IS", field)
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
                'r_accretion': racc,
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
    #print "XPOS =", xpos
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
    args = parse_inputs()
    prev_args = args
    print("Starting mosaic_mod_script on rank", rank)

    # Read in directories:
    input_file = args.input_file
    #save_dir = args.save_directory
    #if os.path.exists(save_dir) == False:
    #    os.makedirs(save_dir)

    # Read in input file
    print("Reading in input mosaic file on rank", rank)
    positions = []
    paths = []
    args_dict = []
    with open(input_file, 'rU') as mosaic_file:
        reader = csv.reader(mosaic_file)
        for row in reader:
            if row[0] == 'Grid_inputs:':
                glr = float(row[1])
                grl = float(row[2])
                glw = float(row[3])
                ghspace = float(row[4])
            elif row[0][0] != '#':
                positions.append((int(row[0]), int(row[1])))
                paths.append(row[2])
                dict = ""
                for col in row[3:]:
                    dict = dict + col
                    if col != row[-1]:
                        dict = dict + ','
                dict = ast.literal_eval(dict)
                args_temp = argparse.Namespace(**vars(args))
                for key in list(dict.keys()):
                    if key in args:
                        exec("args_temp."+ key + " = " + "str(dict[key])")
                args_dict.append(args_temp)
                del args_temp
                args = prev_args
                
    import pdb
    pdb.set_trace()

    positions = np.array(positions)

    c = define_constants()
    mym.set_global_font_size(args.text_font)
    files = []
    simfo = []
    X = []
    Y = []
    X_vel = []
    Y_vel = []
    sim_files = []
    L = None
    for pit in range(len(paths)):
        fs = get_files(paths[pit], args_dict[pit])
        files.append(fs)

        #print "paths =", paths
        #print "fs =", fs
        #print "args_dict =", args_dict
        sfo = sim_info(paths[pit], fs[-1], args_dict[pit])
        simfo.append(sfo)

        if args_dict[pit].yt_proj == False:
            x, y, x_vel, y_vel, cl = mym.initialise_grid(files[pit][-1], zoom_times=args_dict[pit].zoom_times)
            X.append(x)
            Y.append(y)
            X_vel.append(x_vel)
            Y_vel.append(y_vel)
        else:
            x = np.linspace(sfo['xmin'], sfo['xmax'], sfo['dimension'])
            y = np.linspace(sfo['ymin'], sfo['ymax'], sfo['dimension'])
            x, y  = np.meshgrid(x, y)
            
            annotate_space = (simfo[pit]['xmax'] - simfo[pit]['xmin'])/31.
            x_ind = []
            y_ind = []
            counter = 0
            while counter < 31:
                val = annotate_space*counter + annotate_space/2. + simfo[pit]['xmin']
                x_ind.append(int(val))
                y_ind.append(int(val))
                counter = counter + 1
            x_vel, y_vel = np.meshgrid(x_ind, y_ind)
            if args_dict[pit].projection_orientation != None:
                y_val = 1./np.tan(np.deg2rad(float(args_dict[pit].projection_orientation)))
                if np.isinf(y_val):
                    y_val = 0.0
                L = [1.0, y_val, 0.0]
            else:
                if has_particles == False or len(dd['particle_posx']) == 1:
                    L = [0.0, 1.0, 0.0]
                else:
                    pos_vec = [np.diff(dd['particle_posx'].value)[0], np.diff(dd['particle_posy'].value)[0]]
                    L = [-1*pos_vec[-1], pos_vec[0]]
                    L.append(0.0)
                    if L[0] > 0.0:
                        L = [-1.0*L[0], -1.0*L[1], 0.0]
            print("SET PROJECTION ORIENTATION L=", L)
            L = np.array(L)
            X.append(x)
            Y.append(y)
            X_vel.append(x_vel)
            Y_vel.append(y_vel)
        if rank == 0:
            print("shape of x, y", np.shape(x), np.shape(y))

        if args_dict[pit].yt_proj == False and args_dict[pit].image_center != 0:
            sim_fs = sorted(glob.glob(paths[pit] + 'WIND_hdf5_plt_cnt*'))
        elif args_dict[pit].yt_proj != False and args_dict[pit].image_center != 0:
            sim_fs = files
        else:
            sim_fs = []
        sim_files.append(sim_fs)
    #myf.set_normal(L)
    #print "SET PROJECTION ORIENTATION L=", myf.get_normal()

    # Initialise Grid and build lists
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files[0], args.time_step, presink_frames=args.presink_frames, end_time=args.end_time)
    no_frames = len(m_times)
    m_times = m_times[args.start_frame:]
    sys.stdout.flush()
    CW.Barrier()

    usable_files = []
    usable_sim_files = []
    for pit in range(len(paths)):
        usable_fs = mym.find_files(m_times, files[pit])
        usable_files.append(usable_fs)
        if args_dict[pit].image_center != 0 and args_dict[pit].yt_proj == False:
            usable_sfs = mym.find_files(m_times, sim_files[pit])
            usable_sim_files.append(usable_fs)
            del sim_files[pit]
        else:
            usable_sim_files.append([])
    sys.stdout.flush()
    CW.Barrier()
    frames = list(range(args.start_frame, no_frames))

    sink_form_time = []
    for pit in range(len(paths)):
        sink_form = mym.find_sink_formation_time(files[pit])
        print("sink_form_time", sink_form_time)
        sink_form_time.append(sink_form)
    del files

    # Define colourbar bounds
    cbar_max = args.colourbar_max
    cbar_min = args.colourbar_min

    if L is None:
        if args.axis == 'xy':
            L = [0.0, 0.0, 1.0]
        else:
            L = [1.0, 0.0, 0.0]
        L = np.array(L)
    if args.axis == 'xy':
        y_int = 1
    else:
        y_int = 2

    sys.stdout.flush()
    CW.Barrier()
    rit = args.working_rank
    for frame_val in range(len(frames)):
        if rank == rit:
            time_val = m_times[frame_val]
            plt.clf()
            columns = np.max(positions[:,0])
            rows = np.max(positions[:,1])

            width = float(columns)*(14.5/3.)
            height = float(rows)*(17./4.)
            fig =plt.figure(figsize=(width, height))
            
            gs_left = gridspec.GridSpec(rows, columns-1)
            gs_right = gridspec.GridSpec(rows, 1)

            gs_left.update(right=glr, wspace=glw, hspace=ghspace)
            gs_right.update(left=grl, hspace=ghspace)
            
            axes_dict = {}
            counter = 1

            for pit in range(len(paths)):
                
                try:
                    title_parts = args_dict[pit].title
                except:
                    title_parts = args_dict[pit]['title']
                title = ''
                for part in title_parts:
                    if part != title_parts[-1]:
                        title = title + part + ' '
                    else:
                        title = title + part
            
                ax_label = 'ax' + str(counter)
                yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                if positions[pit][0] == 1 and positions[pit][1] == 1:
                    if columns > 1:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[0,0])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    else:
                        axes_dict.update({ax_label:fig.add_subplot(gs_right[0,0])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                elif positions[pit][0] != columns:
                    if args.share_x and args.share_y:
                        if yit >= len(axes_dict):
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'])})
                            #print "ADDED SUBPLOT:", counter, "on rank", rank
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                            #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_x:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[it][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_y and positions[pit][0]!=1:
                        yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_y:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    else:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                else:
                    if args.share_x and args.share_y:
                        yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                        axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0], sharex=axes_dict['ax1'], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_x:
                        axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0], sharex=axes_dict['ax1'])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_y:
                        yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                        axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    else:
                        axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank

                counter = counter + 1
                axes_dict[ax_label].set(adjustable='box-forced', aspect='equal')
                

                if args.yt_proj and args.plot_time==None and os.path.isfile(paths[pit] + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"):
                    pickle_file = paths[pit] + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
                    print("USING PICKLED FILE:", pickle_file)
                    file = open(pickle_file, 'r')
                    #weight_fieldstuff = pickle.load(file)
                    X[pit], Y[pit], image, magx, magy, X_vel[pit], Y_vel[pit], velx, vely, part_info, args_dict[pit], simfo[pit] = pickle.load(file)

                    #file_time = stuff[17]
                    file.close()

                else:
                    time_val = m_times[frame_val]
                    print("FILE =", usable_files[pit][frame_val])
                    has_particles = has_sinks(usable_files[pit][frame_val])
                    if has_particles:
                        part_info = mym.get_particle_data(usable_files[pit][frame_val], args_dict[pit].axis, proj_or=L)
                    else:
                        part_info = {}
                    center_vel = [0.0, 0.0, 0.0]
                    if args.image_center != 0 and has_particles:
                        original_positions = [X[pit], Y[pit], X_vel[pit], y_vel[pit]]
                        x_pos = np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                        y_pos = np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                        pos = np.array([part_info['particle_position'][0][args.image_center - 1], part_info['particle_position'][1][args.image_center - 1]])
                        X[pit] = X[pit] + x_pos
                        Y[pit] = Y[pit] + y_pos
                        X_vel[pit] = X_vel[pit] + x_pos
                        Y_vel[pit] = Y_vel[pit] + y_pos
                        if args.yt_proj == False:
                            sim_file = usable_sim_files[frame_val][:-12] + 'part' + usable_sim_files[frame_val][-5:]
                        else:
                            sim_file = part_file
                        if len(part_info['particle_mass']) == 1:
                            part_ind = 0
                        else:
                            min_dist = 1000.0
                            for part in range(len(part_info['particle_mass'])):
                                f = h5py.File(sim_file, 'r')
                                temp_pos = np.array([f[list(f.keys())[11]][part][13]/c['au'], f[list(f.keys())[11]][part][13+y_int]/c['au']])
                                f.close()
                                dist = np.sqrt(np.abs(np.diff((temp_pos - pos)**2)))[0]
                                if dist < min_dist:
                                    min_dist = dist
                                    part_ind = part
                        f = h5py.File(sim_file, 'r')
                        center_vel = [f[list(f.keys())[11]][part_ind][18], f[list(f.keys())[11]][part_ind][19], f[list(f.keys())[11]][part_ind][20]]
                        f.close()
                    xabel, yabel, xlim, ylim = image_properties(X[pit], Y[pit], args_dict[pit], simfo[pit])
                    if args_dict[pit].axis == 'xy':
                        center_vel=center_vel[:2]
                    else:
                        center_vel=center_vel[::2]
                    
                    if args_dict[pit].ax_lim != None:
                        if has_particles and args_dict[pit].image_center != 0:
                            xlim = [-1*args_dict[pit].ax_lim + part_info['particle_position'][0][args_dict[pit].image_center - 1], args_dict[pit].ax_lim + part_info['particle_position'][0][args_dict[pit].image_center - 1]]
                            ylim = [-1*args_dict[pit].ax_lim + part_info['particle_position'][1][args_dict[pit].image_center - 1], args_dict[pit].ax_lim + part_info['particle_position'][1][args_dict[pit].image_center - 1]]
                        else:
                            xlim = [-1*args_dict[pit].ax_lim, args_dict[pit].ax_lim]
                            ylim = [-1*args_dict[pit].ax_lim, args_dict[pit].ax_lim]

                    if args.yt_proj == False:
                        f = h5py.File(usable_files[pit][frame_val], 'r')
                        image = get_image_arrays(f, simfo[pit]['field'], simfo[pit], args_dict[pit], X[pit], Y[pit])
                        magx = get_image_arrays(f, 'mag'+args.axis[0]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis, simfo[pit], args_dict[pit], X[pit], Y[pit])
                        magy = get_image_arrays(f, 'mag'+args.axis[1]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis, simfo[pit], args_dict[pit], X[pit], Y[pit])
                        x_pos_min = int(np.round(np.min(X[pit]) - simfo[pit]['xmin_full'])/simfo[pit]['cell_length'])
                        y_pos_min = int(np.round(np.min(Y[pit]) - simfo[pit]['xmin_full'])/simfo[pit]['cell_length'])
                        if np.shape(f['vel'+args.axis[0]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis]) == (2048, 2048):
                            velocity_data = [f['vel'+args.axis[0]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis], f['vel'+args.axis[1]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis]]
                        elif args.axis == 'xy':
                            velocity_data = [f['vel'+args.axis[0]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis][:,:,0], f['vel'+args.axis[1]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis][:,:,0]]
                        else:
                            velocity_data = [f['vel'+args.axis[0]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis][:,0,:], f['vel'+args.axis[1]+'_'+simfo[pit]['movie_file_type']+'_'+args.axis][:,0,:]]
                        velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X[pit], velocity_data[0], velocity_data[1], center_vel=center_vel)
                    else:
                        if args_dict[pit].image_center == 0 or has_particles == False:
                            center_pos = np.array([0.0, 0.0, 0.0])
                        else:
                            dd = f.all_data()
                            center_pos = np.array([dd['particle_posx'][args.image_center-1].in_units('AU'), dd['particle_posy'][args.image_center-1].in_units('AU'), dd['particle_posz'][args.image_center-1].in_units('AU')])
                        x_width = (xlim[1] -xlim[0])
                        y_width = (ylim[1] -ylim[0])
                        thickness = yt.YTArray(args.slice_thickness, 'AU')
                        
                        proj = yt.OffAxisProjectionPlot(f, L, [simfo[pit]['field'], 'cell_mass', 'velz_mw', 'magz_mw', 'Projected_Magnetic_Field_mw', 'Projected_Velocity_mw'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))
                        image = (proj.frb.data[simfo[pit]['field']]/thickness.in_units('cm')).value
                        velx_full = (proj.frb.data[('gas', 'Projected_Velocity_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).value
                        vely_full = (proj.frb.data[('gas', 'velz_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).value
                        magx = (proj.frb.data[('gas', 'Projected_Magnetic_Field_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).value
                        magy = (proj.frb.data[('gas', 'magz_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).value
                        mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm')).value
                        
                        velx_full = velx_full/mass
                        vely_full = vely_full/mass
                        magx = magx/mass
                        magy = magy/mass
                        del mass

                        velx, vely = mym.get_quiver_arrays(0.0, 0.0, X[pit], velx_full, vely_full, center_vel=center_vel)
                        del velx_full
                        del vely_full

                        if len(frames) == 1:
                            if rank == 0:
                                pickle_file = paths[pit] + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
                                file = open(pickle_file, 'w+')
                                pickle.dump((X[pit], Y[pit], image, magx, magy, X_vel[pit], Y_vel[pit], velx, vely, xlim, ylim, has_particles, part_info, simfo[pit], time_val,xabel, yabel), file)
                                file.close()
                                print("Created Pickle:", pickle_file, "for  file:", usable_files[pit][frame_val])
                        else:
                            pickle_file = paths[pit] + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
                            file = open(pickle_file, 'w+')
                            pickle.dump((X[pit], Y[pit], image, magx, magy, X_vel[pit], Y_vel[pit], velx, vely, xlim, ylim, has_particles, part_info, simfo[pit], time_val,xabel, yabel), file)
                            file.close()
                            print("Created Pickle:", pickle_file, "for  file:", usable_files[pit][frame_val])
                    
                    f.close()

                plot = axes_dict[ax_label].pcolormesh(X[pit], Y[pit], image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
                plt.gca().set_aspect('equal')
                if frame_val > 0 or time_val > -1.0:
                    axes_dict[ax_label].streamplot(X[pit], Y[pit], magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                else:
                    axes_dict[ax_label].streamplot(X[pit], Y[pit], magx, magy, density=4, linewidth=0.25, minlength=0.5)

                xlim = args_dict[pit]['xlim']
                ylim = args_dict[pit]['ylim']
                mym.my_own_quiver_function(axes_dict[ax_label], X_vel[pit], Y_vel[pit], velx, vely, plot_velocity_legend=bool(args_dict[pit]['annotate_velocity']), limits=[xlim, ylim], standard_vel=args.standard_vel)
                if args_dict[pit]['has_particles']:
                    if args.annotate_particles_mass == True:
                        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'])
                    else:
                        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)
                if args.plot_lref == True:
                    r_acc = np.round(part_info['accretion_rad'])
                    axes_dict[ax_label].annotate('$r_{acc}$='+str(r_acc)+'AU', xy=(0.98*simfo[pit]['xmax'], 0.93*simfo[pit]['ymax']), va="center", ha="right", color='w', fontsize=args_dict[pit].text_font)
                if args.annotate_time == "True" and pit == 0:
                    print("ANNONTATING TIME:", str(int(time_val))+'yr')
                    time_text = axes_dict[ax_label].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(time_val))+'yr', va="center", ha="left", color='w', fontsize=args.text_font)
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                    #ax.annotate('$t$='+str(int(time_val))+'yr', xy=(xlim[0]+0.01*(xlim[1]-xlim[0]), ylim[1]-0.03*(ylim[1]-ylim[0])), va="center", ha="left", color='w', fontsize=args.text_font)
                title_text = axes_dict[ax_label].text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+2))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                if positions[pit][0] == columns:
                    cbar = plt.colorbar(plot, pad=0.0, ax=axes_dict[ax_label])
                    cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                axes_dict[ax_label].set_xlabel(args_dict[pit]['xabel'], labelpad=-1, fontsize=args.text_font)
                if positions[pit][0] == 1:
                    axes_dict[ax_label].set_ylabel(args_dict[pit]['yabel'], labelpad=-20, fontsize=args.text_font)
                axes_dict[ax_label].set_xlim(xlim)
                axes_dict[ax_label].set_ylim(ylim)
                for line in axes_dict[ax_label].xaxis.get_ticklines():
                    line.set_color('white')
                for line in axes_dict[ax_label].yaxis.get_ticklines():
                    line.set_color('white')

                plt.tick_params(axis='both', which='major', labelsize=16)
                for line in axes_dict[ax_label].xaxis.get_ticklines():
                    line.set_color('white')
                for line in axes_dict[ax_label].yaxis.get_ticklines():
                    line.set_color('white')

                if positions[pit][0] != 1:
                    yticklabels = axes_dict[ax_label].get_yticklabels()
                    plt.setp(yticklabels, visible=False)

                if positions[pit][0] == 1:
                    axes_dict[ax_label].tick_params(axis='y', which='major', labelsize=args.text_font)
                if positions[pit][1] == rows:
                    axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
                    if positions[pit][0] != 1:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels[0], visible=False)

                if len(usable_files[pit]) > 1:
                    if args.output_filename == None:
                        import pdb
                        pdb.set_trace()
                        file_name = save_dir + "movie_frame_" + ("%06d" % frames[frame_val])
                    else:
                        file_name = args.output_filename + "_" + str(int(time_val))
                else:
                    if args.output_filename != None:
                        file_name = args.output_filename
                    else:
                        import pdb
                        pdb.set_trace()
                        file_name = save_dir + "time_" + str(args.plot_time)

                plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight')
                #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                
                #plt.savefig(file_name + ".jpg", format='jpeg', bbox_inches='tight')
                call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
                os.remove(file_name + '.eps')

                del image
                del magx
                del magy
                del velx
                del vely
                
                if args.image_center != 0 and has_particles:
                    X[pit], Y[pit], X_vel[pit], Y_vel[pit] = original_positions
            print('Created frame', (frames[frame_val]), 'of', str(frames[-1]), 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.eps')

        rit = rit +1
        if rit == size:
            rit = 0

    print("completed making movie frames on rank", rank)

if __name__ == '__main__': main()

#Read in inputs:


