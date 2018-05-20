#!/usr/bin/env python
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
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.patheffects as path_effects
import yt

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

def define_constants():
    constants = {'year':31557600.0, 'au':1.496e13, 'Msun':1.98841586e+33}
    return constants

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-zt", "--zoom_times", help="0 is default zoom", default=0)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    #parser.add_argument("-st", "--start_time", help="What time woudl you like to start calculating times from?", type=float, default=0.0)
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
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 100AU", type=float, default=100.)
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

def has_sinks(f):
    try:
        if "particlepositions" in f.keys():
            return True
        else:
            return False
    except:
        if ('io', u'particle_mass') in f.field_list:
            return True
        else:
            return False

def sim_info(path, file, args):
    c = define_constants()
    type = get_files(path, args)[1]
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
        if 'r_accretion' in f.keys():
            racc = f['r_accretion'][0]/yt.units.AU.in_units('cm').value
        else:
            racc = 0.0
        for key in f.keys():
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
        part_file = file[:-12] + 'part' + file[-5:]
        f = yt.load(file, particle_filename=part_file)
        dd = f.all_data()
        if has_sinks(f):
            racc = np.min(dd['dx'].in_units('au').value)*2.5
        else:
            racc = 0.0
        for key in f.field_list:
            if args.field in key:
                field = key
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
    
    c = define_constants()
    args = parse_inputs()
    mym.set_global_font_size(args.text_font)
    files = get_files(path, args)
    simfo = sim_info(path, files[-1], args)
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
            L = [1, y_val, 0]
        else:
            if args.axis == 'xz':
                L = [1.0, 0.0, 0.0]
            else:
                L = [0.0, 0.0, 1.0]
            '''
            else:
                pos_vec = [np.diff(dd['particle_posx'].value)[0], np.diff(dd['particle_posy'].value)[0]]
                L = [-1*pos_vec[-1], pos_vec[0]]
                L.append(0.0)
                if L[0] > 0:
                    L = [-1*L[0], -1*L[1], 0.0]
            '''
        print "SET PROJECTION ORIENTATION L=", L
    if args.yt_proj == False and args.image_center != 0:
        sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
    elif args.yt_proj != False and args.image_center != 0:
        sim_files = files
    
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
    if args.image_center != 0 and args.yt_proj == False:
        usable_sim_files = mym.find_files(m_times, sim_files)
        del sim_files
    sys.stdout.flush()
    CW.Barrier()
    frames = range(args.start_frame, no_frames)

    sink_form_time = mym.find_sink_formation_time(files)
    print "sink_form_time", sink_form_time, "on rank", rank
    del files

    # Define colourbar bounds
    cbar_max = args.colourbar_max
    cbar_min = args.colourbar_min

    if args.axis == 'xy':
        y_int = 1
    else:
        y_int = 2
    
    sys.stdout.flush()
    CW.Barrier()
    rit = args.working_rank
    for frame_val in range(len(frames)):
        if rank == rit:
            print "creating frame", frames[frame_val], "on rank", rank
            if args.yt_proj and args.plot_time==None and os.path.isfile(path + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"):
                pickle_file = path + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
                print "USING PICKLED FILE:", pickle_file
                file = open(pickle_file, 'r')
                X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                file.close()
            else:
                print "NOT USING PICKLE ON RANK", rank
                time_val = m_times[frame_val]
                if args.yt_proj != False:
                    file = usable_files[frame_val]
                    part_file = file[:-12] + 'part' + file[-5:]
                    f = yt.load(file, particle_filename=part_file)
                    dd = f.all_data()
                    #file_time = f.current_time.in_units('yr').value - sink_form_time
                else:
                    f = h5py.File(usable_files[frame_val], 'r')
                    #file_time = (f['time'][0]/yt.units.yr.in_units('s').value)-sink_form_time
                print "FILE =", usable_files[frame_val]
                has_particles = has_sinks(f)
                if has_particles:
                    part_info = mym.get_particle_data(usable_files[frame_val], args.axis, proj_or=L)
                else:
                    part_info = {}
                center_vel = [0.0, 0.0, 0.0]
                if args.image_center != 0 and has_particles:
                    
                    original_positions = [X, Y, X_vel, Y_vel]
                    x_pos = np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                    y_pos = np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                    pos = np.array([part_info['particle_position'][0][args.image_center - 1], part_info['particle_position'][1][args.image_center - 1]])
                    X = X + x_pos
                    Y = Y + y_pos
                    X_vel = X_vel + x_pos
                    Y_vel = Y_vel + y_pos
                    if args.yt_proj == False:
                        sim_file = usable_sim_files[frame_val][:-12] + 'part' + usable_sim_files[frame_val][-5:]
                    else:
                        sim_file = part_file
                    print "SIM_FILE =", sim_file
                    f.close()
                    f = h5py.File(sim_file, 'r')
                    if len(part_info['particle_mass']) == 1:
                        part_ind = 0
                    else:
                        min_dist = 1000.0
                        for part in range(len(part_info['particle_mass'])):
                            temp_pos = np.array([f[f.keys()[11]][part][13]/c['au'], f[f.keys()[11]][part][13+y_int]/c['au']])
                            dist = np.sqrt(np.abs(np.diff((temp_pos - pos)**2)))[0]
                            if dist < min_dist:
                                min_dist = dist
                                part_ind = part
                    
                    center_vel = [f[f.keys()[11]][part_ind][18], f[f.keys()[11]][part_ind][19], f[f.keys()[11]][part_ind][20]]
                    print "CENTER_VEL=", center_vel
                    f.close()
                    f = h5py.File(usable_files[frame_val], 'r')
                xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
                if args.axis == 'xy':
                    center_vel=center_vel[:2]
                else:
                    center_vel=center_vel[::2]
                
                if args.ax_lim != None:
                    if has_particles and args.image_center != 0:
                        xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center - 1], args.ax_lim + part_info['particle_position'][0][args.image_center - 1]]
                        ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center - 1], args.ax_lim + part_info['particle_position'][1][args.image_center - 1]]
                    else:
                        xlim = [-1*args.ax_lim, args.ax_lim]
                        ylim = [-1*args.ax_lim, args.ax_lim]

                if args.yt_proj == False:
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
                else:
                    if args.image_center == 0 or has_particles == False:
                        center_pos = np.array([0.0, 0.0, 0.0])
                    else:
                        center_pos = np.array([dd['particle_posx'][args.image_center-1].in_units('AU'), dd['particle_posy'][args.image_center-1].in_units('AU'), dd['particle_posz'][args.image_center-1].in_units('AU')])
                    #CALCULATE PARTICLE PROJECTION
                    if has_particles:
                        L_orth = np.array([L[1], -1*L[0]])
                        L_len = np.sqrt(L_orth[0]**2. + L_orth[1]**2.)
                        r = np.array([part_info['particle_position'][0], part_info['particle_position'][1]])
                        r_pos = r * (L_orth/L_len)

                        #part_plane_position = np.array([dd['particle_posx'].in_units('AU'), dd['particle_posy'].in_units('AU')])
                        #part_info['particle_position'][0] = np.sign(part_plane_position[0])*np.sqrt((part_plane_position[0])**2. + (part_plane_position[1])**2.)
                    x_width = (xlim[1] -xlim[0])
                    y_width = (ylim[1] -ylim[0])
                    thickness = yt.YTArray(args.slice_thickness, 'AU')

                    temp = dd['velx']
                    temp = dd['vely']
                    temp = dd['velz']
                    if has_particles:
                        temp = dd['particle_posx']
                        temp = dd['particle_posy']
                    temp = dd['velocity_magnitude']
                        
                    del temp
                    
                    proj = yt.OffAxisProjectionPlot(f, L, [simfo['field'], 'Projected_Velocity_mw', 'velz_mw', 'Projected_Magnetic_Field_mw', 'magz_mw', 'cell_mass'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))
                    image = (proj.frb.data[('flash', 'dens')]/thickness.in_units('cm')).T.value
                    velx_full = (proj.frb.data[('gas', 'Projected_Velocity_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).T.value
                    vely_full = (proj.frb.data[('gas', 'velz_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).T.value
                    magx = (proj.frb.data[('gas', 'Projected_Magnetic_Field_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).T.value
                    magy = (proj.frb.data[('gas', 'magz_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).T.value
                    mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm')).T.value
                    image = image.T
                    velx_full = velx_full.T
                    vely_full = vely_full.T
                    magx = magx.T
                    magy = magy.T
                    mass = mass.T
                    '''
                    if np.median(image[200:600,400]) > 5.e-15:
                        image = image.T
                        velx_full = velx_full.T
                        vely_full = vely_full.T
                        magx = magx.T
                        magy = magy.T
                        mass = mass.T
                    '''
            
                    velx_full = velx_full/mass
                    vely_full = vely_full/mass
                    magx = magx/mass
                    magy = magy/mass
                    del mass
                        
                    velx, vely = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=center_vel)
                    del velx_full
                    del vely_full
    
                    pickle_file = path + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
                    file = open(pickle_file, 'w+')
                    pickle.dump((X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel), file)
                    file.close()
                    print "Created Pickle:", pickle_file, "for  file:", usable_files[frame_val]

                f.close()
                    
            title_parts = args.title.split('_')
            title = ''
            for part in title_parts:
                if part != title_parts[-1]:
                    title = title + part + ' '
                else:
                    title = title + part
            
            if args.pickle_dump == False:
                plt.clf()
                fig, ax = plt.subplots()
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
                    #ax.annotate('$t$='+str(int(time_val))+'yr', xy=(xlim[0]+0.01*(xlim[1]-xlim[0]), ylim[1]-0.03*(ylim[1]-ylim[0])), va="center", ha="left", color='w', fontsize=args.text_font)
                    
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                cbar = plt.colorbar(plot, pad=0.0)
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
                ax.set_ylabel(yabel, labelpad=-20, fontsize=args.text_font)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

                plt.tick_params(axis='both', which='major', labelsize=16)
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
                #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    
                #plt.savefig(file_name + ".jpg", format='jpeg', bbox_inches='tight')
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
                args_dict.update({'annotate_velocity': args.plot_velocity_legend})
                args_dict.update({'time_val': time_val})
                args_dict.update({'cbar_min': cbar_min})
                args_dict.update({'cbar_max': cbar_max})
                args_dict.update({'title': title})
                args_dict.update({'yabel': yabel})
                args_dict.update({'axlim':args.ax_lim})
                args_dict.update({'xlim':xlim})
                args_dict.update({'ylim':ylim})
                print "Built dictionary"
                pickle_file = save_dir + 'movie_pickle.pkl'
                print "Got pickle file name"
                file = open(pickle_file, 'w+')
                print "Opened pickle file"
                pickle.dump((usable_files[frame_val], X, Y, X_vel, Y_vel, image, velx, vely, part_info, args_dict, simfo, args), file)
                print "Dumped data into pickle"
                file.close()
                print "Created pickle"
        
        rit = rit +1
        if rit == size:
            rit = 0

    print "completed making movie frames on rank", rank

if __name__ == '__main__': main()

#Read in inputs:


