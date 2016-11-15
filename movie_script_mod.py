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

def define_constants():
    constants = {'year':31557600.0, 'au':1.496e13, 'Msun':1.98841586e+33}
    return constants

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--zoom", help="Will movie be zoomed in?", default=False, type=bool)
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default = 4)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 2, type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    parser.add_argument("-st", "--start_time", help="What time woudl you like to start calculating times from?", type=int)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", default = 25)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-plr", "--plot_lref", help="would you like to annotate the refinement level?", default=False)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", default=True)
    parser.add_argument("-sc", "--smooth_cells", help="how many cells would you like the smooth the velocities over? If not defined it is set to half the annotatation frequency")
    parser.add_argument("-ms", "--magnetic_smoothing", help="how many cells would you like the smooth the magnetic field over? If not defined it is set to 1", default=1 )
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", default=True)
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-mt", "--movie_times", help="What movies times would you like plotted?", type=list, default=[])
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-16)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def get_files(path, args):
    if args.axis == "xz":
        type = "proj"
    else:
        type = "slice"
    source_directory = sorted(glob.glob(path + 'WIND_' + type + '*'))
    return source_directory

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
    lref =  path.split('lref')[-1].split('/')[0].split('_')[-1]
    if '0.' in  path.split('lref')[0].split('/')[-1]:
        den_pert = path.split('lref')[0].split('/')[-1].split('_')[-2]
    else:
        den_pert = '0.0'
    movie_type = path.split('/')[-2].split('_')[0]
    f = h5py.File(file, 'r')
    if 'r_accretion' in f.keys():
        racc = f['r_accretion'][0]/c['au']
    else:
        racc = 0.0
    for key in f.keys():
        if args.field in key:
            field = key
    dim = np.shape(f[field])[0]
    if args.zoom:
        zoom_cell = (dim - dim/float(args.zoom_times))/2.
    else:
        zoom_cell = 0
    xmin = f['minmax_xyz'][0][0]/c['au']
    xmax = f['minmax_xyz'][0][1]/c['au']
    cl = (xmax-xmin)/dim
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
                }
    f.close()
    return sim_info

#load initial file:
def initialise_grid(args, sim_info):
    c = define_constants()
    cl = sim_info['cell_length']
    xmin = sim_info['xmin']
    xmax = sim_info['xmax']
    ymin = sim_info['ymin']
    ymax = sim_info['ymax']
    zoom_cell = sim_info['zoom_cell']
    annotate_vel_freq = sim_info['annotate_freq']
    x = np.arange(xmin, xmax, cl)
    y = np.arange(ymin, ymax, cl)
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 31:
        #val = xmin + (zoom_cell*cl) + cl*(annotate_vel_freq/2.) + cl*(annotate_vel_freq*counter)
        val = annotate_vel_freq*counter + annotate_vel_freq/2.
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    x_vel = []
    y_vel = []
    for x_counter in x_ind:
        val = xmin + x_counter*cl
        x_vel.append(val)
        y_vel.append(val)
    x_vel = np.array(x_vel)
    y_vel = np.array(y_vel)
    #print x_vel
    #x_vel = np.arange(xmin+((zoom_cell)*cl), xmax-(zoom_cell)*cl, annotate_vel_freq*cl)
    #y_vel = np.arange(ymin+((zoom_cell)*cl), ymax-(zoom_cell)*cl, annotate_vel_freq*cl)
    X, Y = np.meshgrid(x, y)
    X_vel, Y_vel = np.meshgrid(x_vel, y_vel)
    #print "created meshs"
    return X, Y, X_vel, Y_vel

def has_sinks(f):
    if "particlepositions" in f.keys():
        return True
    else:
        return False

def get_image_arrays(f, field, simfo, args):
    dim = int(simfo['dimension'])
    image = []
    for x in range(int(simfo['zoom_cell']), int(simfo['dimension']-simfo['zoom_cell'])):
        image_val = f[field][x]
        if np.shape(image_val)[0] == 1:
            image_val = image_val.transpose()
        image_val = image_val[simfo['zoom_cell']: simfo['dimension']-simfo['zoom_cell']]
        if simfo['movie_file_type'] == "proj":
            image_val = image_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
        image.append(image_val)
    image = np.array(image)
    if np.shape(image)[-1] == 1:
        image = image[:, :, 0]
    return image

def get_quiver_arrays(f, field, sim_info, args):
    smoothing_val = int(sim_info['smoothing'])
    x_pos = []
    y_pos = []
    counter = 0
    while counter < 31.:
        val = sim_info['zoom_cell'] + sim_info['annotate_freq']*counter + sim_info['annotate_freq']/2.
        x_pos.append(int(val))
        y_pos.append(int(val))
        counter = counter + 1
    #print x_pos
    velx = []
    vely = []
    '''
    x_test = f["vel" + args.axis[0] + "_" + sim_info['type'] + "_" + args.axis][0]
    if np.shape(x_test) == (sim_info['dimension'], 1):
        Transpose = True
    else:
        Transpose = False
    '''
    for x in x_pos:
        x = int(x)
        xarr = []
        yarr = []
        for y in y_pos:
            y = int(y)
            x_vel = 0
            y_vel = 0
            average_num = len(range(x-smoothing_val, x+smoothing_val-1))
            for curx in range(x-smoothing_val, x+smoothing_val-1):
                #for cury in range(y-smoothing_val, y+smoothing_val):
                tempx = f["vel" + args.axis[0] + "_" + sim_info['type'] + "_" + args.axis][curx]
                tempy = f["vel" + args.axis[1] + "_" + sim_info['type'] + "_" + args.axis][curx]
                '''
                if Transpose:
                    tempx = tempx.transpose()
                    tempy = tempy.transpose()
                '''
                x_vel = x_vel + sum(tempx[y-smoothing_val:y+smoothing_val-1])#[cury]
                y_vel = y_vel + sum(tempy[y-smoothing_val:y+smoothing_val-1])#[cury]
            x_vel = x_vel/np.square(average_num)
            y_vel = y_vel/np.square(average_num)
            xarr.append(x_vel)
            yarr.append(y_vel)
        velx.append(xarr)
        vely.append(yarr)
    velx = np.array(velx)
    vely = np.array(vely)
    return velx, vely

def get_particle_data(f, args, sim_info):
    c = define_constants()
    if args.axis == "xy":
        part_pos_x = f["particlepositions"][0]/c['au']
        part_pos_y = f["particlepositions"][1]/c['au']
    else:
        part_pos_x = f["particlepositions"][0]/c['au']
        part_pos_y = f["particlepositions"][2]/c['au']
    part_mass = np.array(f["particlemasses"])/c['Msun']
    return part_mass, [part_pos_x, part_pos_y]

def image_properties(X, Y, args, sim_info):
    if args.axis == "xy":
        ylabel = '$y$ (AU)'
    else:
        ylabel = '$z$ (AU)'
    xlim = [sim_info['xmin'], sim_info['xmax']-sim_info['cell_length']]
    ylim = [sim_info['ymin'], sim_info['ymax']-sim_info['cell_length']]
    return ylabel, xlim, ylim

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

def my_own_quiver_function(axis, X_pos, Y_pos, X_val, Y_val, sim_info, args):
    standard_vel = 500000.
    scale_factor = 100000.
    vel_scale_factor = 90000.
    #len_scale = 4000.
    len_scale = standard_vel/(0.04*(sim_info['xmax']-sim_info['xmin']))
    vels = np.hypot(X_val, Y_val)
    for xp in range(len(X_pos[0])):
        for yp in range(len(Y_pos[0])):
            xvel = X_val[xp][yp]/len_scale
            yvel = Y_val[xp][yp]/len_scale
            #width_val = widths[xp][yp]
            width_val = np.sqrt(X_val[xp][yp]**2. + Y_val[xp][yp]**2.)/standard_vel
            if width_val > 1.0:
                width_val = 1.0
            axis.add_patch(mpatches.FancyArrowPatch((X_pos[xp][yp], Y_pos[xp][yp]), (X_pos[xp][yp]+xvel, Y_pos[xp][yp]+yvel), color='w', linewidth=1.*width_val, arrowstyle='->', mutation_scale=15.*width_val, shrinkA=0.0, shrinkB=0.0))
    #print "done plotting arrows"
    if args.plot_velocity_legend == True:
        if args.zoom_times and float(args.zoom_times) != 1.28:
            pos_start = [0.77*sim_info['xmax'], 0.87*sim_info['ymin']]
        else:
            pos_start = [0.77*sim_info['xmax'], 0.87*sim_info['ymin']]
        xvel = standard_vel/len_scale
        yvel = 0.0
        width_val = 1.0
        axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=1.*width_val, mutation_scale=15.*width_val))
        axis.annotate("5kms$^{-1}$", xy=(0.98*sim_info['xmax'], 0.95*sim_info['ymin']), va="center", ha="right", color='w', fontsize=16)
    #print "plotted arrow legend"

#=======MAIN=======
def main():
    
    rank = CW.Get_rank()
    size = CW.Get_size()
    
    # Read in directories:
    path = sys.argv[1]
    save_dir = sys.argv[2]
    if os.path.exists(save_dir) == False:
        os.makedirs(save_dir)
    
    c = define_constants()
    args = parse_inputs()
    files = get_files(path, args)
    simfo = sim_info(path, files[-1], args)
    sink_form_time = mym.find_sink_formation_time(files)
    X, Y, X_vel, Y_vel = initialise_grid(args, simfo)
    yabel, xlim, ylim = image_properties(X, Y, args, simfo)
    deltax = (simfo['dimension'])/64.
    x_stream = np.arange(0, simfo['dimension'], deltax)
    start_stream = []
    for x in x_stream:
        start_stream.append([x, 0])
    start_stream = np.array(start_stream)

    # Define colourbar bounds
    cbar_max = args.colourbar_max
    cbar_min = args.colourbar_min
    
    # Initialise Grid and build lists
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step)
    no_frames = len(m_times)
    sys.stdout.flush()
    CW.Barrier()
    
    #if rank == 0:
    usable_files = mym.find_files(m_times, files)
    sys.stdout.flush()
    CW.Barrier()
    frames = range(args.start_frame, no_frames)

    sys.stdout.flush()
    CW.Barrier()
    rit = args.working_rank
    for frame_val in frames:
        if rank == rit:
            f = h5py.File(usable_files[frame_val], 'r')
            file_time = (f['time'][0]/c['year'])-sink_form_time
            has_particles = has_sinks(f)
            image = get_image_arrays(f, simfo['field'], simfo, args)
            magx = get_image_arrays(f, 'mag'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args)
            magy = get_image_arrays(f, 'mag'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args)
            #magx = smooth_array_2D(magx, args)
            #magy = smooth_array_2D(magy, args)
            velx, vely = get_quiver_arrays(f, 'vel', simfo, args)
            plt.clf()
            fig, ax = plt.subplots()
            plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
            cbar = plt.colorbar(plot, pad=0.0)
            cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=11, size=14)
            if frame_val > 0 or file_time > -1.0:
                plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)#, start_points=start_stream)
            #my_own_streamplot_function(X, Y, magx, magy, start_stream, sim_info, args)
            else:
                plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5)#, start_points=start_stream)
            #my_own_streamplot_function(X, Y, magx, magy, start_stream, sim_info, args)
#Q = plt.quiver(X_vel, Y_vel, velx, vely, width=0.00003, linewidths=proj_mag_lw, edgecolors='w', headwidth=50, headlength=75, pivot='mid', scale=1.e7,minlength=30)#linewidths=proj_mag_lw, edgecolors='w', minlength=0.0, headwidth=10, headlength=12) #scale=1.e7

#qk = quiverkey(Q, 0.9, 0.05, 500000, r'$5\rm{km}\rm{s}^{-1}$', labelpos='N', color='w', labelcolor='w',fontproperties={'weight': 'bold'}, labelsep=0.05)
            my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, simfo, args)
            if has_particles:
                part_mass, part_pos = get_particle_data(f, args, simfo)
                part_color = ['lime','cyan','r','c','y','w','k']
                if args.zoom and float(args.zoom_times) != 1.28:
                    line_rad = simfo['r_accretion']
                elif simfo['refinement_level'] < 12:
                    line_rad = simfo['r_accretion']
                elif float(args.zoom_times) == 1.28:
                    line_rad = 0.025*simfo['xmax']
                else:
                    line_rad = 0.025*simfo['xmax']
                #print "line_rad=", line_rad
                particle_text = []
                p_t = ''
                for pos_it in range(len(part_pos[0])):
                    #ax.plot((part_pos[0][pos_it]-(line_rad), part_pos[0][pos_it]+(line_rad)), (part_pos[1][pos_it], part_pos[1][pos_it]), lw=1, color='k')
                    #ax.plot((part_pos[0][pos_it], part_pos[0][pos_it]), (part_pos[1][pos_it]-(line_rad), part_pos[1][pos_it]+(line_rad)), lw=1, color='k')
                    ax.plot((part_pos[0][pos_it]-(line_rad), part_pos[0][pos_it]+(line_rad)), (part_pos[1][pos_it], part_pos[1][pos_it]), lw=2., c='k')
                    ax.plot((part_pos[0][pos_it], part_pos[0][pos_it]), (part_pos[1][pos_it]-(line_rad), part_pos[1][pos_it]+(line_rad)), lw=2., c='k')
                    ax.plot((part_pos[0][pos_it]-(line_rad), part_pos[0][pos_it]+(line_rad)), (part_pos[1][pos_it], part_pos[1][pos_it]), lw=1., c=part_color[pos_it])
                    ax.plot((part_pos[0][pos_it], part_pos[0][pos_it]), (part_pos[1][pos_it]-(line_rad), part_pos[1][pos_it]+(line_rad)), lw=1., c=part_color[pos_it])
                    circle = mpatches.Circle([part_pos[0][pos_it], part_pos[1][pos_it]], simfo['r_accretion'], fill=False, lw=1, edgecolor='k')
                    ax.add_patch(circle)
                    P_msun = str(np.round(part_mass[pos_it], decimals=2))
                    if p_t == '':
                        p_t = '$M_'+str(pos_it+1)+'$='+P_msun+'M$_\odot$'
                    else:
                        p_t = p_t+', $M_'+str(pos_it+1)+'$='+P_msun+'M$_\odot$'
                particle_text = p_t.split(',')
                ax.annotate(p_t, xy=(0.98*simfo['xmin'], 0.95*simfo['ymin']), va="center", ha="left", color='w', fontsize=16)
#plt.text((1./2.)*(simfo['xmax']-((simfo['zoom_cell'])*simfo['cell_length'])), (simfo['xmin']+((simfo['zoom_cell'])*simfo['cell_length']))-((1./5.)*(simfo['xmax']-((simfo['zoom_cell'])*simfo['cell_length']))), r''+particle_text)
#rainbow_text((1./2.)*(simfo['xmax']-((simfo['zoom_cell'])*simfo['cell_length'])), (simfo['xmin']+((simfo['zoom_cell'])*simfo['cell_length']))-((1./5.)*(simfo['xmax']-((simfo['zoom_cell'])*simfo['cell_length']))),particle_text, part_color)
            f.close()
            time_val = 10.0*(np.floor(np.round(file_time)/10.0))
            time_val = m_times[frame_val]
            if time_val == -0.0:
                time_val = 0.0
            if args.plot_lref == True:
                r_acc = np.round(simfo['r_accretion'])
                ax.annotate('$r_{acc}$='+str(r_acc)+'AU', xy=(0.98*simfo['xmax'], 0.93*simfo['ymax']), va="center", ha="right", color='w', fontsize=16)
            if args.annotate_time == True:
                ax.annotate('$t$='+str(int(time_val))+'yr', xy=(0.98*simfo['xmin'], 0.93*simfo['ymax']), va="center", ha="left", color='w', fontsize=16)
            ax.annotate(args.title, xy=(0.0, 0.93*simfo['ymax']), va="center", ha="center", color='w', fontsize=20)
            ax.set_xlabel('$x$ (AU)', labelpad=-1)
            ax.set_ylabel(yabel, labelpad=-20)
                #if args.zoom:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
                #else:
                #ax.set_xlim([-1000.,1000])
                #ax.set_ylim([-1000.,1000])
            #plt.show()
            if len(usable_files) > 1:
                if args.output_filename == None:
                    file_name = save_dir + "movie_frame_" + ("%06d" % frame_val)
                else:
                    file_name = args.output_filename + "_" + str(int(time_val))
            else:
                if args.output_filename != None:
                    file_name = args.output_filename
                else:
                    file_name = save_dir + "time_" + str(args.plot_time)
            plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight', pad_inches = 0.02)
            plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight', pad_inches = 0.02)
                
                #plt.savefig(file_name + ".jpg", format='jpeg', bbox_inches='tight')
            call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
            #os.remove(file_name + '.eps')
            print 'Created frame', (frame_val+1), 'of', str(len(usable_files)), 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', save_dir+file_name
        rit = rit +1
        if rit == size:
            rit = 0

    print "completed making movie frames on rank", rank

if __name__ == '__main__': main()

#Read in inputs:


