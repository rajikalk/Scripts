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

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

def define_constants():
    constants = {'year':31557600.0, 'au':1.496e13, 'Msun':1.98841586e+33}
    return constants

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default = 0)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 2, type=float)
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
    parser.add_argument("-mt", "--movie_times", help="What movies times would you like plotted?", type=list, default=[])
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-16)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("-ic", "--image_center", help="where would you like to center the image?", type=int, default=0)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=None)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
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
    lref = path.split('lref')[-1].split('/')[0].split('_')[-1]
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

def has_sinks(f):
    if "particlepositions" in f.keys():
        return True
    else:
        return False

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
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
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
    mym.set_global_font_size(args.text_font)
    files = get_files(path, args)
    sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
    
    # Initialise Grid and build lists
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=args.presink_frames)
    no_frames = len(m_times)
    sys.stdout.flush()
    CW.Barrier()

    #if rank == 0:
    usable_files = mym.find_files(m_times, files)
    usable_sim_files = mym.find_files(m_times, sim_files)
    sys.stdout.flush()
    CW.Barrier()
    frames = range(args.start_frame, no_frames)

    sink_form_time = mym.find_sink_formation_time(files)
    print "sink_form_time", sink_form_time

    # Define colourbar bounds
    cbar_max = args.colourbar_max
    cbar_min = args.colourbar_min

    sys.stdout.flush()
    CW.Barrier()
    rit = args.working_rank
    for frame_val in frames:
        if rank == rit:
            f = h5py.File(usable_files[frame_val], 'r')
            print "FILE =", usable_files[frame_val]
            file_time = (f['time'][0]/c['year'])-sink_form_time
            simfo = sim_info(path, usable_files[frame_val], args)
            has_particles = has_sinks(f)
            if has_particles:
                part_info = mym.get_particle_data(usable_files[frame_val], args.axis)
            X, Y, X_vel, Y_vel, cl = mym.initialise_grid(usable_files[frame_val], zoom_times=args.zoom_times)
            center_vel = [0.0, 0.0, 0.0]
            if args.image_center != 0:
                sim_file = usable_sim_files[frame_val][:-12] + 'part' + usable_sim_files[frame_val][-5:]
                print "SIM_FILE =", sim_file
                f.close()
                f = h5py.File(sim_file, 'r')
                #center_vel = [f[f.keys()[11]][args.image_center-1][18], f[f.keys()[11]][args.image_center-1][19], f[f.keys()[11]][args.image_center - 1][20]]
                center_vel = [f[f.keys()[11]][0][18], f[f.keys()[11]][0][19], f[f.keys()[11]][0][20]]

                import pdb
                pdb.set_trace()
                print "CENTER_VEL=", center_vel
                f.close()
                f = h5py.File(usable_files[frame_val], 'r')
                x_pos = np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                y_pos = np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                X = X + x_pos
                Y = Y + y_pos
                X_vel = X_vel + x_pos
                Y_vel = Y_vel + y_pos
            yabel, xlim, ylim = image_properties(X, Y, args, simfo)
            
            if args.ax_lim != None:
                if args.image_center == 0:
                    xlim = [-1*args.ax_lim, args.ax_lim]
                    ylim = [-1*args.ax_lim, args.ax_lim]
                else:
                    xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center - 1], args.ax_lim + part_info['particle_position'][0][args.image_center - 1]]
                    ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center - 1], args.ax_lim + part_info['particle_position'][1][args.image_center - 1]]

            image = get_image_arrays(f, simfo['field'], simfo, args, X, Y)
            print "image shape=", np.shape(image)
            print "grid shape=", np.shape(X)
            magx = get_image_arrays(f, 'mag'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args, X, Y)
            magy = get_image_arrays(f, 'mag'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis, simfo, args, X, Y)
            x_pos_min = int(np.round(np.min(X) - simfo['xmin_full'])/simfo['cell_length'])
            y_pos_min = int(np.round(np.min(Y) - simfo['xmin_full'])/simfo['cell_length'])
            if args.axis == 'xy':
                if args.image_center != 0:
                    velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X, f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0], center_vel=center_vel[:2])
                else:
                    velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X, f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,:,0])
            else:
                if args.image_center != 0:
                    velx, vely = mym.get_quiver_arrays(y_pos_min, y_pos_min, X, f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,0,:], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,0,:], center_vel=center_vel[::2])
                else:
                    velx, vely = mym.get_quiver_arrays(y_pos_min, y_pos_min, X, f['vel'+args.axis[0]+'_'+simfo['movie_file_type']+'_'+args.axis][:,0,:], f['vel'+args.axis[1]+'_'+simfo['movie_file_type']+'_'+args.axis][:,0,:])
            time_val = 10.0*(np.floor(np.round(file_time)/10.0))
            time_val = m_times[frame_val]
            if time_val == -0.0:
                time_val = 0.0
            title_parts = args.title.split('_')
            title = ''
            for part in title_parts:
                if part != title_parts[-1]:
                    title = title + part + ' '
                else:
                    title = title + part
            f.close()
            
            if args.pickle_dump == False:
                plt.clf()
                fig, ax = plt.subplots()
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
                plt.gca().set_aspect('equal')
                if frame_val > 0 or file_time > -1.0:
                    plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                else:
                    plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5)
                
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)
                if has_particles:
                    if args.annotate_particles_mass == True:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'])
                    else:
                        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)
                if args.plot_lref == True:
                    r_acc = np.round(part_info['accretion_rad'])
                    ax.annotate('$r_{acc}$='+str(r_acc)+'AU', xy=(0.98*simfo['xmax'], 0.93*simfo['ymax']), va="center", ha="right", color='w', fontsize=args.text_font)
                if args.annotate_time == "True":
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(time_val))+'yr', va="center", ha="left", color='w', fontsize=args.text_font)
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                    #ax.annotate('$t$='+str(int(time_val))+'yr', xy=(xlim[0]+0.01*(xlim[1]-xlim[0]), ylim[1]-0.03*(ylim[1]-ylim[0])), va="center", ha="left", color='w', fontsize=args.text_font)
                ax.annotate(title, xy=(np.mean(xlim), ylim[1]-0.03*(ylim[1]-ylim[0])), va="center", ha="center", color='w', fontsize=(args.text_font+4))

                cbar = plt.colorbar(plot, pad=0.0)
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                ax.set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
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
                        file_name = save_dir + "movie_frame_" + ("%06d" % frame_val)
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
                print 'Created frame', (frame_val+1), 'of', str(len(usable_files)), 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.eps'
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


