import matplotlib.pyplot as plt
import pickle
import csv
import numpy as np
import os
from matplotlib.colors import LogNorm
import h5py
import my_module as mym
import matplotlib.gridspec as gridspec
import glob
import matplotlib.patheffects as path_effects
import matplotlib
import yt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", "--in_file", help="input file with images, positions and crops")
    parser.add_argument("-sx", "--share_x", help="do you want to share the x axis?", default=False)
    parser.add_argument("-sy", "--share_y", help="do you want to share the y axis?", default=True)
    parser.add_argument("-sa", "--share_ax", help="do you want to share axes?", default=True)
    parser.add_argument("-sd", "--save_dir", help="Where do you want to save the plot? defaults to same as input file directory", default='./')
    parser.add_argument("-sc", "--share_colourbar", help="Do you want to share colour bars over rows?", default=False)
    parser.add_argument("-tf", "--text_font", help="what font do you want the text to have?", type=int, default=12)
    parser.add_argument("-sn", "--savename", help="what do you want to save the plot as?", type=str, default='multiplot')
    parser.add_argument("-ylp", "--y_label_pad", help="y lable pad", default=-20, type=float)
    parser.add_argument("-title", "--multiplot_title", help="Do you want to title the multiplot", type=str, default="")
    args = parser.parse_args()
    return args

def get_image_arrays(f, field, simfo, args, X, Y):
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


#=============================MAIN==============================
args = parse_inputs()
if args.share_x == 'True':
    args.share_x = True
if args.share_y == 'False':
    args.share_y = False
elif args.share_y == 'True':
    args.share_y = True
if args.share_ax == 'False':
    args.share_ax = False
if args.share_colourbar == 'False':
    args.share_colourbar = False
if args.save_dir == None:
    path = os.path.abspath(args.in_file).split('/')[:-1]
    save_dir = ''
    for dir in path:
        save_dir = save_dir + dir + '/'
else:
    save_dir = args.save_dir
savename = save_dir + args.savename

mym.set_global_font_size(args.text_font)

positions = []
plot_type = []
file_dir = []
input_args = []

with open(args.in_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0].lower() == 'grid_inputs:':
            glr = float(row[1])
            grl = float(row[2])
            glw = float(row[3])
            ghspace = float(row[4])
        elif row[0][0]!= '#':
            positions.append((int(row[0]), int(row[1])))
            plot_type.append(row[2])
            file_dir.append(row[3])
            input_args.append(row[4])
positions = np.array(positions)

columns = np.max(positions[:,0])
rows = np.max(positions[:,1])

width = float(columns)*(14.5/3.)
height = float(rows)*(17./4.)
if plot_type[0] == "outflow":
    f = plt.figure(figsize=(6, 10))
else:
    f = plt.figure(figsize=(width, height))
gs_left = gridspec.GridSpec(rows, columns-1)
gs_right = gridspec.GridSpec(rows, 1)

gs_left.update(right=glr, wspace=glw, hspace=ghspace)
gs_right.update(left=grl, hspace=ghspace)

'''
#single and tight binary slices
gs_left.update(right=0.60)
gs_right.update(left=0.3775)
gs_left.update(wspace=-0.0625)
gs_left.update(hspace=0.05)
gs_right.update(hspace=0.05)
'''
'''
#profile plots
gs_left.update(right=0.60)
gs_right.update(left=0.3475)
gs_left.update(wspace=0.01)
gs_left.update(hspace=0.05)
gs_right.update(hspace=0.05)
'''
'''
#wide slices
gs_left.update(right=0.60)
gs_right.update(left=0.3775)
gs_left.update(wspace=-0.0625)
gs_left.update(hspace=0.05)
gs_right.update(hspace=0.05)
cbar_plotted = False
'''
axes_dict = {}

for it in range(len(positions)):
    #ax_label = 'ax' + str(counter)
    ax_label = 'ax' + str(positions[it][0])+str(positions[it][1])
    if positions[it][0] == 1 and positions[it][1] == 1:
        if columns > 1:
            axes_dict.update({ax_label:f.add_subplot(gs_left[0,0])})
        else:
            axes_dict.update({ax_label:f.add_subplot(gs_right[0,0])})
    elif positions[it][0] != columns:
        if args.share_x and args.share_y == True:
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'], sharey=axes_dict['ax11'])})
        elif args.share_x and args.share_y == 'row':
            if positions[it][0] == 1:
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'])})
            else:
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'], sharey=axes_dict['ax1'+str(positions[it][1])])})
        elif args.share_x:
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'])})
        else:
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])})
    else:
        if args.share_x and args.share_y == True:
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'], sharey=axes_dict['ax11'])})
        elif args.share_x and args.share_y == 'row':
            if positions[it][0] == 1:
                axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'])})
            else:
                axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'], sharey=axes_dict['ax1'+str(positions[it][1])])})
        elif args.share_x:
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'])})
        else:
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0])})
    '''
    yit = np.where(positions[:,1] == positions[it][1])[0][0]
    if positions[it][0] == 1 and positions[it][1] == 1:
        if columns > 1:
            axes_dict.update({ax_label:f.add_subplot(gs_left[0,0])})
        else:
            axes_dict.update({ax_label:f.add_subplot(gs_right[0,0])})
    elif positions[it][0] != columns:
        if args.share_x and args.share_y:
            if yit >= len(axes_dict):
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax1'])})
            else:
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax1'], sharey=axes_dict[axes_dict.keys()[yit]])})
        elif args.share_x:
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax1'])})
        elif args.share_y and positions[it][0]!=1:
            yit = np.where(positions[:,1] == positions[it][1])[0][0]
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharey=axes_dict[axes_dict.keys()[yit]])})
        elif args.share_y:
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])})
        elif args.share_y == 'row':
            yit = np.where(positions[:,1] == positions[it][1])[0][0]
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharey=axes_dict[axes_dict.keys()[yit]])})
        else:
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])})
    else:
        if args.share_x and args.share_y:
            yit = np.where(positions[:,1] == positions[it][1])[0][0]
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax1'], sharey=axes_dict[axes_dict.keys()[yit]])})
        elif args.share_x:
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax1'])})
        elif args.share_y:
            yit = np.where(positions[:,1] == positions[it][1])[0][0]
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharey=axes_dict[axes_dict.keys()[yit-1]])})
        else:
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0])})
    '''
    if 'amb' in plot_type[it]:
        file = file_dir[it]
        times = []
        ang_labels = []
        ang_arrays = []
        header = True
        with open(file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if header == False:
                    times.append(float(row[0]))
                    for col_it in range(len(row[1:])):
                        ang_arrays[col_it].append(float(row[col_it+1]))
                if header == True:
                    for col in row[1:]:
                        ang_labels.append(col)
                        ang_arrays.append([])
                    header = False
        for arr in range(len(ang_arrays)):
            plt.semilogy(times, ang_arrays[arr], label=ang_labels[arr])
    if 'movie' in plot_type[it]:
        #axes_dict[ax_label].set(adjustable='box-forced', aspect='equal')
        mov_args = input_args[it].split(' ')
        arg_list = ['python', '/groups/astro/rlk/Scripts/movie_script_mod.py', file_dir[it], save_dir, '-pd', 'True', '-tf', str(args.text_font)] #['srun', '-n', '20', 'python', '/groups/astro/rlk/Scripts/movie_script_mod.py', file_dir[it], save_dir, '-pd', 'True', '-tf', str(args.text_font)]
        get_stdv = False
        standard_vel = 5.
        get_plot_time = False
        plot_time = 0
        field_str = 'dens'
        get_field = False
        weight_field = 'dens'
        get_weight_field = False
        cbar_lims = [None, None]
        get_cbar_lim = (False, None)
        ax_string = 'xz'
        get_axis_string = False
        ax_limits = None
        get_axis_limits = False
        title = ''
        get_title = False
        thickness = 300
        get_thickness = False
        for mov_arg in mov_args:
            if get_stdv:
                standard_vel = float(mov_arg)
                get_stdv = False
            if mov_arg == '-stdv':
                get_stdv = True
            if get_plot_time:
                plot_time = str(float(mov_arg))
                get_plot_time = False
            if mov_arg == '-pt':
                get_plot_time = True
            if get_field:
                field_str = mov_arg
                get_field = False
            if mov_arg == '-f':
                get_field = True
            if get_weight_field:
                weight_field = mov_arg
                get_weight_field = False
            if mov_arg == '-wf':
                get_weight_field = True
            if get_cbar_lim[0] == True:
                cbar_lims[get_cbar_lim[1]] = float(mov_arg)
                get_cbar_lim = (False, None)
            if mov_arg == '-cmin':
                get_cbar_lim = (True, 0)
            if mov_arg == '-cmax':
                get_cbar_lim = (True, 1)
            if get_axis_string:
                ax_string = mov_arg
                get_axis_string = False
            if mov_arg == '-ax':
                get_axis_string = True
            if get_axis_limits:
                ax_limits = float(mov_arg)
                get_axis_limits = False
            if mov_arg == '-al':
                get_axis_limits = True
            if get_title:
                title = mov_arg
                get_title = False
            if mov_arg == '-t':
                get_title = True
            if get_thickness:
                thickness = int(mov_arg)
                get_thickness = False
            if mov_arg == '-thickness':
                get_thickness = True
            arg_list.append(mov_arg)

        if weight_field == 'None':
            pickle_file = file_dir[it] + ax_string + '_' + field_str + '_thickness_' + str(thickness) + '_AU_movie_time_'+plot_time+'_unweighted.pkl'
        else:
            pickle_file = file_dir[it] + ax_string + '_' + field_str + '_thickness_' + str(thickness) + '_AU_movie_time_'+plot_time+'.pkl'
        if os.path.isfile(pickle_file) == False:
            os.system(" ".join(arg_list))

        file = open(pickle_file, 'rb')
        print("pickle file:", pickle_file)
        X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo = pickle.load(file)
        file.close()
        if cbar_lims != [None, None]:
            args_dict['cbar_min'] = cbar_lims[0]
            args_dict['cbar_max'] = cbar_lims[1]
        if 0.0 in (args_dict['cbar_min'], args_dict['cbar_max']):
            #plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max']), rasterized=True)
            if 'Relative_Keplerian_Velocity' in field_str or 'B_angle' in field_str:
                plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.brg, rasterized=True, vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max'])
            else:
                plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.RdBu, rasterized=True, vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max'])
        elif weight_field == 'None':
            plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max']), rasterized=True)
            #plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(), rasterized=True)
        else:
            print("plotted with log scale")
            plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max']), rasterized=True)
        axes_dict[ax_label].streamplot(X, Y, magx, magy, density=3, linewidth=0.5, minlength=0.5, arrowstyle='-', color='royalblue')
        mym.my_own_quiver_function(axes_dict[ax_label], X_vel, Y_vel, velx, vely, plot_velocity_legend=args_dict['annotate_velocity'], limits=[args_dict['xlim'], args_dict['ylim']], standard_vel=standard_vel)
        part_info['particle_mass'] = np.sort(part_info['particle_mass'])[::-1]
        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], [args_dict['xlim'], args_dict['ylim']], annotate_field=part_info['particle_mass'])
        if 'annotate_time' in list(args_dict.keys()):
            time_string  = args_dict['annotate_time'].split('=')[-1].split('yr')[0]
            time_string = str(int(np.round(float(time_string)/(10**(len(str(int(time_string)))-1)))*(10**(len(str(int(time_string)))-1))))
            time_text = axes_dict[ax_label].text((args_dict['xlim'][0]+0.01*(args_dict['xlim'][1]-args_dict['xlim'][0])), (args_dict['ylim'][1]-0.03*(args_dict['ylim'][1]-args_dict['ylim'][0])), "$t$="+time_string+"yr", va="center", ha="left", color='w', fontsize=args.text_font)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        if args.multiplot_title == "":
            if title != args_dict['title']:
                args_dict[title] = title
            title = axes_dict[ax_label].text(np.mean(args_dict['xlim']), (args_dict['ylim'][1]-0.04*(args_dict['ylim'][1]-args_dict['ylim'][0])), args_dict['title'], va="center", ha="center", color='w', fontsize=(args.text_font+2))
            title.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        axes_dict[ax_label].set_xlim(args_dict['xlim'])
        axes_dict[ax_label].set_ylim(args_dict['ylim'])
        for line in axes_dict[ax_label].xaxis.get_ticklines():
            line.set_color('white')
        for line in axes_dict[ax_label].yaxis.get_ticklines():
            line.set_color('white')
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel(args_dict['yabel'], labelpad=args.y_label_pad, fontsize=args.text_font)
        if args.share_colourbar == False:
            if positions[it][0] == columns:
                cbar = plt.colorbar(plot, pad=0.0, ax=axes_dict[ax_label])
                if field_str == 'dens':
                    cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'magnetic_field_poloidal':
                    cbar.set_label('$B_\mathrm{Pol}$ (gauss)', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'magnetic_field_toroidal':
                    cbar.set_label('$B_\mathrm{Tor}$ (gauss)', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'B_Pol_to_B_Tor':
                    cbar.set_label('$B_\mathrm{Pol}/B_\mathrm{Tor}$', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'B_Tor_to_B_mag':
                    cbar.set_label('$B_\mathrm{Tor}/B_\mathrm{mag}$', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'B_Pol_to_B_mag':
                    cbar.set_label('$B_\mathrm{Pol}/B_\mathrm{mag}$', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'B_angle':
                    cbar.set_label(r"$\theta$ ($^{\circ}$)", rotation=270, labelpad=15, size=args.text_font)
                else:
                    cbar.set_label('Relative Keplerian Velocity ($v_{\phi}/v_{\mathrm{kep}}$)', rotation=270, labelpad=15, size=args.text_font)
        else:
            if cbar_plotted == False:
                cax = f.add_axes([0.9, 0.1, 0.02, 0.8])
                f.colorbar(plot, pad=0.0, cax=cax)
                if field_str == 'dens':
                    cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'magnetic_field_poloidal':
                    cbar.set_label('$B_\mathrm{Pol}$ (gauss)', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'magnetic_field_toroidal':
                    cbar.set_label('$B_\mathrm{Tor}$ (gauss)', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'Pol_to_Tor_Ratio':
                    cbar.set_label('$B_\mathrm{Pol}/B_\mathrm{Tor}$', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'B_Tor_to_B_mag':
                    cbar.set_label('$B_\mathrm{Tor}/B_\mathrm{mag}$', rotation=270, labelpad=15, size=args.text_font)
                elif field_str == 'B_Tor_to_B_mag':
                    cbar.set_label('$B_\mathrm{Pol}/B_\mathrm{mag}$', rotation=270, labelpad=15, size=args.text_font)
                else:
                    cbar.set_label('Relative Keplerian Velocity ($v_{\phi}/v_{\mathrm{kep}}$)', rotation=270, labelpad=15, size=args.text_font)
                cbar_plotted = True
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('$x$ (AU)', fontsize=args.text_font)
        axes_dict[ax_label].set_aspect('equal')
        '''
        if positions[it][0] != 1 and ax_string!='xz':
            xticklabels = axes_dict[ax_label].get_xticklabels()
            plt.setp(xticklabels[1], visible=False)
        '''
        #axes_dict[ax_label].set_adjustable('box', share=True)
        print("added movie segment")
    if 'yt_proj' in plot_type[it]:
        axes_dict[ax_label].set(adjustable='box-forced', aspect='equal')
        yt_args = input_args[it].split(' ')
        arg_list = ['python', '/home/100/rlk100/Scripts/paper_plots.py', file_dir[it], save_dir, '-pd', 'True']
        for yt_arg in yt_args:
            arg_list.append(yt_arg)
        call(arg_list)

        pickle_file = save_dir + 'yt_proj_pickle.pkl'
        file = open(pickle_file, 'r')
        X, Y, X_vel, Y_vel, image, velx, vely, magx, magy, part_info, args_dict = pickle.load(file)
        file.close()

        plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max']), rasterized=True)
        axes_dict[ax_label].streamplot(X, Y, magx, magy, density=3, linewidth=0.5, minlength=0.5, arrowstyle='-', color='royalblue')
        mym.my_own_quiver_function(axes_dict[ax_label], X_vel, Y_vel, velx, vely, plot_velocity_legend=args_dict['annotate_velocity'], limits=[args_dict['xlim'], args_dict['ylim']], standard_vel=args_dict['standard_vel'])
        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], [args_dict['xlim'], args_dict['ylim']], annotate_field=part_info['particle_mass'])
        if 'annotate_time' in list(args_dict.keys()):
            time_text = axes_dict[ax_label].text((args_dict['xlim'][0]+0.01*(args_dict['xlim'][1]-args_dict['xlim'][0])), (args_dict['ylim'][1]-0.03*(args_dict['ylim'][1]-args_dict['ylim'][0])), args_dict['annotate_time'], va="center", ha="left", color='w', fontsize=args.text_font)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        title = axes_dict[ax_label].text(np.mean(args_dict['xlim']), (args_dict['ylim'][1]-0.04*(args_dict['ylim'][1]-args_dict['ylim'][0])), args_dict['title'], va="center", ha="center", color='w', fontsize=(args.text_font+2))
        title.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        axes_dict[ax_label].set_xlim(args_dict['xlim'])
        axes_dict[ax_label].set_ylim(args_dict['ylim'])
        for line in axes_dict[ax_label].xaxis.get_ticklines():
            line.set_color('white')
        for line in axes_dict[ax_label].yaxis.get_ticklines():
            line.set_color('white')
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel(args_dict['yabel'], labelpad=args.y_label_pad, fontsize=args.text_font)
        if args.share_colourbar == False:
            if positions[it][0] == columns:
                cbar = plt.colorbar(plot, pad=0.0, ax=axes_dict[ax_label])
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=15, size=args.text_font)
        else:
            if cbar_plotted == False:
                cax = f.add_axes([0.9, 0.1, 0.02, 0.8])
                f.colorbar(plot, pad=0.0, cax=cax)
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=15, size=args.text_font)
                cbar_plotted = True
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('Distance from center (AU)', fontsize=args.text_font)
    if 'slice' in plot_type[it]:
        axes_dict[ax_label].set(adjustable='box-forced', aspect='equal')
        slice_args = input_args[it].split(' ')
        arg_list = ['python', '/home/100/rlk100/Scripts/paper_plots.py', file_dir[it], save_dir, '-pd', 'True']
        plot_velocity_legend = False
        get_stdv = False
        standard_vel = 5.
        check_annotate_time = False
        annotate_time = False
        title = ''
        check_title = False
        update_slice_field = False
        slice_field = 'vkep'
        for slice_arg in slice_args:
            if get_stdv:
                standard_vel = float(slice_arg)
                get_stdv = False
            if slice_arg == '-stdv':
                get_stdv = True
            arg_list.append(slice_arg)
            if slice_arg == '-pvl':
                plot_velocity_legend = True
            if check_annotate_time and slice_arg == 'True':
                annotate_time = True
                check_annotate_time = False
            if slice_arg == '-at':
                check_annotate_time = True
            if check_title:
                title = slice_arg
                check_title = False
            if slice_arg == '-t':
                check_title = True
            if update_slice_field:
                slice_field = slice_arg
                update_slice_field = False
            if slice_arg == '-f':
                update_slice_field = True
        call(arg_list)

        pickle_file = save_dir + 'slice_pickle.pkl'
        file = open(pickle_file, 'r')
        X_vel, Y_vel, xy, field_grid, weight_field, velx, vely, magx_grid, magy_grid, part_info, limits, time_string = pickle.load(file)
        file.close()
        if slice_field == 'dens':
            plot = axes_dict[ax_label].pcolormesh(xy[0], xy[1], field_grid.value, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=1.e-15, vmax=1.e-12), rasterized=True)#cmap=plt.cm.get_cmap('bwr_r')
            axes_dict[ax_label].streamplot(xy[0], xy[1], magx_grid.value, magy_grid.value, density=4, linewidth=0.25, minlength=0.5, arrowstyle='-')
            if annotate_time:
                time_text = axes_dict[ax_label].text((np.min(xy[0])+0.01*(np.max(xy[0])-np.min(xy[0]))), (np.max(xy[1])-0.04*(np.max(xy[1])-np.min(xy[1]))), time_string, va="center", ha="left", color='w', fontsize=args.text_font)
                time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            title = axes_dict[ax_label].text(np.mean([np.min(xy[0]), np.max(xy[0])]), (np.max(xy[1])-0.08*(np.max(xy[1]))), title, va="center", ha="center", color='w', fontsize=(args.text_font+2))
            title.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], [[np.min(xy[0]), np.max(xy[0])], [np.min(xy[1]), np.max(xy[1])]], annotate_field=part_info['particle_mass'])
        else:
            img_array = plt.get_cmap('brg')(field_grid/np.max(field_grid))
            alphas = np.clip(weight_field/(0.01*np.max(weight_field)), 0.0, 1.0)**(1./10.)
            img_array[..., 3] = alphas
            background = np.ones(np.shape(field_grid))
            back_img_array = plt.get_cmap('Greys')(background)
            extent = np.min(xy[0]), np.max(xy[0]), np.min(xy[1]), np.max(xy[1])
            #plot = axes_dict[ax_label].pcolormesh(xy[0], xy[1], field_grid, cmap=plt.cm.get_cmap('brg'), rasterized=True, vmin=0.0, vmax=2.0)#cmap=plt.cm.get_cmap('bwr_r')
            axes_dict[ax_label].imshow(back_img_array,origin='lower', extent=extent)
            plot = axes_dict[ax_label].imshow(img_array,origin='lower', vmin=0.0, vmax=2.0, extent=extent, cmap=plt.cm.brg)
            mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits)
        mym.my_own_quiver_function(axes_dict[ax_label], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend, limits=limits, standard_vel=standard_vel)
        axes_dict[ax_label].set_xlim(limits[0])
        axes_dict[ax_label].set_ylim(limits[1])
        #plot.set_clim([0.0, 2.0])
        for line in axes_dict[ax_label].xaxis.get_ticklines():
            line.set_color('white')
        for line in axes_dict[ax_label].yaxis.get_ticklines():
            line.set_color('white')
        if positions[it][0] == columns and args.share_colourbar == False:
            cbar = plt.colorbar(plot, pad=0.0)
            if slice_field == 'dens':
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=22, size=args.text_font)
            else:
                cbar.set_label("Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)", rotation=270, labelpad=22, size=args.text_font)
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel('$y$ (AU)', labelpad=args.y_label_pad, fontsize=args.text_font)
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('$x$ (AU)', fontsize=args.text_font)
        axes_dict[ax_label].set_aspect((limits[1][1] - limits[1][0])/(limits[0][1] - limits[0][0]))
    if 'profile' in plot_type[it]:
        prof_args = input_args[it].split(' ')
        arg_list = ['python', '/home/100/rlk100/Scripts/paper_plots.py', file_dir[it], save_dir, '-pd', 'True']
        get_rmax = False
        r_max = 500.
        get_field = False
        field_str = 'Relative_Keplerian_Velocity'
        get_plot_time = False
        plot_time = 0
        is_log = False
        get_log = False
        title = ''
        get_title = False
        cal_col_dens = False
        get_col_mean = False
        for prof_arg in prof_args:
            arg_list.append(prof_arg)
            if get_rmax == True:
                r_max = float(prof_arg)
                get_rmax = False
            if 'rmax' in prof_arg:
                get_rmax = True
            if get_field:
                field_str = prof_arg
                get_field = False
            if prof_arg == '-f':
                get_field = True
            if get_plot_time:
                plot_time = str(float(prof_arg))
                get_plot_time = False
            if prof_arg == '-pt':
                get_plot_time = True
            if get_log:
                if prof_arg == "True":
                    is_log = True
                get_log = False
            if prof_arg == '-log':
                get_log = True
            if get_col_mean:
                if prof_arg == "True":
                    cal_col_dens = True
                get_col_mean = False
            if prof_arg == '-col_mean':
                get_col_mean = True

        if cal_col_dens:
            pickle_file = file_dir[it] + 'column_density_profile_pickle_' + str(int(float(plot_time))) + '.pkl'
        else:
            pickle_file = file_dir[it] + field_str + '_profile_pickle_' + str(int(float(plot_time))) + '.pkl'
        if os.path.isfile(pickle_file) == False:
            call(arg_list)

        file = open(pickle_file, 'r')
        prof_x, prof_y, sampled_points = pickle.load(file)
        file.close()

        if is_log:
            axes_dict[ax_label].semilogy(prof_x, prof_y, 'k-', linewidth=2.)
        else:
            axes_dict[ax_label].plot(prof_x, prof_y, 'k-', linewidth=2.)
        '''
        for sep in separation:
            axes_dict[ax_label].axvline(x=sep, alpha=0.5)
        '''
        axes_dict[ax_label].set_xlim([0.0, r_max])
        if field_str == 'Relative_Keplerian_Velocity':
            axes_dict[ax_label].set_ylim([0.0, 3.5])
            axes_dict[ax_label].axhline(y=1.0, color='k', linestyle='--')
            cm = plt.cm.get_cmap('RdYlBu')
            plot = axes_dict[ax_label].scatter(sampled_points[1], sampled_points[2], c=sampled_points[0], alpha=0.4, cmap=cm, edgecolors='none', vmin=0, vmax=100)
            if positions[it][0] == columns and args.share_colourbar == False:
                cbar = plt.colorbar(plot, pad=0.0)
                cbar.set_label('|$z$ position| (AU)', rotation=270, labelpad=20, size=args.text_font)
        if positions[it][0] == 1:
            if field_str != 'Relative_Keplerian_Velocity':
                axes_dict[ax_label].set_ylabel(field_str)
            else:
                axes_dict[ax_label].set_ylabel('Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)', size=args.text_font)
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('Radius (AU)', size=args.text_font)
        print("added profile segment")
    if 'multi' in plot_type[it]:
        prof_args = input_args[it].split(' ')
        get_field = False
        field_str = 'Tangential_Velocity'
        cal_col_dens = False
        get_col_mean = False
        title = ''
        get_title = False
        log_scale = True
        get_log_scale = False
        for prof_arg in prof_args:
            if get_field:
                field_str = prof_arg
                get_field = False
            if prof_arg == '-f':
                get_field = True
            if get_col_mean:
                if prof_arg == "True":
                    cal_col_dens = True
                get_col_mean = False
            if prof_arg == '-col_mean':
                get_col_mean = True
            if get_title:
                title = ' '.join(prof_arg.split('_'))
                get_title = False
            if prof_arg == '-t':
                get_title = True
            if get_log_scale:
                if prof_arg == 'False':
                    log_scale = False
                get_title = False
            if prof_arg == '-log':
                get_log_scale = True
        if cal_col_dens == True:
            files = sorted(glob.glob(file_dir[it] + 'column_density_profile_pickle_*.pkl'))
        else:
            files = sorted(glob.glob(file_dir[it] + field_str + '_profile_pickle_*.pkl'))
        times = []
        colors = ['k', 'b', 'c', 'g', 'r', 'm']
        dash_list =  [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
        for f_it, file in enumerate(files):
            time_val = int(file.split('_profile_pickle_')[-1].split('.pkl')[0])
            times.append(time_val)
            try:
                open_file = open(file, 'r')
                prof_x, prof_y, sampled_points, separation = pickle.load(open_file)
                open_file.close()
            except:
                open_file = open(file, 'r')
                prof_x, prof_y, sampled_points = pickle.load(open_file)
                open_file.close()
            if log_scale:
                #axes_dict[ax_label].semilogy(prof_x, prof_y, c=colors[-len(files) + f_it], dashes=dash_list[-len(files) + f_it], label=str(time_val)+"yr")
                axes_dict[ax_label].loglog(prof_x, prof_y, c=colors[-len(files) + f_it], dashes=dash_list[-len(files) + f_it], label=str(time_val)+"yr")
            else:
                axes_dict[ax_label].plot(prof_x, prof_y, c=colors[-len(files) + f_it], dashes=dash_list[-len(files) + f_it], label=str(time_val)+"yr")
            if time_val == 2000:
                if positions[it][1] == 1:
                    prof_y_2 = 150*8.e-14*(1/(prof_x**2))
                    #axes_dict[ax_label].semilogy(prof_x, prof_y_2, 'k:', alpha = 0.5)
                    axes_dict[ax_label].loglog(prof_x, prof_y_2, 'k:', alpha = 0.5)
                    axes_dict[ax_label].set_ylim([5.e-16, 8.e-14])
                if positions[it][1] == 2:
                    prof_y_2 = 150*3.2e24*(1/prof_x)
                    #axes_dict[ax_label].semilogy(prof_x, prof_y_2, 'k:', alpha = 0.5)
                    axes_dict[ax_label].loglog(prof_x, prof_y_2, 'k:', alpha = 0.5)
                    axes_dict[ax_label].set_ylim([1.e24, 8.e24])

            if time_val == 5000:
                if positions[it][1] == 3:
                    if positions[it][0] == 3:
                        rad1_ind = np.argmin(np.abs(np.array(prof_x)-70))
                        rad2_ind = np.argmin(np.abs(np.array(prof_x)-82.5))
                        mass_1 = prof_y[rad1_ind]
                        mass_2 = prof_y[rad2_ind-1]
                        print("disk mass is between", mass_1, mass_2)
                        y1 = np.zeros(np.shape(prof_y))
                        axes_dict[ax_label].fill_between(prof_x[rad1_ind:rad2_ind], np.zeros(np.shape(prof_x[rad1_ind:rad2_ind])), np.ones(np.shape(prof_x[rad1_ind:rad2_ind]))*mass_1, facecolor='grey', alpha=0.5)
                        mass_1_arr = np.ones(np.shape(prof_x[:rad2_ind]))*mass_1
                        mass_2_arr = np.ones(np.shape(prof_x[:rad2_ind]))*mass_2
                        axes_dict[ax_label].fill_between(prof_x[:rad2_ind], mass_1_arr, mass_2_arr, facecolor='grey', alpha=0.5)
                    else:
                        mass_1, mass_2 = 0.11395946798614512, 0.12237884009714563
                        mass_1_arr = np.ones(np.shape(prof_x))*mass_1
                        mass_2_arr = np.ones(np.shape(prof_x))*mass_2
                        axes_dict[ax_label].fill_between(prof_x, mass_1_arr, mass_2_arr, facecolor='grey', alpha=0.5)
            
            #axes_dict[ax_label].set_ylim([1.e24, 1.e25])
        if title != '':
            axes_dict[ax_label].set_title(title)
        if positions[it][0] == columns and args.share_colourbar == False and positions[it][1] == 1:
            axes_dict[ax_label].legend(loc='best', fontsize=args.text_font)
        if field_str == 'Tangential_Velocity':
            if positions[it][0] == 1:
                axes_dict[ax_label].set_ylabel('Velocity Dispersion (km/s)', size=args.text_font)
                axes_dict[ax_label].set_ylim(bottom=0.0)
        elif field_str == 'cell_mass':
            if positions[it][0] == 1:
                axes_dict[ax_label].set_ylabel('Cumulative Mass ($M_\odot$)', size=args.text_font)
            axes_dict[ax_label].set_ylim(bottom=0.0)
            '''
            if positions[it][0] == columns:
                axes_dict[ax_label].axvspan(65, 75, alpha=0.5, color='grey')
            '''
        elif cal_col_dens:
            if positions[it][0] == 1:
                axes_dict[ax_label].set_ylabel('H2 Column Density (cm$^{-2}$)', size=args.text_font)
        else:
            if positions[it][0] == 1:
                axes_dict[ax_label].set_ylabel('Density (g/cm$^3$)', size=args.text_font)
        '''
        if positions[it][1] == 1:
            axes_dict[ax_label].axhline(y=5.e-15, linestyle='--', color='k', alpha=0.5)
        if positions[it][1] == 2:
            axes_dict[ax_label].axhline(y=3.2e24, linestyle='--', color='k', alpha=0.5)
        '''
        #axes_dict[ax_label].set_ylim([0.0, 60.0])
        axes_dict[ax_label].set_xlim([np.min(prof_x), np.max(prof_x)])
        axes_dict[ax_label].set_xlim([0.0, 150.0])
        axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('Radius (AU)', size=args.text_font)
            '''
            if positions[it][0] != columns:
                xticklabels = axes_dict[ax_label].get_xticklabels()
                plt.setp(xticklabels[-1], visible=True)
            '''
    if 'force' in plot_type[it]:
        force_args = input_args[it].split(' ')
        arg_list = ['python', '/home/100/rlk100/Scripts/paper_plots.py', file_dir[it], save_dir, '-pd', 'True']
        for force_args in force_args:
            arg_list.append(force_args)
    
        pickle_file = file_dir[it] + 'force_comp_pickle.pkl'
        if os.path.isfile(pickle_file) == False:
            call(arg_list)
        
        file = open(pickle_file, 'r')
        x, y, times, y_label = pickle.load(file)
        file.close()
        colors = ['k', 'b', 'c', 'g', 'r', 'm']
        dash_list =  [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
        plt.axhline(y=1.0, color='k', linestyle='--')
        for time in range(len(x)):
            axes_dict[ax_label].loglog(x[time], y[time], c=colors[-len(x) + time], dashes=dash_list[-len(x) + time], label=str(times[time])+"yr")#, linewidth=1.5)
            if positions[it][0] == columns:
                plt.legend(loc='best')
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel(y_label, fontsize=args.text_font+2)
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('Z-distance (AU)', fontsize=args.text_font)
        axes_dict[ax_label].set_xlim([1.0, 1000.0])
        #data_aspect = (np.log(plt.ylim()[1])-np.log(plt.ylim()[0]))/(np.log(1000.0)-np.log(1.0))
        #axes_dict[ax_label].set_aspect(1./data_aspect)
    if 'outflow' in plot_type[it]:
        legend_labels = ['NT', 'T1', 'T2']
        linestyles = ['k-', 'b--', 'r-.']
        #csv_files = sorted(glob.glob(file_dir[it] + "mach*.csv"))
        csv_files = [file_dir[it]+"mach_0.csv", file_dir[it]+"mach_1.csv", file_dir[it]+"mach_2.csv"]
        panel_label_x_pos = 250
        if 'mass' in input_args[it]:
            time = []
            mass = []
            for file_it in range(len(csv_files)):
                header = 0
                with open(csv_files[file_it], 'r') as file:
                    time.append([])
                    mass.append([])
                    reader = csv.reader(file)
                    for row in reader:
                        if header != 0:
                            time_val = float(row[0])
                            if time_val not in time[-1]:
                                time[-1].append(time_val)
                                m = float(row[1])
                                mass[-1].append(m)
                        if header == 0:
                            header = 1
                    time_it = np.argmin(abs(np.array(time[-1])-5000))
                    mass[-1] = np.nan_to_num(mass[-1])
                    dm = np.array(mass[-1][:time_it+1])[1:] - np.array(mass[-1][:time_it+1])[:-1]
                    dt = np.array(time[-1][:time_it+1])[1:] - np.array(time[-1][:time_it+1])[:-1]
                    time_averaged_mass = np.sum(dm*dt)/np.sum(dt)
                    print("time averaged outflow mass for", csv_files[file_it], "is", time_averaged_mass)
            for t in range(len(time)):
                axes_dict[ax_label].semilogy(time[t], mass[t], linestyles[t], label=legend_labels[t])
            axes_dict[ax_label].text(panel_label_x_pos, 5.5e-2, 'a', fontsize=args.text_font+2)
            axes_dict[ax_label].set_ylabel('Outflow Mass (M$_\odot$)')
            axes_dict[ax_label].set_xlim([0, 5000])
            axes_dict[ax_label].set_ylim([1.e-4, 1.5e-1])
            #axes_dict[ax_label].set_ylim([1.e-3, 1.e-1])
            axes_dict[ax_label].legend(loc='lower right')
        elif 'max' in input_args[it]:
            time = []
            speed = []
            for file_it in range(len(csv_files)):
                header = 0
                with open(csv_files[file_it], 'r') as file:
                    time.append([])
                    speed.append([])
                    reader = csv.reader(file)
                    for row in reader:
                        if header != 0:
                            time_val = float(row[0])
                            if time_val not in time[-1]:
                                time[-1].append(time_val)
                                s = float(row[4])
                                speed[-1].append(s)
                        if header == 0:
                            header = 1
            for t in range(len(time)):
                axes_dict[ax_label].semilogy(time[t], speed[t], linestyles[t], label=legend_labels[t])
            axes_dict[ax_label].set_ylabel('Max. Speed (km/s)')
            axes_dict[ax_label].set_xlim([0, 5000])
            axes_dict[ax_label].set_ylim([1.e0, 1.e2])
            axes_dict[ax_label].set_xlabel('Time since protostar formation (yr)')
            axes_dict[ax_label].text(panel_label_x_pos, 0.5*1.e2, 'd', fontsize=args.text_font+2)
            #axes_dict[ax_label].set_ylim([1.e-3, 1.e-1])
            #axes_dict[ax_label].legend(loc='best')
        elif 'mean' in input_args[it]:
            time = []
            speed = []
            for file_it in range(len(csv_files)):
                header = 0
                with open(csv_files[file_it], 'r') as file:
                    time.append([])
                    speed.append([])
                    reader = csv.reader(file)
                    for row in reader:
                        if header != 0:
                            time_val = float(row[0])
                            if time_val not in time[-1]:
                                time[-1].append(time_val)
                                s = float(row[-1])
                                speed[-1].append(s)
                        if header == 0:
                            header = 1
            for t in range(len(time)):
                axes_dict[ax_label].semilogy(time[t], speed[t], linestyles[t], label=legend_labels[t])
            axes_dict[ax_label].set_ylabel('Mean. Speed (km/s)')
            axes_dict[ax_label].set_xlim([0, 5000])
            axes_dict[ax_label].set_ylim([1.e0, 1.e2])
            axes_dict[ax_label].set_xlabel('Time since protostar formation (yr)')
            axes_dict[ax_label].text(panel_label_x_pos, 0.5*1.e2, 'd', fontsize=args.text_font+2)
            #axes_dict[ax_label].set_ylim([1.e-3, 1.e-1])
            #axes_dict[ax_label].legend(loc='best')
        elif 'ang' in input_args[it]:
            time = []
            angular_momentum = []
            for file_it in range(len(csv_files)):
                header = 0
                with open(csv_files[file_it], 'r') as file:
                    time.append([])
                    angular_momentum.append([])
                    reader = csv.reader(file)
                    for row in reader:
                        if header != 0:
                            time_val = float(row[0])
                            if time_val not in time[-1]:
                                time[-1].append(time_val)
                                L = float(row[3])
                                if 'specific' in input_args[it]:
                                    L = L/float(row[1])
                                #if np.isnan(L):
                                #    L = 1.e-6
                                angular_momentum[-1].append(L)
                        if header == 0:
                            header = 1
                    time_it = np.argmin(abs(np.array(time[-1])-5000))
                    angular_momentum[-1] = np.nan_to_num(angular_momentum[-1])
                    dl = np.array(angular_momentum[-1][:time_it+1])[1:] - np.array(angular_momentum[-1][:time_it+1])[:-1]
                    dt = np.array(time[-1][:time_it+1])[1:] - np.array(time[-1][:time_it+1])[:-1]
                    time_averaged_angular_momentum = np.sum(dl*dt)/np.sum(dt)
                    print("time averaged outflow angular for", csv_files[file_it], "is", time_averaged_angular_momentum)
            for t in range(len(time)):
                axes_dict[ax_label].semilogy(time[t], angular_momentum[t], linestyles[t])
            if 'specific' in input_args[it]:
                axes_dict[ax_label].set_ylabel('Sp. Ang. Mom. (km$^2\,$s$^{-1}$)')
                #axes_dict[ax_label].set_ylim([1.e7, 1.e11])
            else:
                axes_dict[ax_label].set_ylabel('Ang. Mom. (M$_\odot$km$^2\,$s$^{-1}$)')
                axes_dict[ax_label].set_ylim([1.e6, 5.e9])
                axes_dict[ax_label].text(panel_label_x_pos, 1.5e9, 'c', fontsize=args.text_font+2)
                #axes_dict[ax_label].set_ylim([1.e5, 1.e10])
            #axes_dict[ax_label].set_xlabel('Time since protostar formation (yr)')
            axes_dict[ax_label].set_xlim([0, 5000])
        else:
            time = []
            momentum = []
            for file_it in range(len(csv_files)):
                header = 0
                with open(csv_files[file_it], 'r') as file:
                    time.append([])
                    momentum.append([])
                    reader = csv.reader(file)
                    for row in reader:
                        if header != 0:
                            time_val = float(row[0])
                            if time_val not in time[-1]:
                                time[-1].append(time_val)
                                p = float(row[2])
                                if 'specific' in input_args[it]:
                                    p = p/float(row[1])
                                #if np.isnan(p):
                                #    p = 1.e-6
                                momentum[-1].append(p)
                        if header == 0:
                            header = 1
                    time_it = np.argmin(abs(np.array(time[-1])-5000))
                    momentum[-1] = np.nan_to_num(momentum[-1])
                    dp = np.array(momentum[-1][:time_it+1])[1:] - np.array(momentum[-1][:time_it+1])[:-1]
                    dt = np.array(time[-1][:time_it+1])[1:] - np.array(time[-1][:time_it+1])[:-1]
                    time_averaged_linear_momentum = np.sum(dp*dt)/np.sum(dt)
                    print("time averaged outflow linear for", csv_files[file_it], "is", time_averaged_linear_momentum)
            for t in range(len(time)):
                axes_dict[ax_label].semilogy(time[t], momentum[t], linestyles[t])
            if 'specific' in input_args[it]:
                axes_dict[ax_label].set_ylabel('Specific Mom. (km$\,$s$^{-1}$)')
                #axes_dict[ax_label].set_ylim([1.e-1, 1.e1])
            else:
                axes_dict[ax_label].set_ylabel('Momentum (M$_\odot$km$\,$s$^{-1}$)')
                axes_dict[ax_label].set_ylim([1.e-4, 2.5e-1])
                axes_dict[ax_label].text(panel_label_x_pos, 8.e-2, 'b', fontsize=args.text_font+2)
                #axes_dict[ax_label].set_ylim([1.e-3, 1.e0])
            axes_dict[ax_label].set_xlim([0, 5000])
            print("CREATED OUTFLOWS PLOT")
    if 'appendix' in plot_type[it]:
        lref_labels = ['10', '11', '12', '13', '14', '15']
        colors = ['k', 'b', 'c', 'g', 'r', 'm']
        dash_list =  [[5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None), [1,1]]
        Time = []
        Mass = []
        Momentum = []
        Angular_Momentum = []
        Max_speed = []
        csv_files = sorted(glob.glob(file_dir[it] + "out*.csv"))
        for csv_file in csv_files:
            Time.append([])
            Mass.append([])
            Momentum.append([])
            Angular_Momentum.append([])
            Max_speed.append([])
            with open(csv_file, 'rU') as file:
                reader = csv.reader(file)
                header = True
                for row in reader:
                    if header == False:
                        Time[-1].append(float(row[0]))
                        mass_val = float(row[1])
                        if np.isnan(mass_val):
                            mass_val = 1.e-5
                        Mass[-1].append(mass_val)
                        mom_val = float(row[2])
                        if np.isnan(mom_val):
                            mom_val = 1.e-5
                        Momentum[-1].append(mom_val)
                        ang_mom_val = float(row[3])
                        if np.isnan(ang_mom_val):
                            ang_mom_val = 1.e5
                        Angular_Momentum[-1].append(ang_mom_val)
                        max_val = float(row[4])
                        if np.isnan(max_val):
                            max_val = 1.e-1
                        Max_speed[-1].append(max_val)
                    else:
                        header = False
        if input_args[it] == 'Mass':
            plot_field = Mass
            y_label = 'Mass (M$_\odot$)'
            y_limits = [1.e-5, 1.e0]
        elif input_args[it] == 'Momentum':
            plot_field = Momentum
            y_label = 'Momentum (M$_\odot$km$\,$s$^{-1}$)'
            y_limits = [1.e-5, 5.e0]
        elif input_args[it] == 'Angular_Momentum':
            plot_field = Angular_Momentum
            y_label = 'Angular Momentum (M$_\odot$km$^2\,$s$^{-1}$)'
            y_limits = [1.e5, 5.e10]
        else:
            plot_field = Max_speed
            y_label = 'Maximum Outflow Speed (km$\,$s$^{-1}$)'
            y_limits = [1.e-1, 1.e2]

        for sim in range(len(Time)):
            axes_dict[ax_label].semilogy(Time[sim], plot_field[sim], color=colors[sim], dashes=dash_list[sim], label='$L_\mathrm{ref}=$'+lref_labels[sim])
            print("plotted lref=", lref_labels[sim])
        axes_dict[ax_label].set_xlim([0,3000])
        axes_dict[ax_label].set_ylim(y_limits)
            
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel(y_label)
        if positions[it][1] == rows and positions[it][0] == 2:
            axes_dict[ax_label].set_xlabel('Time since protostar formation (yr)')
        if positions[it][0] == columns:
            print("plotting legend")
            axes_dict[ax_label].legend(loc='best', ncol=3) #prop={'size':16},
                
    if 'phasefolded' in plot_type[it]:
        pickle_file = file_dir[it] + 'particle_data.pkl'
        file_open = open(pickle_file, 'r')
        particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
        file_open.close()
        window = 10
        moving_index = window
        moving_average_time = []
        moving_average_sep = []
        moving_average_accretion = [[],[]]
        while moving_index < len(particle_data['time']):
            moving_average_time.append(np.mean(particle_data['time'][moving_index-window: moving_index]))
            moving_average_sep.append(np.mean(particle_data['separation'][moving_index-window: moving_index]))
            moving_average_accretion[0].append(np.mean(particle_data['mdot'][0][moving_index-window: moving_index]))
            moving_average_accretion[1].append(np.mean(particle_data['mdot'][1][moving_index-window: moving_index]))
            moving_index = moving_index + 1
        particle_data['time'] = np.array(moving_average_time)
        particle_data['separation'] = np.array(moving_average_sep)
        particle_data['mdot'] = np.array(moving_average_accretion)
        periastron_inds = [np.argmin(particle_data['separation'])]
        apastron_inds = []
        index_interval = 500
        index = periastron_inds[0] + index_interval
        passed_apastron = False
        while index < len(particle_data['time']):
            if np.min(particle_data['separation'][index-index_interval:index]) != np.min(particle_data['separation'][np.array([index-index_interval, index-index_interval+1, index-2, index-1])]) and passed_apastron == True:
                periastron_ind = np.argmin(particle_data['separation'][index-index_interval:index]) + index-index_interval
                periastron_inds.append(periastron_ind)
                print("found periastron at time", particle_data['time'][periastron_ind], "of separation", particle_data['separation'][periastron_ind])
                passed_apastron = False
            elif np.max(particle_data['separation'][index-index_interval:index]) != np.max(particle_data['separation'][np.array([index-index_interval, index-index_interval+1, index-2, index-1])]) and passed_apastron == False:
                apastron_ind = np.argmax(particle_data['separation'][index-index_interval:index]) + index-index_interval
                apastron_inds.append(apastron_ind)
                print("found apastron at time", particle_data['time'][apastron_ind], "of separation", particle_data['separation'][apastron_ind])
                passed_apastron = True
            index = index + index_interval/2
    
        no_of_folded_orbits = 15
        plot_start_ind = int(input_args[it])
        plot_end_ind = plot_start_ind + no_of_folded_orbits
        '''
        plt.plot(particle_data['time'][periastron_inds[plot_start_ind]:periastron_inds[plot_end_ind]], particle_data['separation'][periastron_inds[plot_start_ind]:periastron_inds[plot_end_ind]], 'k-')
        plt.xlabel('time (yr)')
        plt.ylabel('separation (au)')
        for periastron in periastron_inds[plot_start_ind:plot_end_ind]:
            plt.axvline(x=particle_data['time'][periastron], alpha=0.5)
        for apastron in apastron_inds[plot_start_ind:plot_end_ind-1]:
            plt.axvline(x=particle_data['time'][apastron], color='r', alpha=0.5)
        plt.savefig(save_dir + 'separation_start_orbit_'+ str(plot_start_ind) +'.eps')
        plt.clf()
        print "Created separation evolution plot"
        '''

        hist_ind = 1
        ap_ind = 0
        phase = np.linspace(0, 1, 20)
        #FIND BINNED DATA
        averaged_binned_accretion = [[],[]]
        averaged_total_accretion = []
        while hist_ind < len(periastron_inds[plot_start_ind:plot_end_ind]):
            binned_accretion = [[],[]]
            total_accretion = []
            time_bins_1 = np.linspace(particle_data['time'][periastron_inds[hist_ind-1]], particle_data['time'][apastron_inds[ap_ind]],11)
            time_bins_2 = np.linspace(particle_data['time'][apastron_inds[ap_ind]], particle_data['time'][periastron_inds[hist_ind]],11)
            bin_ind = 1
            while bin_ind < len(time_bins_1):
                time_bin_inds_1 = np.where((particle_data['time'] > time_bins_1[bin_ind-1]) & (particle_data['time'] < time_bins_1[bin_ind]))[0]
                binned_accretion[0].append(np.mean(particle_data['mdot'][0][time_bin_inds_1]))
                binned_accretion[1].append(np.mean(particle_data['mdot'][1][time_bin_inds_1]))
                total_accretion.append(np.mean(particle_data['mdot'][:,time_bin_inds_1]))
                bin_ind = bin_ind + 1
            bin_ind = 1
            while bin_ind < len(time_bins_2):
                time_bin_inds_2 = np.where((particle_data['time'] > time_bins_2[bin_ind-1]) & (particle_data['time'] < time_bins_2[bin_ind]))[0]
                binned_accretion[0].append(np.mean(particle_data['mdot'][0][time_bin_inds_2]))
                binned_accretion[1].append(np.mean(particle_data['mdot'][1][time_bin_inds_2]))
                total_accretion.append(np.mean(particle_data['mdot'][:,time_bin_inds_2]))
                bin_ind = bin_ind + 1
            hist_ind = hist_ind + 1
            ap_ind = ap_ind + 1
            averaged_binned_accretion[0].append(binned_accretion[0])
            averaged_binned_accretion[1].append(binned_accretion[1])
            averaged_total_accretion.append(total_accretion)
        phase_2 = np.linspace(1, 2, 20)
        phase_2 = phase.tolist() + phase_2[1:].tolist()
        median_accretion = []
        standard_deviation = []
        median_accretion.append(np.median(averaged_binned_accretion[0], axis=0))
        median_accretion.append(np.median(averaged_binned_accretion[1], axis=0))
        standard_deviation.append(np.std(averaged_binned_accretion[0], axis=0))
        standard_deviation.append(np.std(averaged_binned_accretion[1], axis=0))
        median_total = np.median(averaged_total_accretion, axis=0)
        standard_deviation_total = np.std(averaged_total_accretion, axis=0)
        #plt.plot(phase_2, median_accretion[0].tolist()+median_accretion[0][1:].tolist(), ls='steps', alpha=0.5, label='Primary')
        #plt.plot(phase_2, median_accretion[1].tolist()+median_accretion[1][1:].tolist(), ls='steps', alpha=0.5, label='Secondary')
        #plt.plot(phase_2, median_total.tolist() + median_total[1:].tolist(),ls='steps', label='Total')
        axes_dict[ax_label].errorbar(phase_2, median_accretion[0].tolist()+median_accretion[0][1:].tolist(), yerr=standard_deviation[0].tolist()+standard_deviation[0][1:].tolist(), ls='steps-mid', alpha=0.5, label='Primary')
        axes_dict[ax_label].errorbar(phase_2, median_accretion[1].tolist()+median_accretion[1][1:].tolist(), yerr=standard_deviation[1].tolist()+standard_deviation[1][1:].tolist(), ls='steps-mid', alpha=0.5, label='Secondary')
        axes_dict[ax_label].errorbar(phase_2, median_total.tolist() + median_total[1:].tolist(), yerr=standard_deviation_total.tolist()+standard_deviation_total[1:].tolist(), ls='steps-mid', label='Total')
        
        #axes_dict[ax_label].legend(loc='best')
        #plt.xlabel("Orbital Phase")
        axes_dict[ax_label].set_ylabel("Normalised accretion")
        #plt.title("Start orbit:" + str(start_ind))
        axes_dict[ax_label].set_xlim([0.0, 1.3])
        axes_dict[ax_label].set_ylim(bottom=0.0)

    if positions[it][0] != 1:
        yticklabels = axes_dict[ax_label].get_yticklabels()
        plt.setp(yticklabels, visible=False)
        yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
        plt.setp(yticklabels, visible=False)
    if positions[it][1] != rows:
        xticklabels = axes_dict[ax_label].get_xticklabels()
        plt.setp(xticklabels, visible=False)
    if positions[it][0] == 1:
        axes_dict[ax_label].tick_params(axis='y', which='major', labelsize=args.text_font)
    if positions[it][1] == rows:
        axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
        if positions[it][0] != 1:
            xticklabels = axes_dict[ax_label].get_xticklabels()
            if 'force' in plot_type[it]:
                plt.setp(xticklabels[1], visible=False)
            else:
                plt.setp(xticklabels[0], visible=False)
                if positions[it][0] == 3:
                    plt.setp(xticklabels[1], visible=False)

    if 'multi' in plot_type[it]:
        if positions[it][1] == rows and positions[it][0] == columns:
            yticklabels = axes_dict['ax13'].get_yticklabels()
            plt.setp(yticklabels[-1], visible=False)
    if args.multiplot_title != "":
        plt.suptitle(args.multiplot_title, x=0.45, y=0.91, fontsize=18)

    #f.savefig(savename + '.pdf', format='pdf')
    #f.savefig(savename + '.eps', format='eps')
    f.savefig(savename + '.pdf', format='pdf', bbox_inches='tight')
    f.savefig(savename + '.eps', format='eps', bbox_inches='tight')
