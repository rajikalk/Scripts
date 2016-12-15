import matplotlib.pyplot as plt
import pickle
import csv
import numpy as np
import os
from subprocess import call
from matplotlib.colors import LogNorm
import h5py
import my_module as mym
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.gridspec as gridspec

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-col", "--columns", help="Number of colums?", type=int)
    parser.add_argument("-row", "--rows", help="Number of rows?", type=int)
    parser.add_argument("-in", "--in_file", help="input file with images, positions and crops")
    parser.add_argument("-sx", "--share_x", help="do you want to share the x axis?", default=False)
    parser.add_argument("-sy", "--share_y", help="do you want to share the y axis?", default=True)
    parser.add_argument("-sa", "--share_ax", help="do you want to share axes?", default=True)
    parser.add_argument("-sd", "--save_dir", help="Where do you want to save the plot? defaults to same as input file directory")
    parser.add_argument("-sc", "--share_colourbar", help="Do you want to share colour bars over rows?", default=False)
    parser.add_argument("-tf", "--text_font", help="what font do you want the text to have?", type=int, default=12)
    parser.add_argument("-sn", "--savename", help="what do you want to save the plot as?", type=str, default='multiplot')
    args = parser.parse_args()
    return args

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

#=============================MAIN==============================
args = parse_inputs()
if args.share_x == 'True':
    args.share_x = True
if args.share_y == 'False':
    args.share_y = False
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

with open(args.in_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0][0]!= '#':
            positions.append((int(row[0]), int(row[1])))
            plot_type.append(row[2])
            file_dir.append(row[3])
            input_args.append(row[4])
'''
f = plt.figure(figsize=(14.5, 17))

if args.share_colourbar:
    grid = ImageGrid(f, 111,          # as in plt.subplot(111)
                 nrows_ncols=(args.rows,args.columns),
                 axes_pad=0.02,
                 share_all=args.share_ax,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="5%",
                 cbar_pad=0.0,
                 )
else:
    grid = ImageGrid(f, 111,          # as in plt.subplot(111)
                     nrows_ncols=(args.rows,args.columns),
                     axes_pad=0.02,
                     share_all=args.share_ax,
                     cbar_location="right",
                     cbar_mode="edge",
                     cbar_size="5%",
                     cbar_pad=0.0,
                     )


if args.share_x and args.share_y:
    f, axarr = plt.subplots(args.rows, args.columns, sharex='col', sharey='row', figsize=(14.5, 17))
elif args.share_x:
    f, axarr = plt.subplots(args.rows, args.columns, sharex='col', figsize=(14.5, 17))
elif args.share_y:
    f, axarr = plt.subplots(args.rows, args.columns, sharey='row', figsize=(14.5, 17))
else:
    f, axarr = plt.subplots(args.rows, args.columns, figsize=(14.5, 17))
for ax1 in range(len(axarr)):
    for ax2 in range(len(axarr[ax1])):
        ax = axarr[ax1][ax2]
        ax.set(adjustable='box-forced', aspect='equal')
plt.subplots_adjust(wspace=-0.3, hspace=0.01)
'''
f = plt.figure(figsize=(14.5, 17))
gs_left = gridspec.GridSpec(args.rows, args.columns-1)
gs_right = gridspec.GridSpec(args.rows, 1)
gs_left.update(right=0.60)
gs_right.update(left=0.3475)
gs_left.update(wspace=-0.08)
gs_left.update(hspace=0.05)
gs_right.update(hspace=0.05)

cbar_plotted = False

axes_dict = {}
counter = 1

for it in range(len(positions)):
    #ax = grid[it]
    #ax = axarr[positions[it][1]-1][positions[it][0]-1]
    ax_label = 'ax' + str(counter)
    counter = counter + 1
    if positions[it][0] != args.columns:
        if positions[it] == (1,1):
            axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])})
            #ax = f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])
            #ax1 = ax
        else:
            if args.share_x and args.share_y:
                #ax = f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=ax1, sharey=ax1)
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax1'], sharey=axes_dict['ax1'])})
            elif args.share_x:
                #ax = f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=ax1)
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax1'])})
            elif args.share_y:
                #ax = f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharey=ax1)
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharey=axes_dict['ax1'])})
            else:
                #ax = f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])})
    else:
        if args.share_x and args.share_y:
            #ax = f.add_subplot(gs_right[positions[it][1]-1,0], sharex=ax1, sharey=ax1)
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax1'], sharey=axes_dict['ax1'])})
        elif args.share_x:
            #ax = f.add_subplot(gs_right[positions[it][1]-1,0], sharex=ax1)
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax1'])})
        elif args.share_y:
            #ax = f.add_subplot(gs_right[positions[it][1]-1,0], sharey=ax1)
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharey=axes_dict['ax1'])})
        else:
            #ax = f.add_subplot(gs_right[positions[it][1]-1,0])
            axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0])})
    axes_dict[ax_label].set(adjustable='box-forced', aspect='equal')
    #ax..set(adjustable='box-forced', aspect='equal')
    if 'movie' in plot_type[it]:
        mov_args = input_args[it].split(' ')
        arg_list = ['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/movie_script_mod.py', file_dir[it], save_dir, '-pd', 'True', '-tf', str(args.text_font)]
        concat_string = ''
        adding = False
        for mov_arg in mov_args:
            arg_list.append(mov_arg)
        call(arg_list)
        pickle_file = save_dir + 'movie_pickle.pkl'
        file = open(pickle_file, 'r')
        movie_file, X, Y, X_vel, Y_vel, image, velx, vely, part_info, args_dict, simfo, margs = pickle.load(file)
        file.close()
        movf = h5py.File(movie_file, 'r')
        #plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max']), rasterized=True)
        plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args_dict['cbar_min'], vmax=args_dict['cbar_max']), rasterized=True)
        #ax.set_aspect('equal')
        axes_dict[ax_label].set_aspect('equal')
        magx = get_image_arrays(movf, 'mag'+margs.axis[0]+'_'+simfo['movie_file_type']+'_'+margs.axis, simfo, margs)
        magy = get_image_arrays(movf, 'mag'+margs.axis[1]+'_'+simfo['movie_file_type']+'_'+margs.axis, simfo, margs)
        #ax.streamplot(X, Y, magx, magy, density=3, linewidth=0.5, minlength=0.5, arrowstyle='-')
        axes_dict[ax_label].streamplot(X, Y, magx, magy, density=3, linewidth=0.5, minlength=0.5, arrowstyle='-')
        #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args_dict['annotate_velocity'])
        mym.my_own_quiver_function(axes_dict[ax_label], X_vel, Y_vel, velx, vely, plot_velocity_legend=args_dict['annotate_velocity'])
        #mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], np.max(X), annotate_field=part_info['particle_mass'])
        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], np.max(X), annotate_field=part_info['particle_mass'])
        if 'annotate_time' in args_dict.keys():
            #ax.annotate(args_dict['annotate_time'], xy=(0.98*simfo['xmin'], 0.93*simfo['ymax']), va="center", ha="left", color='w', fontsize=args.text_font)
            axes_dict[ax_label].annotate(args_dict['annotate_time'], xy=(0.98*simfo['xmin'], 0.93*simfo['ymax']), va="center", ha="left", color='w', fontsize=args.text_font)
        #ax.annotate(args_dict['title'], xy=(0.0, 0.93*simfo['ymax']), va="center", ha="center", color='w', fontsize=(args.text_font+2))
        axes_dict[ax_label].annotate(args_dict['title'], xy=(0.0, 0.93*simfo['ymax']), va="center", ha="center", color='w', fontsize=(args.text_font+2))
        #ax.set_xlim([np.min(X), np.max(X)])
        axes_dict[ax_label].set_xlim([np.min(X), np.max(X)])
        #ax.set_ylim([np.min(Y), np.max(Y)])
        axes_dict[ax_label].set_ylim([np.min(Y), np.max(Y)])
        if positions[it][0] == 1:
            #ax.tick_params(axis='y', which='major', labelsize=args.text_font)
            axes_dict[ax_label].tick_params(axis='y', which='major', labelsize=args.text_font)
            #ax.set_ylabel(args_dict['yabel'], labelpad=-20, fontsize=args.text_font)
            axes_dict[ax_label].set_ylabel(args_dict['yabel'], labelpad=-20, fontsize=args.text_font)
        if args.share_colourbar == False:
            if positions[it][0] == args.columns:
                #cbar = plt.colorbar(plot, pad=0.0, ax=ax)
                cbar = plt.colorbar(plot, pad=0.0, ax=axes_dict[ax_label])
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                '''
                ax.cax.colorbar(plot)
                ax.cax.toggle_label(True)
                ax.cax.set_label('Density (gcm$^{-3}$)')
                grid[it].cax.colorbar(plot)
                cax = grid.cbar_axes[it]
                axis = cax.axis[cax.orientation]
                axis.label.set_text('Density (gcm$^{-3}$)')
                '''
        else:
            if cbar_plotted == False:
                '''
                for cax in grid.cbar_axes:
                    cax.toggle_label(False)
                grid[0].cax.colorbar(plot)
                cax = grid.cbar_axes[0]
                cax.toggle_label(True)
                axis = cax.axis[cax.orientation]
                axis.label.set_text('Density (gcm$^{-3}$)')
                axis.label.set_fontsize(args.text_font)
                '''
                #cbar = grid.cbar_axes[0].colorbar(plot)
                #cbar.set_label_text('Density (gcm$^{-3}$)', fontdict={'size':args.text_font})
                cax = f.add_axes([0.9, 0.1, 0.02, 0.8])
                f.colorbar(plot, pad=0.0, cax=cax)
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
                cbar_plotted = True
        if positions[it][1] == args.rows:
            '''
            if positions[it][0] == args.columns:
                xticks = ax.xaxis.get_major_ticks()
                xticks[1].label1.set_visible(False)
            '''
            #ax.set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
            axes_dict[ax_label].set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
            #ax.tick_params(axis='x', which='major', labelsize=args.text_font)
            axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
        if positions[it][0] != 1:
            #axes_dict[ax_label].set_yticklabels([])
            yticklabels = axes_dict[ax_label].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        if positions[it][1] != args.rows:
            xticklabels = axes_dict[ax_label].get_xticklabels()
            plt.setp(xticklabels, visible=False)
    if 'slice' in plot_type[it]:
        slice_args = input_args[it].split(' ')
        arg_list = ['python', '/home/100/rlk100/Scripts/paper_plots.py', file_dir[it], save_dir, '-pd', 'True']
        for slice_arg in slice_args:
            arg_list.append(slice_arg)
        call(arg_list)
        pickle_file = save_dir + 'slice_pickle.pkl'
        file = open(pickle_file, 'r')
        X, Y, X_vel, Y_vel, field_grid, velx, vely, part_info = pickle.load(file)
        file.close()
        plot = ax.pcolormesh(X, Y, field_grid, cmap=cmap, rasterized=True)
        mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend)
        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], np.max(X))
        ax.set_xlim([np.min(X), np.max(X)])
        ax.set_ylim([np.min(Y), np.max(Y)])
        if positions[it][0] == args.columns and args.share_colourbar == False:
            cbar = plt.colorbar(plot, pad=0.0)
            cbar.set_label("Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)", rotation=270, labelpad=14, size=args.text_font)
            plot.set_clim([0.0, 2.0])
        if positions[it][0] != 1:
            ax.set_ylabel('$y$ (AU)', labelpad=-20, fontsize=args.text_font)
        if positions[it][1] != args.rows:
            ax.set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
    if 'profile' in plot_type[it]:
        prof_args = input_args[it].split(' ')
        arg_list = ['python', '/home/100/rlk100/Scripts/paper_plots.py', file_dir[it], save_dir, '-pd', 'True']
        for prof_arg in prof_args:
            arg_list.append(prof_arg)
        call(arg_list)
        pickle_file = save_dir + 'slice_pickle.pkl'
        file = open(pickle_file, 'r')
        prof_x, prof_y, sampled_points = pickle.load(file)
        file.close()
        cm = plt.cm.get_cmap('RdYlBu')
        plot = ax.scatter(sampled_points[1], sampled_points[2], c=sampled_points[0], alpha=0.4, cmap=cm)
        ax.plot(prof_x, prof_y[args.y_field], 'k-', linewidth=2.)
        ax.set_xlim([0.0, args.r_max])
        ax.set_ylim([0.0, 2.0])
        ax.axhline(y=1.0, color='k', linestyle='--')
        if positions[it][0] == args.columns and args.share_colourbar == False:
            cbar = plt.colorbar(plot, pad=0.0)
            cbar.set_label('|z position| (AU)', rotation=270, labelpad=13, size=args.text_font)
        if positions[it][0] == 1:
            ax.set_ylabel('Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)', labelpad=-3)
        if positions[it][1] == args.rows:
            ax.set_xlabel('Cyclindral Radius (AU)', labelpad=-1)
    f.savefig(savename + '.pdf', format='pdf', bbox_inches='tight', pad_inches=0)
    f.savefig(savename + '.eps', format='eps', bbox_inches='tight', pad_inches=0)


