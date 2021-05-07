import matplotlib.pyplot as plt
import pickle
import csv
import numpy as np
import os
from matplotlib.colors import LogNorm
import h5py
import my_ramses_module as mym
import matplotlib.gridspec as gridspec
import glob
import matplotlib.patheffects as path_effects
import matplotlib
import matplotlib.colors as colors
import yt
from scipy.ndimage import zoom
from matplotlib import ticker
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
    parser.add_argument("-plot_cbar", "--plot_colourbar", help="Do you need to plot a colorbar?", default=True)
    parser.add_argument("-tf", "--text_font", help="what font do you want the text to have?", type=int, default=12)
    parser.add_argument("-sn", "--savename", help="what do you want to save the plot as?", type=str, default='multiplot')
    parser.add_argument("-ylp", "--y_label_pad", help="y lable pad", default=-20, type=float)
    parser.add_argument("-title", "--multiplot_title", help="Do you want to title the multiplot", type=str, default="")
    args = parser.parse_args()
    return args


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
if args.plot_colourbar == 'False':
    args.plot_colourbar = False
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
            try:
                height_ratios = []
                for ratio in row[5].split(':'):
                    height_ratios.append(float(ratio))
            except:
                continue
        elif row[0][0]!= '#':
            positions.append((int(row[0]), int(row[1])))
            plot_type.append(row[2])
            file_dir.append(row[3])
            input_args.append([])
            for input_arg in row[4].split(' '):
                if len(input_arg) > 1:
                    input_args[-1].append(input_arg)
positions = np.array(positions)

if 0 in positions or args.plot_colourbar == False:
    make_right_gs = False
else:
    make_right_gs = True

columns = np.max(positions[:,0])
rows = np.max(positions[:,1])

width = float(columns)*(14.5/3.)
if args.plot_colourbar == True:
    height = float(rows)*(17./4.)
else:
    height = float(rows)*(13./4.)

f = plt.figure(figsize=(width, height))

if make_right_gs:
    gs_left = gridspec.GridSpec(rows, columns-1)
    gs_right = gridspec.GridSpec(rows, 1)

    gs_left.update(right=glr, wspace=glw, hspace=ghspace)
    gs_right.update(left=grl, hspace=ghspace)
else:
    try:
        gs_left = gridspec.GridSpec(rows, columns, height_ratios=height_ratios)
    except:
        gs_left = gridspec.GridSpec(rows, columns)
    gs_left.update(right=glr, wspace=glw, hspace=ghspace)

axes_dict = {}

for it in range(len(positions)):
    #ax_label = 'ax' + str(counter)
    ax_label = 'ax' + str(positions[it][0])+str(positions[it][1])
    if positions[it][0] == 0:
        axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,:])})
    elif positions[it][1] == 0:
        axes_dict.update({ax_label:f.add_subplot(gs_left[:,positions[it][0]-1])})
    elif positions[it][0] == 1 and positions[it][1] == 1:
        if columns > 1:
            axes_dict.update({ax_label:f.add_subplot(gs_left[0,0])})
        else:
            if make_right_gs:
                axes_dict.update({ax_label:f.add_subplot(gs_right[0,0])})
            else:
                axes_dict.update({ax_label:f.add_subplot(gs_left[0,0])})
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
            if make_right_gs:
                axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'], sharey=axes_dict['ax11'])})
            else:
                #spoofing this line. Not sure what to use.
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'], sharey=axes_dict['ax11'])})
        elif args.share_x and args.share_y == 'row':
            if positions[it][0] == 1:
                if make_right_gs:
                    axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'])})
                else:
                    #spoofing this line. Not sure what to use.
                    axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'])})
            else:
                if make_right_gs:
                    axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'], sharey=axes_dict['ax1'+str(positions[it][1])])})
                else:
                    #spoofing this line. Not sure what to use.
                    axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'], sharey=axes_dict['ax1'+str(positions[it][1])])})
        elif args.share_x:
            if make_right_gs:
                axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0], sharex=axes_dict['ax11'])})
            else:
                #spoofing this line. Not sure what to use.
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1], sharex=axes_dict['ax11'])})
        else:
            if make_right_gs:
                axes_dict.update({ax_label:f.add_subplot(gs_right[positions[it][1]-1,0])})
            else:
                #spoofing this line. Not sure what to use.
                axes_dict.update({ax_label:f.add_subplot(gs_left[positions[it][1]-1,positions[it][0]-1])})

    if plot_type[it] == 'RV_proj':
        file = open(file_dir[it], 'rb')
        X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
        file.close()
        
        if int(np.max(part_info['particle_position'][1])) == 200:
            image = np.rot90(np.rot90(image))
            vel_rad = np.rot90(np.rot90(vel_rad))
            velx = np.rot90(np.rot90(velx))
            vely = np.rot90(np.rot90(vely))
            part_info['particle_position'][1][0] = part_info['particle_position'][1][0]*-1
        
        if '-cmin' in input_args[it]:
            cmin_ind = input_args[it].index('-cmin')
            cbar_min = eval(input_args[it][cmin_ind+1])
        else:
            cbar_min = None
        if '-cmax' in input_args[it]:
            cmax_ind = input_args[it].index('-cmax')
            cbar_max = eval(input_args[it][cmax_ind+1])
        else:
            cbar_max = None
        
        if '-threshold' in input_args[it]:
            thres_ind = input_args[it].index('-threshold')
            density_threshold = eval(input_args[it][thres_ind+1])
        else:
            density_threshold = 0.0
            
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
        axes_dict[ax_label].set_xlim(xlim)
        axes_dict[ax_label].set_ylim(ylim)
    
        bool_den_array = image>density_threshold
        vel_rad = vel_rad*bool_den_array #bool_den_array*np.nan*vel_rad
        vel_rad[vel_rad == 0] = np.nan
        
        v_std = np.std(vel_rad/10000)
        v_cbar_min = -10 #center_vel_rv.in_units('km/s').value - 5
        v_cbar_max = 10
        
        plot = axes_dict[ax_label].pcolormesh(X, Y, vel_rad/10000, cmap='idl06_r', rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
        fmt = ticker.LogFormatterSciNotation()
        fmt.create_dummy_axis()
            
        if density_threshold != 0.0:
            exp_min = np.log10(density_threshold)
        else:
            exp_min = np.log10(cbar_min)
        exp_max = np.log10(cbar_max)
        n_level = (exp_max-exp_min)*2 + 1
        contour_levels = np.logspace(exp_min, exp_max, int(n_level))
        CS = axes_dict[ax_label].contour(X,Y,image, locator=plt.LogLocator(), linewidths=0.5, colors='k', levels=contour_levels)
        def func(x):
            s = "%.0g" % x
            if "e" in s:
                tup = s.split('e')
                significand = tup[0].rstrip('0').rstrip('.')
                sign = tup[1][0].replace('+', '')
                exponent = tup[1][1:].lstrip('0')
                s = ('%se%s%s' % (significand, sign, exponent)).rstrip('e')
            return s
        axes_dict[ax_label].clabel(CS,CS.levels,fmt=func)
        
        if has_particles:
            mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
            
        if positions[it][0] == columns and positions[it][1] == rows:
            cax = f.add_axes([0.8, 0.11, 0.02, 0.77])
            cbar =f.colorbar(plot, pad=0.0, cax=cax)#, label='Radial Velocity (km/s)')
            cbar.set_label('Radial Velocity (km/s)', rotation=270, labelpad=7, size=args.text_font)
            cbar.ax.tick_params(labelsize=args.text_font)
            #cax.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=args.text_font)
        axes_dict[ax_label].set_aspect('equal')
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('$x$ (AU)', fontsize=args.text_font, labelpad=-2)
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel('$y$ (AU)', fontsize=args.text_font, labelpad=-7)
    elif plot_type[it] == 'CF_hist':
        file = open("Tobin_CF.pkl", 'rb')
        S_bins, CF_per_bin_Tobin = pickle.load(file)
        file.close()
        
        file = open(file_dir[it], 'rb')
        Times, CF_Array_Full, N_sys_total = pickle.load(file)
        file.close()
        
        S_bins = np.logspace(0.75,4,14)
        bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
        CF_median = []
        CF_err = []
        
        for bit in range(len(S_bins)-1):
            median = np.median(CF_Array_Full[:,bit])
            mean = np.mean(CF_Array_Full[:,bit])
            std = np.std(CF_Array_Full[:,bit], ddof=1)
            
            standard_deviation = [median-(mean-std), (mean+std)-median]
            CF_median.append(median)
            CF_err.append(standard_deviation)

        CF_err = np.array(CF_err)

        try:
            axes_dict[ax_label].bar(bin_centers, CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
        except:
            axes_dict[ax_label].bar(bin_centers[1:], CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
        axes_dict[ax_label].bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
        if positions[it][0] == 1 and positions[it][1] == 1:
            axes_dict[ax_label].legend(loc='best')
        if positions[it][1] == rows:
            axes_dict[ax_label].xlabel('Log Separation (AU)')
        axes_dict[ax_label].ylabel('Companion Frequency')
        axes_dict[ax_label].xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
        axes_dict[ax_label].ylim(bottom=0.0)
    elif plot_type[it] == 'Comp_hist':
        file = open("Tobin_CF.pkl", 'rb')
        S_bins, CF_per_bin_Tobin = pickle.load(file)
        file.close()
        CF_per_bin_Tobin = [0.0, 0.0] + CF_per_bin_Tobin
        
        if '-bound_comp' in input_args[it]:
            bound_comp_ind = input_args[it].index('-bound_comp')
            bound_comp = eval(input_args[it][bound_comp_ind+1])
        else:
            bound_comp = False
        if '-plot_Tobin' in input_args[it]:
            plot_Tobin_ind = input_args[it].index('-plot_Tobin')
            plot_Tobin = eval(input_args[it][plot_Tobin_ind+1])
        else:
            plot_Tobin = True
        
        dir_split = file_dir[it].split('*')
        if bound_comp == False:
            label_string = ['G400', 'G200', 'G100', 'G50']
        else:
            label_string = ['Unbound', 'Bound']
            subplot_title = file_dir[it].split('Global/')[1].split('/')[0]
        
        pickle_files = []
        for label in label_string:
            if label == 'G100':
                pickle_files.append(dir_split[0] + label + '/256' + dir_split[1])
            else:
                pickle_files.append(dir_split[0] + label + dir_split[1])
        hist_colors = ['b', 'orange', 'g', 'magenta']
        if dir_split[1].split('/')[1] == 'SFE':
            subplot_title = 'SFE'
        elif dir_split[1].split('/')[1] == 'SFE_t_ff':
            subplot_title = 'SFE$_n$'
        elif dir_split[1].split('/')[1] == 'Number_Stars':
            subplot_title = 'N$_{stars}$'
        elif dir_split[1].split('/')[1] == 'M_tot_150':
            subplot_title = 'M$_{Tot}$=150M$_\odot$'
        
        
        S_bins = np.logspace(0.75,4,14)
        bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
        
        if plot_Tobin == True:
            axes_dict[ax_label].bar(bin_centers, CF_per_bin_Tobin, edgecolor='k', label="Tobin", width=0.25, fill=False)
        #axes_dict[ax_label].text(2, 0.3, "")
                
        max_y_lim = 0
        elinewidths = [1.0, 0.5]
        N_instances = []
        
        for pit in range(len(pickle_files)):
            try:
                file = open(pickle_files[pit], 'rb')
                Times, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
                file.close()
                N_instances.append(np.shape(CF_Array_Full)[0])
        
                CF_median = []
                CF_err = []
                
                for bit in range(len(S_bins)-1):
                    median = np.median(CF_Array_Full[:,bit])
                    mean = np.mean(CF_Array_Full[:,bit])
                    std = np.std(CF_Array_Full[:,bit], ddof=1)
                    standard_deviation = [median-(mean-std), (mean+std)-median]
                    CF_median.append(median)
                    CF_err.append(standard_deviation)
                
                CF_err = np.array(CF_err)
                CF_median = np.array(CF_median)
                axes_dict[ax_label].bar(bin_centers, CF_median, edgecolor=hist_colors[pit], label=label_string[pit], width=0.25, color=hist_colors[pit], ecolor=hist_colors[pit], fill=False, alpha=0.5)#, yerr=CF_err.T, error_kw={'elinewidth':elinewidths[pit]})
                if np.max(CF_median) > max_y_lim:
                    max_y_lim = np.max(CF_median)
            except:
                print("Skipping", pickle_files[pit], "Because it doesn't exist")
        
        if positions[it][0] == 1 and positions[it][1] == 1:
            axes_dict[ax_label].legend(loc='best')
        if positions[it][0] == 1:
            subplot_title = subplot_title + ':' +str(N_instances[0])+ ' instances over '+str(Times[-1]-Times[0])+'yr'
            axes_dict[ax_label].text(0.75, 0.15, subplot_title)
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('Log Separation (AU)')
            max_y_lim = np.ceil(max_y_lim*10)/10
            axes_dict[ax_label].set_ylim([0.0, max_y_lim])
            axes_dict[ax_label].set_xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
        if bound_comp == False:
            if positions[it][0] == 1:
                if positions[it][1] == 1:
                    axes_dict[ax_label].set_ylabel('CF ($L$ limit)')
                else:
                    axes_dict[ax_label].set_ylabel('CF (No $L$ limit)')
        elif bound_comp == True:
            axes_dict[ax_label].set_ylabel('CF')
        #axes_dict[ax_label].set_ylim(bottom=0.0)
        #axes_dict[ax_label].set_xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
    
    ax_r = axes_dict[ax_label].secondary_yaxis('right')
    ax_t = axes_dict[ax_label].secondary_xaxis('top')
    xticklabels = ax_t.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    yticklabels = ax_r.get_yticklabels()
    plt.setp(yticklabels, visible=False)
    ax_r.tick_params(axis='y', direction='inout')
    ax_t.tick_params(axis='x', direction='inout')
    
    if positions[it][0] != 1:
        yticklabels = axes_dict[ax_label].get_yticklabels()
        plt.setp(yticklabels, visible=False)
        yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
        plt.setp(yticklabels, visible=False)
    if positions[it][0] == 1:
        axes_dict[ax_label].tick_params(axis='y', which='major', labelsize=args.text_font)
        if positions[it][1] > 1:
            yticklabels = axes_dict[ax_label].get_yticklabels()
            plt.setp(yticklabels[-1], visible=False)
    if positions[it][1] != rows:
        xticklabels = axes_dict[ax_label].get_xticklabels()
        plt.setp(xticklabels, visible=False)
    if positions[it][1] == rows:
        axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
        if positions[it][0] != 1:
            xticklabels = axes_dict[ax_label].get_xticklabels()
            plt.setp(xticklabels[0], visible=False)
    
    axes_dict[ax_label].tick_params(direction='inout')

    if args.multiplot_title != "":
        plt.suptitle(args.multiplot_title, x=0.45, y=0.9, fontsize=18)

    #f.savefig(savename + '.pdf', format='pdf')
    #f.savefig(savename + '.eps', format='eps')
    f.savefig(savename + '.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)
    #f.savefig(savename + '.eps', format='eps', bbox_inches='tight', pad_inches = 0.02)
