import matplotlib.pyplot as plt
import pickle
import csv
import numpy as np
import os
from matplotlib.colors import LogNorm
import my_ramses_module as mym
import matplotlib.gridspec as gridspec
import glob
import matplotlib.patheffects as path_effects
import matplotlib
import matplotlib.colors as colors
import yt
from scipy.ndimage import zoom
from matplotlib import ticker
from scipy.ndimage import gaussian_filter
import matplotlib.patches as patches

#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['mathtext.fontset'] = 'custom'
#matplotlib.rcParams['mathtext.fontset'] = 'sans-serif'
matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

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
    parser.add_argument("-conv", "--convolve", help="Do you want to convolve the image with a gaussian beam?", type=str, default='True')
    parser.add_argument("-fig_w", "--figure_width", help="Do you want to set the image width", type=float, default=None)
    parser.add_argument("-fig_h", "--figure_height", help="Do you want to set the image height", type=float, default=None)
    parser.add_argument("-arrow_width", "--arrow_width_scale", help="How do you want to scale the width of the quiver arrows?", type=float, default=1)
    parser.add_argument("-len_scale", "--length_scale", help="How do you want to scale the length of the quiver arrows?", type=float, default=1)
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

#fontdict = {'family': 'sans-serif',
#        'size': args.text_font,
#        }

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
                if len(input_arg) > 0:
                    input_args[-1].append(input_arg)
positions = np.array(positions)

if 0 in positions or args.plot_colourbar == False:
    make_right_gs = False
else:
    make_right_gs = True

columns = np.max(positions[:,0])
rows = np.max(positions[:,1])

if args.figure_width == None:
    width = float(columns)*(14.5/3.)
else:
    width = args.figure_width

if args.figure_height != None:
    height = args.figure_height
else:
    height = float(rows)*(17./4.)
#if args.plot_colourbar == True:
#    height = float(rows)*(17./4.)
    
plt.clf()
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
        if file_dir[it].split('.')[-1] == 'pkl':
    
            file = open(file_dir[it], 'rb')
            X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
            file.close()
            
            if int(np.round(np.max(part_info['particle_position'][1]))) == 200:
                image = np.flip(np.flip(image, axis=0), axis=1)
                vel_rad = np.flip(np.flip(vel_rad, axis=0), axis=1)
                velx = np.flip(np.flip(velx, axis=0), axis=1)
                vely = np.flip(np.flip(vely, axis=0), axis=1)
                part_info['particle_position'][1][0] = part_info['particle_position'][1][0]*-1
                #image = np.rot90(np.rot90(image))
                #vel_rad = np.rot90(np.rot90(vel_rad))
                #velx = np.rot90(np.rot90(velx))
                #vely = np.rot90(np.rot90(vely))
                #part_info['particle_position'][1][0] = part_info['particle_position'][1][0]*-1
                
            if args.convolve == 'True':
                res = (args_dict['xlim'][1] - args_dict['xlim'][0])/np.shape(vel_rad)[0]
                beam_rad = np.sqrt(12*27)/2.355
                beam_rad_pixel = beam_rad/res
                image = gaussian_filter(image,sigma=beam_rad_pixel)
                vel_rad = gaussian_filter(vel_rad,sigma=beam_rad_pixel)
            
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
            
            if '-title' in input_args[it]:
                title_ind = input_args[it].index('-title')
                title_str = input_args[it][title_ind+1]
                #title_comps = title_str.split("_")
                #title_str = ' '.join(title_comps)
            else:
                title_str = ""
                
            if '-at' in input_args[it]:
                at_ind = input_args[it].index('-at')
                annotate_time = input_args[it][at_ind+1]
            else:
                annotate_time = "True"
            
            if '-axlim' in input_args[it]:
                ax_ind = input_args[it].index('-axlim')
                axlim = int(input_args[it][ax_ind+1])
                print("updating axis limits")
            else:
                axlim = np.nan
                
            if np.isnan(axlim):
                if np.round(np.mean(args_dict['xlim'])) == np.round(np.mean(X)):
                    xlim = args_dict['xlim']
                    ylim = args_dict['ylim']
                else:
                    xlim = args_dict['xlim'] + np.mean(X)
                    ylim = args_dict['ylim'] + np.mean(Y)
            else:
                xlim = [-1*axlim, axlim]
                ylim = [-1*axlim, axlim]
            has_particles = args_dict['has_particles']
            time_val = int(args_dict['time_val'])
            time_string = f'{time_val:,}'
            xabel = 'x (au)'
            yabel = 'y (au)'
            #xabel = args_dict['xabel']
            #yabel = args_dict['yabel']
        
            bool_den_array = image>density_threshold
            vel_rad = vel_rad*bool_den_array #bool_den_array*np.nan*vel_rad
            vel_rad[vel_rad == 0] = np.nan
            
            v_std = np.std(vel_rad/100000)
            v_cbar_min = -1 #center_vel_rv.in_units('km/s').value - 5
            v_cbar_max = 1
            
            #plot = axes_dict[ax_label].pcolormesh(X, Y, vel_rad/100000, cmap='idl06_r', rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
            plot = axes_dict[ax_label].pcolormesh(X, Y, vel_rad/100000, cmap='RdYlBu_r', rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
            fmt = ticker.LogFormatterSciNotation()
            fmt.create_dummy_axis()
            print('plotted image')
                
            if density_threshold != 0.0:
                exp_min = np.log10(density_threshold)
            else:
                exp_min = np.log10(cbar_min)
            exp_max = np.log10(cbar_max)
            n_level = (exp_max-exp_min)*2 + 1
            contour_levels = np.logspace(exp_min, exp_max, int(n_level))
            CS = axes_dict[ax_label].contour(X,Y,image, locator=plt.LogLocator(), linewidths=0.5, colors='k', levels=contour_levels)
            print('plotted contours')
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
                mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])#, fontdict=fontdict)
                print('plotted particles')
                
            if annotate_time == "True":
                print("ANNONTATING TIME:", str(int(time_val))+'$\,$yr')
                time_text = axes_dict[ax_label].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), 't$\,=\,$'+time_string+'$\,$yr', va="center", ha="left", color='w', fontsize=args.text_font-0.5, zorder=10)#, fontdict=fontdict)
                time_text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'), path_effects.Normal()])
                print('plotted time')
                
            if title_str != "":
                title_text = axes_dict[ax_label].text((np.mean(xlim)), (ylim[1]-0.06*(ylim[1]-ylim[0])), title_str, va="center", ha="center", color='w', fontsize=(args.text_font+2), zorder=11)#, fontdict=fontdict)
                title_text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'), path_effects.Normal()])
                print('plotted title')
            
            if positions[it][0] == 1 and positions[it][1] == 1 and args.convolve == 'True':
                #plot convolution beam
                theta = np.linspace(0, 2*np.pi, 100)
                x1 = beam_rad*np.cos(theta) + 0.9*xlim[1]
                x2 = beam_rad*np.sin(theta) + 0.9*ylim[1]
                axes_dict[ax_label].plot(x1, x2, c='m')
            
            
            from matplotlib.ticker import AutoMinorLocator
            minor_locator = AutoMinorLocator(2)
            axes_dict[ax_label].xaxis.set_minor_locator(minor_locator)
            axes_dict[ax_label].tick_params(which='minor', color='k', direction='in', axis='x', length=3, top=True)
            axes_dict[ax_label].tick_params(which='major', color='k', direction='in', axis='x', length=3, top=True)
            axes_dict[ax_label].tick_params(which='major', color='k', direction='in', axis='y', length=3, right=True)
            '''
            axes_dict[ax_label].yaxis.set_minor_locator(minor_locator)
            axes_dict[ax_label].tick_params(which='minor', color='k', direction='in', axis='y', length=3, top=True)
            '''
        else:
            import scipy.io
            from matplotlib.patches import Ellipse
            if '-axlim' in input_args[it]:
                ax_ind = input_args[it].index('-axlim')
                axlim = int(input_args[it][ax_ind+1])
                print("updating axis limits")
            else:
                axlim = np.nan
            obs_rv = scipy.io.readsav(file_dir[it])
            image = obs_rv['mom1']
            contour_arr = obs_rv['mom0']
            con_max_it = int(np.max(obs_rv['mom0'])/obs_rv['rms'])+5 + 5
            contour_levels = np.arange(10, con_max_it, 5) * obs_rv['rms']
            x = obs_rv['xaxis']
            y = obs_rv['yaxis']
            X, Y = np.meshgrid(x, y)
            if np.isnan(axlim):
                xlim = [x[0], x[-1]]
                ylim = [y[0], y[-1]]
            else:
                xlim = [axlim, -1*axlim]
                ylim = [-1*axlim, axlim]
            vla_pos = [obs_rv['vla1xpos'], obs_rv['vla1ypos']]
            angle = 113.9
            line_length = 100
            line_pos_x = [-1*line_length*np.sin(np.radians(angle))+vla_pos[0], line_length*np.sin(np.radians(angle))+vla_pos[0]]
            line_pos_y = [-1*line_length*np.cos(np.radians(angle))+vla_pos[1], line_length*np.cos(np.radians(angle))+vla_pos[1]]
            title = str(obs_rv['molname'])[2:-1]
            if it == 3:
                title = 'HDO'
            title = r'{}'.format(title)
            
            axes_dict[ax_label].set_xlim(xlim)
            axes_dict[ax_label].set_ylim(ylim)
            
            v_cbar_min = 6 #center_vel_rv.in_units('km/s').value - 5
            v_cbar_max = 8
            
            plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap='RdYlBu_r', rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
            axes_dict[ax_label].contour(X,Y,contour_arr, locator=plt.LogLocator(), linewidths=0.5, colors='k', levels=contour_levels)

            plt.gca().set_aspect('equal')

            plt.plot(line_pos_x, line_pos_y, 'm--', lw=2)
            '''
            if '_' in title:
                title_list = title.split('$')
                title_list[1] = '_\mathregular{'.join(title_list[1].split('_')) + '}'
                title = '$'.join(title_list)
            elif '^' in title:
                title_list = title.split('$')
                if '{' in title_list[1]:
                    title_list[1] = title_list[1].split('{')[-1].split('}')[0]
                    #title_list[1] = '^\mathregular'.join(title_list[1].split('^'))
                else:
                    title_list[1] = '^\mathregular{'.join(title_list[1].split('^')) + '}'
                    title = '$'.join(title_list)
            #params = {'mathtext.default': 'regular' }
            #plt.rcParams.update(params)
            '''
            title = r'{}'.format(title)
            
            title_text = axes_dict[ax_label].text((np.mean(xlim)), (ylim[1]-0.1*(ylim[1]-ylim[0])), title, va="center", ha="center", color='k', fontsize=(args.text_font+2))
            
            from matplotlib.ticker import AutoMinorLocator
            minor_locator = AutoMinorLocator(2)
            axes_dict[ax_label].xaxis.set_minor_locator(minor_locator)
            axes_dict[ax_label].tick_params(which='minor', color='k', direction='in', axis='x', length=3, top=True)
            axes_dict[ax_label].tick_params(which='major', color='k', direction='in', axis='x', length=3, top=True)
            axes_dict[ax_label].tick_params(which='major', color='k', direction='in', axis='y', length=3, right=True)
            '''
            axes_dict[ax_label].yaxis.set_minor_locator(minor_locator)
            axes_dict[ax_label].tick_params(which='minor', color='k', direction='in', axis='y', length=3, top=True)
            '''
            if positions[it][0] == 1 and positions[it][1] == 1 and args.convolve == 'True':
                ellipse = Ellipse((0.85*xlim[1], 0.85*ylim[1]), width=13, height=28,angle=-15, edgecolor='m', facecolor='w')
                axes_dict[ax_label].add_patch(ellipse)
            
            axes_dict[ax_label].set_xticklabels(axes_dict[ax_label].get_xticks(), fontname = "Arial")
            axes_dict[ax_label].set_yticklabels(axes_dict[ax_label].get_yticks(), fontname = "Arial")
            
        if positions[it][0] == columns and positions[it][1] == rows:
            top_ax_label = ax_label[:-1]+'1'
            cax_x = np.array(axes_dict[top_ax_label].get_position())[1][0]
            cax_y = np.array(axes_dict[ax_label].get_position())[0][1]
            bar_height = np.array(axes_dict[top_ax_label].get_position())[1][1] - np.array(axes_dict[ax_label].get_position())[0][1]
            cax = f.add_axes([cax_x, cax_y, 0.02, bar_height])
            cbar = f.colorbar(plot, pad=0.0, cax=cax)
            #cax = f.add_axes([0.645, 0.11, 0.02, 0.77])
            #cbar =f.colorbar(plot, pad=0.0, cax=cax)#, label='Radial Velocity (km/s)')
            if file_dir[it].split('.')[-1] == 'pkl':
                cbar.set_label('Radial Velocity (km$\,$s$^{-1}$)', rotation=270, labelpad=12, size=args.text_font+1)#fontdict=fontdict)#
            else:
                cbar.set_label('Radial Velocity (km$\,$s$^{-1}$)', rotation=270, labelpad=12, size=args.text_font+1)
            cbar.ax.tick_params(labelsize=args.text_font)
            #import pdb
            #pdb.set_trace()
            cbar.ax.set_yticklabels(cbar.get_ticks(), fontname = "Arial")
            #cbar.ax.set_yticklabels(cbar.ax.yaxis.get_ticklabels(), fontname = "Arial")
            print('plotted colourbar')
            #cax.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=args.text_font)
        axes_dict[ax_label].set_aspect('equal')
        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel('x (au)', fontsize=args.text_font, labelpad=-1)#, fontdict=fontdict)
            print('plotted xlabel')
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel('y (au)', fontsize=args.text_font, labelpad=-10)#, fontdict=fontdict)
            print('plotted ylabel')
            
        axes_dict[ax_label].set_xlim(xlim)
        axes_dict[ax_label].set_ylim(ylim)
        
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
    elif plot_type[it] == 'global_frame':
        try:
            file = open(file_dir[it], 'rb')
            X, Y, image, particle_x_pos, particle_y_pos = pickle.load(file)
            file.close()
            
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
            
            if '-title' in input_args[it]:
                title_ind = input_args[it].index('-title')
                title = input_args[it][title_ind+1]
                #title_comps = title_str.split("_")
                #title_str = ' '.join(title_comps)
            else:
                title = ""
                
            if '-sfe' in input_args[it]:
                sfe_ind = input_args[it].index('-sfe')
                sfe = input_args[it][sfe_ind+1]
            else:
                sfe = 2.5*(positions[it][0]-1)
                
            if '-cbar_range' in input_args[it]:
                cbar_range_ind = input_args[it].index('-cbar_range')
                cbar_range = float(input_args[it][cbar_range_ind+1])
            else:
                cbar_range = 1.5
                
            xlim = [np.min(X), np.max(X)]
            ylim = [np.min(Y), np.max(Y)]
                
            M_gas = [1500, 3000, 3750, 4500, 6000, 12000]
            volume = yt.YTQuantity(4**3, 'pc**3')
            mass = yt.YTQuantity(M_gas[positions[it][1]-1], 'Msun')
            mean_dens = mass/volume
            mean_log = np.log10(mean_dens.in_units('g/cm**3'))
            cbar_min = 10**(mean_log-cbar_range)
            cbar_max = 10**(mean_log+cbar_range)
            plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
            plt.gca().set_aspect('equal')
            if positions[it][0] == columns:
                cbar = plt.colorbar(plot, pad=0.0)
                cbar.ax.tick_params(labelsize=args.text_font)
                cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=20, size=args.text_font)
            
            title_text = axes_dict[ax_label].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[0]+0.05*(ylim[1]-ylim[0])), '$SFE$='+str(sfe)+'\%', va="center", ha="left", color='w', fontsize=(args.text_font))
            title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            time_text = axes_dict[ax_label].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1] - 0.05*(ylim[1]-ylim[0])), title, va="center", ha="left", color='w', fontsize=args.text_font+2, zorder=4)#, fontdict=fontdict)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            axes_dict[ax_label].scatter((particle_x_pos.value - 2), (particle_y_pos.value - 2), color='c', s=0.5)
            
            if positions[it][1] == rows:
                axes_dict[ax_label].set_xlabel("X (pc)", labelpad=-1, fontsize=args.text_font)
            if positions[it][0] == 1:
                axes_dict[ax_label].set_ylabel("Y (pc)", fontsize=args.text_font, labelpad=-9)
        except:
            plt.gca().set_aspect('equal')
            print("pickle doesn't exist")

    elif plot_type[it] == 'movie_frame':
        if positions[it][1] > 1:
            prev_ax_label = ax_label[:-1]+str(int(ax_label[-1])-1)
            hsep = np.array(axes_dict[prev_ax_label].get_position())[0][1] - np.array(axes_dict[ax_label].get_position())[1][1]
            if hsep != ghspace:
                d_sep = hsep - ghspace
                new_pos = np.array(axes_dict[ax_label].get_position())
                new_pos[:,1] = new_pos[:,1] + d_sep
                left = new_pos[0][0]
                bottom = new_pos[0][1]
                width = new_pos[1][0] - left
                height = new_pos[1][1] - bottom
                axes_dict[ax_label].set_position([left, bottom, width, height])
    
        file = open(file_dir[it], 'rb')
        X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
        file.close()
        
        if '-pvl' in input_args[it]:
            pvl_ind = input_args[it].index('-pvl')
            plot_vel_legend = eval(input_args[it][pvl_ind+1])
        else:
            plot_vel_legend = False
        
        if '-stdv' in input_args[it]:
            stdv_ind = input_args[it].index('-stdv')
            standard_vel = eval(input_args[it][stdv_ind+1])
        else:
            standard_vel = 2
            
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
        
        if '-title' in input_args[it]:
            title_ind = input_args[it].index('-title')
            title = input_args[it][title_ind+1]
            #title_comps = title_str.split("_")
            #title_str = ' '.join(title_comps)
        else:
            title = ""
            
        if '-at' in input_args[it]:
            at_ind = input_args[it].index('-at')
            annotate_time = eval(input_args[it][at_ind+1])
        else:
            annotate_time = False
        
        if '-axlim' in input_args[it]:
            ax_ind = input_args[it].index('-axlim')
            axlim = int(input_args[it][ax_ind+1])
            print("updating axis limits")
        else:
            axlim = np.nan
        
        if '-apm' in input_args[it]:
            apm_ind = input_args[it].index('-apm')
            plot_particle_mass = eval(input_args[it][apm_ind+1])
        else:
            plot_particle_mass = False
        
        if '-col_title' in input_args[it]:
            col_it = input_args[it].index('-col_title')
            col_title = input_args[it][col_it+1]
        else:
            col_title = ""
        
        xlim = args_dict['xlim']
        ylim = args_dict['ylim']
        if np.round(np.mean(args_dict['xlim'])) != np.round(np.mean(X)):
            shift_x = np.mean(X)
            shift_y = np.mean(Y)
            X = X - shift_x
            Y = Y - shift_y
            X_vel = X_vel - shift_x
            Y_vel = Y_vel - shift_y
            part_info['particle_position'][0] = part_info['particle_position'][0] - shift_x
            part_info['particle_position'][1] = part_info['particle_position'][1] - shift_y
              
        has_particles = args_dict['has_particles']
        #xabel = args_dict['xabel']
        xabel = 'x (au)'
        yabel = args_dict['yabel']

        if positions[it][1] == rows:
            axes_dict[ax_label].set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
        #if positions[it][0] == 1:
        #    axes_dict[ax_label].set_ylabel(yabel, fontsize=args.text_font, labelpad=-10)
        axes_dict[ax_label].set_xlim(xlim)
        axes_dict[ax_label].set_ylim(ylim)
        
        cbar_min = args_dict['cbar_min']
        cbar_max = args_dict['cbar_max']
        
        if 0.0 in (cbar_min, cbar_max) or len(np.where(np.array([cbar_min, cbar_max]) < 0)[0]) > 0 :
            plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.bwr, rasterized=True, vmin=cbar_min, vmax=cbar_max)
        else:
            plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
        plt.gca().set_aspect('equal')
        
        del image
        
        del X
        del Y
        del magx
        del magy
        
        if positions[it][0] == columns:
            cax_x = np.array(axes_dict[ax_label].get_position())[1][0]
            cax_y = np.array(axes_dict[ax_label].get_position())[0][1]
            bar_height = np.array(axes_dict[ax_label].get_position())[1][1] - np.array(axes_dict[ax_label].get_position())[0][1]
            cax = f.add_axes([cax_x, cax_y, 0.01, bar_height])
            cbar =f.colorbar(plot, pad=0.0, cax=cax)#, label='Radial Velocity (km/s)')
            #cbar.set_label('Radial Velocity (km$\,$s$^{-1}$)', rotation=270, labelpad=20, size=args.text_font)#fontdict=fontdict)#
            #print('plotted colourbar')
            #cbar = plt.colorbar(plot, pad=0.0, ax=axes_dict[ax_label])
            cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
            cbar.ax.tick_params(labelsize=args.text_font)
        
        mym.my_own_quiver_function(axes_dict[ax_label], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_vel_legend, standard_vel=standard_vel, limits=[xlim, ylim])
        
        del X_vel
        del Y_vel
        del velx
        del vely
        
        title_text = axes_dict[ax_label].text((np.mean(xlim)), (ylim[1]-0.05*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+2))
        title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        
        if col_title != "":
            col_title = " ".join(col_title.split("_"))
            axes_dict[ax_label].set_title(col_title, fontsize=(args.text_font+2))
            
        if annotate_time:
            time_val = int(args_dict['time_val'])
            time_string = f'{time_val:,}'
            time_string = "t$\,$=$\,$"+time_string+"$\,$yr"
            time_string_raw = r"{}".format(time_string)
            #import pdb
            #pdb.set_trace()
            time_text = axes_dict[ax_label].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        
        if plot_particle_mass:
            annotate_field = part_info['particle_mass']
        else:
            annotate_field = None

        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=annotate_field, particle_tags=part_info['particle_tag'])
        
        if it == 0:
            rect = patches.Rectangle((750, 750), 150, 150, linewidth=2, edgecolor='w', facecolor='none')
            axes_dict[ax_label].add_patch(rect)
            
        del part_info
        
        
        if positions[it][0] == 1:
            axes_dict[ax_label].set_ylabel('y (au)', labelpad=-20, fontsize=args.text_font)
        
        if positions[it][0] != 1:
            yticklabels = axes_dict[ax_label].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        if positions[it][1] == rows:
            axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
            #if positions[it][0] != 1:
            #    xticklabels = axes_dict[ax_label].get_xticklabels()
            #    plt.setp(xticklabels[0], visible=False)
        plt.tick_params(axis='both', which='major', labelsize=args.text_font)
        axes_dict[ax_label].tick_params(which='both', direction='in', color='white')
        for line in axes_dict[ax_label].xaxis.get_ticklines():
            line.set_color('white')
        for line in axes_dict[ax_label].yaxis.get_ticklines():
            line.set_color('white')
            
        axes_dict[ax_label].tick_params(which='major', color='w', direction='in', axis='x', length=3, top=True)
        axes_dict[ax_label].tick_params(which='major', color='w', direction='in', axis='y', length=3, right=True)
    
    if positions[it][0] != 1:
        yticklabels = axes_dict[ax_label].get_yticklabels()
        plt.setp(yticklabels, visible=False)
        yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
        plt.setp(yticklabels, visible=False)
    if positions[it][0] == 1:
        axes_dict[ax_label].tick_params(axis='y', which='major', labelsize=args.text_font)
        if positions[it][1] > 1 and plot_type[it] != 'movie_frame':
            yticklabels = axes_dict[ax_label].get_yticklabels()
            plt.setp(yticklabels[-1], visible=False)
    if positions[it][1] != rows:
        xticklabels = axes_dict[ax_label].get_xticklabels()
        plt.setp(xticklabels, visible=False)
    if positions[it][1] == rows:
        axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
        if positions[it][0] != 1:
            xticklabels = axes_dict[ax_label].get_xticklabels()
            if file_dir[it].split('.')[-1] == 'idl':
                plt.setp(xticklabels[-1], visible=False)
            else:
                plt.setp(xticklabels[0], visible=False)
                
    ax_r = axes_dict[ax_label].secondary_yaxis('right')
    ax_t = axes_dict[ax_label].secondary_xaxis('top')
    xticklabels = ax_t.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    yticklabels = ax_r.get_yticklabels()
    plt.setp(yticklabels, visible=False)
    ax_r.tick_params(axis='y', direction='in')
    ax_t.tick_params(axis='x', direction='in')
    
    axes_dict[ax_label].tick_params(axis='x', direction='in')
    axes_dict[ax_label].tick_params(axis='y', direction='in')

    if args.multiplot_title != "":
        plt.suptitle(args.multiplot_title, x=0.45, y=0.9, fontsize=18)

    #f.savefig(savename + '.pdf', format='pdf')
    #f.savefig(savename + '.eps', format='eps')
    print("saving", savename + '.pdf')
    #plt.savefig(savename + '.eps', format='eps', bbox_inches='tight', pad_inches = 0.02)
    #if it == (len(positions) - 1):
    plt.savefig(savename + '.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)
    plt.savefig(savename + '.png', format='png', bbox_inches='tight', pad_inches = 0.02, dpi=500)
    print("saved", savename + '.pdf')
    #f.savefig(savename + '.eps', format='eps', bbox_inches='tight', pad_inches = 0.02)
