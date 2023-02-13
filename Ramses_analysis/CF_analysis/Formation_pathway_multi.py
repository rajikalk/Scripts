#!/usr/bin/env python
import numpy as np
import sys
import pickle
import csv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
import matplotlib
import collections


def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

#matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
#matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
#matplotlib.rcParams['mathtext.rm'] = 'Arial'
#matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
#matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
#matplotlib.rcParams['mathtext.rm'] = 'Arial'
#matplotlib.rcParams['mathtext.sf'] = 'Arial'
#matplotlib.rcParams['mathtext.default'] = 'regular'
#matplotlib.rcParams['font.sans-serif'] = 'Arial'
#matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
#       r'\usepackage{helvet}',    # set the normal font here
#       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
#]

plot_pickles = ['/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/bound_core_frag_(106_77)_1_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/bound_core_frag_(106_77)_2_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/bound_core_frag_(106_77)_3_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/unbound_core_frag_121_104_1_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/unbound_core_frag_121_104_2_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/unbound_core_frag_121_104_3_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/dynamical_capt_(101_(13_[77_106]))_1_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/dynamical_capt_(101_(13_[77_106]))_2_all.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/dynamical_capt_(101_(13_[77_106]))_3_all.pkl']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(two_col_width,two_col_width))

for pick_it in range(len(plot_pickles)):
    pickle_file = plot_pickles[pick_it]
    file = open(pickle_file, 'rb')
    try:
        system, particle_x_pos, particle_y_pos, particle_masses, center_pos, thickness, X, Y, image, existing_sinks, time_val = pickle.load(file)
    except:
        file.close()
        file = open(pickle_file, 'rb')
        system, particle_x_pos, particle_y_pos, particle_masses, center_pos, thickness, X, Y, image, time_val = pickle.load(file)
    file.close()
    
    xlim = [-1*thickness, thickness]
    ylim = [-1*thickness, thickness]
    x = np.linspace(xlim[0], xlim[1], 800)
    y = np.linspace(ylim[0], ylim[1], 800)
    X, Y = np.meshgrid(x, y)
    
    X = X + center_pos[0]
    Y = Y + center_pos[1]
    
    X = X/100000
    Y = Y/100000
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    
    cmin = 10**(np.log10(np.mean(image))-1.5)
    cmax = 10**(np.log10(np.mean(image))+1.5)
    
    plot = axs.flatten()[pick_it].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cmin, vmax=cmax), rasterized=True)
    #plt.savefig('formation_pathways.png', format='png', bbox_inches='tight')
    plt.gca().set_aspect('equal')
    
    #cbar = plt.colorbar(plot, pad=0.0)
    #cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
    
    time_string = "$t$="+str(int(time_val))+"yr"
    time_string_raw = r"{}".format(time_string)
    time_text = axs.flatten()[pick_it].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    
    #xabel = 'X (AU)'
    #yabel = 'Y (AU)'
    
    if '/bound_' in pickle_file:
        if len(particle_x_pos) > 1:
            axs.flatten()[pick_it].plot(particle_x_pos, particle_y_pos, linestyle='-', color='grey')
            
        #part_color = ['cyan','magenta','r','b','y','w','k']
        axs.flatten()[pick_it].scatter(particle_x_pos, particle_y_pos, c='y', marker='*', s=100, linewidth=1.5, edgecolor="k", zorder=11)
    elif '/unbound_' in pickle_file:
        if len(particle_x_pos) == 2 and '1_all.pkl' not in pickle_file:
            if '3_all.pkl' in pickle_file:
                linestyle = '-'
            elif '2_all.pkl' in pickle_file:
                linestyle = ':'
            axs.flatten()[pick_it].plot(particle_x_pos, particle_y_pos, linestyle=linestyle, color='grey')
        elif len(particle_x_pos) > 2:
            if '1_all.pkl' not in pickle_file:
                #plot lines between system:
                sys_string = str(system[0][1])
                reduced = False
                replace_int = 999
                while reduced == False:
                    open_bracket_pos = []
                    for char_it in range(len(sys_string)):
                        if sys_string[char_it] == '[':
                            open_bracket_pos.append(char_it)
                        if sys_string[char_it] == ']':
                            open_ind = open_bracket_pos.pop()
                            sub_sys = sys_string[open_ind:char_it+1]
                            first_ind = existing_sinks.index(eval(sub_sys)[0])
                            second_ind = existing_sinks.index(eval(sub_sys)[1])
                            sub_inds = np.array([first_ind, second_ind])
                            axs.flatten()[pick_it].plot(particle_x_pos[sub_inds], particle_y_pos[sub_inds], 'b-', alpha=0.5)
                            
                            x_com = np.sum(particle_x_pos[sub_inds] * particle_masses[sub_inds])/np.sum(particle_masses[sub_inds])
                            y_com = np.sum(particle_y_pos[sub_inds] * particle_masses[sub_inds])/np.sum(particle_masses[sub_inds])
                            com_mass = np.sum(particle_masses[sub_inds])
                            replace_string = str(replace_int)
                            existing_sinks.append(replace_int)
                            replace_int = replace_int - 1

                            particle_masses = yt.YTArray(list(particle_masses.value)+[com_mass.value], 'Msun')
                            particle_x_pos = yt.YTArray(list(particle_x_pos.value)+[x_com.value], 'au')
                            particle_y_pos = yt.YTArray(list(particle_y_pos.value)+[y_com.value], 'au')
                            
                            sys_string = sys_string[:open_ind] + replace_string + sys_string[char_it+1:]
                            if '[' not in sys_string:
                                reduced = True
                            break
                    target_sink_ind = existing_sinks.index(system[0][0])
                    if '3_all.pkl' in pickle_file:
                        linestyle = '-'
                    elif '2_all.pkl' in pickle_file:
                        linestyle = ':'
                    axs.flatten()[pick_it].plot(particle_x_pos[np.array([target_sink_ind, -1])], particle_y_pos[np.array([target_sink_ind, -1])], linestyle=linestyle, color='grey')
        
        axs.flatten()[pick_it].scatter(particle_x_pos, particle_y_pos, c='y', marker='*', s=100, linewidth=1.5, edgecolor="k", zorder=11)
    else:
        if len(particle_x_pos) == 2 and '1_all.pkl' not in pickle_file:
            if '3_all.pkl' in pickle_file:
                linestyle = 'b-'
            elif '2_all.pkl' in pickle_file:
                linestyle = 'b:'
            axs.flatten()[pick_it].plot(particle_x_pos, particle_y_pos, linestyle)
        elif len(particle_x_pos) > 2:
            if '1_all.pkl' not in pickle_file:
                #plot lines between system:
                sys_string = str(system[0][1])
                reduced = False
                replace_int = 999
                while reduced == False:
                    open_bracket_pos = []
                    for char_it in range(len(sys_string)):
                        if sys_string[char_it] == '[':
                            open_bracket_pos.append(char_it)
                        if sys_string[char_it] == ']':
                            open_ind = open_bracket_pos.pop()
                            sub_sys = sys_string[open_ind:char_it+1]
                            first_ind = existing_sinks.index(eval(sub_sys)[0])
                            second_ind = existing_sinks.index(eval(sub_sys)[1])
                            sub_inds = np.array([first_ind, second_ind])
                            axs.flatten()[pick_it].plot(particle_x_pos[sub_inds], particle_y_pos[sub_inds], 'b-', alpha=0.5)
                            
                            x_com = np.sum(particle_x_pos[sub_inds] * particle_masses[sub_inds])/np.sum(particle_masses[sub_inds])
                            y_com = np.sum(particle_y_pos[sub_inds] * particle_masses[sub_inds])/np.sum(particle_masses[sub_inds])
                            com_mass = np.sum(particle_masses[sub_inds])
                            replace_string = str(replace_int)
                            existing_sinks.append(replace_int)
                            replace_int = replace_int - 1

                            particle_masses = yt.YTArray(list(particle_masses.value)+[com_mass.value], 'Msun')
                            particle_x_pos = yt.YTArray(list(particle_x_pos.value)+[x_com.value], 'au')
                            particle_y_pos = yt.YTArray(list(particle_y_pos.value)+[y_com.value], 'au')
                            
                            sys_string = sys_string[:open_ind] + replace_string + sys_string[char_it+1:]
                            if '[' not in sys_string:
                                reduced = True
                            break
                    target_sink_ind = existing_sinks.index(system[0][0])
                    if '3_all.pkl' in pickle_file:
                        linestyle = 'b-'
                    elif '2_all.pkl' in pickle_file:
                        linestyle = 'b:'
                    axs.flatten()[pick_it].plot(particle_x_pos[np.array([target_sink_ind, -1])], particle_y_pos[np.array([target_sink_ind, -1])], linestyle)
                    
        other_ind = system[0][1][1]
        if '[' in other_ind:
            other_ind = flatten(eval(other_ind))
        else:
            other_ind = [eval(other_ind)]
        plot_other_inds = list(set(other_ind).intersection(existing_sinks))
        if len(plot_other_inds) > 0:
            for plot_other_ind in plot_other_inds:
                axs.flatten()[pick_it].scatter(particle_x_pos[existing_sinks.index(plot_other_ind)], particle_y_pos[existing_sinks.index(plot_other_ind)], c='r', marker='*', s=100, linewidth=1.5, edgecolor="k", zorder=11)
        birth_sys = [system[0][0], system[0][1][0]]
        birth_inds = list(set(birth_sys).intersection(existing_sinks))
        for birth_ind in birth_inds:
            axs.flatten()[pick_it].scatter(particle_x_pos[existing_sinks.index(birth_ind)], particle_y_pos[existing_sinks.index(birth_ind)], c='y', marker='*', s=100, linewidth=1.5, edgecolor="k", zorder=11)
        
    plt.tick_params(axis='both', which='major')# labelsize=16)
    for line in axs.flatten()[pick_it].xaxis.get_ticklines():
        line.set_color('white')
    for line in axs.flatten()[pick_it].yaxis.get_ticklines():
        line.set_color('white')
    axs.flatten()[pick_it].tick_params(direction='in', color='white')
    if pick_it == 0:
        axs.flatten()[pick_it].annotate(r'$\times$10$^5$', xy=(np.min(X).value-1, np.max(Y).value))
    if pick_it == 8:
        axs.flatten()[pick_it].annotate(r'$\times$10$^5$', xy=(np.max(X).value, np.min(Y).value-1))
    #axs.flatten()[pick_it].ticklabel_format(axis='both', style='sci', scilimits=(4,4))
        
    plt.savefig("formation_pathways.png", format='png', bbox_inches='tight')
    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
    print('updated "formation_pathways.png')
        

