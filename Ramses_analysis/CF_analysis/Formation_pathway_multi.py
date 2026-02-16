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
import yt

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

plot_pickles = ['/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/bound_core_frag_(106_77)_1_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/bound_core_frag_(106_77)_2_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/bound_core_frag_(106_77)_3_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/unbound_core_frag_121_104_1_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/unbound_core_frag_121_104_2_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/unbound_core_frag_121_104_3_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/dynamical_capt_(101_(13_[77_106]))_1_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/dynamical_capt_(101_(13_[77_106]))_2_all.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/Formation_pathways/Plot_pickles/dynamical_capt_(101_(13_[77_106]))_3_all.pkl']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

pathway_label = ['Bound Core Fragmentation', 'Unbound Core Fragmentation', 'Dynamical Capture']

plt.clf()
fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(two_col_width,two_col_width))
#plt.subplots_adjust(wspace=0.12, hspace=0.3)
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
    
    #get values for scale bar
    upper_x = np.max(X) - 0.05*(np.max(X) - np.min(X))
    lower_x = upper_x - yt.YTQuantity(1000, 'au')
    y_scale = np.min(Y) + 0.1*(np.max(Y) - np.min(Y))
    x_scale = np.array([lower_x, upper_x])/10000
    y_scale = y_scale/10000
    
    X = X/10000
    Y = Y/10000
    particle_x_pos = particle_x_pos/10000
    particle_y_pos = particle_y_pos/10000
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    
    cmin = 10**(np.log10(np.mean(image))-1.5)
    cmax = 10**(np.log10(np.mean(image))+1.5)
    
    plot = axs.flatten()[pick_it].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cmin, vmax=cmax), rasterized=True)
    
    time_string = "$t$="+str(int(time_val))+"yr"
    time_string_raw = r"{}".format(time_string)
    time_text = axs.flatten()[pick_it].text((xlim[0]+0.02*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    
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
                            axs.flatten()[pick_it].plot(particle_x_pos[sub_inds], particle_y_pos[sub_inds], color='grey', alpha=0.5)
                            
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
                            axs.flatten()[pick_it].plot(particle_x_pos[sub_inds], particle_y_pos[sub_inds], color='grey', alpha=0.5)
                            
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
    #if pick_it == 0:
    #    axs.flatten()[pick_it].text((xlim[0]-0.16*(xlim[1]-xlim[0])), (ylim[1]+0.01*(ylim[1]-ylim[0])), r"$\times$10$^5$", va="center", ha="left", color='k', fontsize=10)
    #if pick_it == 8:
    #    axs.flatten()[pick_it].text((xlim[1]-0.04*(xlim[1]-xlim[0])), (ylim[0]-0.06*(ylim[1]-ylim[0])), r"$\times$10$^5$", va="center", ha="left", color='k', fontsize=10)
        
    if np.remainder(pick_it, 3) == 1:
        axs.flatten()[pick_it].set_title(pathway_label[int(pick_it/3)], pad=-2)
    #axs.flatten()[pick_it].ticklabel_format(axis='both', style='sci', scilimits=(4,4))
    axs.flatten()[pick_it].set_xlim(xlim)
    axs.flatten()[pick_it].set_ylim(ylim)
    
    xabel = r"X (10$^4$AU)"
    yabel = r"Y (10$^4$AU)"
    if np.remainder(pick_it,3)==0:
        axs.flatten()[pick_it].set_ylabel(yabel, fontsize=10)
    if pick_it > 5:
        axs.flatten()[pick_it].set_xlabel(xabel, fontsize=10)
        
    if pick_it == 2:
        yticklabels = axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[1], visible=False)
    
    plt.gca().set_aspect('equal')
    
    axs.flatten()[pick_it].tick_params(axis='y', rotation=90)
    
    axs.flatten()[pick_it].axhline(y=y_scale, xmin=x_scale[0], xmax=x_scale[1], color='grey', linewidth=2, zorder=11)
    axs.flatten()[pick_it].plot(x_scale, np.array([y_scale, y_scale]), color='grey', linewidth=3, zorder=11)
    
        
    plt.savefig("formation_pathways.pdf", format='pdf', bbox_inches='tight', pad=0.02)
    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
    print('updated "formation_pathways.png')
        

