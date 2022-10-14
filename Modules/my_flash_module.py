#!/usr/bin/env python
import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as path_effects

def my_own_quiver_function(axis, X_pos, Y_pos, X_val, Y_val, plot_velocity_legend='False', standard_vel=5, limits=None, Z_val=None, width_ceil = 0.8, zorder=3):
    global fontsize_global
    if plot_velocity_legend == 'False' or plot_velocity_legend == False:
        plot_velocity_legend = False
    elif plot_velocity_legend == 'True' or plot_velocity_legend == True:
        plot_velocity_legend = True
        if standard_vel > 1.0:
            legend_text=str(int(standard_vel)) + "km$\,$s$^{-1}$"
        else:
            legend_text=str(int(standard_vel*10.)/10.) + "km$\,$s$^{-1}$"
    standard_vel = yt.units.km.in_units('cm').value * standard_vel
    if limits is None:
        xmin = np.min(X_pos)
        xmax = np.max(X_pos)
        ymin = np.min(Y_pos)
        ymax = np.max(Y_pos)
    else:
        xmin = limits[0][0]
        xmax = limits[0][1]
        ymin = limits[1][0]
        ymax = limits[1][1]
    #rv_colors = np.linspace(-1, 1, 256)
    rv_colors = np.linspace(0, 1, 256)
    rv_cmap = plt.cm.get_cmap('RdGy_r')#plt.cm.get_cmap('bwr')
    
    len_scale = (0.07*(xmax - xmin))
    #len_scale = length_scale*standard_vel/(0.07*(xmax - xmin))
    #vels = np.hypot(X_val, Y_val)
    for xp in range(len(X_pos[0])):
        for yp in range(len(Y_pos[0])):
            xvel = len_scale*(X_val[xp][yp]/standard_vel)
            yvel = len_scale*(Y_val[xp][yp]/standard_vel)
            #xvel = length_scale*X_val[xp][yp]/len_scale
            #yvel = length_scale*Y_val[xp][yp]/len_scale
            #width_val = (np.sqrt(X_val[xp][yp]**2. + Y_val[xp][yp]**2.)/standard_vel)**2.
            width_val = np.sqrt(X_val[xp][yp]**2. + Y_val[xp][yp]**2.)/standard_vel
            if width_val > width_ceil:
                width_val = width_ceil
            try:
                if Z_val == None:
                    color = 'w'
                    #color = 'k'
            except:
                #cmap = 'idl06_r'
                zvel = Z_val[xp][yp]/len_scale
                cit = np.argmin(abs(rv_colors - zvel))
                color = rv_cmap(cit)
            axis.add_patch(mpatches.FancyArrowPatch((X_pos[xp][yp], Y_pos[xp][yp]), (X_pos[xp][yp]+xvel, Y_pos[xp][yp]+yvel), color=color, linewidth=width_val, arrowstyle='->', mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0, alpha=width_val/width_ceil))
    if plot_velocity_legend:
        #print("plotting quiver legend")
        #pos_start = [xmax - 0.15*(xmax-xmin), ymin + (fontsize_global/100)*(ymax-ymin)]
        pos_start = [xmax - 0.2*(xmax-xmin), ymin + (fontsize_global/100)*(ymax-ymin)]
        xvel = len_scale*(standard_vel/standard_vel)
        yvel = 0.0
        width_val = width_ceil
        annotate_text = axis.text((xmax - 0.01*(xmax-xmin)), (ymin + 0.03*(ymax-ymin)), legend_text, va="center", ha="right", color='w', fontsize=fontsize_global)
        annotate_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=width_val, mutation_scale=10.*width_val, alpha=width_val/width_ceil))
    return axis

def annotate_particles(axis, particle_position, accretion_rad, limits, annotate_field=None, field_symbol="M", units=None, particle_tags=None, lw=1.5, zorder=4):
    global fontsize_global
    if annotate_field is not None and units is not None:
        annotate_field = annotate_field.in_units(units)
    #part_color = ['cyan','magenta','r','b','y','w','k']
    part_color = ['r','r','r','b','y','w','k']
    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    box_size = xmax - xmin
    if accretion_rad/box_size < 0.05:
        line_rad = 0.005*box_size
    else:
        line_rad = accretion_rad
    try:
        units = str(annotate_field.unit_quantity).split(' ')[-1]
        s1 = units.split('**')
        unit_string = ""
        for s in s1:
            if s == s1[0]:
                unit_string = unit_string + s
            else:
                unit_string = unit_string + "$^" + s[0] + "$" + s[1:]

    except:
        unit_string = "$\,$M$_\odot$"
    field_symbol = "$" + field_symbol + "_"
    p_t = ""
    rainbow_text_colors = []
    string_pos = []
    if particle_tags == None:
        particle_tags = np.arange(len(particle_position))
    for pos_it in np.argsort(particle_tags):
        axis.scatter(particle_position[0][pos_it], particle_position[1][pos_it], c=part_color[pos_it], s=1, zorder=zorder)
        #axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=lw, c='k')
        #axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=lw, c='k')
        #axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=lw/2, c=part_color[pos_it])
        #axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=lw/2, c=part_color[pos_it])
        #circle = mpatches.Circle([particle_position[0][pos_it], particle_position[1][pos_it]], accretion_rad, fill=False, lw=lw/2, edgecolor='k')
        #axis.add_patch(circle)
        if annotate_field is not None:
            if units is not None:
                annotate_field = annotate_field.in_units(units)
            if unit_string == "$\,$M$_\odot$":
                P_msun = str(np.round(annotate_field[pos_it], 2))
                if len(P_msun.split('.')[-1]) == 1:
                    P_msun = P_msun +"0"
            else:
                if annotate_field[pos_it] == 0.0:
                    P_msun = "0.0"
                else:
                    P_msun = "{:0.1f}".format(annotate_field[pos_it])
            if p_t == "":
                p_t = field_symbol+str(pos_it+1)+"$ =$\,$"+P_msun+unit_string
            else:
                p_t = p_t+", "+field_symbol+str(pos_it+1)+"$ =$\,$"+P_msun+unit_string
            rainbow_text_colors.append(part_color[pos_it])
            rainbow_text_colors.append('white')
    if annotate_field is not None:
        if len(particle_tags) > 3:
            string_l = p_t[:68]
            string_2 = p_t[69:]
            colors_1 = rainbow_text_colors[:9]
            colors_2 = rainbow_text_colors[9:]
            rainbow_text((xmin + 0.01*(box_size)), (ymin + 0.025*(ymax-ymin)*3), string_l.split(' '), colors_1, size=fontsize_global, zorder=10)
            rainbow_text((xmin + 0.01*(box_size)), (ymin + 0.025*(ymax-ymin)), string_2.split(' '), colors_2, size=fontsize_global, zorder=10)
        else:
            #try:
            rainbow_text((xmin + 0.01*(box_size)), (ymin + 0.025*(ymax-ymin)), p_t.split(' '), rainbow_text_colors, size=fontsize_global, zorder=10)
            #except:
            #    print("couldn't annotate particle masses")
    return axis
