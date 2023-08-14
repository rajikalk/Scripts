#!/usr/bin/env python
import yt
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as path_effects
from matplotlib import transforms
from matplotlib.transforms import Affine2D, offset_copy

fontsize_global=12

def set_global_font_size(x):
    global fontsize_global
    fontsize_global = x
    return fontsize_global

def get_global_font_size():
    global fontsize_global
    return fontsize_global

def find_sink_formation_time(files):
    file = files[-2]
    part_file=file[:-12] + 'part' + file[-5:]

    f = h5py.File(part_file, 'r')
    if f[list(f.keys())[1]][-1][-1] > 0:
        sink_form = np.min(f[list(f.keys())[11]][:,5])*yt.units.s
        sink_form = sink_form.in_units('yr').value
    else:
        sink_form = 0.0
    f.close()
    return sink_form

def generate_frame_times(files, dt, start_time=0, presink_frames=25, end_time=None, form_time=None):
    if form_time != None:
        sink_form_time = form_time
    else:
        sink_form_time = find_sink_formation_time(files)
        
    if end_time != None:
        max_time = end_time
    else:
        file = files[-1]
        part_file=file[:-12] + 'part' + file[-5:]
        f = h5py.File(part_file, 'r')
        max_time = (f[list(f.keys())[7]][0][-1]*yt.units.s).in_units('yr').value - sink_form_time
        f.close()

    if presink_frames != 0:
        m_times = np.logspace(0.0, np.log10(sink_form_time), presink_frames) - sink_form_time
    else:
        m_times = np.array([])
    m_times = m_times.tolist()

    postsink = 0.0
    while postsink <= max_time:
        if start_time is not None:
            if postsink >= start_time:
                m_times.append(postsink)
        else:
            m_times.append(postsink)
        postsink = postsink + dt
    return m_times

def find_files(m_times, files):
    sink_form_time = find_sink_formation_time(files)
    if m_times[0] > 10e5:
        m_times[0] = m_times[0]/yt.units.yr.in_units('s').value
    usable_files = []
    mit = 0
    min = 0
    max = len(files)-1
    pit = 0
    prev_diff = np.nan
    while mit < len(m_times):
        it = int(np.round(min + ((max - min)/2.)))
        #print 'search iterator =', it
        file = files[it]
        part_file=file[:-12] + 'part' + file[-5:]
        f = h5py.File(part_file, 'r')
        time = (f['real scalars'][0][1]*yt.units.s).in_units('yr').value-sink_form_time
        f.close()
        print("Current file time =", time, "for interator =", it)
        if pit == it or time == m_times[mit]:
            if it == (len(files)-1):
                pot_files = files[it-1:it]
            elif it == 0:
                pot_files = files[it:it+1]
            else:
                pot_files = files[it-2:it+2]
            diff_arr = []
            for pfile in pot_files:
                part_file=pfile[:-12] + 'part' + pfile[-5:]
                f = h5py.File(part_file, 'r')
                time = (f['real scalars'][0][1]*yt.units.s).in_units('yr').value-sink_form_time
                f.close()
                diff_val = abs(time - m_times[mit])
                diff_arr.append(diff_val)
            append_file = pot_files[np.argmin(diff_arr)]
            if m_times[mit] == 0.0:
                try:
                    f = h5py.File(append_file, 'r')
                    print("found file", append_file, ", checking if particles exist")
                    if f[list(f.keys())[9]][-1][-1] == 0:
                        found_first_particle_file = False
                        app_ind = files.index(append_file) + 1
                        while found_first_particle_file == False:
                            append_file = files[app_ind]
                            f = h5py.File(append_file, 'r')
                            if f[list(f.keys())[9]][-1][-1] > 0:
                                found_first_particle_file = True
                                print("found particles in file", append_file)
                            else:
                                app_ind = app_ind + 1
                except:
                    f = h5py.File(append_file, 'r')
                    if 'particlemasses' not in list(f.keys()):
                        append_file = pot_files[np.argmin(diff_arr)+1]
                    f.close()
            usable_files.append(append_file)
            part_file=append_file[:-12] + 'part' + append_file[-5:]
            f = h5py.File(part_file, 'r')
            time = (f['real scalars'][0][1]*yt.units.s).in_units('yr').value-sink_form_time
            f.close()
            print("found time", time, "for m_time", m_times[mit], "with file:", usable_files[-1])
            #append_file = files[it]
            #usable_files.append(append_file)
            mit = mit + 1
            min = it
            max = len(files)-1
            pit = it
            prev_diff = np.nan
        elif time > m_times[mit]:
            max = it
            pit = it
        elif time < m_times[mit]:
            min = it
            pit = it
    return usable_files

"""
def rainbow_text(x,y,ls,lc,**kw):
    t = plt.gca().transData
    figlocal = plt.gcf()
    space_size = 1.45*kw['size']
            
    #horizontal version
    for string,c in zip(ls,lc):
        string_raw = r'{}'.format(string)
        #string = str_text[1:-1]
        #string = string.encode('unicode_escape')
        text = plt.text(x,y,string_raw,color=c, transform=t, **kw)
        text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'), path_effects.Normal()])
        text.draw(figlocal.canvas.get_renderer())
        ex = text.get_window_extent(renderer=figlocal.canvas.get_renderer())
        if "odot" in string:
            #import pdb
            #pdb.set_trace()
            #t = transforms.offset_copy(text._transform, x=ex.width, units='dots')
            t = transforms.offset_copy(text._transform, x=0.75*ex.width, units='dots')ßß
        else:
            #import pdb
            #pdb.set_trace()
            #t = transforms.offset_copy(text._transform, x=space_size, units='dots')
            t = transforms.offset_copy(text._transform, x=0.75*ex.width, units='dots')
"""
def rainbow_text(x, y, strings, colors, orientation='horizontal',
                 ax=None, **kwargs):
    """
    Take a list of *strings* and *colors* and place them next to each
    other, with text strings[i] being shown in colors[i].

    Parameters
    ----------
    x, y : float
        Text position in data coordinates.
    strings : list of str
        The strings to draw.
    colors : list of color
        The colors to use.
    orientation : {'horizontal', 'vertical'}
    ax : Axes, optional
        The Axes to draw into. If None, the current axes will be used.
    **kwargs
        All other keyword arguments are passed to plt.text(), so you can
        set the font size, family, etc.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transData
    fig = ax.figure
    canvas = fig.canvas

    assert orientation in ['horizontal', 'vertical']
    if orientation == 'vertical':
        kwargs.update(rotation=90, verticalalignment='bottom')

    for s, c in zip(strings, colors):
        text = ax.text(x, y, s + "", color=c, transform=t, **kwargs)
        text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        # Convert window extent from pixels to inches
        # to avoid issues displaying at different dpi
        ex = fig.dpi_scale_trans.inverted().transform_bbox(ex)

        if orientation == 'horizontal':
            t = text.get_transform() + \
                offset_copy(Affine2D(), fig=fig, x=ex.width, y=0)
        else:
            t = text.get_transform() + \
                offset_copy(Affine2D(), fig=fig, x=0, y=ex.height)
                
def get_quiver_arrays(x_pos_min, y_pos_min, image_array, velx_full, vely_full, no_of_quivers=31., smooth_cells=None, center_vel=None, velz_full=None, axis = None):
    if axis == 'xy':
        center_vel_plane = np.array([center_vel[0], center_vel[1]])
        center_vel_perp = center_vel[2]
    elif axis == 'xz':
        center_vel_plane = np.array([center_vel[0], center_vel[2]])
        center_vel_perp = center_vel[1]
    elif axis == 'yz':
        center_vel_plane = np.array([center_vel[1], center_vel[2]])
        center_vel_perp = center_vel[0]
    else:
        center_vel_plane = center_vel
        center_vel_perp = 0
    annotate_freq = float(np.shape(image_array)[0])/float(no_of_quivers)
    if smooth_cells is None:
        smoothing_val = int(annotate_freq/2)
    else:
        smoothing_val = int(smooth_cells)
    x_ind = []
    y_ind = []
    counter = 0
    while counter < (no_of_quivers):
        valx = int(x_pos_min + annotate_freq*counter + annotate_freq/2.)
        valy = int(y_pos_min + annotate_freq*counter + annotate_freq/2.)
        x_ind.append(int(valx))
        y_ind.append(int(valy))
        counter = counter + 1
    
    try:
        if velz_full is None:
            velz_full = np.zeros(np.shape(velx_full))
    except:
        velz_full = velz_full
    
    velx = []
    vely = []
    velz = []
    for x in x_ind:
        x = int(x)
        xarr = []
        yarr = []
        zarr = []
        for y in y_ind:
            y = int(y)
            x_vel = []
            y_vel = []
            z_vel = []
            for curx in range(x-smoothing_val, x+smoothing_val):
                for cury in range(y-smoothing_val, y+smoothing_val):
                    x_vel.append(velx_full[curx][cury])
                    y_vel.append(vely_full[curx][cury])
                    z_vel.append(velz_full[curx][cury])
            x_vel = np.mean(x_vel)
            y_vel = np.mean(y_vel)
            z_vel = np.mean(z_vel)
            xarr.append(x_vel)
            yarr.append(y_vel)
            zarr.append(z_vel)
        velx.append(xarr)
        vely.append(yarr)
        velz.append(zarr)
    velx = np.array(velx)
    vely = np.array(vely)
    velz = np.array(velz)
    try:
        velx = velx - center_vel_plane[0]
        vely = vely - center_vel_plane[1]
        velz = velz - center_vel_perp
    except:
        velx = velx
        vely = vely
        velz = velz
    return velx, vely, velz

def my_own_quiver_function(axis, X_pos, Y_pos, X_val, Y_val, plot_velocity_legend='False', standard_vel=5, limits=None, Z_val=None, width_ceil = 0.8, zorder=3, renderer=None):
    global fontsize_global
    if plot_velocity_legend == 'False' or plot_velocity_legend == False:
        plot_velocity_legend = False
    elif plot_velocity_legend == 'True' or plot_velocity_legend == True:
        plot_velocity_legend = True
        if standard_vel > 1.0:
            legend_text=str(int(standard_vel)) + "km$\,$s$^{-1}$"
        else:
            legend_text=str(int(standard_vel*10.)/10.) + "km$\,$s$^{-1}$"
    standard_vel = yt.YTQuantity(1, 'km').in_units('cm').value * standard_vel
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
    max_length = 0
    
    for xp in range(len(X_pos[0])):
        for yp in range(len(Y_pos[0])):
            xvel = len_scale*(X_val[xp][yp]/standard_vel)
            yvel = len_scale*(Y_val[xp][yp]/standard_vel)
            length = np.sqrt(X_val[xp][yp]**2 + Y_val[xp][yp]**2)
            if length > max_length:
                max_length = length
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
            #import pdb
            #pdb.set_trace()
            axis.add_patch(mpatches.FancyArrowPatch((X_pos[xp][yp], Y_pos[xp][yp]), (X_pos[xp][yp]+xvel, Y_pos[xp][yp]+yvel), color=color, linewidth=width_val, arrowstyle='->', mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0)) #alpha=width_val/width_ceil
    print("Max velocity =", max_length)
    if plot_velocity_legend:
        #print("plotting quiver legend")
        pos_start = [xmax - 0.2*(xmax-xmin), ymin + 0.12*(ymax-ymin)]
        #pos_start = [xmax - 0.25*(xmax-xmin), ymin + (fontsize_global/100)*0.70*(ymax-ymin)]
        xvel = len_scale*(standard_vel/standard_vel)
        yvel = 0.0
        width_val = width_ceil
        #annotate_text = axis.text((xmax - 0.01*(xmax-xmin)), (ymin + 0.05*(ymax-ymin)), legend_text, va="center", ha="right", color='w', fontsize=fontsize_global)
        annotate_text = axis.text(xmax - 0.01*xmax, ymin+0.01*ymax, legend_text, va="bottom", ha="right", color='w', fontsize=fontsize_global)
        annotate_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        pf = plt.gcf()
        renderer = pf.canvas.get_renderer()
        bbox_text = annotate_text.get_window_extent(renderer=renderer)
        bbox_text = annotate_text.get_window_extent(renderer=renderer)
        bbox_text = annotate_text.get_window_extent(renderer=renderer)
        text_height = bbox_text.height
        #set arrow pos to be 10% higher than text pos
        
        import pdb
        pdb.set_trace()
        axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=width_val, edgecolor = 'k', mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0))
        #axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=width_val, edgecolor = 'k', mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0))
    return axis

def annotate_particles(axis, particle_position, accretion_rad, limits, annotate_field=None, field_symbol="M", units=None, particle_tags=None, lw=1.5, zorder=4, split_threshold = 4):
    global fontsize_global
    if annotate_field is not None and units is not None:
        annotate_field = annotate_field.in_units(units)
    #part_color = ['cyan','magenta','r','b','y','w','k']
    part_color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    box_size = xmax - xmin
    if accretion_rad/box_size < 0.05:
        line_rad = yt.YTQuantity(0.005*box_size, 'AU')
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
    for pos_it in range(len(particle_tags)): #np.argsort(particle_tags):
        #axis.scatter(particle_position[0][pos_it], particle_position[1][pos_it], c=part_color[pos_it], s=1, zorder=zorder)
        axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=lw, c='k')
        axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=lw, c='k')
        axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=lw/2, c=part_color[pos_it])
        axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=lw/2, c=part_color[pos_it])
        circle = mpatches.Circle([particle_position[0][pos_it], particle_position[1][pos_it]], accretion_rad, fill=False, lw=lw/2, edgecolor='k')
        axis.add_patch(circle)
        
        if annotate_field is not None:
            if units is not None:
                annotate_field = annotate_field.in_units(units)
                if units == "Msun":
                    unit_string = "$\,$M$_\odot$"
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
                p_t = "M=[ "
                rainbow_text_colors.append('white')
            p_t = p_t + P_msun + ' '
            rainbow_text_colors.append(part_color[pos_it])
            if pos_it != len(particle_tags)-1:
                p_t = p_t + ', '
            else:
                p_t = p_t + ']$\,$M$_\odot$'
            rainbow_text_colors.append('white')
        '''
        if annotate_field is not None:
            if units is not None:
                annotate_field = annotate_field.in_units(units)
                if units == "Msun":
                    unit_string = "$\,$M$_\odot$"
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
            '''
    if annotate_field is not None:
        if len(particle_tags) > split_threshold:
            string_1 = p_t[:31]
            string_2 = p_t[32:]
            colors_1 = rainbow_text_colors[:9]
            colors_2 = rainbow_text_colors[9:]
            rainbow_text((xmin + 0.01*(box_size)), (ymin + 0.029*(ymax-ymin)*4), string_1.split(' '), colors_1, size=fontsize_global, zorder=10, ax=axis)
            rainbow_text((xmin + 0.1*(box_size)), (ymin + 0.029*(ymax-ymin)), string_2.split(' '), colors_2, size=fontsize_global, zorder=10, ax=axis)
        else:
            #try:
            rainbow_text((xmin + 0.01*(box_size)), (ymin + 0.029*(ymax-ymin)), p_t.split(' '), rainbow_text_colors, size=fontsize_global, zorder=10, ax=axis)
            #except:
            #    print("couldn't annotate particle masses")
    return axis
