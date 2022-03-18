#!/usr/bin/env python
import yt
from yt.utilities.exceptions import YTOutputNotIdentified
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import my_ramses_fields as myf
import scipy.spatial as spatial
import pickle
import matplotlib.patheffects as path_effects
from matplotlib import transforms
import csv
import math

fontsize_global=12

units = {"length_unit":(1,"cm"), "mass_unit":(1,"g"), "velocity_unit":(1, "cm/s"), "time_unit":(1, "s"), "density_unit":(1, "g/cm**3")}
length_unit = yt.YTQuantity(units["length_unit"][0],units["length_unit"][1])
velocity_unit = yt.YTQuantity(units["velocity_unit"][0],units["velocity_unit"][1])
mass_unit = yt.YTQuantity(units["mass_unit"][0],units["mass_unit"][1])
time_unit = yt.YTQuantity(units["time_unit"][0],units["time_unit"][1])
density_unit = yt.YTQuantity(units["density_unit"][0],units["density_unit"][1])

def set_units(input_units):
    """
    sets units used for all calculations
    """
    global length_unit
    global velocity_unit
    global mass_unit
    global time_unit
    for key in input_units.keys():
        if key in units.keys():
            units[key] = input_units[key]
    length_unit = yt.YTQuantity(units["length_unit"][0],units["length_unit"][1])
    velocity_unit = yt.YTQuantity(units["velocity_unit"][0],units["velocity_unit"][1])
    mass_unit = yt.YTQuantity(units["mass_unit"][0],units["mass_unit"][1])
    time_unit = yt.YTQuantity(units["time_unit"][0],units["time_unit"][1])
    density_unit = yt.YTQuantity(units["density_unit"][0],units["density_unit"][1])

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
            t = transforms.offset_copy(text._transform, x=0.75*ex.width, units='dots')
        else:
            #import pdb
            #pdb.set_trace()
            #t = transforms.offset_copy(text._transform, x=space_size, units='dots')
            t = transforms.offset_copy(text._transform, x=0.75*ex.width, units='dots')
        
def set_global_font_size(x):
    global fontsize_global
    fontsize_global = x
    return fontsize_global

def get_global_font_size():
    global fontsize_global
    return fontsize_global

def find_sink_formation_time(files,sink_number=164):
    global time_unit
    csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
    part_file = files[0].split("info")[0]+"stars_output.snktxt"
    with open(part_file, 'r') as f:
        reader = csv.reader(f, dialect='dat')
        for row in reader:
            try:
                if int(row[0])==sink_number:
                    sink_form = float(row[10])*time_unit.in_units('yr').value
            except:
                continue
    return sink_form

def generate_frame_times(files, dt, start_time=0, presink_frames=25, end_time=None, form_time=None):
    if form_time != None:
        sink_form_time = form_time
    else:
        sink_form_time = find_sink_formation_time(files)
        
    if end_time != None:
        max_time = end_time
    else:
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        with open(files[-1].split("info")[0]+'stars_output.snktxt', 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                time = yt.YTQuantity(float(row[1])*time_unit.in_units('yr').value, 'yr')
                break
        if form_time != None:
            max_time = time.in_units('yr') - form_time
        else:
            max_time = time.in_units('yr')

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

def find_files(m_times, files, sink_form_time, sink_number, verbatim=True):
    global time_unit
    csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
    usable_files = []
    mit = 0
    min = 0
    max = len(files)-1
    pit = 0
    while mit < len(m_times):
        #define it as half way between the min and max indexes
        it = int(np.round(min + ((max - min)/2.)))
        #Check how close it is to your value:
        with open(files[it].split("info")[0]+'stars_output.snktxt', 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                time = float(row[1])*time_unit.in_units('yr').value - sink_form_time.value
                break
        f.close()
        if verbatim == True:
            print("Current file time =", time, "for interator =", it)
        #check if it is honed in on the closes file, or times match:
        if pit == it or time == m_times[mit]:
            #Checking which file is actually closest to the movie time
            potential_files = files[it-1:it+2]
            diff_arr = []
            for pot_file in potential_files:
                with open(pot_file.split("info")[0]+'stars_output.snktxt', 'r') as f:
                    reader = csv.reader(f, dialect='dat')
                    for row in reader:
                        time = float(row[1])*time_unit.in_units('yr').value - sink_form_time.value
                        break
                f.close()
                diff_val = time - m_times[mit]
                diff_arr.append(diff_val)
            try:
                append_it = np.argmin(abs(np.array(diff_arr)))
            
                if m_times[mit] == 0 and diff_arr[append_it] < 0:
                    append_it = append_it + 1
                usable_files.append(potential_files[append_it])
            except:
                usable_files.append(files[it])
            if verbatim == True:
                print("found time", time, "for m_time", m_times[mit], "with file:", usable_files[-1])
            mit = mit + 1
            min = it
            max = len(files)-1
            pit = 0
        elif time < m_times[mit]:
            min = it
            pit = it
        elif time > m_times[mit]:
            max = it
            pit = it
    return usable_files

def get_particle_data(ds, axis='xy', sink_id=None, region=None):
    """
    Retrieve particle data for plotting. NOTE: CANNOT RETURN PARTICLE VELOCITIES AS THESES ARE NOT STORED IN THE MOVIE FILES.
    """
    dd = ds.all_data()
    if sink_id != None:
        sink_id = sink_id
    else:
        sink_id = np.argmin(dd['sink_particle_speed'])
    myf.set_centred_sink_id(sink_id)
    
    #TODO: Possibly update to set the usuable sink particle IDs to custom values (f.x. [80, 90])
    if region != None:
        usable_sinks = np.argwhere((dd['sink_particle_posx'].in_units('au')>region.left_edge[0])&(dd['sink_particle_posx'].in_units('au')<region.right_edge[0])&(dd['sink_particle_posy'].in_units('au')>region.left_edge[1])&(dd['sink_particle_posy'].in_units('au')<region.right_edge[1])&(dd['sink_particle_posz'].in_units('au')>region.left_edge[2])&(dd['sink_particle_posz'].in_units('au')<region.right_edge[2])).T[0]
    else:
        usable_sinks = np.arange(sink_id, len(dd['sink_particle_tag']))
    
    part_tags = dd['sink_particle_tag'][usable_sinks]
    part_mass = dd['sink_particle_mass'][usable_sinks].in_units('msun').value
    
    accretion_rad = 4.*np.min(dd['dx']).in_units('au').value
    #center_pos = myf.get_center_pos()
    if axis == 'xy':
        part_pos_x = dd['sink_particle_posx'][usable_sinks].in_units('au').value
        part_pos_y = dd['sink_particle_posy'][usable_sinks].in_units('au').value
    elif axis == 'xz':
        part_pos_x = dd['sink_particle_posx'][usable_sinks].in_units('au').value
        part_pos_y = dd['sink_particle_posz'][usable_sinks].in_units('au').value
    elif axis == 'yz':
        part_pos_x = dd['sink_particle_posy'][usable_sinks].in_units('au').value
        part_pos_y = dd['sink_particle_posz'][usable_sinks].in_units('au').value
    positions = np.array([part_pos_x,part_pos_y])
    part_info = {'particle_mass':part_mass,
                 'particle_position':positions,
                 'accretion_rad':accretion_rad,
                 'particle_tag':part_tags}
    return part_info

def initialise_grid(file, zoom_times=0, num_of_vectors=31.):#, center=0):
    f = h5py.File(file, 'r')
    for key in list(f.keys()):
        if 'dens' in key:
            field = key
    dim = np.shape(f[field])[0]
    if zoom_times != 0.0:
        zoom_cell = int((dim - dim/float(zoom_times))/2.)
    else:
        zoom_cell = 0
    xmin = f['minmax_xyz'][0][0]/yt.units.au.in_units('cm').value
    xmax = f['minmax_xyz'][0][1]/yt.units.au.in_units('cm').value
    cl = (xmax-xmin)/dim
    xmin = f['minmax_xyz'][0][0]/yt.units.au.in_units('cm').value + zoom_cell*cl
    xmax = f['minmax_xyz'][0][1]/yt.units.au.in_units('cm').value - zoom_cell*cl
    annotate_freq = ((xmax/cl) - (xmin/cl))/num_of_vectors
    x = np.arange(xmin+cl/2., xmax+cl/2., cl)
    #x = np.arange(xmin, xmax, cl)[:-1]
    '''
    if (len(x) % 2 != 0):
        x = np.arange(int(xmin), int(xmax), cl)
    '''
    x_ind = []
    y_ind = []
    counter = 0
    while counter < num_of_vectors:
        val = annotate_freq*counter + annotate_freq/2.
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    x_vel = []
    y_vel = []
    for x_counter in x_ind:
        val = xmin+cl/2. + x_counter*cl
        #val = xmin + x_counter*cl
        x_vel.append(val)
        y_vel.append(val)
    x_vel = np.array(x_vel)
    y_vel = np.array(y_vel)
    X, Y = np.meshgrid(x, x)
    X_vel, Y_vel = np.meshgrid(x_vel, y_vel)
    '''
    if center != 0:
        x_pos = np.round(f['particlepositions'][0][center-1]/yt.units.au.in_units('cm').value/cl)*cl
        y_pos = np.round(f['particlepositions'][1][center-1]/yt.units.au.in_units('cm').value/cl)*cl
        X = X + x_pos
        Y = Y + y_pos
        X_vel = X_vel + x_pos
        Y_vel = Y_vel + y_pos
    '''
    #print "created meshs"
    return X, Y, X_vel, Y_vel, cl

def get_quiver_arrays(x_pos_min, y_pos_min, image_array, velx_full, vely_full, no_of_quivers=32., smooth_cells=None, center_vel=None, velz_full=None, axis = None):
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
    annotate_freq = float(np.shape(image_array)[0])/float(no_of_quivers-1)
    if smooth_cells is None:
        smoothing_val = int(annotate_freq/2)
    else:
        smoothing_val = int(smooth_cells)
    x_ind = []
    y_ind = []
    counter = 0
    while counter < (no_of_quivers-1):
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
    rv_colors = np.linspace(-1, 1, 256)
    rv_cmap = plt.cm.get_cmap('bwr')
    
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
                    #color = 'w'
                    color = 'k'
            except:
                #cmap = 'idl06_r'
                zvel = Z_val[xp][yp]/len_scale
                cit = np.argmin(abs(rv_colors - zvel))
                color = rv_cmap(cit)
            axis.add_patch(mpatches.FancyArrowPatch((X_pos[xp][yp], Y_pos[xp][yp]), (X_pos[xp][yp]+xvel, Y_pos[xp][yp]+yvel), color=color, linewidth=width_val, arrowstyle='->', mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0, rasterized=False))#, alpha=width_val/width_ceil))
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
        axis.scatter(particle_position[0][pos_it], particle_position[1][pos_it], c=part_color[pos_it], s=5, zorder=zorder)
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

def profile_plot(x, y, weight=None, n_bins=None, log=False, bin_data=None, bin_min=None, bin_max=None, calc_vel_dispersion=False, cumulative=False):
    unit_string = str(x.unit_quantity).split(' ')[-1]
    if bin_data is None:
        bin_data = x
    if bin_min is None:
        bin_min = np.min(bin_data)
    if bin_max is None:
        try:
            bin_max = np.max(bin_data.value)
        except:
            bin_max = np.max(bin_data)
    if n_bins is None:
        unit_string = str(bin_data.unit_quantity).split(' ')[-1]
        bins = [0] + list(set(bin_data.value))
        bins = np.sort(bins)
        #bins = bins[::]
        bins = np.append(bins, bins[-1] + (bins[-1]-bins[-2]))
        bins = yt.YTArray(bins, unit_string)
    elif log == True:
        if bin_min is None:
            bins = np.logspace(np.log10(np.min(bin_data)), np.log10(np.max(x)), n_bins+1)
        else:
            bins = np.logspace(np.log10(bin_min), np.log10(np.max(x)), n_bins+1)
    else:
        bins = np.linspace(bin_min, bin_max, n_bins+1)
    print(("No of bins:", len(bins)))
    neg_inds = np.where(bins < 0)[0]
    #bins[neg_inds[-1]] = -1*(bins[neg_inds[-1]+1] - yt.YTArray(0.01, unit_string))

    x_array = []
    y_array = []
    
    prev_bin = bins[0]
    cum_val = yt.YTQuantity(0.0, str(y.unit_quantity).split(' ')[-1])
    for bin in bins[1:]:
        ind = np.where((x >= prev_bin) & (x < bin))[0]
        mid_x = (bin+prev_bin)/2.
        x_array.append(mid_x)
        if len(ind) != 0:
            if calc_vel_dispersion:
                mean_v_phi = np.mean(y[ind])
                y[ind] = y[ind] - mean_v_phi
                bin_val = np.std(y[ind].in_units('km/s'))
            elif weight is not None:
                bin_val = (np.sum(y[ind]*weight[ind]))/np.sum(weight[ind])
            else:
                bin_val = np.mean(y[ind])
        else:
            bin_val = np.nan
            print("FOUND EMPTY BIN")
        prev_bin = bin
        if cumulative:
            cum_val = cum_val + np.sum(y[ind])
            bin_val = cum_val
        y_array.append(bin_val)
        print("Value =", bin_val, ", at radius=", mid_x)

    x_array = np.array(x_array)
    y_array = np.array(y_array)
    inds = np.where(np.isnan(y_array)==False)
    x_array = x_array[inds]
    y_array = y_array[inds]

    return x_array, y_array

def sample_points(x_field, y_field, z, bin_no=2., no_of_points=2000, weight_arr=None):
    big_array = np.array([z, x_field, y_field])
    plot_array = []
    z_bins = np.linspace(np.min(np.abs(z)), np.max(np.abs(z)), bin_no+1)
    z_bins = z_bins[::-1]
    import random
    prev_z = z_bins[0]
    counter_max = no_of_points/bin_no
    for z_bin in z_bins[1:]:
        filtered_z = (z < prev_z)*(z > z_bin)
        temp_array = filtered_z*big_array
        temp_array = temp_array.transpose()
        counter = 0
        while counter < counter_max:
            data_point = random.choice(temp_array)
            if 0.0 not in data_point:
                plot_array.append(data_point)
                counter = counter + 1
                #print "selected random data point", data_point
        prev_z = z_bin
    #print "selected random data points"
    plot_array = np.array(plot_array)
    plot_array = plot_array.transpose()
    return plot_array

def sliceplot(ds, X, Y, field, cmap=plt.cm.get_cmap('brg'), log=False, resolution=1024, center=0, units=None, cbar_label="Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)", weight=None, vmin=None, vmax=None): #cmap=plt.cm.get_cmap('bwr_r')
    """
       Creates a slice plot along the xy-plane for a data cube. This is done by interpolating  the simulation output onto a grid.
    """
    #print "SLICE PLOT CENTER =", center
    #print("image resolution =", resolution)
    x = np.linspace(np.min(X), np.max(X), resolution)
    y = np.linspace(np.min(Y), np.max(Y), resolution)
    #print("len(x) =", len(x))
    X = yt.YTArray(X, 'AU')
    Y = yt.YTArray(Y, 'AU')
    xy = np.meshgrid(x,y)
    dd = ds.all_data()
    no_of_particles = len(dd['particle_mass'])
    cell_max = np.max(dd['dz'].in_units('AU'))
    dd = ds.region([0.0,0.0,0.0], [np.min(X), np.min(Y), -cell_max], [np.max(X), np.max(Y), cell_max])
    xyz = yt.YTArray([dd['x'].in_units('AU'),dd['y'].in_units('AU'),dd['z'].in_units('AU')]).T
    print("Computing tree...")
    tree = spatial.cKDTree(xyz)
    del xyz
    print("Done computing tree...")
    #create grid of the pixel positions and quiver positions
    xyz_grid = yt.YTArray([xy[0].flatten(), xy[1].flatten(), np.zeros_like(xy[0].flatten())]).T

    #This line also takes a while
    #print("Finding nearest points...")
    nearest_points = tree.query(xyz_grid, k=2, distance_upper_bound=cell_max, n_jobs=-1)[1]
    
    #print("len(nearest_points) =", len(nearest_points))
    del xyz_grid
    #print("Done finding nearest points...")

    if center == 3:
        field_grid = np.zeros_like(xy[0])
        for cen in range(len(dd['particle_mass'])+1):
            myf.set_center(cen)
            print(("Doing center =", cen))
            del dd
            dd = ds.region([0.0,0.0,0.0], [np.min(X), np.min(Y), -cell_max], [np.max(X), np.max(Y), cell_max])
            field_grid = field_grid + (dd[field][nearest_points[:,0]].reshape(xy[0].shape) + dd[field][nearest_points[:,1]].reshape(xy[0].shape))/2.
        field_grid = field_grid/(len(dd['particle_mass'])+1)
    else:
        myf.set_center(center)
        del dd
        dd = ds.region([0.0,0.0,0.0], [np.min(X), np.min(Y), -cell_max], [np.max(X), np.max(Y), cell_max])
        field_grid = (dd[field][nearest_points[:,0]].reshape(xy[0].shape) + dd[field][nearest_points[:,1]].reshape(xy[0].shape))/2.

    if weight is not None:
        weight_field = (dd[weight][nearest_points[:,0]].reshape(xy[0].shape) + dd[weight][nearest_points[:,1]].reshape(xy[0].shape))/2.
    else:
        weight_field = np.ones(np.shape(field_grid))
    weight_field = weight_field/np.max(weight_field)

    if units is not None:
        field_grid = field_grid.in_units(units)
    field_grid = field_grid

    if field == 'Relative_Keplerian_Velocity':
        vmin=0.0
        vmax=2.0
    
    plt.clf()
    fig, ax = plt.subplots()
    #print("len(xy[0]) =", len(xy[0]))
    if log:
        if vmin is not None:
            plot = ax.pcolormesh(xy[0], xy[1], field_grid.value, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
            for i,j in zip(plot.get_facecolors(),weight_field.flatten()):
                i[3] = j
        else:
            plot = ax.pcolormesh(xy[0], xy[1], field_grid.value, cmap=cmap, norm=LogNorm(), rasterized=True)
            for i,j in zip(plot.get_facecolors(),weight_field.flatten()):
                i[3] = j
    else:
        if vmin is not None:
            plot = ax.pcolormesh(xy[0], xy[1], field_grid.value, cmap=cmap, rasterized=True, vmin=vmin, vmax=vmax)
            for i,j in zip(plot.get_facecolors(),weight_field.flatten()):
                i[3] = j
        else:
            plot = ax.pcolormesh(xy[0], xy[1], field_grid.value, cmap=cmap, rasterized=True)
            for i,j in zip(plot.get_facecolors(),weight_field.flatten()):
                i[3] = j
    plt.gca().set_aspect('equal')
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label(cbar_label, rotation=270, labelpad=13, size=14)

    ax.set_xlim([np.min(X), np.max(X)])
    ax.set_ylim([np.min(Y), np.max(Y)])
    ax.set_xlabel('$x$ (AU)', labelpad=-1)
    ax.set_ylabel('$y$ (AU)', labelpad=-20)
    return fig, ax, xy, field_grid, weight_field
