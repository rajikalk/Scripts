#!/usr/bin/env python
import yt
from yt.utilities.exceptions import YTOutputNotIdentified
import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import my_fields as myf
import scipy.spatial as spatial
import pickle
import matplotlib.patheffects as path_effects

fontgize_global=12

def set_global_font_size(x):
    global fontgize_global
    fontgize_global = x
    return fontgize_global

def get_global_font_size():
    global fontgize_global
    return fontgize_global

def yt_file_type(file):
    try:
        part_file=file[:-13] + 'part' + file[-6:]
        ds = yt.load(file, particle_filename=part_file)
        yt_file = True
    except YTOutputNotIdentified:
        yt_file = False
    return yt_file

def find_sink_formation_time(files):
    try:
        file = files[-2]
        part_file=file[:-13] + 'part' + file[-6:]
        '''
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        if ('io', u'particle_creation_time') in ds.field_list:
            sink_form = float(np.min(dd['particle_creation_time']/yt.units.yr.in_units('s')).value)
        else:
            sink_form = 0.0
        '''
        f = h5py.File(part_file, 'r')
        if f[list(f.keys())[1]][-1][-1] > 0:
            sink_form = np.min(f[list(f.keys())[11]][:,5])/yt.units.yr.in_units('s').value
        else:
            sink_form = 0.0
    except:
        sink_form = None
        for source in range(len(files)):
            f = h5py.File(files[source], 'r')
            if 'particlemasses' in list(f.keys()):
                sink_form = f['time'][0]/yt.units.yr.in_units('s').value
                break
        if sink_form is None:
            sink_form = -1*f['time'][0]/yt.units.yr.in_units('s').value
    return sink_form

def get_image_mesh(file, zoom_times):
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
    x = np.arange(xmin, xmax, cl)
    X, Y = np.meshgrid(x, x)
    return X, Y

def generate_frame_times(files, dt, start_time=0, presink_frames=25, end_time=2000.):
    try:
        file = files[-1]
        part_file=file[:-13] + 'part' + file[-6:]
        #ds = yt.load(file, particle_filename=part_file)
        #dd = ds.all_data()
        f = h5py.File(part_file, 'r')
        sink_form_time = find_sink_formation_time(files)
        max_time = f[list(f.keys())[7]][0][-1]/yt.units.yr.in_units('s').value - sink_form_time
    except:
        f = h5py.File(files[-1], 'r')
        sink_form_time = find_sink_formation_time(files)
        max_time = f['time'][0]/yt.units.yr.in_units('s').value - sink_form_time
        f.close()
    #if end_time is not None and end_time < max_time:
    max_time = end_time

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
        try:
            file = files[it]
            part_file=file[:-13] + 'part' + file[-6:]
            f = h5py.File(part_file, 'r')
            time = f['real scalars'][0][1]/yt.units.year.in_units('s').value-sink_form_time
        except:
            f = h5py.File(files[it], 'r')
            time = f['time'][0]/yt.units.yr.in_units('s').value-sink_form_time
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
                try:
                    part_file=pfile[:-13] + 'part' + pfile[-6:]
                    f = h5py.File(part_file, 'r')
                    time = f['real scalars'][0][1]/yt.units.year.in_units('s').value-sink_form_time
                except:
                    f = h5py.File(pfile, 'r')
                    time = f['time'][0]/yt.units.yr.in_units('s').value-sink_form_time
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
            usable_files.append(append_file)
            try:
                part_file=append_file[:-13] + 'part' + append_file[-6:]
                f = h5py.File(part_file, 'r')
                time = f['real scalars'][0][1]/yt.units.year.in_units('s').value-sink_form_time
            except:
                f = h5py.File(append_file, 'r')
                time = f['time'][0]/yt.units.yr.in_units('s').value-sink_form_time
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

def get_particle_data(file, axis='xz', proj_or=None):
    """
    Retrieve particle data for plotting. NOTE: CANNOT RETURN PARTICLE VELOCITIES AS THESES ARE NOT STORED IN THE MOVIE FILES.
    """
    part_mass = []
    part_pos_x = []
    part_pos_y = []
    accretion_rad = []
    #center_pos = myf.get_center_pos()
    try:
        part_file=file[:-13] + 'part' + file[-6:]
        f = h5py.File(part_file, 'r')
        tag_ind = np.where(f['particle names'][:] == b'tag                     ')[0][0]
        ordered_inds = np.argsort(f['tracer particles'][:,tag_ind])
        mass_ind = np.where(f['particle names'][:] == b'mass                    ')[0][0]
        part_mass = f['tracer particles'][:,mass_ind][ordered_inds]/yt.units.msun.in_units('g').value
        posx_ind = np.where(f['particle names'][:] == b'posx                    ')[0][0]
        posy_ind = np.where(f['particle names'][:] == b'posy                    ')[0][0]
        if axis == 'xy':
            part_pos_x = f['tracer particles'][:,posx_ind][ordered_inds]/yt.units.au.in_units('cm').value# - center_pos[0]/yt.units.au.in_units('cm').value
            part_pos_y = f['tracer particles'][:,posy_ind][ordered_inds]/yt.units.au.in_units('cm').value# - center_pos[1]/yt.units.au.in_units('cm').value
            depth_pos = list(range(len(part_mass)))
        else:
            if proj_or is not None:
                L = np.array([proj_or[0],proj_or[1]])
                L_orth = np.array([[proj_or[1]], [-1*proj_or[0]]])
                L_len = np.sqrt(L_orth[0]**2. + L_orth[1]**2.)
                part_pos_x = f['tracer particles'][:,posx_ind][ordered_inds]/yt.units.au.in_units('cm').value
                part_pos_y = f['tracer particles'][:,posy_ind][ordered_inds]/yt.units.au.in_units('cm').value
                r = np.array([part_pos_x, part_pos_y])
                part_pos_x = -1*np.dot(r.T,(L_orth/L_len))
                part_pos_x = part_pos_x.T[0]
                depth_pos = -1*np.dot(r.T,(L/L_len))
                depth_pos = np.argsort(depth_pos)[::-1]
            else:
                part_pos_x = f['tracer particles'][:,posx_ind][ordered_inds]/yt.units.au.in_units('cm').value
                depth_pos = f['tracer particles'][:,posy_ind]/yt.units.au.in_units('cm').value
                depth_pos = depth_pos[::-1]
            posz_ind = np.where(f['particle names'][:] == b'posz                    ')[0][0]
            part_pos_y = f['tracer particles'][:,posz_ind][ordered_inds]/yt.units.au.in_units('cm').value
        radius_ind = [i for i, v in enumerate(f['real runtime parameters'][:]) if v[0] == b'sink_accretion_radius                                                           '][0]
        accretion_rad = f['real runtime parameters'][radius_ind][1]/yt.units.au.in_units('cm').value
    except:
        f = h5py.File(file, 'r')
        part_mass = np.array(f["particlemasses"])/yt.units.msun.in_units('g').value
        ordered_inds = np.argsort(part_mass)[::-1]
        part_mass = np.array(f["particlemasses"][:][ordered_inds])/yt.units.msun.in_units('g').value
        if axis == 'xy':
            part_pos_x = f["particlepositions"][0][ordered_inds]/yt.units.au.in_units('cm').value
            part_pos_y = f["particlepositions"][1][ordered_inds]/yt.units.au.in_units('cm').value
            depth_pos = list(range(len(part_mass)))
        else:
            part_pos_x = f["particlepositions"][0][ordered_inds]/yt.units.au.in_units('cm').value
            depth_pos = f["particlepositions"][1][ordered_inds]/yt.units.au.in_units('cm').value
            depth_pos = depth_pos[::-1]
            part_pos_y = f["particlepositions"][2][ordered_inds]/yt.units.au.in_units('cm').value
        accretion_rad = f['r_accretion'][0]/yt.units.au.in_units('cm').value
    positions = np.array([part_pos_x,part_pos_y])
    part_info = {'particle_mass':part_mass,
                 'particle_position':positions,
                 'accretion_rad':accretion_rad,
                 'depth_position':depth_pos}
    return part_info

def initialise_grid(file, zoom_times=0):#, center=0):
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
    annotate_freq = ((xmax/cl) - (xmin/cl))/31.
    x = np.arange(xmin+cl/2., xmax+cl/2., cl)
    #x = np.arange(xmin, xmax, cl)[:-1]
    '''
    if (len(x) % 2 != 0):
        x = np.arange(int(xmin), int(xmax), cl)
    '''
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 31:
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

def get_quiver_arrays(x_pos_min, y_pos_min, image_array, velx_full, vely_full, no_of_quivers=32., smooth_cells=None, center_vel=None):
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
    velx = []
    vely = []
    for x in x_ind:
        x = int(x)
        xarr = []
        yarr = []
        for y in y_ind:
            y = int(y)
            x_vel =[]
            y_vel = []
            for curx in range(x-smoothing_val, x+smoothing_val):
                for cury in range(y-smoothing_val, y+smoothing_val):
                    x_vel.append(velx_full[curx][cury])
                    y_vel.append(vely_full[curx][cury])
            x_vel = np.mean(x_vel)
            y_vel = np.mean(y_vel)
            xarr.append(x_vel)
            yarr.append(y_vel)
        velx.append(xarr)
        vely.append(yarr)
    velx = np.array(velx)
    vely = np.array(vely)
    try:
        velx = velx - center_vel[0]
        vely = vely - center_vel[1]
    except:
        velx = velx
        vely = vely
    return velx, vely

def my_own_quiver_function(axis, X_pos, Y_pos, X_val, Y_val, plot_velocity_legend='False', standard_vel=5, limits=None):
    legend_text=str(int(standard_vel)) + "kms$^{-1}$"
    global fontgize_global
    if plot_velocity_legend == 'False':
        plot_velocity_legend = False
    elif plot_velocity_legend == 'True':
        plot_velocity_legend = True
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
    len_scale = standard_vel/(0.07*(xmax - xmin))
    vels = np.hypot(X_val, Y_val)
    for xp in range(len(X_pos[0])):
        for yp in range(len(Y_pos[0])):
            xvel = X_val[xp][yp]/len_scale
            yvel = Y_val[xp][yp]/len_scale
            width_val = np.sqrt(X_val[xp][yp]**2. + Y_val[xp][yp]**2.)/standard_vel
            if width_val > 0.8:
                width_val = 0.8
            axis.add_patch(mpatches.FancyArrowPatch((X_pos[xp][yp], Y_pos[xp][yp]), (X_pos[xp][yp]+xvel, Y_pos[xp][yp]+yvel), color='w', linewidth=1.*width_val, arrowstyle='->', mutation_scale=15.*width_val, shrinkA=0.0, shrinkB=0.0))
    if plot_velocity_legend:
        #print("plotting quiver legend")
        pos_start = [xmax - 0.15*(xmax-xmin), ymin + 0.07*(ymax-ymin)]
        xvel = standard_vel/len_scale
        yvel = 0.0
        width_val = 1.0
        axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=1.*width_val, mutation_scale=15.*width_val))
        annotate_text = axis.text((xmax - 0.01*(xmax-xmin)), (ymin + 0.03*(ymax-ymin)), legend_text, va="center", ha="right", color='w', fontsize=fontgize_global)
        annotate_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    return axis

def annotate_particles(axis, particle_position, accretion_rad, limits, annotate_field=None, field_symbol='M', units=None,depth_array=None):
    global fontgize_global
    if depth_array is None:
        depth_array = np.arange(len(particle_position[0]))
    if np.max(depth_array) > len(depth_array):
        depth_array = np.argsort(depth_array)
    if annotate_field is not None and units is not None:
        annotate_field = annotate_field.in_units(units)
    part_color = ['cyan','cyan','r','c','y','w','k']
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
        unit_string = ''
        for s in s1:
            if s == s1[0]:
                unit_string = unit_string + s
            else:
                unit_string = unit_string + "$^" + s[0] + "$" + s[1:]

    except:
        unit_string = 'M$_\odot$'
    field_symbol = '$' + field_symbol + '_'
    p_t = ''
    for pos_it in depth_array:
        axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=2., c='k')
        axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=2., c='k')
        axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=1., c=part_color[pos_it])
        axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=1., c=part_color[pos_it])
        circle = mpatches.Circle([particle_position[0][pos_it], particle_position[1][pos_it]], accretion_rad, fill=False, lw=1, edgecolor='k')
        axis.add_patch(circle)
        if annotate_field is not None:
            if units is not None:
                annotate_field = annotate_field.in_units(units)
            if unit_string == 'M$_\odot$':
                P_msun = str(np.round(annotate_field[pos_it], 2))
            else:
                if annotate_field[pos_it] == 0.0:
                    P_msun = '0.0'
                else:
                    P_msun = '{:0.1e}'.format(annotate_field[pos_it])
            if p_t == '':
                p_t = field_symbol+str(pos_it+1)+'$='+P_msun+unit_string
            else:
                p_t = p_t+', ' +field_symbol+str(pos_it+1)+'$='+P_msun+unit_string
            print("plotted particle at position =", particle_position[0][pos_it], particle_position[1][pos_it])
    if annotate_field is not None:
        part_text = axis.text((xmin + 0.01*(box_size)), (ymin + 0.03*(ymax-ymin)), p_t, va="center", ha="left", color='w', fontsize=fontgize_global)
        part_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        #axis.annotate(p_t, xy=(xmin + 0.01*(box_size), ymin + 0.03*(ymax-ymin)), va="center", ha="left", color='w', fontsize=fontgize_global)
        #print("Annotated particle field")
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
