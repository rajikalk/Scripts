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
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        yt_file = True
    except YTOutputNotIdentified:
        yt_file = False
    return yt_file

def find_sink_formation_time(files):
    sink_form = None
    for source in range(len(files)):
        f = h5py.File(files[source], 'r')
        if 'particlemasses' in f.keys():
            sink_form = f['time'][0]/yt.units.yr.in_units('s').value
            break
    if sink_form == None:
        sink_form = -1*f['time'][0]/yt.units.yr.in_units('s').value
    return sink_form

def get_image_mesh(file, zoom_times):
    f = h5py.File(file, 'r')
    for key in f.keys():
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

def generate_frame_times(files, dt, start_time=None, presink_frames=25):
    try:
        file = files[-1]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        sink_form_time = np.min(dd['particle_creation_time'].value/yt.units.yr.in_units('s').value)
        max_time = ds.current_time.in_units('yr').value - sink_form_time
    except YTOutputNotIdentified:
        f = h5py.File(files[-1], 'r')
        sink_form_time = find_sink_formation_time(files)
        max_time = f['time'][0]/yt.units.yr.in_units('s').value - sink_form_time
        f.close()

    if start_time == None:
        m_times = np.logspace(0.0, np.log10(sink_form_time), presink_frames) - sink_form_time
        m_times = m_times.tolist()
    elif start_time < 0.0:
        m_times = np.logspace(np.log10(start_time), np.log10(sink_form_time), presink_frames) - sink_form_time
        m_times = m_times.tolist()
    else:
        m_times = []

    postsink = 0.0
    while postsink < (max_time):
        if start_time != None:
            if postsink >= start_time:
                m_times.append(postsink)
        else:
            m_times.append(postsink)
        postsink = postsink + dt
    return m_times

def find_files(m_times, files):
    try:
        file = files[-1]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        sink_form_time = np.min(dd['particle_creation_time'].value/yt.units.yr.in_units('s').value)
        yt_file = True
    except YTOutputNotIdentified:
        sink_form_time = find_sink_formation_time(files)
        yt_file = False
    usable_files = []
    mit = 0
    min = 0
    max = len(files)-1
    pit = 0
    while mit < len(m_times):
        it = int(np.round(min + ((max - min)/2.)))
        #print 'search iterator =', it
        if yt_file:
            file = files[it]
            part_file=file[:-12] + 'part' + file[-5:]
            ds = yt.load(file, particle_filename=part_file)
            time = ds.current_time.in_units('yr').value-sink_form_time
        else:
            f = h5py.File(files[it], 'r')
            time = f['time'][0]/yt.units.yr.in_units('s').value-sink_form_time
        if pit == it or time == m_times[mit]:
            if time < 0:
                usable_files.append(files[it+1])
            else:
                usable_files.append(files[it])
            print "found time", time, "for m_time", m_times[mit], "with file:", usable_files[-1]
            mit = mit + 1
            min = it
            max = len(files)-1
            pit = it
        elif time > m_times[mit]:
            max = it
            pit = it
        elif time < m_times[mit]:
            min = it
            pit = it
    return usable_files

def get_particle_data(file, axis='xz'):
    part_mass = []
    part_pos_x = []
    part_pos_y = []
    accretion_rad = []
    try:
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        if axis == 'xy':
            part_pos_x = dd['particle_posx'].in_units('AU').value
            part_pos_y = dd['particle_posy'].in_units('AU').value
        else:
            part_pos_x = dd['particle_posx'].in_units('AU').value
            part_pos_y = dd['particle_posz'].in_units('AU').value
        part_mass = dd['particle_mass'].in_units('msun').value
        accretion_rad = np.min(dd['dx'].in_units('au').value) * 2.5
    except YTOutputNotIdentified:
        f = h5py.File(file, 'r')
        if axis == 'xy':
            part_pos_x = f["particlepositions"][0]/yt.units.au.in_units('cm').value
            part_pos_y = f["particlepositions"][1]/yt.units.au.in_units('cm').value
        else:
            part_pos_x = f["particlepositions"][0]/yt.units.au.in_units('cm').value
            part_pos_y = f["particlepositions"][2]/yt.units.au.in_units('cm').value
        part_mass = np.array(f["particlemasses"])/yt.units.msun.in_units('g').value
        accretion_rad = f['r_accretion'][0]/yt.units.au.in_units('cm').value
    part_info = {'particle_mass':part_mass,
                 'particle_position':[part_pos_x, part_pos_y],
                 'accretion_rad':accretion_rad}
    return part_info

def initialise_grid(file, zoom_times=0):#, center=0):
    f = h5py.File(file, 'r')
    for key in f.keys():
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

def get_quiver_arrays(velx_full, vely_full, no_of_quivers=32., smooth_cells=None):
    annotate_freq = float(np.shape(velx_full)[0])/float(no_of_quivers-1)
    if smooth_cells == None:
        smoothing_val = int(annotate_freq/2)
    else:
        smoothing_val = int(smooth_cells)
    x_ind = []
    y_ind = []
    counter = 0
    while counter < (no_of_quivers-1):
        val = int(annotate_freq*counter + annotate_freq/2.)
        x_ind.append(int(val))
        y_ind.append(int(val))
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
    return velx, vely

def my_own_quiver_function(axis, X_pos, Y_pos, X_val, Y_val, plot_velocity_legend='False', standard_vel=5, legend_text="5kms$^{-1}$", limits=None):
    global fontgize_global
    if plot_velocity_legend == 'False':
        plot_velocity_legend = False
    elif plot_velocity_legend == 'True':
        plot_velocity_legend = True
    standard_vel = yt.units.km.in_units('cm').value * standard_vel
    if limits == None:
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
        print("plotting quiver legend")
        pos_start = [xmax - 0.15*(xmax-xmin), ymin + 0.07*(ymax-ymin)]
        xvel = standard_vel/len_scale
        yvel = 0.0
        width_val = 1.0
        axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=1.*width_val, mutation_scale=15.*width_val))
        axis.annotate(legend_text, xy=(xmax - 0.01*(xmax-xmin), ymin + 0.03*(ymax-ymin)), va="center", ha="right", color='w', fontsize=fontgize_global)
    return axis

def annotate_particles(axis, particle_position, accretion_rad, limits, annotate_field=None, field_symbol='M', units=None):
    global fontgize_global
    if annotate_field != None and units != None:
        annotate_field = annotate_field.in_units(units)
    part_color = ['cyan','cyan','r','c','y','w','k']
    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    box_size = xmax - xmin
    if accretion_rad/box_size < 0.025:
        line_rad = 0.025*box_size
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
    for pos_it in range(len(particle_position[0])):
        axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=2., c='k')
        axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=2., c='k')
        axis.plot((particle_position[0][pos_it]-(line_rad), particle_position[0][pos_it]+(line_rad)), (particle_position[1][pos_it], particle_position[1][pos_it]), lw=1., c=part_color[pos_it])
        axis.plot((particle_position[0][pos_it], particle_position[0][pos_it]), (particle_position[1][pos_it]-(line_rad), particle_position[1][pos_it]+(line_rad)), lw=1., c=part_color[pos_it])
        circle = mpatches.Circle([particle_position[0][pos_it], particle_position[1][pos_it]], accretion_rad, fill=False, lw=1, edgecolor='k')
        axis.add_patch(circle)
        if annotate_field != None:
            if units != None:
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
            print("p_t =", p_t)
    if annotate_field != None:
        axis.annotate(p_t, xy=(xmin + 0.01*(box_size), ymin + 0.03*(ymax-ymin)), va="center", ha="left", color='w', fontsize=fontgize_global)
        print("Annotated particle field")
    return axis

def profile_plot(data, x_field, y_fields, weight_field=None, n_bins=100, log=False, x_units='AU', y_units=None, abs=False, center=0):
    myf.set_center(center)
    if x_units == None:
        x = data[x_field]
    else:
        x = data[x_field].in_units(x_units)
    if abs:
        x = np.abs(x)
    print("Got x values")
    if log == 'True':
        bins = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), n_bins+1)
    elif n_bins == None:
        unit_string = str(x.unit_quantity).split(' ')[-1]
        bins = list(set(x.value))
        bins = np.sort(bins)
        bins = yt.YTArray(bins, unit_string)
        bins = bins[::2]
    else:
        bins = np.linspace(np.min(x), np.max(x), n_bins+1)
    print("No of bins:", len(bins))


    if weight_field != None:
        weight = data[weight_field]
    x_array = []
    y_array = {}

    for y_field in y_fields:
        y_temp = []
        if y_units == None:
            y = data[y_field]
        else:
            y = data[y_field].in_units(y_units)
        y_units = y.units
        print("Got y values")
        prev_bin = bins[0]
        for bin in bins[1:]:
            ind = np.where((x >= prev_bin) & (x < bin))[0]
            mid_x = (bin+prev_bin)/2.
            x_array.append(mid_x)
            if len(ind) != 0:
                y_vals = y[ind]
                if weight_field != None:
                    bin_val = (np.sum(y[ind]*weight[ind]))/np.sum(weight[ind])
                else:
                    bin_val = np.mean(y[ind])
            else:
                bin_val = np.nan
            prev_bin = bin
            y_temp.append(bin_val)
            print "Value =", bin_val, ", at radius=", mid_x

    y_array.update({y_field:y_temp})

    return x_array, y_array

def sample_points(data, x_field, y_field, bin_no=2., no_of_points=2000, x_units='AU', center=0):
    myf.set_center(center)
    z = data['z'].in_units('AU')
    if x_units == None:
        x = data[x_field]
    else:
        x = data[x_field].in_units(x_units)
    big_array = np.array([z, x, data[y_field]])
    plot_array = []
    z_bins = np.linspace(np.min(np.abs(z)), np.max(np.abs(z)), bin_no+1)
    import random
    prev_z = z_bins[0]
    counter_max = no_of_points/bin_no
    for z_bin in z_bins[1:]:
        filtered_z = (z > prev_z)*(z < z_bin)
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

def sliceplot(ds, X, Y, field, cmap=plt.cm.get_cmap('brg'), log=False, resolution=1024, center=0, units=None, cbar_label="Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)"): #cmap=plt.cm.get_cmap('bwr_r')
    """
       Creates a slice plot along the xy-plane for a data cube. This is done by interpolating  the simulation output onto a grid.
    """
    print "SLICE PLOT CENTER =", center
    print("image resolution =", resolution)
    x = np.linspace(np.min(X), np.max(X), resolution)
    y = np.linspace(np.min(Y), np.max(Y), resolution)
    print("len(x) =", len(x))
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
    print("Finding nearest points...")
    nearest_points = tree.query(xyz_grid, k=2, distance_upper_bound=cell_max, n_jobs=-1)[1]
    
    print("len(nearest_points) =", len(nearest_points))
    del xyz_grid
    print("Done finding nearest points...")

    if center == 3:
        field_grid = np.zeros_like(xy[0].shape)
        for cen in range(len(dd['particle_mass'])+1):
            myf.set_center(cen)
            print("Doing center =", cen)
            del dd
            dd = ds.region([0.0,0.0,0.0], [np.min(X), np.min(Y), -cell_max], [np.max(X), np.max(Y), cell_max])
            field_grid = field_grid + (dd[field][nearest_points[:,0]].reshape(xy[0].shape) + dd[field][nearest_points[:,1]].reshape(xy[0].shape))/2.
        field_grid = field_grid/(len(dd['particle_mass'])+1)
    else:
        '''
        if no_of_particles == 2:
            if center == 1:
                myf.set_center(2)
            else:
                myf.set_center(1)
        else:
            myf.set_center(center)
        '''
        myf.set_center(center)
        del dd
        dd = ds.region([0.0,0.0,0.0], [np.min(X), np.min(Y), -cell_max], [np.max(X), np.max(Y), cell_max])
        field_grid = (dd[field][nearest_points[:,0]].reshape(xy[0].shape) + dd[field][nearest_points[:,1]].reshape(xy[0].shape))/2.

    if units != None:
        field_grid = field_grid.in_units(units)
    field_grid = field_grid
    
    plt.clf()
    fig, ax = plt.subplots()
    print("len(xy[0]) =", len(xy[0]))
    if log:
        plot = ax.pcolormesh(xy[0], xy[1], field_grid.value, cmap=cmap, norm=LogNorm(), rasterized=True)
    else:
        plot = ax.pcolormesh(xy[0], xy[1], field_grid.value, cmap=cmap, rasterized=True)
    plt.gca().set_aspect('equal')
    #cbar = plt.colorbar(plot, pad=0.0)
    #cbar.set_label(cbar_label, rotation=270, labelpad=13, size=14)
    '''
    if field == "Relative_Keplerian_Velocity":
       cbar.set_clim([0.0, 2.0])
    '''
    ax.set_xlim([np.min(X), np.max(X)])
    ax.set_ylim([np.min(Y), np.max(Y)])
    ax.set_xlabel('$x$ (AU)', labelpad=-1)
    ax.set_ylabel('$y$ (AU)', labelpad=-20)
    return fig, ax, xy, field_grid
