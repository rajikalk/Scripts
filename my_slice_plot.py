from __future__ import division, print_function
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import numpy as np
import scipy.spatial as spatial
import yt
import argparse
import os
import glob
plt.ion()

temp_yarr = yt.YTArray([0,1,2], 'g')
field_dictionary = {'temp_yarr':temp_yarr}

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="What is the input files. It can be an individual file or a series. If indivivdual, this will be plotted.")
    parser.add_argument("-sd", "--save_dir", help="Where will the output data be saved?")
    parser.add_argument("-zt", "--zoom_times", help="How much do you want to zoom in?", default=18.5, type=float)
    parser.add_argument("-pt", "--plot_time", help="would you like to plot a specific time?", default=None)
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default=2, type=float)
    parser.add_argument("-res", "--resolution", help="how many pixels by how many pixels do you want?", default=512, type=int)
    parser.add_argument("-qf", "--quiver_frequency", help="How many quivers along one edge do you want?", default=31, type=int)
    parser.add_argument("-c", "--center", help="would you like to center the profile on the center of mass, or tone of the particles?, 0=CoM, 1=particle 1, 2=particle 2", default=0, type=int)
    parser.add_argument("-fi", "--field", help="What field would you like to plot?", default="Relative_Keplerian_Velocity_Cyl", nargs='+')
    parser.add_argument("-ip", "--image_prefix", help="do you want your saved images to all have the same prefix?", default=None, type=str)
    parser.add_argument("-log", "--log_colour", help="When you plot the slice do you want the color on a log scale?", default=False, type=bool)
    parser.add_argument("-clim", "--colour_limits", help="colourbar limits", nargs="+", type=float, default=None)
    parser.add_argument("-tmax", "--time_max", help="what is the max time you would like to go up to?", default=None, type=int)
    parser.add_argument("-cl", "--colour_bar_label", help="what do you want as the colour bar label?", default=None, type=str)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", default=False, type=bool)
    parser.add_argument("-sc", "--smooth_cells", help="how many cells would you like the smooth the velocities over? If not defined it is set to half the annotatation frequency", default=None, type=int)
    args = parser.parse_args()
    return args

args = parse_inputs()

def generate_frame_times(files,sink_form_time,args):
    if args.time_max == None:
        file = files[-1]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        max_time = ds.current_time.in_units('yr').value - sink_form_time
    else:
        max_time = args.time_max
    m_times = []
    postsink = 0.0
    while postsink < max_time:
        m_times.append(postsink)
        postsink = postsink + args.time_step
    no_frames = len(m_times)
    return m_times

def find_files_bin(m_times, files, sink_form_time):
    usable_files = []
    mit = 0
    min = 0
    max = len(files)-1
    pit = 0
    while mit < len(m_times):
        it = int(np.round(min + ((max - min)/2.)))
        file = files[it]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        time = ds.current_time.in_units('yr').value-sink_form_time
        if pit == it or time == m_times[mit]:
            if m_times[mit] == 0.0:
                if time < 0.0:
                    it = it + 1
                    while time < 0.0:
                        file = files[it]
                        part_file=file[:-12] + 'part' + file[-5:]
                        ds = yt.load(file, particle_filename=part_file)
                        time = ds.current_time.in_units('yr').value-sink_form_time
                        if time > 0.0:
                            usable_files.append(file)
                        else:
                            it = it + 1
                else:
                    usable_files.append(file)
            else:
                usable_files.append(file)
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

#First find the time that you want to plot! This uses a simple binary search
def find_time(plot_time, files, sink_form_time):
    min = 0
    max = len(files)-1
    pit = 0
    it = 100
    found = False
    while found == False:
        it = int(np.round(min + ((max - min)/2.)))
        ds = yt.load(files[it])
        time_val = (ds.current_time.in_units('yr').value)-sink_form_time
        if pit == it:
            if plot_time == 0.0 and time_val < 0.0:
                it = it + 1
                ds = yt.load(files[it])
                time_val = (ds.current_time.in_units('yr').value)-sink_form_time
            found = True
            break
        elif time_val > plot_time:
            max = it
            pit = it
        elif time_val < plot_time:
            min = it
            pit = it
    return it

def _CoM(field, data):
    TM = np.sum(data['cell_mass'].in_units('g'))
    x_top = yt.YTArray(0.0, 'cm*g')
    y_top = yt.YTArray(0.0, 'cm*g')
    z_top = yt.YTArray(0.0, 'cm*g')
    if ('all', u'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(data['particle_mass'].in_units('g'))
        for part in range(len(data['particle_mass'])):
            x_top = x_top + data['particle_mass'][part].in_units('g')*data['particle_posx'][part].in_units('cm')
            y_top = y_top + data['particle_mass'][part].in_units('g')*data['particle_posy'][part].in_units('cm')
            z_top = z_top + data['particle_mass'][part].in_units('g')*data['particle_posz'][part].in_units('cm')
    x_top = x_top + np.sum(data['cell_mass'].in_units('g')*data['x'].in_units('cm'))
    y_top = y_top + np.sum(data['cell_mass'].in_units('g')*data['y'].in_units('cm'))
    z_top = z_top + np.sum(data['cell_mass'].in_units('g')*data['z'].in_units('cm'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    com = yt.YTArray(com, 'cm')
    return com

yt.add_field("CoM", function=_CoM, units=r"cm")

def _My_Bulk_Velocity(field, data):
    x_vel = np.sum(data['velx'].in_units('cm/s')*data['cell_mass'].in_units('g'))
    y_vel = np.sum(data['vely'].in_units('cm/s')*data['cell_mass'].in_units('g'))
    z_vel = np.sum(data['velz'].in_units('cm/s')*data['cell_mass'].in_units('g'))
    TM = np.sum(data['cell_mass'].in_units('g'))
    if ('all', u'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(data['particle_mass'].in_units('g'))
        for part in range(len(data['particle_mass'])):
            x_vel = x_vel + data['particle_mass'][part].in_units('g')*data['particle_velx'][part].in_units('cm/s')
            y_vel = y_vel + data['particle_mass'][part].in_units('g')*data['particle_vely'][part].in_units('cm/s')
            z_vel = z_vel + data['particle_mass'][part].in_units('g')*data['particle_velz'][part].in_units('cm/s')
    bv = [(x_vel/TM), (y_vel/TM), (z_vel/TM)]
    bv = yt.YTArray(bv, 'cm/s')
    return bv

yt.add_field("My_Bulk_Velocity", function=_My_Bulk_Velocity, units=r"cm/s")

def get_quiver_arrays(velx_full, vely_full):
    global args
    annotate_freq = float(args.resolution)/float(args.quiver_frequency)
    if args.smooth_cells == None:
        smoothing_val = int(annotate_freq/2)
    else:
        smoothing_val = int(args.smooth_cells)
    x_pos = []
    y_pos = []
    counter = 0
    while counter < 31:
        val = int(annotate_freq*counter + annotate_freq/2.)
        x_pos.append(int(val))
        y_pos.append(int(val))
        counter = counter + 1
    velx = []
    vely = []
    for x in x_pos:
        x = int(x)
        xarr = []
        yarr = []
        for y in y_pos:
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
    velx = yt.YTArray(velx, "cm/s")
    vely = yt.YTArray(vely, "cm/s")
    return velx, vely

def my_own_quiver_function(axis, X_pos, Y_pos, X_val, Y_val):
    global field_dictionary
    global args
    standard_vel = 500000.
    scale_factor = 100000.
    vel_scale_factor = 90000.
    xmin = np.min(field_dictionary['X'][0].in_units('AU').value)
    xmax = np.max(field_dictionary['X'][0].in_units('AU').value)
    ymin = np.min(field_dictionary['Y'][0].in_units('AU').value)
    len_scale = standard_vel/(0.04*(xmax - xmin))
    vels = np.hypot(X_val, Y_val)
    for xp in range(len(X_pos[0])):
        for yp in range(len(Y_pos[0])):
            xvel = X_val[xp][yp]/len_scale
            yvel = Y_val[xp][yp]/len_scale
            width_val = np.sqrt(X_val[xp][yp]**2. + Y_val[xp][yp]**2.)/standard_vel
            if width_val > 1.0:
                width_val = 1.0
            axis.add_patch(mpatches.FancyArrowPatch((X_pos[xp][yp], Y_pos[xp][yp]), (X_pos[xp][yp]+xvel, Y_pos[xp][yp]+yvel), color='w', linewidth=1.*width_val, arrowstyle='->', mutation_scale=15.*width_val, shrinkA=0.0, shrinkB=0.0))
    if args.plot_velocity_legend == True:
        print("plotting quiver legend")
        pos_start = [0.77*xmax, 0.87*ymin]
        xvel = standard_vel/len_scale
        yvel = 0.0
        width_val = 1.0
        axis.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=1.*width_val, mutation_scale=15.*width_val))
        axis.annotate("5kms$^{-1}$", xy=(0.98*xmax, 0.95*ymin), va="center", ha="right", color='w', fontsize=16)


def My_plotting_function(dd, field, colormap=None, save_name=None, image_prefix=None, log_colour=False, clim=None, units=None):
    global field_dictionary
    global args
    bx = field_dictionary['box_bound'].value
    if args.colour_bar_label != None:
        cbar_label = args.colour_bar_label
    else:
        cbar_label = field
        if units == None:
            field_data = field_dictionary[field]
        else:
            field_data = field_dictionary[field].in_units(units)
        if str(field_dictionary[field].unit_quantity).split(' ')[-1] != '':
            cbar_label = cbar_label + ' (' + str(field_dictionary[field].unit_quantity).split(' ')[-1] + ')'
    
    plt.clf()
    fig, ax = plt.subplots()
    if log_colour:
        plot = plt.pcolormesh(field_dictionary['X'].value, field_dictionary['Y'].value, field_dictionary[field].value, cmap='brg', norm=LogNorm(), rasterized=True)
    else:
        plot = plt.pcolormesh(field_dictionary['X'].value, field_dictionary['Y'].value, field_dictionary[field].value, cmap='brg', rasterized=True)
    print("plotted image")
    if clim != None:
        plt.clim(clim)
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label(cbar_label, rotation=270, labelpad=12, size=14)
    #print("plotted color bar")
    print("plotting quiver plot")
    my_own_quiver_function(ax, field_dictionary['X_quiver'].in_units('AU').value, field_dictionary['Y_quiver'].in_units('AU').value, field_dictionary['velx_quiver'].value, field_dictionary['vely_quiver'].value)
    if ('all', u'particle_mass') in dd.ds.field_list:
        part_color = ['lime','cyan','r','c','y','w','k']
        line_rad = np.min(dd['dx'].in_units('AU'))*2.5
        p_t = ''
        for pos_it in range(len(dd['particle_posx'])):
            ax.plot((dd['particle_posx'][pos_it].in_units('AU')-(line_rad), dd['particle_posx'][pos_it].in_units('AU')+(line_rad)), (dd['particle_posy'][pos_it].in_units('AU'), dd['particle_posy'][pos_it].in_units('AU')), lw=2., c='k')
            ax.plot((dd['particle_posx'][pos_it].in_units('AU'), dd['particle_posx'][pos_it].in_units('AU')), (dd['particle_posy'][pos_it].in_units('AU')-(line_rad), dd['particle_posy'][pos_it].in_units('AU')+(line_rad)), lw=2., c='k')
            ax.plot((dd['particle_posx'][pos_it].in_units('AU')-(line_rad), dd['particle_posx'][pos_it].in_units('AU')+(line_rad)), (dd['particle_posy'][pos_it].in_units('AU'), dd['particle_posy'][pos_it].in_units('AU')), lw=1., c=part_color[pos_it])
            ax.plot((dd['particle_posx'][pos_it].in_units('AU'), dd['particle_posx'][pos_it].in_units('AU')), (dd['particle_posy'][pos_it].in_units('AU')-(line_rad), dd['particle_posy'][pos_it].in_units('AU')+(line_rad)), lw=1., c=part_color[pos_it])
            circle = mpatches.Circle([dd['particle_posx'][pos_it].in_units('AU').value, dd['particle_posy'][pos_it].in_units('AU').value], line_rad.value, fill=False, lw=1, edgecolor='k')
            ax.add_patch(circle)
    plt.xlim([-bx, bx])
    plt.ylim([-bx, bx])
    plt.xlabel('$x$ (AU)', labelpad=-1)
    plt.ylabel('$y$ (AU)', labelpad=-20)

    if save_name == None:
        save_name = field +'_sliceplot.eps'
    if image_prefix != None:
        save_name = image_prefix + save_name
    plt.savefig(save_name, bbox_inches='tight', pad_inches = 0.02)
    print('Saved file', save_name)
    return fig

#========================================================================================

path = args.files
if args.save_dir == None:
    split_cur = path.split('/O')
    save_dir = split_cur[0] + '/YT_O' + split_cur[1]
else:
    save_dir = args.save_dir

if save_dir[-1] != '/':
    save_dir = save_dir + '/'

print("save_dir =", save_dir)
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

#Get sink formation time
data_files = sorted(glob.glob(args.files+'WIND_hdf5_plt_cnt*'))
file = data_files[-1]
init_ds = yt.load(file, particle_filename=file[:-12] + 'part' + file[-5:])
init_dd = init_ds.all_data()
yr = yt.units.yr.in_units('s')
sink_formation_time_s = yt.YTArray(init_dd['particle_creation_time'][0].value, 's')
sink_formation_time = sink_formation_time_s.in_units('yr').value
print("sink formation time =", sink_formation_time)

#Find relavent files
if args.plot_time != None:
    plot_time_val = float(args.plot_time)
    file_it = find_time(plot_time_val, data_files, sink_formation_time)
    usable_files = [data_files[file_it]]
    m_times = [plot_time_val]
else:
    m_times = generate_frame_times(data_files, sink_formation_time, args)
    usable_files = find_files_bin(m_times, data_files, sink_formation_time)
print("FOUND USABLE FILES")

#Initialise meshs
bx = (init_ds.domain_right_edge[0].in_units('AU')/args.zoom_times)/2.
nx = args.resolution       #number of points in the grid
nx_quiver = args.quiver_frequency #number of points for a quiver plot in the grid
x = np.linspace(-bx, bx, nx)
xy = np.meshgrid(x,x)
delta_quiver = (bx*2)/nx_quiver
annotate_freq = args.resolution/args.quiver_frequency
x_quiver = np.linspace(-bx+(delta_quiver/2.), bx-(delta_quiver/2.), (nx_quiver))
xy_quiver = np.meshgrid(x_quiver,x_quiver)
print("Created plot and quiver meshs")

#Start going through the file
fit = 1
for file in usable_files:
    print("doing file", fit, "of", len(usable_files))
    time = m_times[fit-1]
    part_file=file[:-12] + 'part' + file[-5:]
    #load file
    ds = yt.load(file, particle_filename=part_file)
    dd = ds.box(left_edge=[-bx.in_units('cm'),-bx.in_units('cm'),-bx.in_units('cm')], right_edge=[bx.in_units('cm'),bx.in_units('cm'),bx.in_units('cm')])
    dz = dd['z'][np.argmax(dd['dens'])].in_units('AU')
    
    #Set the center of the system for calculations
    if args.center == 0:
        center = dd['CoM'].in_units('AU')
        center_vel = dd['My_Bulk_Velocity'].in_units('cm/s')
    else:
        center = [dd['particle_posx'][args.center-1].in_units('AU'), dd['particle_posy'][args.center-1].in_units('AU'), dd['particle_posz'][args.center-1].in_units('AU')]
        center_vel = [dd['particle_velx'][args.center-1].in_units('cm/s'), dd['particle_vely'][args.center-1].in_units('cm/s'), dd['particle_velz'][args.center-1].in_units('cm/s')]

    #get xyz positions relative to the chosen center
    xyz = yt.YTArray([dd['x'].in_units('AU'),dd['y'].in_units('AU'),dd['z'].in_units('AU')]).T

    #This line takes a while
    print("Computing tree...")
    tree = spatial.cKDTree(xyz)
    print("Done computing tree...")

    #create grid of the pixel positions and quiver positions
    xyz_grid = yt.YTArray([xy[0].flatten(), xy[1].flatten(), dz*np.ones_like(xy[0].flatten().value)]).T

    #This line also takes a while
    print("Finding nearest points...")
    nearest_points = tree.query(xyz_grid, distance_upper_bound=np.max(dd['dx'].in_units('AU')), n_jobs=-1)
    print("Done finding nearest points...")

    #correct velocities:
    dens_grid = dd['dens'][nearest_points[1]].reshape(xy[0].shape)
    vx_grid = dd['velx'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s') - center_vel[0]
    vy_grid = dd['vely'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s') - center_vel[1]
    vz_grid = dd['velz'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s') - center_vel[2]
    vx_grid_quiver_full = dd['velx'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s')
    vy_grid_quiver_full = dd['vely'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s')
    vx_grid_quiver, vy_grid_quiver = get_quiver_arrays(vx_grid_quiver_full, vy_grid_quiver_full)


    #Try computing mass versus radius
    radii_cyl = np.sqrt((xyz.T[0].in_units('AU')-center[0])**2 + (xyz.T[1].in_units('AU')-center[1])**2)
    radii_sph = np.sqrt((xyz.T[0].in_units('AU')-center[0])**2 + (xyz.T[1].in_units('AU')-center[1])**2 + (xyz.T[2].in_units('AU')-center[2])**2)
    r_grid = np.sqrt(xy[0].in_units('AU')**2 + xy[1].in_units('AU')**2)
    r_corr_grid = np.sqrt((xy[0].in_units('AU')-center[0].in_units('AU'))**2. + (xy[1].in_units('AU')-center[1].in_units('AU'))**2.)

    bin_size = np.sqrt((x[-1]-x[-2])**2.+(x[-1]-x[-2])**2.)
    rs = np.arange(np.min(r_corr_grid), np.max(r_corr_grid)+bin_size, bin_size)
    mass_versus_r = np.ones_like(rs)
    masses = dd['cell_mass']

    #Work out the semimajor axis of the system, if there is one
    if ('all', u'particle_mass') in ds.field_list:
        if len(dd['particle_mass']) == 2:
            pos1 = [dd['particle_posx'][0].in_units('cm'), dd['particle_posy'][0].in_units('cm'), dd['particle_posz'][0].in_units('cm')]
            pos2 = [dd['particle_posx'][1].in_units('cm'), dd['particle_posy'][1].in_units('cm'), dd['particle_posz'][1].in_units('cm')]
            a = np.sqrt((pos1[0] - pos2[0])**2. + (pos1[1] - pos2[1])**2. + (pos1[2] - pos2[2])**2.)
        else:
            a = yt.YTArray(0.0, 'cm')
    else:
        a = yt.YTArray(0.0, 'cm')

    #Calculate enclosed mass in cylindrical coordinates
    print("Calculating cylindrical enclosed mass")
    enclosed_mass_grid_cyl = yt.YTArray(np.zeros(np.shape(r_corr_grid)), 'g')
    prev_r = rs[0]
    prev_enc_mass_cyl = yt.YTArray(0.0, 'g')
    included_particles = False
    for r in rs[1:]:
        ind = np.where((radii_cyl >= prev_r) & (radii_cyl < r))
        grid_ind = np.where((r_corr_grid >= prev_r) & (r_corr_grid < r))
        if len(ind) != 0:
            enclosed_mass_cyl = prev_enc_mass_cyl + np.sum(masses[ind])
            if r == rs[1] and args.center != 0:
                enclosed_mass_cyl = enclosed_mass_cyl + dd['particle_mass'][args.center-1]
            if included_particles == False:
                if args.center != 0 and r > a.in_units('AU'):
                    if args.center == 1:
                        enclosed_mass_cyl = enclosed_mass_cyl + dd['particle_mass'][1]
                    else:
                        enclosed_mass_cyl = enclosed_mass_cyl + dd['particle_mass'][0]
                    included_particles = True
                elif args.center == 0 and r > a.in_units('AU')/2. and ('all', u'particle_mass') in ds.field_list:
                    enclosed_mass_cyl = enclosed_mass_cyl + np.sum(dd['particle_mass'])
                    included_particles = True
            enclosed_mass_grid_cyl[grid_ind] = enclosed_mass_cyl
            prev_enc_mass_cyl = enclosed_mass_cyl
            prev_r = r

    #Then again in spherical coordinates
    #Calculate enclosed mass in cylindrical coordinates
    print("Calculating spherical enclosed mass")
    enclosed_mass_grid_sph = yt.YTArray(np.zeros(np.shape(r_corr_grid)), 'g')
    prev_r = rs[0]
    prev_enc_mass_sph = yt.YTArray(0.0, 'g')
    included_particles = False
    for r in rs[1:]:
        ind = np.where((radii_sph >= prev_r) & (radii_sph < r))
        grid_ind = np.where((r_corr_grid >= prev_r) & (r_corr_grid < r))
        if len(ind) != 0:
            enclosed_mass_sph = prev_enc_mass_sph + np.sum(masses[ind])
            if r == rs[1] and args.center != 0:
                enclosed_mass_sph = enclosed_mass_sph + dd['particle_mass'][args.center-1]
            if included_particles == False:
                if args.center != 0 and r > a.in_units('AU'):
                    if args.center == 1:
                        enclosed_mass_sph = enclosed_mass_sph + dd['particle_mass'][1]
                    else:
                        enclosed_mass_sph = enclosed_mass_sph + dd['particle_mass'][0]
                    included_particles = True
                elif args.center == 0 and r > a.in_units('AU')/2. and ('all', u'particle_mass') in ds.field_list:
                    enclosed_mass_sph = enclosed_mass_sph + np.sum(dd['particle_mass'])
                    included_particles = True
            enclosed_mass_grid_sph[grid_ind] = enclosed_mass_sph
            prev_enc_mass_sph = enclosed_mass_sph
            prev_r = r

    #Work out the velocities
    vkep_grid_cyl = np.sqrt(yt.utilities.physical_constants.G*enclosed_mass_grid_cyl.in_units('g')/r_corr_grid.in_units('cm'))
    vkep_grid_sph = np.sqrt(yt.utilities.physical_constants.G*enclosed_mass_grid_sph.in_units('g')/r_corr_grid.in_units('cm'))

    vrad_grid_cyl = (xy[0].in_units('cm')*vx_grid.in_units('cm/s') + xy[1]*vy_grid.in_units('cm/s'))/r_grid.in_units('cm')
    vmag_grid_cyl = np.sqrt(vx_grid.in_units('cm/s')**2 + vy_grid.in_units('cm/s')**2)

    vrad_grid_sph = (xy[0]*vx_grid.in_units('cm/s') + xy[1]*vy_grid.in_units('cm/s') + dz*vz_grid.in_units('cm/s'))/r_grid.in_units('cm')
    vmag_grid_sph = np.sqrt(vx_grid.in_units('cm/s')**2 + vy_grid.in_units('cm/s')**2 + vz_grid.in_units('cm/s')**2)

    vtan_grid_cyl = np.sqrt(vmag_grid_cyl.in_units('cm/s')**2 - vrad_grid_cyl.in_units('cm/s')**2)
    vtan_grid_sph = np.sqrt(vmag_grid_sph.in_units('cm/s')**2 - vrad_grid_sph.in_units('cm/s')**2)

    vrel_grid_cyl = vtan_grid_cyl.in_units('cm/s')/vkep_grid_cyl.in_units('cm/s')
    vrel_grid_sph = vtan_grid_sph.in_units('cm/s')/vkep_grid_sph.in_units('cm/s')

    #Save fields to dictionaries. Should make plotting easier
    field_dictionary.pop('temp_yarr', None)
    field_update = {'Enclosed_Mass_Cyl':enclosed_mass_grid_cyl,
                    'Keplerian_Velocity_Cyl':vkep_grid_cyl,
                    'Radial_Velocity_Cyl':vrad_grid_cyl,
                    'Velocity_Magnitude_Cyl':vmag_grid_cyl,
                    'Tangential_Velocity_Cyl':vtan_grid_cyl,
                    'Relative_Keplerian_Velocity_Cyl':vrel_grid_cyl,
                    'Enclosed_Mass_Sph':enclosed_mass_grid_sph,
                    'Keplerian_Velocity_Sph':vkep_grid_sph,
                    'Radial_Velocity_Sph':vrad_grid_sph,
                    'Velocity_Magnitude_Sph':vmag_grid_sph,
                    'Tangential_Velocity_Sph':vtan_grid_sph,
                    'Relative_Keplerian_Velocity_Sph':vrel_grid_sph,
                    'quiver_grids':xy_quiver,
                    'velx_quiver':vx_grid_quiver,
                    'vely_quiver':vy_grid_quiver,
                    'X':xy[0],
                    'Y':xy[1],
                    'X_quiver':xy_quiver[0],
                    'Y_quiver':xy_quiver[1],
                    'box_bound':bx,
                    'Density':dens_grid
                    }
    field_dictionary.update(field_update)

    #Now plot
    for field_val in args.field:
        My_plotting_function(dd, field_val, image_prefix=args.image_prefix, save_name=field_val+'_time_'+str(int(time))+'.eps', log_colour=args.log_colour, clim=args.colour_limits)

    fit = fit + 1
