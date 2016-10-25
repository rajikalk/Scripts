from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as spatial
import yt
import argparse
import os
import glob
plt.ion()

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="What is the input files. It can be an individual file or a series. If indivivdual, this will be plotted.")
    parser.add_argument("-sd", "--save_dir", help="Where will the output data be saved?")
    parser.add_argument("-zt", "--zoom_times", help="How much do you want to zoom in?", default=18.0, type=float)
    parser.add_argument("-pt", "--plot_time", help="would you like to plot a specific time?", default=None)
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default=2, type=float)
    parser.add_argument("-res", "--resolution", help="how many pixels by how many pixels do you want?", default=512, type=int)
    parser.add_argument("-qf", "--quiver_frequency", help="How many pixles do you want between quivers?", default=32, type=int)
    parser.add_argument("-c", "--center", help="would you like to center the profile on the center of mass, or tone of the particles?, 0=CoM, 1=particle 1, 2=particle 2", default=0, type=int)
    parser.add_argument("-coords", "--coordinates", help="Do you want to use cyclindrical or spherical coordinates?", default="cyl", type=str)
    parser.add_argument("-fi", "--field", help="What field would you like to plot?", default="Relative_Keplerian_Velocity", type=str)
    parser.add_argument
    args = parser.parse_args()
    return args

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
        #print 'search iterator =', it
        file = files[it]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        time = ds.current_time.in_units('yr').value-sink_form_time
        #print "it =", it
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
                            print("found time", time, "for m_time", m_times[mit])
                        else:
                            it = it + 1
                else:
                    usable_files.append(file)
                    print("found time", time, "for m_time", m_times[mit])
            else:
                usable_files.append(file)
                print("found time", time, "for m_time", m_times[mit])
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
            print("found time", time_val)
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
    #Center of mass calculation:
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

#================================================================================================

args = parse_inputs()

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

data_files = sorted(glob.glob(args.files+'WIND_hdf5_plt_cnt*'))
file = data_files[-1]
init_ds = yt.load(file, particle_filename=file[:-12] + 'part' + file[-5:])
init_dd = init_ds.all_data()
yr = yt.units.yr.in_units('s')
sink_formation_time_s = yt.YTArray(init_dd['particle_creation_time'][0].value, 's')
sink_formation_time = sink_formation_time_s.in_units('yr').value
print("sink formation time =", sink_formation_time)

if args.plot_time != None:
    plot_time_val = float(args.plot_time)
    file_it = find_time(plot_time_val, data_files, sink_formation_time)
    usable_files = [data_files[file_it]]
    m_times = [plot_time_val]
else:
    m_times = generate_frame_times(data_files, sink_formation_time, args)
    usable_files = find_files_bin(m_times, data_files, sink_formation_time)
print("FOUND USABLE FILES")

if "cyl" in args.coordinates:
    coord_sys = "cyl"
else:
    coord_sys = "sph"

bx = (init_ds.domain_right_edge[0].in_units('AU')/args.zoom_times)/2.
nx = args.resolution       #number of points in the grid
nx_quiver = args.quiver_frequency #number of points for a quiver plot in the grid
x = np.linspace(-bx, bx, nx)
xy = np.meshgrid(x,x)
x_quiver = np.linspace(-bx, bx, nx_quiver)
xy_quiver = np.meshgrid(x_quiver,x_quiver)
print("Created plot and quiver meshs")

fit = 1
for file in usable_files:
    print("doing file", fit, "of", len(usable_files))
    time = m_times[fit-1]
    part_file=file[:-12] + 'part' + file[-5:]
    #load file
    ds = yt.load(file, particle_filename=part_file)
    dd = ds.box(left_edge=[-bx.in_units('cm'),-bx.in_units('cm'),-bx.in_units('cm')], right_edge=[bx.in_units('cm'),bx.in_units('cm'),bx.in_units('cm')])
    dz = dd['z'][np.argmax(dd['dens'])].in_units('AU')
    
    if args.center == 0:
        center = dd['CoM'].in_units('AU')
    else:
        center = [dd['particle_posx'][args.center-1].in_units('AU'), dd['particle_posy'][args.center-1].in_units('AU'), dd['particle_posz'][args.center-1].in_units('AU')]

    xyz = yt.YTArray([dd['x'].in_units('AU')-center[0],dd['y'].in_units('AU')-center[1],dd['z'].in_units('AU')-center[2]]).T

    #This line takes a while
    print("Computing tree...")
    tree = spatial.KDTree(xyz)
    print("Done computing tree...")

    xyz_grid = yt.YTArray([xy[0].flatten(), xy[1].flatten(), dz*np.ones_like(xy[0].flatten().value)]).T
    xyz_grid_quiver = yt.YTArray([xy_quiver[0].flatten(), xy_quiver[1].flatten(), np.zeros_like(xy_quiver[0].flatten())]).T

    #This line also takes a while
    print("Finding nearest points...")

    nearest_points_quiver = tree.query(xyz_grid_quiver)
    nearest_points = tree.query(xyz_grid)
    print("Done finding nearest points...")

    #correct velocities:
    vx_grid = dd['velx'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s') - dd['My_Bulk_Velocity'][0].in_units('cm/s')
    vy_grid = dd['vely'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s') - dd['My_Bulk_Velocity'][1].in_units('cm/s')
    vz_grid = dd['velz'][nearest_points[1]].reshape(xy[0].shape).in_units('cm/s') - dd['My_Bulk_Velocity'][2].in_units('cm/s')

    #Try computing mass versus radius
    if coord_sys == "cyl":
        radii = np.sqrt(xyz.T[0].in_units('AU')**2 + xyz.T[1].in_units('AU')**2)
    else:
        radii = np.sqrt(xyz.T[0].in_units('AU')**2 + xyz.T[1].in_units('AU')**2 + xyz.T[2].in_units('AU')**2)
    r_grid = np.sqrt(xy[0].in_units('AU')**2 + xy[1].in_units('AU')**2)

    rs = x[len(x)//2:]
    mass_versus_r = np.ones_like(rs)
    masses = dd['cell_mass']

    if ('all', u'particle_mass') in ds.field_list:
        if len(dd['particle_mass']) == 2:
            pos1 = [dd['particle_posx'][0].in_units('cm'), dd['particle_posy'][0].in_units('cm'), dd['particle_posz'][0].in_units('cm')]
            pos2 = [dd['particle_posx'][1].in_units('cm'), dd['particle_posy'][1].in_units('cm'), dd['particle_posz'][1].in_units('cm')]
            a = np.sqrt((pos1[0] - pos2[0])**2. + (pos1[1] - pos2[1])**2. + (pos1[2] - pos2[2])**2.)
        else:
            a = yt.YTArray(0.0, 'cm')
    else:
        a = yt.YTArray(0.0, 'cm')

    for ix, r in enumerate(rs):
        enclosed_mass = np.sum(masses[radii<r])
        print("enclosed_mass(1) =", enclosed_mass)
        if args.center != 0 and r > a.in_units('AU'):
            enclosed_mass = enclosed_mass + np.sum(dd['particle_mass'])
        elif args.center != 0:
            enclosed_mass = enclosed_mass + dd['particle_mass'][args.center-1]
        elif args.center == 0 and r > a.in_units('AU')/2. and ('all', u'particle_mass') in ds.field_list:
            enclosed_mass = enclosed_mass + np.sum(dd['particle_mass'])
        print("enclosed_mass(2) =", enclosed_mass)
        mass_versus_r[ix] = enclosed_mass

    enclosed_mass_grid = yt.YTArray(np.interp(r_grid, rs, mass_versus_r))

    vkep_grid = np.sqrt(yt.utilities.physical_constants.G*enclosed_mass_grid/r_grid.in_units('cm'))

    if coord_sys == "cyl":
        vrad_grid = (xy[0].in_units('cm')*vx_grid + xy[1]*vy_grid)/r_grid
        vmag_grid = np.sqrt(vx_grid**2 + vy_grid**2)
    else:
        vrad_grid = (xy[0]*vx_grid + xy[1]*vy_grid + dz*vz_grid)/r_grid
        vmag_grid = np.sqrt(vx_grid**2 + vy_grid**2 + vz_grid**2)

    vtan_grid = np.sqrt(vmag_grid**2 - vrad_grid**2)
    vrel_grid = vtan_grid/vkep_grid

    #Save fields to dictionaries. Should make plotting easier
    field_dictionary = {'Enclosed_Mass':enclosed_mass_grid,
                        'Keplerian_Velocity':vkep_grid,
                        'Radial_Velocity':vrad_grid,
                        'Velocity_Magnitude':vmag_grid,
                        'Tangential_Velocity':vtan_grid,
                        'Relative_Keplerian_Velocity':vrel_grid
                        }

    #Now plot
    plt.clf()
    plt.imshow(field_dictionary[args.field], extent=[-bx, bx, -bx, bx], cmap='brg')
    cbar = plt.colorbar(pad=0.0)
    cbar.set_label(str(field_dictionary[args.field].unit_quantity).split(' ')[-1], rotation=270, size=14, labelpad=11)
    plt.quiver(xy_quiver[0], xy_quiver[1], vx_grid_quiver, vy_grid_quiver)
    if ('all', u'particle_mass') in ds.field_list:
        part_color = ['lime','cyan','r','c','y','w','k']
        line_rad = np.min(dd['dx'].value)*2.5
        for pos_it in range(len(dd['particle_posx'])):
            plt.plot((dd['particle_posx'][pos_it].value-(line_rad), dd['particle_posx'][pos_it].value+(line_rad)), (dd['particle_posy'][pos_it].value, dd['particle_posy'][pos_it].value), lw=2., c='k')
            plt.plot((dd['particle_posx'][pos_it].value, dd['particle_posx'][pos_it].value), (dd['particle_posy'][pos_it].value-(line_rad), dd['particle_posy'][pos_it].value+(line_rad)), lw=2., c='k')
            plt.plot((dd['particle_posx'][pos_it].value-(line_rad), dd['particle_posx'][pos_it].value+(line_rad)), (dd['particle_posy'][pos_it].value, dd['particle_posy'][pos_it].value), lw=1., c=part_color[pos_it])
            plt.plot((dd['particle_posx'][pos_it].value, dd['particle_posx'][pos_it].value), (dd['particle_posy'][pos_it].value-(line_rad), dd['particle_posy'][pos_it].value+(line_rad)), lw=1., c=part_color[pos_it])
    plt.savefig('kep_velocity.png')
    fit = fit + 1
