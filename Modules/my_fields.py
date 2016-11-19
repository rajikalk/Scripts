#!/usr/bin/env python
import yt
import numpy as np

center = 0
bin_size = yt.YTArray(100., 'AU')
image_resolution = 512
coordinates = 'cylindrical'

def set_center(x):
    """
    Sets the center used when calculateing fields.
    
    Type: int
    Default: 0
    Options:0=center of mass, 1=particle 1, 2=particle 2.
    """
    global center
    center = x
    return center

def set_bin_size(x):
    """
    sets the bin size used when calculating the enclosed mass
    Type: int or float
    """
    global bin_size
    bin_size = yt.YTArray(x, 'AU')
    return bin_size

def set_image_resolution(x):
    """
    Sets image resolution (assuming square images) of slices of fields with dependancies on the enclosed mass.
    Type: int
    """
    global image_resolution
    image_resolution = x
    return image_resolution

def set_coordinate_system(x):
    """
    Sets the coordinates system used for some fields.
    
    Default: 'cylindrical'
    Options: 'cylindrical', 'spherical'
    """
    global coordinates
    coordinates = x
    return coordinates

def get_center():
    """
    returns the currently set center
    """
    global center
    return center

def get_bin_size():
    """
    returns the currently set bin_size
    """
    global bin_size
    return bin_size

def get_image_resolution():
    """
    returns the currently set image resolution
    """
    global image_resolution
    return image_resolution

def get_coordinate_system(x):
    """
    returns the currently set coordinate system
    """
    global coordinates
    return coordinates

def _CoM(field, data):
    """
    Calculate the center of mass. Always includes particles where possible.
    """
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

def _Semimajor_Axis(field, data):
    """
    Calculates the semimajor axis of the binary. If there are no stars, or only one star, the semimajor axis is set to zero.
    """
    if ('all', u'particle_mass') in data.ds.field_list:
        if len(data['particle_mass']) == 2:
            pos1 = [data['particle_posx'][0].in_units('cm'), data['particle_posy'][0].in_units('cm'), data['particle_posz'][0].in_units('cm')]
            pos2 = [data['particle_posx'][1].in_units('cm'), data['particle_posy'][1].in_units('cm'), data['particle_posz'][1].in_units('cm')]
            a = np.sqrt((pos1[0] - pos2[0])**2. + (pos1[1] - pos2[1])**2. + (pos1[2] - pos2[2])**2.)
        else:
            a = yt.YTArray(0.0, 'cm')
    else:
        a = yt.YTArray(0.0, 'cm')
    return a

yt.add_field("Semimajor_Axis", function=_Semimajor_Axis, units=r"cm")

def _My_Bulk_Velocity(field, data):
    """
    Calculates the bulk velocity. Always includes particles where possible.
    """
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

def _Center_Position(field, data):
    """
    Returns the center position for the current set center.
    """
    global center
    if center == 0:
        center_pos = data['CoM'].in_units('cm')
    else:
        center_pos = [data['particle_posx'][center-1].in_units('cm'), data['particle_posy'][center-1].in_units('cm'), data['particle_posz'][center-1].in_units('cm')]
    print("center_pos =", center_pos)
    return center_pos

yt.add_field("Center_Position", function=_Center_Position, units=r"cm")

def _Center_Velocity(field, data):
    """
    Returns the center velocity for the current set center.
    """
    global center
    if center == 0:
        center_vel = data['My_Bulk_Velocity'].in_units('cm/s')
    else:
        center_vel = [dd['particle_velx'][center-1].in_units('cm/s'), dd['particle_vely'][center-1].in_units('cm/s'), dd['particle_velz'][center-1].in_units('cm/s')]
    print("center_vel", center_vel)
    return center_vel

yt.add_field("Center_Velocity", function=_Center_Velocity, units=r"cm/s")

#Lets also create a new field with the distance of cells from the CoM
def _dx_from_Center(field, data):
    """
    Calculates the change in x position from the current set center.
    """
    dx = data['x'].in_units('cm')-data['Center_Position'][0]
    return dx

yt.add_field("dx_from_Center", function=_dx_from_Center, units=r"cm")

def _dy_from_Center(field, data):
    """
    Calculates the change in y position from the current set center.
    """
    dy = data['y'].in_units('cm')-data['Center_Position'][1]
    return dy

yt.add_field("dy_from_Center", function=_dy_from_Center, units=r"cm")

def _dz_from_Center(field, data):
    """
    Calculates the change in z position from the current set center.
    """
    dz = data['z'].in_units('cm')-data['Center_Position'][2]
    return dz

yt.add_field("dz_from_Center", function=_dz_from_Center, units=r"cm")

def _Distance_from_Center(field, data):
    """
    Calculates the distance from the current set center.
    """
    distance = np.sqrt((data['dx_from_Center'])**2. + (data['dy_from_Center'])**2. + (data['dz_from_Center'])**2.)
    return distance

yt.add_field("Distance_from_Center", function=_Distance_from_Center, units=r"cm")

def _Cylindrical_Distance_from_Center(field, data):
    """
    Calculates the cylindrical distance from the current set center.
    """
    distance = np.sqrt((data['dx_from_Center'])**2. + (data['dy_from_Center'])**2.)
    return distance

yt.add_field("Cylindrical_Distance_from_Center", function=_Cylindrical_Distance_from_Center, units=r"cm")

def _Corrected_velx(field, data):
    """
    Calculates the x-velocity correcnted for the bulk velocity.
    """
    dvx = data['velx'].in_units('cm/s') - data['My_Bulk_Velocity'][0].in_units('cm/s')
    return dvx

yt.add_field("Corrected_velx", function=_Corrected_velx, units=r"cm/s")

def _Corrected_vely(field, data):
    """
    Calculates the y-velocity correcnted for the bulk velocity.
    """
    dvy = data['vely'].in_units('cm/s') - data['My_Bulk_Velocity'][1].in_units('cm/s')
    return dvy

yt.add_field("Corrected_vely", function=_Corrected_vely, units=r"cm/s")

def _Corrected_velz(field, data):
    """
    Calculates the z-velocity correcnted for the bulk velocity.
    """
    dvz = data['velz'].in_units('cm/s') - data['My_Bulk_Velocity'][2].in_units('cm/s')
    return dvz

yt.add_field("Corrected_velz", function=_Corrected_velz, units=r"cm/s")

def _Corrected_vel_mag(field, data):
    """
    Calculates the velocity magnitude corrected for the bulk velocity
    """
    distance = np.sqrt((data['Corrected_velx'])**2. + (data['Corrected_vely'])**2. + (data['Corrected_velz'])**2.)
    return distance

yt.add_field("Corrected_vel_mag", function=_Corrected_vel_mag, units=r"cm/s")

def Enclosed_Mass(data):
    global center
    global bin_size
    global coordinates
    import pdb
    pdb.set_trace()
    a = data['Semimajor_Axis']
    center_pos = data['Center_Position']
    center_vel = data['Center_Velocity']

    if 'cyl' in coordinates.lower():
        distance = data['Cylindrical_Distance_from_Center']
    else:
        distance = data['Distance_from_Center']
    rs = np.arange(np.min(distance), np.max(distance)+bin_size.in_units('cm'), bin_size.in_units('cm'))
    enclosed_mass = yt.YTArray(np.zeros(np.shape(distance)), 'g')

    prev_r = rs[0]
    prev_enc_mass = yt.YTArray(0.0, 'g')
    included_particles = False
    for r in rs[1:]:
        ind = np.where((distance >= prev_r) & (distance < r))
        if len(ind) != 0:
            enclosed_mass_val = prev_enc_mass + np.sum(data['cell_mass'][ind])
            if r == rs[1] and center != 0:
                enclosed_mass_val = enclosed_mass_val + data['particle_mass'][center-1]
            if included_particles == False:
                if center != 0 and r > a.in_units('AU'):
                    if center == 1:
                        enclosed_mass_val = enclosed_mass_val + data['particle_mass'][0]
                    else:
                        enclosed_mass_val = enclosed_mass_val + data['particle_mass'][1]
                    included_particles = True
                elif center == 0 and r > a.in_units('AU')/2. and ('all', u'particle_mass') in data.ds.field_list:
                    enclosed_mass_val = enclosed_mass_val + np.sum(data['particle_mass'])
                    included_particles = True
            enclosed_mass[ind] = enclosed_mass_val
            prev_enc_mass = enclosed_mass_val
            prev_r = r
    return enclosed_mass

def _Enclosed_Mass(field, data):
    dd = data.ds.all_data()
    indices_x = np.where(np.in1d(dd[('index', 'x')], data['x']))[0]
    indices_y = np.where(np.in1d(dd[('index', 'y')], data['y']))[0]
    indices_z = np.where(np.in1d(dd[('index', 'z')], data['z']))[0]
    inds = list(set(indices_x).intersection(indices_y).intersection(indices_z))
    enclosed_mass = Enclosed_Mass(dd)
    #enclosed_mass = enclosed_mass[inds]
    return enclosed_mass

yt.add_field("Enclosed_Mass", function=_Enclosed_Mass, units=r"")





'''
    #because orbital mechanics get a bit complicated in binaries, we try to simplify the problem.
    #I try to define regions, which tells us which center to use when calculating the keplerian velocity.
    #Region 1: closer to particle 1 than particle 2, but <a
    #Region 2: closer to particle 2 than particle 1, but <a
    #Region 3: r>a
    def _Defined_regions(field, data):
    distance = data['Distance_from_CoM'].in_units('cm')
    field_array = np.zeros(len(distance))
    if ('all', u'particle_mass') in data.ds.field_list:
    if len(data['particle_mass']) == 2:
    semimajor_axis = data['Semimajor_Axis']
    distance_from_a = np.sqrt((data['x'].in_units('cm') - data['particle_posx'][0].in_units('cm'))**2. + (data['y'].in_units('cm') - data['particle_posy'][0].in_units('cm'))**2. + (data['z'].in_units('cm') - data['particle_posz'][0].in_units('cm'))**2.)
    distance_from_b = np.sqrt((data['x'].in_units('cm') - data['particle_posx'][1].in_units('cm'))**2. + (data['y'].in_units('cm') - data['particle_posy'][1].in_units('cm'))**2. + (data['z'].in_units('cm') - data['particle_posz'][1].in_units('cm'))**2.)
    region_c = (distance > semimajor_axis)
    region_a = (distance < semimajor_axis)*np.less(distance_from_a, distance_from_b)
    region_b = (distance < semimajor_axis)*np.less(distance_from_b, distance_from_a)
    field_array = field_array+region_c*3.
    field_array = field_array+region_b*2.
    field_array = field_array+region_a*1.
    else:
    field_array = field_array+3.
    else:
    field_array = field_array + 3.
    regions = yt.YTArray(field_array, '')
    return regions
    
    yt.add_field("Defined_regions", function=_Defined_regions, units=r"")
    
    #Here I calculate the enclosed mass.
    def _Enclosed_Mass(field, data):
    distance = data['Distance_from_CoM'].in_units('cm')
    field_array = np.zeros(len(distance))
    AU = yt.units.astronomical_unit.in_units('cm')
    #For the single star case:
    if ('all', u'particle_mass') in data.ds.field_list:
    if len(data['particle_mass']) == 1:
    #I use radial bins to calculate the enclosed mass.
    #These are currently linearly spaced, not sure I should use like a logarithmic bins
    #bin size is currently 100 Au.. bit big, but for testing!
    r_bins = range(0*AU, np.max(distance), 10.*AU)
    r_bins.append(np.max(distance) + 10.*AU)
    r_bins = yt.YTArray(r_bins, 'cm')
    prev_r = r_bins[0]
    for r in r_bins[1:]:
    filtered_dist = (distance > prev_r)*(distance < r) #finds cells in radial shell
    enc_dist = distance < r #gets enclosed cells
    enclosed_mass = np.sum(enc_dist*data['cell_mass'].in_units('g'))
    enclosed_mass = enclosed_mass + data['particle_mass'][0]
    field_array = field_array+(filtered_dist*enclosed_mass.value)
    prev_r = r
    #but gets messy with binaries:
    elif len(data['particle_mass']) == 2:
    #calculates distance from stars a and b
    distance_from_a = np.sqrt((data['x'].in_units('cm') - data['particle_posx'][0].in_units('cm'))**2. + (data['y'].in_units('cm') - data['particle_posy'][0].in_units('cm'))**2. + (data['z'].in_units('cm') - data['particle_posz'][0].in_units('cm'))**2.)
    distance_from_b = np.sqrt((data['x'].in_units('cm') - data['particle_posx'][1].in_units('cm'))**2. + (data['y'].in_units('cm') - data['particle_posy'][1].in_units('cm'))**2. + (data['z'].in_units('cm') - data['particle_posz'][1].in_units('cm'))**2.)
    regions = data['Defined_regions']
    semimajor_axis = data['Semimajor_Axis']
    #first calculate for region around first star:
    #refer to drawing in log book for 09/08/2016 for max dist justification
    max_dist = np.sqrt((semimajor_axis**2.) + (semimajor_axis/2.)**2.)
    r_bins = range(0*AU, max_dist, 2.*AU)
    r_bins = yt.YTArray(r_bins, 'cm')
    prev_r = r_bins[0]
    for r in r_bins[1:]:
    filtered_dist = (distance_from_a > prev_r)*(distance_from_a < r)*(regions==1)
    enc_dist = distance_from_a < r
    enclosed_mass = np.sum(enc_dist*data['cell_mass'].in_units('g'))
    enclosed_mass = enclosed_mass + data['particle_mass'][0]
    field_array = field_array+(filtered_dist*enclosed_mass.value)
    prev_r = r
    #Do again for the region around the second star:
    prev_r = r_bins[0]
    for r in r_bins[1:]:
    filtered_dist = (distance_from_b > prev_r)*(distance_from_b < r)*(regions==2)
    enc_dist = distance_from_b < r
    enclosed_mass = np.sum(enc_dist*data['cell_mass'].in_units('g'))
    enclosed_mass = enclosed_mass + data['particle_mass'][1]
    field_array = field_array+(filtered_dist*enclosed_mass.value)
    prev_r = r
    #Then for the circumbinary area:
    r_bins = range(semimajor_axis, np.max(distance), 10.*AU)
    r_bins.append(np.max(distance) + 10.*AU)
    r_bins = yt.YTArray(r_bins, 'cm')
    prev_r = r_bins[0]
    for r in r_bins[1:]:
    filtered_dist = (distance > prev_r)*(distance < r)*(regions==3)
    enc_dist = distance < r
    enclosed_mass = np.sum(enc_dist*data['cell_mass'].in_units('g'))
    enclosed_mass = enclosed_mass + np.sum(data['particle_mass'])
    field_array = field_array+(filtered_dist*enclosed_mass.value)
    prev_r = r
    else:
    r_bins = range(0*AU, np.max(distance), 10.*AU)
    r_bins.append(np.max(distance) + 10.*AU)
    r_bins = yt.YTArray(r_bins, 'cm')
    prev_r = r_bins[0]
    for r in r_bins[1:]:
    filtered_dist = (distance > prev_r)*(distance < r) #finds cells in radial shell
    enc_dist = distance < r #gets enclosed cells
    enclosed_mass = np.sum(enc_dist*data['cell_mass'].in_units('g'))
    field_array = field_array+(filtered_dist*enclosed_mass.value)
    prev_r = r
    enclosed_mass_field = yt.YTArray(field_array, 'g')
    return enclosed_mass_field
    
    yt.add_field("Enclosed_Mass", function=_Enclosed_Mass, units=r"g")
    '''

