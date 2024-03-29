#!/usr/bin/env python
import yt
yt.unit_system_registry["cgs"]
import numpy as np
from subprocess import call
import pickle
import os
import time

center = 0
n_bins = 100
adaptive_bins = True
coordinates = 'cylindrical'
has_run = False
global_enc_mass = []
normal = [1.0, 0.0, 0.0]

def set_center(x):
    """
    Sets the center used when calculateing fields.
    
    Type: int
    Default: 0
    Options:0=center of mass, 1=particle 1, 2=particle 2.
    """
    global center
    global global_enc_mass
    center = x
    global_enc_mass = []
    return center

def set_n_bins(x):
    """
    sets the number of bins used when calculating the enclosed mass
    Type: int
    """
    global n_bins
    n_bins = int(x)
    return x

def set_adaptive_bins(x):
    """
    sets the whether adaptive bins are used when calculating used when calculating the enclosed mass
    Type: bool
    """
    global adaptive_bins
    adaptive_bins = x
    return x

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
    global global_enc_mass
    global coordinates
    coordinates = x
    global_enc_mass = []
    return coordinates

def set_normal(x):
    """
    Sets the normal used for projected fields.
        
    Default: [x,y,z] = [1,0,0]
    """
    global normal
    if isinstance(x,list) == False:
        x = x.tolist()
    normal = x
    return normal

def get_center():
    """
    returns the currently set center
    """
    global center
    return center

def get_n_bins():
    """
    returns the currently set number of bins
    """
    global n_bins
    return n_bins

def get_adaptive_bins():
    """
    returns the current set adaptive_bins bool
    """
    global adaptive_bins
    return adaptive_bins

def get_image_resolution():
    """
    returns the currently set image resolution
    """
    global image_resolution
    return image_resolution

def get_coordinate_system():
    """
    returns the currently set coordinate system
    """
    global coordinates
    return coordinates

def get_normal():
    """
    returns the currently set normal
    """
    global normal
    return normal

def _CoM(field, data):
    """
    Calculate the center of mass. Always includes particles where possible.
    """
    TM = np.sum(data['cell_mass'].in_units('g'))
    x_top = yt.YTArray(0.0, 'cm*g')
    y_top = yt.YTArray(0.0, 'cm*g')
    z_top = yt.YTArray(0.0, 'cm*g')
    if ('all', 'particle_mass') in data.ds.field_list:
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
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM", function=_CoM, units=r"cm")

def _Semimajor_Axis(field, data):
    """
    Calculates the semimajor axis of the binary. If there are no stars, or only one star, the semimajor axis is set to zero.
    """
    if ('all', 'particle_mass') in data.ds.field_list:
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
    if ('all', 'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(data['particle_mass'].in_units('g'))
        for part in range(len(data['particle_mass'])):
            x_vel = x_vel + data['particle_mass'][part].in_units('g')*data['particle_velx'][part].in_units('cm/s')
            y_vel = y_vel + data['particle_mass'][part].in_units('g')*data['particle_vely'][part].in_units('cm/s')
            z_vel = z_vel + data['particle_mass'][part].in_units('g')*data['particle_velz'][part].in_units('cm/s')
    bv = [(x_vel/TM), (y_vel/TM), (z_vel/TM)]
    bv = yt.YTArray(bv, 'cm/s')
    del x_vel
    del y_vel
    del z_vel
    del TM
    return bv

yt.add_field("My_Bulk_Velocity", function=_My_Bulk_Velocity, units=r"cm/s")

def _Center_Position(field, data):
    """
    Returns the center position for the current set center.
    """
    global center
    dd = data.ds.all_data()
    if np.shape(data)[0] == 16:
        center_pos = data['CoM'].in_units('cm')
    elif center == 0:
        center_pos = dd['CoM'].in_units('cm')
    else:
        center_pos = [dd['particle_posx'][center-1].in_units('cm').value, dd['particle_posy'][center-1].in_units('cm').value, dd['particle_posz'][center-1].in_units('cm').value]
        center_pos = yt.YTArray(center_pos, 'cm')
    del dd
    return center_pos

yt.add_field("Center_Position", function=_Center_Position, units=r"cm")

def _Center_Velocity(field, data):
    """
    Returns the center velocity for the current set center.
    """
    global center
    dd = data.ds.all_data()
    center_vel = yt.YTArray(np.array([0.0,0.0,0.0]), 'cm/s')
    if center == 0 and ('gas', 'velocity_x') in data.ds.derived_field_list:
        if ('all', 'particle_mass') in data.ds.field_list:
            center_vel = dd.quantities.bulk_velocity(use_particles=True)
        else:
            center_vel = dd.quantities.bulk_velocity(use_particles=False)
    elif center == 0 and ('gas', 'velocity_x') not in data.ds.derived_field_list:
        center_vel = yt.YTArray([np.sum(dd['velx'].in_units('cm/s').value), np.sum(dd['vely'].in_units('cm/s').value), np.sum(dd['velz'].in_units('cm/s').value)], 'cm/s')
    else:
        center_vel = yt.YTArray([dd['particle_velx'][center-1].in_units('cm/s').value, dd['particle_vely'][center-1].in_units('cm/s').value, dd['particle_velz'][center-1].in_units('cm/s').value], 'cm/s')
    del dd
    return center_vel

yt.add_field("Center_Velocity", function=_Center_Velocity, units=r"cm/s")

#Lets also create a new field with the distance of cells from the CoM
def _dx_from_Center(field, data):
    """
    Calculates the change in x position from the current set center.
    """
    if ('gas', 'Center_Position') in data.ds.derived_field_list:
        dd = data.ds.all_data()
        center_pos = dd['Center_Position'][0].in_units('cm')
    else:
        center_pos = yt.YTArray(0.0, 'cm')
    dx = data['x'].in_units('cm')-center_pos
    return dx

yt.add_field("dx_from_Center", function=_dx_from_Center, units=r"cm")

def _dy_from_Center(field, data):
    """
    Calculates the change in y position from the current set center.
    """
    if ('gas', 'Center_Position') in data.ds.derived_field_list:
        dd = data.ds.all_data()
        center_pos = dd['Center_Position'][1].in_units('cm')
    else:
        center_pos = yt.YTArray(0.0, 'cm')
    dy = data['y'].in_units('cm')-center_pos
    return dy

yt.add_field("dy_from_Center", function=_dy_from_Center, units=r"cm")

def _dz_from_Center(field, data):
    """
    Calculates the change in z position from the current set center.
    """
    if ('gas', 'Center_Position') in data.ds.derived_field_list:
        dd = data.ds.all_data()
        center_pos = dd['Center_Position'][2].in_units('cm')
    else:
        center_pos = yt.YTArray(0.0, 'cm')
    dz = data['z'].in_units('cm')-center_pos
    return dz

yt.add_field("dz_from_Center", function=_dz_from_Center, units=r"cm")

def _Distance_from_Center(field, data):
    """
    Calculates the distance from the current set center.
    """
    global coordinates
    if 'cyl' in coordinates.lower():
        distance = np.sqrt((data['dx_from_Center'])**2. + (data['dy_from_Center'])**2.)
    else:
        distance = np.sqrt((data['dx_from_Center'])**2. + (data['dy_from_Center'])**2. + (data['dz_from_Center'])**2.)
    return distance

yt.add_field("Distance_from_Center", function=_Distance_from_Center, units=r"cm")

def _Corrected_velx(field, data):
    """
    Calculates the x-velocity correcnted for the bulk velocity.
    """
    if ('gas', 'Center_Velocity') in data.ds.derived_field_list:
        dd = data.ds.all_data()
        center_vel = dd['Center_Velocity'][0].in_units('cm/s')
    else:
        center_vel = yt.YTArray(0.0, 'cm/s')
    dvx = data['velx'].in_units('cm/s') - center_vel
    return dvx

yt.add_field("Corrected_velx", function=_Corrected_velx, units=r"cm/s")

def _Corrected_vely(field, data):
    """
    Calculates the y-velocity correcnted for the bulk velocity.
    """
    if ('gas', 'Center_Velocity') in data.ds.derived_field_list:
        dd = data.ds.all_data()
        center_vel = dd['Center_Velocity'][1].in_units('cm/s')
    else:
        center_vel = yt.YTArray(0.0, 'cm/s')
    dvy = data['vely'].in_units('cm/s') - center_vel
    return dvy

yt.add_field("Corrected_vely", function=_Corrected_vely, units=r"cm/s")

def _Corrected_velz(field, data):
    """
    Calculates the z-velocity correcnted for the bulk velocity.
    """
    if ('gas', 'Center_Velocity') in data.ds.derived_field_list:
        dd = data.ds.all_data()
        center_vel = dd['Center_Velocity'][2].in_units('cm/s')
    else:
        center_vel = yt.YTArray(0.0, 'cm/s')
    dvz = data['velz'].in_units('cm/s') - center_vel
    return dvz

yt.add_field("Corrected_velz", function=_Corrected_velz, units=r"cm/s")

def _Corrected_vel_mag(field, data):
    """
    Calculates the velocity magnitude corrected for the bulk velocity in the current coordinate system.
    """
    global coordinates
    if 'cyl' in coordinates.lower():
        vmag = np.sqrt((data['Corrected_velx'])**2. + (data['Corrected_vely'])**2.)
    else:
        vmag = np.sqrt((data['Corrected_velx'])**2. + (data['Corrected_vely'])**2. + (data['Corrected_velz'])**2.)
    return vmag

yt.add_field("Corrected_vel_mag", function=_Corrected_vel_mag, units=r"cm/s")

def _Projected_Velocity(field, data):
    """
    If the normal is [1,0,0] and there are two particles this calculates the projected velocity of the gas in the plane of their separation axis. If normal is set to something else, it calculates the projected
    """
    global normal
    if np.shape(data)[0] != 16:
        if normal != [1.0,0.0,0.0]:
            pos_vec = np.array([normal[1], -1*normal[0]])
        elif ('all', 'particle_mass') in data.ds.field_list:
            dd = data.ds.all_data()
            if len(dd['particle_mass']) == 2:
                pos_vec = np.array([np.diff(dd[('all', 'particle_posx')].value)[0], np.diff(dd[('all', 'particle_posy')].value)[0]])
            else:
                pos_vec = np.array([0.0, -1.0])
        else:
            pos_vec = np.array([0.0, -1.0])
        import pdb
        pdb.set_trace()
        #pos_mag = np.sqrt(pos_vec[0]**2. + pos_vec[1]**2.)
        vels = np.array([data[('velx')].in_units('cm/s').value, data[('vely')].in_units('cm/s').value])
        vels = vels.T
        c = ((np.dot(vels,pos_vec))/(np.dot(pos_vec,pos_vec)))* pos_vec
        #velx = pos_vec[0] * c
        #vely = pos_vec[1] * c
        #del c
        #vels = np.sqrt(velx**2. + vely**2.)
        #del velx
        #del vely
        #vels = c*pos_mag
        vels = yt.YTArray(vels, 'cm/s')
    else:
        vels = data['velx']
    return vels

yt.add_field("Projected_Velocity", function=_Projected_Velocity, units=r"cm/s")

def _Projected_Magnetic_Field(field, data):
    """
    If the normal is [1,0,0] and there are two particles this calculates the projected magnetic field of the gas in the plane of their separation axis. If normal is set to something else, it calculates the projected
    """
    global normal
    if normal != [1.0,0.0,0.0]:
        pos_vec = np.array([normal[1], -1*normal[0]])
    elif ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        if len(dd['particle_mass']) == 2:
            pos_vec = np.array([np.diff(dd[('all', 'particle_posx')].value)[0], np.diff(dd[('all', 'particle_posy')].value)[0]])
        else:
            pos_vec = np.array([0.0, -1.0])
        del dd
    else:
        pos_vec = np.array([0.0, -1.0])
    #pos_mag = np.sqrt(pos_vec[0]**2. + pos_vec[1]**2.)
    mags = np.array([data['magx'].value, data['magy'].value])
    mags = mags.T
    c = ((np.dot(mags,pos_vec))/(np.dot(pos_vec,pos_vec)))
    magx = pos_vec[0] * c
    magy = pos_vec[1] * c
    del c
    mags = np.sqrt(magx**2. + magy**2.)
    del magx
    del magy
    mags = yt.YTArray(mags, 'gauss')
    return mags

yt.add_field("Projected_Magnetic_Field", function=_Projected_Magnetic_Field, units=r"gauss")
'''
def Enclosed_Mass(file, max_radius,inds):
    global center
    global n_bins
    global adaptive_bins
    global global_enc_mass
    part_file = part_file=file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    data = ds.all_data()
    if ('all', u'particle_mass') in data.ds.field_list:
        particle_mass = data['particle_mass']
        if len(data['particle_mass']) == 2:
            pos1 = [data['particle_posx'][0].in_units('cm'), data['particle_posy'][0].in_units('cm'), data['particle_posz'][0].in_units('cm')]
            pos2 = [data['particle_posx'][1].in_units('cm'), data['particle_posy'][1].in_units('cm'), data['particle_posz'][1].in_units('cm')]
            a = np.sqrt((pos1[0] - pos2[0])**2. + (pos1[1] - pos2[1])**2. + (pos1[2] - pos2[2])**2.)
        else:
            a = yt.YTArray(0.0, 'cm')
    else:
        particle_mass = yt.YTArray([], 'g')
        a = yt.YTArray(0.0, 'cm')
    
    distance = data['Distance_from_Center']
    cell_mass = data['cell_mass']
    time_str = str(int(time.time()))
    pickle_file_enc = '/short/ek9/rlk100/temp_pickles/enclosed_mass_temp_' + time_str + '.pkl'
    pickle_file = '/short/ek9/rlk100/temp_pickles/pickle_temp_' + time_str + '.pkl'
    temp_file = open(pickle_file, 'w')
    pickle.dump((distance.in_units('cm').value, cell_mass.in_units('g').value, particle_mass.in_units('g').value), temp_file)
    temp_file.close()
    call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/Modules/enclosed_mass.py', '-f', pickle_file, '-sf', pickle_file_enc, '-c', str(center), '-a', str(a.in_units('cm').value), '-mr', str(max_radius.in_units('cm').value), '-bins', str(int(n_bins)), '-ab', str(adaptive_bins)])
    #call(['mpirun', '-np', '8', 'python', '/Users/rajikak/Scripts/Modules/enclosed_mass.py', '-f', pickle_file, '-c', str(center), '-bs', str(n_bins), '-co', coordinates, '-a', str(a.in_units('cm').value), '-mr', str(max_radius.in_units('cm').value)])

    file = open(pickle_file_enc, 'r')
    enclosed_mass = pickle.load(file)
    enc_mass = yt.YTArray(enclosed_mass, 'g')
    global_enc_mass = enc_mass[inds]
    os.remove(pickle_file_enc)
    return global_enc_mass

def _Enclosed_Mass(field, data):
    """
    Calculates the enclosed mass for the set center and in the set coordinate system, with the current set bin size
    """
    global global_enc_mass
    
    if np.shape(data)[0] == 16:
        enclosed_mass = yt.YTArray(np.zeros(np.shape(data)), 'g')
    elif np.shape(data) != np.shape(global_enc_mass):
        file = data.ds.fullpath +'/'+data.ds.basename
        max_radius = np.max(data['Distance_from_Center'])
        dd = data.ds.all_data()
        comb = dd['x'].value + dd['y'].in_units('km').value + dd['z'].in_units('au').value
        data_comb = data['x'].value + data['y'].in_units('km').value + data['z'].in_units('au').value
        inds = np.where(np.in1d(comb, data_comb))[0]
        enclosed_mass = Enclosed_Mass(file, max_radius,inds)
    elif np.shape(data) == np.shape(global_enc_mass):
        enclosed_mass = global_enc_mass
    return enclosed_mass

yt.add_field("Enclosed_Mass", function=_Enclosed_Mass, units=r"g")
'''
def _Radial_Velocity(field, data):
    """
    Calculates the radial velocity from the current center, in the current coordinate system, corrected for the velocity of the center.
    """
    global coordinates
    if 'cyl' in coordinates.lower():
        rad_vel = (data['Corrected_velx'].in_units('cm/s')*data['dx_from_Center'].in_units('cm') + data['Corrected_vely'].in_units('cm/s')*data['dy_from_Center'].in_units('cm'))/data['Distance_from_Center'].in_units('cm')
    else:
        rad_vel = (data['Corrected_velx'].in_units('cm/s')*data['dx_from_Center'].in_units('cm') + data['Corrected_vely'].in_units('cm/s')*data['dy_from_Center'].in_units('cm') + data['Corrected_velz'].in_units('cm/s')*data['dz_from_Center'].in_units('cm'))/data['Distance_from_Center'].in_units('cm')
    return rad_vel

yt.add_field("Radial_Velocity", function=_Radial_Velocity, units=r"cm/s")

def _Tangential_Velocity(field, data):
    """
    Calculates the Tangetial velocity for the current center, in the current coordinate system, corrected for the velocity of the center.
    """
    v_mag = data['Corrected_vel_mag']
    v_r = data['Radial_Velocity']
    v_t = np.sqrt(v_mag**2. - v_r**2.)
    del v_mag
    del v_r
    return v_t

yt.add_field("Tangential_Velocity", function=_Tangential_Velocity, units=r"cm/s")

def _Total_Potential(field, data):
    """
        Gives the total potential inclusing contribution from the gas and the sink particles.
        """
    G_pot_total = data['gpot'] + data['Particle_Potential']
    return G_pot_total

yt.add_field("Total_Potential", function=_Total_Potential, units=r"cm**2/s**2")

def _Keplerian_Velocity(field, data):
    """
    Keplerian velocity calculated from the total potential energy (Sum of the potential from the sinks and the gas)
    """
    keplerian_field = np.sqrt(-1*data['Total_Potential'])
    return keplerian_field

yt.add_field("Keplerian_Velocity", function=_Keplerian_Velocity, units=r"cm/s")
'''
def _Keplerian_Velocity(field, data):
    """
    Calculates the keplerian velocity for the enclosed mass calculated from the current center, in the current coordinate system, corrected for the velocity of the center.
    """
    G = yt.utilities.physical_constants.G
    keplerian_field = np.sqrt((G*data['Enclosed_Mass'])/data['Distance_from_Center'])
    return keplerian_field

yt.add_field("Keplerian_Velocity", function=_Keplerian_Velocity, units=r"cm/s")
'''
def _Relative_Keplerian_Velocity(field, data):
    """
    Calculates the Relative Keplerian Velocity.
    """
    v_mag = data['Tangential_Velocity']
    #v_mag = data['Corrected_vel_mag']
    v_kep = data['Keplerian_Velocity']
    rel_kep = v_mag/v_kep
    del v_mag
    del v_kep
    return rel_kep

yt.add_field("Relative_Keplerian_Velocity", function=_Relative_Keplerian_Velocity, units=r"")
'''
def Gradient(x_field, y_field, bin_data):
    global n_bins
    global adaptive_bins
    
    print "X_FIELD =", x_field
    print "Y_FIELD =", y_field
    
    time_str = str(time.time()).split('.')[0] + '_' + str(time.time()).split('.')[1]
    gradient_input = '/short/ek9/rlk100/temp_pickles/gradient_input_' + time_str + '.pkl'
    gradient_output = '/short/ek9/rlk100/temp_pickles/gradient_' + time_str + '.pkl'
    temp_file = open(gradient_input, 'w')
    pickle.dump((x_field.value, y_field.value, bin_data.value), temp_file)
    temp_file.close()
    call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/Modules/gradient.py', '-f', gradient_input, '-sf', gradient_output, '-bins', str(n_bins), '-ab', str(adaptive_bins)])

    file = open(gradient_output, 'r')
    gradient = pickle.load(file)
    os.remove(gradient_output)
    return gradient

def _Squared_B_Mag(field, data):
    """
    Just calculates the dquare of the magnetic fields magnitude
    """
    if ('gas', 'magnetic_field_magnitude') in data.ds.derived_field_list:
        B = data['magnetic_field_magnitude']**2.
    else:
        B = yt.YTArray(np.zeros(np.shape(data)), 'G**2')
    return B

yt.add_field("Squared_B_Mag", function=_Squared_B_Mag, units=r"G**2")

def _B_gradient(field, data):
    """
    Calculates the magnetic field gradient in the z direction
    """
    if np.shape(data)[0] != 16:
        y = data['Squared_B_Mag']
        #bin_data = np.abs(data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2.))
        #x = np.abs(data['dz_from_Center'].in_units('cm'))
        bin_data = data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2)
        x = data['dz_from_Center'].in_units('cm')
        gradient = Gradient(x, y, bin_data)
        dB = yt.YTArray(np.abs(gradient), 'G**2/cm')
    else:
        dB = yt.YTArray(np.zeros(np.shape(data)), 'G**2/cm')
    return dB

yt.add_field("B_gradient", function=_B_gradient, units=r"G**2/cm")

def _Magnetic_Acceleration(field, data):
    """
    Calculates the magnetic pressure in sphere
    """
    magnetic_acceleration = 1/(data['dens'].in_units('g/cm**3')*8.*np.pi) * data['B_gradient']
    magnetic_acceleration = np.sign(data['dz_from_Center'])*(np.abs(magnetic_acceleration.in_units('cm/s**2')))
    return magnetic_acceleration

yt.add_field("Magnetic_Acceleration", function=_Magnetic_Acceleration, units=r"cm/s**2")

def _Gravitational_Acceleration(field, data):
    """
    Calculates gravitational pressure, from gravitational force.
    """
    gravitational_acceleration = -1*np.sign(data['dz_from_Center'])*(yt.physical_constants.G*data['Enclosed_Mass'].in_units('g'))/(data['dz_from_Center'].in_units('cm')**2.)
    return gravitational_acceleration

yt.add_field("Gravitational_Acceleration", function=_Gravitational_Acceleration, units=r"cm/s**2")
'''

def _Gravitational_Acceleration(field, data):
    """
    Calculates gravitational pressure, from gravitational force.
    """
    gravitational_acceleration = (data['Total_Potential']*data['cell_mass'])/data['Distance_from_Center']
    return gravitational_acceleration

yt.add_field("Gravitational_Acceleration", function=_Gravitational_Acceleration, units=r"cm/s**2")

def _Particle_Potential(field, data):
    """
    Calculates the potential from the praticles.
    """
    if np.shape(data)[0] == 16:
        Part_gpot = yt.YTArray(np.zeros(np.shape(data)), 'cm**2/s**2')
    else:
        dd = data.ds.all_data()
        comb = dd['x'].value + dd['y'].in_units('km').value + dd['z'].in_units('au').value
        data_comb = data['x'].value + data['y'].in_units('km').value + data['z'].in_units('au').value
        inds = np.where(np.in1d(comb, data_comb))[0]
        Part_gpot = yt.YTArray(np.zeros(np.shape(dd)), 'cm**2/s**2')
        if ('all', 'particle_mass') in data.ds.field_list:
            for part in range(len(dd['particle_mass'])):
                gpot = -(yt.physical_constants.G*dd['particle_mass'][part].in_units('g'))/(dd['Distance_from_Center'].in_units('cm'))
                Part_gpot = Part_gpot + gpot
        Part_gpot = Part_gpot[inds]
    return Part_gpot

yt.add_field("Particle_Potential", function=_Particle_Potential, units=r"cm**2/s**2")
'''
def _Gravitational_Acceleration_z(field, data):
    """
    Calculates gravitational pressure, from gravitational force.
    """
    y = data['Total_Potential']
    if np.shape(data)[0] != 16:
        #bin_data = np.abs(data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2.))
        #x = np.abs(data['dz_from_Center'].in_units('cm'))
        bin_data = data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2.)
        x = data['dz_from_Center'].in_units('cm')
        gradient = Gradient(x, y, bin_data)
        gravitational_acceleration = yt.YTArray(-1*np.sign(data['dz_from_Center'])*np.abs(gradient), 'cm/s**2')
    else:
        gravitational_acceleration = yt.YTArray(np.zeros(np.shape(data)), 'cm/s**2')
    del y
    return gravitational_acceleration

yt.add_field("Gravitational_Acceleration_z", function=_Gravitational_Acceleration_z, units=r"cm/s**2")

def _Acceleration_Ratio(field, data):
    """
    Ratio of the magnetic acceleration to the gravitational acceleration calculated from enclosed_mass
    """
    A_ratio = np.abs(data['Magnetic_Acceleration']/data['Gravitational_Acceleration'])
    return A_ratio

yt.add_field("Acceleration_Ratio", function=_Acceleration_Ratio, units=r"")

def _Acceleration_Ratio_z(field, data):
    """
    Ratio of the magnetic acceleration to the gravitational acceleration calculated from the potential gradient
    """
    A_ratio = np.abs(data['Magnetic_Acceleration']/data['Gravitational_Acceleration_z'])
    return A_ratio

yt.add_field("Acceleration_Ratio_z", function=_Acceleration_Ratio_z, units=r"")

def _Magnetic_Force(field, data):
    """
    Calculates force from magnetic gradient
    """
    F_mag = data['cell_mass']*data['Magnetic_Acceleration']
    return F_mag

yt.add_field("Magnetic_Force", function=_Magnetic_Force, units=r"g*cm/s**2")

def _Gravitational_Force(field, data):
    """
    Calculates force from gravity
    """
    F_grav = np.abs(data['cell_mass']*data['Gravitational_Acceleration'])
    return F_grav

yt.add_field("Gravitational_Force", function=_Gravitational_Force, units=r"g*cm/s**2")

def _Gravitational_Force_z(field, data):
    """
    Calculates force from gravity
    """
    F_grav = np.abs(data['cell_mass']*data['Gravitational_Acceleration_z'])
    return F_grav

yt.add_field("Gravitational_Force_z", function=_Gravitational_Force_z, units=r"g*cm/s**2")

def _Force_Ratio(field, data):
    """
    Calculates the ratio of the magnetic pressure to the gravitational pressure
    """
    force_ratio = np.abs((data['Magnetic_Force']))/(np.abs(data['Gravitational_Force']))
    return force_ratio

yt.add_field("Force_Ratio", function=_Force_Ratio, units=r"")

def _Force_Ratio_z(field, data):
    """
    Calculates the ratio of the magnetic pressure to the gravitational pressure
    """
    force_ratio = np.abs((data['Magnetic_Force']))/(np.abs(data['Gravitational_Force_z']))
    return force_ratio

yt.add_field("Force_Ratio_z", function=_Force_Ratio_z, units=r"")
'''

def _Angular_Momentum_x(field, data):
    """
    Calculates the angular momentum in the x_direction about current set center.
    """
    L_x = data['cell_mass']*(data['Corrected_velx']*data['dy_from_Center']- data['Corrected_vely']*data['dz_from_Center'])
    return L_x

yt.add_field("Angular_Momentum_x", function=_Angular_Momentum_x, units=r"g*cm**2/s")

def _Angular_Momentum_y(field, data):
    """
    Calculates the angular momentum in the y_direction about current set center.
    """
    L_y = data['cell_mass']*(data['Corrected_velz']*data['dx_from_Center'] - data['Corrected_velx']*data['dz_from_Center'])
    return L_y

yt.add_field("Angular_Momentum_y", function=_Angular_Momentum_y, units=r"g*cm**2/s")

def _Angular_Momentum_z(field, data):
    """
    Calculates the angular momentum in the z_direction about current set center.
    """
    L_z = data['cell_mass']*(data['Corrected_vely']*data['dx_from_Center'] - data['Corrected_velx']*data['dy_from_Center'])
    return L_z

yt.add_field("Angular_Momentum_z", function=_Angular_Momentum_z, units=r"g*cm**2/s")

def _Angular_Momentum(field, data):
    """
    Calculates the angular momentum about current set center.
    """
    L = np.sqrt(data['Angular_Momentum_x']**2. + data['Angular_Momentum_y']**2. + data['Angular_Momentum_z']**2.)
    return L

yt.add_field("Angular_Momentum", function=_Angular_Momentum, units=r"g*cm**2/s")

def _Specific_Angular_Momentum(field, data):
    """
    Calculates the specific angular momentum about current set center.
    """
    l = data['Angular_Momentum']/data['cell_mass']
    return l

yt.add_field("Specific_Angular_Momentum", function=_Specific_Angular_Momentum, units=r"cm**2/s")

def _Particle_dx_from_Center(field, data):
    """
    Calculates the change in x position from the current set center.
    """
    if ('all', 'particle_posx') in data.ds.field_list:
        dx = data['particle_posx'].in_units('cm')-data['Center_Position'][0]
    else:
        dx = yt.YTArray(np.zeros(np.shape(data)), 'cm')
    return dx

yt.add_field("Particle_dx_from_Center", function=_Particle_dx_from_Center, units=r"cm")

def _Particle_dy_from_Center(field, data):
    """
    Calculates the change in y position from the current set center.
    """
    if ('all', 'particle_posy') in data.ds.field_list:
        dy = data['particle_posy'].in_units('cm')-data['Center_Position'][1]
    else:
        dy = yt.YTArray(np.zeros(np.shape(data)), 'cm')
    return dy

yt.add_field("Particle_dy_from_Center", function=_Particle_dy_from_Center, units=r"cm")

def _Particle_dz_from_Center(field, data):
    """
    Calculates the change in z position from the current set center.
    """
    if ('all', 'particle_posz') in data.ds.field_list:
        dz = data['particle_posz'].in_units('cm')-data['Center_Position'][2]
    else:
        dz = yt.YTArray(np.zeros(np.shape(data)), 'cm')
    return dz

yt.add_field("Particle_dz_from_Center", function=_Particle_dz_from_Center, units=r"cm")

def _Particle_Distance_from_Center(field, data):
    """
    Calculates the distance from the current set center.
    """
    global coordinates
    if 'cyl' in coordinates.lower():
        distance = np.sqrt((data['Particle_dx_from_Center'])**2. + (data['Particle_dy_from_Center'])**2.)
    else:
        distance = np.sqrt((data['Particle_dx_from_Center'])**2. + (data['Particle_dy_from_Center'])**2. + (data['Particle_dz_from_Center'])**2.)
    return distance

yt.add_field("Particle_Distance_from_Center", function=_Particle_Distance_from_Center, units=r"cm")

def _Particle_Angular_Momentum_x(field, data):
    """
    Calculates the angular momentum in the x_direction about current set center.
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        L_x = data['particle_mass']*(data['particle_vely']*data['Particle_dz_from_Center'] - data['particle_velz']*data['Particle_dy_from_Center'])
    else:
        L_x = yt.YTArray(np.zeros(np.shape(data)), 'g*cm**2/s')
    return L_x

yt.add_field("Particle_Angular_Momentum_x", function=_Particle_Angular_Momentum_x, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_y(field, data):
    """
    Calculates the angular momentum in the y_direction about current set center.
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        L_y = data['particle_mass']*(data['particle_velx']*data['Particle_dz_from_Center'] - data['particle_velz']*data['Particle_dx_from_Center'])
    else:
        L_y = yt.YTArray(np.zeros(np.shape(data)), 'g*cm**2/s')
    return L_y

yt.add_field("Particle_Angular_Momentum_y", function=_Particle_Angular_Momentum_y, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_z(field, data):
    """
    Calculates the angular momentum in the z_direction about current set center.
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        L_z = data['particle_mass']*(data['particle_velx']*data['Particle_dy_from_Center'] - data['particle_vely']*data['Particle_dx_from_Center'])
    else:
        L_z = yt.YTArray(np.zeros(np.shape(data)), 'g*cm**2/s')
    return L_z

yt.add_field("Particle_Angular_Momentum_z", function=_Particle_Angular_Momentum_z, units=r"g*cm**2/s")

def _Particle_Angular_Momentum(field, data):
    """
    Calculates the angular momentum about current set center.
    """
    
    L = np.sqrt(data['Particle_Angular_Momentum_x']**2. + data['Particle_Angular_Momentum_y']**2. + data['Particle_Angular_Momentum_z']**2.)
    return L

yt.add_field("Particle_Angular_Momentum", function=_Particle_Angular_Momentum, units=r"g*cm**2/s")

def _Particle_Specific_Angular_Momentum(field, data):
    """
    Calculates the specific angular momentum about current set center.
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        l = data['Particle_Angular_Momentum']/data['particle_mass']
    else:
        l = yt.YTArray(np.zeros(np.shape(data)), 'cm**2/s')
    return l

yt.add_field("Particle_Specific_Angular_Momentum", function=_Particle_Specific_Angular_Momentum, units=r"cm**2/s")
'''
def _B_angle(field, data):
    """
    Calculates the angle of the magnetic field from the xy-plane
    """
    if ('flash', u'magz') in data.ds.field_list:
        B_z = data['magnetic_field_z']
    else:
        B_z = yt.YTArray(np.ones(np.shape(data['cell_mass'])), 'gauss')
    B_h = np.sqrt(data['magnetic_field_x']**2. + data['magnetic_field_y']**2.)
    angle = np.arctan(B_z/B_h)
    angle = np.rad2deg(angle)
    angle = yt.YTArray(angle, 'deg')
    del B_h
    return angle

yt.add_field("B_angle", function=_B_angle, units=r"deg")
'''
def _Projected_Velocity_mw(field, data):
    """
    field to be able to created a mass weighted projection of Projected_Velocity
    """
    return data['Projected_Velocity']*data['cell_mass']

yt.add_field("Projected_Velocity_mw", function=_Projected_Velocity_mw, units=r"cm*g/s")

def _velx_mw(field, data):
    """
    field to be able to created a mass weighted projection of velx
    """
    return data['velx']*data['cell_mass']

yt.add_field("velx_mw", function=_velx_mw, units=r"cm*g/s")

def _vely_mw(field, data):
    """
    field to be able to created a mass weighted projection of vely
    """
    return data['vely']*data['cell_mass']

yt.add_field("vely_mw", function=_vely_mw, units=r"cm*g/s")

def _velz_mw(field, data):
    """
    field to be able to created a mass weighted projection of velz
    """
    return data['velz']*data['cell_mass']

yt.add_field("velz_mw", function=_velz_mw, units=r"cm*g/s")

def _Projected_Magnetic_Field_mw(field, data):
    """
    field to be able to created a mass weighted projection of Projected_Magnetic_Field
    """
    return data['Projected_Magnetic_Field']*data['cell_mass']

yt.add_field("Projected_Magnetic_Field_mw", function=_Projected_Magnetic_Field_mw, units=r"gauss*g")

def _magx_mw(field, data):
    """
    field to be able to created a mass weighted projection of magz
    """
    if ('flash', 'magx') in data.ds.field_list:
        return data['magx']*data['cell_mass']
    else:
        return yt.YTArray(np.ones(np.shape(data['cell_mass'])), 'gauss*g')

yt.add_field("magx_mw", function=_magx_mw, units=r"gauss*g")

def _magy_mw(field, data):
    """
    field to be able to created a mass weighted projection of magz
    """
    if ('flash', 'magy') in data.ds.field_list:
        return data['magy']*data['cell_mass']
    else:
        return yt.YTArray(np.ones(np.shape(data['cell_mass'])), 'gauss*g')

yt.add_field("magy_mw", function=_magy_mw, units=r"gauss*g")

def _magz_mw(field, data):
    """
    field to be able to created a mass weighted projection of magz
    """
    if ('flash', 'magz') in data.ds.field_list:
        return data['magz']*data['cell_mass']
    else:
        return yt.YTArray(np.ones(np.shape(data['cell_mass'])), 'gauss*g')

yt.add_field("magz_mw", function=_magz_mw, units=r"gauss*g")

'''
def _Is_Unbound(field, data):
    """
    returned boolean array about whether the gas is bound or unbound using the total potential energy and the kinetic energy
    """
    Total_Energy = (data['Total_Potential']*data['cell_mass'] + data['kinetic_energy']*data['cell_volume']).in_units('erg').value
    unbound_inds = np.where(Total_Energy > 0.0)[0]
    bound_inds = np.where(Total_Energy < 0.0)[0]
    Unbound = Total_Energy
    Unbound[unbound_inds] = True
    Unbound[bound_inds] = False
    del Total_Energy
    del unbound_inds
    del bound_inds
    return Unbound

yt.add_field("Is_Unbound", function=_Is_Unbound, units=r"")
'''
