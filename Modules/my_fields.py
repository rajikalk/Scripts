#!/usr/bin/env python
import yt
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

def set_n_bins(x):
    """
    sets the number of bins used when calculating the enclosed mass
    Type: int
    """
    global n_bins
    n_bins = x
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
    global coordinates
    coordinates = x
    return coordinates

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

def set_adaptive_bins():
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
    dd = data.ds.all_data()
    if np.shape(data) == (16, 16, 16):
        center_pos = data['CoM'].in_units('cm')
    elif center == 0:
        center_pos = dd.quantities.center_of_mass(use_particles=True).in_units('cm')
    else:
        #center_pos = [data['particle_posx'][center-1].in_units('cm').value, data['particle_posy'][center-1].in_units('cm').value, data['particle_posz'][center-1].in_units('cm').value]
        center_pos = [dd['particle_posx'][center-1].in_units('cm').value, dd['particle_posy'][center-1].in_units('cm').value, dd['particle_posz'][center-1].in_units('cm').value]
        center_pos = yt.YTArray(center_pos, 'cm')
    return center_pos

yt.add_field("Center_Position", function=_Center_Position, units=r"cm")

def _Center_Velocity(field, data):
    """
    Returns the center velocity for the current set center.
    """
    global center
    dd = data.ds.all_data()
    if np.shape(data) == (16, 16, 16):
        center_vel = data['My_Bulk_Velocity'].in_units('cm/s')
    elif center == 0:
        center_vel = dd.quantities.bulk_velocity(use_particles=True)
    else:
        #center_vel = [data['particle_velx'][center-1].in_units('cm/s'), data['particle_vely'][center-1].in_units('cm/s'), data['particle_velz'][center-1].in_units('cm/s')]
        center_vel = [dd['particle_velx'][center-1].in_units('cm/s'), dd['particle_vely'][center-1].in_units('cm/s'), dd['particle_velz'][center-1].in_units('cm/s')]
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
    Calculates the velocity magnitude corrected for the bulk velocity in the current coordinate system.
    """
    global coordinates
    if 'cyl' in coordinates.lower():
        vmag = np.sqrt((data['Corrected_velx'])**2. + (data['Corrected_vely'])**2.)
    else:
        vmag = np.sqrt((data['Corrected_velx'])**2. + (data['Corrected_vely'])**2. + (data['Corrected_velz'])**2.)
    return vmag

yt.add_field("Corrected_vel_mag", function=_Corrected_vel_mag, units=r"cm/s")

def Enclosed_Mass(file, max_radius):
    global center
    global n_bins
    global adaptive_bins
    global coordinates
    global has_run
    global global_enc_mass
    if has_run == False:
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
        
        prev_coordinates = coordinates
        coordinates = 'cylindrical'
        distance = data['Distance_from_Center']
        cell_mass = data['cell_mass']
        time_str = str(int(time.time()))
        pickle_file_enc = '/short/ek9/rlk100/temp_pickles/enclosed_mass_temp_' + time_str + '.pkl'
        pickle_file = '/short/ek9/rlk100/temp_pickles/pickle_temp_' + time_str + '.pkl'
        temp_file = open(pickle_file, 'w')
        pickle.dump((distance.in_units('cm').value, cell_mass.in_units('g').value, particle_mass.in_units('g').value), temp_file)
        temp_file.close()
        call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/Modules/enclosed_mass.py', '-f', pickle_file, '-sf', pickle_file_enc, '-c', str(center), '-co', coordinates, '-a', str(a.in_units('cm').value), '-mr', str(max_radius.in_units('cm').value), '-bins', str(int(n_bins)), '-ab', str(adaptive_bins)])
        #call(['mpirun', '-np', '8', 'python', '/Users/rajikak/Scripts/Modules/enclosed_mass.py', '-f', pickle_file, '-c', str(center), '-bs', str(n_bins), '-co', coordinates, '-a', str(a.in_units('cm').value), '-mr', str(max_radius.in_units('cm').value)])

        file = open(pickle_file_enc, 'r')
        enclosed_mass = pickle.load(file)
        global_enc_mass = yt.YTArray(enclosed_mass, 'g')
        os.remove(pickle_file_enc)
        has_run = True
        coordinates = prev_coordinates
    return global_enc_mass

def _Enclosed_Mass(field, data):
    """
    Calculates the enclosed mass for the set center and in the set coordinate system, with the current set bin size
    """
    global coordinates
    global global_enc_mass
    global has_run
    
    import pdb
    pdb.set_trace()
    if np.shape(data) != np.shape(global_enc_mass):
        has_run = False

    if np.shape(data) != (16, 16, 16):
        file = data.ds.fullpath +'/'+data.ds.basename
        max_radius = np.max(data['Distance_from_Center'])
        enclosed_mass = Enclosed_Mass(file, max_radius)
        dd = data.ds.all_data()
        comb = dd['x'].value + dd['y'].in_units('km').value + dd['z'].in_units('au').value
        data_comb = data['x'].value + data['y'].in_units('km').value + data['z'].in_units('au').value
        inds = np.where(np.in1d(comb, data_comb))[0]
        if len(global_enc_mass) == 0:
            enclosed_mass = yt.YTArray(np.zeros(np.shape(data)), 'g')
        else:
            enclosed_mass = global_enc_mass[inds]
    else:
        enclosed_mass = yt.YTArray(np.zeros(np.shape(data)), 'g')
    return enclosed_mass

yt.add_field("Enclosed_Mass", function=_Enclosed_Mass, units=r"g")

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
    return v_t

yt.add_field("Tangential_Velocity", function=_Tangential_Velocity, units=r"cm/s")

def _Keplerian_Velocity(field, data):
    """
    Calculates the keplerian velocity for the enclosed mass calculated from the current center, in the current coordinate system, corrected for the velocity of the center.
    """
    G = yt.utilities.physical_constants.G
    keplerian_field = np.sqrt((G*data['Enclosed_Mass'])/data['Distance_from_Center'])
    return keplerian_field

yt.add_field("Keplerian_Velocity", function=_Keplerian_Velocity, units=r"cm/s")

def _Relative_Keplerian_Velocity(field, data):
    """
    Calculates the Relative Keplerian Velocity.
    """
    v_mag = data['Tangential_Velocity']
    v_kep = data['Keplerian_Velocity']
    rel_kep = v_mag/v_kep
    return rel_kep

yt.add_field("Relative_Keplerian_Velocity", function=_Relative_Keplerian_Velocity, units=r"")

def delta_B(file, data, field):
    global n_bins
    
    part_file = part_file=file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    dd = ds.all_data()
    
    z = data['dz_from_Center']
    comb = dd['x'].value + dd['y'].in_units('km').value + dd['z'].in_units('au').value
    data_comb = data['x'].value + data['y'].in_units('km').value + data['z'].in_units('au').value
    inds = np.where(np.in1d(comb, data_comb))[0]
    B = dd[field][inds]
    print "SHAPE(z) =", np.shape(z)
    print "SHAPE(B) =", np.shape(B)
    time_str = str(int(time.time()))
    pickle_file_B_grad = '/short/ek9/rlk100/temp_pickles/B_grad_temp_' + time_str + '.pkl'
    pickle_file = '/short/ek9/rlk100/temp_pickles/pickle_temp_B_grad_' + time_str + '.pkl'
    '''
    pickle_file_B_grad = '/Users/rajikak/Scripts/Modules/B_grad_temp.pkl'
    pickle_file = '/Users/rajikak/Scripts/Modules/pickle_temp_B_grad.pkl'
    '''
    file = open(pickle_file, 'w+')
    pickle.dump((np.abs(z.in_units('cm').value), B.value), file)
    file.close()
    
    call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/Modules/magnetic_gradient_z.py', '-f', pickle_file, '-sf', pickle_file_B_grad, '-bs', str(n_bins)])
    
    file = open(pickle_file_B_grad, 'r')
    B_grad = pickle.load(file)
    os.remove(pickle_file_B_grad)
    return B_grad

def _Squared_B_Mag(field, data):
    """
    Just calculates the dquare of the magnetic fields magnitude
    """
    B = data['magnetic_field_magnitude']**2.
    return B

yt.add_field("Squared_B_Mag", function=_Squared_B_Mag, units=r"G**2")

    
def _B_gradient(field, data):
    """
    Calculates the magnetic field gradient in the z direction
    """
    if np.shape(data) != (16, 16, 16):
        file = data.ds.fullpath +'/'+data.ds.basename
        dB = delta_B(file, data, 'Squared_B_Mag')
        dB = yt.YTArray(np.abs(dB), 'G**2')
        dB = dB/data['dz'].in_units('cm')
    else:
        dB = yt.YTArray(np.zeros(np.shape(data)), 'G**2/cm')
    return dB

yt.add_field("B_gradient", function=_B_gradient, units=r"G**2/cm")

def _Magnetic_Acceleration(field, data):
    """
    Calculates the magnetic pressure in sphere
    """
    magnetic_acceleration = 1/(data['dens'].in_units('g/cm**3')*8.*np.pi) * data['B_gradient']
    magnetic_acceleration = -1*magnetic_acceleration.in_units('cm/s**2')
    return magnetic_acceleration

yt.add_field("Magnetic_Acceleration", function=_Magnetic_Acceleration, units=r"cm/s**2")

def _Gravitational_Acceleration(field, data):
    """
    Calculates gravitational pressure, from gravitational force.
    """
    gravitational_acceleration = (yt.physical_constants.G*data['Enclosed_Mass'].in_units('g'))/(data['z'].in_units('cm')**2.)
    return gravitational_acceleration

yt.add_field("Gravitational_Acceleration", function=_Gravitational_Acceleration, units=r"cm/s**2")

def _Gravitational_Acceleration_z(field, data):
    """
    Calculates gravitational pressure, from gravitational force.
    """
    if np.shape(data) != (16, 16, 16):
        file = data.ds.fullpath +'/'+data.ds.basename
        delta_pot = delta_B(file, data, 'gravitational_potential')
        g_acc = delta_pot/dd['dz'].in_units('cm').value
        gravitational_acceleration = yt.YTArray(np.abs(g_acc), 'cm/s**2')
    else:
        gravitational_acceleration = yt.YTArray(np.zeros(np.shape(data)), 'cm/s**2')
    return gravitational_acceleration

yt.add_field("Gravitational_Acceleration_z", function=_Gravitational_Acceleration_z, units=r"cm/s**2")

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
    F_grav = -1*data['cell_mass']*data['Gravitational_Acceleration']
    return F_grav

yt.add_field("Gravitational_Force", function=_Gravitational_Force, units=r"g*cm/s**2")

def _Force_Ratio(field, data):
    """
    Calculates the ratio of the magnetic pressure to the gravitational pressure
    """
    global coordinates
    coordinates = 'sph'
    force_ratio = (data['Magnetic_Force'])/(-1*data['Gravitational_Force'])
    return force_ratio

yt.add_field("Force_Ratio", function=_Force_Ratio, units=r"")

def _Angular_Momentum_x(field, data):
    """
    Calculates the angular momentum in the x_direction about current set center.
    """
    L_x = data['cell_mass']*(data['vely']*data['dz_from_Center'] - data['velz']*data['dy_from_Center'])
    return L_x

yt.add_field("Angular_Momentum_x", function=_Angular_Momentum_x, units=r"g*cm**2/s")

def _Angular_Momentum_y(field, data):
    """
    Calculates the angular momentum in the y_direction about current set center.
    """
    L_y = data['cell_mass']*(data['velx']*data['dz_from_Center'] - data['velz']*data['dx_from_Center'])
    return L_y

yt.add_field("Angular_Momentum_y", function=_Angular_Momentum_y, units=r"g*cm**2/s")

def _Angular_Momentum_z(field, data):
    """
    Calculates the angular momentum in the z_direction about current set center.
    """
    L_z = data['cell_mass']*(data['velx']*data['dy_from_Center'] - data['vely']*data['dx_from_Center'])
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
    dx = data['particle_posx'].in_units('cm')-data['Center_Position'][0]
    return dx

yt.add_field("Particle_dx_from_Center", function=_Particle_dx_from_Center, units=r"cm")

def _Particle_dy_from_Center(field, data):
    """
        Calculates the change in y position from the current set center.
        """
    dy = data['particle_posy'].in_units('cm')-data['Center_Position'][1]
    return dy

yt.add_field("Particle_dy_from_Center", function=_Particle_dy_from_Center, units=r"cm")

def _Particle_dz_from_Center(field, data):
    """
        Calculates the change in z position from the current set center.
        """
    dz = data['particle_posz'].in_units('cm')-data['Center_Position'][2]
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
    L_x = data['particle_mass']*(data['particle_vely']*data['Particle_dz_from_Center'] - data['particle_velz']*data['Particle_dy_from_Center'])
    return L_x

yt.add_field("Particle_Angular_Momentum_x", function=_Particle_Angular_Momentum_x, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_y(field, data):
    """
        Calculates the angular momentum in the y_direction about current set center.
        """
    L_y = data['particle_mass']*(data['particle_velx']*data['Particle_dz_from_Center'] - data['particle_velz']*data['Particle_dx_from_Center'])
    return L_y

yt.add_field("Particle_Angular_Momentum_y", function=_Particle_Angular_Momentum_y, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_z(field, data):
    """
        Calculates the angular momentum in the z_direction about current set center.
        """
    L_z = data['particle_mass']*(data['particle_velx']*data['Particle_dy_from_Center'] - data['particle_vely']*data['Particle_dx_from_Center'])
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
    l = data['Particle_Angular_Momentum']/data['particle_mass']
    return l

yt.add_field("Particle_Specific_Angular_Momentum", function=_Particle_Specific_Angular_Momentum, units=r"cm**2/s")
