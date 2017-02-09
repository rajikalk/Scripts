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
    global global_enc_mass
    global coordinates
    coordinates = x
    global_enc_mass = []
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
    
    if np.shape(data) == (16, 16, 16):
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

def Gradient(x_field, y_field, max_radius, bin_data):
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
    call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/Modules/gradient.py', '-f', gradient_input, '-sf', gradient_output, '-mr', str(max_radius.value), '-bins', str(int(n_bins)), '-ab', str(adaptive_bins)])

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
    if np.shape(data) != (16, 16, 16):
        y = data['Squared_B_Mag']
        bin_data = np.abs(data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2.))
        x = np.abs(data['dz_from_Center'].in_units('cm'))
        max_radius = np.max(x)
        gradient = Gradient(x, y, max_radius, bin_data)
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

def _Particle_Potential(field, data):
    """
    Calculates the potential from the praticles.
    """
    if np.shape(data) == (16, 16, 16):
        Part_gpot = yt.YTArray(np.zeros(np.shape(data)), 'cm**2/s**2')
    else:
        dd = data.ds.all_data()
        comb = dd['x'].value + dd['y'].in_units('km').value + dd['z'].in_units('au').value
        data_comb = data['x'].value + data['y'].in_units('km').value + data['z'].in_units('au').value
        inds = np.where(np.in1d(comb, data_comb))[0]
        Part_gpot = yt.YTArray(np.zeros(np.shape(dd)), 'cm**2/s**2')
        if ('all', u'particle_mass') in data.ds.field_list:
            for part in range(len(dd['particle_mass'])):
                gpot = -(yt.physical_constants.G*dd['particle_mass'][part].in_units('g'))/(dd['Distance_from_Center'].in_units('cm'))
                Part_gpot = Part_gpot + gpot
        Part_gpot = Part_gpot[inds]
    return Part_gpot

yt.add_field("Particle_Potential", function=_Particle_Potential, units=r"cm**2/s**2")

def _Total_Potential(field, data):
    """
    Gives the total potential inclusing contribution from the gas and the sink particles.
    """
    G_pot_total = data['gpot'] + data['Particle_Potential']
    return G_pot_total

yt.add_field("Total_Potential", function=_Total_Potential, units=r"cm**2/s**2")

def _Gravitational_Acceleration_z(field, data):
    """
    Calculates gravitational pressure, from gravitational force.
    """
    y = data['Total_Potential']
    if np.shape(data) != (16, 16, 16):
        #bin_data = np.abs(data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2.))
        #x = np.abs(data['dz_from_Center'].in_units('cm'))
        bin_data = data['dz_from_Center'].in_units('cm') - (data['dz'].in_units('cm')/2.)
        x = data['dz_from_Center'].in_units('cm')
        max_radius = np.max(x)
        gradient = Gradient(x, y, max_radius, bin_data)
        gravitational_acceleration = yt.YTArray(-1*np.sign(data['dz_from_Center'])*np.abs(gradient), 'cm/s**2')
    else:
        gravitational_acceleration = yt.YTArray(np.zeros(np.shape(data)), 'cm/s**2')
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
    if ('all', u'particle_posx') in data.ds.field_list:
        dx = data['particle_posx'].in_units('cm')-data['Center_Position'][0]
    else:
        dx = yt.YTArray(np.zeros(np.shape(data)), 'cm')
    return dx

yt.add_field("Particle_dx_from_Center", function=_Particle_dx_from_Center, units=r"cm")

def _Particle_dy_from_Center(field, data):
    """
    Calculates the change in y position from the current set center.
    """
    if ('all', u'particle_posy') in data.ds.field_list:
        dy = data['particle_posy'].in_units('cm')-data['Center_Position'][1]
    else:
        dy = yt.YTArray(np.zeros(np.shape(data)), 'cm')
    return dy

yt.add_field("Particle_dy_from_Center", function=_Particle_dy_from_Center, units=r"cm")

def _Particle_dz_from_Center(field, data):
    """
    Calculates the change in z position from the current set center.
    """
    if ('all', u'particle_posz') in data.ds.field_list:
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
    if ('all', u'particle_mass') in data.ds.field_list:
        L_x = data['particle_mass']*(data['particle_vely']*data['Particle_dz_from_Center'] - data['particle_velz']*data['Particle_dy_from_Center'])
    else:
        L_x = yt.YTArray(np.zeros(np.shape(data)), 'g*cm**2/s')
    return L_x

yt.add_field("Particle_Angular_Momentum_x", function=_Particle_Angular_Momentum_x, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_y(field, data):
    """
    Calculates the angular momentum in the y_direction about current set center.
    """
    if ('all', u'particle_mass') in data.ds.field_list:
        L_y = data['particle_mass']*(data['particle_velx']*data['Particle_dz_from_Center'] - data['particle_velz']*data['Particle_dx_from_Center'])
    else:
        L_y = yt.YTArray(np.zeros(np.shape(data)), 'g*cm**2/s')
    return L_y

yt.add_field("Particle_Angular_Momentum_y", function=_Particle_Angular_Momentum_y, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_z(field, data):
    """
    Calculates the angular momentum in the z_direction about current set center.
    """
    if ('all', u'particle_mass') in data.ds.field_list:
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
    if ('all', u'particle_mass') in data.ds.field_list:
        l = data['Particle_Angular_Momentum']/data['particle_mass']
    else:
        l = yt.YTArray(np.zeros(np.shape(data)), 'cm**2/s')
    return l

yt.add_field("Particle_Specific_Angular_Momentum", function=_Particle_Specific_Angular_Momentum, units=r"cm**2/s")
