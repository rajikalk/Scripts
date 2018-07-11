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
coordinates = 'spherical'
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
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM", function=_CoM, units=r"cm")

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
        if ('all', u'particle_mass') in data.ds.field_list:
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
    del center_pos
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
    del center_pos
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
    del center_pos
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
    velx = data['velx'].in_units('cm/s')
    dvx = velx - center_vel
    del velx
    del center_vel
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
    vely = data['vely'].in_units('cm/s')
    dvy = vely - center_vel
    del vely
    del center_vel
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
    velz = data['velz'].in_units('cm/s')
    dvz = velz - center_vel
    del velz
    del center_vel
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
    Calculates the projected velocity on the plane perpendicular to the normal vector L.
    """
    global normal
    #velx = data[('Corrected_velx')].in_units('cm/s').value
    #vely = data[('Corrected_vely')].in_units('cm/s').value
    velx = data['velx'].in_units('cm/s').value
    vely = data['vely'].in_units('cm/s').value
    vels = np.array([velx, vely])
    del velx
    del vely
    pos_vec = np.array([-1*normal[1], normal[0]])
    if np.shape(data)[0] != 16:
        vels = vels.T
        c = (np.dot(vels,pos_vec))/(np.dot(pos_vec,pos_vec))
        pos_mag = np.sqrt(pos_vec[0]**2. + pos_vec[1]**2.)
        del pos_vec
        vels = c*pos_mag
        del c
        vels = yt.YTArray(vels, 'cm/s')
    else:
        #vels = data['Corrected_vely']
        vels = data['vely'].in_units('cm/s').value
    return vels

yt.add_field("Projected_Velocity", function=_Projected_Velocity, units=r"cm/s")

def _Projected_Magnetic_Field(field, data):
    """
    Calculates the projected magnetic field on the plane perpendicular to the normal vector L.
    """
    global normal
    magx = data[('magx')].value
    magy = data[('magy')].value
    mags = np.array([magx, magy])
    del magx
    del magy
    pos_vec = np.array([-1*normal[1], normal[0]])
    if np.shape(data)[0] != 16:
        mags = mags.T
        c = (np.dot(mags,pos_vec))/(np.dot(pos_vec,pos_vec))
        pos_mag = np.sqrt(pos_vec[0]**2. + pos_vec[1]**2.)
        del pos_vec
        mags = c*pos_mag
        del c
        mags = yt.YTArray(mags, 'gauss')
    else:
        mags = data['magy']
    return mags

yt.add_field("Projected_Magnetic_Field", function=_Projected_Magnetic_Field, units=r"gauss")

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

def _Particle_Potential(field, data):
    """
    Calculates the potential from the praticles.
    """
    Part_gpot = yt.YTArray(np.zeros(np.shape(data)), 'cm**2/s**2')
    if ('all', u'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        for part in range(len(dd['particle_mass'])):
            dx = data['x'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
            dy = data['y'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
            dz = data['z'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
            r = np.sqrt(dx**2. + dy**2. + dz**2.)
            gpot = -(yt.physical_constants.G*dd['particle_mass'][part].in_units('g'))/r
            Part_gpot = Part_gpot + gpot
        del dx
        del dy
        del dz
        del r
        del gpot
        del dd
    return Part_gpot

yt.add_field("Particle_Potential", function=_Particle_Potential, units=r"cm**2/s**2")

def _Total_Potential(field, data):
    """
    Gives the total potential inclusing contribution from the gas and the sink particles.
    """
    gas_pot = data['gpot']
    part_pot = data['Particle_Potential']
    G_pot_total = gas_pot + part_pot
    del gas_pot
    del part_pot
    return G_pot_total

yt.add_field("Total_Potential", function=_Total_Potential, units=r"cm**2/s**2")

def _Keplerian_Velocity(field, data):
    """
    Keplerian velocity calculated from the total potential energy (Sum of the potential from the sinks and the gas)
    """
    total_pot = data['Total_Potential']
    keplerian_field = np.sqrt(-1*total_pot)
    del total_pot
    return keplerian_field

yt.add_field("Keplerian_Velocity", function=_Keplerian_Velocity, units=r"cm/s")

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

def _Projected_Velocity_mw(field, data):
    """
    field to be able to created a mass weighted projection of Projected_Velocity
    """
    proj_vel = data['Projected_Velocity']
    cell_mass = data['cell_mass']
    vel_mw = proj_vel*cell_mass
    del proj_vel
    del cell_mass
    return vel_mw

yt.add_field("Projected_Velocity_mw", function=_Projected_Velocity_mw, units=r"cm*g/s")

def _velz_mw(field, data):
    """
    field to be able to created a mass weighted projection of velz
    """
    velz = data['velz']
    cell_mass = data['cell_mass']
    vel_mw = velz*cell_mass
    del velz
    del cell_mass
    return vel_mw

yt.add_field("velz_mw", function=_velz_mw, units=r"cm*g/s")

def _Projected_Magnetic_Field_mw(field, data):
    """
    field to be able to created a mass weighted projection of Projected_Magnetic_Field
    """
    proj_mag = data['Projected_Magnetic_Field']
    cell_mass = data['cell_mass']
    mag_mw = proj_mag*cell_mass
    del proj_mag
    del cell_mass
    return mag_mw

yt.add_field("Projected_Magnetic_Field_mw", function=_Projected_Magnetic_Field_mw, units=r"gauss*g")

def _magz_mw(field, data):
    """
    field to be able to created a mass weighted projection of magz
    """
    magz = data['magz']
    cell_mass = data['cell_mass']
    mag_mw = magz*cell_mass
    del magz
    del cell_mass
    return mag_mw

yt.add_field("magz_mw", function=_magz_mw, units=r"gauss*g")
