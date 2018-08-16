#!/usr/bin/env python
import yt
yt.enable_parallelism()
import numpy as np

center = 0
center_pos = [0.0, 0.0, 0.0]
center_vel = [0.0, 0.0, 0.0]
coordinates = 'spherical'
normal = [1.0, 0.0, 0.0]
part_pos = yt.YTArray([], 'cm')
part_mass = yt.YTArray([], 'g')
part_vel = yt.YTArray([], 'cm/s')

def set_center(x):
    """
    Sets the center used when calculating fields.
    
    Type: int
    Default: 0
    Options:0=center of mass, 1=particle 1, 2=particle 2.
    """
    global center
    global global_enc_mass
    center = x
    global_enc_mass = []
    return center

def set_part_pos(x):
    """
    Saves the particles positions.
    """
    global part_pos
    part_pos = x
    return part_pos

def set_part_vel(x):
    """
    Saves the particles velocities.
    """
    global part_vel
    part_vel = x
    return part_vel

def set_part_mass(x):
    """
    Saves the particles positions.
    """
    global part_mass
    part_mass = x
    return part_mass

def set_center_vel(x):
    """
    Sets the center velocity when calculating fields.
    """
    global center_vel
    center_vel = x
    return center_vel

def set_center_pos(x):
    """
    Sets the center position when calculating fields.
    """
    global center_pos
    center_pos = x
    return center_pos

def set_coordinate_system(x):
    """
    Sets the coordinates system used for some fields.
    
    Default: 'cylindrical'
    Options: 'cylindrical', 'spherical'
    """
    global coordinates
    coordinates = x
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

def get_part_pos():
    """
    returns the particles positions.
    """
    global part_pos
    return part_pos

def get_part_vel():
    """
    returns the particles velocities.
    """
    global part_vel
    return part_vel

def get_part_mass():
    """
    returns the particles positions.
    """
    global part_mass
    return part_mass

def get_center_vel():
    """
    Gets the center velocity when calculating fields.
    """
    global center_vel
    return center_vel

def get_center_pos():
    """
    Gets the center position when calculating fields.
    """
    global center_pos
    return center_pos

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

def _Center_Position(field, data):
    """
    Returns the center position for the current set center.
    """
    center = get_center()
    if center == 0:
        try:
            dd = data.ds.all_data()
            try:
                center_pos = dd.quantities.center_of_mass(use_particles=True)
            except:
                center_pos = dd.quantities.center_of_mass(use_particles=False)
        except:
            center_pos = yt.YTArray([0.0, 0.0, 0.0], 'cm')
    else:
        dd = data.ds.all_data()
        center_pos = [dd['particle_posx'][center-1].in_units('cm').value, dd['particle_posy'][center-1].in_units('cm').value, dd['particle_posz'][center-1].in_units('cm').value]
        center_pos = yt.YTArray(center_pos, 'cm')
    set_center_pos(center_pos)
    return center_pos

yt.add_field("Center_Position", function=_Center_Position, units=r"cm")

def _Center_Velocity(field, data):
    """
    Returns the center velocity for the current set center.
    """
    center = get_center()
    if center == 0:
        try:
            dd = data.ds.all_data()
            try:
                center_vel = dd.quantities.bulk_velocity(use_particles=True)
            except:
                center_vel = dd.quantities.bulk_velocity(use_particles=False)
        except:
            center_vel = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
    else:
        dd = data.ds.all_data()
        center_vel = yt.YTArray([dd['particle_velx'][center-1].in_units('cm/s').value, dd['particle_vely'][center-1].in_units('cm/s').value, dd['particle_velz'][center-1].in_units('cm/s').value], 'cm/s')
    set_center_vel(center_vel)
    return center_vel

yt.add_field("Center_Velocity", function=_Center_Velocity, units=r"cm/s")

def _All_Particle_Positions(field, data):
    """
    Saves all the particle positions
    """
    if ('all', u'particle_posx') in data.ds.field_list:
        dd = data.ds.all_data()
        pos = np.array([dd['particle_posx'].in_units('cm').value, dd['particle_posy'].in_units('cm').value, dd['particle_posz'].in_units('cm').value])
        pos = yt.YTArray(pos.T, 'cm')
    else:
        pos = yt.YTArray([], 'cm')
    set_part_pos(pos)
    return pos

yt.add_field("All_Particle_Positions", function=_All_Particle_Positions, units=r"cm")

def _All_Particle_Velocities(field, data):
    """
    Saves all the particle velocities
    """
    if ('all', u'particle_velx') in data.ds.field_list:
        dd = data.ds.all_data()
        vel = np.array([dd['particle_velx'].in_units('cm/s').value, dd['particle_vely'].in_units('cm/s').value, dd['particle_velz'].in_units('cm/s').value])
        vel = yt.YTArray(vel.T, 'cm/s')
    else:
        vel = yt.YTArray([], 'cm/s')
    set_part_vel(vel)
    return vel

yt.add_field("All_Particle_Velocities", function=_All_Particle_Velocities, units=r"cm/s")

def _All_Particle_Masses(field, data):
    """
    Saves all the particle masses
    """
    if ('all', u'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        mass = dd['particle_mass']
    else:
        mass = yt.YTArray([], 'g')
    set_part_mass(mass)
    return mass

yt.add_field("All_Particle_Masses", function=_All_Particle_Masses, units=r"g")

def _dx_from_Center(field, data):
    """
    Calculates the change in x position from the current set center.
    """
    center_pos = get_center_pos()
    dx = data['x'].in_units('cm')-center_pos[0]
    return dx

yt.add_field("dx_from_Center", function=_dx_from_Center, units=r"cm")

def _dy_from_Center(field, data):
    """
    Calculates the change in y position from the current set center.
    """
    center_pos = get_center_pos()
    dy = data['y'].in_units('cm')-center_pos[1]
    return dy

yt.add_field("dy_from_Center", function=_dy_from_Center, units=r"cm")

def _dz_from_Center(field, data):
    """
    Calculates the change in z position from the current set center.
    """
    center_pos = get_center_pos()
    dz = data['z'].in_units('cm')-center_pos[2]
    return dz

yt.add_field("dz_from_Center", function=_dz_from_Center, units=r"cm")

def _Distance_from_Center(field, data):
    """
    Calculates the distance from the current set center.
    """
    coordinates = get_coordinate_system()
    if 'cyl' in coordinates.lower():
        distance = np.sqrt((data['dx_from_Center'])**2. + (data['dy_from_Center'])**2.)
    else:
        distance = np.sqrt((data['dx_from_Center'])**2. + (data['dy_from_Center'])**2. + (data['dz_from_Center'])**2.)
    return distance

yt.add_field("Distance_from_Center", function=_Distance_from_Center, units=r"cm")

def _Corrected_velx(field, data):
    """
    Calculates the x-velocity corrected for the bulk velocity.
    """
    center_vel = get_center_vel()
    dvx = data['velx'].in_units('cm/s') - center_vel[0]
    return dvx

yt.add_field("Corrected_velx", function=_Corrected_velx, units=r"cm/s")

def _Corrected_vely(field, data):
    """
    Calculates the y-velocity corrected for the bulk velocity.
    """
    center_vel = get_center_vel()
    dvy = data['vely'].in_units('cm/s') - center_vel[1]
    return dvy

yt.add_field("Corrected_vely", function=_Corrected_vely, units=r"cm/s")

def _Corrected_velz(field, data):
    """
    Calculates the z-velocity corrected for the bulk velocity.
    """
    center_vel = get_center_vel()
    dvz = data['velz'].in_units('cm/s') - center_vel[2]
    return dvz

yt.add_field("Corrected_velz", function=_Corrected_velz, units=r"cm/s")

def _Corrected_vel_mag(field, data):
    """
    Calculates the velocity magnitude corrected for the bulk velocity in the current coordinate system.
    """
    coordinates = get_coordinate_system()
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
    normal = get_normal()
    velx = data[('Corrected_velx')].in_units('cm/s').value
    vely = data[('Corrected_vely')].in_units('cm/s').value
    vels = np.array([velx, vely])
    pos_vec = np.array([-1*normal[1], normal[0]])
    vels = vels.T
    c = (np.dot(vels,pos_vec))/(np.dot(pos_vec,pos_vec))
    pos_mag = np.sqrt(pos_vec[0]**2. + pos_vec[1]**2.)
    vels = c*pos_mag
    vels = yt.YTArray(vels, 'cm/s')
    return vels

yt.add_field("Projected_Velocity", function=_Projected_Velocity, units=r"cm/s")

def _Projected_Magnetic_Field(field, data):
    """
    Calculates the projected magnetic field on the plane perpendicular to the normal vector L.
    """
    normal = get_normal()
    magx = data[('magx')].value
    magy = data[('magy')].value
    mags = np.array([magx, magy])
    pos_vec = np.array([-1*normal[1], normal[0]])
    mags = mags.T
    c = (np.dot(mags,pos_vec))/(np.dot(pos_vec,pos_vec))
    pos_mag = np.sqrt(pos_vec[0]**2. + pos_vec[1]**2.)
    mags = c*pos_mag
    mags = yt.YTArray(mags, 'gauss')
    return mags

yt.add_field("Projected_Magnetic_Field", function=_Projected_Magnetic_Field, units=r"gauss")

def _Radial_Velocity(field, data):
    """
    Calculates the radial velocity from the current center, in the current coordinate system, corrected for the velocity of the center.
    """
    coordinates = get_coordinate_system()
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

def _Particle_Potential(field, data):
    """
    Calculates the potential from the praticles.
    """
    part_pos = get_part_pos()
    part_mass = get_part_mass()
    if len(part_pos) > 0 and len(part_mass) == 0:
        data['All_Particle_Masses']
        part_mass = get_part_mass()
    if len(part_pos) == 0 and len(part_mass) > 0:
        data['All_Particle_Positions']
        part_pos = get_part_pos()
    Part_gpot = yt.YTArray(np.zeros(np.shape(data['dens'])), 'cm**2/s**2')
    if len(part_mass) > 0:
        for part in range(len(part_mass)):
            dx = data['x'].in_units('cm') - part_pos[part][0].in_units('cm')
            dy = data['y'].in_units('cm') - part_pos[part][1].in_units('cm')
            dz = data['z'].in_units('cm') - part_pos[part][2].in_units('cm')
            r = np.sqrt(dx**2. + dy**2. + dz**2.)
            gpot = -(yt.physical_constants.G*part_mass[part])/r
            Part_gpot = Part_gpot + gpot
    return Part_gpot

yt.add_field("Particle_Potential", function=_Particle_Potential, units=r"cm**2/s**2", particle_type=True)

def _Total_Potential(field, data):
    """
    Gives the total potential inclusing contribution from the gas and the sink particles.
    """
    gas_pot = data['gpot']
    part_pot = data['Particle_Potential']
    G_pot_total = gas_pot + part_pot
    return G_pot_total

yt.add_field("Total_Potential", function=_Total_Potential, units=r"cm**2/s**2")

def _Keplerian_Velocity(field, data):
    """
    Keplerian velocity calculated from the total potential energy (Sum of the potential from the sinks and the gas)
    """
    total_pot = data['Total_Potential']
    keplerian_field = np.sqrt(-1*total_pot)
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
'''
def _Projected_Velocity_mw(field, data):
    """
    field to be able to created a mass weighted projection of Projected_Velocity
    """
    proj_vel = data['Projected_Velocity']
    cell_mass = data['cell_mass']
    vel_mw = proj_vel*cell_mass
    return vel_mw

yt.add_field("Projected_Velocity_mw", function=_Projected_Velocity_mw, units=r"cm*g/s")

def _velx_mw(field, data):
    """
    field to be able to created a mass weighted projection of velx
    """
    velx = data['velx']
    cell_mass = data['cell_mass']
    vel_mw = velx*cell_mass
    return vel_mw

yt.add_field("velx_mw", function=_velx_mw, units=r"cm*g/s")

def _vely_mw(field, data):
    """
    field to be able to created a mass weighted projection of vely
    """
    vely = data['vely']
    cell_mass = data['cell_mass']
    vel_mw = vely*cell_mass
    return vel_mw

yt.add_field("vely_mw", function=_vely_mw, units=r"cm*g/s")

def _velz_mw(field, data):
    """
    field to be able to created a mass weighted projection of velz
    """
    velz = data['velz']
    cell_mass = data['cell_mass']
    vel_mw = velz*cell_mass
    return vel_mw

yt.add_field("velz_mw", function=_velz_mw, units=r"cm*g/s")

def _Projected_Magnetic_Field_mw(field, data):
    """
    field to be able to created a mass weighted projection of Projected_Magnetic_Field
    """
    proj_mag = data['Projected_Magnetic_Field']
    cell_mass = data['cell_mass']
    mag_mw = proj_mag*cell_mass
    return mag_mw

yt.add_field("Projected_Magnetic_Field_mw", function=_Projected_Magnetic_Field_mw, units=r"gauss*g")

def _magx_mw(field, data):
    """
    field to be able to created a mass weighted projection of magx
    """
    magx = data['magx']
    cell_mass = data['cell_mass']
    mag_mw = magx*cell_mass
    return mag_mw

yt.add_field("magx_mw", function=_magx_mw, units=r"gauss*g")

def _magy_mw(field, data):
    """
    field to be able to created a mass weighted projection of magx
    """
    magy = data['magy']
    cell_mass = data['cell_mass']
    mag_mw = magy*cell_mass
    return mag_mw

yt.add_field("magy_mw", function=_magy_mw, units=r"gauss*g")

def _magz_mw(field, data):
    """
    field to be able to created a mass weighted projection of magz
    """
    magz = data['magz']
    cell_mass = data['cell_mass']
    mag_mw = magz*cell_mass
    return mag_mw

yt.add_field("magz_mw", function=_magz_mw, units=r"gauss*g")
'''
def _Gravitational_Force_on_particles_x(field, data):
    """
    Calculates the x component of the gravitational force on the sink particles
    """
    try:
        dd = data.ds.all_data()
        F_x = yt.YTArray(np.zeros(np.shape(dd['particle_mass'])), 'cm*g/s**2')
        cell_mass = dd['cell_mass'].in_units('g')
        for part in range(len(dd['particle_mass'])):
            dx = dd['x'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
            dy = dd['y'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
            dz = dd['z'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
            r = np.sqrt(dx**2. + dy**2. + dz**2.)
            F_x_arr = ((yt.physical_constants.G*dd['particle_mass'][part]*cell_mass)/r**3.)*dx
            F_x_tot = np.sum(F_x_arr)
            if len(dd['particle_mass']) > 1:
                dx = dd['particle_posx'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
                dy = dd['particle_posy'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
                dz = dd['particle_posz'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
                r = np.sqrt(dx**2. + dy**2. + dz**2.)
                inds = np.argwhere(r != 0.0)[0]
                F_part = ((yt.physical_constants.G*dd['particle_mass']*dd['particle_mass'][part])/r**3.)*dx
                F_x_tot = F_x_tot + np.sum(F_part[inds])
            F_x[part] = F_x_tot
    except:
        F_x = yt.YTArray([], 'cm*g/s**2')
    return F_x

yt.add_field("Gravitational_Force_on_particles_x", function=_Gravitational_Force_on_particles_x, units=r"cm*g/s**2", particle_type=True)

def _Gravitational_Force_on_particles_y(field, data):
    """
    Calculates the y component of the gravitational force on the sink particles
    """
    try:
        dd = data.ds.all_data()
        F_y = yt.YTArray(np.zeros(np.shape(dd['particle_mass'])), 'cm*g/s**2')
        cell_mass = dd['cell_mass'].in_units('g')
        for part in range(len(dd['particle_mass'])):
            dx = dd['x'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
            dy = dd['y'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
            dz = dd['z'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
            r = np.sqrt(dx**2. + dy**2. + dz**2.)
            F_y_arr = ((yt.physical_constants.G*dd['particle_mass'][part]*cell_mass)/r**3.)*dy
            F_y_tot = np.sum(F_y_arr)
            if len(dd['particle_mass']) > 1:
                dx = dd['particle_posx'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
                dy = dd['particle_posy'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
                dz = dd['particle_posz'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
                r = np.sqrt(dx**2. + dy**2. + dz**2.)
                inds = np.argwhere(r != 0.0)[0]
                F_part = ((yt.physical_constants.G*dd['particle_mass']*dd['particle_mass'][part])/r**3.)*dy
                F_y_tot = F_y_tot + np.sum(F_part[inds])
            F_y[part] = F_y_tot
    except:
        F_y = yt.YTArray([], 'cm*g/s**2')
    return F_y

yt.add_field("Gravitational_Force_on_particles_y", function=_Gravitational_Force_on_particles_y, units=r"cm*g/s**2", particle_type=True)

def _Gravitational_Force_on_particles_z(field, data):
    """
    Calculates the z component of the gravitational force on the sink particles
    """
    try:
        dd = data.ds.all_data()
        F_z = yt.YTArray(np.zeros(np.shape(dd['particle_mass'])), 'cm*g/s**2')
        cell_mass = dd['cell_mass'].in_units('g')
        for part in range(len(dd['particle_mass'])):
            dx = dd['x'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
            dy = dd['y'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
            dz = dd['z'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
            r = np.sqrt(dx**2. + dy**2. + dz**2.)
            F_z_arr = ((yt.physical_constants.G*dd['particle_mass'][part]*cell_mass)/r**3.)*dz
            F_z_tot = np.sum(F_z_arr)
            if len(dd['particle_mass']) > 1:
                dx = dd['particle_posx'].in_units('cm') - dd['particle_posx'][part].in_units('cm')
                dy = dd['particle_posy'].in_units('cm') - dd['particle_posy'][part].in_units('cm')
                dz = dd['particle_posz'].in_units('cm') - dd['particle_posz'][part].in_units('cm')
                r = np.sqrt(dx**2. + dy**2. + dz**2.)
                inds = np.argwhere(r != 0.0)[0]
                F_part = ((yt.physical_constants.G*dd['particle_mass']*dd['particle_mass'][part])/r**3.)*dz
                F_z_tot = F_z_tot + np.sum(F_part[inds])
            F_z[part] = F_z_tot
    except:
        F_z = yt.YTArray([], 'cm*g/s**2')
    return F_z

yt.add_field("Gravitational_Force_on_particles_z", function=_Gravitational_Force_on_particles_z, units=r"cm*g/s**2", particle_type=True)

def _Gravitational_Force_on_particles_mag(field, data):
    """
    Calculates the magnitude of the gravitational force on the sink particles
    """
    F_mag = np.sqrt(data['Gravitational_Force_on_particles_x']**2. + data['Gravitational_Force_on_particles_y']**2. + data['Gravitational_Force_on_particles_z']**2.)
    return F_mag

yt.add_field("Gravitational_Force_on_particles_mag", function=_Gravitational_Force_on_particles_mag, units=r"cm*g/s**2", particle_type=True)

def _Gravitational_Force_on_particles_Rad(field, data):
    """
    The component of the gravitational force on the particles in the radial direction from the current center.
    """
    try:
        dd = data.ds.all_data()
        center_pos = dd['Center_Position']
        F_rad = yt.YTArray(np.zeros(np.shape(dd['particle_mass'])), 'cm*g/s**2')
        for part in range(len(dd['particle_mass'])):
            dx = center_pos[0].in_units('cm') - dd['particle_posx'][part].in_units('cm')
            dy = center_pos[1].in_units('cm') - dd['particle_posy'][part].in_units('cm')
            dz = center_pos[2].in_units('cm') - dd['particle_posz'][part].in_units('cm')
            r = np.array([dx.value, dy.value, dz.value])
            F_x = dd['Gravitational_Force_on_particles_x'][part].value
            F_y = dd['Gravitational_Force_on_particles_y'][part].value
            F_z = dd['Gravitational_Force_on_particles_z'][part].value
            F = np.array([F_x, F_y, F_z])
            F_proj = (np.dot(F, r)/np.dot(r,r))*r
            F_mag = np.sqrt(F_proj[0]**2. + F_proj[1]**2. + F_proj[2]**2.)
            F_rad[part] = F_mag
    except:
        F_rad = yt.YTArray([], 'cm*g/s**2')
    return F_rad

yt.add_field("Gravitational_Force_on_particles_Rad", function=_Gravitational_Force_on_particles_Rad, units=r"cm*g/s**2", particle_type=True)

def _Gravitational_Force_on_particles_Tan(field, data):
    """
    The component of the gravitational force on the particles in the radial direction from the current center.
    """
    F_mag = data['Gravitational_Force_on_particles_mag']
    F_rad = data['Gravitational_Force_on_particles_Rad']
    F_tan = np.sqrt(F_mag**2. - F_rad**2.)
    return F_tan

yt.add_field("Gravitational_Force_on_particles_Tan", function=_Gravitational_Force_on_particles_Tan, units=r"cm*g/s**2", particle_type=True)

def _Angular_Momentum_x(field, data):
    """
    Calculates the angular momentum in the x_direction about current set center.
    """
    L_x = data['cell_mass']*(data['Corrected_velz']*data['dy_from_Center']- data['Corrected_vely']*data['dz_from_Center'])
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

def _Magnetic_To_Gas_Pressure(field, data):
    """
    Calculates the specific angular momentum about current set center.
    """
    mag_pres = data['magnetic_pressure']
    gas_pres = data['pressure']
    ratio = mag_pres/gas_pres
    return ratio

yt.add_field("Magnetic_To_Gas_Pressure", function=_Magnetic_To_Gas_Pressure, units=r"")

'''
def _Particle_Angular_Momentum_x(field, data):
    """
    Calculates the angular momentum of the particles in the x_direction about current set center.
    """
    global part_mass
    global part_vel
    global center_vel
    global part_pos
    global center_pos
    L_x = part_mass*((part_vel.T[2] - center_vel[2])*(part_pos.T[1] - center_pos[1]) - (part_vel.T[1] - center_vel[1])*(part_pos.T[2] - center_pos[2]))
    L_x = yt.YTArray(L_x.value, 'g*cm**2/s')
    return L_x

yt.add_field("Particle_Angular_Momentum_x", function=_Particle_Angular_Momentum_x, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_y(field, data):
    """
    Calculates the angular momentum of the particles in the y_direction about current set center.
    """
    global part_mass
    global part_vel
    global center_vel
    global part_pos
    global center_pos
    L_y = part_mass*((part_vel.T[2] - center_vel[2])*(part_pos.T[0] - center_pos[0]) - (part_vel.T[0] - center_vel[0])*(part_pos.T[2] - center_pos[2]))
    L_y = yt.YTArray(L_y.value, 'g*cm**2/s')
    return L_y

yt.add_field("Particle_Angular_Momentum_y", function=_Particle_Angular_Momentum_y, units=r"g*cm**2/s")

def _Particle_Angular_Momentum_z(field, data):
    """
    Calculates the angular momentum of the particles in the z_direction about current set center.
    """
    global part_mass
    global part_vel
    global center_vel
    global part_pos
    global center_pos
    L_z = part_mass*((part_vel.T[1] - center_vel[1])*(part_pos.T[0] - center_pos[0]) - (part_vel.T[0] - center_vel[0])*(part_pos.T[1] - center_pos[1]))
    L_z = yt.YTArray(L_z.value, 'g*cm**2/s')
    return L_z

yt.add_field("Particle_Angular_Momentum_z", function=_Particle_Angular_Momentum_z, units=r"g*cm**2/s")

def _Particle_Angular_Momentum(field, data):
    """
    Calculates the angular momentum about current set center for the particles
    """
    L = np.sqrt(data['Particle_Angular_Momentum_x']**2. + data['Particle_Angular_Momentum_y']**2. + data['Particle_Angular_Momentum_z']**2.)
    return L

yt.add_field("Particle_Angular_Momentum", function=_Particle_Angular_Momentum, units=r"g*cm**2/s")

def _Particle_Specific_Angular_Momentum(field, data):
    """
    Calculates the specific angular momentum for the particles about current set center.
    """
    global part_mass
    l = data['Particle_Angular_Momentum']/part_mass
    return l

yt.add_field("Particle_Specific_Angular_Momentum", function=_Particle_Specific_Angular_Momentum, units=r"cm**2/s")
'''
