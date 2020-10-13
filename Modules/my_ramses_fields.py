#!/usr/bin/env python
#Created by Rajika Kuruwita, 2019
import yt
yt.enable_parallelism()
import numpy as np
import csv
import glob
from pyramses import rsink

center = 0
center_pos_gas = yt.YTArray([0.0, 0.0, 0.0], 'cm')
center_pos_part = yt.YTArray([0.0, 0.0, 0.0], 'cm')
center_pos = yt.YTArray([0.0, 0.0, 0.0], 'cm')
center_vel_gas = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
center_vel_part = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
center_vel = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
coordinates = 'spherical'
normal = [1.0, 0.0, 0.0]
part_pos = yt.YTArray([], 'cm')
part_mass = yt.YTArray([], 'g')
part_vel = yt.YTArray([], 'cm/s')
com_pos_use_gas = True
com_pos_use_part = True
com_vel_use_gas = True
com_vel_use_part = True
centred_sink_id = 0
active_radius = yt.YTArray(np.nan, 'au')

def set_centred_sink_id(x):
    """
    sets centred sink
    """
    global centred_sink_id
    centred_sink_id = x
    return centred_sink_id
    
def get_centred_sink_id():
    """
    returned centred sink id
    """
    global centred_sink_id
    return centred_sink_id

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

def set_com_pos_use_gas(x):
    """
    Sets whether to use gas when calculate center position
    
    Default: True
    """
    global com_pos_use_gas
    com_pos_use_gas = x
    return com_pos_use_gas
    
def set_com_vel_use_gas(x):
    """
    Sets whether to use gas when calculate center velocity

    Default: True
    """
    global com_vel_use_gas
    com_vel_use_gas = x
    return com_vel_use_gas
    
def set_com_pos_use_part(x):
    """
    Sets whether to use particles when calculate center position

    Default: True
    """
    global com_pos_use_part
    com_pos_use_part = x
    return com_pos_use_part
    
def set_com_vel_use_part(x):
    """
    Sets whether to use particles when calculate center velocity

    Default: True
    """
    global com_vel_use_part
    com_vel_use_part = x
    return com_vel_use_part
    
def set_active_radius(x):
    """
    Sets the active radius to consider when selecting which particles to use

    Default: 10000 AU
    """
    global active_radius
    active_radius = yt.YTArray(x, 'au')
    return active_radius

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
    
def get_com_pos_use_gas():
    """
    returns whether to use gas when calculate center position
    """
    global com_pos_use_gas
    return com_pos_use_gas

def get_com_vel_use_gas():
    """
    returns whether to use gas when calculate center velocity
    """
    global com_vel_use_gas
    return com_vel_use_gas

def get_com_pos_use_part():
    """
    returns whether to use particles when calculate center position
    """
    global com_pos_use_part
    return com_pos_use_part
    
def get_com_vel_use_part():
    """
    returns whether to use particles when calculate center velocity
    """
    global com_vel_use_part
    return com_vel_use_part
    
def get_active_radius():
    """
    returns the active radius to consider when selecting which particles to use
    """
    global active_radius
    return active_radius

#===========================OVERWRITING DENSITY FIELD BECAUSE DENSITY UNIT DOESN'T GET OVERWRITTEN======================================

def _Density(field,data):
    """
    Overwrites density field
    """
    density_unit = (data.ds.mass_unit/data.ds.length_unit**3).in_cgs()
    density = data[('ramses', 'Density')].value*density_unit
    density_arr = yt.YTArray(density, 'g/cm**3')
    del density
    del density_unit
    return density_arr

yt.add_field("Density", function=_Density, units=r"g/cm**3")
    
#===========================CALCULATING RAMSES CENTERED MAGNETIC FIELDS======================================

def _magx(field,data):
    """
    Calculated the centred x-component of the megnatic field
    """
    try:
        magx = yt.YTArray((data['x-Bfield-left'].value + data['x-Bfield-right'].value)/2., "G")
    except:
        magx = yt.YTArray((data['B_x_left'].value + data['B_x_right'].value)/2., "G")
    return magx
    
yt.add_field("magx", function=_magx, units=r"gauss")

def _magy(field,data):
    """
    Calculated the centred y-component of the megnatic field
    """
    try:
        magy = yt.YTArray((data['y-Bfield-left'].value + data['y-Bfield-right'].value)/2., "G")
    except:
        magy = yt.YTArray((data['B_y_left'].value + data['B_y_right'].value)/2., "G")
    return magy
    
yt.add_field("magy", function=_magy, units=r"gauss")

def _magz(field,data):
    """
    Calculated the centred z-component of the megnatic field
    """
    try:
        magz = yt.YTArray((data['z-Bfield-left'].value + data['z-Bfield-right'].value)/2., "G")
    except:
        magz = yt.YTArray((data['B_z_left'].value + data['B_z_right'].value)/2., "G")
    return magz
    
yt.add_field("magz", function=_magz, units=r"gauss")

#===========================READING RAMSES PARTICLE DATA======================================

def _sink_particle_tag(field, data):
    """
    Retrieve particle tags from .snktxt file
    """
    particle_tag = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_tag = yt.YTArray(np.array(particle_tag), "")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_tag = np.arange(float(len(loaded_sink_data['x'])))
        '''
        #particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        #csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        #time_val = None
        #with open(particle_file[0], 'r') as f:
        #    reader = csv.reader(f, dialect='dat')
        #    for row in reader:
        #        if time_val == None:
        #            time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
        #            continue
        #        try:
        #            dummy = float(row[0])
        #            particle_tag.append(float(row[0]))
        #        except:
        #            continue
        '''
    particle_tag = yt.YTArray(np.array(particle_tag), "")
    return particle_tag

yt.add_field("sink_particle_tag", function=_sink_particle_tag, units=r"")

def _sink_particle_posx(field, data):
    """
    Retrieve particle x position from .snktxt file
    """
    particle_posx = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_posx = yt.YTArray(np.array(particle_posx), "pc")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_posx = loaded_sink_data['x']*data.ds.length_unit.in_units("pc").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_posx.append(float(row[1])*data.ds.length_unit.value)
                except:
                    continue
        '''
        particle_posx = yt.YTArray(np.array(particle_posx), "pc")
    return particle_posx

yt.add_field("sink_particle_posx", function=_sink_particle_posx, units=r"pc")

def _sink_particle_posy(field, data):
    """
    Retrieve particle y position from .snktxt file
    """
    particle_posy = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_posy = yt.YTArray(np.array(particle_posy), "pc")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_posy = loaded_sink_data['y']*data.ds.length_unit.in_units("pc").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_posy.append(float(row[2])*data.ds.length_unit.value)
                except:
                    continue
        '''
        particle_posy = yt.YTArray(np.array(particle_posy), "pc")
    return particle_posy

yt.add_field("sink_particle_posy", function=_sink_particle_posy, units=r"pc")

def _sink_particle_posz(field, data):
    """
    Retrieve particle z position from .snktxt file
    """
    particle_posz = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_posz = yt.YTArray(np.array(particle_posz), "pc")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_posz = loaded_sink_data['z']*data.ds.length_unit.in_units("pc").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_posz.append(float(row[3])*data.ds.length_unit.value)
                except:
                    continue
        '''
        particle_posz = yt.YTArray(np.array(particle_posz), "pc")
    return particle_posz

yt.add_field("sink_particle_posz", function=_sink_particle_posz, units=r"pc")

def _sink_particle_velx(field, data):
    """
    Retrieve particle x velocity from .snktxt file
    """
    particle_velx = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_velx = yt.YTArray(np.array(particle_velx), "km/s")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_velx = loaded_sink_data['ux']*data.ds.velocity_unit.in_units("km/s").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_velx.append(float(row[4])*data.ds.velocity_unit.value)
                except:
                    continue
        '''
        particle_velx = yt.YTArray(np.array(particle_velx), "km/s")
    return particle_velx

yt.add_field("sink_particle_velx", function=_sink_particle_velx, units=r"km/s")

def _sink_particle_vely(field, data):
    """
    Retrieve particle y velocity from .snktxt file
    """
    particle_vely = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_vely = yt.YTArray(np.array(particle_vely), "km/s")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_vely = loaded_sink_data['uy']*data.ds.velocity_unit.in_units("km/s").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_vely.append(float(row[5])*data.ds.velocity_unit.value)
                except:
                    continue
        '''
        particle_vely = yt.YTArray(np.array(particle_vely), "km/s")
    return particle_vely

yt.add_field("sink_particle_vely", function=_sink_particle_vely, units=r"km/s")

def _sink_particle_velz(field, data):
    """
    Retrieve particle z velocity from .snktxt file
    """
    particle_velz = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_velz = yt.YTArray(np.array(particle_velz), "km/s")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_velz = loaded_sink_data['uz']*data.ds.velocity_unit.in_units("km/s").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_velz.append(float(row[6])*data.ds.velocity_unit.value)
                except:
                    continue
        '''
        particle_velz = yt.YTArray(np.array(particle_velz), "km/s")
    return particle_velz

yt.add_field("sink_particle_velz", function=_sink_particle_velz, units=r"km/s")

def _sink_particle_speed(field, data):
    """
    Retrieve particle speed from sink particle file
    """
    particle_speed = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_speed = yt.YTArray(np.array(particle_speed), "km/s")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_speed = loaded_sink_data['u']*data.ds.velocity_unit.in_units("km/s").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_velz.append(float(row[6])*data.ds.velocity_unit.value)
                except:
                    continue
        '''
        particle_speed = yt.YTArray(np.array(particle_speed), "km/s")
    return particle_speed

yt.add_field("sink_particle_speed", function=_sink_particle_speed, units=r"km/s")

def _sink_particle_mass(field, data):
    """
    Retrieve particle mass from .snktxt file
    """
    particle_mass = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_mass = yt.YTArray(np.array(particle_mass), "Msun")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_mass = loaded_sink_data['m']*data.ds.mass_unit.in_units("Msun").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_mass.append(float(row[7])*data.ds.mass_unit.value)
                except:
                    continue
        '''
        particle_mass = yt.YTArray(np.array(particle_mass), "Msun")
    return particle_mass

yt.add_field("sink_particle_mass", function=_sink_particle_mass, units=r"Msun")

def _sink_particle_accretion_rate(field, data):
    """
    Retrieve particle accretion rate from sink file
    """
    particle_mdot = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_mdot = yt.YTArray(np.array(particle_mdot), "Msun/yr")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        numerator = loaded_sink_data['dm']*data.ds.mass_unit.in_units("Msun").value
        denominator = (loaded_sink_data['snapshot_time'] - loaded_sink_data['tflush'])*data.ds.time_unit.in_units("yr").value
        particle_mdot = numerator/denominator
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_mass.append(float(row[7])*data.ds.mass_unit.value)
                except:
                    continue
        '''
        particle_mdot = yt.YTArray(np.array(particle_mdot), "Msun/yr")
    return particle_mdot

yt.add_field("sink_particle_accretion_rate", function=_sink_particle_accretion_rate, units=r"Msun/yr")

def _sink_particle_form_time(field, data):
    """
    Retrieve particle formation time from .snktxt file
    """
    particle_form_time = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_form_time = yt.YTArray(np.array(particle_form_time), "yr")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_form_time = loaded_sink_data['tcreate']*data.ds.time_unit.in_units("yr").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_form_time.append(float(row[10])*data.ds.time_unit.value)
                except:
                    continue
        '''
        particle_form_time = yt.YTArray(np.array(particle_form_time), "yr")
    return particle_form_time

yt.add_field("sink_particle_form_time", function=_sink_particle_form_time, units=r"yr")

def _sink_particle_age(field, data):
    """
    Retrieve particle age from .snktxt file
    """
    particle_age = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_age = yt.YTArray(np.array(particle_age), "yr")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        particle_age = (loaded_sink_data['snapshot_time']-loaded_sink_data['tcreate'])*data.ds.time_unit.in_units("yr").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    age = time_val.in_units('s').value - float(row[10])*data.ds.time_unit.value
                    particle_age.append(age)
                except:
                    continue
        '''
        particle_age = yt.YTArray(np.array(particle_age), "yr")
    return particle_age

yt.add_field("sink_particle_age", function=_sink_particle_age, units=r"yr")

#===========================REPLACING YT PARTICLE FIELDS WITH RAMSES ONES======================================
'''
def _particle_identity(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_tag']
    
yt.add_field("particle_tag", function=_particle_tag, units=r"")

def _particle_position_x(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_posx'].in_units('cm')
    
yt.add_field("particle_position_x", function=_particle_position_x, units=r"cm")

def _particle_position_y(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_posy'].in_units('cm')
    
yt.add_field("particle_position_y", function=_particle_position_y, units=r"cm")

def _particle_position_z(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_posz'].in_units('cm')
    
yt.add_field("particle_position_z", function=_particle_position_z, units=r"cm")

def _particle_velocity_x(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_velx'].in_units('cm/s')
    
yt.add_field("particle_velocity_x", function=_particle_velocity_x, units=r"cm/s")

def _particle_velocity_y(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_vely'].in_units('cm/s')
    
yt.add_field("particle_velocity_y", function=_particle_velocity_y, units=r"cm/s")

def _particle_velocity_z(field, data):
    """
    replacing yt particle fields
    """
    return data['particle_velz'].in_units('cm/s')
    
yt.add_field("particle_velocity_z", function=_particle_velocity_z, units=r"cm/s")
'''
def _Center_Position_Gas(field, data):
    """
    Calculates the CoM of gas
    """
    try:
        dd = data.ds.all_data()
        TM = np.sum(dd['cell_mass'].in_units('g'))
        x_top = np.sum(dd['cell_mass'].in_units('g')*dd[('index','x')].in_units('cm'))
        y_top = np.sum(dd['cell_mass'].in_units('g')*dd[('index','y')].in_units('cm'))
        z_top = np.sum(dd['cell_mass'].in_units('g')*dd[('index','z')].in_units('cm'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm')
    center_pos_gas = com
    return com
    
yt.add_field("Center_Position_Gas", function=_Center_Position_Gas, units=r"cm")

def _Center_Position_Particle(field, data):
    """
    Calculates the CoM of gas
    """
    global centred_sink_id
    global active_radius
    try:
        dd = data.ds.all_data()
        if np.isnan(active_radius):
            usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            usable_tags = np.array(usable_tags)
        else:
            centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
            dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
            dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
            dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
            dist = np.sqrt(dx**2+dy**2+dz**2)
            usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
        TM = np.sum(dd['sink_particle_mass'][np.array(usable_tags)].in_units('g'))
        x_top = np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_posx'][usable_tags].in_units('cm'))
        y_top = np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_posy'][usable_tags].in_units('cm'))
        z_top = np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_posz'][usable_tags].in_units('cm'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm')
    center_pos_part = com
    return com
    
yt.add_field("Center_Position_Particle", function=_Center_Position_Particle, units=r"cm")

def _Center_Position(field, data):
    """
    Returns the center position which is derived from the sink particles with tags equal to or greater the set centred particle
    """
    global com_pos_use_gas
    global com_pos_use_part
    global center
    global active_radius
    try:
        dd = data.ds.all_data()
        if center == 0:
            TM = yt.YTArray(0.0, 'g')
            x_top = yt.YTArray(0.0, 'cm*g')
            y_top = yt.YTArray(0.0, 'cm*g')
            z_top = yt.YTArray(0.0, 'cm*g')
            if com_pos_use_part == True:
                try:
                    global centred_sink_id
                    if np.isnan(active_radius):
                        usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
                        usable_tags = np.array(usable_tags)
                    else:
                        centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
                        dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
                        dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
                        dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
                        dist = np.sqrt(dx**2+dy**2+dz**2)
                        usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
                    M_part = np.sum(dd['sink_particle_mass'][np.array(usable_tags)].in_units('g'))
                    TM = TM + M_part
                    x_top = x_top + np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_posx'][usable_tags].in_units('cm'))
                    y_top = y_top + np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_posy'][usable_tags].in_units('cm'))
                    z_top = z_top + np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_posz'][usable_tags].in_units('cm'))
                except:
                    TM = TM
                    x_top = x_top
                    y_top = y_top
                    z_top = z_top
            if com_pos_use_gas == True:
                M_gas = np.sum(dd['cell_mass'].in_units('g'))
                TM = TM + M_gas
                x_top = x_top + np.sum(dd['cell_mass'].in_units('g')*dd[('index','x')].in_units('cm'))
                y_top = y_top + np.sum(dd['cell_mass'].in_units('g')*dd[('index','y')].in_units('cm'))
                z_top = z_top + np.sum(dd['cell_mass'].in_units('g')*dd[('index','z')].in_units('cm'))
            com = [(x_top/TM), (y_top/TM), (z_top/TM)]
            center_pos = yt.YTArray(com, 'cm')
        else:
            particle_tag = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            center_tag = int(particle_tag[center-1])
            center_pos = yt.YTArray([dd['sink_particle_posx'][center_tag].in_units('cm').value, dd['sink_particle_posy'][center_tag].in_units('cm').value, dd['sink_particle_posz'][center_tag].in_units('cm').value], 'cm')
    except:
        center_pos = data.ds.domain_center
    set_center_pos(center_pos)
    return center_pos

yt.add_field("Center_Position", function=_Center_Position, units=r"cm")

def _Center_Velocity_Gas(field, data):
    """
    Calculates the mass weighted bulk velocity of the gas
    """
    try:
        dd = data.ds.all_data()
        TM = np.sum(dd['cell_mass'].in_units('g'))
        x_top = np.sum(dd['cell_mass'].in_units('g')*dd['x-velocity'].in_units('cm/s'))
        y_top = np.sum(dd['cell_mass'].in_units('g')*dd['y-velocity'].in_units('cm/s'))
        z_top = np.sum(dd['cell_mass'].in_units('g')*dd['z-velocity'].in_units('cm/s'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
    center_vel_gas = com
    return com
    
yt.add_field("Center_Velocity_Gas", function=_Center_Velocity_Gas, units=r"cm/s")

def _Center_Velocity_Particle(field, data):
    """
    Calculates the mass weighted bulk velocity of the particles
    """
    global centred_sink_id
    global active_radius
    try:
        dd = data.ds.all_data()
        if np.isnan(active_radius):
            usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            usable_tags = np.array(usable_tags)
        else:
            centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
            dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
            dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
            dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
            dist = np.sqrt(dx**2+dy**2+dz**2)
            usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
        TM = np.sum(dd['sink_particle_mass'][np.array(usable_tags)].in_units('g'))
        x_top = np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_velx'][usable_tags].in_units('cm/s'))
        y_top = np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_vely'][usable_tags].in_units('cm/s'))
        z_top = np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_velz'][usable_tags].in_units('cm/s'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
    center_vel_part = com
    return com

yt.add_field("Center_Velocity_Particle", function=_Center_Velocity_Particle, units=r"cm/s")

def _Center_Velocity(field, data):
    """
    Returns the center velocity for the current set center.
    """
    global com_vel_use_gas
    global com_vel_use_part
    global centred_sink_id
    global center
    global active_radius
    try:
        if center == 0:
            TM = yt.YTArray(0.0, 'g')
            x_top = yt.YTArray(0.0, 'cm*g/s')
            y_top = yt.YTArray(0.0, 'cm*g/s')
            z_top = yt.YTArray(0.0, 'cm*g/s')
            if com_vel_use_part == True:
                try:
                    global centred_sink_id
                    if np.isnan(active_radius):
                        usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
                        usable_tags = np.array(usable_tags)
                    else:
                        centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
                        dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
                        dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
                        dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
                        dist = np.sqrt(dx**2+dy**2+dz**2)
                        usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
                    M_part = np.sum(dd['sink_particle_mass'][np.array(usable_tags)].in_units('g'))
                    TM = TM + M_part
                    x_top = x_top + np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_velx'][usable_tags].in_units('cm/s'))
                    y_top = y_top + np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_vely'][usable_tags].in_units('cm/s'))
                    z_top = z_top + np.sum(dd['sink_particle_mass'][usable_tags].in_units('g')*dd['sink_particle_velz'][usable_tags].in_units('cm/s'))
                except:
                    TM = TM
                    x_top = x_top
                    y_top = y_top
                    z_top = z_top
            if com_vel_use_gas == True:
                M_gas = np.sum(data['cell_mass'].in_units('g'))
                TM = TM + M_gas
                x_top = x_top + np.sum(data['cell_mass'].in_units('g')*data['x-velocity'].in_units('cm/s'))
                y_top = y_top + np.sum(data['cell_mass'].in_units('g')*data['y-velocity'].in_units('cm/s'))
                z_top = z_top + np.sum(data['cell_mass'].in_units('g')*data['z-velocity'].in_units('cm/s'))
            com_vel = [(x_top/TM), (y_top/TM), (z_top/TM)]
            center_vel = yt.YTArray(com_vel, 'cm')
        else:
            dd = data.ds.all_data()
            particle_tag = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            center_tag = int(particle_tag[center-1])
            center_vel = yt.YTArray([dd['sink_particle_velx'][center_tag].in_units('cm/s').value, dd['sink_particle_vely'][center_tag].in_units('cm/s').value, dd['sink_particle_velz'][center_tag].in_units('cm/s').value], 'cm/s')
    except:
        center_vel = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
    set_center_vel(center_vel)
    return center_vel

yt.add_field("Center_Velocity", function=_Center_Velocity, units=r"cm/s")

def _All_Particle_Positions(field, data):
    """
    Saves all the particle positions
    """
    if ('all', 'sink_particle_posx') in data.ds.field_list:
        dd = data.ds.all_data()
        if len(dd['sink_particle_posx'].in_units('cm').value) > 1:
            pos = np.array([dd['sink_particle_posx'].in_units('cm').value, dd['sink_particle_posy'].in_units('cm').value, dd['sink_particle_posz'].in_units('cm').value])
            pos = yt.YTArray(pos.T, 'cm')
        else:
            pos = yt.YTArray([[dd['sink_particle_posx'][0].in_units('cm').value, dd['sink_particle_posy'][0].in_units('cm').value, dd['sink_particle_posz'][0].in_units('cm').value]], 'cm')
    else:
        pos = yt.YTArray([], 'cm')
    set_part_pos(pos)
    return pos

yt.add_field("All_Particle_Positions", function=_All_Particle_Positions, units=r"cm")

def _All_Particle_Velocities(field, data):
    """
    Saves all the particle velocities
    """
    if ('all', 'sink_particle_velx') in data.ds.field_list:
        dd = data.ds.all_data()
        if len(dd['sink_particle_velx'].in_units('cm/s').value) > 1:
            vel = np.array([dd['sink_particle_velx'].in_units('cm/s').value, dd['sink_particle_vely'].in_units('cm/s').value, dd['sink_particle_velz'].in_units('cm/s').value])
            vel = yt.YTArray(vel.T, 'cm/s')
        else:
            vel = yt.YTArray([[dd['sink_particle_velx'][0].in_units('cm/s').value, dd['sink_particle_vely'][0].in_units('cm/s').value, dd['sink_particle_velz'][0].in_units('cm/s').value]], 'cm/s')
    else:
        vel = yt.YTArray([], 'cm/s')
    set_part_vel(vel)
    return vel

yt.add_field("All_Particle_Velocities", function=_All_Particle_Velocities, units=r"cm/s")

def _All_Particle_Masses(field, data):
    """
    Saves all the particle masses
    """
    if ('all', 'sink_particle_mass') in data.ds.field_list:
        if np.shape(data['x']) == (16, 16, 16):
            mass = data['sink_particle_mass']
        else:
            dd = data.ds.all_data()
            mass = dd['sink_particle_mass']
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
    dx = data['x'].in_units('cm')-center_pos[0].in_units('cm')
    return dx

yt.add_field("dx_from_Center", function=_dx_from_Center, units=r"cm")

def _dy_from_Center(field, data):
    """
    Calculates the change in y position from the current set center.
    """
    center_pos = get_center_pos()
    dy = data['y'].in_units('cm')-center_pos[1].in_units('cm')
    return dy

yt.add_field("dy_from_Center", function=_dy_from_Center, units=r"cm")

def _dz_from_Center(field, data):
    """
    Calculates the change in z position from the current set center.
    """
    center_pos = get_center_pos()
    dz = data['z'].in_units('cm')-center_pos[2].in_units('cm')
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
    dvx = data['x-velocity'].in_units('cm/s') - center_vel[0].in_units('cm/s')
    del center_vel
    return dvx

yt.add_field("Corrected_velx", function=_Corrected_velx, units=r"cm/s")

def _Corrected_vely(field, data):
    """
    Calculates the y-velocity corrected for the bulk velocity.
    """
    center_vel = get_center_vel()
    dvy = data['y-velocity'].in_units('cm/s') - center_vel[1].in_units('cm/s')
    del center_vel
    return dvy

yt.add_field("Corrected_vely", function=_Corrected_vely, units=r"cm/s")

def _Corrected_velz(field, data):
    """
    Calculates the z-velocity corrected for the bulk velocity.
    """
    center_vel = get_center_vel()
    dvz = data['z-velocity'].in_units('cm/s') - center_vel[2].in_units('cm/s')
    del center_vel
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
    
def _magnetic_field_magnitude(field, data):
    """
    Calculates the magnitude of the magnetic field
    """
    B = np.sqrt(data['magx']**2. + data['magy']**2. + data['magz']**2.)
    return B
    
yt.add_field("magnetic_field_magnitude", function=_magnetic_field_magnitude, units=r"G")

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

def _Orbital_Angular_Momentum_x(field, data):
    """
    Calculates the orbital angular momentume of all the necessaru sink particles
    """
    global centred_sink_id
    global center_vel_part
    global center_pos_part
    global active_radius
    try:
        dd = data.ds.all_data()
        if np.isnan(active_radius):
            usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            usable_tags = np.array(usable_tags)
        else:
            centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
            dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
            dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
            dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
            dist = np.sqrt(dx**2+dy**2+dz**2)
            usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
        M = dd['sink_particle_mass'][np.array(usable_tags)].in_units('g')
        V_z = dd['sink_particle_velz'][usable_tags].in_units('cm/s') - center_vel_part[2].in_units('cm/s')
        V_y = dd['sink_particle_vely'][usable_tags].in_units('cm/s') - center_vel_part[1].in_units('cm/s')
        P_y = dd['sink_particle_posy'][usable_tags].in_units('cm') - center_pos_part[1].in_units('cm')
        P_z = dd['sink_particle_posz'][usable_tags].in_units('cm') - center_pos_part[2].in_units('cm')
        L_x = M *((V_z*P_y) - (V_y*P_z))
    except:
        L_x = yt.YTArray(0.0, 'g*cm**2/s')
    return L_x
    
yt.add_field("Orbital_Angular_Momentum_x", function=_Orbital_Angular_Momentum_x, units=r"g*cm**2/s")

def _Orbital_Angular_Momentum_y(field, data):
    """
    Calculates the orbital angular momentume of all the necessaru sink particles
    """
    global centred_sink_id
    global center_vel_part
    global center_pos_part
    global active_radius
    try:
        dd = data.ds.all_data()
        if np.isnan(active_radius):
            usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            usable_tags = np.array(usable_tags)
        else:
            centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
            dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
            dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
            dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
            dist = np.sqrt(dx**2+dy**2+dz**2)
            usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
        M = dd['sink_particle_mass'][np.array(usable_tags)].in_units('g')
        V_z = dd['sink_particle_velz'][usable_tags].in_units('cm/s') - center_vel_part[2].in_units('cm/s')
        V_x = dd['sink_particle_velx'][usable_tags].in_units('cm/s') - center_vel_part[0].in_units('cm/s')
        P_x = dd['sink_particle_posx'][usable_tags].in_units('cm') - center_pos_part[0].in_units('cm')
        P_z = dd['sink_particle_posz'][usable_tags].in_units('cm') - center_pos_part[2].in_units('cm')
        L_y = M *((V_z*P_x) - (V_x*P_z))
    except:
        L_y = yt.YTArray(0.0, 'g*cm**2/s')
    return L_y
    
yt.add_field("Orbital_Angular_Momentum_y", function=_Orbital_Angular_Momentum_y, units=r"g*cm**2/s")

def _Orbital_Angular_Momentum_z(field, data):
    """
    Calculates the orbital angular momentume of all the necessaru sink particles
    """
    global centred_sink_id
    global center_vel_part
    global center_pos_part
    global active_radius
    try:
        dd = data.ds.all_data()
        if np.isnan(active_radius):
            usable_tags = dd['sink_particle_tag'][centred_sink_id:].astype(int)
            usable_tags = np.array(usable_tags)
        else:
            centered_sink_pos = yt.YTArray([dd['sink_particle_posx'][sink_ind].in_units('au').value, dd['sink_particle_posy'][sink_ind].in_units('au').value, dd['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
            dx = dd['sink_particle_posx'].in_units('au') - center_pos[0]
            dy = dd['sink_particle_posy'].in_units('au') - center_pos[1]
            dz = dd['sink_particle_posz'].in_units('au') - center_pos[2]
            dist = np.sqrt(dx**2+dy**2+dz**2)
            usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
        M = dd['sink_particle_mass'][np.array(usable_tags)].in_units('g')
        V_y = dd['sink_particle_vely'][usable_tags].in_units('cm/s') - center_vel_part[1].in_units('cm/s')
        V_x = dd['sink_particle_velx'][usable_tags].in_units('cm/s') - center_vel_part[0].in_units('cm/s')
        P_x = dd['sink_particle_posx'][usable_tags].in_units('cm') - center_pos_part[0].in_units('cm')
        P_y = dd['sink_particle_posy'][usable_tags].in_units('cm') - center_pos_part[1].in_units('cm')
        L_z = M *((V_y*P_x) - (V_x*P_y))
    except:
        L_z = yt.YTArray(0.0, 'g*cm**2/s')
    return L_z
    
yt.add_field("Orbital_Angular_Momentum_z", function=_Orbital_Angular_Momentum_z, units=r"g*cm**2/s")

def _Orbital_Angular_Momentum(field, data):
    """
    returns the total orbital angular momentum
    """
    L = np.sqrt(data['Orbital_Angular_Momentum_x']**2 + data['Orbital_Angular_Momentum_y']**2 + data['Orbital_Angular_Momentum_z']**2)
    return L
    
yt.add_field("Orbital_Angular_Momentum", function=_Orbital_Angular_Momentum, units=r"g*cm**2/s")

def _Bulk_Velocity_Gas(field, data):
    """
    Calculates the mass weighted bulk velocity of the gas
    """
    TM = np.sum(data['cell_mass'].in_units('g'))
    x_top = np.sum(data['cell_mass'].in_units('g')*data['x-velocity'].in_units('cm/s'))
    y_top = np.sum(data['cell_mass'].in_units('g')*data['y-velocity'].in_units('cm/s'))
    z_top = np.sum(data['cell_mass'].in_units('g')*data['z-velocity'].in_units('cm/s'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    center_vel_gas = com
    return com

yt.add_field("Bulk_Velocity_Gas", function=_Bulk_Velocity_Gas, units=r"cm/s")

def _Radial_Velocity_x(field, data):
    """
    returns the velocity along the light of sight
    """
    global normal
    v = yt.YTArray([data['x-velocity'].in_units('cm/s').value, data['y-velocity'].in_units('cm/s').value, data['z-velocity'].in_units('cm/s').value], 'cm/s')
    proj_v_onto_L_x = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[0]
    proj_v_onto_L_y = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[1]
    proj_v_onto_L_z = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[2]
    proj_v_onto_L = yt.YTArray([proj_v_onto_L_x.T, proj_v_onto_L_y.T, proj_v_onto_L_z.T], 'cm/s')[0]
    return proj_v_onto_L

yt.add_field("Radial_Velocity_x", function=_Radial_Velocity_x, units=r"cm/s")

def _Radial_Velocity_y(field, data):
    """
    returns the x-component of the velocity along the light of sight
    """
    global normal
    v = yt.YTArray([data['x-velocity'].in_units('cm/s').value, data['y-velocity'].in_units('cm/s').value, data['z-velocity'].in_units('cm/s').value], 'cm/s')
    proj_v_onto_L_x = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[0]
    proj_v_onto_L_y = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[1]
    proj_v_onto_L_z = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[2]
    proj_v_onto_L = yt.YTArray([proj_v_onto_L_x.T, proj_v_onto_L_y.T, proj_v_onto_L_z.T], 'cm/s')[1]
    return proj_v_onto_L

yt.add_field("Radial_Velocity_y", function=_Radial_Velocity_y, units=r"cm/s")

def _Radial_Velocity_z(field, data):
    """
    returns the y-component of the velocity along the light of sight
    """
    global normal
    v = yt.YTArray([data['x-velocity'].in_units('cm/s').value, data['y-velocity'].in_units('cm/s').value, data['z-velocity'].in_units('cm/s').value], 'cm/s')
    proj_v_onto_L_x = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[0]
    proj_v_onto_L_y = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[1]
    proj_v_onto_L_z = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[2]
    proj_v_onto_L = yt.YTArray([proj_v_onto_L_x.T, proj_v_onto_L_y.T, proj_v_onto_L_z.T], 'cm/s')[1]
    return proj_v_onto_L

yt.add_field("Radial_Velocity_z", function=_Radial_Velocity_z, units=r"cm/s")

def _Projected_Velocity_x(field, data):
    """
    returns the z-component of the velocity along the light of sight
    """
    projected_velocity = data['x-velocity'].in_units('cm/s') - data['Radial_Velocity_x'].in_units('cm/s')
    return projected_velocity
    
yt.add_field("Projected_Velocity_x", function=_Projected_Velocity_x, units=r"cm/s")

def _Projected_Velocity_y(field, data):
    """
    returns the y-component of the projected velocity onto the planet defined by the normal
    """
    projected_velocity = data['y-velocity'].in_units('cm/s') - data['Radial_Velocity_y'].in_units('cm/s')
    return projected_velocity
    
yt.add_field("Projected_Velocity_y", function=_Projected_Velocity_y, units=r"cm/s")

def _Projected_Velocity_z(field, data):
    """
    returns the z-component of the projected velocity onto the planet defined by the normal
    """
    projected_velocity = data['z-velocity'].in_units('cm/s') - data['Radial_Velocity_z'].in_units('cm/s')
    return projected_velocity
    
yt.add_field("Projected_Velocity_z", function=_Projected_Velocity_z, units=r"cm/s")

def _Projected_Magnetic_Field_x(field, data):
    """
    returns the projected velocity
    """
    global normal
    v = yt.YTArray([data['magx'].value, data['magy'].value, data['magz'].value], 'G')
    proj_v_onto_L_x = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[0]
    proj_v_onto_L_y = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[1]
    proj_v_onto_L_z = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[2]
    proj_v_onto_L = yt.YTArray([proj_v_onto_L_x.T, proj_v_onto_L_y.T, proj_v_onto_L_z.T], 'G')
    projected_B = (v - proj_v_onto_L)[0]
    return projected_B
    
yt.add_field("Projected_Magnetic_Field_x", function=_Projected_Magnetic_Field_x, units=r"gauss")

def _Projected_Magnetic_Field_y(field, data):
    """
    returns the projected velocity
    """
    global normal
    v = yt.YTArray([data['magx'].value, data['magy'].value, data['magz'].value], 'G')
    proj_v_onto_L_x = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[0]
    proj_v_onto_L_y = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[1]
    proj_v_onto_L_z = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[2]
    proj_v_onto_L = yt.YTArray([proj_v_onto_L_x.T, proj_v_onto_L_y.T, proj_v_onto_L_z.T], 'G')
    projected_B = (v - proj_v_onto_L)[1]
    return projected_B
    
yt.add_field("Projected_Magnetic_Field_y", function=_Projected_Magnetic_Field_y, units=r"gauss")

def _Projected_Magnetic_Field_z(field, data):
    """
    returns the projected velocity
    """
    global normal
    v = yt.YTArray([data['magx'].value, data['magy'].value, data['magz'].value], 'G')
    proj_v_onto_L_x = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[0]
    proj_v_onto_L_y = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[1]
    proj_v_onto_L_z = (np.dot(v.T, normal)/np.dot(normal,normal))*normal[2]
    proj_v_onto_L = yt.YTArray([proj_v_onto_L_x.T, proj_v_onto_L_y.T, proj_v_onto_L_z.T], 'G')
    projected_B = (v - proj_v_onto_L)[2]
    return projected_B
    
yt.add_field("Projected_Magnetic_Field_z", function=_Projected_Magnetic_Field_z, units=r"gauss")

def _Projected_Particle_Posx(field, data):
    """
    returns the projected x position of particles
    """
    global centred_sink_id
    global normal
    global center_pos
    global active_radius
    if np.isnan(active_radius):
        usable_tags = data['sink_particle_tag'][centred_sink_id:].astype(int)
        usable_tags = np.array(usable_tags)
    else:
        centered_sink_pos = yt.YTArray([data['sink_particle_posx'][sink_ind].in_units('au').value, data['sink_particle_posy'][sink_ind].in_units('au').value, data['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
        dx = data['sink_particle_posx'].in_units('au') - center_pos[0]
        dy = data['sink_particle_posy'].in_units('au') - center_pos[1]
        dz = data['sink_particle_posz'].in_units('au') - center_pos[2]
        dist = np.sqrt(dx**2+dy**2+dz**2)
        usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
    try:
        dd = data.ds.all_data()
        pos = yt.YTArray([dd['sink_particle_posx'][usable_tags].in_units('cm').value - center_pos[0].in_units('cm').value, dd['sink_particle_posy'][usable_tags].in_units('cm').value - center_pos[1].in_units('cm').value, dd['sink_particle_posz'][usable_tags].in_units('cm').value - center_pos[2].in_units('cm').value], 'cm')
    except:
        pos = yt.YTArray([data['sink_particle_posx'][usable_tags].in_units('cm').value - center_pos[0].in_units('cm').value, data['sink_particle_posy'][usable_tags].in_units('cm').value - center_pos[1].in_units('cm').value, data['sink_particle_posz'][usable_tags].in_units('cm').value - center_pos[2].in_units('cm').value], 'cm')
    proj_pos_onto_L_x = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[0]
    proj_pos_onto_L_y = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[1]
    proj_pos_onto_L_z = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[2]
    proj_pos_onto_L = yt.YTArray([proj_pos_onto_L_x.T, proj_pos_onto_L_y.T, proj_pos_onto_L_z.T], 'cm')
    projected_pos = (pos - proj_pos_onto_L)[0]
    return projected_pos
    
yt.add_field("Projected_Particle_Posx", function=_Projected_Particle_Posx, units=r"cm")

def _Projected_Particle_Posy(field, data):
    """
    returns the projected y position of particles
    """
    global centred_sink_id
    global normal
    global center_pos
    global active_radius
    if np.isnan(active_radius):
        usable_tags = data['sink_particle_tag'][centred_sink_id:].astype(int)
        usable_tags = np.array(usable_tags)
    else:
        centered_sink_pos = yt.YTArray([data['sink_particle_posx'][sink_ind].in_units('au').value, data['sink_particle_posy'][sink_ind].in_units('au').value, data['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
        dx = data['sink_particle_posx'].in_units('au') - center_pos[0]
        dy = data['sink_particle_posy'].in_units('au') - center_pos[1]
        dz = data['sink_particle_posz'].in_units('au') - center_pos[2]
        dist = np.sqrt(dx**2+dy**2+dz**2)
        usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
    try:
        dd = data.ds.all_data()
        pos = yt.YTArray([dd['sink_particle_posx'][usable_tags].in_units('cm').value - center_pos[0].in_units('cm').value, dd['sink_particle_posy'][usable_tags].in_units('cm').value - center_pos[1].in_units('cm').value, dd['sink_particle_posz'][usable_tags].in_units('cm').value - center_pos[2].in_units('cm').value], 'cm')
    except:
        pos = yt.YTArray([data['sink_particle_posx'][usable_tags].in_units('cm').value - center_pos[0].in_units('cm').value, data['sink_particle_posy'][usable_tags].in_units('cm').value - center_pos[1].in_units('cm').value, data['sink_particle_posz'][usable_tags].in_units('cm').value - center_pos[2].in_units('cm').value], 'cm')
    proj_pos_onto_L_x = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[0]
    proj_pos_onto_L_y = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[1]
    proj_pos_onto_L_z = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[2]
    proj_pos_onto_L = yt.YTArray([proj_pos_onto_L_x.T, proj_pos_onto_L_y.T, proj_pos_onto_L_z.T], 'cm')
    projected_pos = (pos - proj_pos_onto_L)[1]
    return projected_pos
    
yt.add_field("Projected_Particle_Posy", function=_Projected_Particle_Posy, units=r"cm")

def _Projected_Particle_Posz(field, data):
    """
    returns the projected z position of particles
    """
    global centred_sink_id
    global normal
    global center_pos
    global active_radius
    if np.isnan(active_radius):
        usable_tags = data['sink_particle_tag'][centred_sink_id:].astype(int)
        usable_tags = np.array(usable_tags)
    else:
        centered_sink_pos = yt.YTArray([data['sink_particle_posx'][sink_ind].in_units('au').value, data['sink_particle_posy'][sink_ind].in_units('au').value, data['sink_particle_posz'][sink_ind].in_units('au').value], 'au')
        dx = data['sink_particle_posx'].in_units('au') - center_pos[0]
        dy = data['sink_particle_posy'].in_units('au') - center_pos[1]
        dz = data['sink_particle_posz'].in_units('au') - center_pos[2]
        dist = np.sqrt(dx**2+dy**2+dz**2)
        usable_tags = np.argwhere(dist.value < active_radius.value).T[0]
    try:
        dd = data.ds.all_data()
        pos = yt.YTArray([dd['sink_particle_posx'][usable_tags].in_units('cm').value - center_pos[0].in_units('cm').value, dd['sink_particle_posy'][usable_tags].in_units('cm').value - center_pos[1].in_units('cm').value, dd['sink_particle_posz'][usable_tags].in_units('cm').value - center_pos[2].in_units('cm').value], 'cm')
    except:
        pos = yt.YTArray([data['sink_particle_posx'][usable_tags].in_units('cm').value - center_pos[0].in_units('cm').value, data['sink_particle_posy'][usable_tags].in_units('cm').value - center_pos[1].in_units('cm').value, data['sink_particle_posz'][usable_tags].in_units('cm').value - center_pos[2].in_units('cm').value], 'cm')
    proj_pos_onto_L_x = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[0]
    proj_pos_onto_L_y = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[1]
    proj_pos_onto_L_z = (np.dot(pos.T, normal)/np.dot(normal,normal))*normal[2]
    proj_pos_onto_L = yt.YTArray([proj_pos_onto_L_x.T, proj_pos_onto_L_y.T, proj_pos_onto_L_z.T], 'cm')
    projected_pos = (pos - proj_pos_onto_L)[2]
    return projected_pos
    
yt.add_field("Projected_Particle_Posz", function=_Projected_Particle_Posz, units=r"cm")
