#!/usr/bin/env python
#Created by Rajika Kuruwita, 2019
import yt
yt.enable_parallelism()
import numpy as np
#import csv
#import glob
from pyramses import rsink

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

yt.add_field("Density", function=_Density, units=r"g/cm**3", sampling_type="local", force_override=True)

def _cell_mass(field,data):
    """
    Overwrites cell mass field
    """
    cell_mass = data[('gas', 'Density')].in_units('g/cm**3')*data['cell_volume'].in_units('cm**3')
    return cell_mass

yt.add_field("cell_mass", function=_cell_mass, units=r"g", sampling_type="local", force_override=True)
    
#===========================CALCULATING RAMSES CENTERED MAGNETIC FIELDS======================================

def _magx(field,data):
    """
    Calculated the centred x-component of the megnatic field
    """
    try:
        magx = yt.YTArray((data['x-Bfield-left'].value + data['x-Bfield-right'].value)/2., "G")
    except:
        try:
            magx = yt.YTArray((data['B_x_left'].value + data['B_x_right'].value)/2., "G")
        except:
            magx = yt.YTArray((data['hydro_B_left_x'].value + data['hydro_B_right_x'].value)/2., "G")
    return magx
    
yt.add_field("magx", function=_magx, units=r"gauss", sampling_type="local", force_override=True)

def _magy(field,data):
    """
    Calculated the centred y-component of the megnatic field
    """
    try:
        magy = yt.YTArray((data['y-Bfield-left'].value + data['y-Bfield-right'].value)/2., "G")
    except:
        try:
            magy = yt.YTArray((data['B_y_left'].value + data['B_y_right'].value)/2., "G")
        except:
            magy = yt.YTArray((data['hydro_B_left_y'].value + data['hydro_B_right_y'].value)/2., "G")
    return magy
    
yt.add_field("magy", function=_magy, units=r"gauss", sampling_type="local", force_override=True)

def _magz(field,data):
    """
    Calculated the centred z-component of the megnatic field
    """
    try:
        magz = yt.YTArray((data['z-Bfield-left'].value + data['z-Bfield-right'].value)/2., "G")
    except:
        try:
            magz = yt.YTArray((data['B_z_left'].value + data['B_z_right'].value)/2., "G")
        except:
            magz = yt.YTArray((data['hydro_B_left_z'].value + data['hydro_B_right_z'].value)/2., "G")
    return magz
    
yt.add_field("magz", function=_magz, units=r"gauss", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_tag", function=_sink_particle_tag, units=r"", sampling_type="local", force_override=True)

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
        import pdb
        pdb.set_trace()
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

yt.add_field("sink_particle_posx", function=_sink_particle_posx, units=r"pc", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_posy", function=_sink_particle_posy, units=r"pc", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_posz", function=_sink_particle_posz, units=r"pc", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_velx", function=_sink_particle_velx, units=r"km/s", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_vely", function=_sink_particle_vely, units=r"km/s", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_velz", function=_sink_particle_velz, units=r"km/s", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_speed", function=_sink_particle_speed, units=r"km/s", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_mass", function=_sink_particle_mass, units=r"Msun", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_accretion_rate", function=_sink_particle_accretion_rate, units=r"Msun/yr", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_form_time", function=_sink_particle_form_time, units=r"yr", sampling_type="local", force_override=True)

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

yt.add_field("sink_particle_age", function=_sink_particle_age, units=r"yr", sampling_type="local", force_override=True)
