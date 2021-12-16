import yt
from pyramses import rsink
import numpy as np

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
        particle_age = yt.YTArray(np.array(particle_age), "yr")
    return particle_age

yt.add_field("sink_particle_age", function=_sink_particle_age, units=r"yr")
