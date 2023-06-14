#!/usr/bin/env python
#Created by Rajika Kuruwita, 2019
import yt
yt.enable_parallelism()
import numpy as np

def _Neg_z(field, data):
    """
    returns the negative of the z-positions
    """
    return -1*data['z']

yt.add_field("Neg_z", function=_Neg_z, units=r"cm", sampling_type="local")
    
def _Neg_dz(field, data):
    """
    returns the negative of the dz
    """
    return -1*data['dz']
    
yt.add_field("Neg_dz", function=_Neg_z, units=r"cm", sampling_type="local")

def _CoM_full(field, data):
    """
    Returns center of mass using both the gas and the sink particles
    """
    dd = data.ds.all_data()
    TM = np.sum(dd['cell_mass'].in_units('g'))
    x_top = np.sum(dd['cell_mass'].in_units('g')*dd['x'].in_units('cm'))
    y_top = np.sum(dd['cell_mass'].in_units('g')*dd['y'].in_units('cm'))
    z_top = np.sum(dd['cell_mass'].in_units('g')*dd['z'].in_units('cm'))
    if ('all', 'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(dd['particle_mass'].in_units('g'))
        x_top = x_top + np.sum(dd['particle_mass'].in_units('g')*dd['particle_posx'].in_units('cm'))
        y_top = y_top + np.sum(dd['particle_mass'].in_units('g')*dd['particle_posy'].in_units('cm'))
        z_top = z_top + np.sum(dd['particle_mass'].in_units('g')*dd['particle_posz'].in_units('cm'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    com = yt.YTArray(com, 'cm')
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM_full", function=_CoM, units=r"cm", sampling_type="local")

def _CoM_Velocity_full(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    dd = data.ds.all_data()
    TM = np.sum(dd['cell_mass'].in_units('g'))
    x_top = np.sum(dd['cell_mass'].in_units('g')*dd['velx'].in_units('cm/s'))
    y_top = np.sum(dd['cell_mass'].in_units('g')*dd['vely'].in_units('cm/s'))
    z_top = np.sum(dd['cell_mass'].in_units('g')*dd['velz'].in_units('cm/s'))
    if ('all', 'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(dd['particle_mass'].in_units('g'))
        x_top = x_top + np.sum(dd['particle_mass'].in_units('g')*dd['particle_velx'].in_units('cm/s'))
        y_top = y_top + np.sum(dd['particle_mass'].in_units('g')*dd['particle_vely'].in_units('cm/s'))
        z_top = z_top + np.sum(dd['particle_mass'].in_units('g')*dd['particle_velz'].in_units('cm/s'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    com = yt.YTArray(com, 'cm/s')
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM_Velocity_full", function=_CoM_Velocity, units=r"cm/s", sampling_type="local")

def _CoM(field, data):
    """
    Returns center of mass using both the gas and the sink particles
    """
    TM = np.sum(data['cell_mass'].in_units('g'))
    x_top = np.sum(data['cell_mass'].in_units('g')*data['x'].in_units('cm'))
    y_top = np.sum(data['cell_mass'].in_units('g')*data['y'].in_units('cm'))
    z_top = np.sum(data['cell_mass'].in_units('g')*data['z'].in_units('cm'))
    if ('all', 'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(data['particle_mass'].in_units('g'))
        x_top = x_top + np.sum(data['particle_mass'].in_units('g')*data['particle_posx'].in_units('cm'))
        y_top = y_top + np.sum(data['particle_mass'].in_units('g')*data['particle_posy'].in_units('cm'))
        z_top = z_top + np.sum(data['particle_mass'].in_units('g')*data['particle_posz'].in_units('cm'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    com = yt.YTArray(com, 'cm')
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM", function=_CoM, units=r"cm", sampling_type="local")

def _CoM_Velocity(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    TM = np.sum(data['cell_mass'].in_units('g'))
    x_top = np.sum(data['cell_mass'].in_units('g')*data['velx'].in_units('cm/s'))
    y_top = np.sum(data['cell_mass'].in_units('g')*data['vely'].in_units('cm/s'))
    z_top = np.sum(data['cell_mass'].in_units('g')*data['velz'].in_units('cm/s'))
    if ('all', 'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(data['particle_mass'].in_units('g'))
        x_top = x_top + np.sum(data['particle_mass'].in_units('g')*data['particle_velx'].in_units('cm/s'))
        y_top = y_top + np.sum(data['particle_mass'].in_units('g')*data['particle_vely'].in_units('cm/s'))
        z_top = z_top + np.sum(data['particle_mass'].in_units('g')*data['particle_velz'].in_units('cm/s'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    com = yt.YTArray(com, 'cm/s')
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM_Velocity", function=_CoM_Velocity, units=r"cm/s", sampling_type="local")

def _L_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    TM = np.sum(data['cell_mass'].in_units('g'))
    x_top = np.sum(data['cell_mass'].in_units('g')*data['velx'].in_units('cm/s'))
    y_top = np.sum(data['cell_mass'].in_units('g')*data['vely'].in_units('cm/s'))
    z_top = np.sum(data['cell_mass'].in_units('g')*data['velz'].in_units('cm/s'))
    if ('all', 'particle_mass') in data.ds.field_list:
        TM = TM + np.sum(data['particle_mass'].in_units('g'))
        x_top = x_top + np.sum(data['particle_mass'].in_units('g')*data['particle_velx'].in_units('cm/s'))
        y_top = y_top + np.sum(data['particle_mass'].in_units('g')*data['particle_vely'].in_units('cm/s'))
        z_top = z_top + np.sum(data['particle_mass'].in_units('g')*data['particle_velz'].in_units('cm/s'))
    com = [(x_top/TM), (y_top/TM), (z_top/TM)]
    com = yt.YTArray(com, 'cm/s')
    del x_top
    del y_top
    del z_top
    del TM
    return com

yt.add_field("CoM_Velocity", function=_CoM_Velocity, units=r"cm/s", sampling_type="local")
