#!/usr/bin/env python
#Created by Rajika Kuruwita, 2019
import yt
from yt.fields.api import ValidateParameter
yt.enable_parallelism()
import numpy as np

def _CoM_full(field, data):
    """
    Returns center of mass using both the gas and the sink particles
    """
    try:
        dd = data.ds.all_data()
        TM = np.sum(dd['gas', 'mass'].in_units('g'))
        x_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['gas', 'x'].in_units('cm'))
        y_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['gas', 'y'].in_units('cm'))
        z_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['gas', 'z'].in_units('cm'))
        if ('all', 'particle_mass') in data.ds.field_list:
            TM = TM + np.sum(dd['particle_mass'].in_units('g'))
            x_top = x_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_posx'].in_units('cm'))
            y_top = y_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_posy'].in_units('cm'))
            z_top = z_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_posz'].in_units('cm'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
        com = yt.YTArray(com, 'cm')
        del x_top
        del y_top
        del z_top
        del TM
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm')
    return com

yt.add_field("CoM_full", function=_CoM_full, units=r"cm", sampling_type="local")

def _CoM_Velocity_full(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    try:
        dd = data.ds.all_data()
        TM = np.sum(dd['gas', 'mass'].in_units('g'))
        x_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['flash','velx'].in_units('cm/s'))
        y_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['flash','vely'].in_units('cm/s'))
        z_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['flash','velz'].in_units('cm/s'))
        if ('all', 'particle_mass') in data.ds.field_list:
            TM = TM + np.sum(dd['all', 'particle_mass'].in_units('g'))
            x_top = x_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_velx'].in_units('cm/s'))
            y_top = y_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_vely'].in_units('cm/s'))
            z_top = z_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_velz'].in_units('cm/s'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
        com = yt.YTArray(com, 'cm/s')
        del x_top
        del y_top
        del z_top
        del TM
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
    return com

yt.add_field("CoM_Velocity_full", function=_CoM_Velocity_full, units=r"cm/s", sampling_type="local")

def _CoM(field, data):
    """
    Returns center of mass using both the gas and the sink particles
    """
    try:
        TM = np.sum(data['gas', 'mass'].in_units('g'))
        x_top = np.sum(data['gas', 'mass'].in_units('g')*data['gas', 'x'].in_units('cm'))
        y_top = np.sum(data['gas', 'mass'].in_units('g')*data['gas', 'y'].in_units('cm'))
        z_top = np.sum(data['gas', 'mass'].in_units('g')*data['gas', 'z'].in_units('cm'))
        if ('all', 'particle_mass') in data.ds.field_list:
            TM = TM + np.sum(data['all', 'particle_mass'].in_units('g'))
            x_top = x_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_posx'].in_units('cm'))
            y_top = y_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_posy'].in_units('cm'))
            z_top = z_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_posz'].in_units('cm'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
        com = yt.YTArray(com, 'cm')
        del x_top
        del y_top
        del z_top
        del TM
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm')
    return com

yt.add_field("CoM", function=_CoM, units=r"cm", sampling_type="local")

def _CoM_Velocity(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    try:
        TM = np.sum(data['gas', 'mass'].in_units('g'))
        x_top = np.sum(data['gas', 'mass'].in_units('g')*data['flash','velx'].in_units('cm/s'))
        y_top = np.sum(data['gas', 'mass'].in_units('g')*data['flash','vely'].in_units('cm/s'))
        z_top = np.sum(data['gas', 'mass'].in_units('g')*data['flash','velz'].in_units('cm/s'))
        if ('all', 'particle_mass') in data.ds.field_list:
            TM = TM + np.sum(data['particle_mass'].in_units('g'))
            x_top = x_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_velx'].in_units('cm/s'))
            y_top = y_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_vely'].in_units('cm/s'))
            z_top = z_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_velz'].in_units('cm/s'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
        com = yt.YTArray(com, 'cm/s')
        del x_top
        del y_top
        del z_top
        del TM
    except:
        com = yt.YTArray([0.0, 0.0, 0.0], 'cm/s')
    return com

yt.add_field("CoM_Velocity", function=_CoM_Velocity, units=r"cm/s", sampling_type="local")
'''
def _L_gas_x_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    CoM_pos = data['CoM_full'].in_units('cm')
    CoM_vel = data['CoM_Velocity_full'].in_units('cm/s')
    d_gas = data['x'].in_units('cm') - CoM_pos[0]
    
    dv_gas = data['velx'].in_units('cm/s') - CoM_vel[0]
    
    L_gas = data['mass'].value * np.cross(dv_gas, d_gas).T
    L_gas = yt.YTQuantity(L_gas, 'g*cm**2/s')
    return L_gas

yt.add_field("L_gas_x_wrt_CoM", function=_L_gas_x_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _L_gas_y_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    CoM_pos = data['CoM_full'].in_units('cm')
    CoM_vel = data['CoM_Velocity_full'].in_units('cm/s')
    d_gas = data['y'].in_units('cm') - CoM_pos[1]
    
    dv_gas = data['vely'].in_units('cm/s') - CoM_vel[1]
    
    L_gas = data['mass'].value * np.cross(dv_gas, d_gas).T
    L_gas = yt.YTQuantity(L_gas, 'g*cm**2/s')
    return L_gas

yt.add_field("L_gas_y_wrt_CoM", function=_L_gas_y_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _L_gas_z_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    CoM_pos = data['CoM_full'].in_units('cm')
    CoM_vel = data['CoM_Velocity_full'].in_units('cm/s')
    d_gas = data['z'].in_units('cm') - CoM_pos[2]
    
    dv_gas = data['velz'].in_units('cm/s') - CoM_vel[2]
    
    L_gas = data['mass'].value * np.cross(dv_gas, d_gas).T
    L_gas = yt.YTQuantity(L_gas, 'g*cm**2/s')
    return L_gas

yt.add_field("L_gas_z_wrt_CoM", function=_L_gas_z_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")
'''
def _L_gas_wrt_CoM(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    CoM_pos = data['gas', 'CoM_full'].in_units('cm')
    CoM_vel = data['gas', 'CoM_Velocity_full'].in_units('cm/s')
    
    dd = data.ds.all_data()
    dx_gas = CoM_pos[0] - data['gas', 'x'].in_units('cm')
    dy_gas = CoM_pos[1] - data['gas', 'y'].in_units('cm')
    dz_gas = CoM_pos[2] - data['gas', 'z'].in_units('cm')
    d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T

    dvx_gas = CoM_vel[0] - data['flash','velx'].in_units('cm/s')
    dvy_gas = CoM_vel[1] - data['flash','vely'].in_units('cm/s')
    dvz_gas = CoM_vel[2] - data['flash','velz'].in_units('cm/s')
    d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
    
    L_gas = data['gas', 'mass'].value * np.cross(d_vel_gas, d_pos_gas).T
    L_gas_wrt_CoM = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
    return L_gas_wrt_CoM


yt.add_field("L_gas_wrt_CoM", function=_L_gas_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _nearest_particle_index(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        d_all = []
        dd = data.ds.all_data()
        for part_pos_it in range(len(dd['particle_tag'])):
            dx_gas = dd['all', 'particle_posx'][part_pos_it].in_units('cm') - data['gas', 'x'].in_units('cm')
            dy_gas = dd['all', 'particle_posy'][part_pos_it].in_units('cm') - data['gas', 'y'].in_units('cm')
            dz_gas = dd['all', 'particle_posz'][part_pos_it].in_units('cm') - data['gas', 'z'].in_units('cm')
            d_gas = np.sqrt(dx_gas**2 + dy_gas**2 + dz_gas**2)
            d_all.append(d_gas)
        #Nearest_tag = data['particle_tag'][np.argmin(d_all, axis=0)]
        Nearest_tag_ind = yt.YTArray(np.argmin(d_all, axis=0), '')
        
    else:
        Nearest_tag_ind = yt.YTArray(np.nan*np.ones(np.shape(data['gas', 'x'])), '')
    return Nearest_tag_ind

yt.add_field("nearest_particle_index", function=_nearest_particle_index, units=r"", sampling_type="local")

def _L_gas_wrt_nearest_sink(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        Nearest_tag_ind = data['gas', 'nearest_particle_index'].value.astype(int)
        dx_gas = dd['all', 'particle_posx'][Nearest_tag_ind].in_units('cm') - data['gas', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][Nearest_tag_ind].in_units('cm') - data['gas', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][Nearest_tag_ind].in_units('cm') - data['gas', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
    
        dvx_gas = dd['all', 'particle_velx'][Nearest_tag_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][Nearest_tag_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][Nearest_tag_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        
        L_gas = data['gas', 'mass'].value * np.cross(d_vel_gas, d_pos_gas).T
        L_wrt_nearest = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
    else:
        L_wrt_nearest = yt.YTArray(np.nan*np.ones(np.shape(data['x'])), 'g*cm**2/s')
    return L_wrt_nearest

yt.add_field("L_gas_wrt_nearest_sink", function=_L_gas_wrt_nearest_sink, units=r"g*cm**2/s", sampling_type="local")

'''
def _L_gas_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    L_tot = np.sqrt(data['L_gas_x_wrt_CoM']**2 + data['L_gas_y_wrt_CoM']**2 + data['L_gas_z_wrt_CoM']**2)
    return L_tot

yt.add_field("L_gas_wrt_CoM", function=_L_gas_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")


def _L_part_x_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    CoM_pos = data['CoM_full'].in_units('cm')
    CoM_vel = data['CoM_Velocity_full'].in_units('cm/s')
    d_part = data['particle_posx'].in_units('cm') - CoM_pos[0]
    
    dv_part = data['particle_velx'].in_units('cm/s') - CoM_vel[0]
    
    L_part = data['particle_mass'].value * np.cross(dv_part, d_part).T
    L_part = yt.YTQuantity(L_part, 'g*cm**2/s')
    return L_gas

yt.add_field("L_part_x_wrt_CoM", function=_L_part_x_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _L_part_y_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    CoM_pos = data['CoM_full'].in_units('cm')
    CoM_vel = data['CoM_Velocity_full'].in_units('cm/s')
    d_part = data['particle_posy'].in_units('cm') - CoM_pos[0]
    
    dv_part = data['particle_vely'].in_units('cm/s') - CoM_vel[0]
    
    L_part = data['particle_mass'].value * np.cross(dv_part, d_part).T
    L_part = yt.YTQuantity(L_part, 'g*cm**2/s')
    return L_gas

yt.add_field("L_part_y_wrt_CoM", function=_L_part_y_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _L_part_z_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    CoM_pos = data['CoM_full'].in_units('cm')
    CoM_vel = data['CoM_Velocity_full'].in_units('cm/s')
    d_part = data['particle_posz'].in_units('cm') - CoM_pos[0]
    
    dv_part = data['particle_velz'].in_units('cm/s') - CoM_vel[0]
    
    L_part = data['particle_mass'].value * np.cross(dv_part, d_part).T
    L_part = yt.YTQuantity(L_part, 'g*cm**2/s')
    return L_gas

yt.add_field("L_part_z_wrt_CoM", function=_L_part_z_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _L_part_wrt_CoM(field, data):
    """
    Calculates the velcoity fo the CoM
    """
    L_tot = np.sqrt(data['L_part_x_wrt_CoM']**2 + data['L_part_y_wrt_CoM']**2 + data['L_part_z_wrt_CoM']**2)
    return L_tot

yt.add_field("L_part_wrt_CoM", function=_L_part_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

'''
