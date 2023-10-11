#!/usr/bin/env python
#Created by Rajika Kuruwita, 2019
import yt
from yt.fields.api import ValidateParameter
yt.enable_parallelism()
import numpy as np

def projected_vector(vector, proj_vector):
    """
    Calculates the projection of vecter projected onto vector
    """
    vector_units = vector.units
    if len(proj_vector)>3:
        #Calc vector.proj
        v_dot_pv = vector.T[0]*proj_vector.T[0] + vector.T[1]*proj_vector.T[1] + vector.T[2]*proj_vector.T[2]
        pv_dot_pv = proj_vector.T[0]**2 + proj_vector.T[1]**2 + proj_vector.T[2]**2
        proj_v_x = (v_dot_pv/pv_dot_pv)*proj_vector.T[0]
        proj_v_y = (v_dot_pv/pv_dot_pv)*proj_vector.T[1]
        proj_v_z = (v_dot_pv/pv_dot_pv)*proj_vector.T[2]
        proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T)
        del v_dot_pv, pv_dot_pv, proj_v_x, proj_v_y, proj_v_z
    else:
        proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
        proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
        proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
        proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
        del proj_v_x, proj_v_y, proj_v_z
    return proj_v

def _CoM_full(field, data):
    """
    Returns center of mass using both the gas and the sink particles
    """
    try:
        dd = data.ds.all_data()
        TM = np.sum(dd['gas', 'mass'].in_units('g'))
        x_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['flash', 'x'].in_units('cm'))
        y_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['flash', 'y'].in_units('cm'))
        z_top = np.sum(dd['gas', 'mass'].in_units('g')*dd['flash', 'z'].in_units('cm'))
        if ('all', 'particle_mass') in data.ds.field_list:
            TM = TM + np.sum(dd['particle_mass'].in_units('g'))
            x_top = x_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_posx'].in_units('cm'))
            y_top = y_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_posy'].in_units('cm'))
            z_top = z_top + np.sum(dd['all', 'particle_mass'].in_units('g')*dd['all', 'particle_posz'].in_units('cm'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
        com = yt.YTArray(com, 'cm')
        del x_top, y_top, z_top, TM
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
        del x_top, y_top, z_top, TM
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
        x_top = np.sum(data['gas', 'mass'].in_units('g')*data['flash', 'x'].in_units('cm'))
        y_top = np.sum(data['gas', 'mass'].in_units('g')*data['flash', 'y'].in_units('cm'))
        z_top = np.sum(data['gas', 'mass'].in_units('g')*data['flash', 'z'].in_units('cm'))
        if ('all', 'particle_mass') in data.ds.field_list:
            TM = TM + np.sum(data['all', 'particle_mass'].in_units('g'))
            x_top = x_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_posx'].in_units('cm'))
            y_top = y_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_posy'].in_units('cm'))
            z_top = z_top + np.sum(data['all', 'particle_mass'].in_units('g')*data['all', 'particle_posz'].in_units('cm'))
        com = [(x_top/TM), (y_top/TM), (z_top/TM)]
        com = yt.YTArray(com, 'cm')
        del x_top, y_top, z_top, TM
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
        del x_top, y_top, z_top, TM
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
    
    L_gas = data['mass'].value * np.cross(d_gas, dv_gas).T
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
    
    L_gas = data['mass'].value * np.cross(d_gas, dv_gas).T
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
    
    L_gas = data['mass'].value * np.cross(d_gas, dv_gas).T
    L_gas = yt.YTQuantity(L_gas, 'g*cm**2/s')
    return L_gas

yt.add_field("L_gas_z_wrt_CoM", function=_L_gas_z_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")
'''

def _L_gas_wrt_CoM_spec(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    CoM_pos = data['gas', 'CoM_full'].in_units('cm')
    CoM_vel = data['gas', 'CoM_Velocity_full'].in_units('cm/s')
    
    dx_gas = CoM_pos[0] - data['flash', 'x'].in_units('cm')
    dy_gas = CoM_pos[1] - data['flash', 'y'].in_units('cm')
    dz_gas = CoM_pos[2] - data['flash', 'z'].in_units('cm')
    d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
    del CoM_pos, dx_gas, dy_gas, dz_gas

    dvx_gas = CoM_vel[0] - data['flash','velx'].in_units('cm/s')
    dvy_gas = CoM_vel[1] - data['flash','vely'].in_units('cm/s')
    dvz_gas = CoM_vel[2] - data['flash','velz'].in_units('cm/s')
    d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
    del CoM_vel, dvx_gas, dvy_gas, dvz_gas
    
    L_gas = np.cross(d_pos_gas, d_vel_gas).T
    L_gas_wrt_CoM = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'cm**2/s')
    del L_gas
    return L_gas_wrt_CoM

yt.add_field("L_gas_wrt_CoM_spec", function=_L_gas_wrt_CoM_spec, units=r"cm**2/s", sampling_type="local")

def _L_gas_wrt_CoM_spec_cyl(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    CoM_pos = data['gas', 'CoM_full'].in_units('cm')
    CoM_vel = data['gas', 'CoM_Velocity_full'].in_units('cm/s')
    
    dx_gas = CoM_pos[0] - data['flash', 'x'].in_units('cm')
    dy_gas = CoM_pos[1] - data['flash', 'y'].in_units('cm')
    dz_gas = CoM_pos[2] - data['flash', 'z'].in_units('cm')
    d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas*0]).T
    del CoM_pos, dx_gas, dy_gas, dz_gas

    dvx_gas = CoM_vel[0] - data['flash','velx'].in_units('cm/s')
    dvy_gas = CoM_vel[1] - data['flash','vely'].in_units('cm/s')
    dvz_gas = CoM_vel[2] - data['flash','velz'].in_units('cm/s')
    d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas*0]).T
    del CoM_vel, dvx_gas, dvy_gas, dvz_gas
    
    L_gas = np.cross(d_pos_gas, d_vel_gas).T
    L_gas_wrt_CoM = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'cm**2/s')
    del L_gas
    return L_gas_wrt_CoM


yt.add_field("L_gas_wrt_CoM_spec_cyl", function=_L_gas_wrt_CoM_spec_cyl, units=r"cm**2/s", sampling_type="local")

def _L_gas_wrt_CoM(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    CoM_pos = data['gas', 'CoM_full'].in_units('cm')
    CoM_vel = data['gas', 'CoM_Velocity_full'].in_units('cm/s')
    
    dx_gas = CoM_pos[0] - data['flash', 'x'].in_units('cm')
    dy_gas = CoM_pos[1] - data['flash', 'y'].in_units('cm')
    dz_gas = CoM_pos[2] - data['flash', 'z'].in_units('cm')
    d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
    del CoM_pos, dx_gas, dy_gas, dz_gas

    dvx_gas = CoM_vel[0] - data['flash','velx'].in_units('cm/s')
    dvy_gas = CoM_vel[1] - data['flash','vely'].in_units('cm/s')
    dvz_gas = CoM_vel[2] - data['flash','velz'].in_units('cm/s')
    d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
    del CoM_vel, dvx_gas, dvy_gas, dvz_gas
    
    L_gas = data['gas', 'mass'].value * np.cross(d_pos_gas, d_vel_gas).T
    L_gas_wrt_CoM = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
    del L_gas
    return L_gas_wrt_CoM


yt.add_field("L_gas_wrt_CoM", function=_L_gas_wrt_CoM, units=r"g*cm**2/s", sampling_type="local")

def _L_gas_wrt_CoM_density(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    CoM_pos = data['gas', 'CoM_full'].in_units('cm')
    CoM_vel = data['gas', 'CoM_Velocity_full'].in_units('cm/s')
    
    dx_gas = CoM_pos[0] - data['flash', 'x'].in_units('cm')
    dy_gas = CoM_pos[1] - data['flash', 'y'].in_units('cm')
    dz_gas = CoM_pos[2] - data['flash', 'z'].in_units('cm')
    d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
    del CoM_pos, dx_gas, dy_gas, dz_gas

    dvx_gas = CoM_vel[0] - data['flash','velx'].in_units('cm/s')
    dvy_gas = CoM_vel[1] - data['flash','vely'].in_units('cm/s')
    dvz_gas = CoM_vel[2] - data['flash','velz'].in_units('cm/s')
    d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
    del CoM_vel, dvx_gas, dvy_gas, dvz_gas
    
    L_gas = data['flash', 'dens'].value * np.cross(d_pos_gas, d_vel_gas).T
    L_gas_wrt_CoM = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g/cm*s')
    del L_gas
    return L_gas_wrt_CoM

yt.add_field("L_gas_wrt_CoM_density", function=_L_gas_wrt_CoM_density, units=r"g/cm*s", sampling_type="local")

def _L_gas_wrt_CoM_cyl(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    CoM_pos = data['gas', 'CoM_full'].in_units('cm')
    CoM_vel = data['gas', 'CoM_Velocity_full'].in_units('cm/s')
    
    dx_gas = CoM_pos[0] - data['flash', 'x'].in_units('cm')
    dy_gas = CoM_pos[1] - data['flash', 'y'].in_units('cm')
    dz_gas = CoM_pos[2] - data['flash', 'z'].in_units('cm')
    d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas*0]).T
    del CoM_pos, dx_gas, dy_gas, dz_gas

    dvx_gas = CoM_vel[0] - data['flash','velx'].in_units('cm/s')
    dvy_gas = CoM_vel[1] - data['flash','vely'].in_units('cm/s')
    dvz_gas = CoM_vel[2] - data['flash','velz'].in_units('cm/s')
    d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas*0]).T
    del CoM_vel, dvx_gas, dvy_gas, dvz_gas
    
    L_gas = data['gas', 'mass'].value * np.cross(d_pos_gas, d_vel_gas).T
    L_gas_wrt_CoM = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
    del L_gas
    return L_gas_wrt_CoM


yt.add_field("L_gas_wrt_CoM_cyl", function=_L_gas_wrt_CoM_cyl, units=r"g*cm**2/s", sampling_type="local")

def _nearest_particle_index(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        d_all = []
        dd = data.ds.all_data()
        for part_pos_it in range(len(dd['particle_tag'])):
            dx_gas = dd['all', 'particle_posx'][part_pos_it].in_units('cm') - data['flash', 'x'].in_units('cm')
            dy_gas = dd['all', 'particle_posy'][part_pos_it].in_units('cm') - data['flash', 'y'].in_units('cm')
            dz_gas = dd['all', 'particle_posz'][part_pos_it].in_units('cm') - data['flash', 'z'].in_units('cm')
            d_gas = np.sqrt(dx_gas**2 + dy_gas**2 + dz_gas**2)
            del dx_gas, dy_gas, dz_gas
            d_all.append(d_gas)
            del d_gas
        del dd
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
        dx_gas = dd['all', 'particle_posx'][Nearest_tag_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][Nearest_tag_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][Nearest_tag_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
        del dx_gas, dy_gas, dz_gas
    
        dvx_gas = dd['all', 'particle_velx'][Nearest_tag_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][Nearest_tag_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][Nearest_tag_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        del dvx_gas, dvy_gas, dvz_gas, Nearest_tag_ind, dd
        
        L_gas = data['gas', 'mass'].value * np.cross(d_pos_gas, d_vel_gas).T
        L_wrt_nearest = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
        del L_gas
    else:
        L_wrt_nearest = data['gas', 'L_gas_wrt_CoM']
    return L_wrt_nearest

yt.add_field("L_gas_wrt_nearest_sink", function=_L_gas_wrt_nearest_sink, units=r"g*cm**2/s", sampling_type="local")

def _velx_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        del dd, primary_ind
    else:
        dvx_gas = data['flash','velx'].in_units('cm/s')
    return dvx_gas

yt.add_field("velx_wrt_primary", function=_velx_wrt_primary, units=r"cm/s", sampling_type="local")

def _vely_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dvx_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        del dd, primary_ind
    else:
        dvx_gas = data['flash','vely'].in_units('cm/s')
    return dvx_gas

yt.add_field("vely_wrt_primary", function=_vely_wrt_primary, units=r"cm/s", sampling_type="local")

def _velz_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dvx_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        del dd, primary_ind
    else:
        dvx_gas = data['flash','velz'].in_units('cm/s')
    return dvx_gas

yt.add_field("velz_wrt_primary", function=_velz_wrt_primary, units=r"cm/s", sampling_type="local")

def _L_gas_wrt_primary_spec(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
        del primary_ind, dx_gas, dy_gas, dz_gas
    
        dvx_gas = data['gas','velx_wrt_primary'].in_units('cm/s')
        dvy_gas = data['gas','vely_wrt_primary'].in_units('cm/s')
        dvz_gas = data['gas','velz_wrt_primary'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        del dvx_gas, dvy_gas, dvz_gas, dd
        
        L_gas = np.cross(d_pos_gas, d_vel_gas).T
        L_wrt_primary = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'cm**2/s')
        del L_gas
    else:
        L_wrt_primary = data['gas', 'L_gas_wrt_CoM_spec']
    return L_wrt_primary

yt.add_field("L_gas_wrt_primary_spec", function=_L_gas_wrt_primary_spec, units=r"cm**2/s", sampling_type="local")

def _L_gas_wrt_primary_spec_cyl(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas*0]).T
        del dx_gas, dy_gas, dz_gas
    
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas*0]).T
        del dvx_gas, dvy_gas, dvz_gas, primary_ind, dd
        
        L_gas = np.cross(d_pos_gas, d_vel_gas).T
        L_wrt_primary = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'cm**2/s')
        del L_gas
    else:
        L_wrt_primary = data['gas', 'L_gas_wrt_CoM_spec_cyl']
    return L_wrt_primary

yt.add_field("L_gas_wrt_primary_spec_cyl", function=_L_gas_wrt_primary_spec_cyl, units=r"cm**2/s", sampling_type="local")

def _L_gas_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
        del dx_gas, dy_gas, dz_gas
    
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        del dvx_gas, dvy_gas, dvz_gas, primary_ind, dd
        
        L_gas = data['gas', 'mass'].value * np.cross(d_pos_gas, d_vel_gas).T
        L_wrt_primary = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
        del L_gas
    else:
        L_wrt_primary = data['gas', 'L_gas_wrt_CoM']
    return L_wrt_primary

yt.add_field("L_gas_wrt_primary", function=_L_gas_wrt_primary, units=r"g*cm**2/s", sampling_type="local")

def _L_gas_wrt_primary_density(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
        del dx_gas, dy_gas, dz_gas
    
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        del dvx_gas, dvy_gas, dvz_gas, primary_ind, dd
        
        L_gas = data['flash', 'dens'].value * np.cross(d_pos_gas, d_vel_gas).T
        L_wrt_primary = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g/cm*s')
        del L_gas
    else:
        L_wrt_primary = data['gas', 'L_gas_wrt_CoM_density']
    return L_wrt_primary

yt.add_field("L_gas_wrt_primary_density", function=_L_gas_wrt_primary_density, units=r"g/cm*s", sampling_type="local")

def _L_gas_wrt_primary_cyl(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas*0]).T
        del dx_gas, dy_gas, dz_gas
    
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas*0]).T
        del dvx_gas, dvy_gas, dvz_gas, primary_ind, dd
        
        L_gas = data['gas', 'mass'].value * np.cross(d_pos_gas, d_vel_gas).T
        L_wrt_primary = yt.YTArray(np.sqrt(np.sum(L_gas**2, axis=0)), 'g*cm**2/s')
        del L_gas
    else:
        L_wrt_primary = data['gas', 'L_gas_wrt_CoM_cyl']
    return L_wrt_primary

yt.add_field("L_gas_wrt_primary_cyl", function=_L_gas_wrt_primary_cyl, units=r"g*cm**2/s", sampling_type="local")

def _Distance_from_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')

        radius = (dx_gas**2 + dy_gas**2 + dz_gas**2)**0.5
        del dd, primary_ind, dx_gas, dy_gas, dz_gas
    else:
        radius = yt.YTArray(np.zeros(np.shape(data['flash', 'x'])), 'cm')
    return radius

yt.add_field("Distance_from_primary", function=_Distance_from_primary, units=r"cm", sampling_type="local")

def _Distance_from_domain_centre(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    dx_gas = data['flash', 'x'].in_units('cm')
    dy_gas = data['flash', 'y'].in_units('cm')
    dz_gas = data['flash', 'z'].in_units('cm')

    radius = (dx_gas**2 + dy_gas**2 + dz_gas**2)**0.5
    del dx_gas, dy_gas, dz_gas
    return radius

yt.add_field("Distance_from_domain_centre", function=_Distance_from_domain_centre, units=r"cm", sampling_type="local")

def _cell_mass(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    volume = data['flash', 'dx'].in_units('cm') * data['flash', 'dy'].in_units('cm') * data['flash', 'dz'].in_units('cm')
    density = data['flash', 'dens'].in_units('g/cm**3')

    mass = density * volume
    del volume, density
    return mass

yt.add_field("cell_mass", function=_cell_mass, units=r"g", sampling_type="local")



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
    
    L_part = data['particle_mass'].value * np.cross(d_part, dv_part).T
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
    
    L_part = data['particle_mass'].value * np.cross(d_part, dv_part).T
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
    
    L_part = data['particle_mass'].value * np.cross(d_part, dv_part).T
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

def _Keplerian_velocity_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        v_kep = np.sqrt(abs(data['flash', 'gpot'].in_cgs())).in_units('cm/s')
    else:
        v_kep = yt.YTArray(np.ones(np.shape(data['gas', 'mass']))*np.nan, 'cm/s')
    return v_kep

yt.add_field("Keplerian_velocity_wrt_primary", function=_Keplerian_velocity_wrt_primary, units=r"cm/s", sampling_type="local")

def _Position_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        r_vec = yt.YTArray([dx_gas, dy_gas, dz_gas])
        del dx_gas, dy_gas, dz_gas, primary_ind, dd
    else:
        r_vec = yt.YTArray(np.ones(np.shape(data['gas', 'mass']))*np.nan, 'cm')
    return r_vec

yt.add_field("Position_wrt_primary", function=_Position_wrt_primary, units=r"cm", sampling_type="local")

def _Velocity_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        v_vec = yt.YTArray([dvx_gas, dvy_gas, dvz_gas])
        del dvx_gas, dvy_gas, dvz_gas, primary_ind, dd
    else:
        v_vec = yt.YTArray(np.ones(np.shape(data['gas', 'mass']))*np.nan, 'cm/s')
    return v_vec

yt.add_field("Velocity_wrt_primary", function=_Velocity_wrt_primary, units=r"cm/s", sampling_type="local")

def _Radial_velocity_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    #if np.shape(data['flash','velx']) == (16, 16, 16):
    #    rad_vel = yt.YTArray(np.ones(np.shape(data['flash','velx']))*np.nan, 'cm/s')
    #else:
    if ('all', 'particle_mass') in data.ds.field_list:
        #dd = data.ds.all_data()
        #primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        #dx_gas = dd['all', 'particle_posx'][primary_ind].in_units('cm') - data['flash', 'x'].in_units('cm')
        #dy_gas = dd['all', 'particle_posy'][primary_ind].in_units('cm') - data['flash', 'y'].in_units('cm')
        #dz_gas = dd['all', 'particle_posz'][primary_ind].in_units('cm') - data['flash', 'z'].in_units('cm')
        #r_vec = yt.YTArray([dx_gas, dy_gas, dz_gas])
        #del dx_gas, dy_gas, dz_gas
        #distance = np.sqrt(np.sum(r_vec**2, axis=0))
        #r_unit = r_vec/distance
        #del r_vec, distance
        
        #dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        #dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        #dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        #v_vec = yt.YTArray([dvx_gas, dvy_gas, dvz_gas])
        #del dvx_gas, dvy_gas, dvz_gas
        
        #rad_vel = v_vec[0]*r_unit[0] + v_vec[1]*r_unit[1] + v_vec[2]*r_unit[2]
        #del v_vec, r_unit
        #import pdb
        #pdb.set_trace()
        
        r_vec = data['Position_wrt_primary']
        v_vec = data['Velocity_wrt_primary']

        rad_vel = projected_vector(v_vec.T, r_vec.T)
        if np.shape(rad_vel) == (16, 16, 16, 3):
            rad_vel_mag = yt.YTArray(np.sqrt(np.sum(rad_vel**2, axis=3)).value, 'cm/s')
        else:
            rad_vel_mag = yt.YTArray(np.sqrt(np.sum(rad_vel**2, axis=1)).value, 'cm/s')
        del r_vec, v_vec, rad_vel
    else:
        rad_vel_mag = yt.YTArray(np.ones(np.shape(data['flash','velx']))*np.nan, 'cm/s')
    return rad_vel_mag

yt.add_field("Radial_velocity_wrt_primary", function=_Radial_velocity_wrt_primary, units=r"cm/s", sampling_type="local")

'''
def _Radial_velocity_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if np.shape(data['flash','velx']) == (16, 16, 16):
        rad_vel = yt.YTArray(np.ones(np.shape(data['flash','velx']))*np.nan, 'cm/s')
    else:
        if ('all', 'particle_mass') in data.ds.field_list:
            r_vec = data['Position_wrt_primary']
            v_vec = data['Velocity_wrt_primary']
            
            rad_vel = projected_vector(v_vec.T, r_vec.T)
            rad_vel = yt.YTArray(np.sqrt(np.sum(rad_vel**2, axis=1)).value, 'cm/s')
            del r_vec, v_vec
        else:
            rad_vel = yt.YTArray(np.ones(np.shape(data['flash','velx']))*np.nan, 'cm/s')
    return rad_vel

yt.add_field("Radial_velocity_wrt_primary", function=_Radial_velocity_wrt_primary, units=r"cm/s", sampling_type="local")
'''
def _Tangential_velocity_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        '''
        dd = data.ds.all_data()
        primary_ind = np.argmin(dd['all', 'particle_creation_time'])
        dvx_gas = dd['all', 'particle_velx'][primary_ind].in_units('cm/s') - data['flash','velx'].in_units('cm/s')
        dvy_gas = dd['all', 'particle_vely'][primary_ind].in_units('cm/s') - data['flash','vely'].in_units('cm/s')
        dvz_gas = dd['all', 'particle_velz'][primary_ind].in_units('cm/s') - data['flash','velz'].in_units('cm/s')
        v_mag_sq = dvx_gas**2 + dvy_gas**2 + dvz_gas**2
        del dd, primary_ind, dvx_gas, dvy_gas, dvz_gas
        '''
        
        v_vec = data['Velocity_wrt_primary']
        v_mag_sq = v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2
        
        rad_vel = data['Radial_velocity_wrt_primary'].in_units('cm/s')
        tang_vel = np.sqrt(v_mag_sq - rad_vel**2)
        del v_mag_sq, rad_vel
    else:
        tang_vel = yt.YTArray(np.ones(np.shape(data['flash','velx']))*np.nan, 'cm/s')
    return tang_vel

yt.add_field("Tangential_velocity_wrt_primary", function=_Tangential_velocity_wrt_primary, units=r"cm/s", sampling_type="local")

def _Relative_keplerian_velocity_wrt_primary(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        v_kep = np.sqrt(abs(data['flash', 'gpot'].in_cgs())).in_units('cm/s')
        #v_kep = data['Keplerian_velocity_wrt_primary']
        vel = data['Tangential_velocity_wrt_primary']
        rel_kep = vel/v_kep
        del v_kep, vel
    else:
        rel_kep = yt.YTArray(np.ones(np.shape(data['gas', 'mass']))*np.nan, '')
    return rel_kep

yt.add_field("Relative_keplerian_velocity_wrt_primary", function=_Relative_keplerian_velocity_wrt_primary, units=r"", sampling_type="local")

def _Relative_keplerian_velocity_wrt_primary_full_v(field, data):
    """
    Calculates the angular momentum w.r.t to the CoM
    """
    if ('all', 'particle_mass') in data.ds.field_list:
        v_kep = np.sqrt(abs(data['flash', 'gpot'].in_cgs())).in_units('cm/s')
        #v_kep = data['Keplerian_velocity_wrt_primary']
        vel_mag = np.sqrt(np.sum(data['Velocity_wrt_primary'][0].in_units('cm/s')**2 + data['Velocity_wrt_primary'][1].in_units('cm/s')**2 + data['Velocity_wrt_primary'][2].in_units('cm/s')**2))
        rel_kep = vel_mag/v_kep
        del v_kep, vel_mag
    else:
        rel_kep = yt.YTArray(np.ones(np.shape(data['gas', 'mass']))*np.nan, '')
    return rel_kep

yt.add_field("Relative_keplerian_velocity_wrt_primary_full_v", function=_Relative_keplerian_velocity_wrt_primary_full_v, units=r"", sampling_type="local")

def _Neg_z(field, data):
    """
    returns the negative of the z-positions
    """
    return -1*data['flash', 'z']

yt.add_field("Neg_z", function=_Neg_z, units=r"cm", sampling_type="local")
    
def _Neg_dz(field, data):
    """
    returns the negative of the dz
    """
    return -1*data['flash', 'dz']

yt.add_field("Neg_dz", function=_Neg_z, units=r"cm", sampling_type="local")

def _N_cells(field, data):
    """
    returns the negative of the dz
    """
    return np.ones(np.shape(data['flash', 'x']))

yt.add_field("N_cells", function=_N_cells, units=r"", sampling_type="local")
