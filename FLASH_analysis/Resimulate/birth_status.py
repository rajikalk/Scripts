import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
sink_evol_pickle = sys.argv[1]

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

Birth_data = {}

first_sink = True
for key in sink_data.keys():
    if first_sink == False:
        form_time = sink_data[key]['time'][0]
        other_keys = []
        other_pos = []
        other_vel = []
        other_mass = []
        for other_key in sink_data.keys():
            try:
                form_ind = np.where((sink_data[other_key]['time'] - form_time)==0)[0][0]
                other_keys.append(other_key)
                other_pos.append([sink_data[other_key]['posx'][form_ind], sink_data[other_key]['posy'][form_ind], sink_data[other_key]['posz'][form_ind]])
                other_vel.append([sink_data[other_key]['velx'][form_ind], sink_data[other_key]['vely'][form_ind], sink_data[other_key]['velz'][form_ind]])
                other_mass.append(sink_data[other_key]['mass'][form_ind])
            except:
                pass
        other_keys.pop()
        other_pos.pop()
        other_vel.pop()
        other_mass.pop()
        
        #Calculate closest sink
        curr_pos = yt.YTArray([sink_data[key]['posx'][0], sink_data[key]['posy'][0], sink_data[key]['posz'][0]], 'cm')
        curr_vel = yt.YTArray([sink_data[key]['velx'][0], sink_data[key]['vely'][0], sink_data[key]['velz'][0]], 'cm/s')
        curr_mass = yt.YTQuantity(sink_data[key]['mass'][0], 'g')
        other_keys = yt.YTArray(other_keys, '')
        other_pos = yt.YTArray(other_pos, 'cm')
        other_vel = yt.YTArray(other_vel, 'cm/s')
        other_mass = yt.YTArray(other_mass, 'g')
        
        dx = other_pos.T[0] - curr_pos[0]
        dy = other_pos.T[1] - curr_pos[1]
        dz = other_pos.T[2] - curr_pos[2]
        separations = np.sqrt(dx**2 + dy**2 + dz**2)
        min_sep = np.min(separations.in_units('au'))
        closest_sink_id = other_keys[np.argmin(separations)]
        #CALCULATE TIME BETWEEN FORMATION
        
        #calculate energy
        GPE = -1*(yt.units.gravitational_constant*curr_mass*other_mass)/separations
        dvx = other_vel.T[0] - curr_vel[0]
        dvy = other_vel.T[1] - curr_vel[1]
        dvz = other_vel.T[2] - curr_vel[2]
        rel_vel = np.sqrt(dvx**2 + dvy**2 + dvz**2)
        reduced_mass = (curr_mass*other_mass)/(curr_mass + other_mass)
        KE = 1/2 * reduced_mass * rel_vel**2
        E_tot = GPE + KE
        Bound_inds = np.argwhere(E_tot<0)
        if len(Bound_inds)>0:
            bound_sink_ids = other_keys[Bound_inds]
            bound_sink_ids = bound_sink_ids.value.tolist()
        else:
            bound_sink_ids = []
        
        birth_data = {'closest_separation': float(min_sep.value), 'closest_sink': str(closest_sink_id.value), 'bound_sink_ids': bound_sink_ids}
        
        Birth_data.update({key: birth_data})
    else:
        first_sink = False


