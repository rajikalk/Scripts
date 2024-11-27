import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
sink_evol_pickle = sys.argv[1]

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

final_keys = []
final_pos = []
final_vel = []
final_mass = []

for final_key in sink_data.keys():
    final_keys.append(final_key)
    final_pos.append([sink_data[final_key]['posx'][-1], sink_data[final_key]['posy'][-1], sink_data[final_key]['posz'][-1]])
    final_vel.append([sink_data[final_key]['velx'][-1], sink_data[final_key]['vely'][-1], sink_data[final_key]['velz'][-1]])
    final_mass.append(sink_data[final_key]['mass'][-1])
    
final_pos = yt.YTArray(final_pos, 'cm')
final_vel = yt.YTArray(final_vel, 'cm/s')
final_mass = yt.YTArray(final_mass, 'g')

Final_data = {}

for key in sink_data.keys():
    curr_ind = final_keys.index(key)
    curr_mass = final_mass[curr_ind]
    curr_pos = final_pos[curr_ind]
    curr_vel = final_vel[curr_ind]
    
    dx = final_pos.T[0] - curr_pos[0]
    dy = final_pos.T[1] - curr_pos[1]
    dz = final_pos.T[2] - curr_pos[2]
    separations = np.sqrt(dx**2 + dy**2 + dz**2)
    if key == '65910':
        import pdb
        pdb.set_trace()
    min_ind = np.argsort(separations)[1]
    min_sep = separations.in_units('au')[min_ind]
    closest_sink_id = final_keys[min_ind]
    #CALCULATE TIME BETWEEN FORMATION
    
    #calculate energy
    GPE = -1*(yt.units.gravitational_constant*curr_mass*final_mass)/separations
    dvx = final_vel.T[0] - curr_vel[0]
    dvy = final_vel.T[1] - curr_vel[1]
    dvz = final_vel.T[2] - curr_vel[2]
    rel_vel = np.sqrt(dvx**2 + dvy**2 + dvz**2)
    reduced_mass = (curr_mass*final_mass)/(curr_mass + final_mass)
    KE = 1/2 * reduced_mass * rel_vel**2
    E_tot = GPE + KE
    E_tot[np.where(np.isinf(E_tot))[0]] = yt.YTArray([np.inf], 'cm**2*g/s**2')
    Bound_inds = np.argwhere(E_tot<0)
    if len(Bound_inds)>0:
        bound_sink_ids = yt.YTArray(final_keys, '')[Bound_inds]
        bound_sink_ids = bound_sink_ids.value.tolist()
    else:
        bound_sink_ids = []
    
    final_data = {'closest_separation': float(min_sep.value), 'closest_sink': closest_sink_id, 'bound_sink_ids': bound_sink_ids}
        
    Final_data.update({key: final_data})
