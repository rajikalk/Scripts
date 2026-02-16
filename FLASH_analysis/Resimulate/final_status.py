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
    #if key == '65910':
    #    import pdb
    #    pdb.set_trace()
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

Masses = []
Spec_angs = []
Final_Accretion = []
averaging_time = yt.YTQuantity(100, 'yr')

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=False)
for key in sink_data.keys():
    sim_time = yt.YTArray(sink_data[key]['time'], 's').in_units('yr')
    time = sink_data[key]['time'] - sink_data[key]['time'][0]
    time = yt.YTArray(time, 's').in_units('yr')
    mass = yt.YTArray(sink_data[key]['mass'], 'g')
    ang = np.sqrt(sink_data[key]['anglx']**2 + sink_data[key]['angly']**2 + sink_data[key]['anglz']**2)
    spec_ang = ang/mass.in_units('g').value
    axs[0].plot(sim_time, spec_ang)
    axs[1].plot(sim_time, mass.in_units('msun'))

axs[1].set_xlabel('Time since formation (yr)')
axs[0].set_ylabel('specific angular momentum')
axs[1].set_ylabel('mass (msun)')
#axs[0].set_xlim(left=0)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
plt.savefig('all_spec_ang.png')


for key in sink_data.keys():
    ang_total = np.sqrt(sink_data[key]['anglx'][-1]**2 + sink_data[key]['angly'][-1]**2 + sink_data[key]['anglz'][-1]**2)
    Mass_array = yt.YTArray(sink_data[key]['mass'], 'g')
    Masses.append(Mass_array[-1].in_units('msun'))
    Time_array = yt.YTArray(sink_data[key]['time'], 's')
    average_start = Time_array[-1] - averaging_time.in_units('s')
    average_start_ind = np.argmin(abs(Time_array - average_start))
    dM = Mass_array[-1] - Mass_array[average_start_ind]
    dt = Time_array[-1] - Time_array[average_start_ind]
    accretion = dM.in_units('Msun')/dt.in_units('yr')
    ang_spec = ang_total/sink_data[key]['mass'][-1]
    Spec_angs.append(ang_spec)
    Final_Accretion.append(accretion)
    
plt.clf()
plt.scatter(Masses, Spec_angs)
plt.xlabel('Mass (msun)')
plt.ylabel('Specific angular momentum')
plt.savefig('spin_distribution.png')

plt.clf()
plt.scatter(Final_Accretion, Spec_angs)
plt.xlabel('Final_Accretion (Msun/yr)')
plt.ylabel('Specific angular momentum')
plt.savefig('accretion_vs_spin.png')
    
