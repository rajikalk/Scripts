import pickle
import numpy as np
import matplotlib.pyplot as plt
import yt
import gc
from mpi4py.MPI import COMM_WORLD as CW
import glob
import sys

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

#Set units
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('.')[0][1:]
try:
    simulation_density_id_int = int(simulation_density_id)
except:
    simulation_density_id = args.global_data_pickle_file.split('_G')[-1].split('.pkl')[0]

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/(units_override['length_unit'][0]**3), "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
gc.collect()

file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()


sys.stdout.flush()
CW.Barrier()

window = yt.YTQuantity(100, 'yr')
Sink_masses = {}
Sink_sigma_v = {}
Sink_delta_v = {}
rit = -1
for sink_id in range(len(global_data['ux'].T)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        Sink_masses.update({str(sink_id):[]})
        Sink_sigma_v.update({str(sink_id):[]})
        Sink_delta_v.update({str(sink_id):[]})
        print('Doing sink', sink_id, 'on rank', rank)
        for time_it in range(len(global_data['time'])):
            curr_mass = global_data['m'].T[sink_id][time_it]*units['mass_unit']
            Sink_masses[str(sink_id)].append(curr_mass)
        
            curr_time = global_data['time'][time_it]*scale_t.in_units('yr')
            start_time = curr_time - window/2
            end_time = curr_time + window/2
            
            start_ind = np.argmin(abs(global_data['time']*scale_t.in_units('yr') - start_time))
            end_ind = np.argmin(abs(global_data['time']*scale_t.in_units('yr') - end_time))
            
            if end_ind == start_ind == 0:
                std = yt.YTQuantity(np.nan, 'km/s')
                dv = yt.YTQuantity(np.nan, 'km/s')
            elif end_ind == start_ind:
                start_ind = start_ind -1
                std = np.std(global_data['ux'].T[sink_id][start_ind:end_ind+1]*scale_v.in_units('km/s'))
                dv = (np.max(global_data['ux'].T[sink_id][start_ind:end_ind+1]) - np.min(global_data['ux'].T[sink_id][start_ind:end_ind+1]))*scale_v.in_units('km/s')
            else:
                std = np.std(global_data['ux'].T[sink_id][start_ind:end_ind+1]*scale_v.in_units('km/s'))
                dv = (np.max(global_data['ux'].T[sink_id][start_ind:end_ind+1]) - np.min(global_data['ux'].T[sink_id][start_ind:end_ind+1]))*scale_v.in_units('km/s')
            Sink_sigma_v[str(sink_id)].append(std)
            Sink_delta_v[str(sink_id)].append(dv)
            if np.remainder(time_it, 1000) == 0:
                print('Calculate v_spread of sink', sink_id, 'for time_it', time_it, 'out of', len(global_data['time']))
        print('Calculated sigma_v evolution for sink', sink_id, 'on rank', rank)
        plt.clf()
        plt.plot((global_data['time']-global_data['time'][0])*scale_t.in_units('Myr'), Sink_delta_v[str(sink_id)], 'b-', label='Delta V')
        plt.plot((global_data['time']-global_data['time'][0])*scale_t.in_units('Myr'), Sink_sigma_v[str(sink_id)], 'r-', label='Sigma V')
        plt.xlabel('Time (Myr)')
        plt.ylabel('V spread (km/s)')
        plt.legend()
        #plt.xlim(left=0)
        #plt.ylim(bottom=0)
        plt.savefig('v_spread_vs_time_'+str(sink_id)+'.png')
        
#save pickle of v_spread over time for each sink
pickle_file_rank = 'V_spread_'+str(rank)+'.pkl'
file = open(pickle_file_rank, 'wb')
pickle.dump((Sink_masses, Sink_sigma_v, Sink_delta_v),file)
file.close()

#compile pickles
if rank == 0:
    pickle_files = glob.glob('V_spread_*.pkl')
    Sink_masses_all = {}
    Sink_sigma_v_all = {}
    Sink_delta_v_all = {}
    for pickle_file in pickle_files:
        file = open(pickle_file, 'rb')
        Sink_masses, Sink_sigma_v, Sink_delta_v = pickle.load(file)
        file.close()
        
        for key in Sink_masses.keys():
            Sink_masses_all.update({key:Sink_masses[key]})
            Sink_sigma_v_all.update({key:Sink_sigma_v[key]})
            Sink_delta_v_all.update({key:Sink_delta_v[key]})
        
    file = open('V_spread.pkl', 'wb')
    pickle.dump((Sink_masses_all, Sink_sigma_v_all, Sink_delta_v_all),file)
    file.close()
