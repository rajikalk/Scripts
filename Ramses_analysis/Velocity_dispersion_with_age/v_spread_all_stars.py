import pickle
import numpy as np
import yt
import gc
from mpi4py.MPI import COMM_WORLD as CW
import glob
import sys

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-window", "--intergration_window", default=100, type=float)
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

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
gc.collect()

if rank == 0:
    file_open = open(args.global_data_pickle_file, 'rb')
    try:
        global_data = pickle.load(file_open,encoding="latin1")
    except:
        file_open.close()
        import pickle5 as pickle
        file_open = open(args.global_data_pickle_file, 'rb')
        global_data = pickle.load(file_open,encoding="latin1")
    file_open.close()

    print('read in global data')
    del global_data['x']
    del global_data['y']
    del global_data['z']
    del global_data['uy']
    del global_data['uz']
    gc.collect()
    
    file = open('global_data_reduced.pkl', 'wb')
    pickle.dump((global_data), file)
    file.close()

sys.stdout.flush()
CW.Barrier()

file = open('global_data_reduced.pkl', 'rb')
global_data = pickle.load(file)
file.close()

sys.stdout.flush()
CW.Barrier()

convective_boundary = 0.2
intermediate_mass = 5
high_mass = 8

window = yt.YTQuantity(args.intergration_window, 'yr')
Time_arr = global_data['time']*units['time_unit'].in_units('yr')
del global_data['time']
gc.collect()
V_std_all = []
print('About to start calcualting RV spread')
rit = -1
for time_it in range(len(Time_arr)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        #get indexes of integration window
        curr_time = Time_arr[time_it]
        start_time = curr_time - window/2
        end_time = curr_time + window/2
        
        start_it = np.argmin(abs(Time_arr - start_time))
        end_it = np.argmin(abs(Time_arr - end_time))
        
        #iterate over the stars to measure RV dispersion over window
        
        #First find usable stars
        exisitng_stars = np.argwhere(global_data['m'][time_it]>0).T[0]
        V_spread_arr = []
        for star_it in exisitng_stars:
            vx_vals = global_data['ux'].T[star_it][start_it:end_it+1]*units['velocity_unit'].in_units('km/s')
            dv = np.max(vx_vals) - np.min(vx_vals)
            V_spread_arr.append(dv)
            
        #now that you have the Delta v, lets calculate the spread. Is the spread jsut the standard deviation?
        dv_std_all = np.mean(V_spread_arr)
        #I can also filter by mass
        Mass_arr = global_data['m'][time_it]*units['mass_unit'].in_units('Msun')
        Mass_convective_inds = np.argwhere(Mass_arr > convective_boundary)
        if len(Mass_convective_inds) == 0:
            dv_std_conv = np.nan
        else:
            dv_std_conv = np.mean(np.array(V_spread_arr)[Mass_convective_inds.T[0]])
            
        Mass_intermediate_inds = np.argwhere(Mass_arr > intermediate_mass)
        if len(Mass_intermediate_inds) == 0:
            dv_std_inter = np.nan
        else:
            dv_std_inter = np.mean(np.array(V_spread_arr)[Mass_intermediate_inds.T[0]])
            
        Mass_high_inds = np.argwhere(Mass_arr > high_mass)
        if len(Mass_high_inds) == 0:
            dv_std_high = np.nan
        else:
            dv_std_high = np.mean(np.array(V_spread_arr)[Mass_high_inds.T[0]])
            
        #let's save all these spreads
        V_std_all.append([curr_time, dv_std_all, dv_std_conv, dv_std_inter, dv_std_high])
        
        pickle_file = 'v_spread_'+str(rank)+'.pkl'
        file = open(pickle_file, 'wb')
        pickle.dump((V_std_all), file)
        file.close()
        print('Calculated v_spread for time_it', time_it, 'of', len(Time_arr), 'on rank', rank)
        
sys.stdout.flush()
CW.Barrier()

#Compile together pickles
print('collecting pickles')
pickle_files = glob.glob('v_spread_*.pkl')
All_V_spread = []
for pickle_file in pickle_files:
    file = open(pickle_file, 'rb')
    V_std_all = pickle.load(file)
    file.close()
    
    if len(All_V_spread) == 0:
        All_V_spread = V_std_all
    else:
        All_V_spread = All_V_spread + V_std_all
    
sort_inds = np.argsort(np.array(All_V_spread).T[0])
All_V_spread = np.array(All_V_spread)[sort_inds]

file = open('v_spread.pkl', 'wb')
pickle.dump((V_std_all), file)
file.close()

import matplotlib.pyplot as plt
#Let's try plotting these!
T_arr = All_V_spread.T[0] - All_V_spread.T[0][0]
plt.clf()
plt.plot(T_arr, All_V_spread.T[1], label='all stars')
plt.plot(T_arr, All_V_spread.T[2], label='>0.2Msun')
plt.plot(T_arr, All_V_spread.T[3], label='>5Msun')
plt.plot(T_arr, All_V_spread.T[4], label='>8Msun')
plt.legend()
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel('Time since first sink formation (yr)')
plt.ylabel('Mean Delta RV (km/s)')
plt.savefig('sigma_v_vs_time.png')
print('created plot')

