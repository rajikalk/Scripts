import pickle
import numpy as np
import yt
import gc
from mpi4py.MPI import COMM_WORLD as CW
import glob
import sys
import matplotlib.pyplot as plt

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
Time_arr = Time_arr-Time_arr[0]
del global_data['time']
gc.collect()
V_std_all = []
print('About to start calcualting RV spread')
rit = -1
for sink_id in range(len(global_data['m'].T)):
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=2)#, figsize=(two_col_width,page_height))
    plt.subplots_adjust(wspace=0.0)
    plt.subplots_adjust(hspace=0.07)
    
    axs[0].plot(Time_arr, global_data['m'].T[sink_id]*units['mass_unit'].in_units('msun'))
    axs[0].set_xlim(left=0)
    axs[0].set_ylim(bottom=0)
    axs[0].set_ylabel('Mass (Msun)')
    axs[0].axhline(y=5)
    
    axs[1].plot(Time_arr, global_data['ux'].T[sink_id]*units['velocity_unit'].in_units('km/s'))
    axs[1].set_xlim(left=0)
    axs[1].set_xlabel('Time (yr)')
    axs[1].set_ylabel('V_x (km/s)')
    plt.savefig('Sink_'+str(sink_id)+'_v_x_evol.png')
    print('plotted sink', sink_id)
    
    

