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

convective_boundary = 0.2
intermediate_mass = 5
high_mass = 8

sys.stdout.flush()
CW.Barrier()

window = yt.YTQuantity(100, 'yr')
Time_arr = global_data['time']*units['time_unit'].in_units('yr')
V_std_all = []
for time_it in range(len(Time_arr)):
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
        vx_vals = global_data['ux'].T[star_it][start_it:end_it]*units['velocity_unit'].in_units('km/s')
        dv = np.max(vx_vals) - np.min(vx_vals)
        V_spread_arr.append(dv)
        
    #now that you have the Delta v, lets calculate the spread. Is the spread jsut the standard deviation?
    dv_std_all = np.std(V_spread_arr)
    #I can also filter by mass
    Mass_arr = global_data['m'][time_it]*units['mass_unit'].in_units('Msun')
    Mass_convective_inds = np.argwhere(Mass_arr > convective_boundary)
    if len(Mass_convective_inds) == 0:
        dv_std_conv = np.nan
    else:
        import pdb
        pdb.set_trace()
        
    Mass_intermediate_inds = np.argwhere(Mass_arr > intermediate_mass)
    if len(Mass_intermediate_inds) == 0:
        dv_std_inter = np.nan
    else:
        import pdb
        pdb.set_trace()
        
    Mass_high_inds = np.argwhere(Mass_arr > high_mass)
    if len(Mass_high_inds) == 0:
        dv_std_high = np.nan
    else:
        import pdb
        pdb.set_trace()
        
    #let's save all these spreads
    V_std_all.append([dv_std_all, dv_std_conv, dv_std_inter, dv_std_high])
        
import pdb
pdb.set_trace()


