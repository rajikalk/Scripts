import pickle
import numpy as np
import matplotlib.pyplot as plt
import yt
import gc

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

#Set units
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('_')[0][1:]
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
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

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

Cluster_age = []
V_spread_all = []
V_spread_1 = []
V_spread_5 = []
V_spread_8 = []
for time_it in range(len(global_data['ux'])):
    #get indices fo stars that exist
    all_stars = np.argwhere(global_data['m'][time_it]>0).T[0]
    v_std_all_stars = np.std(global_data['ux'][time_it][all_stars]*scale_v.in_units('km/s'))
    V_spread_all.append(v_std_all_stars)
    
    M_1_stars = np.argwhere(global_data['m'][time_it]>1).T[0]
    v_std_1 = np.std(global_data['ux'][time_it][M_1_stars]*scale_v.in_units('km/s'))
    V_spread_1.append(v_std_1)
    
    M_5_stars = np.argwhere(global_data['m'][time_it]>1).T[0]
    v_std_5 = np.std(global_data['ux'][time_it][M_5_stars]*scale_v.in_units('km/s'))
    V_spread_5.append(v_std_5)
    
    M_8_stars = np.argwhere(global_data['m'][time_it]>1).T[0]
    v_std_8 = np.std(global_data['ux'][time_it][M_8_stars]*scale_v.in_units('km/s'))
    V_spread_8.append(v_std_8)
    
    age = (global_data['time'][time_it][0] - global_data['time'][0][0])*scale_t.in_units('Myr')
    Cluster_age.append(age)
    
    if np.remainder(time_it, 1000) == 0:
        print('calculated v_spread for time_it', time_it, 'of', len(global_data['ux']))
        

plt.clf()
plt.plot(Cluster_age, V_spread_all, label='all stars', alpha=0.5)
plt.plot(Cluster_age, V_spread_1, label='> 1Msun', alpha=0.5)
plt.plot(Cluster_age, V_spread_5, label='> 5Msun', alpha=0.5)
plt.plot(Cluster_age, V_spread_8, label='> 8Msun', alpha=0.5)
plt.legend()
plt.xlabel('Time since first star formation (Myr)')
plt.ylabel('Velocity spread (km/s)')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('v_spread_vs_time_'+simulation_density_id+'.png')
