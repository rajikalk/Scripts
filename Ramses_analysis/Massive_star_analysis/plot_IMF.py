import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
import sys
import pickle
import os

global_pickles = ['/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G50/stars_imf_G50.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G100/256/stars_imf_G100.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G125/stars_imf_G125.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G150/stars_imf_G150.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G200/stars_imf_G200.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G400/stars_imf_G400.pkl']

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
plt.clf()
for global_pickle in global_pickles:
    simulation_density_id = global_pickle.split('/G')[-1].split('/')[0]

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
        

    units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
        
    scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
    scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
    scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
    scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

    units={}
    for key in units_override.keys():
        units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

    #====================================================================================================================================================
    file_open = open(global_pickle, 'rb')
    try:
        global_data = pickle.load(file_open,encoding="latin1")
    except:
        file_open.close()
        import pickle5 as pickle
        file_open = open(global_pickle, 'rb')
        global_data = pickle.load(file_open,encoding="latin1")
    file_open.close()

    SFE_5_time_ind = np.argmin(abs(np.sum(global_data['m'], axis=1)-0.05))
    masses_SFE_5 = global_data['m'][SFE_5_time_ind]*units['mass_unit']

    masses_end = global_data['m'][-1]*units['mass_unit']

    bins = np.logspace(-2,2,13)
    IMF_SFE_5, bins =np.histogram(masses_SFE_5, bins=bins)
    IMF_end, bins =np.histogram(masses_end, bins=bins)
    bins_centers = (bins[1:] + bins[:-1])/2

    plt.semilogx(bins_centers, IMF_end, label='G'+str(simulation_density_id))
plt.xlabel('Mass (M$_\odot$)')
plt.ylabel('Number')
plt.legend()
plt.ylim(bottom=0)
plt.savefig('IMF.png')

