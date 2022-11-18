import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
from mpi4py.MPI import COMM_WORLD as CW
import sys
import collections
import os
import pickle
import gc

f_acc= 0.5

def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def accretion(sink_inds, time_ind):
    """
    Calculates the accretion of the given indeexes
    """
    global Accretion_array
    M_dot = Accretion_array[time_ind, sink_inds]
    return M_dot
    
def luminosity(global_data, sink_inds, global_ind):
    """
    Calculates the luminosity of the given indexes
    """
    global f_acc
    radius = yt.YTQuantity(2.0, 'rsun')
    M_dot = accretion(sink_inds, global_ind)
    M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    return L_tot

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-superplot", "--superplot_pickle", help="what pickle will I reade the multiplicity analysis from?", type=str)

    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()


#Set units
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/G')[-1].split('.')[0]

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

sys.stdout.flush()
CW.Barrier()

#loading global data
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
print('loaded global data')

sys.stdout.flush()
CW.Barrier()

#Iterate through systems from superplot pickle, and retrieve masses from global data. I guess I should have lifetime limits?

#load superplot data
file = open(args.superplot_pickle, 'rb')
superplot_dict, Sink_bound_birth, Sink_formation_times = pickle.load(file)
file.close()
candidate_systems = []
final_masses = []

for sys_key in superplot_dict['System_times'].keys():
    if len(flatten(eval(sys_key))) > 2:
        sep_arr = np.array(superplot_dict['System_seps'][sys_key]).T[0]
        non_nan_inds = np.where(np.isnan(sep_arr)==False)[0]
        last_sys_time = np.array(superplot_dict['System_times'][sys_key])[non_nan_inds][-1]
        t_ind = np.argmin(abs((global_data['time']*units['time_unit'].in_units('yr').value) - last_sys_time))
        #lifetime = last_sys_time - np.array(superplot_dict['System_times'][sys_key])[non_nan_inds][0]
        #if lifetime > 10000:
        masses = global_data['m'][t_ind][flatten(eval(sys_key))]*units['mass_unit'].in_units('msun')
        if len(flatten(eval(sys_key))) == 3:
            inner_mass_max = np.max(masses)
        else:
            #check if inner trinary exists
            #If n is odd, there is a trinary, but quadruples could be [[],[]]
            #need to searhc for [[,],] or [,[,]]
            structure = ''.join([i for i in sys_key if not i.isdigit()])
            if '[, [, ]]' in structure or '[[, ], ]' in structure:
                #strip down to inner triple
                if '[, [, ]]' in structure and '[[, ], ]' in structure:
                    import pdb
                    pdb.set_trace()
                    
                if '[, [, ]]' in structure:
                    rm_bracket_structure = structure.split('[, [, ]]')
                    fitted_struct = '[, [, ]]'
                else:
                    rm_bracket_structure = structure.split('[[, ], ]')
                    fitted_struct = '[[, ], ]'

                if len(rm_bracket_structure) == 2:
                    stripped_string = sys_key
                    for char_it in range(len(stripped_string)):
                        sub_struct = ''.join([i for i in stripped_string[:char_it] if not i.isdigit()])
                        if sub_struct == rm_bracket_structure[0]:
                            trun_ind = char_it
                            break
                    stripped_string = stripped_string[trun_ind:]
                        #trim from end
                    for char_it in range(len(stripped_string)):
                        sub_struct = ''.join([i for i in stripped_string[-1*char_it:] if not i.isdigit()])
                        if sub_struct == rm_bracket_structure[1]:
                            trun_ind = char_it
                            break
                    stripped_string = stripped_string[:-1*char_it]
                    
                    inner_mass_max = global_data['m'][t_ind][flatten(eval(stripped_string))]*units['mass_unit'].in_units('msun')
                elif len(rm_bracket_structure)>2:
                    #break down in the two trinaries
                    #Find midpoint to split
                    stripped_string = sys_key
                    for char_it in range(len(stripped_string)):
                        sub_struct = ''.join([i for i in stripped_string[:char_it] if not i.isdigit()])
                        if sub_struct == rm_bracket_structure[0] + fitted_struct + rm_bracket_structure[1]:
                            trun_ind = char_it
                            break
                    first_sys = stripped_string[:trun_ind][1:-2]
                    second_sys = stripped_string[trun_ind:][:-1]
                    mass_first = global_data['m'][t_ind][flatten(eval(first_sys))]*units['mass_unit'].in_units('msun')
                    mass_second = global_data['m'][t_ind][flatten(eval(second_sys))]*units['mass_unit'].in_units('msun')
                    if np.max(mass_first) > np.max(mass_second):
                        inner_mass_max = mass_first
                    else:
                        inner_mass_max = mass_second
                        
        if np.max(inner_mass_max) > 5:
            candidate_systems.append(sys_key)
            final_masses.append(masses)
            
            if np.max(inner_mass_max) > 8:
                savename = '8msun_candidate_'+sys_key.replace(' ', '')+'.png'
                print('This is one over 8msun!')
            else:
                savename = 'candidate_'+sys_key.replace(' ', '')+'.png'
            
            plt.clf()
            plt.semilogy(superplot_dict['System_times'][sys_key], superplot_dict['System_seps'][sys_key])
            plt.xlabel('time (yr)')
            plt.ylabel('separation (au)')
            plt.title('system = '+str(stripped_string) +', final mass = '+ str(inner_mass_max.value))
            plt.savefig(savename)
            print('plotted a candidate')
        else:
            print("found a multiple, but it's low mass")
