import numpy as np
import yt
import pickle
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import collections

#Define globals
f_acc= 0.5
Accretion_array = []

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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-verbose", "--verbose_printing", help="Would you like to print debug lines?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    

#=====================================================================================================

args = parse_inputs()

#==========================================================================================

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('_')[0][1:]

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
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
    
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
dm = global_data['dm']*units['mass_unit'].in_units('Msun')
dt = (global_data['time'] - global_data['tflush'])*units['time_unit'].in_units('yr')

Sink_bound_birth = []
Mass_plus_blank_row = np.vstack([np.zeros(len(global_data['m'][0])), global_data['m']])
diff_arr =  (Mass_plus_blank_row[1:]-Mass_plus_blank_row[:-1])
zero_inds = np.where(diff_arr == 0)
diff_arr[zero_inds] = 1
formation_inds = np.where(diff_arr == global_data['m'])
n_stars = 0
for sink_id in formation_inds[1]:
    new_sink_pos = np.array([global_data['x'][formation_inds[0][sink_id]][sink_id], global_data['y'][formation_inds[0][sink_id]][sink_id], global_data['z'][formation_inds[0][sink_id]][sink_id]]).T
    new_sink_vel = np.array([global_data['ux'][formation_inds[0][sink_id]][sink_id], global_data['uy'][formation_inds[0][sink_id]][sink_id], global_data['uz'][formation_inds[0][sink_id]][sink_id]]).T
    new_sink_mass = np.array(global_data['m'][formation_inds[0][sink_id]][sink_id])

    abspos = np.array([global_data['x'][formation_inds[0][sink_id]][:sink_id], global_data['y'][formation_inds[0][sink_id]][:sink_id], global_data['z'][formation_inds[0][sink_id]][:sink_id]]).T
    absvel = np.array([global_data['ux'][formation_inds[0][sink_id]][:sink_id], global_data['uy'][formation_inds[0][sink_id]][:sink_id], global_data['uz'][formation_inds[0][sink_id]][:sink_id]]).T
    mass = np.array(global_data['m'][formation_inds[0][sink_id]][:sink_id])
    
    rel_pos = abspos - new_sink_pos
    update_seps = np.argwhere(abs(rel_pos)>0.5)
    for update_sep in update_seps:
        if rel_pos[update_sep[0]][update_sep[1]] < 0:
            rel_pos[update_sep[0]][update_sep[1]] = rel_pos[update_sep[0]][update_sep[1]] + 0.5
        else:
            rel_pos[update_sep[0]][update_sep[1]] = rel_pos[update_sep[0]][update_sep[1]] - 0.5
    rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
    rel_vel = absvel - new_sink_vel
    rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
    mtm = new_sink_mass * mass
    mpm = new_sink_mass + mass
    reducedMass = mtm/mpm
    newtonianPotential = -1./rel_sep
    
    Ekin = 0.5 * mtm/mpm * rel_speed**2
    Epot = Grho * mtm * newtonianPotential
    Etot = Ekin + Epot
    
    #Do multiplicity analysis
    time_it = formation_inds[0][sink_id]
    n_stars = np.where(global_data['m'][time_it]>0)[0]
    if len(n_stars) > 1:
        abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
        absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
        mass = np.array(global_data['m'][time_it][n_stars])
        time = global_data['time'][time_it][n_stars][0]
        S = pr.Sink()
        S._jet_factor = 1.
        S._scale_l = scale_l.value
        S._scale_v = scale_v.value
        S._scale_t = scale_t.value
        S._scale_d = scale_d.value
        S._time = yt.YTArray(time, '')
        S._abspos = yt.YTArray(abspos, '')
        S._absvel = yt.YTArray(absvel, '')
        S._mass = yt.YTArray(mass, '')
        res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho)
        multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
        for multi_ind in multi_inds:
            sys_comps = losi(multi_ind, res)
            sys_string = sorted(flatten(sys_comps))
            if sink_id in sys_string:
                import pdb
                pdb.set_trace()
    
    if True in (Etot<0):
        born_bound = True
        most_bound_sink_id = np.argmin(Etot)
        most_bound_sep = rel_sep[most_bound_sink_id]*scale_l.in_units('au')
    else:
        born_bound = False
        if len(Etot) > 0:
            most_bound_sink_id = np.argmin(Etot)
        else:
            most_bound_sink_id = np.nan
        most_bound_sep = np.nan
    Sink_bound_birth.append([born_bound, most_bound_sink_id, most_bound_sep])
    print("Found birth conditions of sink", sink_id)

file = open("sink_birth_conditions.pkl", 'wb')
pickle.dump((Sink_bound_birth),file)
file.close()
