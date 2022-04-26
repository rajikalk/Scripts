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

try:
    file = open("sink_birth_conditions.pkl", 'wb')
    Sink_bound_birth = pickle.load(file)
    file.close()
except:
    Sink_bound_birth = []
Mass_plus_blank_row = np.vstack([np.zeros(len(global_data['m'][0])), global_data['m']])
diff_arr =  (Mass_plus_blank_row[1:]-Mass_plus_blank_row[:-1])
zero_inds = np.where(diff_arr == 0)
diff_arr[zero_inds] = 1
formation_inds = np.where(diff_arr == global_data['m'])
if len(Sink_bound_birth) > 0:
    import pdb
    pdb.set_trace()
n_stars = 0
for sink_id in formation_inds[1]:
    form_time_it = formation_inds[0][sink_id]
    new_sink_pos = np.array([global_data['x'][form_time_it][sink_id], global_data['y'][form_time_it][sink_id], global_data['z'][form_time_it][sink_id]]).T
    new_sink_vel = np.array([global_data['ux'][form_time_it][sink_id], global_data['uy'][form_time_it][sink_id], global_data['uz'][form_time_it][sink_id]]).T
    new_sink_mass = np.array(global_data['m'][form_time_it][sink_id])

    abspos = np.array([global_data['x'][form_time_it][:sink_id], global_data['y'][form_time_it][:sink_id], global_data['z'][form_time_it][:sink_id]]).T
    absvel = np.array([global_data['ux'][form_time_it][:sink_id], global_data['uy'][form_time_it][:sink_id], global_data['uz'][form_time_it][:sink_id]]).T
    mass = np.array(global_data['m'][form_time_it][:sink_id])
    
    rel_pos = abspos - new_sink_pos
    update_seps_neg = np.argwhere(rel_pos<-0.5)
    update_seps_pos = np.argwhere(rel_pos>0.5)
    
    rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] = rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] + 0.5
    rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] = rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] - 0.5

    rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
    rel_vel = absvel - new_sink_vel
    rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
    mtm = new_sink_mass * mass
    mpm = new_sink_mass + mass
    newtonianPotential = -1./rel_sep
    
    Ekin = 0.5 * mtm/mpm * rel_speed**2
    Epot = Grho * mtm * newtonianPotential
    Etot = Ekin + Epot
    
    sep_below_10000 = np.where((units['length_unit'].in_units('au')*rel_sep)<10000)[0]
    
    if True in (Etot[sep_below_10000]<0):
        born_bound = True
        most_bound_sink_id = np.argmin(Etot)
        lowest_Etot = np.nanmin(Etot)
        #Do multiplicity analysis
        time_it = formation_inds[0][sink_id]
        n_stars = np.where(global_data['m'][time_it]>0)[0]

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
        if sink_id in res['index1']:
            sys_id = np.argwhere(res['index1'] == sink_id)[0][0]
            first_bound_sink = res['index2'][sys_id]
        elif sink_id in res['index2']:
            sys_id = np.argwhere(res['index2'] == sink_id)[0][0]
            first_bound_sink = res['index1'][sys_id]
        else:
            sys_id = np.nan
        if np.isnan(sys_id) == False:
            first_bound_sink = losi(first_bound_sink, res)
            lowest_Etot = res['epot'][sys_id] + res['ekin'][sys_id]
            most_bound_sep = res['separation'][sys_id]
            break
    else:
        born_bound = False
        if len(Etot) > 0:
            most_bound_sink_id = np.argmin(Etot)
        else:
            most_bound_sink_id = np.nan
        most_bound_sep = np.nan
        first_bound_sink = np.nan
        lowest_Etot = np.nan
        
        time_it = formation_inds[0][sink_id]

        new_sink_pos_x = global_data['x'][time_it:,sink_id]
        new_sink_pos_y = global_data['y'][time_it:,sink_id]
        new_sink_pos_z = global_data['z'][time_it:,sink_id]
        
        new_sink_vel_x = global_data['ux'][time_it:,sink_id]
        new_sink_vel_y = global_data['uy'][time_it:,sink_id]
        new_sink_vel_z = global_data['uz'][time_it:,sink_id]
        
        new_sink_mass = np.array(global_data['m'][time_it:,sink_id])
        
        abspos_x = global_data['x'][time_it:]
        abspos_y = global_data['y'][time_it:]
        abspos_z = global_data['z'][time_it:]
        
        absvel_x = global_data['ux'][time_it:]
        absvel_y = global_data['uy'][time_it:]
        absvel_z = global_data['uz'][time_it:]
        
        mass = np.array(global_data['m'][time_it:])
        
        rel_pos_x = abspos_x.T - new_sink_pos_x
        rel_pos_y = abspos_y.T - new_sink_pos_y
        rel_pos_z = abspos_z.T - new_sink_pos_z
        
        update_seps_x_neg = np.argwhere(rel_pos_x<-0.5)
        update_seps_x_pos = np.argwhere(rel_pos_x>0.5)
        rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] = rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] + 0.5
        rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] = rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] - 0.5
        
        update_seps_y_neg = np.argwhere(rel_pos_y<-0.5)
        update_seps_y_pos = np.argwhere(rel_pos_y>0.5)
        rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] = rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] + 0.5
        rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] = rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] - 0.5
        
        update_seps_z_neg = np.argwhere(rel_pos_z<-0.5)
        update_seps_z_pos = np.argwhere(rel_pos_z>0.5)
        rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] = rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] + 0.5
        rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] = rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] - 0.5
        
        rel_sep = np.sqrt(rel_pos_x**2 + rel_pos_y**2 + rel_pos_z**2)
        
        rel_vel_x = absvel_x.T - new_sink_vel_x
        rel_vel_y = absvel_y.T - new_sink_vel_y
        rel_vel_z = absvel_z.T - new_sink_vel_z
        rel_speed = np.sqrt(rel_vel_x**2 + rel_vel_y**2 + rel_vel_z**2)
        
        mtm = new_sink_mass * mass.T
        mpm = new_sink_mass + mass.T
        newtonianPotential = -1./rel_sep
        Ekin = 0.5 * mtm/mpm * rel_speed**2
        Epot = Grho * mtm * newtonianPotential
        Etot = Ekin + Epot
        Etot[Etot == -1*np.inf] = 0
        
        Etot_min = np.min(Etot, axis=0)
        Etot_min_sink_id = np.argmin(Etot, axis=0)
        Etot_bound_inds = np.where(Etot_min<0)[0]
        
        rel_sep[np.where(rel_sep == 0)] = np.inf
        closest_separations = np.min(rel_sep, axis=0)
        closest_sink_id = np.argmin(rel_sep, axis=0)
        sep_below_10000 = np.where((units['length_unit'].in_units('au')*closest_separations)<10000)[0]
        
        test_time_inds = list(set(Etot_bound_inds).intersection(set(sep_below_10000)))
        
        for test_time_ind in test_time_inds:
            if np.isnan(first_bound_sink):
                time_it = formation_inds[0][sink_id] + test_time_ind
                print("testing time_it", time_it)
                n_stars = np.where(global_data['m'][time_it]>0)[0]
                
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
                if sink_id in res['index1']:
                    sys_id = np.argwhere(res['index1'] == sink_id)[0][0]
                    first_bound_sink = res['index2'][sys_id]
                elif sink_id in res['index2']:
                    sys_id = np.argwhere(res['index2'] == sink_id)[0][0]
                    first_bound_sink = res['index1'][sys_id]
                else:
                    sys_id = np.nan
                if np.isnan(sys_id) == False:
                    first_bound_sink = losi(first_bound_sink, res)
                    lowest_Etot = res['epot'][sys_id] + res['ekin'][sys_id]
                    most_bound_sep = res['separation'][sys_id]
                    break

    if np.isnan(most_bound_sep) and lowest_Etot < 0:
        import pdb
        pdb.set_trace()
    Sink_bound_birth.append([born_bound, most_bound_sink_id, first_bound_sink, most_bound_sep, lowest_Etot])
    print("Birth conditions of sink", sink_id, "is", Sink_bound_birth[-1])

    file = open("sink_birth_conditions.pkl", 'wb')
    pickle.dump((Sink_bound_birth),file)
    file.close()
