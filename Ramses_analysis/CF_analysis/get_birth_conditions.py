import numpy as np
import yt
import pickle
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import collections
from mpi4py.MPI import COMM_WORLD as CW
import sys

rank = CW.Get_rank()
size = CW.Get_size()

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
# 136578 time inds for G50
del file_open

sys.stdout.flush()
CW.Barrier()

try:
    file = open("sink_birth_conditions.pkl", 'rb')
    Sink_bound_birth = pickle.load(file)
    file.close()
except:
    Sink_bound_birth = []
Mass_plus_blank_row = np.vstack([np.zeros(len(global_data['m'][0])), global_data['m']])
diff_arr =  (Mass_plus_blank_row[1:]-Mass_plus_blank_row[:-1])
del Mass_plus_blank_row
zero_inds = np.where(diff_arr == 0)
diff_arr[zero_inds] = 1
del zero_inds
formation_inds = np.where(diff_arr == global_data['m'])
del diff_arr
if len(Sink_bound_birth) > 0:
    sink_id = len(Sink_bound_birth)
else:
    sink_id = 0

sys.stdout.flush()
CW.Barrier()

file = open("All_sink_birth_conditions.pkl", 'rb')
Sink_birth_all = pickle.load(file)
file.close()

mismatched_inds = [13, 15, 16, 20, 24, 28, 31, 40, 46, 50, 63, 65]
true_delay = []
true_first_sys = []

for birth_con in Sink_birth_all:
    if birth_con[0] in mismatched_inds:
        true_delay.append(birth_con[-1])
        true_first_sys.append(birth_con[3])

#true_delay = [700973.119076509, 825781.5734443031, 0.0, 656453.2060968727, 0.0, 0.0, 856.0596429631114, 13638.631729342043, 404753.4436713271, 6036.782605849206, 70337.8161722906, 79088.95493660122]

sys.stdout.flush()
CW.Barrier()

rit = -1
while sink_id < len(formation_inds[1]):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        form_time_it = formation_inds[0][sink_id]

        new_sink_pos = np.array([global_data['x'][form_time_it][sink_id], global_data['y'][form_time_it][sink_id], global_data['z'][form_time_it][sink_id]]).T
        abspos = np.array([global_data['x'][form_time_it][:sink_id], global_data['y'][form_time_it][:sink_id], global_data['z'][form_time_it][:sink_id]]).T
        rel_pos = abspos - new_sink_pos
        del new_sink_pos
        del abspos
        update_seps_neg = np.argwhere(rel_pos<-0.5)
        update_seps_pos = np.argwhere(rel_pos>0.5)
        
        rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] = rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] + 1.0
        del update_seps_neg
        rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] = rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] - 1.0
        del update_seps_pos

        rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
        del rel_pos
        
        new_sink_vel = np.array([global_data['ux'][form_time_it][sink_id], global_data['uy'][form_time_it][sink_id], global_data['uz'][form_time_it][sink_id]]).T
        absvel = np.array([global_data['ux'][form_time_it][:sink_id], global_data['uy'][form_time_it][:sink_id], global_data['uz'][form_time_it][:sink_id]]).T
        rel_vel = absvel - new_sink_vel
        del new_sink_vel
        del absvel
        rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
        del rel_vel
        
        new_sink_mass = np.array(global_data['m'][form_time_it][sink_id])
        mass = np.array(global_data['m'][form_time_it][:sink_id])
        mtm = new_sink_mass * mass
        mpm = new_sink_mass + mass
        del new_sink_mass
        del mass
        
        newtonianPotential = -1./rel_sep
        
        Ekin = 0.5 * mtm/mpm * rel_speed**2
        del rel_speed
        del mpm
        Epot = Grho * mtm * newtonianPotential
        del newtonianPotential
        del mtm
        Etot = Ekin + Epot
        try:
            most_bound_sink_id = np.argmin(Etot)
        except:
            most_bound_sink_id = np.nan
        del Ekin
        del Epot
        
        sep_below_10000 = np.where((units['length_unit'].in_units('au')*rel_sep)<10000)[0]
        del rel_sep
        
        sys_id = np.nan
        if True in (Etot[sep_below_10000]<0):
            #del sep_below_10000
            born_bound = True
            lowest_Etot = np.nanmin(Etot)
            delay_time = 0
            #Do multiplicity analysis
            time_it = formation_inds[0][sink_id]
            n_stars = np.where(global_data['m'][time_it]>0)[0]

            abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][time_it][n_stars])
            time = global_data['time'][time_it][n_stars][0]
            del n_stars
            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l.value
            S._scale_v = scale_v.value
            S._scale_t = scale_t.value
            S._scale_d = scale_d.value
            S._time = yt.YTArray(time, '')
            del time
            S._abspos = yt.YTArray(abspos, '')
            del abspos
            S._absvel = yt.YTArray(absvel, '')
            del absvel
            S._mass = yt.YTArray(mass, '')
            del mass
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
            del res
        #if True not in (Etot[sep_below_10000]<0) or np.isnan(sys_id):
        if np.isnan(sys_id):
            born_bound = False
            most_bound_sep = np.nan
            first_bound_sink = np.nan
            lowest_Etot = np.nan
            delay_time = np.nan
            
            time_it = formation_inds[0][sink_id]
            formation_time = global_data['time'][time_it][0]*units['time_unit'].in_units('yr')
            #test_time_inds = range(len(global_data['x'][time_it:,sink_id]))
            
            new_sink_pos_x = global_data['x'][time_it:,sink_id]
            new_sink_pos_y = global_data['y'][time_it:,sink_id]
            new_sink_pos_z = global_data['z'][time_it:,sink_id]
            
            abspos_x = global_data['x'][time_it:]
            abspos_y = global_data['y'][time_it:]
            abspos_z = global_data['z'][time_it:]
            
            rel_pos_x = abspos_x.T - new_sink_pos_x
            rel_pos_y = abspos_y.T - new_sink_pos_y
            rel_pos_z = abspos_z.T - new_sink_pos_z
            del new_sink_pos_x
            del new_sink_pos_y
            del new_sink_pos_z
            del abspos_x
            del abspos_y
            del abspos_z
            
            update_seps_x_neg = np.argwhere(rel_pos_x<-0.5)
            update_seps_x_pos = np.argwhere(rel_pos_x>0.5)
            rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] = rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] + 1.0
            rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] = rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] - 1.0
            del update_seps_x_neg
            del update_seps_x_pos
            
            update_seps_y_neg = np.argwhere(rel_pos_y<-0.5)
            update_seps_y_pos = np.argwhere(rel_pos_y>0.5)
            rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] = rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] + 1.0
            rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] = rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] - 1.0
            del update_seps_y_neg
            del update_seps_y_pos
            
            update_seps_z_neg = np.argwhere(rel_pos_z<-0.5)
            update_seps_z_pos = np.argwhere(rel_pos_z>0.5)
            rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] = rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] + 1.0
            rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] = rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] - 1.0
            del update_seps_z_neg
            del update_seps_z_pos
            
            rel_sep = np.sqrt(rel_pos_x**2 + rel_pos_y**2 + rel_pos_z**2)
            del rel_pos_x
            del rel_pos_y
            del rel_pos_z
            '''
            new_sink_vel_x = global_data['ux'][time_it:,sink_id]
            new_sink_vel_y = global_data['uy'][time_it:,sink_id]
            new_sink_vel_z = global_data['uz'][time_it:,sink_id]
            
            absvel_x = global_data['ux'][time_it:]
            absvel_y = global_data['uy'][time_it:]
            absvel_z = global_data['uz'][time_it:]
            
            rel_vel_x = absvel_x.T - new_sink_vel_x
            rel_vel_y = absvel_y.T - new_sink_vel_y
            rel_vel_z = absvel_z.T - new_sink_vel_z
            del new_sink_vel_x
            del new_sink_vel_y
            del new_sink_vel_z
            del absvel_x
            del absvel_y
            del absvel_z
            
            rel_speed = np.sqrt(rel_vel_x**2 + rel_vel_y**2 + rel_vel_z**2)
            del rel_vel_x
            del rel_vel_y
            del rel_vel_z
            
            new_sink_mass = np.array(global_data['m'][time_it:,sink_id])
            mass = np.array(global_data['m'][time_it:])
            
            mtm = new_sink_mass * mass.T
            mpm = new_sink_mass + mass.T
            del new_sink_mass
            del mass
            newtonianPotential = -1./rel_sep
            Ekin = 0.5 * mtm/mpm * rel_speed**2
            del mpm
            Epot = Grho * mtm * newtonianPotential
            del mtm
            del newtonianPotential
            Etot = Ekin + Epot
            del Ekin
            del Epot
            Etot[Etot == -1*np.inf] = 0
            
            Etot_min = np.min(Etot, axis=0)
            Etot_min_sink_id = np.argmin(Etot, axis=0)
            Etot_bound_inds = np.where(Etot_min<0)[0]
            
            #del Etot
            #del Etot_min
            '''
            rel_sep[np.where(rel_sep == 0)] = np.inf
            closest_separations = np.min(rel_sep, axis=0)
            closest_sink_id = np.argmin(rel_sep, axis=0)
            #del rel_sep
            #import pdb
            #pdb.set_trace()
            sep_below_10000 = np.where((units['length_unit'].in_units('au')*closest_separations)<10000)[0]
            test_time_inds = sep_below_10000
            #del closest_separations
            
            if sink_id in mismatched_inds:
                goal_delay = true_delay[mismatched_inds.index(sink_id)]
                goal_time = formation_time.value + goal_delay
                all_times = global_data['time'][time_it:].T[0]*units['time_unit'].in_units('yr')
                goal_ind = np.argmin(abs(all_times.value - goal_time))
            '''
            if True not in (Etot_min_sink_id == closest_sink_id):
                test_time_inds = []
            else:
                test_time_inds = sorted(list(set(Etot_bound_inds).intersection(set(sep_below_10000))))
            
            del closest_sink_id
            del Etot_bound_inds
            del sep_below_10000
            '''
            for test_time_ind in test_time_inds:
                if np.isnan(first_bound_sink):
                    time_it = formation_inds[0][sink_id] + test_time_ind
                    print("testing time_it", time_it, "on rank", rank)
                    n_stars = np.where(global_data['m'][time_it]>0)[0]
                    if len(n_stars)>1:
                        abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                        absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                        mass = np.array(global_data['m'][time_it][n_stars])
                        time = global_data['time'][time_it][n_stars][0]
                        del n_stars
                        S = pr.Sink()
                        S._jet_factor = 1.
                        S._scale_l = scale_l.value
                        S._scale_v = scale_v.value
                        S._scale_t = scale_t.value
                        S._scale_d = scale_d.value
                        S._time = yt.YTArray(time, '')
                        del time
                        S._abspos = yt.YTArray(abspos, '')
                        del abspos
                        S._absvel = yt.YTArray(absvel, '')
                        del absvel
                        S._mass = yt.YTArray(mass, '')
                        del mass
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
                            bound_time = global_data['time'][time_it][0]*units['time_unit'].in_units('yr')
                            delay_time = float((bound_time - formation_time).value)
                            try:
                                if Sink_birth_all[np.argwhere(Sink_birth_all[:,0]==sink_id)][0][0][4] != most_bound_sep:
                                    import pdb
                                    pdb.set_trace()
                            except:
                                pass
                            break
                        del res

        if np.isnan(most_bound_sep) and lowest_Etot < 0:
            import pdb
            pdb.set_trace()
        Sink_bound_birth.append([sink_id, born_bound, most_bound_sink_id, str(first_bound_sink), most_bound_sep, lowest_Etot, delay_time])
        print("Birth conditions of sink", sink_id, "is", Sink_bound_birth[-1])

        file = open("sink_birth_conditions_"+str(rank)+".pkl", 'wb')
        pickle.dump((Sink_bound_birth),file)
        file.close()
        
    sink_id = sink_id + 1

sys.stdout.flush()
CW.Barrier()

#compile pickles
if rank == 0:
    import glob
    birth_pickles = sorted(glob.glob("sink_birth_conditions_*.pkl"))
    Sink_birth_all = []
    for birth_pick in birth_pickles:
        file = open(birth_pick, 'rb')
        Sink_bound_birth_rank = pickle.load(file)
        file.close()
        Sink_birth_all = Sink_birth_all + Sink_bound_birth_rank
        os.remove(birth_pick)
    
    sorted_inds = np.argsort(list(map(int, np.array(Sink_birth_all)[:,0])))
    Sink_birth_all = np.array(Sink_birth_all)[sorted_inds]
        
    file = open("sink_birth_all.pkl", 'wb')
    pickle.dump((Sink_birth_all), file)
    file.close()
    print("Collected sink birth data into sink_birth_all.pkl" )

'''
for sink_id in range(np.shape(Sink_birth_fast)[0]):
    if sink_id in Sink_birth_all[:,0]:
        match_ind = np.argwhere(Sink_birth_all[:,0]==sink_id)[0]
        Sink_match = Sink_birth_all[match_ind]
        Sink_fast = Sink_birth_fast[sink_id]
        try:
            if False not in (Sink_match[0][1:] == np.array(Sink_fast)):
                print('Birth conditions match for sink', sink_id)
            else:
                False_inds = np.argwhere((Sink_match[0][1:] == np.array(Sink_fast))==False).T[0]
                if np.sum(np.nan_to_num(Sink_match[0][1:][False_inds].tolist())) == 0 and np.sum(np.nan_to_num(np.array(Sink_fast)[False_inds].tolist())) == 0:
                    print("it's all nan")
                else:
                    mismatched_inds.append(sink_id)
                    true_delay.append(Sink_match[-1])
                    true_first_sys.append(Sink_bound_birth[3])
        except:
            mismatched_inds.append(sink_id)
            true_delay.append(Sink_match[-1])
'''
"""
import pickle
import numpy as np
file_open = open("/groups/astro/rlk/rlk/Global_sink_pickles/G100_full.pkl", "rb")
global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
SFE_5_ind = np.argmin(abs(np.sum(global_data['m'],axis=1)-0.05)) + 1
max_sink = np.where(global_data['m'][SFE_5_ind]>0)[0][-1]+1
file_open = open("global_reduced.pkl", "wb")
pickle.dump(({'time':global_data['time'][:SFE_5_ind].T[0],'m': global_data['m'][:SFE_5_ind,:max_sink], 'x': global_data['x'][:SFE_5_ind,:max_sink], 'y': global_data['y'][:SFE_5_ind,:max_sink], 'z': global_data['z'][:SFE_5_ind,:max_sink], 'ux': global_data['ux'][:SFE_5_ind,:max_sink], 'uy': global_data['uy'][:SFE_5_ind,:max_sink], 'uz': global_data['uz'][:SFE_5_ind,:max_sink]}), file_open)
file_open.close()
