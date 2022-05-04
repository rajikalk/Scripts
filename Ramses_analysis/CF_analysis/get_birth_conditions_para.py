import numpy as np
import yt
import pickle
import pyramses as pr
from pyramses import rsink
import multiplicity as m
from mpi4py.MPI import COMM_WORLD as CW
import sys

rank = CW.Get_rank()
size = CW.Get_size()

def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

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

#simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('_')[0][1:]
simulation_density_id = args.global_data_pickle_file.split('/G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50.
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100.
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125.
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150.
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200.
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400.
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
del scale_v

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
    
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open)
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open)
file_open.close()
# 136578 time inds for G50
del file_open

sys.stdout.flush()
CW.Barrier()

formation_inds = [0]
for sink_id in range(1, np.shape(global_data['m'].T)[0]):
    formation_inds.append(np.argwhere(global_data['m'].T[sink_id]>0)[0][0])

formation_inds = np.array(formation_inds)
formation_times = global_data['time'][formation_inds]

sys.stdout.flush()
CW.Barrier()
Sink_bound_birth = []

#Testing accuracy:
try:
    file = open("sink_birth_all_true.pkl", 'rb')
    True_sink_birth_conditions = pickle.load(file)
    file.close()
    mismatched_inds = []
except:
    print("True birth conditions don't exist")


rit = -1
sink_id = 16 #0
while sink_id < len(formation_inds):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        form_time_it = np.where(global_data['time']==formation_times[sink_id])[0][0]# formation_inds[0][sink_id]
        
        #truncate global data
        global_data['time'] = global_data['time'][form_time_it:]
        global_data['m'] = global_data['m'][form_time_it:]
        global_data['x'] = global_data['x'][form_time_it:]
        global_data['y'] = global_data['y'][form_time_it:]
        global_data['z'] = global_data['z'][form_time_it:]
        global_data['ux'] = global_data['ux'][form_time_it:]
        global_data['uy'] = global_data['uy'][form_time_it:]
        global_data['uz'] = global_data['uz'][form_time_it:]
        
        new_sink_pos = np.array([global_data['x'][0][sink_id], global_data['y'][0][sink_id], global_data['z'][0][sink_id]]).T
        abspos = np.array([global_data['x'][0][:sink_id], global_data['y'][0][:sink_id], global_data['z'][0][:sink_id]]).T
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
        
        new_sink_vel = np.array([global_data['ux'][0][sink_id], global_data['uy'][0][sink_id], global_data['uz'][0][sink_id]]).T
        absvel = np.array([global_data['ux'][0][:sink_id], global_data['uy'][0][:sink_id], global_data['uz'][0][:sink_id]]).T
        rel_vel = absvel - new_sink_vel
        del new_sink_vel
        del absvel
        rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
        del rel_vel
        
        new_sink_mass = np.array(global_data['m'][0][sink_id])
        mass = np.array(global_data['m'][0][:sink_id])
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
        
        born_bound = True
        delay_time = 0
        
        n_stars = np.where(global_data['m'][0]>0)[0]
        sys_id = np.nan
        if len(n_stars)>1:
            abspos = np.array([global_data['x'][0][n_stars], global_data['y'][0][n_stars], global_data['z'][0][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][0][n_stars], global_data['uy'][0][n_stars], global_data['uz'][0][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][0][n_stars])
            del n_stars
            time = global_data['time'][0]
            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l.value
            #S._scale_v = scale_v.value
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
            #Find most bound sink!
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
        '''
        form_time_it = 0

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
            delay_time = 0
            #Do multiplicity analysis
            time_it = 0
            n_stars = np.where(global_data['m'][time_it]>0)[0]

            abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][time_it][n_stars])
            time = global_data['time'][time_it]
            del n_stars
            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l.value
            #S._scale_v = scale_v.value
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
        '''
        if np.isnan(sys_id):
            born_bound = False
            most_bound_sep = np.nan
            first_bound_sink = np.nan
            lowest_Etot = np.nan
            delay_time = np.nan
            
            formation_time = formation_times[sink_id]*units['time_unit'].in_units('yr')
            #test_time_inds = range(len(global_data['x'][time_it:,sink_id]))
            
            new_sink_pos_x = global_data['x'][:,sink_id]
            new_sink_pos_y = global_data['y'][:,sink_id]
            new_sink_pos_z = global_data['z'][:,sink_id]
            
            abspos_x = global_data['x'][:]
            abspos_y = global_data['y'][:]
            abspos_z = global_data['z'][:]
            
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
            del rel_sep
            #import pdb
            #pdb.set_trace()
            test_time_inds = np.where((units['length_unit'].in_units('au')*closest_separations)<10000)[0]
            del closest_separations
            
            '''
            if sink_id in mismatched_inds:
                goal_delay = true_delay[mismatched_inds.index(sink_id)]
                goal_time = formation_time.value + goal_delay
                all_times = global_data['time'][time_it:].T[0]*units['time_unit'].in_units('yr')
                goal_ind = np.argmin(abs(all_times.value - goal_time))
                
            if True not in (Etot_min_sink_id == closest_sink_id):
                test_time_inds = []
            else:
                test_time_inds = sorted(list(set(Etot_bound_inds).intersection(set(sep_below_10000))))
            
            del closest_sink_id
            del Etot_bound_inds
            del sep_below_10000
            '''
            #test_time_inds = np.arange(np.shape(global_data['time'])[0])
            
            for test_time_ind in test_time_inds:
                if np.isnan(first_bound_sink):
                    time_it = 0 + test_time_ind
                    if np.remainder(time_it,5000) == 0:
                        print("testing time_it", time_it, "on rank", rank)
                    n_stars = np.where(global_data['m'][time_it]>0)[0]
                    if len(n_stars)>1:
                        abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                        absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                        mass = np.array(global_data['m'][time_it][n_stars])
                        time = global_data['time'][time_it]
                        del n_stars
                        S = pr.Sink()
                        S._jet_factor = 1.
                        S._scale_l = scale_l.value
                        #S._scale_v = scale_v.value
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
                        import pdb
                        pdb.set_trace()
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
                            bound_time = global_data['time'][time_it]*units['time_unit'].in_units('yr')
                            delay_time = float((bound_time - formation_time).value)
                            if delay_time == 0:
                                import pdb
                                pdb.set_trace()
                                born_bound = True
                                most_bound_sink_id = str(first_bound_sink)
                            break
                        del res

        if str(sink_id) in True_sink_birth_conditions.keys():
            if True_sink_birth_conditions[str(sink_id)][3] != most_bound_sep or True_sink_birth_conditions[str(sink_id)][5] != delay_time:
                mismatched_inds.append(sink_id)
                print("SHORT CUT DOESN'T WORK FOR SINK_ID", sink_id)
        Sink_bound_birth.append([sink_id, born_bound, most_bound_sink_id, str(first_bound_sink), most_bound_sep, lowest_Etot, delay_time])
        print("Birth conditions of sink", sink_id, "is", Sink_bound_birth[-1])

        file = open("sink_birth_conditions_"+("%03d" % rank)+".pkl", 'wb')
        pickle.dump((Sink_bound_birth),file)
        file.close()
        
        try:
            file = open("mismatched_inds"+("%03d" % rank)+".pkl", 'wb')
            pickle.dump((mismatched_inds),file)
            file.close()
        except:
            pass
        
    sink_id = sink_id + 1

sys.stdout.flush()
CW.Barrier()

#compile pickles
if rank == 0:
    import glob
    birth_pickles = sorted(glob.glob("sink_birth_conditions_*.pkl"))
    Sink_birth_all = {}
    for birth_pick in birth_pickles:
        file = open(birth_pick, 'rb')
        Sink_bound_birth_rank = pickle.load(file)
        file.close()
        for sink_birth_con in Sink_bound_birth_rank:
            Sink_birth_all.update({str(sink_birth_con[0]):sink_birth_con[1:]})
        #os.remove(birth_pick)
        
    file = open("sink_birth_all.pkl", 'wb')
    pickle.dump((Sink_birth_all), file)
    file.close()
    print("Collected sink birth data into sink_birth_all.pkl" )
    
    try:
        mismatched_pickles = sorted(glob.glob("mismatched_inds_*.pkl"))
        mismatched_inds = []
        for mismatched_pickle in mismatched_pickles:
            file = open(mismatched_pickle, 'rb')
            mismatched_inds_rank = pickle.load(file)
            file.close()
            mismatched_inds = mismatched_inds + mismatched_inds_rank
            
        file = open("mismatched_inds.pkl", 'wb')
        pickle.dump((mismatched_inds), file)
        file.close()
        print("Collected mismatched ind into mismatched_inds.pkl" )

    except:
        pass

'''
import pickle
import numpy as np
file_open = open("/groups/astro/rlk/rlk/Global_sink_pickles/G400_full.pkl", "rb")
global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
SFE_5_ind = np.argmin(abs(np.sum(global_data['m'],axis=1)-0.05)) + 1
file_open = open("global_reduced.pkl", "wb")
pickle.dump(({'time':global_data['time'][:SFE_5_ind].T[0],'m': global_data['m'][:SFE_5_ind], 'x': global_data['x'][:SFE_5_ind], 'y': global_data['y'][:SFE_5_ind], 'z': global_data['z'][:SFE_5_ind], 'ux': global_data['ux'][:SFE_5_ind], 'uy': global_data['uy'][:SFE_5_ind], 'uz': global_data['uz'][:SFE_5_ind]}), file_open)
file_open.close()
'''
