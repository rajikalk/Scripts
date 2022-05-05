import numpy as np
import pickle
import pyramses as pr
#from pyramses import rsink
import multiplicity as m
from mpi4py.MPI import COMM_WORLD as CW
import sys
import gc

rank = CW.Get_rank()
size = CW.Get_size()

def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

'''
def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    args = parser.parse_args()
    return args
'''

#=====================================================================================================

global_data_pickle_file = sys.argv[1]
#args = parse_inputs()

#==========================================================================================

if rank == 0:
    print("creating units")

simulation_density_id = global_data_pickle_file.split('/G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50.
    scale_m = 1500*1.98841586e+33
elif simulation_density_id == '100':
    Grho=100.
    scale_m = 3000*1.98841586e+33
elif simulation_density_id == '125':
    Grho=125.
    scale_m = 3750*1.98841586e+33
elif simulation_density_id == '150':
    Grho=150.
    scale_m = 4500*1.98841586e+33
elif simulation_density_id == '200':
    Grho=200.
    scale_m = 6000*1.98841586e+33
elif simulation_density_id == '400':
    Grho=400.
    scale_m = 12000*1.98841586e+33
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
del simulation_density_id
gc.collect()

scale_l = 1.2342710323849298e+19
scale_l_au = 825059.2245669405
scale_v = 18000.0
scale_t = scale_l/scale_v
scale_t_yr = 21728716.033625457
scale_d = scale_m/(scale_l**3)
del scale_v
gc.collect()

if rank == 0:
    file_open = open(global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open)
    file_open.close()
    del file_open
    gc.collect()

    print("Finding formation inds")

    formation_inds = [0]
    for sink_id in range(1, np.shape(global_data['m'].T)[0]):
        formation_inds.append(np.argwhere(global_data['m'].T[sink_id]>0)[0][0])

    print("Found formation inds")

    formation_inds = np.array(formation_inds)
    formation_times = global_data['time'][formation_inds]
    del formation_inds
    gc.collect()

    print("Found formation times")
    
    for trunc_it in range(size):
        form_time_it = np.where(global_data['time']==formation_times[trunc_it])[0][0]
        
        #truncate global data
        global_data['time'] = global_data['time'][form_time_it:]
        global_data['m'] = global_data['m'][form_time_it:]
        global_data['x'] = global_data['x'][form_time_it:]
        global_data['y'] = global_data['y'][form_time_it:]
        global_data['z'] = global_data['z'][form_time_it:]
        global_data['ux'] = global_data['ux'][form_time_it:]
        global_data['uy'] = global_data['uy'][form_time_it:]
        global_data['uz'] = global_data['uz'][form_time_it:]
        
        file_open = open("global_data_rank_"+str(trunc_it)+".pkl", "wb")
        pickle.dump((formation_times, global_data), file_open)
        file_open.close()
        del form_time_it
        gc.collect()
        
    del formation_times
    del global_data
    gc.collect()
sys.stdout.flush()
CW.Barrier()

file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
formation_times, global_data = pickle.load(file_open)
file_open.close()
del file_open
del global_data
gc.collect()
Sink_bound_birth = []

sys.stdout.flush()
CW.Barrier()

rit = -1
sink_id = 0
while sink_id < len(formation_times):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
        formation_times, global_data = pickle.load(file_open)
        file_open.close()
        del file_open
        gc.collect()
        
        #Calculate energies to find most bound sink
        new_sink_pos = np.array([global_data['x'][0][sink_id], global_data['y'][0][sink_id], global_data['z'][0][sink_id]]).T
        abspos = np.array([global_data['x'][0][:sink_id], global_data['y'][0][:sink_id], global_data['z'][0][:sink_id]]).T
        new_sink_vel = np.array([global_data['ux'][0][sink_id], global_data['uy'][0][sink_id], global_data['uz'][0][sink_id]]).T
        absvel = np.array([global_data['ux'][0][:sink_id], global_data['uy'][0][:sink_id], global_data['uz'][0][:sink_id]]).T
        new_sink_mass = np.array(global_data['m'][0][sink_id])
        mass = np.array(global_data['m'][0][:sink_id])
        n_stars = np.where(global_data['m'][0]>0)[0]
        del global_data
        gc.collect()
        
        rel_pos = abspos - new_sink_pos
        del new_sink_pos
        del abspos
        gc.collect()
        update_seps_neg = np.argwhere(rel_pos<-0.5)
        update_seps_pos = np.argwhere(rel_pos>0.5)
        
        rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] = rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] + 1.0
        del update_seps_neg
        gc.collect()
        rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] = rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] - 1.0
        del update_seps_pos
        gc.collect()

        rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
        del rel_pos
        gc.collect()
        
        rel_vel = absvel - new_sink_vel
        del new_sink_vel
        del absvel
        gc.collect()
        rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
        del rel_vel
        gc.collect()
        
        mtm = new_sink_mass * mass
        mpm = new_sink_mass + mass
        del new_sink_mass
        del mass
        gc.collect()
        
        newtonianPotential = -1./rel_sep
        
        Ekin = 0.5 * mtm/mpm * rel_speed**2
        del rel_speed
        del mpm
        gc.collect()
        Epot = Grho * mtm * newtonianPotential
        del newtonianPotential
        del mtm
        gc.collect()
        Etot = Ekin + Epot
        try:
            most_bound_sink_id = np.argmin(Etot)
        except:
            most_bound_sink_id = np.nan
        del Ekin
        del Epot
        gc.collect()
        
        born_bound = True
        delay_time = 0
        
        sys_id = np.nan
        if len(n_stars)>1:
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = pickle.load(file_open)
            file_open.close()
            del file_open
            gc.collect()
            abspos = np.array([global_data['x'][0][n_stars], global_data['y'][0][n_stars], global_data['z'][0][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][0][n_stars], global_data['uy'][0][n_stars], global_data['uz'][0][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][0][n_stars])
            del n_stars
            
            time = global_data['time'][0]
            del global_data
            gc.collect()
            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l
            #S._scale_v = scale_v.value
            S._scale_t = scale_t
            S._scale_d = scale_d
            S._time = time
            del time
            gc.collect()
            S._abspos = abspos
            del abspos
            gc.collect()
            S._absvel = absvel
            del absvel
            gc.collect()
            S._mass = mass
            del mass
            gc.collect()
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
                if str(most_bound_sink_id) != str(first_bound_sink):
                    most_bound_sink_id = str(first_bound_sink)
            del res
            gc.collect()
        
        if np.isnan(sys_id) == False:
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = pickle.load(file_open)
            file_open.close()
            del file_open
            gc.collect()
            next_id = sink_id + size
            if next_id < len(formation_times):
                form_time_it = np.where(global_data['time']==formation_times[next_id])[0][0]
            
                #truncate global data
                global_data['time'] = global_data['time'][form_time_it:]
                global_data['m'] = global_data['m'][form_time_it:]
                global_data['x'] = global_data['x'][form_time_it:]
                global_data['y'] = global_data['y'][form_time_it:]
                global_data['z'] = global_data['z'][form_time_it:]
                global_data['ux'] = global_data['ux'][form_time_it:]
                global_data['uy'] = global_data['uy'][form_time_it:]
                global_data['uz'] = global_data['uz'][form_time_it:]
                
                file_open = open("global_data_rank_"+str(rank)+".pkl", "wb")
                pickle.dump((formation_times, global_data), file_open)
                file_open.close()
                del form_time_it
            del global_data
            gc.collect()
                

        if np.isnan(sys_id):
            born_bound = False
            most_bound_sep = np.nan
            first_bound_sink = np.nan
            lowest_Etot = np.nan
            delay_time = np.nan
            
            formation_time = formation_times[sink_id]*scale_t_yr#units['time_unit'].in_units('yr')
            #test_time_inds = range(len(global_data['x'][time_it:,sink_id]))
            
            
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = pickle.load(file_open)
            file_open.close()
            del file_open
            gc.collect()
            
            new_sink_pos_x = global_data['x'][:,sink_id]
            new_sink_pos_y = global_data['y'][:,sink_id]
            new_sink_pos_z = global_data['z'][:,sink_id]
            
            abspos_x = global_data['x'][:]
            abspos_y = global_data['y'][:]
            abspos_z = global_data['z'][:]
            
            del global_data
            gc.collect()
            
            rel_pos_x = abspos_x.T - new_sink_pos_x
            rel_pos_y = abspos_y.T - new_sink_pos_y
            rel_pos_z = abspos_z.T - new_sink_pos_z
            del new_sink_pos_x
            del new_sink_pos_y
            del new_sink_pos_z
            del abspos_x
            del abspos_y
            del abspos_z
            gc.collect()
            
            update_seps_x_neg = np.argwhere(rel_pos_x<-0.5)
            update_seps_x_pos = np.argwhere(rel_pos_x>0.5)
            rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] = rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] + 1.0
            rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] = rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] - 1.0
            del update_seps_x_neg
            del update_seps_x_pos
            gc.collect()
            
            update_seps_y_neg = np.argwhere(rel_pos_y<-0.5)
            update_seps_y_pos = np.argwhere(rel_pos_y>0.5)
            rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] = rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] + 1.0
            rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] = rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] - 1.0
            del update_seps_y_neg
            del update_seps_y_pos
            gc.collect()
            
            update_seps_z_neg = np.argwhere(rel_pos_z<-0.5)
            update_seps_z_pos = np.argwhere(rel_pos_z>0.5)
            rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] = rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] + 1.0
            rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] = rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] - 1.0
            del update_seps_z_neg
            del update_seps_z_pos
            gc.collect()
            
            rel_sep = np.sqrt(rel_pos_x**2 + rel_pos_y**2 + rel_pos_z**2)
            del rel_pos_x
            del rel_pos_y
            del rel_pos_z
            gc.collect()
            
            rel_sep[np.where(rel_sep == 0)] = np.inf
            closest_separations = np.min(rel_sep, axis=0)
            del rel_sep
            gc.collect()
            #import pdb
            #pdb.set_trace()
            test_time_inds = np.where((scale_l_au*closest_separations)<10000)[0]
            del closest_separations
            gc.collect()
            
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = pickle.load(file_open)
            file_open.close()
            del file_open
            gc.collect()
            
            global_test_inds = {}
            global_test_inds.update({'time':global_data['time'][test_time_inds]})
            global_test_inds.update({'m':global_data['m'][test_time_inds]})
            global_test_inds.update({'x':global_data['x'][test_time_inds]})
            global_test_inds.update({'y':global_data['y'][test_time_inds]})
            global_test_inds.update({'z':global_data['z'][test_time_inds]})
            global_test_inds.update({'ux':global_data['ux'][test_time_inds]})
            global_test_inds.update({'uy':global_data['uy'][test_time_inds]})
            global_test_inds.update({'uz':global_data['uz'][test_time_inds]})
            
            next_id = sink_id + size
            if next_id < len(formation_times):
                form_time_it = np.where(global_data['time']==formation_times[next_id])[0][0]
            
                #truncate global data
                global_data['time'] = global_data['time'][form_time_it:]
                global_data['m'] = global_data['m'][form_time_it:]
                global_data['x'] = global_data['x'][form_time_it:]
                global_data['y'] = global_data['y'][form_time_it:]
                global_data['z'] = global_data['z'][form_time_it:]
                global_data['ux'] = global_data['ux'][form_time_it:]
                global_data['uy'] = global_data['uy'][form_time_it:]
                global_data['uz'] = global_data['uz'][form_time_it:]
                
                file_open = open("global_data_rank_"+str(rank)+".pkl", "wb")
                pickle.dump((formation_times, global_data), file_open)
                file_open.close()
                del form_time_it
            del global_data
            gc.collect()
            
            for time_it in range(len(global_test_inds['m'])):
                if np.isnan(first_bound_sink):
                    if np.remainder(time_it,5000) == 0:
                        print("testing time_it", time_it, "on rank", rank)
                    n_stars = np.where(global_test_inds['m'][time_it]>0)[0]
                    if len(n_stars)>1:
                        abspos = np.array([global_test_inds['x'][time_it][n_stars], global_test_inds['y'][time_it][n_stars], global_test_inds['z'][time_it][n_stars]]).T#*scale_l
                        absvel = np.array([global_test_inds['ux'][time_it][n_stars], global_test_inds['uy'][time_it][n_stars], global_test_inds['uz'][time_it][n_stars]]).T#*scale_v
                        mass = np.array(global_test_inds['m'][time_it][n_stars])
                        time = global_test_inds['time'][time_it]
                        
                        #Remove global data:
                        del n_stars
                        gc.collect()
                        S = pr.Sink()
                        S._jet_factor = 1.
                        S._scale_l = scale_l
                        #S._scale_v = scale_v.value
                        S._scale_t = scale_t
                        S._scale_d = scale_d
                        S._time = time
                        del time
                        gc.collect()
                        S._abspos = abspos
                        del abspos
                        gc.collect()
                        S._absvel = absvel
                        del absvel
                        gc.collect()
                        S._mass = mass
                        del mass
                        gc.collect()
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
                            bound_time = res['time']*scale_t_yr
                            delay_time = float((bound_time - formation_time).value)
                            if delay_time == 0:
                                import pdb
                                pdb.set_trace()
                                born_bound = True
                                most_bound_sink_id = str(first_bound_sink)
                            break
                        del res
                        gc.collect()

        Sink_bound_birth.append([sink_id, born_bound, most_bound_sink_id, str(first_bound_sink), most_bound_sep, lowest_Etot, delay_time])
        print("Birth conditions of sink", sink_id, "is", Sink_bound_birth[-1])

        file = open("sink_birth_conditions_"+("%03d" % rank)+".pkl", 'wb')
        pickle.dump((Sink_bound_birth),file)
        file.close()

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
