import numpy as np
import pickle
import sys
import gc
from mpi4py.MPI import COMM_WORLD as CW

rank = CW.Get_rank()
size = CW.Get_size()


def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

#=====================================================================================================
if rank == 0:
    print("creating units", flush=True)
    sys.stdout.flush()

global_data_pickle_file = sys.argv[1]
means_pickle = sys.argv[2]

#let's try getting all the systems that were found
file = open(means_pickle, 'rb')
superplot_dict, Sink_bound_birth_means, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
file.close()

sys_times = {}
for key in superplot_dict['System_times'].keys():
    sys_times.update({key:superplot_dict['System_times'][key][0]})
system_keys = list(superplot_dict['System_times'].keys())
del superplot_dict
del Sink_bound_birth_means
del Sink_formation_times
del means_dict
del Lifetimes_sys
del Sep_maxs
del Sep_mins
del Initial_Seps
del Final_seps
gc.collect()


Grho = int(global_data_pickle_file.split('/G')[-1].split('/')[0])

#Low_cadence = False
#if 'Low_cadence' in global_data_pickle_file:
#    Low_cadence = True


if Grho == 50:
    scale_m = 1500*1.98841586e+33
    diff_conds = [5,  9, 12, 13, 14, 15, 17, 18, 20, 21, 23, 29, 31, 33, 34, 35, 36, 37, 39, 40, 41, 43, 44, 46, 48, 49, 50, 51, 53, 54, 57, 58, 59, 62, 63, 64, 65, 66, 67, 68, 72, 76, 77, 78, 79, 80, 81, 83, 85]
elif Grho == 100:
    scale_m = 3000*1.98841586e+33
elif Grho == 125:
    scale_m = 3750*1.98841586e+33
elif Grho == 150:
    scale_m = 4500*1.98841586e+33
elif Grho == 200:
    scale_m = 6000*1.98841586e+33
elif Grho == 400:
    scale_m = 12000*1.98841586e+33
else:
    print("MASS UNIT NOT SET", flush=True)
    sys.stdout.flush()
    import pdb
    pdb.set_trace()

scale_l = 1.2342710323849298e+19
scale_l_au = 825059.2245669405
scale_v = 18000.0
scale_t = scale_l/scale_v
scale_t_yr = 21728716.033625457
scale_d = scale_m/(scale_l**3)
del scale_v
gc.collect()

CW.Barrier()

##Loading birth conditions that have already been completed
Sink_bound_birth = []
if rank == 0 :
    import glob
    birth_con_pickles = sorted(glob.glob("sink_birth_conditions_*.pkl"))
    found_sinks = []
    for birth_con_pickle in birth_con_pickles:
        file_open = open(birth_con_pickle, 'rb')
        Sink_bound_birth_rank = pickle.load(file_open)
        file_open.close()
        for birth_con in Sink_bound_birth_rank:
            found_sinks.append(birth_con[0])
        Sink_bound_birth = Sink_bound_birth + Sink_bound_birth_rank
    del birth_con_pickles
    
CW.Barrier()

if rank == 0:
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    file_open = open(global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open)
    file_open.close()
    del file_open
    gc.collect()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    print("Finding formation inds", flush=True)
    sys.stdout.flush()
    ##print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    
    del global_data['x']
    del global_data['y']
    del global_data['z']
    del global_data['ux']
    del global_data['uy']
    del global_data['uz']
    gc.collect()
    
    sink_ids = np.arange(np.shape(global_data['m'].T)[0])
    if len(Sink_bound_birth) > 0:
        sink_ids = sorted(list(set(found_sinks).symmetric_difference(set(sink_ids))))
    del found_sinks
    gc.collect()
    
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    formation_inds = []
    for sink_id in sink_ids:
        if len(np.argwhere(global_data['m'].T[sink_id]>0)) > 0:
            new_ind = np.argwhere(global_data['m'].T[sink_id]>0)[0][0]
        else:
            new_ind = 0
        global_data['m'] = global_data['m'][new_ind:]
        if len(formation_inds) == 0:
            formation_ind = new_ind
        else:
            formation_ind = formation_inds[-1]+new_ind
        formation_inds.append(formation_ind)
    gc.collect()

    print("Found formation inds", flush=True)
    sys.stdout.flush()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    formation_inds = np.array(formation_inds)
    formation_times = global_data['time'][formation_inds]
    del formation_inds
    del global_data
    gc.collect()
    
    file_open = open(global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open)
    file_open.close()
    del file_open

    print("Found formation times", flush=True)
    sys.stdout.flush()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

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
        pickle.dump((sink_ids, formation_times, global_data), file_open)
        file_open.close()
        del form_time_it
        gc.collect()
        print("Saved global_data_rank_"+str(trunc_it)+".pkl", flush=True)
        sys.stdout.flush()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
        
    del formation_times
    del global_data
    gc.collect()
del global_data_pickle_file
#print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
#sys.stdout.flush()
CW.Barrier()
#print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
sink_ids, formation_times, global_data = pickle.load(file_open)
file_open.close()
del file_open
del formation_times
del global_data
gc.collect()
print("loaded global_data_rank_"+str(rank)+".pkl", flush=True)
sys.stdout.flush()
#print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

import pyramses as pr
import multiplicity as m

rit = -1
sink_id = 0
for sink_id in sink_ids:
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    rit = rit + 1
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    if rit == size:
        rit = 0
        #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    if rank == rit:
        #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
        file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
        sink_ids, formation_times, global_data = pickle.load(file_open)
        file_open.close()
        print("Finding birth conditions for sink", sink_id, "on rank", rank)
        sys.stdout.flush()
        del formation_times
        del file_open
        gc.collect()
        #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
        
        #FIND OUT WHY BIRTH CONDITIONS ARE DIFFERENT:
        if Grho == 50:
            if sink_id in diff_conds:
                import pdb
                pdb.set_trace()

        #Calculate energies to find most bound sink
        abspos = np.array([global_data['x'][0][:sink_id+1], global_data['y'][0][:sink_id+1], global_data['z'][0][:sink_id+1]]).T
        absvel = np.array([global_data['ux'][0][:sink_id+1], global_data['uy'][0][:sink_id+1], global_data['uz'][0][:sink_id+1]]).T
        mass = np.array(global_data['m'][0][:sink_id+1])
        n_stars = np.arange(len(mass))
        del global_data
        gc.collect()
        
        rel_pos = abspos[:-1] - abspos[-1]
        del abspos
        gc.collect()
        
        update_seps_neg = np.argwhere(rel_pos<-0.5)
        rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] = rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] + 1.0
        del update_seps_neg
        gc.collect()
        
        update_seps_pos = np.argwhere(rel_pos>0.5)
        rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] = rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] - 1.0
        del update_seps_pos
        gc.collect()

        rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
        del rel_pos
        gc.collect()
            
        rel_vel = absvel[:-1] - absvel[-1]
        del absvel
        gc.collect()
        
        rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
        del rel_vel
        gc.collect()
        
        mtm = mass[-1] * mass[:-1]
        mpm = mass[-1] + mass[:-1]
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
            sink_ids, formation_times, global_data = pickle.load(file_open)
            file_open.close()
            del file_open
            del formation_times
            gc.collect()
            abspos = np.array([global_data['x'][0][n_stars], global_data['y'][0][n_stars], global_data['z'][0][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][0][n_stars], global_data['uy'][0][n_stars], global_data['uz'][0][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][0][n_stars])
            sfe = np.sum(mass)
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
            #del time
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
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho, verbose=False)
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
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
                sys_form_time = sfe
                #Find initial separation
                if str(most_bound_sink_id) != str(first_bound_sink):
                    most_bound_sink_id = str(first_bound_sink)
            del res
            gc.collect()
    
        if np.isnan(sys_id) == False:
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            sink_ids, formation_times, global_data = pickle.load(file_open)
            file_open.close()
            del file_open
            gc.collect()
            
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            curr_it = np.argwhere(sink_ids == sink_id)[0][0]
            next_id = curr_it + size
            del curr_it
            if next_id < len(sink_ids):#len(formation_times):
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
                pickle.dump((sink_ids, formation_times, global_data), file_open)
                file_open.close()
                del form_time_it
            del global_data
            del next_id
            gc.collect()
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
        
        if np.isnan(sys_id):
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            born_bound = False
            most_bound_sep = np.nan
            first_bound_sink = np.nan
            lowest_Etot = np.nan
            delay_time = np.nan
            sys_form_time = np.nan
            
            if len([s for s in system_keys if ' '+str(sink_id)+']' in s]) == 0 or len([s for s in system_keys if '['+str(sink_id)+',' in s]) == 0:
                if Grho == 50:
                    if sink_id in diff_conds:
                        import pdb
                        pdb.set_trace()
                Sink_bound_birth.append([sink_id, born_bound, most_bound_sink_id, str(first_bound_sink), most_bound_sep, lowest_Etot, delay_time, sys_form_time])
                print("Rank:", rank, "Birth conditions of sink", sink_id, "(of", sink_ids[-1],") is", Sink_bound_birth[-1], flush=True)
                sys.stdout.flush()
            else:
                if Grho == 50:
                    if sink_id in diff_conds:
                        import pdb
                        pdb.set_trace()
                #find which system form first
                if len([s for s in system_keys if ' '+str(sink_id)+']' in s]) > 0 and len([s for s in system_keys if '['+str(sink_id)+',' in s]) == 0:
                    first_sys = [s for s in system_keys if ' '+str(sink_id)+']' in s][0]
                elif len([s for s in system_keys if ' '+str(sink_id)+']' in s]) == 0 or len([s for s in system_keys if '['+str(sink_id)+',' in s]) > 0:
                    first_sys = [s for s in system_keys if '['+str(sink_id)+',' in s][0]
                else:
                    if sys_times[[s for s in system_keys if ' '+str(sink_id)+']' in s][0]] < sys_times[[s for s in system_keys if '['+str(sink_id)+',' in s][0]]:
                        first_sys = [s for s in system_keys if ' '+str(sink_id)+']' in s][0]
                    else:
                        first_sys = [s for s in system_keys if '['+str(sink_id)+',' in s][0]
                
                sys_start_time = sys_times[first_sys]
                for key in system_keys:
                    if key == first_sys:
                        break
                    else:
                        del sys_times[key]
                        system_keys.remove(key)
                
                if size == 0:
                    import pdb
                    pdb.set_trace()
                del first_sys
                
                file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
                sink_ids, formation_times, global_data = pickle.load(file_open)
                file_open.close()
                
                first_test_ind = np.argmin(abs(global_data['time']*scale_t_yr - sys_start_time))
                test_time_inds = np.arange(first_test_ind-500, first_test_ind+1)
                del first_test_ind
            
                #units['time_unit'].in_units('yr')
                #formation_time = formation_times[sink_id]*scale_t_yr#units['time_unit'].in_units('yr')
                #test_time_inds = range(len(global_data['x'][time_it:,sink_id]))
                #file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
                #sink_ids, formation_times, global_data = pickle.load(file_open)
                #file_open.close()
                del file_open
                sink_it = np.argwhere(sink_ids == sink_id)[0][0]
                formation_time = formation_times[sink_it]*scale_t_yr
                gc.collect()
                """
                abspos_x = global_data['x'][:]
                abspos_y = global_data['y'][:]
                abspos_z = global_data['z'][:]
                
                del global_data
                gc.collect()
                
                rel_pos_x = abspos_x.T - abspos_x.T[sink_id]
                rel_pos_y = abspos_y.T - abspos_y.T[sink_id]
                rel_pos_z = abspos_z.T - abspos_z.T[sink_id]
                del abspos_x
                del abspos_y
                del abspos_z
                gc.collect()
                
                update_seps_x_neg = np.argwhere(rel_pos_x<-0.5)
                rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] = rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] + 1.0
                del update_seps_x_neg
                gc.collect()
                
                update_seps_x_pos = np.argwhere(rel_pos_x>0.5)
                rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] = rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] - 1.0
                del update_seps_x_pos
                gc.collect()
                
                update_seps_y_neg = np.argwhere(rel_pos_y<-0.5)
                rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] = rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] + 1.0
                del update_seps_y_neg
                gc.collect()
                
                update_seps_y_pos = np.argwhere(rel_pos_y>0.5)
                rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] = rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] - 1.0
                del update_seps_y_pos
                gc.collect()
                
                update_seps_z_neg = np.argwhere(rel_pos_z<-0.5)
                rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] = rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] + 1.0
                del update_seps_z_neg
                gc.collect()
                
                update_seps_z_pos = np.argwhere(rel_pos_z>0.5)
                rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] = rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] - 1.0
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
                
                test_time_inds = np.where((scale_l_au*closest_separations)<10000)[0]
                import pdb
                pdb.set_trace()
                #if Low_cadence == False:
                test_time_inds = test_time_inds[::10]
                del closest_separations
                gc.collect()
                
                file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
                sink_ids, formation_times, global_data = pickle.load(file_open)
                file_open.close()
                """
                #del file_open
                gc.collect()
                #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
                
                curr_it = np.argwhere(sink_ids == sink_id)[0][0]
                next_id = curr_it + size
                del curr_it
                if next_id < len(sink_ids):
                    form_time_it = np.where(global_data['time']==formation_times[next_id])[0][0]
                else:
                    form_time_it = -2
                '''
                next_id = sink_id + size
                if next_id < len(formation_times):
                    form_time_it = np.where(global_data['time']==formation_times[next_id])[0][0]
                else:
                    form_time_it = -2
                '''
                global_test_inds = {}
                global_test_inds.update({'time':global_data['time'][test_time_inds]})
                global_data['time'] = global_data['time'][form_time_it:]
                global_test_inds.update({'m':global_data['m'][test_time_inds]})
                global_data['m'] = global_data['m'][form_time_it:]
                global_test_inds.update({'x':global_data['x'][test_time_inds]})
                global_data['x'] = global_data['x'][form_time_it:]
                global_test_inds.update({'y':global_data['y'][test_time_inds]})
                global_data['y'] = global_data['y'][form_time_it:]
                global_test_inds.update({'z':global_data['z'][test_time_inds]})
                global_data['z'] = global_data['z'][form_time_it:]
                global_test_inds.update({'ux':global_data['ux'][test_time_inds]})
                global_data['ux'] = global_data['ux'][form_time_it:]
                global_test_inds.update({'uy':global_data['uy'][test_time_inds]})
                global_data['uy'] = global_data['uy'][form_time_it:]
                global_test_inds.update({'uz':global_data['uz'][test_time_inds]})
                global_data['uz'] = global_data['uz'][form_time_it:]
                
                file_open = open("global_data_rank_"+str(rank)+".pkl", "wb")
                pickle.dump((sink_ids, formation_times, global_data), file_open)
                file_open.close()
                del form_time_it
                del formation_times
                del global_data
                gc.collect()
                #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
                
                counter = 0
                while np.isnan(first_bound_sink):
                    counter = counter + 1
                    if np.remainder(counter,5000) == 0:
                        print("trying test ind No.", counter, "on rank", rank, flush=True)
                        sys.stdout.flush()
                    if len(global_test_inds['m']) == 0:
                        break
                    else:
                        n_stars = np.where(global_test_inds['m'][0]>0)[0]
                        if len(n_stars)>1:
                            abspos = np.array([global_test_inds['x'][0][n_stars], global_test_inds['y'][0][n_stars], global_test_inds['z'][0][n_stars]]).T#*scale_l
                            absvel = np.array([global_test_inds['ux'][0][n_stars], global_test_inds['uy'][0][n_stars], global_test_inds['uz'][0][n_stars]]).T#*scale_v
                            mass = np.array(global_test_inds['m'][0][n_stars])
                            sfe = np.sum(mass)
                            time = global_test_inds['time'][0]
                            
                            global_test_inds['time'] = global_test_inds['time'][1:]
                            global_test_inds['m'] = global_test_inds['m'][1:]
                            global_test_inds['x'] = global_test_inds['x'][1:]
                            global_test_inds['y'] = global_test_inds['y'][1:]
                            global_test_inds['z'] = global_test_inds['z'][1:]
                            global_test_inds['ux'] = global_test_inds['ux'][1:]
                            global_test_inds['uy'] = global_test_inds['uy'][1:]
                            global_test_inds['uz'] = global_test_inds['uz'][1:]
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
                            #del time
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
                            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
                            res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho)
                            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
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
                                delay_time = float(bound_time - formation_time)
                                sys_form_time = sfe
                                #Find initial separation
                                if delay_time == 0:
                                    born_bound = True
                                    most_bound_sink_id = str(first_bound_sink)
                                break
                            del res
                            gc.collect()
        
                Sink_bound_birth.append([sink_id, born_bound, most_bound_sink_id, str(first_bound_sink), most_bound_sep, lowest_Etot, delay_time, sys_form_time])
                print("Rank:", rank, "Birth conditions of sink", sink_id, "(of", sink_ids[-1],") is", Sink_bound_birth[-1], flush=True)
                sys.stdout.flush()

        file = open("sink_birth_conditions_"+("%03d" % rank)+".pkl", 'wb')
        pickle.dump((Sink_bound_birth),file)
        file.close()

    sink_id = sink_id + 1

#sys.stdout.flush()
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
    print("Collected sink birth data into sink_birth_all.pkl", flush=True)
    sys.stdout.flush()


if rank==0:
    import collections
    def flatten(x):
        if isinstance(x, collections.Iterable):
            return [a for i in x for a in flatten(i)]
        else:
            return [x]
    file = open('sink_birth_all.pkl', 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()

    for sys_key in Sink_birth_all.keys():
        if Sink_birth_all[sys_key][0]==False and np.isnan(Sink_birth_all[sys_key][3])==False and Sink_birth_all[sys_key][1] in flatten(eval(Sink_birth_all[sys_key][2])):
            if np.sum(np.array(flatten(eval(Sink_birth_all[sys_key][2])))>eval(sys_key)) == 0:
                print('For', sys_key, 'Changing', Sink_birth_all[sys_key][1], 'to', Sink_birth_all[sys_key][2])
                sys.stdout.flush()
                Sink_birth_all[sys_key][1] = Sink_birth_all[sys_key][2]
            else:
                print('delayed core fragmentation system connects to system with new sink IDS. Birth con for', sys_key, ':', Sink_birth_all[sys_key])
                sys.stdout.flush()

    file = open("sink_birth_all_delayed_core_frag_cleaned.pkl", 'wb')
    pickle.dump((Sink_birth_all), file)
    file.close()
