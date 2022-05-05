import numpy as np
#import pickle
from pickle import load, dump
import pyramses as pr
import multiplicity as m
from mpi4py.MPI import COMM_WORLD as CW
from sys import argv, stdout
from gc import collect
from psutil import virtual_memory
from inspect import currentframe, getframeinfo

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
    print("creating units")
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

global_data_pickle_file = argv[1]
Grho = int(global_data_pickle_file.split('/G')[-1].split('/')[0])

if Grho == 50:
    scale_m = 1500*1.98841586e+33
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
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()

scale_l = 1.2342710323849298e+19
scale_l_au = 825059.2245669405
scale_v = 18000.0
scale_t = scale_l/scale_v
scale_t_yr = 21728716.033625457
scale_d = scale_m/(scale_l**3)
del scale_v
collect()

stdout.flush()
CW.Barrier()

if rank == 0:
    file_open = open(global_data_pickle_file, 'rb')
    global_data = load(file_open)
    file_open.close()
    del file_open
    collect()

    print("Finding formation inds")
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    
    del global_data['x']
    del global_data['y']
    del global_data['z']
    del global_data['ux']
    del global_data['uy']
    del global_data['uz']
    collect()

    formation_inds = [0]
    for sink_id in range(1, np.shape(global_data['m'].T)[0]):
        new_ind = np.argwhere(global_data['m'].T[sink_id]>0)[0][0]
        global_data['m'] = global_data['m'][new_ind:]
        formation_ind = formation_inds[-1]+new_ind
        formation_inds.append(formation_ind)

    print("Found formation inds")
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    formation_inds = np.array(formation_inds)
    formation_times = global_data['time'][formation_inds]
    del formation_inds
    del global_data
    collect()
    
    file_open = open(global_data_pickle_file, 'rb')
    global_data = load(file_open)
    file_open.close()
    del file_open

    print("Found formation times")
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
        dump((formation_times, global_data), file_open)
        file_open.close()
        del form_time_it
        collect()
    print("Saved global_data for each rank")
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
        
    del formation_times
    del global_data
    collect()
del global_data_pickle_file
stdout.flush()
CW.Barrier()

file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
formation_times, global_data = load(file_open)
file_open.close()
del file_open
del global_data
collect()
Sink_bound_birth = []
if rank == 0:
    print("loaded formation_times")
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

stdout.flush()
CW.Barrier()

rit = -1
sink_id = 0
while sink_id < len(formation_times):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
        formation_times, global_data = load(file_open)
        file_open.close()
        del file_open
        collect()
        print("loaded global data on rank", rank)
        #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
        
        #Calculate energies to find most bound sink
        abspos = np.array([global_data['x'][0][:sink_id+1], global_data['y'][0][:sink_id+1], global_data['z'][0][:sink_id+1]]).T
        absvel = np.array([global_data['ux'][0][:sink_id+1], global_data['uy'][0][:sink_id+1], global_data['uz'][0][:sink_id+1]]).T
        mass = np.array(global_data['m'][0][:sink_id+1])
        n_stars = np.arange(len(mass))
        del global_data
        collect()
        
        rel_pos = abspos[:-1] - abspos[-1]
        del abspos
        collect()
        
        update_seps_neg = np.argwhere(rel_pos<-0.5)
        rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] = rel_pos[update_seps_neg.T[0], update_seps_neg.T[1]] + 1.0
        del update_seps_neg
        collect()
        
        update_seps_pos = np.argwhere(rel_pos>0.5)
        rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] = rel_pos[update_seps_pos.T[0], update_seps_pos.T[1]] - 1.0
        del update_seps_pos
        collect()

        rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
        del rel_pos
        collect()
            
        rel_vel = absvel[:-1] - absvel[-1]
        del absvel
        collect()
        
        rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
        del rel_vel
        collect()
        
        mtm = mass[-1] * mass[:-1]
        mpm = mass[-1] + mass[:-1]
        del mass
        collect()
        
        newtonianPotential = -1./rel_sep
        
        Ekin = 0.5 * mtm/mpm * rel_speed**2
        del rel_speed
        del mpm
        collect()
        
        Epot = Grho * mtm * newtonianPotential
        del newtonianPotential
        del mtm
        collect()

        Etot = Ekin + Epot
        try:
            most_bound_sink_id = np.argmin(Etot)
        except:
            most_bound_sink_id = np.nan
        del Ekin
        del Epot
        collect()
        
        born_bound = True
        delay_time = 0
        
        sys_id = np.nan
        if len(n_stars)>1:
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = load(file_open)
            file_open.close()
            del file_open
            collect()
            abspos = np.array([global_data['x'][0][n_stars], global_data['y'][0][n_stars], global_data['z'][0][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][0][n_stars], global_data['uy'][0][n_stars], global_data['uz'][0][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][0][n_stars])
            del n_stars
            
            time = global_data['time'][0]
            del global_data
            collect()

            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l
            #S._scale_v = scale_v.value
            S._scale_t = scale_t
            S._scale_d = scale_d
            S._time = time
            del time
            collect()
            S._abspos = abspos
            del abspos
            collect()
            S._absvel = absvel
            del absvel
            collect()
            S._mass = mass
            del mass
            collect()
            print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho, verbose=False)
            print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
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
            collect()
    
        if np.isnan(sys_id) == False:
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = load(file_open)
            file_open.close()
            del file_open
            collect()
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
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
                dump((formation_times, global_data), file_open)
                file_open.close()
                del form_time_it
            del global_data
            collect()
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
                
        if np.isnan(sys_id):
            born_bound = False
            most_bound_sep = np.nan
            first_bound_sink = np.nan
            lowest_Etot = np.nan
            delay_time = np.nan
            
            formation_time = formation_times[sink_id]*scale_t_yr#units['time_unit'].in_units('yr')
            #test_time_inds = range(len(global_data['x'][time_it:,sink_id]))
            
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = load(file_open)
            file_open.close()
            del file_open
            collect()
            
            abspos_x = global_data['x'][:]
            abspos_y = global_data['y'][:]
            abspos_z = global_data['z'][:]
            
            del global_data
            collect()
            
            rel_pos_x = abspos_x.T - abspos_x.T[sink_id]
            rel_pos_y = abspos_y.T - abspos_y.T[sink_id]
            rel_pos_z = abspos_z.T - abspos_z.T[sink_id]
            del abspos_x
            del abspos_y
            del abspos_z
            collect()
            
            update_seps_x_neg = np.argwhere(rel_pos_x<-0.5)
            rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] = rel_pos_x[update_seps_x_neg.T[0], update_seps_x_neg.T[1]] + 1.0
            del update_seps_x_neg
            collect()
            
            update_seps_x_pos = np.argwhere(rel_pos_x>0.5)
            rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] = rel_pos_x[update_seps_x_pos.T[0], update_seps_x_pos.T[1]] - 1.0
            del update_seps_x_pos
            collect()
            
            update_seps_y_neg = np.argwhere(rel_pos_y<-0.5)
            rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] = rel_pos_y[update_seps_y_neg.T[0], update_seps_y_neg.T[1]] + 1.0
            del update_seps_y_neg
            collect()
            
            update_seps_y_pos = np.argwhere(rel_pos_y>0.5)
            rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] = rel_pos_y[update_seps_y_pos.T[0], update_seps_y_pos.T[1]] - 1.0
            del update_seps_y_pos
            collect()
            
            update_seps_z_neg = np.argwhere(rel_pos_z<-0.5)
            rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] = rel_pos_z[update_seps_z_neg.T[0], update_seps_z_neg.T[1]] + 1.0
            del update_seps_z_neg
            collect()
            
            update_seps_z_pos = np.argwhere(rel_pos_z>0.5)
            rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] = rel_pos_z[update_seps_z_pos.T[0], update_seps_z_pos.T[1]] - 1.0
            del update_seps_z_pos
            collect()
            
            rel_sep = np.sqrt(rel_pos_x**2 + rel_pos_y**2 + rel_pos_z**2)
            del rel_pos_x
            del rel_pos_y
            del rel_pos_z
            collect()
            
            rel_sep[np.where(rel_sep == 0)] = np.inf
            closest_separations = np.min(rel_sep, axis=0)
            del rel_sep
            collect()
            
            test_time_inds = np.where((scale_l_au*closest_separations)<10000)[0]
            del closest_separations
            collect()
            
            file_open = open("global_data_rank_"+str(rank)+".pkl", 'rb')
            formation_times, global_data = load(file_open)
            file_open.close()
            del file_open
            collect()
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            
            next_id = sink_id + size
            if next_id < len(formation_times):
                form_time_it = np.where(global_data['time']==formation_times[next_id])[0][0]
            else:
                form_time_it = -2
            
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
            dump((formation_times, global_data), file_open)
            file_open.close()
            del form_time_it
            del global_data
            collect()
            #print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
            
            counter = 0
            while np.isnan(first_bound_sink):
                counter = counter + 1
                if np.remainder(counter,5000) == 0:
                    print("testing time_it", time_it, "on rank", rank)
                n_stars = np.where(global_test_inds['m'][0]>0)[0]
                if len(n_stars)>1:
                    abspos = np.array([global_test_inds['x'][0][n_stars], global_test_inds['y'][0][n_stars], global_test_inds['z'][0][n_stars]]).T#*scale_l
                    absvel = np.array([global_test_inds['ux'][0][n_stars], global_test_inds['uy'][0][n_stars], global_test_inds['uz'][0][n_stars]]).T#*scale_v
                    mass = np.array(global_test_inds['m'][0][n_stars])
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
                    collect()
                    S = pr.Sink()
                    S._jet_factor = 1.
                    S._scale_l = scale_l
                    #S._scale_v = scale_v.value
                    S._scale_t = scale_t
                    S._scale_d = scale_d
                    S._time = time
                    del time
                    collect()
                    S._abspos = abspos
                    del abspos
                    collect()
                    S._absvel = absvel
                    del absvel
                    collect()
                    S._mass = mass
                    del mass
                    collect()
                    print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
                    res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho)
                    print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
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
                        if delay_time == 0:
                            born_bound = True
                            most_bound_sink_id = str(first_bound_sink)
                        break
                    del res
                    collect()
                    print("Memory_useage on rank", rank,":", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

        Sink_bound_birth.append([sink_id, born_bound, most_bound_sink_id, str(first_bound_sink), most_bound_sep, lowest_Etot, delay_time])
        print("Birth conditions of sink", sink_id, "is", Sink_bound_birth[-1])

        file = open("sink_birth_conditions_"+("%03d" % rank)+".pkl", 'wb')
        dump((Sink_bound_birth),file)
        file.close()

    sink_id = sink_id + 1

stdout.flush()
CW.Barrier()

#compile pickles
if rank == 0:
    import glob
    birth_pickles = sorted(glob.glob("sink_birth_conditions_*.pkl"))
    Sink_birth_all = {}
    for birth_pick in birth_pickles:
        file = open(birth_pick, 'rb')
        Sink_bound_birth_rank = load(file)
        file.close()
        for sink_birth_con in Sink_bound_birth_rank:
            Sink_birth_all.update({str(sink_birth_con[0]):sink_birth_con[1:]})
        #os.remove(birth_pick)
        
    file = open("sink_birth_all.pkl", 'wb')
    dump((Sink_birth_all), file)
    file.close()
    print("Collected sink birth data into sink_birth_all.pkl" )
