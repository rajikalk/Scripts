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
    
def is_hierarchical(sys_structure):
    close_braket = False
    for char in str(sys_structure):
        if char == ']':
            close_braket = True
        if char == '[' and close_braket == True:
            is_hierarchical = False
            break
        else:
            is_hierarchical = True
    return is_hierarchical

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-update", "--update_pickles", help="Do you want to remake the pickles?", type=str, default='True')
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=35.6)
    parser.add_argument("-lower_L", "--lower_L_limit", help="What is the upper Luminosity limit?", type=float, default=0.04)
    parser.add_argument("-cyclic", "--cyclic_bool", help="Do you want to set cyclic?", type=str, default='True')
    parser.add_argument("-bound", "--bound_check", help="Do you actually want to analyse bound systems?", type=str, default='True')
    parser.add_argument("-start_time_it", "--start_time_index", help="What time index do you want to start at? mostly for diagnostic reasons.", type=int, default=0)
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifeimteimte threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("-replace_ind_type", "--replace_index_type", help="Replacing subsystems with a primary, do you want to use the tru mass, or the  observe brightness?", type=str, default="L")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/G')[-1].split('/')[0]

if simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

sys.stdout.flush()
CW.Barrier()

Times = []
SFE = []
SFE_n = []
M_tot = []
M_tot_multi = []
N_stars = []
N_vis_stars = []
N_multi_stars = []
System_seps = {}
System_semimajor = {}
System_times = {}

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
Accretion_array = dm/dt
print('loaded global data')

sys.stdout.flush()
CW.Barrier()

pickle_file = args.pickled_file

if rank == 0 and args.update_pickles == 'False':
    try:
        file = open(pickle_file+'.pkl', 'rb')
        Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times = pickle.load(file)
        file.close()
        print('read in completed pickle')
    except:
        pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
        if len(pickle_files) > 0:
            print('reading pickles from ranks')
            Times_full = []
            SFE_full = []
            SFE_n_full = []
            M_tot_full = []
            M_tot_multi_full = []
            N_stars_full = []
            N_vis_stars_full = []
            N_multi_stars_full = []
            System_seps_full = {}
            System_times_full = {}
            for pick_file in pickle_files:
                file = open(pick_file, 'rb')
                Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times = pickle.load(file)
                file.close()
                for time_key in System_times.keys():
                    if time_key not in System_times_full.keys():
                        System_times_full.update({time_key:System_times[time_key]})
                        System_seps_full.update({time_key:System_seps[time_key]})
                    else:
                        System_times_full[time_key] = System_times_full[time_key] + System_times[time_key]
                        System_seps_full[time_key] = System_seps_full[time_key] + System_seps[time_key]
                Times_full = Times_full + Times
                SFE_full = SFE_full + SFE
                SFE_n_full = SFE_n_full + SFE_n
                M_tot_full = M_tot_full + M_tot
                M_tot_multi_full = M_tot_multi_full + M_tot_multi
                N_stars_full = N_stars_full + N_stars
                N_vis_stars_full = N_vis_stars_full + N_vis_stars
                N_multi_stars_full = N_multi_stars_full + N_multi_stars
                os.remove(pick_file)
            
            #Let's sort the data
            for time_key in System_times_full.keys():
                sorted_array = np.array(System_times_full[time_key])
                dict_sort = np.argsort(sorted_array)
                sorted_array = np.array(System_times_full[time_key])[dict_sort]
                System_times[time_key] = sorted_array.tolist()
                sorted_array = np.array(System_seps_full[time_key])[dict_sort]
                System_seps[time_key] = sorted_array.tolist()
            sorted_inds = np.argsort(Times_full)
            Times = np.array(Times_full)[sorted_inds].tolist()
            SFE = np.array(SFE_full)[sorted_inds].tolist()
            SFE_n = np.array(SFE_n_full)[sorted_inds].tolist()
            M_tot = np.array(M_tot_full)[sorted_inds].tolist()
            M_tot_multi = np.array(M_tot_multi_full)[sorted_inds].tolist()
            N_stars = np.array(N_stars_full)[sorted_inds].tolist()
            N_vis_stars = np.array(N_vis_stars_full)[sorted_inds].tolist()
            N_multi_stars = np.array(N_multi_stars_full)[sorted_inds].tolist()
            
            file = open(pickle_file+'.pkl', 'wb')
            pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times),file)
            file.close()
        else:
            print('not pickles exist')
            
sys.stdout.flush()
CW.Barrier()

start_time_it = args.start_time_index

if args.cyclic_bool == 'False':
    cyclic_arg = False
else:
    cyclic_arg = True
if args.bound_check == 'True':
    bound_check = True
else:
    bound_check = False
if args.update_pickles == 'True':
    rit = -1
    for time_it in range(start_time_it, len(global_data['time'].T[0])):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            #preamble
            n_stars = np.where(global_data['m'][time_it]>0)[0]
            if len(n_stars) > 1:
                abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                mass = np.array(global_data['m'][time_it][n_stars])
                time = global_data['time'][time_it][n_stars][0]
                time_yr = yt.YTQuantity(scale_t*time, 's').in_units('yr').value
                Times.append(time_yr)
                
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
                
                L_tot = luminosity(global_data, n_stars, time_it)
                M_dot = accretion(n_stars, time_it)
                vis_inds = np.where((L_tot>=args.lower_L_limit)&(M_dot>=args.accretion_limit)&(L_tot<=args.upper_L_limit))[0]
                N_vis_stars.append(len(vis_inds))
                
                N_stars.append(len(n_stars))
                M_tot_msun = np.sum(mass*units['mass_unit'].in_units('Msun'))
                M_tot.append(M_tot_msun)
                SFE_val = M_tot_msun/units['mass_unit'].in_units('Msun')
                SFE.append(SFE_val)
                SFE_n.append(SFE_val/time)
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=bound_check, nmax=6, cyclic=cyclic_arg)
                if True in (res['separation']>10000):
                    import pdb
                    pdb.set_trace()
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
                n_multi = np.sum(res['n'][multi_inds])
                M_multi = np.sum(res['mass'][multi_inds])
                updated_systems = []
                for multi_ind in multi_inds:
                    sys_comps = losi(multi_ind, res)
                    sys_string = sorted(flatten(sys_comps))
                    
                    visible_stars = list(set(sys_string).intersection(set(vis_inds)))
                    sys_comps = str(sys_comps)
                    if len(visible_stars) == 0:
                        sys_comps = ""
                    elif len(visible_stars) != len(sys_string):
                        invisible_stars = list(set(sys_string).symmetric_difference(set(visible_stars)))
                        for inv_ind in invisible_stars:
                            if ", "+str(inv_ind)+"]" in sys_comps:
                                split_ind = sys_comps.index(", "+str(inv_ind)+"]")
                                open_brackets = []
                                for char_it in range(len(sys_comps)):
                                    if sys_comps[char_it] == '[':
                                        open_brackets.append(char_it)
                                    if sys_comps[char_it] == ']':
                                        if char_it > split_ind:
                                            bracket_ind = open_brackets.pop()
                                            break
                                        else:
                                            op_ind = open_brackets.pop()
                                sys_comps = sys_comps[:bracket_ind] + sys_comps[bracket_ind+1:]
                                sys_comps = "".join(sys_comps.split(", "+str(inv_ind)+"]"))
                            elif "["+str(inv_ind)+", " in sys_comps:
                                split_ind = sys_comps.index("["+str(inv_ind)+", ")
                                open_brackets = []
                                for char_it in range(len(sys_comps)):
                                    if sys_comps[char_it] == '[':
                                        open_brackets.append(char_it)
                                    if sys_comps[char_it] == ']':
                                        op_ind = open_brackets.pop()
                                        if op_ind == split_ind:
                                            bracket_ind = char_it
                                            break
                                sys_comps = sys_comps[:bracket_ind] + sys_comps[bracket_ind+1:]
                                sys_comps = "".join(sys_comps.split("["+str(inv_ind)+", "))
                            else:
                                import pdb
                                pdb.set_trace()
                            
                    sys_comps_str = str(sys_comps)
                        
                    if len(visible_stars) > 1:
                        updated_systems.append(sys_comps_str)
                        sep_arr = []
                        semimajor_arr = []
                        reduced = False
                        while reduced == False:
                            open_braket_ind = []
                            for char_it in range(len(sys_comps)):
                                if sys_comps[char_it] == '[':
                                    open_braket_ind.append(char_it)
                                if sys_comps[char_it] == ']':
                                    open_ind = open_braket_ind.pop()
                                    sub_sys = eval(sys_comps[open_ind:char_it+1])
                                    if len(sub_sys) == 2:
                                        ind_1 = sub_sys[0]
                                        ind_2 = sub_sys[1]
                                        pos_diff = res['abspos'][ind_1] - res['abspos'][ind_2]
                                        sep_value = np.sqrt(np.sum(pos_diff**2))
                                        position = yt.YTArray(res['abspos'][np.array([ind_1,ind_2])].T, 'au')
                                        if sep_value > (scale_l.in_units('AU')/2):
                                            update_inds = np.where(abs(pos_diff)>scale_l.in_units('AU')/2)[0]
                                            for ind in update_inds:
                                                if pos_diff[ind] < 0:
                                                    pos_diff[ind] = pos_diff[ind] + scale_l.in_units('AU').value
                                                else:
                                                    pos_diff[ind] = pos_diff[ind] - scale_l.in_units('AU').value
                                            sep_value = np.sqrt(np.sum(pos_diff**2))
                                            for ind in update_inds:
                                                pos_update_ind = np.where(position[ind]<scale_l.in_units('AU')/2)[0]
                                                position[ind][pos_update_ind] = position[ind][pos_update_ind] + scale_l.in_units('AU')
                                                velocity[ind][pos_update_ind] = -1*velocity[ind][pos_update_ind]# + scale_l.in_units('AU')
                                            if sep_value > 20000:
                                                import pdb
                                                pdb.set_trace()
                                        
                                        
                                        #Calculation of semimajor axis:
                                        mass_sys = yt.YTArray(res['mass'][np.array([ind_1,ind_2])], 'msun')
                                        velocity = yt.YTArray(res['absvel'][np.array([ind_1,ind_2])].T, 'km/s')
                                        CoM_pos = np.sum(position*mass_sys, axis=1)/np.sum(mass_sys)
                                        CoM_vel = np.sum(velocity*mass_sys, axis=1)/np.sum(mass_sys)
                                        pos_rel_to_com = (position.T - CoM_pos).T
                                        distance_from_com = np.sqrt(np.sum(pos_rel_to_com**2, axis=0))
                                        vel_rel_to_com = (velocity.T - CoM_vel).T
                                        separation = np.sum(distance_from_com)
                                        
                                        relative_speed_to_com = np.sqrt(np.sum(vel_rel_to_com**2, axis=0))
                                        reduced_mass = (np.product(mass_sys)/np.sum(mass_sys)).in_units('g')
                                        E_pot = (-1*(yt.units.G*np.product(mass_sys.in_units('g')))/separation.in_units('cm')).in_units('erg')
                                        E_kin = np.sum((0.5*mass_sys.in_units('g')*relative_speed_to_com.in_units('cm/s')**2).in_units('erg'))
                                        epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
                                        m_x_r = yt.YTArray(np.cross(pos_rel_to_com.T.in_units('cm'), vel_rel_to_com.T.in_units('cm/s')).T, 'cm**2/s')
                                        Ang_mom = mass_sys.in_units('g').T*m_x_r
                                        Ang_mom_tot = np.sqrt(np.sum(np.sum(Ang_mom, axis=1)**2, axis=0))
                                        h_val = Ang_mom_tot/reduced_mass.in_units('g')
                                        ecc = np.sqrt(1 + (2.*epsilon*h_val**2.)/((yt.units.G*np.sum(mass_sys.in_units('g')))**2.))
                                        semimajor_a = ((h_val**2)/(yt.units.G*np.sum(mass_sys.in_units('g'))*(1-ecc**2))).in_units('AU')
                                        semimajor_arr.append(semimajor_a.value.tolist())
                                        
                                        sep_arr.append(sep_value)
                                        if args.replace_index_type == "M":
                                            replace_ind = np.argmax(mass[sub_sys])
                                        elif args.replace_index_type == "L":
                                            replace_ind = np.argmax(L_tot[sub_sys])
                                        replace_string = str(sub_sys[replace_ind])
                                        sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                                        if '[' not in sys_comps:
                                            reduced = True
                                        elif '[' in sys_comps:
                                            if len(eval(sys_comps)) == 1:
                                                reduced = True
                                        break
                                    else:
                                        import pdb
                                        pdb.set_trace()

                        '''
                        reduced = False
                        sep_arr = []
                        while reduced == False:
                            bracket_pos = []
                            for char_it in range(len(sys_comps)):
                                if sys_comps[char_it] == '[':
                                    bracket_pos.append(char_it)
                                elif sys_comps[char_it] == ']':
                                    open_ind = bracket_pos.pop()
                                    try:
                                        sub_sys_comps = eval(sys_comps[open_ind:char_it+1])
                                    except:
                                        comma_split = sys_comps[open_ind:char_it+1].split(',')
                                        sub_sys_comps = eval(comma_split[0]+comma_split[1])
                                    binary_ind = np.where((res['index1']==sub_sys_comps[0])&(res['index2']==sub_sys_comps[1]))[0][0]
                                    ind_1 = res['index1'][binary_ind]
                                    ind_2 = res['index2'][binary_ind]
                                    pos_diff = res['abspos'][ind_1] - res['abspos'][ind_2]
                                    sep_value = np.sqrt(np.sum(pos_diff**2))
                                    if sep_value > 10000.:
                                        update_inds = np.where(pos_diff>scale_l.in_units('AU')/2)[0]
                                        pos_diff[update_inds] = pos_diff[update_inds] - scale_l.in_units('AU').value
                                        sep_value = np.sqrt(np.sum(pos_diff**2))
                                    sep_arr.append(sep_value)
                                    replace_string = str(binary_ind)
                                    str_1 = sys_comps[:open_ind]
                                    str_2 = sys_comps[char_it+1:]
                                    sys_comps = str_1 + replace_string + str_2
                                    if '[' not in sys_comps:
                                        reduced = True
                                    break
                            '''
                        if sys_comps_str not in System_times.keys():
                            System_times.update({sys_comps_str:[time_yr]})
                            System_seps.update({sys_comps_str:[sep_arr]})
                            System_semimajor.update({sys_comps_str:[semimajor_arr]})
                        else:
                            System_times[sys_comps_str].append(time_yr)
                            System_seps[sys_comps_str].append(sep_arr)
                            System_semimajor[sys_comps_str].append(semimajor_arr)
                for sys_key in System_times.keys():
                    if sys_key not in updated_systems:
                        System_times[sys_key].append(time_yr)
                        System_seps[sys_key].append(np.ones(np.shape(System_seps[sys_key][-1]))*np.nan)
                        System_semimajor[sys_key].append(np.ones(np.shape(System_semimajor[sys_key][-1]))*np.nan)
                N_multi_stars.append(n_multi)
                M_tot_multi.append(M_multi)
                
                pickle_file_rank = pickle_file.split('.pkl')[0] + "_" + ("%02d" % rank) + ".pkl"
                file = open(pickle_file_rank, 'wb')
                pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times),file)
                file.close()
                
                print("time_it", time_it, "of", len(global_data['time'].T[0]), "Updated pickle file:", pickle_file.split('.pkl')[0] + "_" +str(rank) + ".pkl")
            else:
                print("not enough stars for multi analysis")
    print("finished going through times on rank", rank)
            
sys.stdout.flush()
CW.Barrier()


#compile all the pickles
if rank == 0:
    try:
        try:
            file = open(pickle_file+'_with_means.pkl', 'rb')
            Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_mean_times, System_mean_seps, System_lifetimes, Initial_Seps, Final_seps, Sep_maxs, Sep_mins = pickle.load(file)
            file.close()
        except:
            file = open(pickle_file+'.pkl', 'rb')
            Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times = pickle.load(file)
            file.close()
    except:
        pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
        if len(pickle_files) > 0:
            Times_full = []
            SFE_full = []
            SFE_n_full = []
            M_tot_full = []
            M_tot_multi_full = []
            N_stars_full = []
            N_vis_stars_full = []
            N_multi_stars_full = []
            System_seps_full = {}
            System_semimajor_full = {}
            System_times_full = {}
            for pick_file in pickle_files:
                file = open(pick_file, 'rb')
                Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times = pickle.load(file)
                file.close()
                for time_key in System_times.keys():
                    if time_key not in System_times_full.keys():
                        System_times_full.update({time_key:System_times[time_key]})
                        System_seps_full.update({time_key:System_seps[time_key]})
                        System_semimajor_full.update({time_key:System_semimajor[time_key]})
                    else:
                        System_times_full[time_key] = System_times_full[time_key] + System_times[time_key]
                        System_seps_full[time_key] = System_seps_full[time_key] + System_seps[time_key]
                        System_semimajor_full[time_key] = System_semimajor_full[time_key] + System_semimajor[time_key]
                Times_full = Times_full + Times
                SFE_full = SFE_full + SFE
                SFE_n_full = SFE_n_full + SFE_n
                M_tot_full = M_tot_full + M_tot
                M_tot_multi_full = M_tot_multi_full + M_tot_multi
                N_stars_full = N_stars_full + N_stars
                N_vis_stars_full = N_vis_stars_full + N_vis_stars
                N_multi_stars_full = N_multi_stars_full + N_multi_stars
                os.remove(pick_file)
            
            #Let's sort the data
            for time_key in System_times_full:
                sorted_array = np.array(System_times_full[time_key])
                dict_sort = np.argsort(sorted_array)
                sorted_array = np.array(System_times_full[time_key])[dict_sort]
                System_times[time_key] = sorted_array.tolist()
                sorted_array = np.array(System_seps_full[time_key])[dict_sort]
                System_seps[time_key] = sorted_array.tolist()
                sorted_array = np.array(System_semimajor_full[time_key])[dict_sort]
                System_semimajor[time_key] = sorted_array.tolist()
            sorted_inds = np.argsort(Times_full)
            Times = np.array(Times_full)[sorted_inds].tolist()
            SFE = np.array(SFE_full)[sorted_inds].tolist()
            SFE_n = np.array(SFE_n_full)[sorted_inds].tolist()
            M_tot = np.array(M_tot_full)[sorted_inds].tolist()
            M_tot_multi = np.array(M_tot_multi_full)[sorted_inds].tolist()
            N_stars = np.array(N_stars_full)[sorted_inds].tolist()
            N_vis_stars = np.array(N_vis_stars_full)[sorted_inds].tolist()
            N_multi_stars = np.array(N_multi_stars_full)[sorted_inds].tolist()
            
            file = open(pickle_file+'.pkl', 'wb')
            pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times),file)
            file.close()
    
    #Create plot
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=6, figsize=(8, 8))
    plt.subplots_adjust(hspace=0.0)
    gs = axs[3].get_gridspec()
    for ax in axs[3:]:
        ax.remove()
    axbig = fig.add_subplot(gs[3:])
    
    axs[0].plot(Times, SFE, label='SFE')
    axs[0].plot(Times, SFE_n, label='SFE$_n$')
    axs[0].legend(loc='best')
    axs[0].set_ylabel('SFE')
    axs[0].set_ylim(bottom=0)
    axs[0].set_xlim([Times[0], Times[-1]])
    xticklabels = axs[0].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    
    axs[1].plot(Times, M_tot, label='Total Mass')
    axs[1].plot(Times, M_tot_multi, label='Mass in Multiple systems')
    axs[1].set_ylabel('M$_{tot}$ (M$_\odot$)')
    axs[1].set_ylim(bottom=0)
    axs[1].legend(loc='best')
    axs[1].set_xlim([Times[0], Times[-1]])
    xticklabels = axs[1].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    
    axs[2].plot(Times, N_stars, label='Total stars')
    axs[2].plot(Times, N_multi_stars, label='Stars in Multiple systems')
    axs[2].plot(Times, N_vis_stars, label='Visible stars')
    axs[2].set_ylabel('\# stars')
    axs[2].set_ylim(bottom=0)
    axs[2].legend(loc='best')
    axs[2].set_xlim([Times[0], Times[-1]])
    xticklabels = axs[2].get_xticklabels()
    plt.setp(xticklabels, visible=False)
    
    window = 50000#year
    S_bins = np.logspace(1,4,13)
    if ("System_mean_times" in locals()) == False:
        System_mean_times = {}
        System_mean_seps = {}
        System_lifetimes = []
        Sep_maxs = []
        Sep_mins = []
        Initial_Seps = []
        Final_seps = []
    lifetime_it = -1
    for time_key in System_times.keys():
        axbig.semilogy(System_times[time_key], System_seps[time_key], alpha=0.1, color='k')
        try:
            lifetime_it = lifetime_it + 1
            if System_lifetimes[lifetime_it] > args.sys_lifetime_threshold:
                axbig.semilogy(System_mean_times[time_key], System_mean_seps[time_key], color='k', linewidth=0.5)
        except:
            mean_time = []
            mean_sep = []
            for time in System_times[time_key]:
                start_time = time - window/2.
                end_time = time + window/2.
                start_ind = np.argmin(abs(np.array(System_times[time_key]) - start_time))
                end_ind = np.argmin(abs(np.array(System_times[time_key]) - end_time))
                if end_ind != start_ind:
                    non_nan_inds = np.argwhere(np.isnan(np.array(System_seps[time_key][start_ind:end_ind])[:,0])==False)
                    mean_t = np.mean(np.array(System_times[time_key][start_ind:end_ind])[non_nan_inds])
                    mean_s = np.mean(np.array(System_seps[time_key][start_ind:end_ind])[non_nan_inds], axis=0)
                    #mean_t = np.mean(np.array(System_times[time_key][start_ind:end_ind]))
                    #mean_s = np.mean(np.array(System_seps[time_key][start_ind:end_ind]), axis=0)
                    mean_time.append(mean_t)
                    #mean_sep.append(mean_s)
                    mean_sep.append(mean_s[0])
            non_nan_inds = np.argwhere(np.isnan(np.array(System_seps[time_key]).T[0])==False).T[0]
            sys_life_time = np.array(System_times[time_key])[non_nan_inds][-1] - np.array(System_times[time_key])[non_nan_inds][0]
            System_lifetimes.append(sys_life_time)
            if sys_life_time > args.sys_lifetime_threshold:
                end_sep_time = np.array(System_times[time_key])[non_nan_inds][0] + sys_life_time - 1000
                #start_sep_time = np.array(System_times[time_key])[non_nan_inds][0] + 1000
                end_time_it = np.argmin(abs(np.array(System_times[time_key])[non_nan_inds] - end_sep_time))
                #start_sep_it = np.argmin(abs(np.array(System_times[time_key])[non_nan_inds] - start_sep_time))
                for sys in range(np.shape(np.array(mean_sep).T)[0]):
                    sep_max = np.nanmax(np.array(System_seps[time_key]).T[sys])
                    sep_min = np.nanmin(np.array(System_seps[time_key]).T[sys])
                    #sep_init = np.mean(np.array(System_seps[time_key]).T[sys][:start_sep_it])
                    sep_init = np.array(System_seps[time_key]).T[sys][0]
                    sep_final = np.mean(np.array(System_seps[time_key]).T[sys][non_nan_inds][end_time_it:])
                    Sep_maxs.append(sep_max)
                    Sep_mins.append(sep_min)
                    Initial_Seps.append(sep_init)
                    Final_seps.append(sep_final)
            System_mean_times.update({time_key:mean_time})
            System_mean_seps.update({time_key:np.array(mean_sep)})
            if sys_life_time > args.sys_lifetime_threshold:
                axbig.semilogy(System_mean_times[time_key], System_mean_seps[time_key], color='k', linewidth=0.5)
        
        print("plotted system:", time_key)
    for bin_bound in S_bins:
        axbig.axhline(y=bin_bound, linewidth=0.5)
    axbig.set_ylabel('Separation (AU)')
    axbig.set_xlabel('Time ($yr$)')
    axbig.set_ylim([10, 10000])
    axbig.set_xlim([Times[0], Times[-1]])
    
    file = open(pickle_file+'_with_means.pkl', 'wb')
    pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_mean_times, System_mean_seps, System_lifetimes, Initial_Seps, Final_seps, Sep_maxs, Sep_mins),file)
    file.close()
    
    plt.savefig('superplot.jpg', format='jpg', bbox_inches='tight')
    plt.savefig('superplot.pdf', format='pdf', bbox_inches='tight')
    print('Created superplot.jpg')
    
    d_sep = np.array(Initial_Seps) - np.array(Final_seps)
    shrinkage = 100*(d_sep/np.array(Initial_Seps))
    bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    shrink_hist, bins = np.histogram(shrinkage, bins=bins)
    bin_centers = np.linspace(-5, 95, 11)
    plt.clf()
    plt.bar(bin_centers, shrink_hist, width=10, edgecolor='black', alpha=0.5, label="Simulation")
    plt.xlabel('Orbital Shrinkage (\% of inital separation)')
    plt.ylabel('Number')
    plt.ylim(bottom=0.0)
    plt.savefig('shrinkage.jpg', format='jpg', bbox_inches='tight')
    print('created shrinkage.jpg')
    
