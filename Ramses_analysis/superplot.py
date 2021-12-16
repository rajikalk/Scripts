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
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-update", "--update_pickles", help="Do you want to remake the pickles?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif args.simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif args.simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.key():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

sys.stdout.flush()
CW.Barrier()

Times = []
SFE = []
SFE_n = []
M_tot = []
M_tot_multi = []
N_stars = []
N_multi_stars = []
System_seps = {}
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
print('loaded global data')

sys.stdout.flush()
CW.Barrier()

pickle_file = args.pickled_file

if rank == 0:
    try:
        file = open(pickle_file+'.pkl', 'rb')
        Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times = pickle.load(file)
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
            N_multi_stars_full = []
            System_seps_full = {}
            System_times_full = {}
            for pick_file in pickle_files:
                file = open(pick_file, 'rb')
                Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times = pickle.load(file)
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
            sorted_inds = np.argsort(Times_full)
            Times = np.array(Times_full)[sorted_inds].tolist()
            SFE = np.array(SFE_full)[sorted_inds].tolist()
            SFE_n = np.array(SFE_n_full)[sorted_inds].tolist()
            M_tot = np.array(M_tot_full)[sorted_inds].tolist()
            M_tot_multi = np.array(M_tot_multi_full)[sorted_inds].tolist()
            N_stars = np.array(N_stars_full)[sorted_inds].tolist()
            N_multi_stars = np.array(N_multi_stars_full)[sorted_inds].tolist()
            
            file = open(pickle_file+'.pkl', 'wb')
            pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times),file)
            file.close()
        else:
            print('not pickles exist')
            
sys.stdout.flush()
CW.Barrier()


if args.update_pickles == 'True':
    rit = -1
    for time_it in range(len(global_data['time'].T[0])):
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
                Times.append(time)
                
                S = pr.Sink()
                S._jet_factor = 1.
                S._scale_l = scale_l
                S._scale_v = scale_v
                S._scale_t = scale_t
                S._scale_d = scale_d
                S._time = yt.YTArray(time, '')
                S._abspos = yt.YTArray(abspos, '')
                S._absvel = yt.YTArray(absvel, '')
                S._mass = yt.YTArray(mass, '')
                
                N_stars.append(len(n_stars))
                M_tot_msun = np.sum(mass*units['mass_unit'].in_units('Msun'))
                M_tot.append(M_tot_msun)
                SFE_val = M_tot_msun/units['mass_unit'].in_units('Msun')
                SFE.append(SFE_val)
                SFE_n.append(SFE_val/time)
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=False, Grho=Grho)
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
                n_multi = np.sum(res['n'][multi_inds])
                M_multi = np.sum(res['mass'][multi_inds])
                for multi_ind in multi_inds:
                    sys_comps = losi(multi_ind, res)
                    sys_string = sorted(flatten(sys_comps))
                    sys_string = str(sys_string)
                    sys_comps = str(sys_comps)
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
                                sep_value = np.sqrt(np.sum((res['abspos'][ind_1] - res['abspos'][ind_2])**2))
                                if sep_value > 9990.0:
                                    print("Separation is", sep_value, "at time_it", time_it, "for multi_ind", multi_ind)
                                sep_arr.append(sep_value)
                                replace_string = str(binary_ind)
                                str_1 = sys_comps[:open_ind]
                                str_2 = sys_comps[char_it+1:]
                                sys_comps = str_1 + replace_string + str_2
                                if '[' not in sys_comps:
                                    reduced = True
                                break
                    if sys_string not in System_times.keys():
                        System_times.update({sys_string:[time]})
                        System_seps.update({sys_string:[sep_arr]})
                    else:
                        System_times[sys_string].append(time)
                        System_seps[sys_string].append(sep_arr)
                N_multi_stars.append(n_multi)
                M_tot_multi.append(M_multi)
                
                pickle_file_rank = pickle_file.split('.pkl')[0] + "_" +str(rank) + ".pkl"
                file = open(pickle_file_rank, 'wb')
                pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times),file)
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
        file = open(pickle_file+'.pkl', 'rb')
        Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times = pickle.load(file)
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
            N_multi_stars_full = []
            System_seps_full = {}
            System_times_full = {}
            for pick_file in pickle_files:
                file = open(pick_file, 'rb')
                Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times = pickle.load(file)
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
            sorted_inds = np.argsort(Times_full)
            Times = np.array(Times_full)[sorted_inds].tolist()
            SFE = np.array(SFE_full)[sorted_inds].tolist()
            SFE_n = np.array(SFE_n_full)[sorted_inds].tolist()
            M_tot = np.array(M_tot_full)[sorted_inds].tolist()
            M_tot_multi = np.array(M_tot_multi_full)[sorted_inds].tolist()
            N_stars = np.array(N_stars_full)[sorted_inds].tolist()
            N_multi_stars = np.array(N_multi_stars_full)[sorted_inds].tolist()
            
            file = open(pickle_file+'.pkl', 'wb')
            pickle.dump((Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_multi_stars, System_seps, System_times),file)
            file.close()
    
    #Create plot
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=6, figsize=(5, 8))
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
    
    axs[1].plot(Times, M_tot, label='Total Mass')
    axs[1].plot(Times, M_tot_multi, label='Mass in Multple systems')
    axs[1].set_ylabel('M$_{tot}$ (M$_\odot$)')
    axs[1].set_ylim(bottom=0)
    axs[1].legend(loc='best')
    axs[1].set_xlim([Times[0], Times[-1]])
    
    axs[2].plot(Times, N_stars, label='Total stars')
    axs[2].plot(Times, N_multi_stars, label='Stars in Multiple systems')
    axs[2].set_ylabel('\# stars')
    axs[2].set_ylim(bottom=0)
    axs[2].legend(loc='best')
    axs[2].set_xlim([Times[0], Times[-1]])
    
    for time_key in System_times.keys():
        axbig.semilogy(System_times[time_key], System_seps[time_key], alpha=0.25, color='k')
    S_bins = np.logspace(0.75,4,14)
    for bin_bound in S_bins:
        axbig.axhline(y=bin_bound)
    axbig.set_ylabel('Separation (AU)')
    axbig.set_xlabel('Time (t$_{ff}$)')
    axbig.set_ylim(top=10000)
    axbig.set_xlim([Times[0], Times[-1]])
    
    plt.savefig('superplot.png')
    print('Created superplot.png')
    
