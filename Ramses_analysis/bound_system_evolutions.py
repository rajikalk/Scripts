import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
import sys
import collections
import matplotlib as mpl
import pickle
import os
from mpi4py.MPI import COMM_WORLD as CW
import matplotlib.gridspec as gridspec


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-preffix", "--figure_prefix", help="Do you want to give saved figures a preffix?", default="", type=str)
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/stars_red_512.pkl', type=str)
    parser.add_argument("-time_window", "--time_intergration_window", help="Over what time do you want to intergrate to calculate accretion rate?", default=1000.0, type=float)
    parser.add_argument("-verbose", "--verbose_printing", help="Would you like to print debug lines?", type=str, default='False')
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-plot_only", "--make_plots_only", help="Do you just want to make plots? Not calculate the CF", type=str, default='False')
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=0.01)
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=35.6)
    parser.add_argument("-start_ind", "--start_index", help="Do you want to start at a particular index?", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
def losi(i, res):
    if (res['n'][i]==1):
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
    
def accretion(global_data, sink_inds, max_time_ind, min_time_ind):
    """
    Calculates the accretion of the given indeexes
    """
    dM = (global_data['m'][max_time_ind,sink_inds] - global_data['m'][min_time_ind,sink_inds])*units['mass_unit'].in_units('msun')
    dt = (global_data['time'][max_time_ind,sink_inds] - global_data['time'][min_time_ind,sink_inds])*units['time_unit'].in_units('yr')
    M_dot = dM/dt
    return M_dot
    
def luminosity(global_data, sink_inds, max_time_ind, min_time_ind):
    """
    Calculates the luminosity of the given indexes
    """
    global f_acc
    M_dot = accretion(global_data, sink_inds, max_time_ind, min_time_ind)
    M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    return L_tot
#=====================================================================================================
args = parse_inputs()
f_acc = 0.5
window = yt.YTQuantity(args.time_intergration_window, 'yr')
window = yt.YTQuantity(args.time_intergration_window, 'yr')
luminosity_lower_limit = 0.01
luminosity_upper_limit = args.upper_L_limit
accretion_limit = args.accretion_limit

rank = CW.Get_rank()
size = CW.Get_size()

datadir = sys.argv[1]
savedir = sys.argv[2]

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

#How will I save the data?
Times_Total = {}
Separations_Total = {}
Times = {}
Separations = {}
Visible = {}

sys.stdout.flush()
CW.Barrier()

'''
#low resolution data
datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data_256'
nout = 133
'''

#high resolution data
#datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data'

#Load global pickle data
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
print('Loaded global pickle data')

SFE_ind = np.argmin(np.abs(0.049-(np.sum(global_data['m'], axis=1)*units['mass_unit'].value)/units['mass_unit'].value))
SFE_t_ff = global_data['time'].T[0][SFE_ind]
dt = args.time_spread
time_bounds = [SFE_t_ff-dt,SFE_t_ff+dt]

if args.start_index == None:
    start_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[0]))
else:
    start_time_ind = args.start_index
end_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[1]))

radius = yt.YTQuantity(2.0, 'rsun')
temperature = yt.YTQuantity(3000, 'K')
    
sys.stdout.flush()
CW.Barrier()

if args.make_plots_only == 'False':
    rit = -1
    for time_it in range(start_time_ind, end_time_ind+1):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            print("doing", time_it, "of", end_time_ind)
            non_zero_inds = np.where(global_data['m'][time_it]>0)[0]
            abspos = np.array([global_data['x'][time_it][non_zero_inds], global_data['y'][time_it][non_zero_inds], global_data['z'][time_it][non_zero_inds]]).T#*scale_l
            absvel = np.array([global_data['ux'][time_it][non_zero_inds], global_data['uy'][time_it][non_zero_inds], global_data['uz'][time_it][non_zero_inds]]).T#*scale_v
            mass = np.array(global_data['m'][time_it][non_zero_inds])
            time = global_data['time'][time_it][non_zero_inds][0]

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

            time_yt = yt.YTArray(time*scale_t, 's')
            
            min_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')+window))
            max_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')-window))
            global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')))
            
            res = m.multipleAnalysis(S,cutoff=10000,nmax=6, cyclic=False)
            multi_systems = np.where((res['topSystem']==True) & (res['n']>1))[0]
            
            sink_inds = np.where((res['n']==1))[0]
            sink_inds_total = np.arange(len(res['n']))
            nan_size = len(sink_inds_total) - len(sink_inds)
            nan_array = yt.YTArray(np.ones(nan_size)*np.nan, 'Lsun')
            #sink_inds = np.arange(len(res['n']))
            L_tot = luminosity(global_data, sink_inds, max_time_ind, min_time_ind)
            M_dot = accretion(global_data, sink_inds, max_time_ind, min_time_ind)
            L_tot = np.append(L_tot, nan_array)
            M_dot = np.append(M_dot, nan_array)
            s_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<luminosity_upper_limit))[0]
            visible_stars = sink_inds[s_inds]
            
            for multi_sys in multi_systems:
                sys_comps = losi(multi_sys, res)
                #sys_comps = sorted(flatten(sys_comps))
                sep_array = []
                
                reduced = False
                hierarchical_bool = is_hierarchical(sys_comps)
                sys_string = str(sys_comps)
                if hierarchical_bool == True:
                    sys_inds = sorted(flatten(sys_comps))
                    central_ind = sys_inds[np.argmax(L_tot[sys_inds])]
                    central_pos = res['abspos'][central_ind]
                    positions = res['abspos'][sys_inds]
                    separations = np.sqrt(np.sum((positions - central_pos)**2, axis=1))
                    for sep in separations:
                        if sep > 0:
                            sep_array.append(sep)
                else:
                    reduced = False
                    sys_string = str(sys_comps)
                    while reduced == False:
                        bracket_pos = []
                        for char_it in range(len(sys_string)):
                            if sys_string[char_it] == '[':
                                bracket_pos.append(char_it)
                            elif sys_string[char_it] == ']':
                                open_ind = bracket_pos.pop()
                                sub_sys_comps = eval(sys_string[open_ind:char_it+1])
                                binary_ind = np.where((res['index1']==sub_sys_comps[0])&(res['index2']==sub_sys_comps[1]))[0][0]
                                ind_1 = res['index1'][binary_ind]
                                ind_2 = res['index2'][binary_ind]
                                sep_value = np.sqrt(np.sum((res['abspos'][ind_1] - res['abspos'][ind_2])**2))
                                sep_array.append(sep_value)
                                str_1 = sys_string[:open_ind]
                                str_2 = sys_string[char_it+1:]
                                sys_string = str_1 + str(binary_ind) + str_2
                                if '[' not in sys_string:
                                    reduced = True
                                break
                '''
                reduced = False
                sys_string = str(sys_comps)
                while reduced == False:
                    bracket_pos = []
                    for char_it in range(len(sys_string)):
                        if sys_string[char_it] == '[':
                            bracket_pos.append(char_it)
                        elif sys_string[char_it] == ']':
                            open_ind = bracket_pos.pop()
                            sub_sys_comps = eval(sys_string[open_ind:char_it+1])
                            binary_ind = np.where((res['index1']==sub_sys_comps[0])&(res['index2']==sub_sys_comps[1]))[0][0]
                            ind_1 = res['index1'][binary_ind]
                            ind_2 = res['index2'][binary_ind]
                            if ind_1 in visible_stars and ind_2 in visible_stars:
                                sep_value = np.sqrt(np.sum((res['abspos'][ind_1] - res['abspos'][ind_2])**2))
                                sep_array.append(sep_value)
                                visible_stars = np.append(visible_stars, binary_ind)
                                if np.isnan(L_tot[ind_1]) or np.isnan(L_tot[ind_2]):
                                    import pdb
                                    pdb.set_trace()
                                
                                L_tot[binary_ind] = L_tot[ind_1] + L_tot[ind_2]
                            elif ind_1 not in visible_stars and ind_2 not in visible_stars:
                                sep_array.append(np.nan)
                            elif ind_1 in visible_stars and ind_2 not in visible_stars:
                                res['abspos'][binary_ind] = res['abspos'][ind_1]
                                sep_array.append(np.nan)
                                visible_stars = np.append(visible_stars, binary_ind)
                                #res['abspos'][binary_ind] = res['abspos'][ind_1]
                                L_tot[binary_ind] =  L_tot[ind_1]
                            else:
                                res['abspos'][binary_ind] = res['abspos'][ind_2]
                                sep_array.append(np.nan)
                                visible_stars = np.append(visible_stars, binary_ind)
                                #res['abspos'][binary_ind] = res['abspos'][ind_2]
                                L_tot[binary_ind] =  L_tot[ind_2]
                            str_1 = sys_string[:open_ind]
                            str_2 = sys_string[char_it+1:]
                            sys_string = str_1 + str(binary_ind) + str_2
                            if '[' not in sys_string:
                                reduced = True
                            break
                '''
                sys_comps = sorted(flatten(sys_comps))
                sys_string = str(sys_comps)
                if sys_string not in Times.keys():
                    Times.update({sys_string: [float(time_yt.in_units('yr').value)]})
                    Separations.update({sys_string: [sep_array]})
                    print("Added system:", sys_string)
                else:
                    Times[sys_string].append(float(time_yt.in_units('yr').value))
                    Separations[sys_string].append(sep_array)
    
    print("Finished going over times on rank", rank)
    #Send Data to Root
    
    file = open("bound_system_data_"+str(rank)+".pkl", 'wb')
    pickle.dump((Times, Separations),file)
    print('updated pickle', "bound_system_data_"+str(rank)+".pkl")
    file.close()

    
    sys.stdout.flush()
    CW.Barrier()
    '''
    if rank > 0:
        print("Sending data to root from rank", rank)
        #data = np.array([Times,Separations])
        data = {'Times': Times, 'Separations': Separations}
        CW.send(data, dest=0, tag=rank)
        print("Sent data to root from rank", rank)
    
    sys.stdout.flush()
    CW.Barrier()
    '''
    if rank == 0:
        pickle_files = sorted(glob.glob("bound_system_data_*.pkl"))
        for pickle_file in pickle_files:
            file = open(pickle_file, 'rb')
            Times, Separations = pickle.load(file)
            file.close()
    
            for key in Times.keys():
                if key not in Times_Total.keys():
                    Times_Total.update({key: Times[key]})
                    Separations_Total.update({key: Separations[key]})
                else:
                    Times_Total[key] = Times_Total[key] + Times[key]
                    Separations_Total[key] = Separations_Total[key] + Separations[key]
            os.remove(pickle_file)
        '''
        for rit_recv in range(1,size):
            data = CW.Recv(source=rit_recv, tag=rit_recv)
            for key in data[0].keys():
                if key not in Times_Total.keys():
                    Times_Total.update({key: data[0][key]})
                    Separations_Total.update({key: data[1][key]})
                else:
                    Times_Total[key] = Times_Total[key] + data[0][key]
                    Separations_Total[key] = Separations_Total[key] + data[1][key]
        
        for key in data['Times'].keys():
            data = CW.recv(source=rit_recv, tag=rit_recv)
            if key not in Times_Total.keys():
                Times_Total.update({key: data['Times'][key]})
                Separations_Total.update({key: data['Separations'][key]})
            else:
                Times_Total[key] = Times_Total[key] + data['Times'][key]
                Separations_Total[key] = Separations_Total[key] + data['Separations'][key]
        '''
        file = open("bound_system_data.pkl", 'wb')
        pickle.dump((Times_Total, Separations_Total),file)
        print('updated pickle', "bound_system_data.pkl")
        file.close()

    sys.stdout.flush()
    CW.Barrier()
    
if rank == 0:
    file = open("bound_system_data.pkl", 'rb')
    Times_Total, Separations_Total = pickle.load(file)
    file.close()

    for key in Times_Total.keys():
        sorted_inds = np.argsort(Times_Total[key])
        Times_Total[key] = np.array(Times_Total[key])[sorted_inds]
        Separations_Total[key] = np.array(Separations_Total[key])[sorted_inds]

    file = open("bound_system_data_sorted.pkl", 'wb')
    pickle.dump((Times_Total, Separations_Total),file)
    print('updated pickle', "bound_system_data_sorted.pkl")
    file.close()

    print("Sorted data, Now time to plot")

    plt.clf()
    window = 50
    S_bins = np.logspace(0.75,4,14)
    n_comp = []
    for key in Times_Total.keys():
        #n_comp.append(len(eval(key)))
        
        if len(Times_Total[key])>2*window:
            n_comp.append(len(eval(key)))
        
    n_comp = np.array(n_comp)
    N = len(np.where(n_comp>2)[0])
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, N)
    c_num = np.arange(N)
    np.random.shuffle(c_num)
    cit = 0
    for key in Times_Total.keys():
        '''
        if len(eval(key)) == 2:
            color = 'grey'
        else:
            color = cmap(c_num[cit])
            cit = cit + 1
            
        for systems_it in range(np.shape(Separations_Total[key])[1]):
            plt.semilogy(Times_Total[key], Separations_Total[key][:, systems_it], color=color)
        '''
        if len(Times_Total[key])>2*window:
            if len(eval(key)) == 2:
                color = 'grey'
            else:
                color = cmap(c_num[cit])
                cit = cit + 1

            moving_time = []
            moving_sep = []
            moving_upper_sep = []
            moving_lower_sep = []
            for smooth_ind in range(window, len(Times_Total[key])-window):
                ave_time = np.mean(Times_Total[key][smooth_ind-window:smooth_ind+window])
                seps = []
                upper_seps = []
                lower_seps = []
                for systems_it in range(np.shape(Separations_Total[key])[1]):
                    non_nan_inds = np.where(np.isnan(Separations_Total[key][:,systems_it][smooth_ind-window:smooth_ind+window])==False)[0]
                    if len(non_nan_inds) == 0:
                        ave_sep = np.nan
                        upper_sep = np.nan
                        lower_sep = np.nan
                    else:
                        ave_sep = np.mean(Separations_Total[key][:,systems_it][smooth_ind-window:smooth_ind+window][non_nan_inds])
                        upper_sep = np.max(Separations_Total[key][:,systems_it][smooth_ind-window:smooth_ind+window][non_nan_inds])
                        lower_sep = np.min(Separations_Total[key][:,systems_it][smooth_ind-window:smooth_ind+window][non_nan_inds])
                    seps.append(ave_sep)
                    upper_seps.append(upper_sep)
                    lower_seps.append(lower_sep)
                moving_time.append(ave_time)
                moving_sep.append(seps)
                moving_upper_sep.append(upper_seps)
                moving_lower_sep.append(lower_seps)
            
            moving_sep = np.array(moving_sep)
            moving_upper_sep = np.array(moving_upper_sep)
            moving_lower_sep = np.array(moving_lower_sep)
            #time to plot
            for systems_it in range(np.shape(Separations_Total[key])[1]):
                plt.semilogy(moving_time, moving_sep[:, systems_it], color=color)
                plt.fill_between(moving_time, moving_lower_sep[:, systems_it], moving_upper_sep[:, systems_it], alpha=0.05, color=color)
        
    for bin_bound in S_bins:
        plt.axhline(y=bin_bound, color='k')
    plt.axhline(y=8, color='k', ls='--')
    plt.xlabel("Time ($t_{ff}$)")
    plt.ylabel("Log Separation (AU)")
    plt.ylim(bottom=1)
    plt.savefig("Separation_evolution.jpg")
    plt.savefig("Separation_evolution.pdf")
    
    
    

