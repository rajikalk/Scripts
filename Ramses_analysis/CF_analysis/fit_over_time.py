import numpy as np
import matplotlib.pyplot as plt
import yt
import glob
import sys
import collections
import matplotlib as mpl
import pickle
import os
from mpi4py.MPI import COMM_WORLD as CW
import matplotlib.gridspec as gridspec

#Define globals
f_acc= 0.5
Accretion_array = []

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-savedir", "--save_directory", help="Where do you want to save the figures?", type=str, default="./")
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=10000)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
#=====================================================================================================

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

#=====================================================================================================
file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Per = pickle.load(file_open)
file_open.close()

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Orion_data.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Ori = pickle.load(file_open)
file_open.close()

#=====================================================================================================

pickle_file = args.pickled_file
try:
    file = open(pickle_file, 'rb')
    Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
    file.close()
except:
    pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
    Times_full = []
    CF_per_bin_full = []
    n_systems_full = []
    for pick_file in pickle_files:
        file = open(pick_file, 'rb')
        CF_arrays, N_sys_total, Times, Sink_Luminosities, Sink_Accretion = pickle.load(file)
        file.close()
        Sink_Luminosities_full = {}
        Sink_Accretion_full = {}
        for L_key in Sink_Luminosities.keys():
            if L_key not in Sink_Luminosities_full.keys():
                Sink_Luminosities_full.update({L_key:Sink_Luminosities[L_key]})
                Sink_Accretion_full.update({L_key:Sink_Accretion[L_key]})
            else:
                Sink_Luminosities_full[L_key] = Sink_Luminosities_full[L_key] + Sink_Luminosities[L_key]
                Sink_Accretion_full[L_key] = Sink_Accretion_full[L_key] + Sink_Accretion[L_key]
        Times_full = Times_full + Times
        CF_per_bin_full = CF_per_bin_full + CF_arrays
        n_systems_full = n_systems_full + N_sys_total
        os.remove(pick_file)
    
    #Let's sort the data
    for L_key in Sink_Luminosities_full:
        sorted_array = np.array(Sink_Luminosities_full[L_key])
        dict_sort = np.argsort(sorted_array.T[0])
        sorted_array.T[0] = np.array(Sink_Luminosities_full[L_key]).T[0][dict_sort]
        sorted_array.T[1] = np.array(Sink_Luminosities_full[L_key]).T[1][dict_sort]
        Sink_Luminosities[L_key] = sorted_array.tolist()
        sorted_array.T[0] = np.array(Sink_Accretion_full[L_key]).T[0][dict_sort]
        sorted_array.T[1] = np.array(Sink_Accretion_full[L_key]).T[1][dict_sort]
        Sink_Accretion[L_key] = sorted_array.tolist()
    sorted_inds = np.argsort(Times_full)
    Times = np.array(Times_full)[sorted_inds]
    CF_Array_Full = np.array(CF_per_bin_full)[sorted_inds]
    N_sys_total = np.array(n_systems_full)[sorted_inds]
    
    file = open(pickle_file+'.pkl', 'wb')
    pickle.dump((Times, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion),file)
    file.close()

#do chi squared
reduced_chi_square_tobin = []
usable_bin_inds = np.array([ 2,  4,  5,  6,  7,  8,  9, 10, 11, 12])
Times = (Times-Times[0])/1e6
Non_bimodal_dist = np.array([0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 ,
       0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 ,
       0.05457945, 0.09212946, 0.12967947])
Bimodal_dist = np.array([0.00780592, 0.01872544, 0.03498372, 0.05090101, 0.05767853,
       0.05090462, 0.03504812, 0.01944561, 0.01284564, 0.02460516,
       0.06113053, 0.10387658, 0.11138321])
CF_it = -1
for CF_hist in CF_Array_Full:
    CF_it = CF_it + 1
    chi_red_tobin = (np.sum(((CF_hist[usable_bin_inds]-CF_per_bin_Tobin_Per[usable_bin_inds])**2)/(CF_errs[usable_bin_inds]**2)))/len(CF_hist[usable_bin_inds])
    reduced_chi_square_tobin.append(chi_red_tobin)
    if chi_red_tobin < 0.5:
        plt.clf()
        plt.bar(bin_centers, CF_hist, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
        plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
        plt.legend(loc='best')
        plt.xlabel('Log Separation (AU)')
        plt.ylabel('Companion Frequency')
        plt.xlim([1, 4])
        plt.ylim(bottom=0.0)
        plt.title("reduced chi$^2$ = " + str(chi_red_tobin))
        plt.savefig(args.save_directory+'Time_'+str(Times[CF_it])+'Myr_chi_tobin.png', format='png', bbox_inches='tight')
        print('saved', args.save_directory+'Time_'+str(Times[CF_it])+'Myr_chi_tobin.png')
    
plt.clf()
plt.semilogy(SFE, reduced_chi_square_tobin)
plt.xlabel("SFE")
plt.ylabel("mean Chi squared (<$\chi^2$>)")
plt.xlim([Times[0], Times[-1]])
plt.ylim([np.min(reduced_chi_square_tobin), np.max(reduced_chi_square_tobin)])
plt.savefig(args.save_directory + "reduced_chi_squared_tobin.png")
