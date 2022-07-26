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

def CF_err(CF, N_sys, z=1):
    Err_lower = (1/(1+(np.square(z)/N_sys)))*(CF + (np.square(z)/(2*N_sys))) - (z/(1+(np.square(z)/N_sys)))*np.sqrt(((CF*(1-CF))/N_sys)+(np.square(z)/(4*(N_sys**2))))
    Err_upper = (1/(1+(np.square(z)/N_sys)))*(CF + (np.square(z)/(2*N_sys))) + (z/(1+(np.square(z)/N_sys)))*np.sqrt(((CF*(1-CF))/N_sys)+(np.square(z)/(4*(N_sys**2))))
    Err_lower = CF - Err_lower
    Err_upper = Err_upper - CF
    return Err_lower, Err_upper

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
bin_centers, CF_per_bin_Tobin_Per, CF_errs_Per = pickle.load(file_open)
file_open.close()

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Orion_data.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Ori, CF_errs_Ori = pickle.load(file_open)
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
usable_bin_inds = np.array([2, 4, 6, 7, 8, 9, 10, 11, 12])

import scipy.stats as stats

#first peak
lower_amp = 0.11
lower_mean = np.log10(120)
lower_std = 0.4
upper_amp = 0.045
upper_mean = np.log10(3500)
upper_std = 0.3
x = np.linspace(1,4,100)
lower_gauss = lower_amp*stats.norm.pdf(bin_centers[1:], lower_mean, lower_std)
upper_gauss = upper_amp*stats.norm.pdf(bin_centers[1:], upper_mean, upper_std)
gauss_total = lower_gauss + upper_gauss


for CF_it in range(len(CF_Array_Full)):
    time_val = Times[CF_it]
    start_time = Times[CF_it] - args.time_spread/2.
    end_time = Times[CF_it] + args.time_spread/2.
    
    start_integration_it = np.argmin(abs(Times - start_time))
    end_integration_it = np.argmin(abs(Times - end_time))
    
    CF_median = np.median(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]

    CF_hist = CF_median#CF_Array_Full[CF_it]
    curr_errs = []
    for cf_val in range(len(CF_hist[1:])):
        if CF_hist[1:][cf_val] < gauss_total[cf_val]:
            curr_errs.append(CF_err[1][cf_val])
        else:
            curr_errs.append(CF_err[0][cf_val])
    #N_sys = np.sum(N_sys_total[CF_it],axis=1)
    #CF_errs = CF_std#np.mean(CF_errs_Per,axis=0)
    chi_red_tobin = (np.median(((CF_hist[1:]-gauss_total)**2)/(np.array(curr_errs)**2)))/len(CF_hist[1:])
    #chi_red_tobin = (np.sum(((CF_hist[usable_bin_inds]-CF_per_bin_Tobin_Per[usable_bin_inds])**2)/(CF_errs[usable_bin_inds]**2)))/len(CF_hist[usable_bin_inds])
    reduced_chi_square_tobin.append(chi_red_tobin)
    if chi_red_tobin < 0.05:
        plt.clf()
        plt.bar(bin_centers[1:], CF_hist[1:], yerr=curr_errs, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
        plt.bar(bin_centers, CF_per_bin_Tobin_Per, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
        plt.legend(loc='best')
        plt.xlabel('Log Separation (AU)')
        plt.ylabel('Companion Frequency')
        plt.xlim([1, 4])
        plt.ylim(bottom=0.0)
        plt.title("reduced chi$^2$ = " + str(chi_red_tobin))
        plt.savefig(args.save_directory+'SFE_'+str(SFE[CF_it])+'Myr_chi_tobin.png', format='png', bbox_inches='tight')
        print('saved', args.save_directory+'SFE_'+str(SFE[CF_it])+'Myr_chi_tobin.png')
    
plt.clf()
plt.semilogy(SFE, reduced_chi_square_tobin)
plt.xlabel("SFE")
plt.ylabel("mean Chi squared (<$\chi^2$>)")
plt.xlim([0, 0.05])
#plt.ylim([np.argmin(reduced_chi_square_tobin), np.argmax(reduced_chi_square_tobin)])
plt.savefig(args.save_directory + "reduced_chi_squared_tobin.png")
