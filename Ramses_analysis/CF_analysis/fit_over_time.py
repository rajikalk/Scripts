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
import matplotlib
from scipy import stats

matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

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
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=None)
    parser.add_argument("-SFE_spread", "--SFE_spread_val", help="how much SFE around the central SFE do you want to intergrate over?", type=float, default=0.001)
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

S_bins = np.logspace(0.75,4,14)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_CDF.pkl", "rb")
Perseus_sep = pickle.load(file_open)
file_open.close()

Perseus_log_sep = np.log10(np.sort(Perseus_sep))
Perseus_frequency = np.cumsum(np.sort(Perseus_sep))/np.cumsum(np.sort(Perseus_sep))[-1]
#Perseus_log_sep = bin_centers[2:]
Perseus_frequency_CF = np.cumsum(CF_per_bin_Tobin_Per[2:])/np.cumsum(CF_per_bin_Tobin_Per[2:])[-1]

#=====================================================================================================

pickle_file = args.pickled_file
try:
    file = open(pickle_file, 'rb')
    try:
        Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion, All_separations = pickle.load(file)
        file.close()
    except:
        file.close()
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
reduced_chi_square_selected = []
reduced_chi_square_sim = []
KS_test_sep = []
D_crits_sep = []
KS_test_CF = []
D_crits_CF = []
selected_bins = np.array([2, 4, 6, 7, 8, 9, 10, 11, 12])
tobin_bins = np.argwhere(CF_per_bin_Tobin_Per>0).T[0]
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

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

SFE_1_ind = np.argmin(abs(SFE-0.01))
CF_Array_Full = CF_Array_Full[SFE_1_ind:]
alpha = 0.99

for CF_it in range(len(CF_Array_Full)):
    time_val = Times[CF_it]
    if args.time_spread != None:
        start_time = Times[CF_it] - args.time_spread/2.
        end_time = Times[CF_it] + args.time_spread/2.
    
        start_integration_it = np.argmin(abs(Times - start_time))
        end_integration_it = np.argmin(abs(Times - end_time))
    else:
        SFE_val = SFE[CF_it]
        start_SFE = SFE[CF_it] - args.SFE_spread_val
        end_SFE = SFE[CF_it] + args.SFE_spread_val
        
        start_integration_it = np.argmin(abs(SFE - start_SFE))
        end_integration_it = np.argmin(abs(SFE - end_SFE))
    
    CF_median = np.median(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]

    CF_hist = CF_median#CF_Array_Full[CF_it]
    curr_errs = []
    for cf_val in range(len(CF_hist)):
        if CF_hist[cf_val] < CF_per_bin_Tobin_Per[cf_val]:
            curr_errs.append(CF_err[1][cf_val])
        else:
            curr_errs.append(CF_err[0][cf_val])
            
    Tobin_errs = []
    for cf_val in range(len(CF_Array_Full[CF_it])):
        if CF_per_bin_Tobin_Per[cf_val] < CF_Array_Full[CF_it][cf_val]:
            Tobin_errs.append(CF_errs_Per[1][cf_val])
        else:
            Tobin_errs.append(CF_errs_Per[0][cf_val])
    
    #N_sys = np.sum(N_sys_total[CF_it],axis=1)
    #CF_errs = CF_std#np.mean(CF_errs_Per,axis=0)
    #chi_red_tobin = (np.median(((CF_hist[1:]-gauss_total)**2)/(np.array(curr_errs)**2)))/len(CF_hist[1:])
    chi_red_calc = (np.mean(((CF_hist[selected_bins]-CF_per_bin_Tobin_Per[selected_bins])**2)/(np.array(curr_errs)[selected_bins]**2)))
    chi_red_selected = (np.mean(((CF_hist[selected_bins]-CF_per_bin_Tobin_Per[selected_bins])**2)/(np.array(Tobin_errs)[selected_bins]**2)))
    chi_red_tobin = (np.mean(((CF_hist[tobin_bins]-CF_per_bin_Tobin_Per[tobin_bins])**2)/(np.array(Tobin_errs)[tobin_bins]**2)))
    reduced_chi_square_sim.append(chi_red_calc)
    reduced_chi_square_tobin.append(chi_red_tobin)
    reduced_chi_square_selected.append(chi_red_selected)
    
    usable_seps = np.where(All_separations[CF_it] >= 10**1.25)[0]
    if len(usable_seps) > 0:
        log_sep = np.log10(np.sort(All_separations[CF_it][usable_seps]))
        frequency = np.cumsum(np.sort(All_separations[CF_it][usable_seps]))/np.cumsum(np.sort(All_separations[CF_it][usable_seps]))[-1]
        KS_test_result = stats.ks_2samp(Perseus_frequency, frequency)[0]
        m = len(usable_seps)
        n = len(Perseus_frequency)
        D_crit = np.sqrt((-1*np.log(alpha/2))*((1+(m/n))/(2*m)))
        
        frequency = np.cumsum(CF_hist[2:])/np.cumsum(CF_hist[2:])[-1]
        KS_test_result_CF = np.max(abs(Perseus_frequency_CF - frequency))
        m = len(CF_hist[2:])
        n = len(CF_hist[2:])
        D_crit_CF = np.sqrt((-1*np.log(alpha/2))*((1+(m/n))/(2*m)))
        #log_sep = bin_centers[2:]
        #KS_test_result = np.max(abs(Perseus_frequency - frequency))
    else:
        KS_test_result = np.nan
        D_crit = np.nan
        KS_test_result_CF = np.nan
        D_crit_CF = np.nan
    KS_test_sep.append(KS_test_result)
    D_crits_sep.append(D_crit)
    KS_test_CF.append(KS_test_result_CF)
    D_crits_CF.append(D_crit_CF)
    if chi_red_tobin < 0.01:
        plt.clf()
        plt.bar(bin_centers, CF_hist, yerr=curr_errs, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
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
plt.figure(figsize=(single_col_width,0.7*single_col_width))
plt.semilogy(SFE[SFE_1_ind:], reduced_chi_square_tobin)
plt.xlabel("SFE")
plt.ylabel("Fit (<$\chi^2$>)", labelpad=-0.2)
plt.xlim([0.01, 0.05])
#plt.ylim(top=1000)
#plt.ylim([np.argmin(reduced_chi_square_tobin), np.argmax(reduced_chi_square_tobin)])
plt.savefig(args.save_directory + "reduced_chi_squared_tobin.pdf", format='pdf', bbox_inches='tight', pad_inches = 0.02)

plt.clf()
plt.figure(figsize=(single_col_width,0.7*single_col_width))
plt.semilogy(SFE[SFE_1_ind:], KS_test_sep, label="KS statistic")
plt.semilogy(SFE[SFE_1_ind:], D_crits_sep, label="Critical Value ($\alpha="+str(alpha)+"$)")
plt.xlabel("SFE")
plt.ylabel("KS statistic", labelpad=-0.2)
plt.xlim([0.01, 0.05])
#plt.ylim(top=1000)
#plt.ylim([np.argmin(reduced_chi_square_tobin), np.argmax(reduced_chi_square_tobin)])
plt.savefig(args.save_directory + "KS_test_alpha_"+str(alpha)+".pdf", format='pdf', bbox_inches='tight', pad_inches = 0.02)

file = open(args.save_directory + 'chi_squared_fit.pkl', 'wb')
pickle.dump((SFE[SFE_1_ind:], reduced_chi_square_selected, reduced_chi_square_tobin, KS_test_sep, D_crits_sep, KS_test_CF, D_crits_CF), file)
file.close()
