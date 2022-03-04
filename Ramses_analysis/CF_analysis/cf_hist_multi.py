import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import os

#Tobin data
L_bins = np.logspace(-1.25,3.5,20)
S_bins = np.logspace(0.75,4,14)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
Tobin_luminosities_0_1 = [1.8, 0.50, 0.6, 0.7, 0.4, 0.36, 1.40, 0.80, 0.43, 1.20, 3.70, 1.70, 0.16, 0.05, 1.6, 0.54, 0.3, 1.20, 23.2, 0.16, 4.7, 16.80, 0.54, 0.09, 0.24, 1.80, 1.90, 3.20, 0.69, 1.50, 0.90, 1.30, 4.20, 3.60, 8.40, 0.68, 2.6, 1.80, 0.40, 0.70, 0.30, 0.60, 19.00, 5.30, 1.50, 0.30, 0.50, 0.54, 0.17, 0.32, 0.70, 2.80, 6.9, 1.10, 4.00, 7.00, 0.10, 32.50, 1.00, 8.3, 9.2, 1.4, 0.87, 1.50, 3.20, 9.1, 0.63, 0.16]

Tobin_Luminosities_multiples = np.array([0.9, 1.3, 4.2, 0.87, 1.5, 1.3, 3.6, 3.2, 9.1, 0.04, 9.08, 1.5, 4.4, 1.1, 1.19, 18.9, 10.8, 0.79, 11.1, 24.3, 0.9, 2.44, 35.54, 2.1])

#2018 Class 0+I
Tobin_objects = {
"Per-emb-1 (0)": [],
"Per-emb-3 (0)": [],
"Per-emb-9 (0)": [],
"Per-emb-14 (0)": [],
"Per-emb-15 (0)": [],
"Per-emb-19 (0/I)": [],
"Per-emb-20 (0/I)": [],
"Per-emb-23 (0)": [],
"Per-emb-24 (0/I)": [],
"Per-emb-25 (0/I)": [],
"Per-emb-29 (0)": [],
"Per-emb-30 (0/I)": [],
"Per-emb-31 (0/I)": [],
"L1451-MMS (0)": [],
"Per-emb-34 (I)": [],
"Per-emb-38 (I)": [],
"Per-emb-46 (I)": [],
"Per-emb-47 (I)": [],
"Per-emb-50 (I)": [],
"Per-emb-52 (I)": [],
"Per-emb-53 (I)": [],
"Per-emb-54 (I)": [],
"Per-emb-56 (I)": [],
"Per-emb-57 (I)": [],
"Per-emb-61 (I)": [],
"Per-emb-62 (I)": [],
"Per-emb-63 (I)": [],
"Per-emb-64 (I)": [],
"Per-emb-66 (I)": [],
"IRAS 03363+3207 (I?)": [],
"SVS13C (0)": [],
"Per-emb-2": [24.0],
"Per-emb-5": [29.1],
"Per-emb-17": [83.3],
"Per-emb-22": [225.4],
"Per-emb-26+Per-emb-42": [2431.3],
"Per-emb-8+Per-emb-55": [185.3, 2867.2],
"Per-emb-16+Per-emb-28": [4818.9],
"Per-emb-6+Per-emb-10": [9584.2],
"Per-emb-27+Per-emb-36": [93.4, 186.0, 9426.0],
"Per-emb-11": [885.4, 2840.6],
"Per-emb-32": [1820.0],
"Per-emb-37+EDJ2009+235": [],#"Per-emb-37+EDJ2009+235": [3166.7],
"B1-bS+Per-emb-41+B1-bN": [4187.1, 5218.4],
"Per-emb-18+Per-emb-21+Per-emb-49": [25.6, 93.8, 3975.7, 8242.3],
"Per-emb-13+IRAS4B’+Per-emb-12": [548.9, 3196.2, 8921.7],
"Per-emb-44+SVS13A2+SVS13B": [90.0, 1594.2, 4479.7],
"Per-emb-33+L1448IRS3A+L1448NW": [75.3, 79.2, 238.4, 2195.2, 6450.8],
"Per-emb-48": [103.7],
"Per-emb-40": [117.4],
"EDJ2009-183": [307.6],
"L1448IRS1": [427.0],
"Per-emb-35": [572.3],
"Per-emb-58+Per-emb-65": [8663.3]
}

CF_per_bin_Tobin = []
CF_errs = []
#S_true = 17#14

for bin_it in range(1, len(S_bins)):
    N_comps_in_sys = []
    for key in Tobin_objects.keys():
        N_comps = len(Tobin_objects[key]) + 1
        smaller_seps = len(np.argwhere(np.array(Tobin_objects[key]) < S_bins[bin_it-1]))
        larger_seps = len(np.argwhere(np.array(Tobin_objects[key]) > S_bins[bin_it]))
        binaries = len(np.argwhere((np.array(Tobin_objects[key]) < S_bins[bin_it])&(np.array(Tobin_objects[key]) > S_bins[bin_it-1])))
        if N_comps == 1:
            N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 2:
            if larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 3:
            if larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1]
            elif binaries == 2:
                N_comps_in_sys = N_comps_in_sys + [3]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 4:
            if larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 2 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 5:
            if larger_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1, 1]
            elif smaller_seps == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 3 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 3 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif binaries == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1]
            elif smaller_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 6:
            if larger_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1, 1]
            elif binaries == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 3 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 3 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 4 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 4 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
    N_t = len(np.argwhere(np.array(N_comps_in_sys) == 3))
    N_b = len(np.argwhere(np.array(N_comps_in_sys) == 2))
    N_s = len(np.argwhere(np.array(N_comps_in_sys) == 1))
    N_sys = N_s + N_b + N_t
    cf = (N_b+2*N_t)/N_sys
    N_comp = np.sum(np.array(N_comps_in_sys) - 1)
    CF_err = ((N_comp*(1-(N_comp/N_sys)))**0.5)*(1/N_sys)
    #CF_err = (N_comp*(1-(N_comp/N_sys))**0.5)*(1/N_sys)
    
    CF_per_bin_Tobin.append(cf)
    CF_errs.append(CF_err)
    

CF_per_bin_Tobin = np.array(CF_per_bin_Tobin)
CF_errs = np.array(CF_errs)

#=====================================================
#Create plots below

#pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G125/Replace_with_most_massive_star/CF_plots/CF_movies/Unbound_visible_tobin_stars/cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G125/Replace_with_most_massive_star/CF_plots/CF_movies/Bound_all_stars/cf_hist.pkl']


pickle_files = ['/groups/astro/rlk/Analysis_plots/Ramses/Global_early_2021/G100/512/SFE/Unbound/L_Limits/cf_hist.pkl', '/groups/astro/rlk/Analysis_plots/Ramses/Global_early_2021/G100/512/SFE/Bound/No_Limits/cf_hist.pkl']

xlabels = ['Unbound with L limits', 'Bound without L limits']

plt.clf()
fig, axs = plt.subplots(2, 1,sharex=True, sharey=True, figsize=(4, 6.5))
plt.subplots_adjust(wspace=0.0, hspace=0.0)

ax_it = -1
for pickle_file in pickle_files:
    ax_it = ax_it + 1
    file = open(pickle_file, 'rb')
    Times, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
    file.close()

    Summed_systems = np.sum(np.array(N_sys_total), axis=0)
    CF_top = Summed_systems[:,1] + Summed_systems[:,2]*2 + Summed_systems[:,3]*3 + Summed_systems[:,4]*4 + Summed_systems[:,5]*5 + Summed_systems[:,6]*6
    CF_bot = np.sum(Summed_systems, axis=1)
    CF_Total = CF_top/CF_bot

    #Create mean CF
    CF_median = []
    CF_err = []
    N = np.shape(N_sys_total)[0]

    for bit in range(len(S_bins)-1):
        median = np.median(CF_Array_Full[:,bit])
        mean = np.mean(CF_Array_Full[:,bit])
        std = np.std(CF_Array_Full[:,bit], ddof=1)
        standard_deviation = [median-(mean-std), (mean+std)-median]
        CF_median.append(median)
        CF_err.append(standard_deviation)

    CF_err = np.array(CF_err)

    try:
        #axs[ax_it].bar(bin_centers, CF_median, edgecolor='k', label="Simulations", width=0.25, alpha=0.5)
        axs[ax_it].bar(bin_centers, CF_median, yerr=CF_err.T, edgecolor='k', label="Simulations", width=0.25, alpha=0.5)
    except:
        axs[ax_it].bar(bin_centers[1:], CF_median, yerr=CF_err.T, edgecolor='k', label="Simulations", width=0.25, alpha=0.5)
        #axs[ax_it].bar(bin_centers[1:], CF_median, edgecolor='k', label="Simulations", width=0.25, alpha=0.5)
    axs[ax_it].bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al. 2018 (Class 0+I)")
    if ax_it == 0:
        axs[ax_it].legend(loc='upper left')
    if ax_it == 1:
        axs[ax_it].set_xlabel('Log Separation (AU)')
    axs[ax_it].set_ylabel('CF '+xlabels[ax_it])
    axs[ax_it].set_xlim([1,4])
    axs[ax_it].set_ylim(bottom=0.0)
    plt.savefig('Multi_Obs_Median_companion_frequency.pdf', format='pdf', bbox_inches='tight')
    print('created Multi_Obs_Median_companion_frequency.pdf')
