import numpy as np
import matplotlib.pyplot as plt
import yt
import glob
import sys
import matplotlib as mpl
import pickle
import os

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-TB", "--Tobin_pickle", type=str)
    parser.add_argument("-pickles", "--pickle_files", type=str)
    parser.add_argument("-save_name", "--save_name_of_image", type=str)
    parser.add_argument("-titles", "--subplot_titles", type=str, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

file = open(args.Tobin_pickle, 'rb')
S_bins, CF_per_bin_Tobin = pickle.load(file)
file.close()
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2

S_bins_lower_sep = np.logspace(0.75,4,14)
bin_centers_lower_sep = (np.log10(S_bins_lower_sep[:-1])+np.log10(S_bins_lower_sep[1:]))/2

pickle_files = eval(args.pickle_files)
fig, axs = plt.subplots(len(pickle_files), 1, figsize=(4, 3*len(pickle_files)), constrained_layout=True, sharex=True, sharey=True)
axs = axs.flatten()
plt.subplots_adjust(wspace=0.0, hspace=0.0)

titles = []#['Bound:X', 'Y', 'Z']
if args.subplot_titles != None:
    titles = eval(args.subplot_titles)

fig_counter = 0
for pickle_f in pickle_files:
    print("opening file", pickle_f)
    file = open(pickle_f, 'rb')
    Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities = pickle.load(file)
    file.close()
    
    #CF calculated from all systems
    Summed_systems = np.sum(np.array(N_sys_total), axis=0)
    CF_top = Summed_systems[:,1] + Summed_systems[:,2]*2 + Summed_systems[:,3]*3 + Summed_systems[:,4]*4 + Summed_systems[:,5]*5 + Summed_systems[:,6]*6
    CF_bot = np.sum(Summed_systems, axis=1)
    CF_Total = CF_top/CF_bot
    
    print("CF_Total =", CF_Total)
    
    #Create mean CF
    CF_median = np.median(np.array(CF_Array_Full), axis=0)
    mean = np.mean(np.array(CF_Array_Full), axis=0)
    std = np.std(np.array(CF_Array_Full), axis=0)
    CF_err = np.array([CF_median-(mean-std), (mean+std)-CF_median])
    
    CF_err = np.array(CF_err)
    print("CF_median =", CF_median)

    axs[fig_counter].bar(bin_centers_lower_sep, CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
    axs[fig_counter].bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
    if fig_counter == 0:
        axs[fig_counter].legend(loc='best')
    if fig_counter == len(pickle_files)-1:
        axs[fig_counter].set_xlabel('Log Separation (AU)')
    if fig_counter != len(pickle_files)-1:
        plt.setp([axs[fig_counter].get_xticklabels()], visible=False)
        plt.setp([axs[fig_counter].get_yticklabels()[0]], visible=False)
    axs[fig_counter].set_xlim([1,4])
    axs[fig_counter].tick_params(axis='y', direction="in")
    #axs[fig_counter].set_aspect('equal')
    if len(titles)!=0:
        axs[fig_counter].set_ylabel(titles[fig_counter] +' Companion Frequency')
    else:
        axs[fig_counter].set_ylabel('Companion Frequency')
    
    fig_counter = fig_counter + 1

#plt.xlabel('Separation (AU)')
plt.savefig(args.save_name_of_image, bbox_inches='tight', pad_inches=0.02)
print('created', args.save_name_of_image)
