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


pickle_files = eval(args.pickle_files)
fig, axs = plt.subplots(len(pickle_files), 1, figsize=(4, 3*len(pickle_files)), constrained_layout=True, sharex=True, sharey=True)
axs = axs.flatten()
plt.subplots_adjust(wspace=0.0, hspace=0.0)

titles = []#['Bound:X', 'Y', 'Z']
if args.subplot_titles != None:
    titles = eval(args.subplot_titles)

fig_counter = 0
for pickle_f in pickle_files:
    file = open(pickle_f, 'rb')
    Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities = pickle.load(file)
    file.close()
    
    sink_inds = []
    system_number = ['[28, 35, 42]', '[84, 109]', '[147, 159]', '[129, 150]', '[119, 120]', '[67, 69, 72]', '[8, 89]', '[12, 97]', '[104, 142, 155]', '[25, 66, 77, 78]', '[124, 156]', '[127, 140, 141, 152]', '[118, 121]', '[146, 151]', '[153, 154]', '[110, 114, 118, 121, 123]', '[64, 116]', '[164, 166]', '[144, 160, 161, 168]', '[139, 143, 157]', '[160, 161, 168]', '[117, 173]', '[173, 174, 178]', '[175, 176, 181]', '[167, 177]', '[48, 92, 111]', '[143, 157]', '[182, 183]', ]
    for key in All_unique_systems_L.keys():
        key_list = eval(key)
        for kit in key_list:
            if kit not in sink_inds:
                sink_inds.append(kit)
    
    #import pdb
    #pdb.set_trace()
    
    CF_median = []
    CF_err = []
    
    for bit in range(11):
        non_zero_inds = np.where(np.array(CF_Array_Full)[:,bit]>0)[0]
        median = np.median(np.array(CF_Array_Full)[:,bit][non_zero_inds])
        mean = np.mean(np.array(CF_Array_Full)[:,bit][non_zero_inds])
        std = np.std(np.array(CF_Array_Full)[:,bit][non_zero_inds])
        standard_deviation = [median-(mean-std), (mean+std)-median]
        CF_median.append(median)
        CF_err.append(standard_deviation)

    CF_err = np.array(CF_err)

    axs[fig_counter].bar(bin_centers, CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
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
