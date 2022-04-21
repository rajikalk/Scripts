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
from scipy.stats import norm
from scipy.optimize import curve_fit
import subprocess

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=10000)
    parser.add_argument("-start_ind", "--starting_ind", help="What frame do you want to start at?", type=int, default=0)
    parser.add_argument("-y_lim", "--y_limit", help="do you want to set a ylim for the frames?", type=float, default=0.2)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
rank = CW.Get_rank()
size = CW.Get_size()
#os.system('export MPLCONFIGDIR="./Process'+str(rank)+'/temp/matpllotlib/"')

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

datadir = sys.argv[1]
savedir = sys.argv[2]

args = parse_inputs()
pickle_file = datadir + args.pickled_file
if pickle_file[-4:] != '.pkl':
    pickle_file = pickle_file + '.pkl'

file = open(pickle_file, 'rb')
Times, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
file.close()

if args.starting_ind == 0:
    start_time_it = np.argmin(abs(Times - args.time_spread))
else:
    start_time_it = args.starting_ind
end_time_it = np.argmin(abs(Times - (Times[-1] - args.time_spread)))
rit = -1

for time_it in range(start_time_it, end_time_it):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        file_name = savedir + "movie_frame_" + ("%06d" % time_it)
        if os.path.isfile(file_name+'.jpg') == False:
            time_val = Times[time_it]
            start_time = Times[time_it] - args.time_spread
            end_time = Times[time_it] + args.time_spread
            
            start_integration_it = np.argmin(abs(Times - start_time))
            end_integration_it = np.argmin(abs(Times - end_time))
            
            CF_median = np.median(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
            CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
            CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
            CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]
            
            plt.clf()
            try:
                plt.bar(bin_centers, CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
            except:
                plt.bar(bin_centers[1:], CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
            plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
            plt.legend(loc='best')
            plt.xlabel('Log Separation (AU)')
            plt.ylabel('Companion Frequency')
            plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
            plt.ylim([0, args.y_limit])
            plt.ylim(bottom=0.0)
            plt.title("Time:"+str(time_val - Times[0])+"yr, Median over " + str(end_integration_it-start_integration_it) + " histograms")
            if size > 1:
                try:
                    plt.savefig(file_name+'.jpg', format='jpg', bbox_inches='tight')
                    print("created "+file_name)
                except:
                    error_type = sys.exc_info()[0]
                    error_str = sys.exc_info()[1]
                    if str(error_type).split("'")[1] == 'TimeoutError':
                        try:
                            lock_file = str(error_str).split('\n')[1][4:]
                            os.remove(lock_file)
                            print("removed lock-file and trying to save again")
                            plt.savefig(file_name+'.jpg', format='jpg', bbox_inches='tight')
                            print("created "+file_name)
                        except:
                            print("not a lock file problem")
                    else:
                        print("FAILED ON TIME_IT", time_it)
            else:
                #/groups/astro/rlk/.cache/matplotlib/tex.cache/bb5fd47356022d30aa7fc0dd54e28e5c.png
                #OSError
                '''
                plt.savefig(file_name+'.jpg', format='jpg', bbox_inches='tight')
                print("created "+file_name)
                '''
                try:
                    plt.savefig(file_name+'.jpg', format='jpg', bbox_inches='tight')
                    print("created "+file_name)
                except:
                    error_type = sys.exc_info()[0]
                    error_str = sys.exc_info()[1]
                    if str(error_str)[:20] == 'error with pngfile: ':
                        pngfile = str(error_str).split(": ")[1].split('.png')[0] + '.png'
                        os.remove(pngfile)
                        print("removed png file on rank "+ str(rank)+ ", trying to save again")
                        try:
                            plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                            print('Created frame of projection', fs, 'on rank', rank, 'at time of', str(int(args_dict['time_val'])), 'to save_dir:', file_name + '.jpg')
                        except:
                            error_type = sys.exc_info()[0]
                            error_str = sys.exc_info()[1]
                            if str(error_type).split("'")[1] == 'TimeoutError':
                                lock_file = str(error_str).split('\n')[1][4:]
                                try:
                                    os.remove(lock_file)
                                    print("removed lock-file on rank "+ str(rank)+ ",and trying save again")
                                except:
                                    print("lock_file already removed. Trying to save again on rank "+ str(rank))
                                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                                print('Created frame of projection', fs, 'on rank', rank, 'at time of', str(int(args_dict['time_val'])), 'to save_dir:', file_name + '.jpg')

                
                
        else:
            print(file_name + " already exists, so skipping it.")

print("Finished making frames")

