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
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=None)
    parser.add_argument("-SFE_spread", "--SFE_spread_val", help="If you don't want to use the t_spread, you can use SFE_spread", type=float, default=0.001)
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
file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Per, CF_errs_Per = pickle.load(file_open)
file_open.close()

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Orion_data.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Ori, CF_errs_Ori = pickle.load(file_open)
file_open.close()

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/CF_corr.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Per, CF_corr, CF_errs_Per = pickle.load(file_open)
file_open.close()

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data_all.pkl", "rb")
CF_per_bin_all, CF_errs, All_separations = pickle.load(file_open)
file_open.close()

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data_66.pkl", "rb")
CF_per_bin_66, CF_errs, All_separations = pickle.load(file_open)
file_open.close()


units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v

import scipy.stats as stats

#first peak
lower_amp = 0.11
lower_mean = np.log10(120)
lower_std = 0.4
upper_amp = 0.045
upper_mean = np.log10(3500)
upper_std = 0.3
x = np.linspace(1,4,100)
lower_gauss = lower_amp*stats.norm.pdf(x, lower_mean, lower_std)
upper_gauss = upper_amp*stats.norm.pdf(x, upper_mean, upper_std)
gauss_total = lower_gauss + upper_gauss

datadir = sys.argv[1]
savedir = sys.argv[2]

args = parse_inputs()
pickle_file = datadir + args.pickled_file
if pickle_file[-4:] != '.pkl':
    pickle_file = pickle_file + '.pkl'

file = open(pickle_file, 'rb')
try:
    Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
    file.close
except:
    file = open(pickle_file, 'rb')
    Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion, All_separations = pickle.load(file)
    file.close()

if args.starting_ind == 0:
    if args.time_spread != None:
        start_time_it = np.argmin(abs(Times - args.time_spread/2.))
    else:
        start_time_it = np.argmin(abs(SFE - args.SFE_spread_val))
else:
    start_time_it = args.starting_ind
if args.time_spread != None:
    end_time_it = np.argmin(abs(Times - (Times[-1] - args.time_spread/2.)))
else:
    end_time_it = np.argmin(abs(SFE - (SFE[-1] - args.SFE_spread_val)))
rit = -1

for time_it in range(start_time_it, end_time_it):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        file_name = savedir + "movie_frame_" + ("%06d" % time_it)
        if os.path.isfile(file_name+'.jpg') == False:
            time_val = Times[time_it]
            if args.time_spread != None:
                start_time = Times[time_it] - args.time_spread/2.
                end_time = Times[time_it] + args.time_spread/2.
            
                start_integration_it = np.argmin(abs(Times - start_time))
                end_integration_it = np.argmin(abs(Times - end_time))
            else:
                SFE_val = SFE[time_it]
                start_SFE = SFE[time_it] - args.SFE_spread_val
                end_SFE = SFE[time_it] + args.SFE_spread_val
                
                start_integration_it = np.argmin(abs(SFE - start_SFE))
                end_integration_it = np.argmin(abs(SFE - end_SFE))
                
            N_stars = np.sum(N_sys_total[time_it]*np.array([1, 2, 3, 4, 5, 6, 7]))
            
            CF_median = np.median(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
            CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
            CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
            CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]
            
            plt.clf()
            try:
                #plt.bar(bin_centers, CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
                plt.bar(bin_centers, CF_median, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
            except:
                #plt.bar(bin_centers[1:], CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
                plt.bar(bin_centers[1:], CF_median, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
            plt.bar(bin_centers, CF_per_bin_Tobin_Per, yerr=CF_errs, width=0.25, edgecolor='black', alpha=0.5, label="Perseus", fill=None, ls='-')
            #plt.bar(bin_centers, CF_per_bin_all, width=0.25, edgecolor='black', alpha=0.5, fill=None, ls='--')
            #plt.bar(bin_centers, CF_per_bin_Tobin_Ori, width=0.25, edgecolor='black', alpha=0.5, label="Orion", fill=None, ls='-.')
            #plt.plot(x,gauss_total)
            plt.legend(loc='upper left')
            plt.xlabel('Log Separation (AU)')
            plt.ylabel('Companion Frequency')
            plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
            plt.ylim([0, args.y_limit])
            plt.ylim(bottom=0.0)
            if args.time_spread != None:
                plt.title("SFE:"+str(np.round(SFE[time_it]*100, decimals=1))+"\% ("+str(int((time_val - Times[0])/1000))+"kyr), Integration window:" + str(args.time_spread) + "yr, N_stars:" + str(N_stars))
            else:
                plt.title("SFE:"+str(np.round(SFE[time_it]*100, decimals=1))+"\% ("+str(int((time_val - Times[0])/1000))+"kyr), Integration window:" + str(args.SFE_spread_val*100) + "\% SFE, Nstars:" + str(N_stars))
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

