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
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
rank = CW.Get_rank()
size = CW.Get_size()
#os.system('export MPLCONFIGDIR="./Process'+str(rank)+'/temp/matpllotlib/"')

S_bins = np.logspace(0.75,4,14)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_CDF.pkl", "rb")
Perseus_sep = pickle.load(file_open)
file_open.close()

Perseus_log_sep = np.log10(np.sort(Perseus_sep))
Perseus_frequency = np.arange(len(Perseus_sep))/np.arange(len(Perseus_sep))[-1]

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v

import scipy.stats as stats
datadir = sys.argv[1]
savedir = sys.argv[2]

args = parse_inputs()
pickle_file = datadir + args.pickled_file
if pickle_file[-4:] != '.pkl':
    pickle_file = pickle_file + '.pkl'

file = open(pickle_file, 'rb')
Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion, All_separations = pickle.load(file)
file.close()

for time_it in range(len(Times)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        file_name = savedir + "movie_frame_" + ("%06d" % time_it)
        if os.path.isfile(file_name+'.jpg') == False:
            log_sep = np.log10(np.sort(All_separations[time_it]))
            frequency = np.arange(len(All_separations[time_it]))/np.arange(len(All_separations[time_it]))[-1]
            
            plt.clf()
            plt.step(Perseus_log_sep, Perseus_frequency, label='Perseus')
            plt.step(log_sep, frequency, label='Simulation')
            plt.legend(loc='upper left')
            plt.xlabel('Separation (Log$_{10}$(AU))')
            plt.ylabel('Frequency')
            plt.xlim([0, 4])
            plt.ylim([0, 1])
            plt.title("SFE:"+str(np.round(SFE[time_it]*100, decimals=1))+"\% ("+str(int((time_val - Times[0])/1000))+"kyr)")
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

