import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import pickle
import os
from mpi4py.MPI import COMM_WORLD as CW
import matplotlib

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

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data_all.pkl", "rb")
CF_per_bin_all, CF_errs_all, Perseus_sep = pickle.load(file_open)
file_open.close()
CF_errs_all[0][np.array([1, 3])] = 0
CF_errs_all[1][np.array([1, 3])] = 0

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data_66.pkl", "rb")
CF_per_bin_66, CF_errs_66, Perseus_sep = pickle.load(file_open)
file_open.close()

CF_errs_66[0][np.array([1, 3])] = 0
CF_errs_66[1][np.array([1, 3])] = 0

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
'''
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
'''
datadir = sys.argv[1]
savedir = sys.argv[2]

args = parse_inputs()
pickle_file = datadir + args.pickled_file
if pickle_file[-4:] != '.pkl':
    pickle_file = pickle_file + '.pkl'

files = glob.glob(pickle_file)
if len(files) == 0:
    file = open(pickle_file, 'rb')
    try:
        Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
        file.close
    except:
        file = open(pickle_file, 'rb')
        Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion, All_separations = pickle.load(file)
        file.close()
else:
    for pickle_file in files:
        file = open(pickle_file, 'rb')
        suffix = pickle_file.split('.pkl')[0].split('_')[-1]
        try:
            Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
            file.close
        except:
            file = open(pickle_file, 'rb')
            Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion, All_separations = pickle.load(file)
            file.close()
        exec("CF_Array_Full_"+suffix+" = CF_Array_Full")
    

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
                
            if len(files)==1:
                CF_median = np.median(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
                CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
                CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
                CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]
                
                plt.clf()
                plt.figure(figsize=(1.5*single_col_width, 1.5*0.4*two_col_width))
                try:
                    #plt.bar(bin_centers, CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
                    plt.bar(bin_centers, CF_median, edgecolor='k', label="Simulation", width=0.25, alpha=0.5)
                except:
                    #plt.bar(bin_centers[1:], CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
                    plt.bar(bin_centers[1:], CF_median, edgecolor='k', label="Simulation", width=0.25, alpha=0.5)
                plt.text((1.03), (0.187), "3D-Full", zorder=11, fontsize=font_size)
            else:
                CF_median_120 = np.median(CF_Array_Full_120[start_integration_it:end_integration_it], axis=0)
                CF_median_55 = np.median(CF_Array_Full_55[start_integration_it:end_integration_it], axis=0)
                #CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
                #CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
                #CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]
                
                plt.clf()
                plt.figure(figsize=(1.5*single_col_width, 1.5*0.4*two_col_width))
                try:
                    #plt.bar(bin_centers, CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
                    plt.bar(bin_centers, CF_median_120, edgecolor='tab:blue', label='$L_{max}=120$L$_\odot$', width=0.25, alpha=0.5)
                    plt.bar(bin_centers, CF_median_55, edgecolor='tab:orange', label='$L_{max}=55$L$_\odot$', width=0.25, alpha=0.5)
                except:
                    #plt.bar(bin_centers[1:], CF_median, yerr=CF_err, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
                    plt.bar(bin_centers[1:], CF_median_120, edgecolor='tab:blue', label='$L_{max}=120$L$_\odot$', width=0.25, alpha=0.5)
                    plt.bar(bin_centers[1:], CF_median_55, edgecolor='tab:orange', label='$L_{max}=55$L$_\odot$', width=0.25, alpha=0.5)
                plt.text((1.1), (0.187), "2D-Bound", zorder=11, fontsize=font_size)
            plt.bar(bin_centers, CF_per_bin_66, yerr=CF_errs_66, width=0.25, edgecolor='black', label="Perseus", fill=None, ls='-')
            plt.bar(bin_centers, CF_per_bin_all, width=0.25, edgecolor='black', fill=None, ls='--')
            #plt.bar(bin_centers, CF_per_bin_all, width=0.25, edgecolor='black', alpha=0.5, fill=None, ls='--')
            #plt.bar(bin_centers, CF_per_bin_Tobin_Ori, width=0.25, edgecolor='black', alpha=0.5, label="Orion", fill=None, ls='-.')
            #plt.plot(x,gauss_total)
            plt.legend(loc='upper right', fontsize=font_size)
            plt.xlabel('Separation (Log$_{10}$(AU))', fontsize=font_size)
            plt.ylabel('Companion Frequency', fontsize=font_size)
            plt.tick_params(axis='both', which='major', labelsize=font_size)
            plt.tick_params(axis='both', which='minor', labelsize=font_size)
            plt.tick_params(axis='x', direction='in')
            plt.tick_params(axis='y', direction='in', right=True)
            plt.xlim([1,4])
            plt.ylim([0, args.y_limit])
            plt.ylim(bottom=0.0)
            
            if args.time_spread != None:
                plt.title('SFE:'+str(np.round(SFE[time_it]*100, decimals=1))+'% ('+str(int((time_val - Times[0])/1000))+'kyr), Integration window:' + str(int(args.time_spread)) + 'yr')
            else:
                plt.title('SFE:'+str(np.round(SFE[time_it]*100, decimals=1))+'% ('+str(int((time_val - Times[0])/1000))+'kyr), Integration window:' + str(args.SFE_spread_val*100) + '% SFE')
            
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
                #/groups/astro/rlk/.cache/matplotlib/tex.cache/bb5fd47356022d30aa7fc0dd54e28e5c.jpg
                #OSError
                '''
                plt.savefig(file_name+'.jpg', format='jpg', bbox_inches='tight')
                print("created "+file_name)
                '''
                try:
                    plt.savefig(file_name+'.jpg', format='jpg', bbox_inches='tight', pad_inches=0.02)
                    print("created "+file_name)
                except:
                    import pdb
                    pdb.set_trace()
                    error_type = sys.exc_info()[0]
                    error_str = sys.exc_info()[1]
                    if str(error_str)[:20] == 'error with jpgfile: ':
                        jpgfile = str(error_str).split(": ")[1].split('.jpg')[0] + '.jpg'
                        os.remove(jpgfile)
                        print("removed jpg file on rank "+ str(rank)+ ", trying to save again")
                        try:
                            plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', pad_inches=0.02)
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
                                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', pad_inches=0.02)
                                print('Created frame of projection', fs, 'on rank', rank, 'at time of', str(int(args_dict['time_val'])), 'to save_dir:', file_name + '.jpg')

                
                
        else:
            print(file_name + " already exists, so skipping it.")

print("Finished making frames")

