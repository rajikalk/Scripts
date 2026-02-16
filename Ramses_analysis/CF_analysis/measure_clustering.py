import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
from mpi4py.MPI import COMM_WORLD as CW
import sys
import os
from scipy.optimize import curve_fit
import pickle

def power_law(x, a, k):
    return a * np.power(x, k)
    
def line(x, m, b):
    return m*x + b

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-uniform_dir", "--uniform_dist_dir", help="what directory do you want to use?", default='/groups/astro/rlk/Analysis_plots/Clustering_measure/Underlying_dist/', type=str)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

#Set units
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

sys.stdout.flush()
CW.Barrier()

#Generate underlying distribution:
n_bin_bound = 11
exp_bins = np.linspace(1, 6, n_bin_bound)
sep_bins = 10**exp_bins
sep_centers = np.linspace(np.mean(exp_bins[:2]), np.mean(exp_bins[-2:]), len(exp_bins)-1)
if os.path.exists(args.uniform_dist_dir + 'all_uniform_dist.pkl')==False:
    uni_hists = []
    pickle_files = glob.glob(args.uniform_dist_dir + 'uniform_dist_*.pkl')
    for pickle_file in pickle_files:
        file_open = open(pickle_file, 'rb')
        sep_his_uni = pickle.load(file_open)
        file_open.close()
        uni_hists.append(sep_his_uni)
else:
    file_open = open(args.uniform_dist_dir + 'all_uniform_dist.pkl', 'rb')
    uni_hists = pickle.load(file_open)
    file_open.close()
    
#Calculate uniform distribution
uni_hist_sum = np.sum(uni_hists,axis=0)
if uni_hist_sum[0] == 0:
    uni_hists[0][0] = 1
    max_bin = np.argmax(uni_hists[0])
    uni_hists[0][max_bin] = uni_hists[0][max_bin] - 1
P_RR = np.mean(uni_hists,axis=0)
P_RR_err = np.std(uni_hists,axis=0)
N_R = int(args.uniform_dist_dir.split('/')[-2].split('_')[0])
RR = P_RR/(N_R * (N_R - 1))
RR_err = P_RR_err/(N_R * (N_R - 1))
RR_rel_err = RR_err/RR

#loading global data
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
print('loaded global data')

sys.stdout.flush()
CW.Barrier()

M_stars = np.sum(global_data['m'],axis=1)*units['mass_unit'].in_units('Msun')
SFE = M_stars/units['mass_unit'].in_units('Msun')
SFE_5_ind = np.argmin(abs(SFE-0.05))

time_end_it = SFE_5_ind
rit = -1
time_it_range = range(0, time_end_it+1)
#time_it_range = [time_end_it]
exp_fits = []
exp_err = []#error on the powerlaw index
saved_t_ind = []
for time_it in time_it_range:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        n_stars = np.where(global_data['m'][time_it]>0)[0]
        N_D = len(n_stars)
        #define random stars
        if N_D > 1:
            n_bin_bound = 11
            all_bins_filled = False
            #while all_bins_filled == False:
            exp_bins = np.linspace(1, 6, n_bin_bound)
            sep_bins = 10**exp_bins
            sep_centers = np.linspace(np.mean(exp_bins[:2]), np.mean(exp_bins[-2:]), len(exp_bins)-1)
            width = (sep_centers[1:] - sep_centers[:-1])[0]
            separations = []
            for n_star in n_stars:
                cur_x = global_data['x'][time_it][n_star]
                cur_y = global_data['y'][time_it][n_star]
                cur_z = global_data['z'][time_it][n_star]
                
                other_x = global_data['x'][time_it][n_star+1:n_stars[-1]+1]
                other_y = global_data['y'][time_it][n_star+1:n_stars[-1]+1]
                other_z = global_data['z'][time_it][n_star+1:n_stars[-1]+1]
                
                dx = other_x - cur_x
                dy = other_y - cur_y
                dz = other_z - cur_z
                
                if np.sum(abs(dx) > 0.5)>0:
                    update_x = np.argwhere(abs(dx) > 0.5).T[0]
                    dx[update_x] = abs(dx[update_x]) - 1.0
                    #print('updated dx')
                if np.sum(abs(dy) > 0.5)>0:
                    update_y = np.argwhere(abs(dy) > 0.5).T[0]
                    dy[update_y] = abs(dy[update_y]) - 1.0
                    #print('updated dy')
                if np.sum(abs(dz) > 0.5)>0:
                    update_z = np.argwhere(abs(dz) > 0.5).T[0]
                    dz[update_z] = abs(dz[update_z]) - 1.0
                    #print('updated dz')
                
                sep = np.sqrt(dx**2 + dy**2 + dz**2)
                #if np.sum(sep > 0.8660254037844386) > 0:
                #    import pdb
                #    pdb.set_trace()
                separations = separations + sep.tolist()
            P_DD, bins = np.histogram(np.array(separations)*units['length_unit'].in_units('AU'), bins=sep_bins)
            P_DD_err = np.sqrt(P_DD)
            DD = P_DD/(N_D*(N_D-1))
            DD_err = P_DD_err/(N_D*(N_D-1))
            DD_rel_err = DD_err/DD
            DD_rel_err = np.nan_to_num(DD_rel_err, nan=1e15)
            TPCF = DD/RR
            TPCF_rel_err = DD_rel_err + RR_rel_err
            TPCF_err = TPCF_rel_err * (TPCF+1e-5)
            TPCF_frac = TPCF #+1
            
            plt.clf()
            plt.errorbar(10**sep_centers, TPCF_frac, yerr=TPCF_err, fmt = 'o')
            #try:
                #power_law_break_ind = np.where(TPCF_frac>30)[0][-1] + 2
                
            power_law_break_ind = 6
            if TPCF_frac[power_law_break_ind-1] == 0:
                power_law_break_ind = 5
            dy = np.log10(TPCF_frac[:power_law_break_ind])[-1]-np.log10(TPCF_frac[:power_law_break_ind])[0]
            dx = sep_centers[:power_law_break_ind][-1] - sep_centers[:power_law_break_ind][0]
            grad_guess = dy/dx
            y_intercept_guess = np.log10(TPCF_frac[:power_law_break_ind][0]) - (grad_guess * sep_centers[:power_law_break_ind][0])
            power_law_break_ind = 6
        
            popt1, pcov1 = curve_fit(line, sep_centers[:power_law_break_ind], np.log10(TPCF_frac[:power_law_break_ind]+1e-5), p0=[grad_guess, y_intercept_guess], sigma=np.log10(TPCF_err/(TPCF_frac+1e-5))[:power_law_break_ind], absolute_sigma=True)
            plt.loglog(10**sep_centers[:power_law_break_ind+1], 10**line(sep_centers[:power_law_break_ind+1], popt1[0], popt1[1]), ls='--', color='k')
            #popt2, pcov2 = curve_fit(line, sep_centers[power_law_break_ind:], np.log10(TPCF_frac[power_law_break_ind:]))
            #plt.loglog(10**sep_centers[power_law_break_ind-1:], 10**line(sep_centers[power_law_break_ind-1:], popt2[0], popt2[1]), ls='--', color='k')
            exp_err.append(np.sqrt(np.diag(pcov1))[0])
            
            exp_fits.append(popt1[0])
            saved_t_ind.append(time_it)
            #except:
            #    print("can't make powerlaw fit")
            
            plt.xscale("log")
            plt.yscale("log")
            plt.ylabel("$1+\\omega_(r)$")
            plt.xlabel('Separation (Log$_{10}$(AU))')
            plt.xlim([10**sep_centers[0], 10**sep_centers[-1]])
            plt.savefig("movie_frame_" + ("%06d" % time_it) + ".jpg")
            print('made ' + "movie_frame_" + ("%06d" % time_it) + ".jpg" +" on rank", rank)

            if time_it == time_end_it:
                SFE_5_pickle = 'SFE_5_TPCF.pkl'
                file = open(SFE_5_pickle, 'wb')
                pickle.dump((sep_centers, TPCF_frac, TPCF_err, power_law_break_ind, popt1), file)
                file.close()

#write pickle
pickle_name = 'TPCF_rank_' + ("%06d" % rank) + '.pkl'
file = open(pickle_name, 'wb')
pickle.dump((saved_t_ind, exp_fits, exp_err), file)
file.close()

print('Finished on rank', rank)

sys.stdout.flush()
CW.Barrier()

if rank == 0:
    #compile pickles
    exp_fits_full = np.zeros(np.shape(time_it_range))
    exp_err_full = np.zeros(np.shape(time_it_range))
    pickles = glob.glob('TPCF_rank_*.pkl')
    for pickle_file in pickles:
        file = open(pickle_file, 'rb')
        saved_t_ind, exp_fits, exp_err = pickle.load(file)
        file.close()
        
        exp_fits_full[saved_t_ind] = exp_fits
        exp_err_full[saved_t_ind] = exp_err
        #os.remove(pickle_file)

    pickle_name = 'TPCF.pkl'
    file = open(pickle_name, 'wb')
    pickle.dump((SFE[:time_end_it+1], exp_fits_full, exp_err_full), file)
    file.close()
    
sys.stdout.flush()
CW.Barrier()
print('Compiled all pickles')
