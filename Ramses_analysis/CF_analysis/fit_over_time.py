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

if rank == 0:
    plt.clf()
    plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
    plt.ylabel("Companion Frequency")
    plt.xlabel("Log (AU)")
    plt.xlim([1,4])
    plt.ylim([0, 0.2])
    
    plt.savefig("Tobin_2018_class_0_I.png")

sys.stdout.flush()
CW.Barrier()
        
file = open("Tobin_CF.pkl", 'wb')
pickle.dump((S_bins, CF_per_bin_Tobin, CF_errs),file)
print('updated pickle', "Tobin_CF.pkl")
file.close()

#=====================================================================================================

pickle_file = args.pickled_file
try:
    file = open(pickle_file+'.pkl', 'rb')
    Times, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
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
usable_bin_inds = np.array([ 2,  4,  5,  6,  7,  8,  9, 10, 11, 12])
Times = (Times-Times[0])/1e6
Non_bimodal_dist = np.array([0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 ,
       0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 , 0.0379096 ,
       0.05457945, 0.09212946, 0.12967947])
Bimodal_dist = np.array([0.00780592, 0.01872544, 0.03498372, 0.05090101, 0.05767853,
       0.05090462, 0.03504812, 0.01944561, 0.01284564, 0.02460516,
       0.06113053, 0.10387658, 0.11138321])
CF_it = -1
for CF_hist in CF_Array_Full:
    CF_it = CF_it + 1
    chi_red_tobin = (np.sum(((CF_hist[usable_bin_inds]-CF_per_bin_Tobin[usable_bin_inds])**2)/(CF_errs[usable_bin_inds]**2)))/len(CF_hist[usable_bin_inds])
    reduced_chi_square_tobin.append(chi_red_tobin)
    if chi_red_tobin < 0.5:
        plt.clf()
        plt.bar(bin_centers, CF_hist, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
        plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
        plt.legend(loc='best')
        plt.xlabel('Log Separation (AU)')
        plt.ylabel('Companion Frequency')
        plt.xlim([1, 4])
        plt.ylim(bottom=0.0)
        plt.title("reduced chi$^2$ = " + str(chi_red_tobin))
        plt.savefig(args.save_directory+'Time_'+str(Times[CF_it])+'Myr_chi_tobin.png', format='png', bbox_inches='tight')
        print('saved', args.save_directory+'Time_'+str(Times[CF_it])+'Myr_chi_tobin.png')
    
plt.clf()
plt.semilogy(Times, reduced_chi_square_tobin)
plt.xlabel("Time (Myr)")
plt.ylabel("mean Chi squared (<$\chi^2$>)")
plt.xlim([Times[0], Times[-1]])
plt.ylim([np.min(reduced_chi_square_tobin), np.max(reduced_chi_square_tobin)])
plt.savefig(args.save_directory + "reduced_chi_squared_tobin.png")
'''
plt.clf()
plt.plot(Times, mean_bi_likelihoods)
plt.xlabel("Time (Myr)")
plt.ylabel("mean likelihood")
plt.axhline(y=0.5, c='k')
plt.xlim([Times[0], Times[-1]])
plt.ylim([0, 1])
plt.savefig(args.save_directory + "mean_likelihood.png")

plt.clf()
plt.plot(Times, likelihood_ratio)
plt.xlabel("Time (Myr)")
plt.ylabel("mean likelihood ratio")
plt.axhline(y=1, c='k')
plt.ylim(bottom=0)
plt.xlim([Times[0], Times[-1]])
plt.savefig(args.save_directory + "mean_likelihood_ratio.png")

#Integrated CF
ts = args.time_spread/1e6
start_time_it = np.argmin(abs(Times - ts))
end_time_it = np.argmin(abs(Times - (Times[-1] - ts)))
chi_squared_int = []
mean_b_likelihoods = []
mean_nb_likelihoods = []
likelihood_ratio = []
for time_it in range(start_time_it, end_time_it):
    time_val = Times[time_it]
    start_time = Times[time_it] - ts
    end_time = Times[time_it] + ts
            
    start_integration_it = np.argmin(abs(Times - start_time))
    end_integration_it = np.argmin(abs(Times - end_time))
    
    CF_median = np.median(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_mean = np.mean(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_std = np.std(CF_Array_Full[start_integration_it:end_integration_it], axis=0)
    CF_err = [CF_median-(CF_mean-CF_std), (CF_mean+CF_std)-CF_median]
    
    chi_per_bin = ((CF_median - Bimodal_dist)**2)/Bimodal_dist
    chi_squared_int.append(np.mean(chi_per_bin))
    b_likelihoods = np.exp(((-1*(CF_mean[usable_bin_inds] - Bimodal_dist[usable_bin_inds])**2)/(2*CF_std[usable_bin_inds]**2)))
    if np.median(b_likelihoods) != 0:
        mean_bl = np.mean(b_likelihoods)
    else:
        mean_bl = np.nan
    mean_b_likelihoods.append(mean_bl)
    nb_likelihoods = np.exp(((-1*(CF_mean[usable_bin_inds] - Non_bimodal_dist[usable_bin_inds])**2)/(2*CF_std[usable_bin_inds]**2)))
    if np.median(nb_likelihoods) != 0:
        mean_nbl = np.mean(nb_likelihoods)
    else:
        mean_nbl = np.nan
    mean_nb_likelihoods.append(mean_nbl)
    ratio = mean_bl/mean_nbl
    likelihood_ratio.append(ratio)

    """
    chi_per_bin = ((CF_median[usable_bin_inds] - CF_per_bin_Tobin[usable_bin_inds])**2)/CF_median[usable_bin_inds]
    chi_squared_int.append(np.mean(chi_per_bin))
    b_likelihoods = np.exp(((-1*(CF_mean[usable_bin_inds] - CF_per_bin_Tobin[usable_bin_inds])**2)/(2*CF_std[usable_bin_inds]**2)))
    mean_b_likelihoods.append(np.mean(b_likelihoods))
    nb_likelihoods = np.exp(((-1*(CF_mean[usable_bin_inds] - Non_bimodal_dist[usable_bin_inds])**2)/(2*CF_std[usable_bin_inds]**2)))
    mean_nb_likelihoods.append(np.mean(nb_likelihoods))
    ratio = np.mean(b_likelihoods)/np.mean(nb_likelihoods)
    likelihood_ratio.append(ratio)
    """
    if ratio > 5.0:
        plt.clf()
        plt.bar(bin_centers, CF_mean, yerr=CF_std, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
        plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
        plt.legend(loc='best')
        plt.xlabel('Log Separation (AU)')
        plt.ylabel('Companion Frequency')
        plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
        plt.ylim(bottom=0.0)
        plt.title("Likelihood ratio = " + str(ratio))
        plt.savefig(args.save_directory+'Time_'+str(time_val)+'Myr_int.png', format='png', bbox_inches='tight')
        print('saved', args.save_directory+'Time_'+str(time_val)+'Myr_int.png')
    #print("done time", time_it, "in", end_time_it)

plt.clf()
plt.plot(Times[start_time_it:end_time_it], chi_squared_int)
plt.xlabel("Time (Myr)")
plt.ylabel("mean Chi squared (<$\chi^2$>)")
plt.xlim([Times[0], Times[-1]])
plt.ylim([np.min(chi_squared_int), np.max(chi_squared_int)])
plt.savefig(args.save_directory + "chi_squared_int_"+str(args.time_spread)+".png")

plt.clf()
plt.plot(Times[start_time_it:end_time_it], mean_b_likelihoods)
plt.xlabel("Time (Myr)")
plt.ylabel("mean likelihood")
plt.axhline(y=0.5, c='k')
plt.xlim([Times[0], Times[-1]])
plt.ylim([0, 1])
plt.savefig(args.save_directory + "mean_likelihood_int_"+str(args.time_spread)+".png")

plt.clf()
plt.plot(Times[start_time_it:end_time_it], likelihood_ratio)
plt.xlabel("Time (Myr)")
plt.ylabel("mean likelihood ratio")
plt.axhline(y=1, c='k')
plt.ylim(bottom=0)
plt.xlim([Times[0], Times[-1]])
plt.savefig(args.save_directory + "mean_likelihood_ratio_int_"+str(args.time_spread)+".png")
'''


[(2M05414580-0154297-B-2M05414580-0154297-A)
(2M05414550-0154286-H384-region-D))
OMC1N-1-C+OMC1N-1-B
H384-region-G+[2M05414325-0154343+(NGC2024-FIR3-A+NGC2024-FIR3-B)]
H76+[[H78-D+(H78-C+H78-B)]+(H78-E+H78-A)]
(H361-G-B+H361-G-A)+H361-H
(H369+[(OMC2-FIR4-ALMA1+VLA16)+VLA15])+(H108+H64)
(H181-A+H181-B)+[H182-C+(H182-B+H182-A)]
OMC1N-6-7-8-H+OMC1N-6-7-8-G
H117+H118
H322+(H389-A+H389-B)
(H86-B+H86-A)+H87
H358-B+H358-A
H384-region-E+2M05414512-0154470
(H361-D+[(H361-F+[(H361-E-A+H361-E-B)+H361-A])+
[H361-B+(H361-C-A+H361-C-B)]])+[(H361-G-B+H361-G-A)+H361-H]
OMC1N-6-7-8-D+OMC1N-6-7-8-E
[H92-B+(H92-A-B+H92-A-A)]+J05351805-050017.98
(H323-B+H323-A)+[H322+(H389-A+H389-B)]
(OMC1N-6-7-8-H+OMC1N-6-7-8-G)+OMC1N-6-7-8-F
H240+H241
H93+H94
H409+(H59-A+H59-B)
(H377+H144)+H143
(2MJ05352746-0509441+H370)+(H66-A+H66-B)
[(2M05414483-0154357+2M05414482-0154251)+
[(2M05414580-0154297-B+2M05414580-0154297-A)+-
(2M05414550-0154286+H384-region-D)]]+-
(2M05414611-0154147+H384-region-H)
H89+H91
H374+H254
(OMC1N-6-7-8-C+OMC1N-6-7-8-B)+(OMC1N-6-7-8-A-B+OMC1N-6-7-8-A-A)
([(H71-A+H71-B)+H394-B]+H72)+H394-A
(H384-region-E+2M05414512-0154470)+(H384-C+[(H384-A+H384-A-B)+H384-B])
H219+H220
H407+H331
(H183-A+H183-B)+[(H181-A+H181-B)+[H182-C+(H182-B+H182-A)]]
H210+H211
J05352074-0515492+[V2358Ori+(H56-B+[(H56-A-A+H56-A-B)+H56-A-C])]
HH270mms2+(HH270mms1-A+HH270mms1-B)
H252+(H255-A+H255-B)
H257+(H261-B+H261-A)

C0+C0
C0+[(C0+C0)+C0]

7615.4±0.7
7629.0±2.7
7766.9±8.6
7946.2±1.9
8144.2±2.7
8199.3±15.8
8623.4±1.1
8968.7±1.0
9130.3±1.8
9138.2±3.3
9212.0±0.9
9415.0±1.3
9459.9±3.9
9999.5±6.0
