import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import pickle
import matplotlib.patches
import collections
from mpi4py.MPI import COMM_WORLD as CW

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifetime threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("-ts", "--timescale", help="Do you want to plot in terms of the of the actual simulation time, or the the time since the formation of the first sink?", type=str, default="first_sink")
    parser.add_argument("-fig_suffix", "--figure_suffix", help="Do you want add a suffix to the figures?", type=str, default="")
    parser.add_argument("-add_hist", "--add_histograms", help="Do you want to add the histograms at the end of the super plots?", type=str, default="True")
    parser.add_argument("-all_sep_evol", "--plot_all_separation_evolution", help="do you want to plot all separation evolution?", type=str, default="True")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

rank = CW.Get_rank()
size = CW.Get_size()

plot_booleans = [[True, True], [True, False], [False, True], [False, False]]

args = parse_inputs()

pickle_files = ["/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G50/superplot_with_means.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G100/256/superplot_with_means.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G125/superplot_with_means.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G150/superplot_with_means.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G200/superplot_with_means.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/G400/superplot_with_means.pkl"]

plt.clf()
if plot_booleans[rank][0] == True:
    fig, axs = plt.subplots(len(pickle_files), 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, len(pickle_files)*3), sharex='col')
    iter_range = range(0, len(pickle_files)*2, 2)
    args.figure_suffix = args.figure_suffix + '_hist'
else:
    fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(12, len(pickle_files)*3))
    iter_range = range(0, len(pickle_files)*2)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.07)

if plot_booleans[rank][1] == True:
    args.figure_suffix = args.figure_suffix + '_thres'

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
GMC_mass_arr = [1500, 3000, 3750, 4500, 6000, 12000]
Shrinkage_hist = []
Shrinkage_bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
Shrinks = []
Lifetime_bins = [0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
Lifetime_hist = []
S_bins = np.logspace(1,4,13)
CF_hist = np.zeros((len(pickle_files),12)).tolist()
Times_all = []
M_tot_all = []
M_tot_multi_all = []
N_stars_all = []
N_multi_stars_all = []
SFE_all = []
t_max = np.nan
G50_t_max = 0

if args.plot_all_separation_evolution == "True":
    for pick_it in iter_range:
        file = open(pickle_files[int(pick_it/2)], 'rb')
        Times, SFE, SFE_n, M_tot, M_tot_vis, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_ecc, Sink_formation_times, System_mean_times, System_mean_seps, System_mean_ecc, System_lifetimes, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
        file.close()
        
        #Add SFE to the big array to plot later
        Times_all.append(Times)
        SFE_all.append(SFE)
        M_tot_all.append(M_tot)
        M_tot_multi_all.append(M_tot_multi)
        N_stars_all.append(N_stars)
        N_multi_stars_all.append(N_multi_stars)
        
        shrink = 100*((np.array(Initial_Seps) - np.array(Final_seps))/np.array(Initial_Seps))
        Shrinks.append(shrink)
        shrink_hist, bins = np.histogram(shrink, bins=Shrinkage_bins)
        Shrinkage_hist.append(shrink_hist/np.sum(shrink_hist))
        life_hist, bins = np.histogram(list(System_lifetimes.values()), bins=Lifetime_bins)
        Lifetime_hist.append(life_hist/np.sum(life_hist))
    
    
        sim_start_time = np.nan
        for time_key in System_times.keys():
            for sit in range(len(S_bins[1:])):
                bin_inst = np.argwhere((System_seps[time_key]>S_bins[sit])&(System_seps[time_key]<S_bins[sit+1]))
                CF_hist[int(pick_it/2)][sit] = CF_hist[int(pick_it/2)][sit] + np.shape(bin_inst)[0]
            key_inds = flatten(eval(time_key))
            if args.timescale == "first_sink":
                first_sink = np.argmin(Sink_formation_times[key_inds])
                Time_adjusted_formation = System_times[time_key] - Sink_formation_times[key_inds[first_sink]].value
            else:
                if np.isnan(sim_start_time):
                    sim_start_time = System_times[time_key][0]
                Time_adjusted_formation = System_times[time_key] - np.array(sim_start_time)
            if np.isnan(t_max):
                if Time_adjusted_formation[-1] > G50_t_max:
                    G50_t_max = Time_adjusted_formation[-1]
            axs.flatten()[pick_it].semilogy(Time_adjusted_formation, System_seps[time_key], alpha=0.01, color='k')
            if plot_booleans[rank][1] == False:
                axs.flatten()[pick_it].scatter(np.ones(np.shape(System_seps[time_key][0]))*Time_adjusted_formation[0], System_seps[time_key][0], color='m', marker='o')
            else:
                if System_lifetimes[time_key] > args.sys_lifetime_threshold:
                    axs.flatten()[pick_it].scatter(Time_adjusted_formation[0], System_seps[time_key][0][np.argmax(System_seps[time_key][0])], color='m', marker='o', s=0.5)
            '''
            if System_lifetimes[time_key] > args.sys_lifetime_threshold:
                window = 50000
                mean_time = []
                mean_sep = []
                for time in Time_adjusted_formation:
                    start_time = time - window/2.
                    end_time = time + window/2.
                    start_ind = np.argmin(abs(np.array(Time_adjusted_formation) - start_time))
                    end_ind = np.argmin(abs(np.array(Time_adjusted_formation) - end_time))
                    if end_ind != start_ind:
                        non_nan_inds = np.argwhere(np.isnan(np.array(System_seps[time_key][start_ind:end_ind]).T[0])==False).T[0]
                        mean_t = np.mean(np.array(Time_adjusted_formation[start_ind:end_ind])[non_nan_inds])
                        mean_s = np.mean(np.array(System_seps[time_key][start_ind:end_ind])[non_nan_inds], axis=0)
                        mean_time.append(mean_t)
                        mean_sep.append(mean_s)
                axs.flatten()[pick_it].semilogy(mean_time, mean_sep, color='k', linewidth=0.25)
            '''
            #axs.flatten()[pick_it].semilogy(Time_adjusted_formation, System_semimajor[time_key], alpha=0.5, color='k')
        print("Systems with lifetime greater than", args.sys_lifetime_threshold, "is", np.sum(np.array(list(System_lifetimes.values()))>args.sys_lifetime_threshold))

        if plot_booleans[rank][0] == True:
            for bin_bound in S_bins:
                axs.flatten()[pick_it].axhline(y=bin_bound, linewidth=0.5)
        #if np.remainder(int(pick_it/2), 2) == 0:
        axs.flatten()[pick_it].set_ylabel('Separation (AU)')
        #else:
        #    yticklabels = axs.flatten()[pick_it].get_yticklabels()
        #    plt.setp(yticklabels, visible=False)
        axs.flatten()[pick_it].set_ylim([10, 10000])
        if np.isnan(t_max):
            t_max = G50_t_max
        axs.flatten()[pick_it].set_xlim([0, t_max])
        #axs.flatten()[pick_it].set_xlim([Times[0], Times[-1]])
        axs.flatten()[pick_it].add_patch(matplotlib.patches.Rectangle((0.85*t_max, 3500), 400000, 5600, facecolor="white", edgecolor="black", zorder=10))
        axs.flatten()[pick_it].text((0.86*t_max), 6500, subplot_titles[int(pick_it/2)], zorder=11)
        axs.flatten()[pick_it].text((0.86*t_max), 4000, "N$_{sys}$="+str(np.sum(np.array(list(System_lifetimes.values()))>args.sys_lifetime_threshold)), zorder=12)
        
        if args.add_histograms == "True":
            CF_norm = np.array(CF_hist[int(pick_it/2)])/np.sum(CF_hist[int(pick_it/2)])
            axs.flatten()[pick_it+1].step([CF_norm[0]] + CF_norm.tolist(), np.log10(S_bins), linewidth=2)
            axs.flatten()[pick_it+1].set_ylim([1,4])
            yticklabels =axs.flatten()[pick_it+1].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        
        plt.savefig('superplot_multi'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
        #plt.savefig('superplot_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
    
    if args.timescale == "first_sink":
        axs.flatten()[pick_it].set_xlabel('Time since formation of first component ($yr$)')
    else:
        axs.flatten()[pick_it].set_xlabel('Time since formation of first multiple star system ($yr$)')
    if args.add_histograms == "True":
        axs.flatten()[pick_it+1].set_xlabel('% of systems')
    plt.savefig('superplot_multi'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
    plt.savefig('superplot_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
    print('Created superplot_multi.jpg')

#================================================================================================

#Plot simulation SFE eovlution
plt.clf()
colors = ['b', 'orange', 'g', 'r', 'c', 'm']
fig, axs = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
plt.subplots_adjust(wspace=0.0)
SFE_5_inds = []
Total_Mass_5 = []
N_Stars_5 = []
for pick_it in range(len(pickle_files)):
    file = open(pickle_files[pick_it], 'rb')
    Times, SFE, SFE_n, M_tot, M_tot_vis, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_ecc, Sink_formation_times, System_mean_times, System_mean_seps, System_mean_ecc, System_lifetimes, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()
    
    SFE_5_ind = np.argmin(abs(np.array(SFE)-0.05))
    SFE_5_inds.append(SFE_5_ind)
    Total_Mass_5.append(M_tot[SFE_5_ind])
    N_Stars_5.append(N_stars[SFE_5_ind])
    
    axs.flatten()[0].plot(Times, M_tot, label=subplot_titles[pick_it], color=colors[pick_it])
    axs.flatten()[0].scatter(Times[SFE_5_ind], M_tot[SFE_5_ind], color=colors[pick_it])
    axs.flatten()[0].plot(Times, M_tot_multi, color=colors[pick_it], linestyle="--", alpha=0.5)
    axs.flatten()[0].scatter(Times[SFE_5_ind], M_tot_multi[SFE_5_ind], color=colors[pick_it], alpha=0.5)
    
    axs.flatten()[1].plot(Times, N_stars, color=colors[pick_it])
    axs.flatten()[1].scatter(Times[SFE_5_ind], N_stars[SFE_5_ind], color=colors[pick_it])
    axs.flatten()[1].plot(Times, N_multi_stars, color=colors[pick_it], linestyle="--", alpha=0.5)
    axs.flatten()[1].scatter(Times[SFE_5_ind], N_multi_stars[SFE_5_ind], color=colors[pick_it], alpha=0.5)
    
    axs.flatten()[2].plot(Times, SFE, color=colors[pick_it])
    axs.flatten()[2].scatter(Times[SFE_5_ind], SFE[SFE_5_ind], color=colors[pick_it])
axs.flatten()[0].legend(loc='best')
plt.xlabel("Time in Simulation")
#plt.xlabel("Time since the first star formed")
axs.flatten()[0].set_ylabel("Total acccreted mass ($M_\\odot$)")
axs.flatten()[1].set_ylabel("Total number of stars ($M_\\odot$)")
axs.flatten()[2].set_ylabel("Star formation efficiency (SFE)")
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylim(bottom=0)
#axs.flatten()[0].set_xlim(left=0)
plt.savefig('SFE_evolution.jpg', format='jpg', bbox_inches='tight')
plt.savefig('SFE_evolution.pdf', format='pdf', bbox_inches='tight')
print('Saved SFE_evolution')

#plot quantities VS cloud mass at SFE=0.05
plt.clf()
fig, axs = plt.subplots(1, 1, figsize=(4, 3), sharex=True)
plt.plot(GMC_mass_arr, N_Stars_5)
plt.scatter(GMC_mass_arr, N_Stars_5)
plt.xlabel('GMC Mass (M$_\odot$)')
plt.ylabel('# stars')
plt.xlim([GMC_mass_arr[0], GMC_mass_arr[-1]])
plt.savefig('N_stars_vs_GMC_mass.jpg', format='jpg', bbox_inches='tight')
plt.savefig('N_stars_vs_GMC_mass.pdf', format='pdf', bbox_inches='tight')
print('Saved N_stars_vs_GMC_mass.jpg')

plot_truncated_super_mult = True
if plot_truncated_super_mult == True:
    plt.clf()
    if plot_booleans[rank][0] == True:
        fig, axs = plt.subplots(len(pickle_files), 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, len(pickle_files)*3), sharex='col')
        iter_range = range(0, len(pickle_files)*2, 2)
        args.figure_suffix = args.figure_suffix + '_hist'
    else:
        fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(12, len(pickle_files)*3))
        iter_range = range(0, len(pickle_files)*2)
    plt.subplots_adjust(wspace=0.0)
    plt.subplots_adjust(hspace=0.07)

    t_max = np.nan
    G50_t_max = 0
    CF_hist = np.zeros((len(pickle_files),12)).tolist()
    
    for pick_it in iter_range:
        file = open(pickle_files[int(pick_it/2)], 'rb')
        Times, SFE, SFE_n, M_tot, M_tot_vis, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_ecc, Sink_formation_times, System_mean_times, System_mean_seps, System_mean_ecc, System_lifetimes, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
        file.close()

        SFE_5_ind = SFE_5_inds[int(pick_it/2)]

        sim_start_time = np.nan
        for time_key in System_times.keys():
            for sit in range(len(S_bins[1:])):
                bin_inst = np.argwhere((System_seps[time_key]>S_bins[sit])&(System_seps[time_key]<S_bins[sit+1]))
                CF_hist[int(pick_it/2)][sit] = CF_hist[int(pick_it/2)][sit] + np.shape(bin_inst)[0]
            key_inds = flatten(eval(time_key))
            if args.timescale == "first_sink":
                first_sink = np.argmin(Sink_formation_times[key_inds])
                Time_adjusted_formation = System_times[time_key] - Sink_formation_times[key_inds[first_sink]].value
            else:
                if np.isnan(sim_start_time):
                    sim_start_time = System_times[time_key][0]
                Time_adjusted_formation = System_times[time_key] - np.array(sim_start_time)
            if np.isnan(t_max):
                if Time_adjusted_formation[-1] > G50_t_max:
                    G50_t_max = Time_adjusted_formation[-1]
            axs.flatten()[pick_it].semilogy(Time_adjusted_formation[:SFE_5_ind], System_seps[time_key][:SFE_5_ind], alpha=0.05, color='k')
            
        for bin_bound in S_bins:
            axs.flatten()[pick_it].axhline(y=bin_bound, linewidth=0.5)

        axs.flatten()[pick_it].set_ylabel('Separation (AU)')
        axs.flatten()[pick_it].set_ylim([10, 10000])
        if np.isnan(t_max):
            t_max = G50_t_max
        axs.flatten()[pick_it].set_xlim([0, t_max])
        
        axs.flatten()[pick_it].add_patch(matplotlib.patches.Rectangle((0.85*t_max, 3500), 400000, 5600, facecolor="white", edgecolor="black", zorder=10))
        axs.flatten()[pick_it].text((0.86*t_max), 6500, subplot_titles[int(pick_it/2)], zorder=11)
        axs.flatten()[pick_it].text((0.86*t_max), 4000, "N$_{sys}$="+str(np.sum(np.array(list(System_lifetimes.values())[:SFE_5_ind])>args.sys_lifetime_threshold)), zorder=12)
        
        if args.add_histograms == "True":
            CF_norm = np.array(CF_hist[int(pick_it/2)])/np.sum(CF_hist[int(pick_it/2)])
            axs.flatten()[pick_it+1].step([CF_norm[0]] + CF_norm.tolist(), np.log10(S_bins), linewidth=2)
            axs.flatten()[pick_it+1].set_ylim([1,4])
            yticklabels =axs.flatten()[pick_it+1].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        
        plt.savefig('superplot_multi_truncated'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
        print('plotted separations for pickle', pickle_files[int(pick_it/2)])
        #plt.savefig('superplot_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
    
    if args.timescale == "first_sink":
        axs.flatten()[pick_it].set_xlabel('Time since formation of first component ($yr$)')
    else:
        axs.flatten()[pick_it].set_xlabel('Time since formation of first multiple star system ($yr$)')
    if args.add_histograms == "True":
        axs.flatten()[pick_it+1].set_xlabel('% of systems')
    plt.savefig('superplot_multi_truncated'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
    plt.savefig('superplot_multi_truncated'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
    print('Created superplot_multi_truncated'+args.figure_suffix+'.jpg')











'''
plt.clf()
colors = ['b', 'orange', 'g', 'r', 'c', 'm']
fig, axs = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
plt.subplots_adjust(wspace=0.0)
for time_it in range(len(Times_all)):
    axs.flatten()[0].plot(Times_all[time_it] - np.array(Times_all[time_it][0]), M_tot_all[time_it], label=subplot_titles[time_it], color=colors[time_it])
    axs.flatten()[0].plot(Times_all[time_it] - np.array(Times_all[time_it][0]), M_tot_multi_all[time_it], color=colors[time_it], linestyle="--", alpha=0.5)
    axs.flatten()[1].plot(Times_all[time_it] - np.array(Times_all[time_it][0]), N_stars_all[time_it], color=colors[time_it])
    axs.flatten()[1].plot(Times_all[time_it] - np.array(Times_all[time_it][0]), N_multi_stars_all[time_it], color=colors[time_it], linestyle="--", alpha=0.5)
    axs.flatten()[2].plot(Times_all[time_it] - np.array(Times_all[time_it][0]), SFE_all[time_it], color=colors[time_it])
axs.flatten()[0].legend(loc='best')
#plt.xlabel("Time in Simulation")
plt.xlabel("Time since the first star formed")
axs.flatten()[0].set_ylabel("Total acccreted mass ($M_\\odot$)")
axs.flatten()[1].set_ylabel("Total number of stars ($M_\\odot$)")
axs.flatten()[2].set_ylabel("Star formation efficiency (SFE)")
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[0].set_xlim(left=0)
plt.savefig('SFE_evolution.jpg', format='jpg', bbox_inches='tight')
plt.savefig('SFE_evolution.pdf', format='pdf', bbox_inches='tight')

plt.clf()
fig = plt.figure(figsize=(8, 6))
linestyles = ["-", "--", "-.", (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1, 1, 1)), ":"]
for shrink_it in range(len(Shrinkage_hist)):
    plt.step(np.linspace(-10, 100, 12), [Shrinkage_hist[shrink_it][0]] + Shrinkage_hist[shrink_it].tolist(), label=subplot_titles[shrink_it], linewidth=2, ls=linestyles[shrink_it])
plt.xlabel('Orbital Shrinkage (% of inital separation)')
plt.ylabel('% of systems')
plt.legend(loc='best')
plt.ylim(bottom=0.0)
plt.xlim([-10, 100])
plt.savefig('shrinkage_multi'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
plt.savefig('shrinkage_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
print('created shrinkage_multi.jpg')

plt.clf()
plt.plot(GMC_mass_arr, np.array(Shrinkage_hist)[:,0], label='ejections')
plt.plot(GMC_mass_arr, np.array(Shrinkage_hist)[:,-1], label='strong inspiral')
plt.xlabel('GMC cloud mass (M$_\odot$)')
plt.ylabel('Fraction of systems')
plt.legend(loc='best')
plt.xlim([GMC_mass_arr[0], GMC_mass_arr[-1]])
plt.savefig('ejection_inspiral'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
plt.savefig('ejection_inspiral'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
print('created ejection_inspiral.jpg')


plt.clf()
fig = plt.figure(figsize=(8, 6))
linestyles = ["-", "--", "-.", (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1, 1, 1)), ":"]
for CF_it in range(len(CF_hist)):
    CF_norm = np.array(CF_hist[CF_it])/np.sum(CF_hist[CF_it])
    plt.step(np.log10(S_bins), [CF_norm[0]] + CF_norm.tolist(), label=subplot_titles[CF_it], linewidth=2, ls=linestyles[CF_it])
plt.xlabel('Log(Separation) (AU)')
plt.ylabel('% of systems')
plt.legend(loc='best')
plt.ylim(bottom=0.0)
plt.xlim([1, 4])
plt.savefig('separation_hist_multi'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
plt.savefig('separation_hist_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
print('created separation_hist_multi.jpg')

plt.clf()
fig = plt.figure(figsize=(8, 6))
linestyles = ["-", "--", "-.", (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1, 1, 1)), ":"]
for life_it in range(len(Lifetime_hist)):
    plt.step([0,1,2,3,4,5,6,7], [Lifetime_hist[life_it][0]] + Lifetime_hist[life_it].tolist(), label=subplot_titles[life_it], linewidth=2, ls=linestyles[life_it])
plt.xlabel('Log(Lifetime) (yr)')
plt.ylabel('% of systems')
plt.legend(loc='best')
plt.ylim(bottom=0.0)
plt.xlim([-10, 100])
plt.savefig('lifetime_multi'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight')
plt.savefig('lifetime_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
print('created lifetime_multi.jpg')
'''
