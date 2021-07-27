import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import pickle
import matplotlib.patches

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifeimteimte threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

pickle_files = ["/groups/astro/rlk/Analysis_plots/Ramses/Superplots/Cyclic/G50/Replace_with_most_massive_star/superplot_with_means.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Superplots/Cyclic/G100/256/Replace_with_most_massive_star/superplot_with_means.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Superplots/Cyclic/G200/Replace_with_most_massive_star/superplot_with_means.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Superplots/Cyclic/G400/Replace_with_most_massive_star/superplot_with_means.pkl"]

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
plt.subplots_adjust(wspace=0.0)

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
Shrinkage_hist = []
Shrinkage_bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
Shrinks = []
Lifetime_bins = [0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
Lifetime_hist = []
S_bins = np.logspace(1,4,13)
CF_hist = [0,0,0,0,0,0,0,0,0,0,0,0]

for pick_it in range(len(pickle_files)):
    file = open(pickle_files[pick_it], 'rb')
    Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_mean_times, System_mean_seps, System_lifetimes, Initial_Seps, Final_seps, Sep_maxs, Sep_mins = pickle.load(file)
    file.close()
    
    shrink = 100*((np.array(Initial_Seps) - np.array(Final_seps))/np.array(Initial_Seps))
    Shrinks.append(shrink)
    shrink_hist, bins = np.histogram(shrink, bins=Shrinkage_bins)
    Shrinkage_hist.append(shrink_hist/np.sum(shrink_hist))
    life_hist, bins = np.histogram(System_lifetimes, bins=Lifetime_bins)
    Lifetime_hist.append(life_hist/np.sum(life_hist))
    
    lifetime_it = -1
    for time_key in System_times.keys():
        for sit in range(len(S_bins[1:])):
            bin_inst = np.argwhere((System_seps[time_key]>S_bins[sit])&(System_seps[time_key]<S_bins[sit+1]))
            CF_hist[sit] = CF_hist[sit] + np.shape(bin_inst)[0]
        axs.flatten()[pick_it].semilogy(System_times[time_key], System_seps[time_key], alpha=0.05, color='k')
        lifetime_it = lifetime_it + 1
        if System_lifetimes[lifetime_it] > args.sys_lifetime_threshold:
            axs.flatten()[pick_it].semilogy(System_mean_times[time_key], System_mean_seps[time_key], color='k', linewidth=0.5)
    print("Systems with lifetime greater than", args.sys_lifetime_threshold, "is", np.sum(np.array(System_lifetimes)>args.sys_lifetime_threshold))

    for bin_bound in S_bins:
        axs.flatten()[pick_it].axhline(y=bin_bound, linewidth=0.5)
    if np.remainder(pick_it, 2) == 0:
        axs.flatten()[pick_it].set_ylabel('Separation (AU)')
    else:
        yticklabels = axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels, visible=False)
    axs.flatten()[pick_it].set_xlabel('Time ($yr$)')
    axs.flatten()[pick_it].set_ylim([10, 10000])
    axs.flatten()[pick_it].set_xlim([Times[0], Times[-1]])
    axs.flatten()[pick_it].add_patch(matplotlib.patches.Rectangle(((0.02*(Times[-1]-Times[0])+Times[0]), 3500), 0.23*(Times[-1]-Times[0]), 5500, facecolor="white", edgecolor="black", zorder=10))
    axs.flatten()[pick_it].text((0.03*(Times[-1]-Times[0])+Times[0]), 6500, subplot_titles[pick_it], zorder=11)
    axs.flatten()[pick_it].text((0.03*(Times[-1]-Times[0])+Times[0]), 4000, "N$_{sys}$="+str(np.sum(np.array(System_lifetimes)>args.sys_lifetime_threshold)), zorder=12)
    
    plt.savefig('superplot_multi.jpg', format='jpg', bbox_inches='tight')
    plt.savefig('superplot_multi.pdf', format='pdf', bbox_inches='tight')
print('Created superplot_multi.jpg')


plt.clf()
linestyles = ["-", "--", "-.", ":"]
for shrink_it in range(len(Shrinkage_hist)):
    plt.step(np.linspace(-10, 100, 12), [Shrinkage_hist[shrink_it][0]] + Shrinkage_hist[shrink_it].tolist(), label=subplot_titles[shrink_it], linewidth=2, ls=linestyles[shrink_it])
plt.xlabel('Orbital Shrinkage (% of inital separation)')
plt.ylabel('% of systems')
plt.legend(loc='best')
plt.ylim(bottom=0.0)
plt.xlim([-10, 100])
plt.savefig('shrinkage_multi.jpg', format='jpg', bbox_inches='tight')
plt.savefig('shrinkage_multi.pdf', format='pdf', bbox_inches='tight')
print('created shrinkage_multi.jpg')

plt.clf()
linestyles = ["-", "--", "-.", ":"]
for shrink_it in range(len(Lifetime_hist)):
    plt.step([0,1,2,3,4,5,6,7], [Lifetime_hist[shrink_it][0]] + Lifetime_hist[shrink_it].tolist(), label=subplot_titles[shrink_it], linewidth=2, ls=linestyles[shrink_it])
plt.xlabel('Log(Lifetime) (yr)')
plt.ylabel('% of systems')
plt.legend(loc='best')
plt.ylim(bottom=0.0)
plt.xlim([-10, 100])
plt.savefig('lifetime_multi.jpg', format='jpg', bbox_inches='tight')
plt.savefig('lifetime_multi.pdf', format='pdf', bbox_inches='tight')
print('created lifetime_multi.jpg')

