import numpy as np
import matplotlib.pyplot as plt
import pickle

Dens = [50, 100, 100, 200, 200, 400]
pickles = ['/groups/astro/rlk/Analysis_plots/Ramses/Superplots/G50/Replace_with_most_massive_star/superplot_with_means.pkl', '/groups/astro/rlk/Analysis_plots/Ramses/Superplots/G100/256/Replace_with_most_massive_star/superplot_with_means.pkl', '/groups/astro/rlk/Analysis_plots/Ramses/Superplots/G100/512/Replace_with_most_massive_star/superplot_with_means.pkl', '/groups/astro/rlk/Analysis_plots/Ramses/Superplots/G200/Replace_with_most_massive_star/superplot_with_means.pkl', '/groups/astro/rlk/Analysis_plots/Ramses/Superplots/G200_new/Replace_with_most_massive_star/superplot_with_means.pkl', '/groups/astro/rlk/Analysis_plots/Ramses/Superplots/G400/Replace_with_most_massive_star/superplot_with_means.pkl']

Vis_stars = []
Ratios = []

for pickle_file in pickles:
    file = open(pickle_file, 'rb')
    Times, SFE, SFE_n, M_tot, M_tot_multi, N_stars, N_vis_stars, N_multi_stars, System_seps, System_semimajor, System_times, System_mean_times, System_mean_seps, System_lifetimes, Initial_Seps, Final_seps, Sep_maxs, Sep_mins, Sink_formation_times = pickle.load(file)
    file.close()
    
    end_time = Times[-1] - 10000
    end_time_it = np.argmin(abs(Times - np.array(end_time)))
    vis_end = np.mean(N_vis_stars[end_time_it:])
    ratio_end = vis_end/np.mean(N_stars[end_time_it:])
    Vis_stars.append(vis_end)
    Ratios.append(ratio_end)


m, b = np.polyfit(Dens, Vis_stars, 1)
x = np.arange(0,200)
y = m*x+b
g_match = np.argmin(abs(y - 115))

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2)
plt.subplots_adjust(hspace=0.0)
axs[0].scatter(Dens, Vis_stars, label="visible stars")
axs[0].plot(Dens, m*np.array(Dens)+b)
axs[0].set_ylabel("visible stars")
axs[0].set_ylim(bottom=0)
axs[0].axhline(y=115)
axs[0].axvline(x=g_match)
axs[1].scatter(Dens, Ratios, label="ratio vis to all stars")
axs[1].set_ylabel("Ratio of vis to all stars")
plt.xlabel("Density")
plt.savefig("vis_vs_dens.png")
