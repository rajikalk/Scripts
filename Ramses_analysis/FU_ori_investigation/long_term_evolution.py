import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt
import sys


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_id", help="which sink?", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}" "\sisetup{detect-all}" r"\usepackage{helvet}" r"\usepackage{sansmath}" "\sansmath"               # <- tricky! -- gotta actually tell tex to use!

args = parse_inputs()
save_dir = sys.argv[1]
if args.sink_id == None:
    global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global.pkl"
else:
    global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_"+str(args.sink_id)+".pkl"
print('global pickle:', global_pickle)
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()
print('successfully read global pickle')

if 'ltot' not in particle_data.keys() or np.shape(particle_data['ltot'])[1] == 2:
    print('Calculating Total Luminosity')
    Baraffe_mass = yt.YTArray([0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.072, 0.075, 0.080, 0.090, 0.100, 0.110, 0.130, 0.150, 0.170, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400], 'msun')
    Baraffe_logL = np.array([-2.469, -2.208, -2.044, -1.783, -1.655, -1.481, -1.399, -1.324, -1.291, -1.261, -1.197, -1.127, -1.154, -1.075, -0.926, -0.795, -0.669, -0.539, -0.199, -0.040, 0.076, 0.171, 0.268, 0.356, 0.436, 0.508, 0.573, 0.634, 0.688, 0.740])
    Baraffe_radius = yt.YTArray([0.341, 0.416, 0.472, 0.603, 0.665, 0.796, 0.846, 0.905, 0.942, 0.972, 1.045, 1.113, 1.033, 1.115, 1.270, 1.412, 1.568, 1.731, 2.215, 2.364, 2.458, 2.552, 2.687, 2.821, 2.960, 3.096, 3.227, 3.362, 3.488, 3.621], 'rsun')

    #Derive a stellar luminosity
    print('Calculating luminosity for candidate')
    facc = 0.5
    lstar_baraffe_prim = []
    rstar_barrafe_prim = []
    Mass_prim = yt.YTArray(particle_data['mass']).T[0]
    Mdot_prim = yt.YTArray(particle_data['mdot']).T[0]
    for mass_val in Mass_prim:
        if mass_val < Baraffe_mass[0]:
            lstar_baraffe_prim.append(10**Baraffe_logL[0])
            rstar_barrafe_prim.append(Baraffe_radius[0])
        else:
            closest_inds = sorted(np.argsort(np.abs(Baraffe_mass - mass_val))[:2])
            gradient = (Baraffe_logL[closest_inds][1] - Baraffe_logL[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
            y_intercept = Baraffe_logL[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
            logL = gradient*mass_val + y_intercept
            lstar_baraffe_prim.append(10**logL)
            
            gradient = (Baraffe_radius[closest_inds][1] - Baraffe_radius[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
            y_intercept = Baraffe_radius[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
            radius = gradient*mass_val + y_intercept
            rstar_barrafe_prim.append(radius)

    lacc_prim = facc * (yt.units.gravitational_constant_cgs * Mass_prim.in_units('g') * Mdot_prim.in_units('g/s'))/yt.YTArray(rstar_barrafe_prim).in_units('cm')
    ltot_prim = lacc_prim.in_units('lsun') + yt.YTArray(np.array(lstar_baraffe_prim), 'lsun')
    del lstar_baraffe_prim, rstar_barrafe_prim, Mass_prim, Mdot_prim, lacc_prim
    print('Candidate luminosity finished')
    
    print('Calculating luminosity for companion')
    lstar_baraffe_comp = []
    rstar_barrafe_comp = []
    Mass_comp = yt.YTArray(particle_data['mass']).T[1]
    Mdot_comp = yt.YTArray(particle_data['mdot']).T[1]
    for mass_val in Mass_comp:
        if mass_val < Baraffe_mass[0]:
            lstar_baraffe_comp.append(10**Baraffe_logL[0])
            rstar_barrafe_comp.append(Baraffe_radius[0])
        else:
            closest_inds = sorted(np.argsort(np.abs(Baraffe_mass - mass_val))[:2])
            gradient = (Baraffe_logL[closest_inds][1] - Baraffe_logL[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
            y_intercept = Baraffe_logL[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
            logL = gradient*mass_val + y_intercept
            lstar_baraffe_comp.append(10**logL)
            
            gradient = (Baraffe_radius[closest_inds][1] - Baraffe_radius[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
            y_intercept = Baraffe_radius[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
            radius = gradient*mass_val + y_intercept
            rstar_barrafe_comp.append(radius)

    lacc_comp = facc * (yt.units.gravitational_constant_cgs * Mass_comp.in_units('g') * Mdot_comp.in_units('g/s'))/yt.YTArray(rstar_barrafe_comp).in_units('cm')
    ltot_comp = lacc_comp.in_units('lsun') + yt.YTArray(np.array(lstar_baraffe_comp), 'lsun')
    del lstar_baraffe_comp, rstar_barrafe_comp, Mass_comp, Mdot_comp, lacc_comp
    print('Companion luminosity finished')

    particle_data.update({'ltot':yt.YTArray([ltot_prim, ltot_comp]).T})

    print('updating pickle')
    file = open(global_pickle, 'wb')
    pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
    file.close()
    print('Finished updating pickle', global_pickle)
print('Got total luminosity')

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(two_col_width, 1.5*two_col_width), sharey=True)#, sharey=True)
#plt.subplots_adjust(wspace=0.0)
#plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].semilogy(particle_data['time'], yt.YTArray(particle_data['ltot']).T[0])
for closest_id in np.unique(particle_data['closest_sink'], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'], return_index=True)[1])]:
    curr_inds = np.argwhere(np.array(particle_data['closest_sink']) == closest_id).T[0]
    diff_inds = np.setdiff1d(np.arange(len(particle_data['time'])), curr_inds)
    ltot_curr = np.copy(yt.YTArray(particle_data['ltot']).T[1])
    ltot_curr[diff_inds] = np.nan
    axs.flatten()[0].semilogy(yt.YTArray(particle_data['time']), ltot_curr, ls=':')
axs.flatten()[0].set_xlim([0, 10000])
axs.flatten()[0].set_ylabel("L$_{tot}$ (L$_\odot$)")
axs.flatten()[0].tick_params(axis='both', direction='in', top=True)
ax0 = axs.flatten()[0].twinx()
ax0.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax0.set_ylabel('Separation (AU)')
ax0.set_ylim([5,1000])
ax0.tick_params(axis='both', direction='in', top=True)
print('plotted time [0, 10000]')

axs.flatten()[1].semilogy(particle_data['time'], yt.YTArray(particle_data['ltot']).T[0])
for closest_id in np.unique(particle_data['closest_sink'], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'], return_index=True)[1])]:
    curr_inds = np.argwhere(np.array(particle_data['closest_sink']) == closest_id).T[0]
    diff_inds = np.setdiff1d(np.arange(len(particle_data['time'])), curr_inds)
    ltot_curr = np.copy(yt.YTArray(particle_data['ltot']).T[1])
    ltot_curr[diff_inds] = np.nan
    axs.flatten()[1].semilogy(yt.YTArray(particle_data['time']), ltot_curr, ls=':')
axs.flatten()[1].set_xlim([10000, 20000])
axs.flatten()[1].set_ylabel("L$_{tot}$ (L$_\odot$)")
axs.flatten()[1].tick_params(axis='both', direction='in', top=True)
ax1 = axs.flatten()[1].twinx()
ax1.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax1.set_ylabel('Separation (AU)')
ax1.set_ylim([5,1000])
ax1.tick_params(axis='both', direction='in', top=True)
print('plotted time [10000, 20000]')

axs.flatten()[2].semilogy(particle_data['time'], yt.YTArray(particle_data['ltot']).T[0])
for closest_id in np.unique(particle_data['closest_sink'], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'], return_index=True)[1])]:
    curr_inds = np.argwhere(np.array(particle_data['closest_sink']) == closest_id).T[0]
    diff_inds = np.setdiff1d(np.arange(len(particle_data['time'])), curr_inds)
    ltot_curr = np.copy(yt.YTArray(particle_data['ltot']).T[1])
    ltot_curr[diff_inds] = np.nan
    axs.flatten()[2].semilogy(yt.YTArray(particle_data['time']), ltot_curr, ls=':')
axs.flatten()[2].set_xlim([20000, 30000])
axs.flatten()[2].set_ylabel("L$_{tot}$ (L$_\odot$)")
axs.flatten()[2].tick_params(axis='both', direction='in', top=True)
ax2 = axs.flatten()[2].twinx()
ax2.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax2.set_ylabel('Separation (AU)')
ax2.set_ylim([5,1000])
ax2.tick_params(axis='both', direction='in', top=True)
print('plotted time [20000, 30000]')

axs.flatten()[3].semilogy(particle_data['time'], yt.YTArray(particle_data['ltot']).T[0])
for closest_id in np.unique(particle_data['closest_sink'], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'], return_index=True)[1])]:
    curr_inds = np.argwhere(np.array(particle_data['closest_sink']) == closest_id).T[0]
    diff_inds = np.setdiff1d(np.arange(len(particle_data['time'])), curr_inds)
    ltot_curr = np.copy(yt.YTArray(particle_data['ltot']).T[1])
    ltot_curr[diff_inds] = np.nan
    axs.flatten()[3].semilogy(yt.YTArray(particle_data['time']), ltot_curr, ls=':')
axs.flatten()[3].set_xlim([30000, 40000])
axs.flatten()[3].set_ylabel("L$_{tot}$ (L$_\odot$)")
axs.flatten()[3].tick_params(axis='both', direction='in', top=True)
ax3 = axs.flatten()[3].twinx()
ax3.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax3.set_ylabel('Separation (AU)')
ax3.set_ylim([5,1000])
ax3.tick_params(axis='both', direction='in', top=True)
print('plotted time [30000, 40000]')

axs.flatten()[4].semilogy(particle_data['time'], yt.YTArray(particle_data['ltot']).T[0])
for closest_id in np.unique(particle_data['closest_sink'], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'], return_index=True)[1])]:
    curr_inds = np.argwhere(np.array(particle_data['closest_sink']) == closest_id).T[0]
    diff_inds = np.setdiff1d(np.arange(len(particle_data['time'])), curr_inds)
    ltot_curr = np.copy(yt.YTArray(particle_data['ltot']).T[1])
    ltot_curr[diff_inds] = np.nan
    axs.flatten()[4].semilogy(yt.YTArray(particle_data['time']), ltot_curr, ls=':')
axs.flatten()[4].set_xlim([40000, 50000])
axs.flatten()[4].set_ylabel("L$_{tot}$ (L$_\odot$)")
axs.flatten()[4].tick_params(axis='both', direction='in', top=True)
ax4 = axs.flatten()[4].twinx()
ax4.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax4.set_ylabel('Separation (AU)')
ax4.set_ylim([5,1000])
ax4.tick_params(axis='both', direction='in', top=True)
axs.flatten()[4].set_xlabel("Time (yr)")
#axs.flatten()[4].set_ylim([5.e-2, 2.5e1])
print('plotted time [40000, 50000]')

plt.savefig('long_term_evolution_ltot_'+str(args.sink_id)+'.pdf', bbox_inches='tight', pad_inches=0.02)

