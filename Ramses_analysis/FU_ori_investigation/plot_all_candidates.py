import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt

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

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global_all_candidates.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_form_time = pickle.load(file_open)
file_open.close()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

for tag in particle_data['tag']:
    tag_it = particle_data['tag'].index(tag)
    
    Baraffe_mass = yt.YTArray([0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.072, 0.075, 0.080, 0.090, 0.100, 0.110, 0.130, 0.150, 0.170, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400], 'msun')
    Baraffe_logL = np.array([-2.469, -2.208, -2.044, -1.783, -1.655, -1.481, -1.399, -1.324, -1.291, -1.261, -1.197, -1.127, -1.154, -1.075, -0.926, -0.795, -0.669, -0.539, -0.199, -0.040, 0.076, 0.171, 0.268, 0.356, 0.436, 0.508, 0.573, 0.634, 0.688, 0.740])
    Baraffe_radius = yt.YTArray([0.341, 0.416, 0.472, 0.603, 0.665, 0.796, 0.846, 0.905, 0.942, 0.972, 1.045, 1.113, 1.033, 1.115, 1.270, 1.412, 1.568, 1.731, 2.215, 2.364, 2.458, 2.552, 2.687, 2.821, 2.960, 3.096, 3.227, 3.362, 3.488, 3.621], 'rsun')

    #Derive a stellar luminosity
    facc = 0.5
    lstar_baraffe = []
    rstar_barrafe = []
    Mass_prim = particle_data['mass'][tag_it]
    Mdot_prim = particle_data['mdot'][tag_it]
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
    ltot = lacc_prim.in_units('lsun') + yt.YTArray(np.array(lstar_baraffe_prim), 'lsun')
    
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(two_col_width, 1.5*two_col_width), sharey=True)
    right_ax = axs.twinx()
    
    axs.flatten()[0].semilogy(particle_data['time'][tag_it], ltot)
    axs.flatten()[0].set_xlim([0, 10000])
    axs.flatten()[0].set_ylabel("L$_{acc}$ (L$_\odot$/)")
    axs.flatten()[0].tick_params(axis='both', direction='in', top=True)
    ax0 = axs.flatten()[0].twinx()
    ax0.semilogy(particle_data['time'][tag_it], particle_data['separation'][tag_it], color='k', ls="--", alpha=0.25)
    ax0.set_ylabel('Separation (AU)')
    #ax0.set_ylim([5,1000])
    ax0.tick_params(axis='both', direction='in', top=True)

    axs.flatten()[1].semilogy(particle_data['time'][tag_it], ltot)
    axs.flatten()[1].set_xlim([10000, 20000])
    axs.flatten()[1].set_ylabel("L$_{acc}$ (L$_\odot$/)")
    axs.flatten()[1].tick_params(axis='both', direction='in', top=True)
    ax1 = axs.flatten()[1].twinx()
    ax1.semilogy(particle_data['time'][tag_it], particle_data['separation'][tag_it], color='k', ls="--", alpha=0.25)
    ax1.set_ylabel('Separation (AU)')
    #ax1.set_ylim([5,1000])
    ax1.tick_params(axis='both', direction='in', top=True)

    axs.flatten()[2].semilogy(particle_data['time'][tag_it], ltot)
    axs.flatten()[2].set_xlim([20000, 30000])
    axs.flatten()[2].set_ylabel("L$_{acc}$ (L$_\odot$/)")
    axs.flatten()[2].tick_params(axis='both', direction='in', top=True)
    ax2 = axs.flatten()[2].twinx()
    ax2.semilogy(particle_data['time'][tag_it], particle_data['separation'][tag_it], color='k', ls="--", alpha=0.25)
    ax2.set_ylabel('Separation (AU)')
    #ax2.set_ylim([5,1000])
    ax2.tick_params(axis='both', direction='in', top=True)

    axs.flatten()[3].semilogy(particle_data['time'][tag_it], ltot)
    axs.flatten()[3].set_xlim([30000, 40000])
    axs.flatten()[3].set_ylabel("L$_{acc}$ (L$_\odot$/)")
    axs.flatten()[3].tick_params(axis='both', direction='in', top=True)
    ax3 = axs.flatten()[3].twinx()
    ax3.semilogy(particle_data['time'][tag_it], particle_data['separation'][tag_it], color='k', ls="--", alpha=0.25)
    ax3.set_ylabel('Separation (AU)')
    #ax3.set_ylim([5,1000])
    ax3.tick_params(axis='both', direction='in', top=True)

    axs.flatten()[4].semilogy(particle_data['time'][tag_it], ltot)
    axs.flatten()[4].set_xlim([40000, 50000])
    axs.flatten()[4].set_ylabel("L$_{acc}$ (L$_\odot$/)")
    axs.flatten()[4].tick_params(axis='both', direction='in', top=True)
    ax4 = axs.flatten()[4].twinx()
    ax4.semilogy(particle_data['time'][tag_it], particle_data['separation'][tag_it], color='k', ls="--", alpha=0.25)
    ax4.set_ylabel('Separation (AU)')
    #ax4.set_ylim([5,1000])
    ax4.tick_params(axis='both', direction='in', top=True)
    axs.flatten()[4].set_xlabel("Time (yr)")
   # axs.flatten()[4].set_ylim([1.e-1, 2.e1])
    
    plt.savefig('Sink_'+str(tag)+'_accretion_evol.pdf', bbox_inches='tight', pad_inches=0.02)
    print('plotted accretion evolution for sink', tag)

