import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt
from matplotlib.ticker import FormatStrFormatter

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

FU_temp = np.concatenate((np.zeros(15), np.ones(85)))
FU_temp_inv = np.concatenate((np.ones(15), np.zeros(85)))
#FU_temp = np.concatenate((np.zeros(7), np.ones(43)))
#FU_temp_inv = np.concatenate((np.ones(7), np.zeros(43)))

time_window = yt.YTQuantity(100, 'yr')

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()

age = yt.YTArray(particle_data['time'])
mass = yt.YTArray(particle_data['mass']).T[1]
mdot = yt.YTArray(particle_data['mdot']).T[1]
facc = 0.5

#M/Ms      log t(yr)    Teff     L/Ls    g    R/Rs  Log(Li/Li0) log Tc  log ROc   Mrad     Rrad       k2conv      k2rad
#0.010     5.695210     2388.  -2.469  3.372  0.341   0.0000   5.417   0.2639  0.0000   0.000E+00  4.535E-01  0.000E+00
#0.015     5.693523     2514.  -2.208  3.376  0.416   0.0000   5.503   0.2014  0.0000   0.000E+00  4.509E-01  0.000E+00
#0.020     5.689884     2594.  -2.044  3.392  0.472   0.0000   5.566   0.1700  0.0000   0.000E+00  4.502E-01  0.000E+00
#0.030     5.695099     2666.  -1.783  3.355  0.603   0.0000   5.634   0.0491  0.0000   0.000E+00  4.459E-01  0.000E+00
#0.040     5.694350     2731.  -1.655  3.394  0.665   0.0000   5.707   0.0471  0.0000   0.000E+00  4.464E-01  0.000E+00
#0.050     5.694123     2761.  -1.481  3.335  0.796   0.0000   5.729  -0.0741  0.0000   0.000E+00  4.424E-01  0.000E+00
#0.060     5.694808     2807.  -1.399  3.361  0.846   0.0000   5.778  -0.0728  0.0000   0.000E+00  4.429E-01  0.000E+00
#0.070     5.692913     2834.  -1.324  3.369  0.905   0.0000   5.813  -0.0918  0.0000   0.000E+00  4.429E-01  0.000E+00
#0.072     5.694756     2831.  -1.291  3.347  0.942   0.0000   5.811  -0.1271  0.0000   0.000E+00  4.420E-01  0.000E+00
#0.075     5.694352     2835.  -1.261  3.338  0.972   0.0000   5.816  -0.1479  0.0000   0.000E+00  4.416E-01  0.000E+00
#0.080     5.694252     2836.  -1.197  3.302  1.045   0.0000   5.814  -0.2109  0.0000   0.000E+00  4.401E-01  0.000E+00
#0.090     5.689450     2861.  -1.127  3.299  1.113   0.0000   5.837  -0.2391  0.0000   0.000E+00  4.400E-01  0.000E+00
#0.100     5.703957     2925.  -1.154  3.410  1.033   0.0000   5.904  -0.1091  0.0000   0.000E+00  4.440E-01  0.000E+00
#0.110     5.709729     2945.  -1.075  3.384  1.115   0.0000   5.914  -0.1628  0.0000   0.000E+00  4.430E-01  0.000E+00
#0.130     5.702708     3006.  -0.926  3.344  1.270   0.0000   5.931  -0.2541  0.0000   0.000E+00  4.416E-01  0.000E+00
#0.150     5.701659     3076.  -0.795  3.314  1.412   0.0000   5.947  -0.3249  0.0000   0.000E+00  4.407E-01  0.000E+00
#0.170     5.693846     3139.  -0.669  3.278  1.568   0.0000   5.955  -0.4031  0.0000   0.000E+00  4.399E-01  0.000E+00
#0.200     5.693823     3220.  -0.539  3.262  1.731   0.0000   5.980  -0.4614  0.0000   0.000E+00  4.400E-01  0.000E+00
#0.300     5.690044     3460.  -0.199  3.224  2.215   0.0000   6.042  -0.6068  0.0000   0.000E+00  4.406E-01  0.000E+00
#0.400     5.689995     3672.  -0.040  3.292  2.364   0.0000   6.133  -0.5733  0.0000   0.000E+00  4.432E-01  0.000E+00
#0.500     5.694573     3849.   0.076  3.355  2.458   0.0000   6.209  -0.5319  0.0000   0.000E+00  4.450E-01  0.000E+00
#0.600     5.690793     3988.   0.171  3.402  2.552   0.0000   6.269  -0.5046  0.0000   0.000E+00  4.462E-01  0.000E+00
#0.700     5.691618     4111.   0.268  3.424  2.687   0.0000   6.312  -0.5062  0.0000   0.000E+00  4.469E-01  0.000E+00
#0.800     5.693354     4221.   0.356  3.440  2.821   0.0000   6.348  -0.5124  0.0000   0.000E+00  4.473E-01  0.000E+00
#0.900     5.692710     4315.   0.436  3.449  2.960   0.0000   6.377  -0.5246  0.0000   0.000E+00  4.477E-01  0.000E+00
#1.000     5.693063     4397.   0.508  3.456  3.096   0.0000   6.402  -0.5378  0.0000   0.000E+00  4.479E-01  0.000E+00
#1.100     5.694315     4471.   0.573  3.461  3.227   0.0000   6.425  -0.5508  0.0000   0.000E+00  4.481E-01  0.000E+00
#1.200     5.693263     4537.   0.634  3.463  3.362   0.0000   6.444  -0.5666  0.0000   0.000E+00  4.483E-01  0.000E+00
#1.300     5.694626     4596.   0.688  3.466  3.488   0.0000   6.462  -0.5799  0.0000   0.000E+00  4.485E-01  0.000E+00
#1.400     5.692977     4647.   0.740  3.466  3.621   0.0000   6.478  -0.5965  0.0000   0.000E+00  4.486E-01  0.000E+00

Baraffe_mass = yt.YTArray([0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.072, 0.075, 0.080, 0.090, 0.100, 0.110, 0.130, 0.150, 0.170, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400], 'msun')
Baraffe_logL = np.array([-2.469, -2.208, -2.044, -1.783, -1.655, -1.481, -1.399, -1.324, -1.291, -1.261, -1.197, -1.127, -1.154, -1.075, -0.926, -0.795, -0.669, -0.539, -0.199, -0.040, 0.076, 0.171, 0.268, 0.356, 0.436, 0.508, 0.573, 0.634, 0.688, 0.740])
Baraffe_radius = yt.YTArray([0.341, 0.416, 0.472, 0.603, 0.665, 0.796, 0.846, 0.905, 0.942, 0.972, 1.045, 1.113, 1.033, 1.115, 1.270, 1.412, 1.568, 1.731, 2.215, 2.364, 2.458, 2.552, 2.687, 2.821, 2.960, 3.096, 3.227, 3.362, 3.488, 3.621], 'rsun')

#Derive a stellar luminosity
lstar_baraffe = []
rstar_barrafe = []
for mass_val in mass:
    if mass_val < Baraffe_mass[0]:
        lstar_baraffe.append(10**Baraffe_logL[0])
        rstar_barrafe.append(Baraffe_radius[0])
    else:
        closest_inds = sorted(np.argsort(np.abs(Baraffe_mass - mass_val))[:2])
        gradient = (Baraffe_logL[closest_inds][1] - Baraffe_logL[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
        y_intercept = Baraffe_logL[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
        logL = gradient*mass_val + y_intercept
        lstar_baraffe.append(10**logL)
        
        gradient = (Baraffe_radius[closest_inds][1] - Baraffe_radius[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
        y_intercept = Baraffe_radius[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
        radius = gradient*mass_val + y_intercept
        rstar_barrafe.append(radius)

lstar_baraffe = yt.YTArray(lstar_baraffe, 'Lsun')
lacc = facc * (yt.units.gravitational_constant_cgs * mass.in_units('g') * mdot.in_units('g/s'))/yt.YTArray(rstar_barrafe).in_units('cm')
ltot = lacc.in_units('lsun') + lstar_baraffe
magnitude = -2.5 * np.log10(ltot.in_units('watt')/3.0128e28)
plt.clf()
plt.semilogy(age, ltot, label='Total', alpha=0.5)
plt.semilogy(age, lacc.in_units('lsun'), label='Acc', alpha=0.5)
plt.semilogy(age, lstar_baraffe, label='Star', alpha=0.5)
plt.legend()
plt.savefig('L_component_evol.png')

plt.clf()
plt.semilogy(age, ltot, label='Total')
plt.legend()
plt.savefig('L_tot_evol.png')

plt.clf()
plt.plot(age, magnitude)
plt.gca().invert_yaxis()
plt.xlim([9000, 11500])
plt.ylim([8,3])
plt.xlabel('Time (yr)')
plt.ylabel('Magnitude')
plt.savefig('Magnitude_evol.png')

plt.clf()
plt.plot(age, ltot)
plt.xlim([9000, 11500])
plt.ylim([0, 4])
plt.xlabel('Time (yr)')
plt.ylabel('Ltot')
plt.savefig('Ltot_evol_zoom.png')

plot_times = [10317.928611457348, 10861.812506761402, 12135.403063911945, 13096.016646496952, 13379.528082296252, 13908.540377408266, 14625.588010121137, 15691.680855810642, 15891.27923379466, 17353.72329491377, 26899.01125465706, 29443.443914979696] #14154.07486982271

plt.clf()
fig, axs = plt.subplots(ncols=3, nrows=4, figsize=(two_col_width, 1.4*two_col_width), sharey=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.1)

rank = CW.Get_rank()
size = CW.Get_size()
L_diff_arr = []
M_diff_arr = []
time_arr = []
cor_arr = []
L_lims = [-1.8, 0.8]
M_lims = [9.1,2.7]
for time_it in range(len(age)):
    end_time = age[time_it] + time_window
    end_it = np.argmin(abs(age - end_time))
    useable_times = age[time_it:end_it]
    useable_M = magnitude[time_it:end_it]
    useable_L = ltot[time_it:end_it]
    if len(useable_L) > 0:
        L_diff = np.max(np.log10(useable_L)) - np.min(np.log10(useable_L))
        L_diff_arr.append(L_diff)
        M_diff = np.max(useable_M) - np.min(useable_M)
        M_diff_arr.append(M_diff)
        time_arr.append(age[time_it])
        scaled_T = useable_times - useable_times[0]
        scaled_L = useable_L - np.min(useable_L)
        scaled_L = scaled_L/np.max(scaled_L)
        scaled_M = useable_M - np.min(useable_M)
        scaled_M = scaled_M/np.max(scaled_M)
        cor_M = np.correlate(scaled_M,FU_temp_inv,'same')
        cor_L = np.correlate(scaled_L,FU_temp,'same')
        median_cor_L = np.median(cor_L[np.where(cor_L>0)[0]])
        median_cor_M = np.median(cor_M[np.where(cor_M>0)[0]])
        cor_arr.append([median_cor_L, median_cor_M])
        if age[time_it] in plot_times:
            plot_it = plot_times.index(age[time_it])
            axs.flatten()[plot_it].plot(useable_times.in_units('kyr'), np.log10(useable_L), label="Total Luminosity")
            axs.flatten()[plot_it].set_xlim([useable_times.in_units('kyr')[0], useable_times.in_units('kyr')[-1]])
            right_ax = axs.flatten()[plot_it].twinx()
            right_ax.invert_yaxis()
            right_ax.plot(useable_times.in_units('kyr'), useable_M, c='orange', ls='--', label="Magnitude")
            if plot_it == 0:
                axs.flatten()[plot_it].legend(loc='lower right')
            if plot_it == 1:
                right_ax.legend(loc='lower right')
            if np.remainder(plot_it, 3) == 0:
                axs.flatten()[plot_it].set_ylabel("Log L (L$_\odot$)")
            if np.remainder(plot_it, 3) == 2:
                right_ax.set_ylabel("Absolute Magnitude")
            else:
                yticklabels = right_ax.get_yticklabels()
                plt.setp(yticklabels, visible=False)
            if plot_it == 9:
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-2], visible=False)
            if plot_it > 8:
                axs.flatten()[plot_it].set_xlabel("Time (kyr)")
            axs.flatten()[plot_it].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            '''
            if np.min(np.log10(useable_L)) < L_lims[0]:
                L_lims[0] = np.min(np.log10(useable_L))
            if np.max(np.log10(useable_L)) > L_lims[1]:
                L_lims[1] = np.max(np.log10(useable_L))
            '''
            '''
            if np.max(useable_M) > M_lims[0]:
                M_lims[0] = np.max(useable_M)
            if np.min(useable_M) < M_lims[1]:
                M_lims[1] = np.min(useable_M)
            '''
            axs.flatten()[plot_it].set_ylim(L_lims)
            right_ax.set_ylim(M_lims)
            axs.flatten()[plot_it].tick_params(axis='both', direction='in', top=True)
            right_ax.tick_params(axis='both', direction='in')
            plt.savefig('multiplot_burst_events.pdf', bbox_inches='tight', pad_inches=0.02)
            print("updated multiplot with time", age[time_it])
            if age[time_it] == plot_times[-1]:
                break
            
            
        '''
        if L_diff>2 and M_diff>5 and median_cor_L>29 and np.min(np.log10(useable_L)) not in np.log10(useable_L)[-10:]:
            #import pdb
            #pdb.set_trace()
            plt.clf()
            fig, left_ax = plt.subplots(ncols=1, nrows=1)
            right_ax = left_ax.twinx()
            left_ax.plot(useable_times, np.log10(useable_L), c='b')
            left_ax.set_ylabel("Log L (L$_\odot$)")
            left_ax.set_xlabel("Time (yr)")
            #left_ax.plot(scaled_T, scaled_L, c='b')
            #left_ax.plot(scaled_T, scaled_M, c='b', ls="--")
            #left_ax.plot(np.linspace(0, scaled_T[-1], len(FU_temp)), FU_temp, c='orange')
            #left_ax.set_ylim([0,1.05])
            right_ax.plot(np.linspace(useable_times[0], useable_times[-1], len(cor_L)), cor_L, c='g')
            #right_ax.plot(np.linspace(0, scaled_T[-1], len(cor_M)), cor_M, c='g', ls="--")
            right_ax.set_ylim(bottom=0)
            right_ax.set_ylabel("Correlation")
            plt.title('Max cor ='+ str(np.nanmax(cor_L)))
            plt.show()
            plt.savefig("Scaled_T_"+str(age[time_it])+".png")
        '''
        '''
        if median_cor_L>30 and median_cor_M>30 and L_diff > 5: #and mass[time_it] > 0.1
            plt.clf()
            fig, ax1 = plt.subplots()

            ax2 = ax1.twinx()
            ax1.invert_yaxis()
            ax1.plot(useable_times, scaled_L, label="scaled L_acc", color='b')
            
            ax2.plot(useable_times[0]+np.linspace(0, scaled_T[-1], len(cor)), cor, label="correlation", color='r')

            ax1.set_xlabel('Time (yr)')
            ax1.set_ylabel('scaled L_acc and correlation')
            ax2.set_ylabel('Total log Luminosity')
            
            ax1.set_ylim([1, 0])
            ax2.set_ylim([np.min(useable_L), np.max(useable_L)])
        
            ax1.legend()
            plt.savefig('Best_match_time_'+str(age[time_it])+'_mass_'+str(np.round(mass[time_it], decimals=2))+'.png',  bbox_inches='tight')
            print("Found potential match at age", age[time_it])
        '''
plt.clf()
plt.plot(time_arr, L_diff_arr)
plt.xlabel('age (yr)')
plt.ylabel('max L diff over 80 yr (log)')
plt.savefig('L_diff.png')
print("plotted L diff history for sink 45 on rank", rank)
