#!/usr/bin/env python
#Plots MESA luminsoties calculated for all sinks int eh G100 512 run
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle

lsun = 3.828e26*1.e7 # solar luminosity in erg
FU_temp = np.concatenate((np.zeros(20), np.ones(80)))
time_window = 80

rank = CW.Get_rank()
size = CW.Get_size()

try:
    sink_files = sorted(glob.glob('/g/data/ek9/rlk100/RAMSES/Global/G100/512_Resolution/Mesa_pickles/Mesa_pickle_0*.pkl'))
    use_pickles = True
    #sink_files = sorted(glob.glob('/data/scratch/troels/IMF_512/mesa/sink_*/LOGS'))
except:
    sink_files = sorted(glob.glob('/home/100/rlk100/rlk/RAMSES/Analysis/MESA_raw/sink_*/LOGS'))
    import mesaPlot
    m=mesaPlot.MESA()
    p=mesaPlot.plot()
    use_pickles = False
#sink_files = sorted(glob.glob('/lustre/astro/troels/IMF_512/mesa/success/sink_*/LOGS'))
#if rank == 0:
#    print('sink_files:', sink_files)
#max_age=150000
rit = -1
sink_it = -1
best_corr_pickle = 'best_cors_'+str(rank)+'.pkl'
if os.path.isfile(best_corr_pickle):
    cors_file = open('best_cors_'+str(rank)+'.pkl', 'rb')
    best_sink, best_time, best_corr = pickle.load(cors_file)
    cors_file.close()
else:
    best_sink = np.array([])
    best_time = np.array([])
    best_corr = np.array([])
for sink_file in sink_files:
    sink_it = sink_it + 1
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        if os.path.isfile('L_diff_Sink_'+str(sink_it)+'.png') == False:
            if use_pickles == False:
                try:
                    m.log_fold=sink_file
                    m.loadHistory()
                    could_read_file = True
                    print('Successfully read', sink_file)
                except:
                    print("Can't read sink_file")
                    could_read_file = False
            else:
                pickle_open = open(sink_file, "rb")
                pickle_data = pickle.load(pickle_open)
                pickle_open.close()
                could_read_file = True
                
            if could_read_file:
                if use_pickles == False:
                    mass = m.hist.star_mass
                    age = m.hist.star_age
                    #idx = np.where(age <= max_age)
                    #age = age[idx]
                    lum = 10.**m.hist.log_L
                    lacc = m.hist.extra_lum / lsun
                    ltot = lum + lacc
                    
                    #save Mesa sink data to pickle
                    pickle_data = {"age":age, "mass":mass, "stellar_luminosity":lum, "accretion_luminosity":lacc, "total_luminosity":ltot}
                    pickle_open = open("Mesa_pickle_"+sink_file.split('sink_')[-1].split('/')[0]+"_full_age.pkl", "wb")
                    pickle.dump((pickle_data), pickle_open)
                    pickle_open.close()
                    print("saved pickle data for sink", sink_file.split('sink_')[-1].split('/')[0])
                else:
                    age = pickle_data['age']
                    mass = pickle_data['mass']
                    lum = pickle_data['stellar_luminosity']
                    lacc = pickle_data['accretion_luminosity']
                    ltot = pickle_data['total_luminosity']
                time_arr = []
                L_diff_arr = []
                
                #have moving window:
                for time_it in range(len(age)):
                    end_time = age[time_it] + time_window
                    end_it = np.argmin(abs(age - end_time))
                    useable_times = age[time_it:end_it]
                    useable_L = ltot[time_it:end_it]
                    if len(useable_L) > 0:
                        L_diff = np.log10(useable_L[-1]) - np.log10(useable_L[0])
                        L_diff_arr.append(L_diff)
                        time_arr.append(age[time_it])
                        useable_L = np.log10(useable_L)
                        #L_diff = np.max(np.log10(useable_L)) - np.min(np.log10(useable_L))
                        scaled_T = useable_times - useable_times[0]
                        scaled_L = useable_L - np.min(useable_L)
                        scaled_L = scaled_L/np.max(scaled_L)
                        cor = np.correlate(scaled_L,FU_temp,'same')
                        if L_diff>1 and np.median(cor)>66.6:# and mass[time_it] > 0.1:
                            plt.clf()
                            fig, ax1 = plt.subplots()

                            ax2 = ax1.twinx()
                            ax1.plot(useable_times, scaled_L, label="scaled Luminosity", color='b')
                            ax1.plot(useable_times, cor[:len(useable_times)]/100., label="correlation", color='r')
                            
                            ax2.plot(useable_times, useable_L, color='b')

                            ax1.set_xlabel('Time (yr)')
                            ax1.set_ylabel('scaled L and correlation')
                            ax2.set_ylabel('Total log Luminosity')
                            
                            ax1.set_ylim([0, 1])
                            ax2.set_ylim([np.min(useable_L), np.max(useable_L)])
                        
                            ax1.legend()
                            plt.savefig('Sink_' + str(sink_it) + '_time_'+str(age[time_it])+'_mass_'+str(np.round(mass[time_it], decimals=2))+'.png',  bbox_inches='tight')
                            print("Found potential match for sink", sink_it, "at age", age[time_it])
                            best_sink = np.append(best_sink, sink_it)
                            best_time = np.append(best_time, age[time_it])
                            best_corr = np.append(best_corr, np.median(cor))
                plt.clf()
                plt.plot(time_arr, L_diff_arr)
                plt.xlabel('age (yr)')
                plt.ylabel('max L diff over 100 yr (log)')
                plt.savefig('L_diff_Sink_'+str(sink_it)+'.png')
                print("plotted L diff history for sink", sink_it, "on rank", rank)
                cors_file = open('best_cors_'+str(rank)+'.pkl', 'wb')
                pickle.dump((best_sink, best_time, best_corr), cors_file)
                cors_file.close()
                
if rank == 0:
    best_sink_all = np.array([])
    best_time_all = np.array([])
    best_corr_all = np.array([])
    
    best_corr_pickle = sorted(glob.glob('best_cors_*.pkl'))
    for best_pick in best_corr_pickle:
        cors_file = open(best_pick, 'rb')
        best_sink, best_time, best_corr = pickle.load(cors_file)
        cors_file.close()
        best_sink_all = np.append(best_sink_all, best_sink)
        best_time_all = np.append(best_time_all, best_time)
        best_corr_all = np.append(best_corr_all, best_corr)
        
    
    top_sinks = []
    top_times = []
    top_corrs = []
    
    for ind in np.argsort(best_corr_all)[::-1]:
        if best_sink_all[ind] not in top_sinks:
            top_sinks.append(best_sink_all[ind])
            top_times.append(best_time_all[ind])
            top_corrs.append(best_corr_all[ind])


    top_times = np.array(top_times)[np.argsort(top_sinks)]
    top_sinks = np.array(top_sinks)[np.argsort(top_sinks)]
    top_clean = np.array([177, 292, 48, 51, 262, 195, 17, 10, 75, 159, 272, 275, 176, 118, 54, 45, 85, 103, 71, 101, 258, 150, 93, 221, 151, 154, 102, 168, 175, 56, 309, 239, 109, 73, 72, 83, 141])
    top_clean = top_clean[np.argsort(top_clean)]
    
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472 #inches
    font_size = 7
    
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['mathtext.fontset'] = 'custom' #stixsans'
    matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
    matplotlib.rcParams['mathtext.rm'] = 'Arial'
    matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
    matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
    matplotlib.rcParams['mathtext.rm'] = 'Arial'
    matplotlib.rcParams['mathtext.sf'] = 'Arial'
    matplotlib.rcParams['mathtext.default'] = 'regular'
    
    plt.cla()
    plt.clf()
    fig, axs = plt.subplots(ncols=5, nrows=2, figsize=(two_col_width, 0.48*two_col_width), sharey=True, linewidth=1)
    plt.subplots_adjust(wspace=0.0)
    plt.subplots_adjust(hspace=0.13)
    
    for plot_it in range(len(top_clean[:10])):
        pickle_open = open('Mesa_pickle_'+("%04d" % top_clean[plot_it])+'_full_age.pkl', "rb")
        pickle_data = pickle.load(pickle_open)
        pickle_open.close()
        
        age = pickle_data['age']
        mass = pickle_data['mass']
        lum = pickle_data['stellar_luminosity']
        lacc = pickle_data['accretion_luminosity']
        ltot = pickle_data['total_luminosity']
        
        top_sink_it = np.where(top_sinks==top_clean[plot_it])[0]
        plot_time = top_times[top_sink_it]
        time_it = np.argmin(abs(age - plot_time))
        end_time = age[time_it] + time_window
        end_it = np.argmin(abs(age - end_time))
        useable_times = age[time_it:end_it]
        useable_L = ltot[time_it:end_it]
        useable_L = np.log10(useable_L)
        scaled_T = useable_times - useable_times[0]
        scaled_L = useable_L - np.min(useable_L)
        scaled_L = scaled_L/np.max(scaled_L)
        cor = np.correlate(scaled_L,FU_temp,'same')
    
        ax1 = axs.flatten()[plot_it]
        #plt.gca().set_aspect('equal')
        ax2 = ax1.twinx()
        ax1.ticklabel_format(useOffset=False)
        ax1.plot(useable_times/1000, scaled_L, label="Scaled Lum.", color='b')
        ax1.plot(useable_times/1000, cor[:len(useable_times)]/100., label="Correlation", color='r')
                            
        ax2.plot(useable_times/1000, useable_L, color='b')

        if plot_it >4:
            ax1.set_xlabel('Time (kyr)', fontsize=font_size, labelpad=-1)
        if np.remainder(plot_it, 5) == 0:
            ax1.set_ylabel('scaled L and correlation', fontsize=font_size, labelpad=0)
        else:
            yticklabels = ax1.get_yticklabels()
            plt.setp(yticklabels, visible=False)
        if plot_it == 4 or plot_it == 9:
            ax2.set_ylabel('Total log Luminosity', fontsize=font_size, labelpad=0)
        else:
            yticklabels = ax2.get_yticklabels()
            plt.setp(yticklabels, visible=False)
            
        if plot_it == 5 or plot_it == 8:
            xticklabels = ax1.get_xticklabels()
            plt.setp(xticklabels[-1], visible=False)
                            
        #ax1.set_xlim([np.min(useable_times), np.max(useable_times)])
        ax1.set_ylim([0, 1])
        ax2.set_ylim([np.min(useable_L), np.max(useable_L)])
        
        ax1.tick_params(axis='x', which='major', direction='in', color='k', top=True, length=2)
        ax1.tick_params(axis='y', which='major', direction='in', color='k', length=0)
        ax2.tick_params(axis='y', which='major', direction='in', color='k', length=3)
        ax1.xaxis.label.set_color('black')
        ax1.yaxis.label.set_color('black')
        ax1.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
        ax2.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
        
        if plot_it == 0:
            ax1.legend(loc="center", fontsize=font_size)
            
        useable_times = useable_times/1000
        Cand_string = "Cand. "+str(plot_it+1)
        Cand_string_raw = r"{}".format(Cand_string)
        Cand_text = ax1.text(np.max(useable_times), 0.15, Cand_string_raw, va="center", ha="right", color='k', fontsize=font_size)
        
        Corr_string = "Med. Corr="+str(np.round(np.median(cor), decimals=2))
        Corr_string_raw = r"{}".format(Corr_string)
        Corr_text = ax1.text(np.max(useable_times), 0.05, Corr_string_raw, va="center", ha="right", color='k', fontsize=font_size)
        
        plt.savefig('Main_body_best_matches.pdf', bbox_inches='tight', pad_inches=0.02)
        print('Updated Main_body_best_matches.pdf with sink', top_clean[plot_it])





