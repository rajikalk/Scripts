#Plots MESA luminsoties calculated for all sinks int eh G100 512 run

import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle

lsun = 3.828e26*1e7 # solar luminosity in erg
FU_temp = np.concatenate((np.zeros(20), np.ones(80)))
time_window = 80

rank = CW.Get_rank()
size = CW.Get_size()

try:
    sink_files = sorted(glob.glob('/data/scratch/troels/IMF_512/mesa/sink_*/LOGS'))
    import mesaPlot
    m=mesaPlot.MESA()
    p=mesaPlot.plot()
    use_pickles = False
except:
    sink_files = sorted(glob.glob('/g/data/ek9/rlk100/RAMSES/Global/G100/512_Resolution/Mesa_pickles/Mesa_pickle_0*.pkl'))
    use_pickles = True
#sink_files = sorted(glob.glob('/lustre/astro/troels/IMF_512/mesa/success/sink_*/LOGS'))
#if rank == 0:
#    print('sink_files:', sink_files)
max_age=150000
rit = -1
sink_it = -1
best_cors = []
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
                    idx = np.where(age <= max_age)
                    age = age[idx]
                    lum = 10.**m.hist.log_L
                    lacc = m.hist.extra_lum / lsun
                    ltot = lum + lacc
                    
                    #save Mesa sink data to pickle
                    pickle_data = {"age":age, "mass":mass, "stellar_luminosity":lum, "accretion_luminosity":lacc, "total_luminosity":ltot}
                    pickle_open = open("Mesa_pickle_"+sink_file.split('sink_')[-1].split('/')[0]+".pkl", "wb")
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
                            best_cors.append([sink_it, age[time_it], np.median(cor)])
                plt.clf()
                plt.plot(time_arr, L_diff_arr)
                plt.xlabel('age (yr)')
                plt.ylabel('max L diff over 100 yr (log)')
                plt.savefig('L_diff_Sink_'+str(sink_it)+'.png')
                print("plotted L diff history for sink", sink_it, "on rank", rank)
                cors_file = open('best_cors.pkl', 'wb')
                pickle.dump((np.array(best_cors)), cors_file)
                cors_file.close()
                

##Data for multiply plot
plot_sink = [17, 17, 17, 17, 17, 51, 74, 260, 270, 273, 290, 290, 290, 290, 290, 174, 175, 54, 84, 45]
plot_time = [69614.81013489112, 77207.60066041496, 77340.74994086033, 83816.50575862321, 88474.8956793223, 34617.36712460182, 30982.481234925883, 36357.680496830486, 45271.959678627885, 23027.750556024763, 20559.836032810344, 26877.25760296158, 26975.78215359961, 37188.631669261405, 50732.408920391805, 73943.3525513, 149484.80417919, 112569.82249152, 47797.44646577, 34786.14689247]
plot_corr = [73.34874538, 78.29124944, 77.15026779, 77.74458688, 78.0603211, 78.96473908, 78.22833953, 78.59386108, 77.8718668, 77.30801646, 78.3453675, 79.30769461, 78.78456987, 78.13137853, 78.41088411, 78.03140985, 76.67806032, 76.23207367, 76.1266242, 76.07275835]


