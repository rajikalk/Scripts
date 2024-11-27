import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
sink_evol_pickle = sys.argv[1]

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

for key in sink_data.keys():
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=False)
    time = sink_data[key]['time'] - sink_data[key]['time'][0]
    time = yt.YTArray(time, 's').in_units('yr')
    mass = yt.YTArray(sink_data[key]['mass'], 'g').in_units('Msun')
    mdot = yt.YTArray(sink_data[key]['mdot'], 'g/s').in_units('Msun/yr')
    time_mid = (time[:-1] + time[1:])/2
    m_dot_calc = (mass[1:]-mass[:-1])/(time[1:]-time[:-1])
    m_dot_calc[np.argwhere(m_dot_calc<0)] = np.nan
    
    axs[0].plot(time, mass)
    axs[0].set_title('Sink id:'+key)
    axs[0].set_ylim(bottom=0)
    axs[0].set_xlim(left=0)
    axs[0].set_ylabel('Mass (M$_\odot$)')
    axs[1].plot(time_mid, m_dot_calc)
    axs[1].set_ylim(bottom=0)
    axs[1].set_ylabel('Mass accretion (M$_\odot$/yr)')
    axs[1].set_xlabel('Time (yr)')
    
    plt.savefig(key+'.png')
    print('saved '+key+'.png')

