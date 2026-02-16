import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
sink_evol_pickle = sys.argv[1]
primary_key = sys.argv[2]
secondary_key = sys.argv[3]

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

Secondary_form_time = sink_data[secondary_key]['time'][0]
form_ind = int(np.where((sink_data[primary_key]['time']-Secondary_form_time)==0)[0])
plot_time = sink_data[primary_key]['time'][form_ind:]
plot_time = plot_time - plot_time[0]
plot_time = yt.YTArray(plot_time, 's').in_units('yr')
 econdary_key]['posx']
dy = sink_data[primary_key]['posy'][form_ind:] - sink_data[secondary_key]['posy']
dz = sink_data[primary_key]['posz'][form_ind:] - sink_data[secondary_key]['posz']
sep = np.sqrt(dx**2 + dy**2 + dz**2)
sep = yt.YTArray(sep, 'cm').in_units('AU')

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=False)
axs[0].plot(plot_time, yt.YTArray(sink_data[primary_key]['mass'][form_ind:], 'g').in_units('msun'), label='Primary Star')
axs[0].plot(plot_time, yt.YTArray(sink_data[secondary_key]['mass'], 'g').in_units('msun'), label='Secondary Star')
axs[0].set_ylabel('Mass (M$_\odot$)')
axs[1].semilogy(plot_time, sep)AVZ
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Separation (au)')
plt.savefig('Separation_'+ primary_key + '_' + secondary_key + '.png')

