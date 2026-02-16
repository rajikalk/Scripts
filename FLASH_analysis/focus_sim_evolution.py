import numpy as np
import pickle5 as pickle
import matplotlib.pyplot as plt
import yt

#pickle_files = ['sink.pkl', 'disk.pkl']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
line_styles = [':', '-.', '--', '-']# ['-', '--', '-.', ':']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
plot_it = -1
end_time = 10000

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=3, figsize=(two_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

file = open('sink.pkl', 'rb')
sink_data, line_counter = pickle.load(file)
file.close()

form_time = np.inf
sink_it = -1
for sink_id in sink_data.keys():
    sink_it = sink_it + 1
    if sink_data[sink_id]['time'][0] < form_time:
        form_time = sink_data[sink_id]['time'][0]

    time = sink_data[sink_id]['time'] - form_time
    time = yt.YTArray(time, 's')
    end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
    plot_time = time.in_units('yr')[:end_ind+1]
    
    L_spec_tot = np.sqrt((sink_data[sink_id]['anglx']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['angly']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['anglz']/sink_data[sink_id]['mass'])**2)
    L_spec_tot = yt.YTArray(L_spec_tot, 'cm**2/s')
    plot_L_spec_tot = L_spec_tot.in_units('cm**2/s')[:end_ind+1]

    L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
    L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
    plot_L_tot = L_tot.in_units('g*cm**2/s')[:end_ind+1]
    
    mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
    plot_mass = mass.in_units('msun')[:end_ind+1]

    axs.flatten()[0].plot(plot_time, plot_mass)
    axs.flatten()[2].plot(plot_time, plot_L_tot)
    axs.flatten()[4].plot(plot_time, plot_L_spec_tot)
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)")
axs.flatten()[0].set_title("Star")
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[2].set_ylabel("L (g$\,$cm$^2$/s)")
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[4].set_ylabel("h (cm$^2$/s)")
axs.flatten()[4].set_xlabel("time (yr)")
axs.flatten()[4].set_xlim([0, 10000])
axs.flatten()[4].set_ylim(bottom=0)
#axs.flatten()[0].legend(loc='best')

#==================================
#Plot disk

file = open('disk.pkl', 'rb')
Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mass_all, Separation = pickle.load(file)
file.close()

axs.flatten()[1].plot(Time_array, Mass_all)
axs.flatten()[3].plot(Time_array, Total_L)
axs.flatten()[5].plot(Time_array, Total_L_spec)
axs.flatten()[5].set_xlabel("time (yr)")
axs.flatten()[1].set_title("Disk")

plt.savefig('inner_sys_evolution.pdf', bbox_inches='tight')
        
    
