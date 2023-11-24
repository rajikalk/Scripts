import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt

pickle_files = ['Lref_8.pkl', 'Lref_9.pkl', 'Lref_10.pkl']#, 'Lref_11.pkl']
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
scaling_factor = [1, 2, 4, 8]

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pickle_file in pickle_files:
    plot_it = plot_it + 1
    file = open(pickle_file, 'rb')
    sink_data, line_counter = pickle.load(file)
    file.close()
    
    form_time = np.inf
    sink_it = -1
    primary_mass = 0
    primary_ind = np.nan
    for sink_id in sink_data.keys():
        mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
        if mass[-1].in_units('msun').value > primary_mass:
            primary_mass = mass[-1].in_units('msun').value
            primary_ind = sink_id
    
    form_time = sink_data[primary_ind]['time'][0]
    
    for sink_id in sink_data.keys():
        sink_it = sink_it + 1

        time = sink_data[sink_id]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        L_spec_tot = np.sqrt((sink_data[sink_id]['anglx']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['angly']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['anglz']/sink_data[sink_id]['mass'])**2)
        L_spec_tot = yt.YTArray(L_spec_tot, 'cm**2/s')
        plot_L_spec_tot = L_spec_tot.in_units('m**2/s')[:end_ind+1]

        L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        plot_L_tot = L_tot.in_units('kg*m**2/s')[:end_ind+1]
        
        mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
        plot_mass = mass.in_units('msun')[:end_ind+1]
    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].plot(plot_time, plot_mass*scaling_factor[plot_it], linestyle=line_styles[plot_it], label=pickle_file.split('.')[0], color=colors[plot_it])
        else:
            axs.flatten()[0].plot(plot_time, plot_mass*scaling_factor[plot_it], linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[1].plot(plot_time, plot_L_tot*scaling_factor[plot_it], linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[2].plot(plot_time, plot_L_spec_tot*scaling_factor[plot_it], linestyle=line_styles[plot_it], color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)")
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylabel("L (kg$\,$m$^2$/s)")
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylabel("h (m$^2$/s)")
axs.flatten()[2].set_xlabel("time (yr)")
axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[0].legend(loc='best')

plt.savefig('resolution_study.pdf', bbox_inches='tight')
        
    
