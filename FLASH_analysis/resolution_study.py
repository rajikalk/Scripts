import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
import sys

r_sink_9 = yt.YTQuantity(4.89593797, 'AU')
L_reference = int(sys.argv[1])
r_scale = (1/(2**(L_reference-9)))*r_sink_9
pickle_files = ['Lref_8.pkl', 'Lref_9.pkl', 'Lref_10.pkl', 'Lref_11.pkl']
R_proto = yt.YTQuantity(2, 'rsun')
r_sink = [2*r_sink_9, r_sink_9, (1/(2))*r_sink_9, (1/(2**2))*r_sink_9, R_proto.in_units('au')]
#          8            9           10                  11 need Lref 18 to resolve protostellar surface
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
L_scale= np.sqrt(2)/16
L_spec_scale = np.sqrt(2)/2
m_scale = 1/8
scaling_factor = [L_reference-8, L_reference-9, L_reference-10, L_reference-11]
comp_time = np.nan
comp_h = []
comp_period = []

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
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
    if np.isnan(comp_time):
        comp_time = yt.YTQuantity((sink_data[primary_ind]['time'] - form_time)[-1], 's')
    
    for sink_id in sink_data.keys():
        sink_it = sink_it + 1

        time = sink_data[sink_id]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        L_spec_tot = np.sqrt((sink_data[sink_id]['anglx']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['angly']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['anglz']/sink_data[sink_id]['mass'])**2)
        L_spec_tot = yt.YTArray(L_spec_tot, 'cm**2/s')
        plot_L_spec_tot = L_spec_tot.in_units('m**2/s')[:end_ind+1]
        
        #Sink_rot = plot_L_spec_tot.in_units('km**2/s') / r_sink[plot_it].in_units('km') #I think this is wrong
        #Sink_circum = 2*np.pi * r_sink[plot_it].in_units('km')
        #Sink_period = Sink_circum/Sink_rot
        
        L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        plot_L_tot = L_tot.in_units('kg*m**2/s')[:end_ind+1]
        
        mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
        plot_mass = mass.in_units('msun')[:end_ind+1]
        
        #===================================================
        #Calculate rotation Period
        #v_scale = plot_L_spec_tot.in_units('km**2/day')*(L_spec_scale**(scaling_factor[plot_it]))/r_sink_9.in_units('km')
        v_scale = plot_L_spec_tot.in_units('km**2/day')*(L_spec_scale**(scaling_factor[plot_it]))/r_scale.in_units('km')
        Sink_circum = 2*np.pi * r_scale.in_units('km')
        Sink_period = Sink_circum/v_scale
        
        #Sink_period = 0.8 * np.pi * plot_mass.in_units('g') * r_sink[plot_it].in_units('cm')**2 / plot_L_tot.in_units('g*cm**2/day')
        
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].plot(plot_time, plot_mass*(m_scale**(scaling_factor[plot_it])), linestyle=line_styles[plot_it], label="R$_{sink}$="+str(np.round(r_sink[plot_it], decimals=2))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].plot(plot_time, plot_mass*(m_scale**(scaling_factor[plot_it])), linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[1].plot(plot_time, plot_L_tot*(L_scale**(scaling_factor[plot_it])), linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[2].plot(plot_time, plot_L_spec_tot*(L_spec_scale**(scaling_factor[plot_it])), linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)")
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylabel("L (kg$\,$m$^2$/s)")
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylabel("h (m$^2$/s)")
axs.flatten()[3].set_ylabel("Period (days)")
axs.flatten()[3].set_xlabel("time (yr)")

if L_reference == 11:
    axs.flatten()[3].set_ylim(top=3.e5)
if L_reference == 19:
    axs.flatten()[3].set_ylim(top=2.e0)
else:
    axs.flatten()[3].set_ylim(top=1.e1)
axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[0].legend(loc='best')

plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.png', bbox_inches='tight')
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('R$_{sink}$ (au)')
plt.ylabel('P$_{rot}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')

#============Not scaled version============================

plot_it = -1
end_time = 10000
L_scale= np.sqrt(2)/16
L_spec_scale = np.sqrt(2)/2
m_scale = 1/8
scaling_factor = [L_reference-8, L_reference-9, L_reference-10, L_reference-11]
comp_time = np.nan
comp_h = []
comp_period = []

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
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
    if np.isnan(comp_time):
        comp_time = yt.YTQuantity((sink_data[primary_ind]['time'] - form_time)[-1], 's')
    
    for sink_id in sink_data.keys():
        sink_it = sink_it + 1

        time = sink_data[sink_id]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        L_spec_tot = np.sqrt((sink_data[sink_id]['anglx']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['angly']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['anglz']/sink_data[sink_id]['mass'])**2)
        L_spec_tot = yt.YTArray(L_spec_tot, 'cm**2/s')
        plot_L_spec_tot = L_spec_tot.in_units('m**2/s')[:end_ind+1]
        
        #Sink_rot = plot_L_spec_tot.in_units('km**2/s') / r_sink[plot_it].in_units('km') #I think this is wrong
        #Sink_circum = 2*np.pi * r_sink[plot_it].in_units('km')
        #Sink_period = Sink_circum/Sink_rot
        
        L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        plot_L_tot = L_tot.in_units('kg*m**2/s')[:end_ind+1]
        
        mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
        plot_mass = mass.in_units('msun')[:end_ind+1]
        
        #===================================================
        #Calculate rotation Period
        #v_scale = plot_L_spec_tot.in_units('km**2/day')*(L_spec_scale**(scaling_factor[plot_it]))/r_sink_9.in_units('km')
        v_scale = plot_L_spec_tot.in_units('km**2/day')/r_scale.in_units('km')
        Sink_circum = 2*np.pi * r_scale.in_units('km')
        Sink_period = Sink_circum/v_scale
        
        #Sink_period = 0.8 * np.pi * plot_mass.in_units('g') * r_sink[plot_it].in_units('cm')**2 / plot_L_tot.in_units('g*cm**2/day')
        
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].plot(plot_time, plot_mass, linestyle=line_styles[plot_it], label="R$_{sink}$="+str(np.round(r_sink[plot_it], decimals=2))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].plot(plot_time, plot_mass, linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[1].plot(plot_time, plot_L_tot/1e45, linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[2].plot(plot_time, plot_L_spec_tot/1e15, linestyle=line_styles[plot_it], color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)")
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylabel("L (kg$\,$m$^2$/s)")
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylabel("h (m$^2$/s)")
axs.flatten()[3].set_ylabel("Period (days)")
axs.flatten()[3].set_xlabel("time (yr)")

.flatten()[3].set_ylim(top=1.e4)
axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[0].legend(loc='best')

plt.savefig('resolution_study_raw.png', bbox_inches='tight')
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('R$_{sink}$ (au)')
plt.ylabel('P$_{rot}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')
