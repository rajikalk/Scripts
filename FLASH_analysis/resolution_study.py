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
'''
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pickle_file in pickle_files:
    plot_it = plot_it + 1
    file = open(pickle_file, 'rb')
    sink_data, line_counter = pickle.load(file)
    file.close()
    L_ref = int(pickle_file.split('_')[-1].split('.')[0])
    
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
        r_sink = (1/(2**(L_ref-9)))*r_sink_9
        ang_vel = (5/2)*(plot_L_spec_tot.in_units('au**2/day')/(r_sink.in_units('au')**2))
        Sink_period = 1/ang_vel
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].plot(plot_time, plot_mass*(m_scale**(scaling_factor[plot_it])), linestyle=line_styles[plot_it], label="$R_{\mathrm{sink}}$="+str(np.round(r_sink[plot_it], decimals=1))+"au", color=colors[plot_it])
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
    axs.flatten()[3].set_ylim(top=5.e4)
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
plt.xlabel('$R_{\mathrm{sink}}$ (au)')
plt.ylabel('$T_{\mathrm{rot}}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')
'''
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

Mass_2500 = []
L_2500 = []
h_2500 = []
T_rot_2500 = []

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pickle_file in pickle_files:
    plot_it = plot_it + 1
    file = open(pickle_file, 'rb')
    sink_data, line_counter = pickle.load(file)
    file.close()
    L_ref = int(pickle_file.split('_')[-1].split('.')[0])
    
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
        r_sink_val = (1/(2**(L_ref-9)))*r_sink_9
        ang_vel = (5/2)*(plot_L_spec_tot.in_units('au**2/day')/(r_sink_val.in_units('au')**2))
        Sink_period = 1/ang_vel
        if sink_id == primary_ind:
            ind_2500 = np.argmin(abs(plot_time.value-2500))
            h_2500.append(plot_L_spec_tot[ind_2500])
            L_2500.append(plot_L_tot[ind_2500])
            Mass_2500.append(plot_mass[ind_2500])
            T_rot_2500.append(Sink_period[ind_2500])
        
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].plot(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1,  label="R_{\mathrm{sink}}$="+str(np.round(r_sink[plot_it], decimals=1))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].plot(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[1].plot(plot_time, plot_L_tot/1e45, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[2].plot(plot_time, plot_L_spec_tot/1e15, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.3)
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylabel("$L$ (10$^{45}$kg$\,$m$^2$/s)")
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylabel("$h$ (10$^{15}$m$^2$/s)")
axs.flatten()[3].set_ylabel("$T_{\mathrm{rot}}$ (days)", labelpad=-0.3)
axs.flatten()[3].set_xlabel("time (yr)", labelpad=-0.1)

axs.flatten()[3].set_ylim([3.e2, 5.e4])
axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[0].legend(loc='upper right')
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)

plt.savefig('resolution_study_raw.pdf', bbox_inches='tight', pad_inches=0.02)
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('$R_{\mathrm{sink}}$ (au)')
plt.ylabel('$T_{\mathrm{rot}}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')

#============Not scaled semilogy version============================

plot_it = -1
end_time = 10000
L_scale= np.sqrt(2)/16
L_spec_scale = np.sqrt(2)/2
m_scale = 1/8
scaling_factor = [L_reference-8, L_reference-9, L_reference-10, L_reference-11]
comp_time = np.nan
comp_h = []
comp_period = []

Mass_2500 = []
L_2500 = []
h_2500 = []
T_rot_2500 = []

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pickle_file in pickle_files:
    plot_it = plot_it + 1
    file = open(pickle_file, 'rb')
    sink_data, line_counter = pickle.load(file)
    file.close()
    L_ref = int(pickle_file.split('_')[-1].split('.')[0])
    
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
        r_sink_val = (1/(2**(L_ref-9)))*r_sink_9
        ang_vel = (5/2)*(plot_L_spec_tot.in_units('au**2/day')/(r_sink_val.in_units('au')**2))
        Sink_period = 1/ang_vel
        if sink_id == primary_ind:
            ind_2500 = np.argmin(abs(plot_time.value-2500))
            h_2500.append(plot_L_spec_tot[ind_2500])
            L_2500.append(plot_L_tot[ind_2500])
            Mass_2500.append(plot_mass[ind_2500])
            T_rot_2500.append(Sink_period[ind_2500])
        
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].semilogy(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, label="$R_{\mathrm{sink}}$="+str(np.round(r_sink[plot_it], decimals=1))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].semilogy(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[1].semilogy(plot_time, plot_L_tot, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[2].semilogy(plot_time, plot_L_spec_tot, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.3)
axs.flatten()[0].set_ylim([1.e-2, 1.e0])
axs.flatten()[1].set_ylabel("$L$ (kg$\,$m$^2$/s)", labelpad=-0.3)
axs.flatten()[1].set_ylim(bottom=1.e42)
axs.flatten()[2].set_ylabel("$h$ (m$^2$/s)", labelpad=-0.3)
axs.flatten()[2].set_ylim(bottom=1.e14)

axs.flatten()[3].set_ylabel("$T_{\mathrm{rot}}$ (days)", labelpad=-0.3)
axs.flatten()[3].set_xlabel("time (yr)", labelpad=-0.1)
axs.flatten()[3].set_ylim([3.e2, 5.e4])

axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[0].legend(loc='lower right')
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)

plt.savefig('resolution_study_raw_semilogy.pdf', bbox_inches='tight', pad_inches=0.02)
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('$R_{\mathrm{sink}}$ (au)')
plt.ylabel('$T_{\mathrm{rot}}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')

#============rough scale version============================

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
    L_ref = int(pickle_file.split('_')[-1].split('.')[0])
    
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
        #WORK OUT SCALED VERSION
    
        r_sink_val = (1/(2**(L_reference-9)))*r_sink_9
        ang_vel = (5/2)*((plot_L_spec_tot.in_units('au**2/day')*(1/(2**(L_reference - L_ref))))/(r_sink_val.in_units('au')**2))
        Sink_period = 1/ang_vel
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].plot(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, label="$R_{\mathrm{sink}}$="+str(np.round(r_sink[plot_it], decimals=1))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].plot(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[1].plot(plot_time, (plot_L_tot*(1/(2**(L_reference - L_ref))))/1.e44, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[2].plot(plot_time, (plot_L_spec_tot*(1/(2**(L_reference - L_ref))))/1.e14, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.3)
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].set_ylabel("$L$ (10$^{44}$kg$\,$m$^2$/s)")
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[2].set_ylabel("$h$ (10$^{14}$m$^2$/s)")
axs.flatten()[3].set_ylabel("$T_{\mathrm{rot}}$ (days)", labelpad=-0.3)
axs.flatten()[3].set_xlabel("time (yr)", labelpad=-0.1)

#axs.flatten()[3].set_ylim([5.e3, 5.e4])
axs.flatten()[3].set_ylim([3.e2, 5.e4])
axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[2].set_ylim(bottom=0)
axs.flatten()[0].legend(loc='upper right')
axs.flatten()[0].tick_params(axis='x', direction='in')
axs.flatten()[0].tick_params(axis='y', direction='in')
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in')
axs.flatten()[1].tick_params(axis='x', direction='in')
axs.flatten()[1].tick_params(axis='y', direction='in')
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in')
axs.flatten()[2].tick_params(axis='x', direction='in')
axs.flatten()[2].tick_params(axis='y', direction='in')
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in')
axs.flatten()[3].tick_params(axis='x', direction='in')
axs.flatten()[3].tick_params(axis='y', direction='in')
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in')

plt.savefig('rough_resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight', pad_inches=0.02)
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('$R_{\mathrm{sink}}$ (au)')
plt.ylabel('$T_{\mathrm{rot}}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')

#============rough scale semilogy version============================

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
    L_ref = int(pickle_file.split('_')[-1].split('.')[0])
    
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
        #WORK OUT SCALED VERSION
    
        r_sink_val = (1/(2**(L_reference-9)))*r_sink_9
        ang_vel = (5/2)*((plot_L_spec_tot.in_units('au**2/day')*(1/(2**(L_reference - L_ref))))/(r_sink_val.in_units('au')**2))
        Sink_period = 1/ang_vel
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].semilogy(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, label="$R_{\mathrm{sink}}$="+str(np.round(r_sink[plot_it], decimals=1))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].semilogy(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[1].semilogy(plot_time, (plot_L_tot*(1/(2**(L_reference - L_ref))))/1.e44, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[2].semilogy(plot_time, (plot_L_spec_tot*(1/(2**(L_reference - L_ref))))/1.e14, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.3)
axs.flatten()[0].set_ylim([1.e-2, 1.e0])
axs.flatten()[1].set_ylabel("$L$ (10$^{44}$kg$\,$m$^2$/s)")
axs.flatten()[1].set_ylim([5.e-3, 5.e0])
axs.flatten()[2].set_ylabel("$h$ (10$^{14}$m$^2$/s)")
axs.flatten()[2].set_ylim([1.e-1, 5.e0])
axs.flatten()[3].set_ylabel("$T_{\mathrm{rot}}$ (days)", labelpad=-0.3)
axs.flatten()[3].set_xlabel("time (yr)", labelpad=-0.1)

#axs.flatten()[3].set_ylim([5.e3, 5.e4])
axs.flatten()[3].set_ylim([3.e2, 5.e4])
axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[0].legend(loc='upper right')
axs.flatten()[0].tick_params(axis='x', direction='in')
axs.flatten()[0].tick_params(axis='y', direction='in')
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in')
axs.flatten()[1].tick_params(axis='x', direction='in')
axs.flatten()[1].tick_params(axis='y', direction='in')
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in')
axs.flatten()[2].tick_params(axis='x', direction='in')
axs.flatten()[2].tick_params(axis='y', direction='in')
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in')
axs.flatten()[3].tick_params(axis='x', direction='in')
axs.flatten()[3].tick_params(axis='y', direction='in')
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in')

plt.savefig('rough_resolution_study_scaled_L'+str(L_reference)+'_semilogy.pdf', bbox_inches='tight', pad_inches=0.02)
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('$R_{\mathrm{sink}}$ (au)')
plt.ylabel('$T_{\mathrm{rot}}$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')

#=============Ploting values (2500) vs resolution
r_sink = [2*r_sink_9, r_sink_9, (1/(2))*r_sink_9]#, (1/(2**2))*r_sink_9]

plt.clf()
plt.cla()
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].loglog(r_sink, Mass_2500[:-1], lw=1, label='Mass')
axs.flatten()[0].scatter(r_sink, Mass_2500[:-1], label='Mass')
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)")

axs.flatten()[1].loglog(r_sink, L_2500[:-1], lw=1, label='Angular momentum')
axs.flatten()[1].scatter(r_sink, L_2500[:-1], label='Angular momentum')
axs.flatten()[1].set_ylabel("$L$ (10$^{44}$kg$\,$m$^2$/s)")

axs.flatten()[2].loglog(r_sink, h_2500[:-1], lw=1, label='Specific Angular momentum')
axs.flatten()[2].scatter(r_sink, h_2500[:-1], label='Specific Angular momentum')
axs.flatten()[2].set_ylabel("$h$ (10$^{14}$m$^2$/s)")

axs.flatten()[3].loglog(r_sink, T_rot_2500[:-1], lw=1, label='Rotation period')
axs.flatten()[3].scatter(r_sink, T_rot_2500[:-1], label='Rotation period')
axs.flatten()[3].set_ylabel("$T_{\mathrm{rot}}$ (days)")

axs.flatten()[3].set_xlabel('Sink particle accretion radius (au)')

plt.savefig('quantity_vs_r_sink.pdf', bbox_inches='tight', pad_inches=0.02)

#=============Ploting values (2500) vs resolution, with Lref 11
r_sink = [2*r_sink_9, r_sink_9, (1/(2))*r_sink_9, (1/(2**2))*r_sink_9]

plt.cla()
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].loglog(r_sink, Mass_2500, lw=1, label='Mass')
axs.flatten()[0].scatter(r_sink, Mass_2500, label='Mass')
axs.flatten()[0].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.01)

axs.flatten()[1].loglog(r_sink, L_2500, lw=1, label='Angular momentum')
axs.flatten()[1].scatter(r_sink, L_2500, label='Angular momentum')
axs.flatten()[1].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[1].set_ylabel("$L$ (kg$\,$m$^2$/s)", labelpad=-0.01)

axs.flatten()[2].loglog(r_sink, h_2500, lw=1, label='Specific Angular momentum')
axs.flatten()[2].scatter(r_sink, h_2500, label='Specific Angular momentum')
axs.flatten()[2].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[2].set_ylabel("$h$ (m$^2$/s)", labelpad=-0.01)

axs.flatten()[3].loglog(r_sink, T_rot_2500, lw=1, label='Rotation period')
axs.flatten()[3].scatter(r_sink, T_rot_2500, label='Rotation period')
axs.flatten()[3].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[3].set_ylabel("$T_{\mathrm{rot}}$ (days)", labelpad=-0.01)

axs.flatten()[3].set_xlabel('$R_\mathrm{sink}$ (au)', labelpad=-0.02)

from scipy.optimize import curve_fit

x_fit = np.linspace(np.log10(R_proto.in_units('au').value), np.log10(r_sink[0].value))
#x_fit = np.linspace(np.log10(r_sink[-1].value), np.log10(r_sink[0].value))
x = 10**x_fit
def power(x, a, b, c):
    return a + b*x**c
def power(x, b, c):
    return b*x**c

c_guess = (np.max(np.log10(Mass_2500[:-1])) - np.min(np.log10(Mass_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(Mass_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], Mass_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]))
axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[0].text(1.e-2, Mass_2500[1], '$y=('+str(np.round(popt[0], decimals=2))+')x^{'+str(np.round(popt[1], decimals=2))+'}$')

c_guess = (np.max(np.log10(L_2500[:-1])) - np.min(np.log10(L_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(L_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], L_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]))
'''
b_err = 0.1
c_err = 0.1
lower_guess = [-1, 10**(np.log10(b_guess)-b_err), 1.1422525230062428]# c_guess-c_err]
higher_guess = [1, 10**(np.log10(b_guess)+b_err), 1.1422525230062428]# c_guess+c_err]
plt.clf()
plt.scatter(r_sink, L_2500, label='L')
plt.loglog(x, power(x, *initial_guess), label="initial")
plt.loglog(x, power(x, *lower_guess), label="lower")
plt.loglog(x, power(x, *higher_guess), label="higher")
plt.loglog(x, power(x, *popt), ls='--')
plt.savefig("testing_fits.png")
'''

axs.flatten()[1].loglog(x, power(x, *popt), lw=1)
axs.flatten()[1].text(1.e-2, L_2500[1], '$y=('+'{:.2e}'.format(popt[0])+')x^{'+str(np.round(popt[1], decimals=2))+'}$')


c_guess = (np.max(np.log10(h_2500[:-1])) - np.min(np.log10(h_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(h_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], h_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]))
'''
b_err = 0.1
c_err = 0.1
lower_guess = [-1, 10**(np.log10(b_guess)-b_err), c_guess-c_err]# c_guess-c_err]
higher_guess = [1, 10**(np.log10(b_guess)+b_err), c_guess+c_err]# c_guess+c_err]
plt.clf()
plt.scatter(r_sink, h_2500, label='L')
plt.loglog(x, power(x, *initial_guess), label="initial")
plt.loglog(x, power(x, *lower_guess), label="lower")
plt.loglog(x, power(x, *higher_guess), label="higher")
plt.loglog(x, power(x, *popt), ls='--')
plt.savefig("testing_fits.png")
'''

axs.flatten()[2].loglog(x, power(x, *popt), lw=1)
axs.flatten()[2].text(1.e-2, h_2500[1], '$y=('+'{:.2e}'.format(popt[0])+')x^{'+str(np.round(popt[1], decimals=2))+'}$')

c_guess = (np.max(np.log10(T_rot_2500[:-1])) - np.min(np.log10(T_rot_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(T_rot_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], T_rot_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]))
axs.flatten()[3].loglog(x, power(x, *popt), lw=1)
axs.flatten()[3].text(1.e-2, T_rot_2500[1], '$y=('+str(np.round(popt[0], decimals=2))+') x^{'+str(np.round(popt[1], decimals=2))+'}$')

axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)

plt.savefig('quantity_vs_r_sink_L11_extrapolation.pdf', bbox_inches='tight', pad_inches=0.02)
